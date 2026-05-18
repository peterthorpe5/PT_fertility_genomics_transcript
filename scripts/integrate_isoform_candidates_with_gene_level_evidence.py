#!/usr/bin/env python3
"""
Integrate GTEx testis isoform candidates with gene-level evidence layers.

This script joins transcript-level candidate isoforms from the GTEx testis
isoform workflow to gene-level fertility, sperm RNA, sperm proteomics,
biochemical accessibility, and Open Targets tractability evidence.

The intended use is to move from transcript-only prioritisation towards a
reviewable evidence table in which each candidate isoform can be assessed in
biological and translational context.

Main outputs
------------
<prefix>.isoform_candidates_with_gene_level_evidence.tsv
    Full integrated candidate table.

<prefix>.tier1_protein_coding_rescue_with_gene_evidence.tsv
    Candidate rows that are tier 1 protein-coding rescue candidates when the
    priority tier is available, or broad-gene protein-coding rescue candidates
    when only raw isoform flags are available.

<prefix>.top_sperm_supported_druggable_isoform_rescue_candidates.tsv
    Review-focused table requiring isoform rescue, sperm/protein evidence, and
    at least one tractability or accessibility signal.

<prefix>.selected_project_genes_isoform_gene_evidence.tsv
    Rows for selected project genes, by default AFG2B, CFAP99, SLC16A7 and
    ABCG4.

<prefix>.summary_by_integrated_evidence_class.tsv
    Counts summarised by integrated evidence class.

<prefix>.merge_diagnostics.tsv
    Join diagnostics for each gene-level evidence source.

Formatted Excel copies of every tabular result are written by default, alongside a combined multi-sheet review workbook.
"""

from __future__ import annotations

import argparse
import logging
import re
import sys
import time
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class Config:
    """
    Workflow configuration.

    Attributes
    ----------
    isoform_candidates_tsv
        Transcript-level candidate or prioritised isoform table.
    biochem_xlsx
        Excel workbook containing the gene-level biochemical/evidence table.
    biochem_sheet
        Sheet name in the biochemical workbook.
    gene_priority_xlsx
        Optional Excel workbook containing gene-level prioritisation and Open
        Targets evidence.
    gene_priority_sheet
        Sheet name in the gene-priority workbook.
    tractability_tsv
        Optional TSV/TSV.GZ/TSV.ZIP file containing full gene-universe Open
        Targets tractability annotations.
    out_dir
        Output directory.
    out_prefix
        Prefix used for all output files.
    selected_genes
        Project genes to extract into a separate review file.
    top_n_review
        Number of top genes to write to convenience gene-list files.
    write_excel_outputs
        Whether to write formatted Excel copies of the main review outputs.
    log_path
        Optional log file.
    log_level
        Logging level.
    """

    isoform_candidates_tsv: Path
    biochem_xlsx: Path
    biochem_sheet: str
    gene_priority_xlsx: Optional[Path]
    gene_priority_sheet: str
    tractability_tsv: Optional[Path]
    out_dir: Path
    out_prefix: str
    selected_genes: Tuple[str, ...]
    top_n_review: int
    write_excel_outputs: bool
    log_path: Optional[Path]
    log_level: str


@contextmanager
def timed(*, logger: logging.Logger, label: str) -> Iterator[None]:
    """
    Log start, completion, and elapsed time for a block.

    Parameters
    ----------
    logger
        Logger instance.
    label
        Human-readable label for the timed block.

    Yields
    ------
    Iterator[None]
        Context manager body.
    """
    start = time.time()
    logger.info("%s: started", label)
    try:
        yield
    except Exception:
        elapsed = time.time() - start
        logger.exception("%s: failed after %.1f seconds", label, elapsed)
        raise
    elapsed = time.time() - start
    logger.info("%s: completed in %.1f seconds", label, elapsed)


def setup_logging(*, log_level: str, log_path: Optional[Path]) -> logging.Logger:
    """
    Configure console and optional file logging.

    Parameters
    ----------
    log_level
        Logging level name.
    log_path
        Optional log-file path.

    Returns
    -------
    logging.Logger
        Configured logger.
    """
    logger = logging.getLogger("integrate_isoform_gene_evidence")
    logger.setLevel(getattr(logging, log_level.upper(), logging.INFO))
    logger.handlers = []
    logger.propagate = False

    formatter = logging.Formatter(
        fmt="%(asctime)s | %(levelname)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    stream_handler = logging.StreamHandler(stream=sys.stdout)
    stream_handler.setFormatter(formatter)
    stream_handler.setLevel(getattr(logging, log_level.upper(), logging.INFO))
    logger.addHandler(stream_handler)

    if log_path is not None:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(filename=log_path)
        file_handler.setFormatter(formatter)
        file_handler.setLevel(getattr(logging, log_level.upper(), logging.INFO))
        logger.addHandler(file_handler)

    return logger


def parse_args(argv: Optional[Sequence[str]] = None) -> Config:
    """
    Parse command-line arguments.

    Parameters
    ----------
    argv
        Optional argument list. If omitted, arguments are read from sys.argv.

    Returns
    -------
    Config
        Parsed workflow configuration.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Integrate GTEx testis isoform candidates with gene-level sperm, "
            "biochemical and tractability evidence."
        )
    )
    parser.add_argument(
        "--isoform_candidates_tsv",
        required=True,
        type=Path,
        help="Prioritised isoform candidate table, TSV or TSV.GZ.",
    )
    parser.add_argument(
        "--biochem_xlsx",
        required=True,
        type=Path,
        help="Biochemical/gene evidence Excel workbook.",
    )
    parser.add_argument(
        "--biochem_sheet",
        default="Genes_Master",
        help="Sheet containing the gene-level biochemical table.",
    )
    parser.add_argument(
        "--gene_priority_xlsx",
        default=None,
        type=Path,
        help=(
            "Optional gene-level priority/Open Targets Excel workbook. "
            "The ranked sheet is used by default."
        ),
    )
    parser.add_argument(
        "--gene_priority_sheet",
        default="ranked",
        help="Sheet name in the gene-level priority workbook.",
    )
    parser.add_argument(
        "--tractability_tsv",
        default=None,
        type=Path,
        help=(
            "Optional full-universe Open Targets tractability TSV, TSV.GZ, "
            "or single-file TSV.ZIP."
        ),
    )
    parser.add_argument(
        "--out_dir",
        required=True,
        type=Path,
        help="Output directory.",
    )
    parser.add_argument(
        "--out_prefix",
        default="gtex_isoform_candidates_gene_level_evidence",
        help="Output filename prefix.",
    )
    parser.add_argument(
        "--selected_genes",
        nargs="+",
        default=["AFG2B", "CFAP99", "SLC16A7", "ABCG4"],
        help="Selected project genes to extract for separate review.",
    )
    parser.add_argument(
        "--top_n_review",
        default=50,
        type=int,
        help="Number of top review genes to write to text gene lists.",
    )
    parser.add_argument(
        "--no_write_excel_outputs",
        action="store_true",
        help="Disable formatted Excel output. Excel output is on by default.",
    )
    parser.add_argument(
        "--log_path",
        default=None,
        type=Path,
        help="Optional log-file path.",
    )
    parser.add_argument(
        "--log_level",
        default="INFO",
        help="Logging level. Default: INFO.",
    )

    args = parser.parse_args(argv)
    return Config(
        isoform_candidates_tsv=args.isoform_candidates_tsv,
        biochem_xlsx=args.biochem_xlsx,
        biochem_sheet=args.biochem_sheet,
        gene_priority_xlsx=args.gene_priority_xlsx,
        gene_priority_sheet=args.gene_priority_sheet,
        tractability_tsv=args.tractability_tsv,
        out_dir=args.out_dir,
        out_prefix=args.out_prefix,
        selected_genes=tuple(args.selected_genes),
        top_n_review=args.top_n_review,
        write_excel_outputs=not args.no_write_excel_outputs,
        log_path=args.log_path,
        log_level=args.log_level,
    )


def read_table(*, path: Path, logger: logging.Logger) -> pd.DataFrame:
    """
    Read a TSV, TSV.GZ or single-file TSV.ZIP table.

    Parameters
    ----------
    path
        Input table path.
    logger
        Logger instance.

    Returns
    -------
    pandas.DataFrame
        Loaded table.
    """
    logger.info("Reading table: %s", path)
    compression = "zip" if str(path).endswith(".zip") else "infer"
    dataframe = pd.read_csv(
        filepath_or_buffer=path,
        sep="\t",
        compression=compression,
        low_memory=False,
    )
    logger.info("Loaded %d rows x %d columns from %s", *dataframe.shape, path)
    return dataframe


def read_excel_sheet(
    *, path: Path, sheet_name: str, logger: logging.Logger
) -> pd.DataFrame:
    """
    Read an Excel worksheet.

    Parameters
    ----------
    path
        Excel workbook path.
    sheet_name
        Worksheet name.
    logger
        Logger instance.

    Returns
    -------
    pandas.DataFrame
        Loaded worksheet.
    """
    logger.info("Reading Excel sheet %s from %s", sheet_name, path)
    dataframe = pd.read_excel(
        io=path,
        sheet_name=sheet_name,
    )
    logger.info(
        "Loaded %d rows x %d columns from %s:%s",
        dataframe.shape[0],
        dataframe.shape[1],
        path,
        sheet_name,
    )
    return dataframe


def strip_ensembl_version(*, identifier: object) -> str:
    """
    Remove an Ensembl dot-version suffix.

    Parameters
    ----------
    identifier
        Identifier value.

    Returns
    -------
    str
        Version-stripped identifier, or empty string for missing values.
    """
    if pd.isna(identifier):
        return ""
    value = str(identifier).strip()
    if not value or value.lower() == "nan":
        return ""
    return value.split(".", maxsplit=1)[0]


def normalise_gene_symbol(*, value: object) -> str:
    """
    Normalise a gene symbol for robust joining.

    Parameters
    ----------
    value
        Gene symbol-like value.

    Returns
    -------
    str
        Upper-case, stripped gene symbol.
    """
    if pd.isna(value):
        return ""
    symbol = str(value).strip()
    if not symbol or symbol.lower() == "nan":
        return ""
    return symbol.upper()


def add_join_keys(
    *,
    dataframe: pd.DataFrame,
    symbol_candidates: Sequence[str],
    gene_id_candidates: Sequence[str],
    logger: Optional[logging.Logger] = None,
) -> pd.DataFrame:
    """
    Add normalised gene-symbol and Ensembl-gene join keys.

    Parameters
    ----------
    dataframe
        Input table.
    symbol_candidates
        Candidate gene-symbol columns, checked in order.
    gene_id_candidates
        Candidate gene-ID columns, checked in order.
    logger
        Optional logger.

    Returns
    -------
    pandas.DataFrame
        Copy of input table with ``join_gene_symbol`` and ``join_gene_id``.
    """
    output = dataframe.copy()
    symbol_col = first_present_column(
        columns=output.columns,
        candidates=symbol_candidates,
        required=False,
    )
    gene_id_col = first_present_column(
        columns=output.columns,
        candidates=gene_id_candidates,
        required=False,
    )

    if symbol_col is not None:
        output["join_gene_symbol"] = output[symbol_col].map(
            lambda value: normalise_gene_symbol(value=value)
        )
    else:
        output["join_gene_symbol"] = ""

    if gene_id_col is not None:
        output["join_gene_id"] = output[gene_id_col].map(
            lambda value: strip_ensembl_version(identifier=value)
        )
    else:
        output["join_gene_id"] = ""

    if logger is not None:
        logger.info(
            "Join keys: symbol_col=%s, gene_id_col=%s, non-empty symbols=%d, "
            "non-empty gene IDs=%d",
            symbol_col,
            gene_id_col,
            int((output["join_gene_symbol"] != "").sum()),
            int((output["join_gene_id"] != "").sum()),
        )

    return output


def first_present_column(
    *, columns: Iterable[str], candidates: Sequence[str], required: bool = True
) -> Optional[str]:
    """
    Return the first candidate column present in a column list.

    Parameters
    ----------
    columns
        Available column names.
    candidates
        Candidate column names to search for.
    required
        Whether to raise an error if no candidate is found.

    Returns
    -------
    Optional[str]
        First present column, or None when not required.
    """
    column_set = set(columns)
    for candidate in candidates:
        if candidate in column_set:
            return candidate
    if required:
        raise ValueError(f"Missing required columns. Tried: {tuple(candidates)}")
    return None


def make_prefixed_evidence_table(
    *,
    dataframe: pd.DataFrame,
    prefix: str,
    keep_columns: Optional[Sequence[str]] = None,
) -> pd.DataFrame:
    """
    Prefix evidence columns while preserving join keys.

    Parameters
    ----------
    dataframe
        Evidence table with join keys already added.
    prefix
        Prefix to add to evidence columns.
    keep_columns
        Optional evidence columns to keep. Missing columns are ignored. If
        omitted, all non-join columns are retained.

    Returns
    -------
    pandas.DataFrame
        Deduplicated evidence table with prefixed evidence columns.
    """
    join_cols = ["join_gene_symbol", "join_gene_id"]
    if keep_columns is None:
        evidence_cols = [column for column in dataframe.columns if column not in join_cols]
    else:
        evidence_cols = [column for column in keep_columns if column in dataframe.columns]

    output = dataframe[join_cols + evidence_cols].copy()
    rename_map = {
        column: f"{prefix}{column}"
        for column in evidence_cols
        if not column.startswith(prefix)
    }
    output = output.rename(columns=rename_map)
    output = output.drop_duplicates(subset=join_cols, keep="first")
    return output


def merge_by_gene_keys(
    *,
    base: pd.DataFrame,
    evidence: pd.DataFrame,
    source_name: str,
    logger: logging.Logger,
) -> Tuple[pd.DataFrame, Dict[str, object]]:
    """
    Merge evidence into base using symbol first, then gene ID fallback.

    Parameters
    ----------
    base
        Base isoform table with join keys.
    evidence
        Evidence table with join keys.
    source_name
        Source label for diagnostics.
    logger
        Logger instance.

    Returns
    -------
    tuple[pandas.DataFrame, dict[str, object]]
        Merged table and diagnostic record.
    """
    base_count = base.shape[0]
    evidence_cols = [
        column
        for column in evidence.columns
        if column not in {"join_gene_symbol", "join_gene_id"}
    ]
    if not evidence_cols:
        logger.warning("Evidence source %s has no evidence columns", source_name)
        return base.copy(), {
            "source": source_name,
            "base_rows": base_count,
            "matched_by_symbol": 0,
            "matched_by_gene_id_fallback": 0,
            "unmatched_rows": base_count,
            "evidence_columns_added": 0,
        }

    symbol_evidence = evidence.loc[evidence["join_gene_symbol"] != ""].copy()
    symbol_evidence = symbol_evidence.drop_duplicates(
        subset=["join_gene_symbol"], keep="first"
    )
    merged = base.merge(
        symbol_evidence.drop(columns=["join_gene_id"]),
        on="join_gene_symbol",
        how="left",
        validate="many_to_one",
    )
    symbol_match_mask = merged[evidence_cols].notna().any(axis=1)

    gene_id_evidence = evidence.loc[evidence["join_gene_id"] != ""].copy()
    gene_id_evidence = gene_id_evidence.drop_duplicates(
        subset=["join_gene_id"], keep="first"
    )
    fallback = base.loc[~symbol_match_mask].merge(
        gene_id_evidence.drop(columns=["join_gene_symbol"]),
        on="join_gene_id",
        how="left",
        validate="many_to_one",
    )

    fallback_match_mask = fallback[evidence_cols].notna().any(axis=1)
    for column in evidence_cols:
        if column in fallback.columns:
            merged.loc[~symbol_match_mask, column] = fallback[column].values

    matched_by_symbol = int(symbol_match_mask.sum())
    matched_by_gene_id = int(fallback_match_mask.sum())
    unmatched = base_count - matched_by_symbol - matched_by_gene_id

    logger.info(
        "%s merge: %d rows matched by symbol, %d by gene ID fallback, %d unmatched",
        source_name,
        matched_by_symbol,
        matched_by_gene_id,
        unmatched,
    )

    diagnostics = {
        "source": source_name,
        "base_rows": base_count,
        "matched_by_symbol": matched_by_symbol,
        "matched_by_gene_id_fallback": matched_by_gene_id,
        "unmatched_rows": unmatched,
        "evidence_columns_added": len(evidence_cols),
    }
    return merged, diagnostics


def truthy_series(*, dataframe: pd.DataFrame, columns: Sequence[str]) -> pd.Series:
    """
    Collapse multiple possible boolean columns into one truthy series.

    Parameters
    ----------
    dataframe
        Input table.
    columns
        Candidate columns to inspect.

    Returns
    -------
    pandas.Series
        Boolean series that is true if any available column is truthy.
    """
    if dataframe.empty:
        return pd.Series(dtype=bool)
    result = pd.Series(False, index=dataframe.index)
    for column in columns:
        if column not in dataframe.columns:
            continue
        result = result | dataframe[column].map(is_truthy_value).fillna(False)
    return result


def is_truthy_value(value: object) -> bool:
    """
    Interpret common truthy values robustly.

    Parameters
    ----------
    value
        Value to parse.

    Returns
    -------
    bool
        Parsed boolean value.
    """
    if pd.isna(value):
        return False
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (int, float, np.integer, np.floating)):
        return bool(value)
    text = str(value).strip().lower()
    if text in {"true", "t", "yes", "y", "1", "present", "detected", "strong"}:
        return True
    return False


def contains_any_series(
    *, dataframe: pd.DataFrame, columns: Sequence[str], patterns: Sequence[str]
) -> pd.Series:
    """
    Return whether any candidate text column contains one of the patterns.

    Parameters
    ----------
    dataframe
        Input table.
    columns
        Candidate text columns.
    patterns
        Case-insensitive text patterns.

    Returns
    -------
    pandas.Series
        Boolean series.
    """
    result = pd.Series(False, index=dataframe.index)
    regex = re.compile("|".join(re.escape(pattern) for pattern in patterns), re.I)
    for column in columns:
        if column not in dataframe.columns:
            continue
        result = result | dataframe[column].astype(str).str.contains(regex, na=False)
    return result


def first_numeric_series(
    *, dataframe: pd.DataFrame, columns: Sequence[str], default: float = 0.0
) -> pd.Series:
    """
    Return the first available numeric column as a Series.

    Parameters
    ----------
    dataframe
        Input table.
    columns
        Candidate columns.
    default
        Default value if no column is available.

    Returns
    -------
    pandas.Series
        Numeric series.
    """
    for column in columns:
        if column in dataframe.columns:
            return pd.to_numeric(dataframe[column], errors="coerce").fillna(default)
    return pd.Series(default, index=dataframe.index, dtype=float)


def derive_integrated_evidence(*, dataframe: pd.DataFrame) -> pd.DataFrame:
    """
    Add harmonised evidence flags, classes and review score.

    Parameters
    ----------
    dataframe
        Integrated candidate table.

    Returns
    -------
    pandas.DataFrame
        Table with additional harmonised evidence columns.
    """
    output = dataframe.copy()

    output["gene_level_has_sperm_rna"] = truthy_series(
        dataframe=output,
        columns=(
            "biochem_sperm_present_any",
            "geneprio_sperm_present_any",
            "geneprio_sperm_rnaseq_present",
            "geneprio_sperm_rnaseq_present",
        ),
    )
    output["gene_level_has_any_proteomics"] = truthy_series(
        dataframe=output,
        columns=(
            "biochem_prot_present_any",
            "biochem_public_prot_present_any",
            "geneprio_prot_present_any",
            "geneprio_public_prot_present_any",
            "geneprio_proteomics_present_any_source",
            "geneprio_prot_any_detected_or_strong",
        ),
    )
    output["gene_level_has_strong_public_proteomics"] = contains_any_series(
        dataframe=output,
        columns=(
            "biochem_public_proteomics_evidence_level",
            "biochem_proteomics_evidence_level",
            "geneprio_public_proteomics_evidence_level",
            "geneprio_proteomics_evidence_level",
        ),
        patterns=("strong",),
    )
    output["gene_level_has_known_or_clinical_fertility_support"] = truthy_series(
        dataframe=output,
        columns=(
            "biochem_in_hpo_gene_set",
            "biochem_clinvar_best_pathogenic_present",
            "biochem_clinvar_hc_pathogenic_present",
            "biochem_in_literature_fertility_set",
            "geneprio_in_hpo_gene_set",
            "geneprio_clinvar_best_pathogenic_present_x",
            "geneprio_clinvar_hc_pathogenic_present_x",
            "geneprio_in_literature_fertility_set",
        ),
    )

    output["gene_level_ot_small_molecule"] = truthy_series(
        dataframe=output,
        columns=(
            "geneprio_ot_any_small_molecule_tractable",
            "tract_ot_any_small_molecule_tractable",
        ),
    )
    output["gene_level_ot_antibody"] = truthy_series(
        dataframe=output,
        columns=(
            "geneprio_ot_any_antibody_tractable",
            "tract_ot_any_antibody_tractable",
        ),
    )
    output["gene_level_ot_protac"] = truthy_series(
        dataframe=output,
        columns=(
            "geneprio_ot_any_protac_tractable",
            "tract_ot_any_protac_tractable",
        ),
    )
    output["gene_level_has_any_ot_tractability"] = (
        output["gene_level_ot_small_molecule"]
        | output["gene_level_ot_antibody"]
        | output["gene_level_ot_protac"]
        | truthy_series(dataframe=output, columns=("geneprio_ot_any_tractable",))
    )

    output["gene_level_is_cell_surface"] = truthy_series(
        dataframe=output,
        columns=("biochem_is_cell_surface_candidate",),
    )
    output["gene_level_is_secreted"] = truthy_series(
        dataframe=output,
        columns=("biochem_is_secreted",),
    )
    output["gene_level_is_membrane"] = truthy_series(
        dataframe=output,
        columns=("biochem_is_membrane", "biochem_has_transmembrane"),
    )
    output["gene_level_has_accessibility_signal"] = (
        output["gene_level_is_cell_surface"]
        | output["gene_level_is_secreted"]
        | output["gene_level_is_membrane"]
    )

    output["gene_level_candidate_druggable_sperm_protein"] = truthy_series(
        dataframe=output,
        columns=("geneprio_candidate_druggable_sperm_protein",),
    )

    output["isoform_is_rescue_candidate"] = truthy_series(
        dataframe=output,
        columns=("is_broad_gene_isoform_rescue_candidate",),
    ) | contains_any_series(
        dataframe=output,
        columns=("priority_tier",),
        patterns=("rescue",),
    )
    output["isoform_is_protein_coding"] = truthy_series(
        dataframe=output,
        columns=("gencode_is_protein_coding_transcript", "gencode_has_cds"),
    ) | contains_any_series(
        dataframe=output,
        columns=("gencode_transcript_type",),
        patterns=("protein_coding",),
    )

    output["gene_level_sperm_or_protein_supported"] = (
        output["gene_level_has_sperm_rna"] | output["gene_level_has_any_proteomics"]
    )
    output["gene_level_druggable_or_accessible"] = (
        output["gene_level_has_any_ot_tractability"]
        | output["gene_level_has_accessibility_signal"]
        | output["gene_level_candidate_druggable_sperm_protein"]
    )

    output["integrated_evidence_class"] = classify_integrated_rows(
        dataframe=output
    )
    output["gene_evidence_score"] = calculate_gene_evidence_score(dataframe=output)
    isoform_score = first_numeric_series(
        dataframe=output,
        columns=("priority_score", "isoform_priority_score"),
        default=0.0,
    )
    output["integrated_review_score"] = isoform_score + output["gene_evidence_score"]
    output["integrated_score_components"] = make_score_component_strings(
        dataframe=output
    )
    return output


def classify_integrated_rows(*, dataframe: pd.DataFrame) -> pd.Series:
    """
    Classify integrated candidate rows into review-focused evidence classes.

    Parameters
    ----------
    dataframe
        Candidate table with harmonised evidence flags.

    Returns
    -------
    pandas.Series
        Integrated evidence class per row.
    """
    classes = pd.Series("other_isoform_candidate", index=dataframe.index)

    full_support = (
        dataframe["isoform_is_rescue_candidate"]
        & dataframe["isoform_is_protein_coding"]
        & dataframe["gene_level_sperm_or_protein_supported"]
        & dataframe["gene_level_druggable_or_accessible"]
    )
    sperm_only = (
        dataframe["isoform_is_rescue_candidate"]
        & dataframe["isoform_is_protein_coding"]
        & dataframe["gene_level_sperm_or_protein_supported"]
        & ~dataframe["gene_level_druggable_or_accessible"]
    )
    drug_only = (
        dataframe["isoform_is_rescue_candidate"]
        & dataframe["isoform_is_protein_coding"]
        & ~dataframe["gene_level_sperm_or_protein_supported"]
        & dataframe["gene_level_druggable_or_accessible"]
    )
    rescue_limited = (
        dataframe["isoform_is_rescue_candidate"]
        & dataframe["isoform_is_protein_coding"]
        & ~full_support
        & ~sperm_only
        & ~drug_only
    )
    non_coding_rescue = (
        dataframe["isoform_is_rescue_candidate"]
        & ~dataframe["isoform_is_protein_coding"]
    )

    classes.loc[non_coding_rescue] = "non_coding_or_uncertain_isoform_rescue"
    classes.loc[rescue_limited] = "protein_coding_isoform_rescue_limited_gene_support"
    classes.loc[drug_only] = "protein_coding_isoform_rescue_druggable_no_sperm_support"
    classes.loc[sperm_only] = "protein_coding_isoform_rescue_sperm_supported"
    classes.loc[full_support] = (
        "protein_coding_isoform_rescue_sperm_supported_druggable"
    )
    return classes


def calculate_gene_evidence_score(*, dataframe: pd.DataFrame) -> pd.Series:
    """
    Calculate a transparent gene-level evidence score.

    Parameters
    ----------
    dataframe
        Candidate table with harmonised evidence flags.

    Returns
    -------
    pandas.Series
        Gene-level evidence score.
    """
    score = pd.Series(0.0, index=dataframe.index)
    score += dataframe["gene_level_has_sperm_rna"].astype(float) * 3.0
    score += dataframe["gene_level_has_any_proteomics"].astype(float) * 4.0
    score += dataframe["gene_level_has_strong_public_proteomics"].astype(float) * 2.0
    score += dataframe["gene_level_has_known_or_clinical_fertility_support"].astype(float) * 3.0
    score += dataframe["gene_level_ot_small_molecule"].astype(float) * 4.0
    score += dataframe["gene_level_ot_antibody"].astype(float) * 3.0
    score += dataframe["gene_level_ot_protac"].astype(float) * 3.0
    score += dataframe["gene_level_is_cell_surface"].astype(float) * 3.0
    score += dataframe["gene_level_is_secreted"].astype(float) * 2.0
    score += dataframe["gene_level_is_membrane"].astype(float) * 2.0
    score += dataframe["gene_level_candidate_druggable_sperm_protein"].astype(float) * 5.0
    return score


def make_score_component_strings(*, dataframe: pd.DataFrame) -> pd.Series:
    """
    Create readable score-component strings.

    Parameters
    ----------
    dataframe
        Candidate table with harmonised evidence flags.

    Returns
    -------
    pandas.Series
        Semicolon-separated score component strings.
    """
    component_map = [
        ("gene_level_has_sperm_rna", "sperm_rna:+3"),
        ("gene_level_has_any_proteomics", "proteomics:+4"),
        ("gene_level_has_strong_public_proteomics", "strong_proteomics:+2"),
        ("gene_level_has_known_or_clinical_fertility_support", "known_clinical:+3"),
        ("gene_level_ot_small_molecule", "small_molecule:+4"),
        ("gene_level_ot_antibody", "antibody:+3"),
        ("gene_level_ot_protac", "protac:+3"),
        ("gene_level_is_cell_surface", "cell_surface:+3"),
        ("gene_level_is_secreted", "secreted:+2"),
        ("gene_level_is_membrane", "membrane:+2"),
        ("gene_level_candidate_druggable_sperm_protein", "candidate_druggable:+5"),
    ]
    values = []
    for _, row in dataframe.iterrows():
        components = [
            label
            for column, label in component_map
            if column in dataframe.columns and bool(row[column])
        ]
        values.append(";".join(components))
    return pd.Series(values, index=dataframe.index)


def select_outputs(
    *, dataframe: pd.DataFrame, selected_genes: Sequence[str]
) -> Dict[str, pd.DataFrame]:
    """
    Select review-focused output tables.

    Parameters
    ----------
    dataframe
        Full integrated table.
    selected_genes
        Gene symbols to extract.

    Returns
    -------
    dict[str, pandas.DataFrame]
        Named output tables.
    """
    outputs: Dict[str, pd.DataFrame] = {}
    outputs["full"] = dataframe.copy()

    tier1_mask = contains_any_series(
        dataframe=dataframe,
        columns=("priority_tier", "integrated_evidence_class"),
        patterns=("tier_1_protein_coding_rescue", "protein_coding_isoform_rescue"),
    ) & dataframe["isoform_is_rescue_candidate"] & dataframe["isoform_is_protein_coding"]

    outputs["tier1"] = dataframe.loc[tier1_mask].copy()

    top_mask = (
        dataframe["isoform_is_rescue_candidate"]
        & dataframe["isoform_is_protein_coding"]
        & dataframe["gene_level_sperm_or_protein_supported"]
        & dataframe["gene_level_druggable_or_accessible"]
    )
    outputs["top_sperm_supported_druggable"] = dataframe.loc[top_mask].copy()

    selected_symbols = {
        normalise_gene_symbol(value=gene_symbol) for gene_symbol in selected_genes
    }
    outputs["selected_project_genes"] = dataframe.loc[
        dataframe["join_gene_symbol"].isin(selected_symbols)
    ].copy()

    summary = (
        dataframe.groupby("integrated_evidence_class", dropna=False)
        .agg(
            n_rows=("integrated_evidence_class", "size"),
            n_unique_genes=("join_gene_symbol", "nunique"),
            n_isoform_rescue=("isoform_is_rescue_candidate", "sum"),
            n_protein_coding=("isoform_is_protein_coding", "sum"),
            n_sperm_or_protein_supported=("gene_level_sperm_or_protein_supported", "sum"),
            n_druggable_or_accessible=("gene_level_druggable_or_accessible", "sum"),
            median_integrated_review_score=("integrated_review_score", "median"),
            max_integrated_review_score=("integrated_review_score", "max"),
        )
        .reset_index()
        .sort_values(
            by=["max_integrated_review_score", "n_rows"],
            ascending=[False, False],
        )
    )
    outputs["summary"] = summary

    present_selected = set(outputs["selected_project_genes"]["join_gene_symbol"])
    missing = sorted(selected_symbols - present_selected)
    outputs["selected_missing"] = pd.DataFrame({"gene_symbol": missing})
    return outputs


def sort_review_table(*, dataframe: pd.DataFrame) -> pd.DataFrame:
    """
    Sort an integrated review table.

    Parameters
    ----------
    dataframe
        Candidate table.

    Returns
    -------
    pandas.DataFrame
        Sorted table.
    """
    sort_cols = [
        column
        for column in (
            "integrated_review_score",
            "gene_evidence_score",
            "priority_score",
            "target_median_tpm",
            "target_median_isoform_usage",
        )
        if column in dataframe.columns
    ]
    if not sort_cols:
        return dataframe
    return dataframe.sort_values(by=sort_cols, ascending=[False] * len(sort_cols))


def write_tsv(*, dataframe: pd.DataFrame, path: Path, logger: logging.Logger) -> None:
    """
    Write a DataFrame as tab-separated text.

    Parameters
    ----------
    dataframe
        DataFrame to write.
    path
        Output path. If the suffix is .gz, gzip compression is used.
    logger
        Logger instance.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    compression = "gzip" if str(path).endswith(".gz") else None
    dataframe.to_csv(
        path_or_buf=path,
        sep="\t",
        index=False,
        compression=compression,
    )
    logger.info(
        "Wrote %d rows x %d columns to %s",
        dataframe.shape[0],
        dataframe.shape[1],
        path,
    )


def write_gene_list(
    *, dataframe: pd.DataFrame, path: Path, top_n: int, logger: logging.Logger
) -> None:
    """
    Write a simple gene-symbol list for downstream plotting/review.

    Parameters
    ----------
    dataframe
        Candidate table.
    path
        Output text path.
    top_n
        Maximum number of genes.
    logger
        Logger instance.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    if dataframe.empty or "gene_symbol" not in dataframe.columns:
        genes = []
    else:
        genes = (
            dataframe["gene_symbol"]
            .dropna()
            .astype(str)
            .drop_duplicates()
            .head(top_n)
            .tolist()
        )
    path.write_text("\n".join(genes) + ("\n" if genes else ""))
    logger.info("Wrote %d genes to %s", len(genes), path)


def make_safe_excel_sheet_name(*, sheet_name: str) -> str:
    """
    Make a safe Excel worksheet name.

    Parameters
    ----------
    sheet_name
        Proposed worksheet name.

    Returns
    -------
    str
        Excel-safe worksheet name.
    """
    safe_name = re.sub(r"[\[\]:*?/\\]", "_", sheet_name).strip()
    if not safe_name:
        safe_name = "results"
    return safe_name[:31]


def estimate_excel_column_widths(
    *,
    dataframe: pd.DataFrame,
    max_scan_rows: int = 2000,
    min_width: int = 8,
    max_text_width: int = 48,
    max_numeric_width: int = 18,
) -> Dict[str, int]:
    """
    Estimate sensible Excel column widths.

    Parameters
    ----------
    dataframe
        DataFrame to inspect.
    max_scan_rows
        Maximum number of rows scanned.
    min_width
        Minimum width.
    max_text_width
        Maximum width for text columns.
    max_numeric_width
        Maximum width for numeric columns.

    Returns
    -------
    dict[str, int]
        Column-width mapping.
    """
    widths = {}
    scan_dataframe = dataframe.head(n=max_scan_rows)
    for column in dataframe.columns:
        header_width = len(str(column)) + 2
        if scan_dataframe.empty:
            value_width = 0
        else:
            value_width = (
                scan_dataframe[column]
                .astype(str)
                .replace(to_replace="nan", value="")
                .str.len()
                .max()
            )
        if pd.isna(value_width):
            value_width = 0
        raw_width = max(header_width, int(value_width) + 2)
        if pd.api.types.is_numeric_dtype(dataframe[column]):
            width = min(max(raw_width, min_width), max_numeric_width)
        else:
            width = min(max(raw_width, min_width), max_text_width)
        widths[str(column)] = width
    return widths


def write_formatted_excel(
    *, dataframe: pd.DataFrame, path: Path, sheet_name: str, logger: logging.Logger
) -> None:
    """
    Write a formatted Excel copy of a review table.

    Parameters
    ----------
    dataframe
        Table to write.
    path
        Output XLSX path.
    sheet_name
        Worksheet name.
    logger
        Logger instance.
    """
    excel_max_rows = 1_048_576
    excel_max_cols = 16_384
    n_rows, n_cols = dataframe.shape

    if n_rows + 1 > excel_max_rows or n_cols > excel_max_cols:
        logger.warning("Skipping Excel output for %s because it exceeds Excel limits", path)
        return

    path.parent.mkdir(parents=True, exist_ok=True)
    safe_sheet_name = make_safe_excel_sheet_name(sheet_name=sheet_name)
    try:
        with pd.ExcelWriter(path=path, engine="xlsxwriter") as writer:
            dataframe.to_excel(
                excel_writer=writer,
                sheet_name=safe_sheet_name,
                index=False,
            )
            workbook = writer.book
            worksheet = writer.sheets[safe_sheet_name]
            header_format = workbook.add_format(
                {
                    "bold": True,
                    "text_wrap": True,
                    "valign": "top",
                    "border": 1,
                }
            )
            integer_format = workbook.add_format({"num_format": "0"})
            float_format = workbook.add_format({"num_format": "0.0000"})
            worksheet.freeze_panes(1, 0)
            if n_cols > 0 and n_rows > 0:
                worksheet.add_table(
                    first_row=0,
                    first_col=0,
                    last_row=n_rows,
                    last_col=n_cols - 1,
                    options={
                        "columns": [{"header": str(column)} for column in dataframe.columns],
                        "style": "Table Style Medium 2",
                        "autofilter": True,
                    },
                )
            worksheet.set_row(0, 30, header_format)
            column_widths = estimate_excel_column_widths(dataframe=dataframe)
            for column_index, column in enumerate(dataframe.columns):
                if pd.api.types.is_integer_dtype(dataframe[column]):
                    cell_format = integer_format
                elif pd.api.types.is_float_dtype(dataframe[column]):
                    cell_format = float_format
                else:
                    cell_format = None
                worksheet.set_column(
                    first_col=column_index,
                    last_col=column_index,
                    width=column_widths[str(column)],
                    cell_format=cell_format,
                )
    except ImportError as error:
        logger.warning("Could not write Excel output %s: %s", path, error)
        return
    logger.info(
        "Wrote formatted Excel file with %d rows x %d columns to %s",
        n_rows,
        n_cols,
        path,
    )


def write_excel_output_set(
    *,
    outputs: Dict[str, pd.DataFrame],
    output_paths: Dict[str, Path],
    diagnostics_frame: pd.DataFrame,
    prefix: Path,
    logger: logging.Logger,
) -> None:
    """
    Write formatted Excel copies for every tabular output.

    Individual Excel files are written beside the TSV outputs, using the same
    basename and an ``.xlsx`` suffix. A combined workbook containing the main
    review sheets is also written for easier manual inspection.

    Parameters
    ----------
    outputs
        Named output DataFrames.
    output_paths
        TSV output paths for the named output DataFrames.
    diagnostics_frame
        Merge diagnostics table.
    prefix
        Output path prefix.
    logger
        Logger instance.
    """
    logger.info("Excel output is enabled; writing formatted XLSX files")

    for name, dataframe in outputs.items():
        if name not in output_paths:
            logger.warning(
                "Skipping Excel output for %s because no TSV path is defined",
                name,
            )
            continue
        xlsx_path = output_paths[name].with_suffix(".xlsx")
        write_formatted_excel(
            dataframe=dataframe,
            path=xlsx_path,
            sheet_name=name,
            logger=logger,
        )

    diagnostics_path = Path(f"{prefix}.merge_diagnostics.xlsx")
    write_formatted_excel(
        dataframe=diagnostics_frame,
        path=diagnostics_path,
        sheet_name="merge_diagnostics",
        logger=logger,
    )

    combined_path = Path(f"{prefix}.review_workbook.xlsx")
    review_sheets = {
        "top_review": outputs.get("top_sperm_supported_druggable"),
        "tier1_rescue": outputs.get("tier1"),
        "selected_genes": outputs.get("selected_project_genes"),
        "selected_missing": outputs.get("selected_missing"),
        "summary": outputs.get("summary"),
        "merge_diagnostics": diagnostics_frame,
    }
    write_combined_formatted_excel(
        sheets=review_sheets,
        path=combined_path,
        logger=logger,
    )


def write_combined_formatted_excel(
    *,
    sheets: Dict[str, Optional[pd.DataFrame]],
    path: Path,
    logger: logging.Logger,
) -> None:
    """
    Write several review tables to one formatted Excel workbook.

    Parameters
    ----------
    sheets
        Mapping from sheet name to DataFrame. Missing DataFrames are skipped.
    path
        Output XLSX path.
    logger
        Logger instance.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        with pd.ExcelWriter(path=path, engine="xlsxwriter") as writer:
            workbook = writer.book
            header_format = workbook.add_format(
                {
                    "bold": True,
                    "text_wrap": True,
                    "valign": "top",
                    "border": 1,
                }
            )
            integer_format = workbook.add_format({"num_format": "0"})
            float_format = workbook.add_format({"num_format": "0.0000"})

            written_sheets = 0
            used_sheet_names = set()
            for requested_sheet_name, dataframe in sheets.items():
                if dataframe is None:
                    continue

                sheet_name = make_safe_excel_sheet_name(
                    sheet_name=requested_sheet_name
                )
                original_sheet_name = sheet_name
                suffix = 1
                while sheet_name in used_sheet_names:
                    suffix_text = f"_{suffix}"
                    sheet_name = f"{original_sheet_name[:31 - len(suffix_text)]}{suffix_text}"
                    suffix += 1
                used_sheet_names.add(sheet_name)

                dataframe.to_excel(
                    excel_writer=writer,
                    sheet_name=sheet_name,
                    index=False,
                )
                worksheet = writer.sheets[sheet_name]
                n_rows, n_cols = dataframe.shape
                worksheet.freeze_panes(1, 0)

                if n_cols > 0 and n_rows > 0:
                    worksheet.add_table(
                        first_row=0,
                        first_col=0,
                        last_row=n_rows,
                        last_col=n_cols - 1,
                        options={
                            "columns": [
                                {"header": str(column)}
                                for column in dataframe.columns
                            ],
                            "style": "Table Style Medium 2",
                            "autofilter": True,
                        },
                    )
                worksheet.set_row(0, 30, header_format)
                column_widths = estimate_excel_column_widths(dataframe=dataframe)
                for column_index, column in enumerate(dataframe.columns):
                    if pd.api.types.is_integer_dtype(dataframe[column]):
                        cell_format = integer_format
                    elif pd.api.types.is_float_dtype(dataframe[column]):
                        cell_format = float_format
                    else:
                        cell_format = None
                    worksheet.set_column(
                        first_col=column_index,
                        last_col=column_index,
                        width=column_widths[str(column)],
                        cell_format=cell_format,
                    )
                written_sheets += 1
    except ImportError as error:
        logger.warning("Could not write combined Excel output %s: %s", path, error)
        return

    logger.info(
        "Wrote combined formatted Excel workbook with %d sheets to %s",
        written_sheets,
        path,
    )


def prepare_biochem_table(*, dataframe: pd.DataFrame, logger: logging.Logger) -> pd.DataFrame:
    """
    Prepare the biochemical/evidence table for merging.

    Parameters
    ----------
    dataframe
        Raw biochemical table.
    logger
        Logger instance.

    Returns
    -------
    pandas.DataFrame
        Prefixed evidence table.
    """
    table = add_join_keys(
        dataframe=dataframe,
        symbol_candidates=("gene_key", "hgnc_symbol", "Description", "gene_name"),
        gene_id_candidates=("ensembl_gene_id", "Name", "gene_id"),
        logger=logger,
    )
    keep_columns = [
        "gene_key",
        "in_testis_high_conf_final",
        "in_hpo_gene_set",
        "clinvar_best_present",
        "clinvar_best_pathogenic_present",
        "clinvar_hc_present",
        "clinvar_hc_pathogenic_present",
        "in_literature_fertility_set",
        "sperm_present_any",
        "sperm_present_frac",
        "sperm_tpm_mean",
        "sperm_tpm_median",
        "prot_present_any",
        "prot_present_fraction",
        "proteomics_evidence_level",
        "public_prot_present_any",
        "public_prot_present_fraction",
        "public_proteomics_evidence_level",
        "prot_unique_peptides_max",
        "prot_coverage_pct_max",
        "lit_support_protein",
        "lit_support_sperm_rna",
        "lit_support_testis_rna",
        "lit_mode_of_action",
        "lit_reference_role",
        "target_median_tpm",
        "max_non_target_tissue",
        "max_non_target_median_tpm",
        "log2_fc_target_vs_max_non_target",
        "tau",
        "target_present_fraction",
        "is_target_max_tissue",
        "gene_full_name",
        "gene_summary",
        "gene_type",
        "go_bp",
        "go_mf",
        "go_cc",
        "hpo_reproductive_term_count",
        "hpo_reproductive_ids",
        "hpo_reproductive_terms",
        "clinvar_hc_pathogenic_gene_summary__n_variants",
        "clinvar_hc_pathogenic_gene_summary__top_phenotypes",
        "clinvar_hc_pathogenic_gene_summary__review_statuses",
        "uniprot_acc",
        "uniprot_subcellular_location",
        "uniprot_keywords",
        "uniprot_go_cc",
        "uniprot_protein_name",
        "has_signal_peptide",
        "has_transmembrane",
        "is_extracellular",
        "is_plasma_membrane",
        "is_membrane",
        "is_secreted",
        "is_cell_surface_candidate",
        "predicted_target_class",
        "biochemical_accessibility_score",
    ]
    return make_prefixed_evidence_table(
        dataframe=table,
        prefix="biochem_",
        keep_columns=keep_columns,
    )


def prepare_gene_priority_table(
    *, dataframe: pd.DataFrame, logger: logging.Logger
) -> pd.DataFrame:
    """
    Prepare the gene-level priority/Open Targets table for merging.

    Parameters
    ----------
    dataframe
        Raw gene priority table.
    logger
        Logger instance.

    Returns
    -------
    pandas.DataFrame
        Prefixed evidence table.
    """
    table = add_join_keys(
        dataframe=dataframe,
        symbol_candidates=("gene_key", "ot_approved_symbol", "hgnc_symbol", "Description"),
        gene_id_candidates=("gene_id", "ensembl_gene_id", "Name"),
        logger=logger,
    )
    keep_columns = [
        "gene_key",
        "gene_id",
        "priority_score",
        "score_components",
        "candidate_druggable_sperm_protein",
        "n_memberships",
        "membership_names",
        "prot_class",
        "prot_any_detected_or_strong",
        "ot_any_tractable",
        "ot_any_small_molecule_tractable",
        "ot_any_antibody_tractable",
        "ot_any_protac_tractable",
        "ot_approved_symbol",
        "ot_tractability_summary",
        "biochemical_accessibility_score",
        "biochemical_accessibility_score_norm",
        "in_testis_high_conf_final",
        "in_hpo_gene_set",
        "in_literature_fertility_set",
        "sperm_present_any",
        "sperm_present_frac",
        "sperm_tpm_mean",
        "sperm_tpm_median",
        "prot_present_any",
        "proteomics_evidence_level",
        "public_prot_present_any",
        "public_proteomics_evidence_level",
        "prot_unique_peptides_max",
        "prot_coverage_pct_max",
        "gene_full_name",
        "gene_summary",
        "gene_type",
        "go_bp",
        "go_mf",
        "go_cc",
        "uniprot_acc",
        "uniprot_subcellular_location",
        "uniprot_keywords",
        "uniprot_go_cc",
        "uniprot_protein_name",
        "has_signal_peptide",
        "has_transmembrane",
        "is_extracellular",
        "is_plasma_membrane",
        "is_membrane",
        "is_secreted",
        "is_cell_surface_candidate",
        "predicted_target_class",
        "proteomics_present_internal",
        "proteomics_present_public",
        "proteomics_present_any_source",
        "sperm_rnaseq_present",
    ]
    return make_prefixed_evidence_table(
        dataframe=table,
        prefix="geneprio_",
        keep_columns=keep_columns,
    )


def prepare_tractability_table(
    *, dataframe: pd.DataFrame, logger: logging.Logger
) -> pd.DataFrame:
    """
    Prepare the full gene-universe tractability table for merging.

    Parameters
    ----------
    dataframe
        Raw tractability table.
    logger
        Logger instance.

    Returns
    -------
    pandas.DataFrame
        Prefixed evidence table.
    """
    table = add_join_keys(
        dataframe=dataframe,
        symbol_candidates=("gene_key", "gene_name", "ot_approved_symbol"),
        gene_id_candidates=("gene_id",),
        logger=logger,
    )
    keep_columns = [
        "gene_key",
        "gene_id",
        "gene_name",
        "Chromosome",
        "Start",
        "End",
        "Strand",
        "gene_length_bp",
        "ot_any_small_molecule_tractable",
        "ot_any_antibody_tractable",
        "ot_any_protac_tractable",
        "ot_tractability_summary",
        "ot_approved_symbol",
    ]
    return make_prefixed_evidence_table(
        dataframe=table,
        prefix="tract_",
        keep_columns=keep_columns,
    )


def move_review_columns_to_front(*, dataframe: pd.DataFrame) -> pd.DataFrame:
    """
    Move key review columns to the front of the table.

    Parameters
    ----------
    dataframe
        Integrated table.

    Returns
    -------
    pandas.DataFrame
        Reordered table.
    """
    front_cols = [
        "gene_symbol",
        "transcript_id_with_version",
        "transcript_id",
        "gene_id",
        "priority_tier",
        "integrated_evidence_class",
        "integrated_review_score",
        "gene_evidence_score",
        "priority_score",
        "target_median_tpm",
        "max_non_target_median_tpm",
        "target_median_isoform_usage",
        "max_non_target_isoform_usage",
        "log2_target_vs_max_non_target_isoform_usage",
        "target_gene_median_tpm",
        "max_non_target_gene_median_tpm",
        "log2_target_vs_max_non_target_gene_tpm",
        "candidate_rank_tier",
        "isoform_is_rescue_candidate",
        "isoform_is_protein_coding",
        "gene_level_has_sperm_rna",
        "gene_level_has_any_proteomics",
        "gene_level_has_strong_public_proteomics",
        "gene_level_has_any_ot_tractability",
        "gene_level_ot_small_molecule",
        "gene_level_ot_antibody",
        "gene_level_ot_protac",
        "gene_level_has_accessibility_signal",
        "gene_level_is_cell_surface",
        "gene_level_is_secreted",
        "gene_level_is_membrane",
        "gene_level_candidate_druggable_sperm_protein",
        "gene_level_has_known_or_clinical_fertility_support",
        "integrated_score_components",
        "gencode_transcript_name",
        "gencode_transcript_type",
        "gencode_cds_length_bp",
        "gencode_exon_count",
        "gencode_cds_exon_count",
        "gencode_protein_id",
        "biochem_predicted_target_class",
        "biochem_biochemical_accessibility_score",
        "geneprio_prot_class",
        "geneprio_candidate_druggable_sperm_protein",
        "geneprio_ot_tractability_summary",
        "tract_ot_tractability_summary",
        "biochem_uniprot_subcellular_location",
        "biochem_uniprot_keywords",
        "biochem_gene_summary",
    ]
    available_front = [column for column in front_cols if column in dataframe.columns]
    remaining = [column for column in dataframe.columns if column not in available_front]
    return dataframe[available_front + remaining]


def run(*, config: Config, logger: logging.Logger) -> None:
    """
    Run the full evidence integration workflow.

    Parameters
    ----------
    config
        Workflow configuration.
    logger
        Logger instance.
    """
    config.out_dir.mkdir(parents=True, exist_ok=True)
    diagnostics = []

    with timed(logger=logger, label="Load isoform candidates"):
        isoforms = read_table(path=config.isoform_candidates_tsv, logger=logger)
        isoforms = add_join_keys(
            dataframe=isoforms,
            symbol_candidates=("gene_symbol", "gencode_gene_name", "gene_key"),
            gene_id_candidates=("gene_id", "gene_id_with_version", "gencode_gene_id"),
            logger=logger,
        )

    with timed(logger=logger, label="Load biochemical gene evidence"):
        biochem_raw = read_excel_sheet(
            path=config.biochem_xlsx,
            sheet_name=config.biochem_sheet,
            logger=logger,
        )
        biochem = prepare_biochem_table(dataframe=biochem_raw, logger=logger)
        isoforms, diag = merge_by_gene_keys(
            base=isoforms,
            evidence=biochem,
            source_name="biochemical_gene_evidence",
            logger=logger,
        )
        diagnostics.append(diag)

    if config.gene_priority_xlsx is not None:
        with timed(logger=logger, label="Load gene priority/Open Targets evidence"):
            gene_priority_raw = read_excel_sheet(
                path=config.gene_priority_xlsx,
                sheet_name=config.gene_priority_sheet,
                logger=logger,
            )
            gene_priority = prepare_gene_priority_table(
                dataframe=gene_priority_raw,
                logger=logger,
            )
            isoforms, diag = merge_by_gene_keys(
                base=isoforms,
                evidence=gene_priority,
                source_name="gene_priority_open_targets_evidence",
                logger=logger,
            )
            diagnostics.append(diag)

    if config.tractability_tsv is not None:
        with timed(logger=logger, label="Load full-universe tractability evidence"):
            tractability_raw = read_table(path=config.tractability_tsv, logger=logger)
            tractability = prepare_tractability_table(
                dataframe=tractability_raw,
                logger=logger,
            )
            isoforms, diag = merge_by_gene_keys(
                base=isoforms,
                evidence=tractability,
                source_name="full_universe_tractability_evidence",
                logger=logger,
            )
            diagnostics.append(diag)

    with timed(logger=logger, label="Derive integrated evidence fields"):
        integrated = derive_integrated_evidence(dataframe=isoforms)
        integrated = sort_review_table(dataframe=integrated)
        integrated = move_review_columns_to_front(dataframe=integrated)
        outputs = select_outputs(
            dataframe=integrated,
            selected_genes=config.selected_genes,
        )
        for key in list(outputs):
            outputs[key] = sort_review_table(dataframe=outputs[key])
            outputs[key] = move_review_columns_to_front(dataframe=outputs[key])

    with timed(logger=logger, label="Write outputs"):
        prefix = config.out_dir / config.out_prefix
        output_paths = {
            "full": Path(f"{prefix}.isoform_candidates_with_gene_level_evidence.tsv"),
            "tier1": Path(
                f"{prefix}.tier1_protein_coding_rescue_with_gene_evidence.tsv"
            ),
            "top_sperm_supported_druggable": Path(
                f"{prefix}.top_sperm_supported_druggable_isoform_rescue_candidates.tsv"
            ),
            "selected_project_genes": Path(
                f"{prefix}.selected_project_genes_isoform_gene_evidence.tsv"
            ),
            "selected_missing": Path(
                f"{prefix}.selected_project_genes_missing.tsv"
            ),
            "summary": Path(f"{prefix}.summary_by_integrated_evidence_class.tsv"),
        }
        for name, path in output_paths.items():
            write_tsv(dataframe=outputs[name], path=path, logger=logger)

        diagnostics_frame = pd.DataFrame(diagnostics)
        write_tsv(
            dataframe=diagnostics_frame,
            path=Path(f"{prefix}.merge_diagnostics.tsv"),
            logger=logger,
        )

        write_gene_list(
            dataframe=outputs["top_sperm_supported_druggable"],
            path=Path(f"{prefix}.genes_to_review_top_{config.top_n_review}.txt"),
            top_n=config.top_n_review,
            logger=logger,
        )
        write_gene_list(
            dataframe=outputs["tier1"],
            path=Path(f"{prefix}.tier1_genes_to_review_top_{config.top_n_review}.txt"),
            top_n=config.top_n_review,
            logger=logger,
        )

        if config.write_excel_outputs:
            write_excel_output_set(
                outputs=outputs,
                output_paths=output_paths,
                diagnostics_frame=diagnostics_frame,
                prefix=prefix,
                logger=logger,
            )
        else:
            logger.info("Excel output disabled by --no_write_excel_outputs")

    logger.info(
        "Integration finished: %d total rows, %d tier1 rescue rows, %d top review rows",
        outputs["full"].shape[0],
        outputs["tier1"].shape[0],
        outputs["top_sperm_supported_druggable"].shape[0],
    )


def main() -> None:
    """Run the command-line entry point."""
    config = parse_args()
    logger = setup_logging(log_level=config.log_level, log_path=config.log_path)
    logger.info("GTEx isoform/gene-level evidence integration starting")
    logger.info("Isoform candidate table: %s", config.isoform_candidates_tsv)
    logger.info("Biochemical evidence workbook: %s", config.biochem_xlsx)
    logger.info("Gene priority workbook: %s", config.gene_priority_xlsx)
    logger.info("Tractability table: %s", config.tractability_tsv)
    logger.info("Output directory: %s", config.out_dir)
    run(config=config, logger=logger)


if __name__ == "__main__":
    main()
