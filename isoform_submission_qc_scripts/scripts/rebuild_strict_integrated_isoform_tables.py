#!/usr/bin/env python3
"""
Rebuild strict integrated isoform candidate tables with stable definitions.

This script repairs a common downstream-analysis problem: an integrated table may
contain rows classified as protein-coding rescue candidates by later harmonised
flags even though they were not part of the strict transcript-level tier-1 set.
For manuscript-ready counts, this script keeps definitions explicit.

The default strict tier-1 definition requires:
- priority_tier == tier_1_protein_coding_rescue_candidate
- strict GENCODE protein-coding flag where available
- CDS support where available
- broad-gene isoform-rescue status where available

It then writes strict tier-1 rows, strict sperm/protein-supported and druggable
rows, selected-gene rows, and summaries.
"""

from __future__ import annotations

import argparse
import logging
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Optional, Sequence

import pandas as pd


LOGGER_NAME = "rebuild_strict_integrated_isoform_tables"
STRICT_TIER1 = "tier_1_protein_coding_rescue_candidate"
SUPPORTED_DRUGGABLE_CLASS = (
    "protein_coding_isoform_rescue_sperm_supported_druggable"
)


@dataclass(frozen=True)
class Config:
    """Runtime configuration."""

    integrated_table: Path
    out_dir: Path
    out_prefix: str
    selected_genes: Sequence[str]
    require_priority_tier1: bool
    require_strict_protein_coding: bool
    require_cds: bool
    require_rescue: bool
    write_excel_outputs: bool
    log_level: str
    log_path: Optional[Path]


def setup_logging(*, log_level: str, log_path: Optional[Path]) -> logging.Logger:
    """Configure logging and return a logger."""
    logger = logging.getLogger(LOGGER_NAME)
    logger.handlers.clear()
    logger.setLevel(getattr(logging, log_level.upper(), logging.INFO))
    formatter = logging.Formatter(
        "%(asctime)s | %(levelname)s | %(message)s", "%Y-%m-%d %H:%M:%S"
    )
    stream_handler = logging.StreamHandler(sys.stderr)
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)
    if log_path is not None:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(log_path)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
    return logger


def parse_args(argv: Optional[Sequence[str]] = None) -> Config:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Rebuild strict integrated isoform candidate outputs."
    )
    parser.add_argument("--integrated_table", required=True, type=Path)
    parser.add_argument("--out_dir", required=True, type=Path)
    parser.add_argument(
        "--out_prefix",
        default="gtex_v11_isoform_candidates_gene_level_evidence.strict_rebuilt",
    )
    parser.add_argument(
        "--selected_genes",
        nargs="*",
        default=["AFG2B", "CFAP99", "SLC16A7", "ABCG4"],
        help="Project genes to extract from the rebuilt table.",
    )
    parser.add_argument(
        "--allow_non_tier1",
        action="store_true",
        help="Do not require the original strict tier-1 priority label.",
    )
    parser.add_argument(
        "--allow_non_strict_protein_coding",
        action="store_true",
        help="Do not require strict GENCODE protein-coding flag.",
    )
    parser.add_argument(
        "--allow_no_cds",
        action="store_true",
        help="Do not require a CDS flag/length where available.",
    )
    parser.add_argument(
        "--allow_non_rescue",
        action="store_true",
        help="Do not require broad-gene isoform-rescue status.",
    )
    parser.add_argument("--no_write_excel_outputs", action="store_true")
    parser.add_argument("--log_level", default="INFO")
    parser.add_argument("--log_path", default=None, type=Path)
    args = parser.parse_args(argv)
    return Config(
        integrated_table=args.integrated_table,
        out_dir=args.out_dir,
        out_prefix=args.out_prefix,
        selected_genes=args.selected_genes,
        require_priority_tier1=not args.allow_non_tier1,
        require_strict_protein_coding=not args.allow_non_strict_protein_coding,
        require_cds=not args.allow_no_cds,
        require_rescue=not args.allow_non_rescue,
        write_excel_outputs=not args.no_write_excel_outputs,
        log_level=args.log_level,
        log_path=args.log_path,
    )


def read_table(*, path: Path, logger: logging.Logger) -> pd.DataFrame:
    """Read a TSV/TSV.GZ/XLSX table."""
    logger.info("Reading %s", path)
    if "".join(path.suffixes).lower().endswith(".xlsx"):
        dataframe = pd.read_excel(path)
    else:
        dataframe = pd.read_csv(path, sep="\t", low_memory=False)
    logger.info("Loaded %d rows x %d columns", *dataframe.shape)
    return dataframe


def as_bool(series: pd.Series) -> pd.Series:
    """Convert common bool-like values to Boolean values."""
    return series.astype(str).str.strip().str.lower().isin({"true", "1", "yes", "y"})


def numeric_positive(series: pd.Series) -> pd.Series:
    """Return True for numeric values greater than zero."""
    return pd.to_numeric(series, errors="coerce").fillna(0) > 0


def first_existing_bool(*, dataframe: pd.DataFrame, columns: Iterable[str]) -> pd.Series:
    """Return Boolean OR over available columns, else all-False."""
    output = pd.Series(False, index=dataframe.index)
    for column in columns:
        if column in dataframe.columns:
            output = output | as_bool(dataframe[column])
    return output


def strict_protein_coding_mask(*, dataframe: pd.DataFrame) -> pd.Series:
    """
    Build a strict protein-coding mask from GENCODE-derived fields.

    Returns true for rows with a strict GENCODE protein-coding transcript flag or
    a transcript type exactly equal to ``protein_coding``. It deliberately does
    not count ``protein_coding_CDS_not_defined`` or
    ``nonsense_mediated_decay`` as strict protein-coding targets.
    """
    mask = first_existing_bool(
        dataframe=dataframe,
        columns=("gencode_is_protein_coding_transcript", "flag_protein_coding_transcript"),
    )
    if "gencode_transcript_type" in dataframe.columns:
        mask = mask | (dataframe["gencode_transcript_type"].astype(str) == "protein_coding")
    return mask


def cds_mask(*, dataframe: pd.DataFrame) -> pd.Series:
    """Build a CDS-support mask from available fields."""
    mask = first_existing_bool(dataframe=dataframe, columns=("gencode_has_cds", "flag_has_cds"))
    if "gencode_cds_length_bp" in dataframe.columns:
        mask = mask | numeric_positive(dataframe["gencode_cds_length_bp"])
    return mask


def rescue_mask(*, dataframe: pd.DataFrame) -> pd.Series:
    """Build an isoform-rescue mask from available fields."""
    return first_existing_bool(
        dataframe=dataframe,
        columns=(
            "is_broad_gene_isoform_rescue_candidate",
            "flag_broad_gene_isoform_rescue",
            "isoform_is_rescue_candidate",
        ),
    )


def supported_mask(*, dataframe: pd.DataFrame) -> pd.Series:
    """Return rows with sperm RNA or proteomics support."""
    return first_existing_bool(
        dataframe=dataframe,
        columns=(
            "gene_level_has_sperm_rna",
            "gene_level_has_any_proteomics",
            "gene_level_has_strong_public_proteomics",
        ),
    )


def druggable_or_accessible_mask(*, dataframe: pd.DataFrame) -> pd.Series:
    """Return rows with tractability or accessibility evidence."""
    return first_existing_bool(
        dataframe=dataframe,
        columns=(
            "gene_level_has_any_ot_tractability",
            "gene_level_ot_small_molecule",
            "gene_level_ot_antibody",
            "gene_level_ot_protac",
            "gene_level_has_accessibility_signal",
            "gene_level_candidate_druggable_sperm_protein",
        ),
    )


def build_strict_mask(*, dataframe: pd.DataFrame, config: Config) -> pd.Series:
    """Build the final strict table mask."""
    mask = pd.Series(True, index=dataframe.index)
    if config.require_priority_tier1 and "priority_tier" in dataframe.columns:
        mask = mask & (dataframe["priority_tier"] == STRICT_TIER1)
    if config.require_strict_protein_coding:
        mask = mask & strict_protein_coding_mask(dataframe=dataframe)
    if config.require_cds:
        mask = mask & cds_mask(dataframe=dataframe)
    if config.require_rescue:
        mask = mask & rescue_mask(dataframe=dataframe)
    return mask


def sort_for_review(*, dataframe: pd.DataFrame) -> pd.DataFrame:
    """Sort output by integrated review score and supporting evidence."""
    output = dataframe.copy()
    for column in ("integrated_review_score", "priority_score", "target_median_tpm"):
        if column in output.columns:
            output[column] = pd.to_numeric(output[column], errors="coerce")
    sort_cols = [
        column
        for column in (
            "integrated_review_score",
            "priority_score",
            "target_median_tpm",
        )
        if column in output.columns
    ]
    if sort_cols:
        output = output.sort_values(sort_cols, ascending=[False] * len(sort_cols))
    return output


def summarise_outputs(*, all_rows: pd.DataFrame, strict_rows: pd.DataFrame) -> pd.DataFrame:
    """Create summary counts for rebuilt tables."""
    rows = []
    for label, dataframe in (("input_all_rows", all_rows), ("strict_rebuilt_rows", strict_rows)):
        row = {
            "table": label,
            "n_rows": len(dataframe),
            "n_unique_genes": dataframe.get("gene_symbol", pd.Series(dtype=object)).nunique(),
            "n_unique_transcripts": dataframe.get(
                "transcript_id_with_version", pd.Series(dtype=object)
            ).nunique(),
            "n_supported": int(supported_mask(dataframe=dataframe).sum()),
            "n_druggable_or_accessible": int(druggable_or_accessible_mask(dataframe=dataframe).sum()),
            "n_supported_and_druggable_or_accessible": int(
                (supported_mask(dataframe=dataframe) & druggable_or_accessible_mask(dataframe=dataframe)).sum()
            ),
        }
        rows.append(row)
    return pd.DataFrame(rows)


def write_tsv(*, dataframe: pd.DataFrame, path: Path, logger: logging.Logger) -> None:
    """Write a TSV output."""
    path.parent.mkdir(parents=True, exist_ok=True)
    dataframe.to_csv(path, sep="\t", index=False)
    logger.info("Wrote %d rows x %d columns to %s", *dataframe.shape, path)


def write_excel(*, dataframe: pd.DataFrame, path: Path, sheet_name: str, logger: logging.Logger) -> None:
    """Write a formatted Excel output."""
    if len(dataframe) + 1 > 1_048_576:
        logger.warning("Skipping Excel output for %s because it exceeds Excel limits", path)
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(path, engine="xlsxwriter") as writer:
        dataframe.to_excel(writer, sheet_name=sheet_name[:31], index=False)
        workbook = writer.book
        worksheet = writer.sheets[sheet_name[:31]]
        header_format = workbook.add_format({"bold": True, "text_wrap": True, "border": 1})
        float_format = workbook.add_format({"num_format": "0.0000"})
        worksheet.freeze_panes(1, 0)
        worksheet.set_row(0, 30, header_format)
        if dataframe.shape[1] > 0:
            worksheet.add_table(
                0,
                0,
                len(dataframe),
                dataframe.shape[1] - 1,
                {
                    "columns": [{"header": str(c)} for c in dataframe.columns],
                    "autofilter": True,
                    "style": "Table Style Medium 2",
                },
            )
        for i, column in enumerate(dataframe.columns):
            sample = dataframe[column].head(1000).astype(str)
            value_width = int(sample.str.len().max()) if not sample.empty else 0
            width = min(max(len(str(column)) + 2, value_width + 2, 8), 45)
            fmt = float_format if pd.api.types.is_float_dtype(dataframe[column]) else None
            worksheet.set_column(i, i, width, fmt)
    logger.info("Wrote formatted Excel output to %s", path)


def write_outputs(*, outputs: Dict[str, pd.DataFrame], config: Config, logger: logging.Logger) -> None:
    """Write all outputs as TSV and optionally XLSX."""
    for suffix, dataframe in outputs.items():
        tsv_path = config.out_dir / f"{config.out_prefix}.{suffix}.tsv"
        write_tsv(dataframe=dataframe, path=tsv_path, logger=logger)
        if config.write_excel_outputs:
            xlsx_path = config.out_dir / f"{config.out_prefix}.{suffix}.xlsx"
            write_excel(dataframe=dataframe, path=xlsx_path, sheet_name=suffix, logger=logger)


def run(*, config: Config, logger: logging.Logger) -> Dict[str, pd.DataFrame]:
    """Run strict output rebuilding."""
    integrated = read_table(path=config.integrated_table, logger=logger)
    strict_rows = sort_for_review(dataframe=integrated.loc[build_strict_mask(dataframe=integrated, config=config)].copy())
    top_supported_druggable = strict_rows.loc[
        supported_mask(dataframe=strict_rows) & druggable_or_accessible_mask(dataframe=strict_rows)
    ].copy()
    selected = strict_rows.loc[
        strict_rows.get("gene_symbol", pd.Series(dtype=str)).astype(str).isin(config.selected_genes)
    ].copy()
    missing = pd.DataFrame(
        {
            "gene_symbol": [
                gene
                for gene in config.selected_genes
                if gene not in set(selected.get("gene_symbol", pd.Series(dtype=str)).astype(str))
            ]
        }
    )
    summary = summarise_outputs(all_rows=integrated, strict_rows=strict_rows)
    top_genes = pd.DataFrame(
        {
            "gene_symbol": strict_rows.get("gene_symbol", pd.Series(dtype=str))
            .dropna()
            .astype(str)
            .drop_duplicates()
            .head(50)
        }
    )
    outputs = {
        "strict_tier1_protein_coding_rescue_candidates": strict_rows,
        "strict_top_supported_druggable_isoform_rescue_candidates": top_supported_druggable,
        "strict_selected_project_genes": selected,
        "strict_selected_project_genes_missing": missing,
        "strict_summary": summary,
        "strict_top_50_genes_to_review": top_genes,
    }
    write_outputs(outputs=outputs, config=config, logger=logger)
    return outputs


def main(argv: Optional[Sequence[str]] = None) -> None:
    """Run the command-line program."""
    config = parse_args(argv)
    logger = setup_logging(log_level=config.log_level, log_path=config.log_path)
    logger.info("Starting strict integrated table rebuild")
    run(config=config, logger=logger)
    logger.info("Finished strict integrated table rebuild")


if __name__ == "__main__":
    main()
