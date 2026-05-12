#!/usr/bin/env python3
"""
First-pass GTEx transcript-level isoform-usage analysis for sperm target genes.

This script is designed as the transcript-level companion to the existing
GTEx gene-level tissue-specificity scripts. It extracts selected genes from the
large GTEx RSEM transcript TPM matrix, maps samples to GTEx tissues, computes
transcript-level tissue medians, and computes isoform usage within each gene.

The key biological question is:

    Does a target gene that is not strongly testis-specific at gene level have
    an individual transcript isoform that is preferentially expressed or used
    in testis?

Recommended first-pass use is for selected targets such as:

    AFG2B CFAP99 SLC16A7 ABCG4

Primary outputs are tab-separated files:

1. <prefix>.transcript_subset_tpm_matrix.tsv.gz
   Selected transcript TPM matrix, with transcript and gene metadata.

2. <prefix>.transcript_sample_tpm_long.tsv.gz
   Long-format transcript TPM values for selected transcripts and samples.

3. <prefix>.transcript_tissue_median_tpm.tsv
   Transcript by tissue median TPM matrix.

4. <prefix>.transcript_isoform_usage_median_by_tissue.tsv
   Transcript by tissue median isoform-usage matrix. Usage is calculated
   sample-wise as transcript TPM divided by total TPM for that gene in that
   sample, then summarised by tissue.

5. <prefix>.transcript_testis_isoform_summary.tsv
   Ranked summary of target-tissue transcript expression and isoform usage.

6. <prefix>.transcript_subset_expected_count_matrix.tsv.gz
   Optional selected transcript expected-count matrix, if an expected-count
   file is provided. This is extracted for later formal modelling; TPM remains
   the preferred input for the first-pass tissue/usage summaries.

All generated tabular outputs are TSV/TSV.GZ, not comma-separated files.
"""

from __future__ import annotations

import argparse
import gzip
import logging
import sys
import time
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Set, Tuple

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class Config:
    """
    Configuration for GTEx transcript-level isoform analysis.

    Attributes
    ----------
    transcript_tpm_path
        Path to the GTEx RSEM transcript TPM matrix.
    gene_tpm_gct_path
        Path to the GTEx gene TPM GCT file used to map Ensembl gene IDs to
        gene symbols.
    sample_attributes_path
        Path to GTEx SampleAttributesDS.txt.
    out_dir
        Output directory.
    out_prefix
        Prefix for output files.
    expected_count_path
        Optional path to the GTEx RSEM transcript expected-count matrix.
    target_genes
        Gene symbols and/or Ensembl gene IDs requested on the command line.
    genes_tsv
        Optional TSV containing selected genes.
    gene_symbol_col
        Optional gene-symbol column name in genes_tsv. If omitted, the script
        tries common names.
    ensembl_gene_id_col
        Optional Ensembl-gene-ID column name in genes_tsv. If omitted, the
        script tries common names.
    sample_id_col
        Sample ID column in SampleAttributesDS.txt.
    tissue_col
        Tissue column in SampleAttributesDS.txt.
    target_tissue
        Tissue of interest for scoring.
    min_tpm_present
        TPM threshold used to calculate target-tissue presence fraction.
    chunk_size
        Number of transcript rows to read per chunk from the large matrix.
    log_path
        Optional log file path.
    log_level
        Logging level.
    epsilon
        Small constant used to avoid division by zero.
    """

    transcript_tpm_path: Path
    gene_tpm_gct_path: Path
    sample_attributes_path: Path
    out_dir: Path
    out_prefix: str
    expected_count_path: Optional[Path]
    target_genes: Tuple[str, ...]
    genes_tsv: Optional[Path]
    gene_symbol_col: Optional[str]
    ensembl_gene_id_col: Optional[str]
    sample_id_col: str
    tissue_col: str
    target_tissue: str
    min_tpm_present: float
    chunk_size: int
    log_path: Optional[Path]
    log_level: str
    epsilon: float = 1e-8


@dataclass(frozen=True)
class TranscriptSchema:
    """
    Column schema for a GTEx transcript matrix.

    Attributes
    ----------
    transcript_col
        Column containing transcript IDs.
    gene_col
        Column containing Ensembl gene IDs.
    metadata_cols
        Non-sample columns to preserve from the matrix.
    sample_cols
        Sample expression/count columns.
    """

    transcript_col: str
    gene_col: str
    metadata_cols: Tuple[str, ...]
    sample_cols: Tuple[str, ...]


@dataclass(frozen=True)
class ResolvedTargets:
    """
    Resolved target-gene information.

    Attributes
    ----------
    requested_terms
        User-requested gene symbols and/or Ensembl gene IDs.
    requested_symbols
        User-requested symbols, upper-cased.
    requested_gene_ids
        User-requested Ensembl gene IDs with version suffixes removed.
    resolved_gene_ids
        Ensembl gene IDs to extract from the transcript matrix.
    unresolved_terms
        Requested terms that could not be resolved from the available map.
    """

    requested_terms: Tuple[str, ...]
    requested_symbols: Tuple[str, ...]
    requested_gene_ids: Tuple[str, ...]
    resolved_gene_ids: Tuple[str, ...]
    unresolved_terms: Tuple[str, ...]


@contextmanager
def timed(*, logger: logging.Logger, label: str) -> Iterator[None]:
    """
    Log start, completion, and elapsed time for a block.

    Parameters
    ----------
    logger
        Logger instance.
    label
        Human-readable name for the timed block.

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
    Configure logging to stdout and optionally to a file.

    Parameters
    ----------
    log_level
        Logging level name.
    log_path
        Optional path to a log file.

    Returns
    -------
    logging.Logger
        Configured logger.
    """
    logger = logging.getLogger("gtex_transcript_isoform_first_pass")
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


def parse_args() -> Config:
    """
    Parse command-line arguments.

    Returns
    -------
    Config
        Parsed configuration object.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Extract selected GTEx transcript TPMs and compute testis-focused "
            "isoform-usage summaries."
        )
    )
    parser.add_argument(
        "--transcript_tpm_path",
        required=True,
        type=Path,
        help="Path to GTEx RSEM transcript TPM matrix txt.gz.",
    )
    parser.add_argument(
        "--gene_tpm_gct_path",
        required=True,
        type=Path,
        help="Path to GTEx gene TPM GCT/GCT.GZ for Ensembl-to-symbol mapping.",
    )
    parser.add_argument(
        "--sample_attributes_path",
        required=True,
        type=Path,
        help="Path to GTEx SampleAttributesDS.txt.",
    )
    parser.add_argument(
        "--expected_count_path",
        default=None,
        type=Path,
        help="Optional path to GTEx RSEM transcript expected-count matrix txt.gz.",
    )
    parser.add_argument(
        "--target_genes",
        nargs="*",
        default=(),
        help=(
            "Gene symbols and/or Ensembl gene IDs to extract, for example: "
            "AFG2B CFAP99 SLC16A7 ABCG4."
        ),
    )
    parser.add_argument(
        "--genes_tsv",
        default=None,
        type=Path,
        help=(
            "Optional TSV of target genes. The script auto-detects common "
            "gene-symbol and Ensembl-gene-ID columns unless column names are "
            "given explicitly."
        ),
    )
    parser.add_argument(
        "--gene_symbol_col",
        default=None,
        help="Optional gene-symbol column name in --genes_tsv.",
    )
    parser.add_argument(
        "--ensembl_gene_id_col",
        default=None,
        help="Optional Ensembl-gene-ID column name in --genes_tsv.",
    )
    parser.add_argument(
        "--sample_id_col",
        default="SAMPID",
        help="Sample ID column in SampleAttributesDS.txt. Default: SAMPID.",
    )
    parser.add_argument(
        "--tissue_col",
        default="SMTSD",
        help="Tissue column in SampleAttributesDS.txt. Default: SMTSD.",
    )
    parser.add_argument(
        "--target_tissue",
        default="Testis",
        help="Target tissue for scoring. Default: Testis.",
    )
    parser.add_argument(
        "--min_tpm_present",
        default=1.0,
        type=float,
        help="TPM threshold for target-tissue presence fraction. Default: 1.0.",
    )
    parser.add_argument(
        "--chunk_size",
        default=250,
        type=int,
        help="Rows per chunk when scanning the large transcript matrix. Default: 250.",
    )
    parser.add_argument(
        "--out_dir",
        required=True,
        type=Path,
        help="Output directory.",
    )
    parser.add_argument(
        "--out_prefix",
        default="selected_targets_gtex_v11",
        help="Output prefix. Default: selected_targets_gtex_v11.",
    )
    parser.add_argument(
        "--log_path",
        default=None,
        type=Path,
        help="Optional log file path.",
    )
    parser.add_argument(
        "--log_level",
        default="INFO",
        help="Logging level. Default: INFO.",
    )

    args = parser.parse_args()

    if not args.target_genes and args.genes_tsv is None:
        raise ValueError("Provide --target_genes and/or --genes_tsv.")

    return Config(
        transcript_tpm_path=args.transcript_tpm_path,
        gene_tpm_gct_path=args.gene_tpm_gct_path,
        sample_attributes_path=args.sample_attributes_path,
        out_dir=args.out_dir,
        out_prefix=args.out_prefix,
        expected_count_path=args.expected_count_path,
        target_genes=tuple(args.target_genes),
        genes_tsv=args.genes_tsv,
        gene_symbol_col=args.gene_symbol_col,
        ensembl_gene_id_col=args.ensembl_gene_id_col,
        sample_id_col=args.sample_id_col,
        tissue_col=args.tissue_col,
        target_tissue=args.target_tissue,
        min_tpm_present=args.min_tpm_present,
        chunk_size=args.chunk_size,
        log_path=args.log_path,
        log_level=args.log_level,
    )


def open_text_maybe_gzip(*, path: Path):
    """
    Open a plain-text or gzip-compressed text file.

    Parameters
    ----------
    path
        Input file path.

    Returns
    -------
    TextIO
        Text-mode file handle.
    """
    if str(path).endswith(".gz"):
        return gzip.open(filename=path, mode="rt")
    return open(file=path, mode="rt", encoding="utf-8")


def strip_ensembl_version(*, identifier: object) -> str:
    """
    Remove a version suffix from an Ensembl identifier.

    Parameters
    ----------
    identifier
        Identifier value, usually like ENSG000001234.5.

    Returns
    -------
    str
        Identifier without version suffix.
    """
    text = str(identifier).strip()
    if text in {"", "nan", "None"}:
        return ""
    return text.split(".")[0]


def looks_like_ensembl_gene_id(*, value: object) -> bool:
    """
    Test whether a value resembles an Ensembl gene ID.

    Parameters
    ----------
    value
        Candidate identifier.

    Returns
    -------
    bool
        True if the value resembles an Ensembl gene ID.
    """
    return strip_ensembl_version(identifier=value).upper().startswith("ENSG")


def load_gene_symbol_map(
    *,
    gene_tpm_gct_path: Path,
    logger: logging.Logger,
) -> Tuple[Dict[str, str], Dict[str, Set[str]]]:
    """
    Load Ensembl-gene-ID to gene-symbol mapping from a GTEx gene GCT.

    This function streams only the first two columns from the GCT-style file,
    avoiding loading the large sample-level gene TPM matrix into memory.

    Parameters
    ----------
    gene_tpm_gct_path
        Path to a GTEx gene TPM GCT/GCT.GZ file.
    logger
        Logger instance.

    Returns
    -------
    Tuple[Dict[str, str], Dict[str, Set[str]]]
        First dictionary maps stripped Ensembl gene IDs to gene symbols.
        Second dictionary maps upper-case gene symbols to Ensembl gene IDs.
    """
    gene_id_to_symbol: Dict[str, str] = {}
    symbol_to_gene_ids: Dict[str, Set[str]] = {}

    with timed(logger=logger, label="Loading gene-symbol map from gene GCT"):
        with open_text_maybe_gzip(path=gene_tpm_gct_path) as handle:
            first = handle.readline()
            if first.startswith("#"):
                _ = handle.readline()
                header = handle.readline().rstrip("\n").split("\t")
            else:
                header = first.rstrip("\n").split("\t")

            if "Name" not in header or "Description" not in header:
                raise ValueError(
                    "Expected a GTEx GCT header containing 'Name' and "
                    f"'Description'. Observed first columns: {header[:5]}"
                )

            name_index = header.index("Name")
            description_index = header.index("Description")
            required_index = max(name_index, description_index)

            for line_number, line in enumerate(handle, start=4):
                parts = line.rstrip("\n").split("\t", required_index + 2)
                if len(parts) <= required_index:
                    logger.warning(
                        "Skipping malformed gene GCT line %d with %d fields",
                        line_number,
                        len(parts),
                    )
                    continue
                gene_id = strip_ensembl_version(identifier=parts[name_index])
                symbol = str(parts[description_index]).strip()
                if not gene_id or not symbol:
                    continue
                gene_id_to_symbol[gene_id] = symbol
                symbol_to_gene_ids.setdefault(symbol.upper(), set()).add(gene_id)

        logger.info("Loaded %d Ensembl gene ID to symbol mappings", len(gene_id_to_symbol))
        logger.info("Loaded %d unique gene symbols", len(symbol_to_gene_ids))

    return gene_id_to_symbol, symbol_to_gene_ids


def choose_first_existing_column(
    *,
    columns: Sequence[str],
    requested: Optional[str],
    candidates: Sequence[str],
) -> Optional[str]:
    """
    Choose a column by explicit request or from candidate names.

    Parameters
    ----------
    columns
        Available column names.
    requested
        Explicitly requested column name, or None.
    candidates
        Ordered candidate column names.

    Returns
    -------
    Optional[str]
        Selected column name, or None if no candidate is found.
    """
    if requested is not None:
        if requested not in columns:
            raise ValueError(
                f"Requested column '{requested}' was not found. "
                f"Available columns: {list(columns)}"
            )
        return requested

    lower_to_original = {column.lower(): column for column in columns}
    for candidate in candidates:
        if candidate.lower() in lower_to_original:
            return lower_to_original[candidate.lower()]
    return None


def load_requested_terms(
    *,
    target_genes: Sequence[str],
    genes_tsv: Optional[Path],
    gene_symbol_col: Optional[str],
    ensembl_gene_id_col: Optional[str],
    logger: logging.Logger,
) -> Tuple[str, ...]:
    """
    Load requested target terms from command-line values and/or a TSV.

    Parameters
    ----------
    target_genes
        Gene symbols and/or Ensembl gene IDs provided on the command line.
    genes_tsv
        Optional TSV containing target genes.
    gene_symbol_col
        Optional gene-symbol column in the TSV.
    ensembl_gene_id_col
        Optional Ensembl-gene-ID column in the TSV.
    logger
        Logger instance.

    Returns
    -------
    Tuple[str, ...]
        Deduplicated requested target terms in input order.
    """
    terms: List[str] = []
    terms.extend([str(gene).strip() for gene in target_genes if str(gene).strip()])

    if genes_tsv is not None:
        logger.info("Loading selected genes TSV: %s", genes_tsv)
        genes_df = pd.read_csv(filepath_or_buffer=genes_tsv, sep="\t", dtype=str)
        logger.info(
            "Selected genes TSV loaded: %d rows x %d columns",
            genes_df.shape[0],
            genes_df.shape[1],
        )

        symbol_col = choose_first_existing_column(
            columns=tuple(genes_df.columns),
            requested=gene_symbol_col,
            candidates=("gene_symbol", "gene_label", "symbol", "gene_name", "hgnc_symbol"),
        )
        ensembl_col = choose_first_existing_column(
            columns=tuple(genes_df.columns),
            requested=ensembl_gene_id_col,
            candidates=("ensembl_gene_id", "gene_id", "ensembl_id", "Name"),
        )

        if symbol_col is None and ensembl_col is None:
            raise ValueError(
                "Could not auto-detect a gene-symbol or Ensembl-gene-ID column "
                f"in {genes_tsv}. Columns: {list(genes_df.columns)}"
            )

        for column in [symbol_col, ensembl_col]:
            if column is None:
                continue
            values = genes_df[column].dropna().astype(str).str.strip()
            terms.extend([value for value in values if value])

    deduplicated = tuple(dict.fromkeys(terms))
    logger.info("Requested target terms after deduplication: %d", len(deduplicated))
    logger.info("Requested target terms: %s", "; ".join(deduplicated))
    return deduplicated


def resolve_targets(
    *,
    requested_terms: Sequence[str],
    symbol_to_gene_ids: Dict[str, Set[str]],
    logger: logging.Logger,
) -> ResolvedTargets:
    """
    Resolve requested gene symbols and Ensembl IDs to Ensembl gene IDs.

    Parameters
    ----------
    requested_terms
        Gene symbols and/or Ensembl gene IDs.
    symbol_to_gene_ids
        Mapping from upper-case symbols to stripped Ensembl gene IDs.
    logger
        Logger instance.

    Returns
    -------
    ResolvedTargets
        Resolved target-gene information.
    """
    requested_symbols: List[str] = []
    requested_gene_ids: List[str] = []
    resolved_gene_ids: List[str] = []
    unresolved_terms: List[str] = []

    with timed(logger=logger, label="Resolving target genes"):
        for term in requested_terms:
            cleaned = str(term).strip()
            if not cleaned:
                continue
            if looks_like_ensembl_gene_id(value=cleaned):
                gene_id = strip_ensembl_version(identifier=cleaned)
                requested_gene_ids.append(gene_id)
                resolved_gene_ids.append(gene_id)
                continue

            symbol = cleaned.upper()
            requested_symbols.append(symbol)
            matched_gene_ids = sorted(symbol_to_gene_ids.get(symbol, set()))
            if not matched_gene_ids:
                unresolved_terms.append(cleaned)
                continue
            resolved_gene_ids.extend(matched_gene_ids)

        resolved_gene_ids = list(dict.fromkeys(resolved_gene_ids))
        requested_symbols = list(dict.fromkeys(requested_symbols))
        requested_gene_ids = list(dict.fromkeys(requested_gene_ids))
        unresolved_terms = list(dict.fromkeys(unresolved_terms))

        logger.info("Resolved Ensembl gene IDs: %d", len(resolved_gene_ids))
        if resolved_gene_ids:
            logger.info("Resolved gene IDs: %s", "; ".join(resolved_gene_ids))
        if unresolved_terms:
            logger.warning(
                "Unresolved requested target term(s): %s",
                "; ".join(unresolved_terms),
            )

    if not resolved_gene_ids:
        raise ValueError("No target genes could be resolved to Ensembl gene IDs.")

    return ResolvedTargets(
        requested_terms=tuple(requested_terms),
        requested_symbols=tuple(requested_symbols),
        requested_gene_ids=tuple(requested_gene_ids),
        resolved_gene_ids=tuple(resolved_gene_ids),
        unresolved_terms=tuple(unresolved_terms),
    )


def find_matrix_header_line(*, matrix_path: Path, logger: logging.Logger) -> int:
    """
    Find the zero-based line index containing the transcript-matrix header.

    Parameters
    ----------
    matrix_path
        Matrix file path.
    logger
        Logger instance.

    Returns
    -------
    int
        Zero-based line index to use as the header after skipping previous rows.
    """
    with open_text_maybe_gzip(path=matrix_path) as handle:
        for line_index, line in enumerate(handle):
            fields = line.rstrip("\n").split("\t")
            lower_fields = {field.lower() for field in fields[:10]}
            has_transcript = bool(
                lower_fields.intersection({"transcript_id", "target_id", "name", "transcript"})
            )
            has_gene = bool(
                lower_fields.intersection({"gene_id", "gene", "description", "geneid"})
            )
            if has_transcript and has_gene:
                logger.info("Detected matrix header at zero-based line %d", line_index)
                logger.info("First header columns: %s", fields[:8])
                return line_index

    raise ValueError(f"Could not detect a transcript-matrix header in {matrix_path}")


def detect_transcript_schema(*, columns: Sequence[str]) -> TranscriptSchema:
    """
    Detect transcript ID, gene ID, metadata, and sample columns.

    Parameters
    ----------
    columns
        Matrix column names.

    Returns
    -------
    TranscriptSchema
        Detected transcript matrix schema.
    """
    lower_to_original = {column.lower(): column for column in columns}

    transcript_col = None
    for candidate in ("transcript_id", "target_id", "name", "transcript"):
        if candidate in lower_to_original:
            transcript_col = lower_to_original[candidate]
            break

    gene_col = None
    for candidate in ("gene_id", "gene", "description", "geneid"):
        if candidate in lower_to_original:
            gene_col = lower_to_original[candidate]
            break

    if transcript_col is None or gene_col is None:
        raise ValueError(
            "Could not detect transcript and gene columns from matrix header. "
            f"First columns were: {list(columns)[:10]}"
        )

    candidate_metadata = {
        transcript_col,
        gene_col,
        "transcript_name",
        "gene_name",
        "length",
        "effective_length",
        "Description",
    }
    metadata_cols = tuple(column for column in columns if column in candidate_metadata)
    sample_cols = tuple(column for column in columns if column not in metadata_cols)

    if not sample_cols:
        raise ValueError("No sample columns were detected in the transcript matrix.")

    return TranscriptSchema(
        transcript_col=transcript_col,
        gene_col=gene_col,
        metadata_cols=metadata_cols,
        sample_cols=sample_cols,
    )


def normalise_transcript_subset(
    *,
    dataframe: pd.DataFrame,
    schema: TranscriptSchema,
    gene_id_to_symbol: Dict[str, str],
    matrix_label: str,
) -> pd.DataFrame:
    """
    Add standard metadata columns to an extracted transcript subset.

    Parameters
    ----------
    dataframe
        Extracted transcript matrix rows.
    schema
        Detected transcript schema.
    gene_id_to_symbol
        Mapping from Ensembl gene IDs to symbols.
    matrix_label
        Label used in error messages.

    Returns
    -------
    pd.DataFrame
        Normalised subset with standard metadata columns first.
    """
    if dataframe.empty:
        raise ValueError(f"No rows were extracted from {matrix_label}.")

    raw_transcript = dataframe[schema.transcript_col].astype(str)
    raw_gene = dataframe[schema.gene_col].astype(str)
    gene_id = raw_gene.map(lambda value: strip_ensembl_version(identifier=value))

    metadata = pd.DataFrame(
        {
            "transcript_id_with_version": raw_transcript.values,
            "transcript_id": raw_transcript.map(
                lambda value: strip_ensembl_version(identifier=value)
            ).values,
            "gene_id_with_version": raw_gene.values,
            "gene_id": gene_id.values,
            "gene_symbol": gene_id.map(gene_id_to_symbol).fillna("").values,
        }
    )

    columns_to_drop = [
        column for column in [schema.transcript_col, schema.gene_col]
        if column in dataframe.columns
    ]
    expression_and_extra = dataframe.drop(columns=columns_to_drop).reset_index(drop=True)
    out = pd.concat(objs=[metadata, expression_and_extra], axis=1)
    return out


def extract_transcript_matrix_subset(
    *,
    matrix_path: Path,
    target_gene_ids: Set[str],
    gene_id_to_symbol: Dict[str, str],
    chunk_size: int,
    logger: logging.Logger,
    matrix_label: str,
) -> pd.DataFrame:
    """
    Extract rows for selected genes from a large transcript matrix.

    Parameters
    ----------
    matrix_path
        GTEx RSEM transcript matrix path.
    target_gene_ids
        Stripped Ensembl gene IDs to retain.
    gene_id_to_symbol
        Mapping from stripped Ensembl gene IDs to symbols.
    chunk_size
        Number of transcript rows to read per chunk.
    logger
        Logger instance.
    matrix_label
        Human-readable matrix label used in logs.

    Returns
    -------
    pd.DataFrame
        Extracted and normalised transcript matrix subset.
    """
    extracted_chunks: List[pd.DataFrame] = []
    header_line = find_matrix_header_line(matrix_path=matrix_path, logger=logger)
    schema: Optional[TranscriptSchema] = None
    total_rows = 0
    retained_rows = 0

    with timed(logger=logger, label=f"Scanning {matrix_label} transcript matrix"):
        reader = pd.read_csv(
            filepath_or_buffer=matrix_path,
            sep="\t",
            header=0,
            skiprows=header_line,
            chunksize=chunk_size,
            dtype=str,
            low_memory=False,
        )

        for chunk_index, chunk in enumerate(reader, start=1):
            if schema is None:
                schema = detect_transcript_schema(columns=tuple(chunk.columns))
                logger.info(
                    "%s schema: transcript_col=%s gene_col=%s sample_cols=%d",
                    matrix_label,
                    schema.transcript_col,
                    schema.gene_col,
                    len(schema.sample_cols),
                )

            total_rows += int(chunk.shape[0])
            gene_ids = chunk[schema.gene_col].map(
                lambda value: strip_ensembl_version(identifier=value)
            )
            mask = gene_ids.isin(target_gene_ids)
            if mask.any():
                retained = chunk.loc[mask, :].copy()
                extracted_chunks.append(retained)
                retained_rows += int(retained.shape[0])
                logger.info(
                    "%s chunk %d retained %d rows; total retained=%d",
                    matrix_label,
                    chunk_index,
                    retained.shape[0],
                    retained_rows,
                )

            if chunk_index % 100 == 0:
                logger.info(
                    "%s scan progress: chunks=%d rows=%d retained=%d",
                    matrix_label,
                    chunk_index,
                    total_rows,
                    retained_rows,
                )

    if schema is None:
        raise ValueError(f"No rows were read from {matrix_path}")
    if not extracted_chunks:
        raise ValueError(
            f"No transcript rows matched the requested gene IDs in {matrix_path}. "
            "Check gene symbols, Ensembl versions, and file annotation version."
        )

    subset = pd.concat(objs=extracted_chunks, axis=0, ignore_index=True)
    subset = normalise_transcript_subset(
        dataframe=subset,
        schema=schema,
        gene_id_to_symbol=gene_id_to_symbol,
        matrix_label=matrix_label,
    )

    logger.info(
        "%s subset extracted: %d transcripts x %d columns",
        matrix_label,
        subset.shape[0],
        subset.shape[1],
    )
    return subset


def load_sample_attributes(
    *,
    sample_attributes_path: Path,
    sample_id_col: str,
    tissue_col: str,
    logger: logging.Logger,
) -> pd.DataFrame:
    """
    Load GTEx sample attributes and retain sample-to-tissue mapping.

    Parameters
    ----------
    sample_attributes_path
        Path to SampleAttributesDS.txt.
    sample_id_col
        Sample ID column.
    tissue_col
        Tissue label column.
    logger
        Logger instance.

    Returns
    -------
    pd.DataFrame
        DataFrame containing sample ID and tissue label columns.
    """
    with timed(logger=logger, label="Loading GTEx sample attributes"):
        dataframe = pd.read_csv(
            filepath_or_buffer=sample_attributes_path,
            sep="\t",
            dtype=str,
            low_memory=False,
        )
        missing = [
            column for column in (sample_id_col, tissue_col)
            if column not in dataframe.columns
        ]
        if missing:
            raise ValueError(
                f"Missing required sample-attribute column(s): {missing}. "
                f"Available columns include: {list(dataframe.columns)[:30]}"
            )

        out = dataframe.loc[:, [sample_id_col, tissue_col]].copy()
        out[sample_id_col] = out[sample_id_col].astype(str)
        out[tissue_col] = out[tissue_col].astype(str)
        before = out.shape[0]
        out = out.drop_duplicates(subset=[sample_id_col], keep="first")
        if out.shape[0] != before:
            logger.warning(
                "Dropped %d duplicate sample-attribute row(s)",
                before - out.shape[0],
            )

        logger.info("Sample attributes loaded: %d samples", out.shape[0])
        logger.info(
            "Top tissue sample counts:\n%s",
            out[tissue_col].value_counts(dropna=False).head(10).to_string(),
        )
    return out


def get_sample_columns(
    *,
    dataframe: pd.DataFrame,
    sample_attributes: pd.DataFrame,
    sample_id_col: str,
    logger: logging.Logger,
) -> List[str]:
    """
    Identify expression columns that overlap GTEx sample IDs.

    Parameters
    ----------
    dataframe
        Transcript subset matrix.
    sample_attributes
        Sample-attribute DataFrame.
    sample_id_col
        Sample ID column in sample_attributes.
    logger
        Logger instance.

    Returns
    -------
    List[str]
        Sample columns in matrix order.
    """
    sample_ids = set(sample_attributes[sample_id_col].astype(str))
    metadata_cols = {
        "transcript_id_with_version",
        "transcript_id",
        "gene_id_with_version",
        "gene_id",
        "gene_symbol",
    }
    sample_cols = [
        column for column in dataframe.columns
        if column not in metadata_cols and column in sample_ids
    ]
    if not sample_cols:
        raise ValueError("No transcript matrix sample columns overlap SampleAttributesDS.")
    logger.info("Overlapping transcript matrix samples: %d", len(sample_cols))
    return sample_cols


def numeric_expression_matrix(
    *,
    dataframe: pd.DataFrame,
    sample_cols: Sequence[str],
    logger: logging.Logger,
) -> pd.DataFrame:
    """
    Convert selected sample columns to numeric expression values.

    Parameters
    ----------
    dataframe
        Transcript subset matrix.
    sample_cols
        Sample columns to convert.
    logger
        Logger instance.

    Returns
    -------
    pd.DataFrame
        Numeric expression matrix indexed by row position.
    """
    with timed(logger=logger, label="Converting selected expression values to numeric"):
        expression = dataframe.loc[:, list(sample_cols)].apply(
            pd.to_numeric,
            errors="coerce",
        )
        expression = expression.fillna(0.0)
        logger.info(
            "Numeric selected expression matrix: %d transcripts x %d samples",
            expression.shape[0],
            expression.shape[1],
        )
    return expression


def make_transcript_row_metadata(*, dataframe: pd.DataFrame) -> pd.DataFrame:
    """
    Extract standard transcript metadata columns from a subset matrix.

    Parameters
    ----------
    dataframe
        Transcript subset matrix.

    Returns
    -------
    pd.DataFrame
        Standard transcript metadata columns.
    """
    metadata_cols = [
        "transcript_id_with_version",
        "transcript_id",
        "gene_id_with_version",
        "gene_id",
        "gene_symbol",
    ]
    return dataframe.loc[:, metadata_cols].copy()


def make_transcript_sample_long(
    *,
    metadata: pd.DataFrame,
    expression: pd.DataFrame,
    sample_to_tissue: pd.Series,
    value_name: str,
    logger: logging.Logger,
) -> pd.DataFrame:
    """
    Convert a selected transcript matrix to sample-level long format.

    Parameters
    ----------
    metadata
        Transcript metadata with one row per transcript.
    expression
        Numeric transcript by sample matrix.
    sample_to_tissue
        Series mapping sample IDs to tissue labels.
    value_name
        Name for the value column.
    logger
        Logger instance.

    Returns
    -------
    pd.DataFrame
        Long-format sample-level transcript values.
    """
    with timed(logger=logger, label=f"Creating sample-level long table for {value_name}"):
        long_df = expression.copy()
        long_df.insert(loc=0, column="row_id", value=np.arange(long_df.shape[0]))
        long_df = long_df.melt(
            id_vars="row_id",
            var_name="sample_id",
            value_name=value_name,
        )
        out = long_df.merge(
            right=metadata.reset_index(drop=True).reset_index().rename(columns={"index": "row_id"}),
            how="left",
            on="row_id",
        )
        out["tissue"] = out["sample_id"].map(sample_to_tissue)
        out = out.drop(columns=["row_id"])
        ordered_cols = [
            "gene_symbol",
            "gene_id",
            "gene_id_with_version",
            "transcript_id",
            "transcript_id_with_version",
            "sample_id",
            "tissue",
            value_name,
        ]
        out = out.loc[:, ordered_cols]
        logger.info("Long table created: %d rows", out.shape[0])
    return out


def compute_tissue_medians(
    *,
    metadata: pd.DataFrame,
    expression: pd.DataFrame,
    sample_to_tissue: pd.Series,
    logger: logging.Logger,
) -> pd.DataFrame:
    """
    Compute median transcript TPM by GTEx tissue.

    Parameters
    ----------
    metadata
        Transcript metadata with one row per transcript.
    expression
        Numeric transcript by sample TPM matrix.
    sample_to_tissue
        Series mapping sample IDs to tissue labels.
    logger
        Logger instance.

    Returns
    -------
    pd.DataFrame
        Transcript rows by tissue columns, with metadata columns first.
    """
    with timed(logger=logger, label="Computing transcript tissue median TPM"):
        tissues = sample_to_tissue.loc[list(expression.columns)].astype(str)
        medians = (
            expression.T.assign(_tissue=tissues.values)
            .groupby("_tissue", dropna=False)
            .median(numeric_only=True)
            .T
        )
        out = pd.concat(
            objs=[metadata.reset_index(drop=True), medians.reset_index(drop=True)],
            axis=1,
        )
        logger.info(
            "Transcript tissue median TPM matrix: %d transcripts x %d tissues",
            medians.shape[0],
            medians.shape[1],
        )
    return out


def compute_isoform_usage(
    *,
    metadata: pd.DataFrame,
    expression: pd.DataFrame,
    epsilon: float,
    logger: logging.Logger,
) -> pd.DataFrame:
    """
    Compute sample-wise isoform usage for each transcript within its gene.

    Isoform usage is calculated for each sample as:

        transcript TPM / sum(TPM of all extracted transcripts for that gene)

    Parameters
    ----------
    metadata
        Transcript metadata with gene IDs.
    expression
        Numeric transcript by sample TPM matrix.
    epsilon
        Small constant to avoid division by zero.
    logger
        Logger instance.

    Returns
    -------
    pd.DataFrame
        Transcript by sample isoform-usage matrix.
    """
    with timed(logger=logger, label="Computing sample-wise isoform usage"):
        usage = pd.DataFrame(
            data=0.0,
            index=expression.index,
            columns=expression.columns,
        )
        for gene_id, row_index in metadata.groupby("gene_id", sort=False).groups.items():
            gene_rows = list(row_index)
            gene_total = expression.loc[gene_rows, :].sum(axis=0)
            usage.loc[gene_rows, :] = expression.loc[gene_rows, :].div(
                gene_total + epsilon,
                axis=1,
            )
            logger.debug(
                "Computed isoform usage for %s with %d transcripts",
                gene_id,
                len(gene_rows),
            )
        logger.info("Isoform usage matrix computed: %d x %d", usage.shape[0], usage.shape[1])
    return usage


def compute_usage_tissue_medians(
    *,
    metadata: pd.DataFrame,
    usage: pd.DataFrame,
    sample_to_tissue: pd.Series,
    logger: logging.Logger,
) -> pd.DataFrame:
    """
    Compute median isoform usage by GTEx tissue.

    Parameters
    ----------
    metadata
        Transcript metadata with one row per transcript.
    usage
        Transcript by sample isoform-usage matrix.
    sample_to_tissue
        Series mapping sample IDs to tissue labels.
    logger
        Logger instance.

    Returns
    -------
    pd.DataFrame
        Transcript rows by tissue columns, with metadata columns first.
    """
    with timed(logger=logger, label="Computing tissue median isoform usage"):
        tissues = sample_to_tissue.loc[list(usage.columns)].astype(str)
        medians = (
            usage.T.assign(_tissue=tissues.values)
            .groupby("_tissue", dropna=False)
            .median(numeric_only=True)
            .T
        )
        out = pd.concat(
            objs=[metadata.reset_index(drop=True), medians.reset_index(drop=True)],
            axis=1,
        )
        logger.info(
            "Isoform usage tissue median matrix: %d transcripts x %d tissues",
            medians.shape[0],
            medians.shape[1],
        )
    return out


def format_top_tissues(*, values: pd.Series, top_n: int = 5) -> str:
    """
    Format the top tissues for a transcript as semicolon-separated entries.

    Parameters
    ----------
    values
        Tissue-level values for one transcript.
    top_n
        Number of tissues to include.

    Returns
    -------
    str
        Semicolon-separated tissue=value entries.
    """
    top = values.sort_values(ascending=False).head(top_n)
    return ";".join([f"{tissue}={value:.4g}" for tissue, value in top.items()])


def compute_summary(
    *,
    metadata: pd.DataFrame,
    expression: pd.DataFrame,
    tissue_medians: pd.DataFrame,
    usage_tissue_medians: pd.DataFrame,
    sample_to_tissue: pd.Series,
    target_tissue: str,
    min_tpm_present: float,
    epsilon: float,
    logger: logging.Logger,
) -> pd.DataFrame:
    """
    Create a target-tissue transcript and isoform-usage summary.

    Parameters
    ----------
    metadata
        Transcript metadata with one row per transcript.
    expression
        Numeric transcript by sample TPM matrix.
    tissue_medians
        Transcript by tissue median TPM matrix with metadata columns.
    usage_tissue_medians
        Transcript by tissue median usage matrix with metadata columns.
    sample_to_tissue
        Series mapping sample IDs to tissue labels.
    target_tissue
        Target tissue for scoring.
    min_tpm_present
        TPM threshold for target-tissue presence fraction.
    epsilon
        Small constant to avoid division by zero.
    logger
        Logger instance.

    Returns
    -------
    pd.DataFrame
        Ranked transcript-level summary table.
    """
    metadata_cols = list(metadata.columns)
    tissue_cols = [column for column in tissue_medians.columns if column not in metadata_cols]
    usage_cols = [column for column in usage_tissue_medians.columns if column not in metadata_cols]

    if target_tissue not in tissue_cols:
        raise ValueError(
            f"Target tissue '{target_tissue}' was not found in transcript TPM tissues. "
            f"Available examples: {tissue_cols[:20]}"
        )
    if target_tissue not in usage_cols:
        raise ValueError(
            f"Target tissue '{target_tissue}' was not found in isoform usage tissues."
        )

    with timed(logger=logger, label="Computing target-tissue isoform summary"):
        tpm_values = tissue_medians.loc[:, tissue_cols].copy()
        usage_values = usage_tissue_medians.loc[:, usage_cols].copy()
        non_target_tpm = tpm_values.drop(columns=[target_tissue])
        non_target_usage = usage_values.drop(columns=[target_tissue])

        summary = metadata.copy().reset_index(drop=True)
        summary["n_transcripts_in_gene"] = summary.groupby("gene_id")[
            "transcript_id"
        ].transform("count")
        summary["target_tissue"] = target_tissue
        summary["target_median_tpm"] = tpm_values[target_tissue].values

        if non_target_tpm.shape[1] > 0:
            summary["max_non_target_tpm"] = non_target_tpm.max(axis=1).values
            summary["max_non_target_tpm_tissue"] = non_target_tpm.idxmax(axis=1).values
        else:
            summary["max_non_target_tpm"] = 0.0
            summary["max_non_target_tpm_tissue"] = ""

        summary["log2_target_vs_max_non_target_tpm"] = np.log2(
            (summary["target_median_tpm"] + epsilon)
            / (summary["max_non_target_tpm"] + epsilon)
        )
        summary["target_is_max_tpm_tissue"] = (
            tpm_values.idxmax(axis=1).values == target_tissue
        ).astype(int)
        summary["top_tissues_by_median_tpm"] = tpm_values.apply(
            lambda row: format_top_tissues(values=row),
            axis=1,
        ).values

        summary["target_median_isoform_usage"] = usage_values[target_tissue].values
        if non_target_usage.shape[1] > 0:
            summary["max_non_target_isoform_usage"] = non_target_usage.max(axis=1).values
            summary["max_non_target_isoform_usage_tissue"] = non_target_usage.idxmax(axis=1).values
        else:
            summary["max_non_target_isoform_usage"] = 0.0
            summary["max_non_target_isoform_usage_tissue"] = ""

        summary["log2_target_vs_max_non_target_isoform_usage"] = np.log2(
            (summary["target_median_isoform_usage"] + epsilon)
            / (summary["max_non_target_isoform_usage"] + epsilon)
        )
        summary["target_is_max_usage_tissue"] = (
            usage_values.idxmax(axis=1).values == target_tissue
        ).astype(int)
        summary["top_tissues_by_median_isoform_usage"] = usage_values.apply(
            lambda row: format_top_tissues(values=row),
            axis=1,
        ).values

        target_samples = [
            sample for sample in expression.columns
            if sample_to_tissue.loc[sample] == target_tissue
        ]
        if not target_samples:
            raise ValueError(f"No samples labelled as target tissue '{target_tissue}'.")
        target_present_fraction = (
            expression.loc[:, target_samples] >= min_tpm_present
        ).sum(axis=1) / float(len(target_samples))
        summary["target_present_fraction_tpm_ge_threshold"] = target_present_fraction.values
        summary["target_present_tpm_threshold"] = min_tpm_present

        summary["target_tpm_rank_within_gene"] = summary.groupby("gene_id")[
            "target_median_tpm"
        ].rank(method="dense", ascending=False).astype(int)
        summary["target_usage_rank_within_gene"] = summary.groupby("gene_id")[
            "target_median_isoform_usage"
        ].rank(method="dense", ascending=False).astype(int)

        summary = summary.sort_values(
            by=[
                "target_is_max_usage_tissue",
                "log2_target_vs_max_non_target_isoform_usage",
                "target_median_isoform_usage",
                "target_median_tpm",
            ],
            ascending=[False, False, False, False],
        )
        logger.info("Summary created for %d transcripts", summary.shape[0])
    return summary


def write_tsv(
    *,
    dataframe: pd.DataFrame,
    path: Path,
    logger: logging.Logger,
    index: bool = False,
) -> None:
    """
    Write a DataFrame as TSV/TSV.GZ with logging.

    Parameters
    ----------
    dataframe
        DataFrame to write.
    path
        Output path.
    logger
        Logger instance.
    index
        Whether to include the DataFrame index.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    logger.info("Writing %d rows x %d columns: %s", dataframe.shape[0], dataframe.shape[1], path)
    dataframe.to_csv(path_or_buf=path, sep="\t", index=index)


def run(*, config: Config, logger: logging.Logger) -> None:
    """
    Run the GTEx transcript-level first-pass isoform workflow.

    Parameters
    ----------
    config
        Pipeline configuration.
    logger
        Logger instance.
    """
    config.out_dir.mkdir(parents=True, exist_ok=True)
    logger.info("Starting GTEx transcript-level isoform first-pass workflow")
    logger.info("Transcript TPM path: %s", config.transcript_tpm_path)
    logger.info("Gene TPM GCT path: %s", config.gene_tpm_gct_path)
    logger.info("Sample attributes path: %s", config.sample_attributes_path)
    logger.info("Output directory: %s", config.out_dir)
    logger.info("Output prefix: %s", config.out_prefix)

    gene_id_to_symbol, symbol_to_gene_ids = load_gene_symbol_map(
        gene_tpm_gct_path=config.gene_tpm_gct_path,
        logger=logger,
    )
    requested_terms = load_requested_terms(
        target_genes=config.target_genes,
        genes_tsv=config.genes_tsv,
        gene_symbol_col=config.gene_symbol_col,
        ensembl_gene_id_col=config.ensembl_gene_id_col,
        logger=logger,
    )
    resolved = resolve_targets(
        requested_terms=requested_terms,
        symbol_to_gene_ids=symbol_to_gene_ids,
        logger=logger,
    )
    target_gene_ids = set(resolved.resolved_gene_ids)

    sample_attributes = load_sample_attributes(
        sample_attributes_path=config.sample_attributes_path,
        sample_id_col=config.sample_id_col,
        tissue_col=config.tissue_col,
        logger=logger,
    )
    sample_to_tissue = sample_attributes.set_index(config.sample_id_col)[config.tissue_col]

    transcript_tpm_subset = extract_transcript_matrix_subset(
        matrix_path=config.transcript_tpm_path,
        target_gene_ids=target_gene_ids,
        gene_id_to_symbol=gene_id_to_symbol,
        chunk_size=config.chunk_size,
        logger=logger,
        matrix_label="TPM",
    )
    tpm_subset_path = config.out_dir / f"{config.out_prefix}.transcript_subset_tpm_matrix.tsv.gz"
    write_tsv(dataframe=transcript_tpm_subset, path=tpm_subset_path, logger=logger)

    sample_cols = get_sample_columns(
        dataframe=transcript_tpm_subset,
        sample_attributes=sample_attributes,
        sample_id_col=config.sample_id_col,
        logger=logger,
    )
    transcript_tpm_subset = transcript_tpm_subset.loc[:, [
        "transcript_id_with_version",
        "transcript_id",
        "gene_id_with_version",
        "gene_id",
        "gene_symbol",
    ] + sample_cols]

    metadata = make_transcript_row_metadata(dataframe=transcript_tpm_subset)
    expression = numeric_expression_matrix(
        dataframe=transcript_tpm_subset,
        sample_cols=sample_cols,
        logger=logger,
    )

    sample_to_tissue = sample_to_tissue.loc[sample_cols]

    sample_long = make_transcript_sample_long(
        metadata=metadata,
        expression=expression,
        sample_to_tissue=sample_to_tissue,
        value_name="transcript_tpm",
        logger=logger,
    )
    sample_long_path = config.out_dir / f"{config.out_prefix}.transcript_sample_tpm_long.tsv.gz"
    write_tsv(dataframe=sample_long, path=sample_long_path, logger=logger)

    tissue_medians = compute_tissue_medians(
        metadata=metadata,
        expression=expression,
        sample_to_tissue=sample_to_tissue,
        logger=logger,
    )
    tissue_medians_path = config.out_dir / f"{config.out_prefix}.transcript_tissue_median_tpm.tsv"
    write_tsv(dataframe=tissue_medians, path=tissue_medians_path, logger=logger)

    usage = compute_isoform_usage(
        metadata=metadata,
        expression=expression,
        epsilon=config.epsilon,
        logger=logger,
    )
    usage_tissue_medians = compute_usage_tissue_medians(
        metadata=metadata,
        usage=usage,
        sample_to_tissue=sample_to_tissue,
        logger=logger,
    )
    usage_path = (
        config.out_dir
        / f"{config.out_prefix}.transcript_isoform_usage_median_by_tissue.tsv"
    )
    write_tsv(dataframe=usage_tissue_medians, path=usage_path, logger=logger)

    summary = compute_summary(
        metadata=metadata,
        expression=expression,
        tissue_medians=tissue_medians,
        usage_tissue_medians=usage_tissue_medians,
        sample_to_tissue=sample_to_tissue,
        target_tissue=config.target_tissue,
        min_tpm_present=config.min_tpm_present,
        epsilon=config.epsilon,
        logger=logger,
    )
    summary_path = config.out_dir / f"{config.out_prefix}.transcript_testis_isoform_summary.tsv"
    write_tsv(dataframe=summary, path=summary_path, logger=logger)

    if config.expected_count_path is not None:
        if not config.expected_count_path.exists():
            logger.warning(
                "Expected-count path was provided but does not exist yet: %s",
                config.expected_count_path,
            )
        else:
            expected_count_subset = extract_transcript_matrix_subset(
                matrix_path=config.expected_count_path,
                target_gene_ids=target_gene_ids,
                gene_id_to_symbol=gene_id_to_symbol,
                chunk_size=config.chunk_size,
                logger=logger,
                matrix_label="expected-count",
            )
            count_subset_path = (
                config.out_dir
                / f"{config.out_prefix}.transcript_subset_expected_count_matrix.tsv.gz"
            )
            write_tsv(dataframe=expected_count_subset, path=count_subset_path, logger=logger)

    logger.info("Workflow finished successfully")
    logger.info("Main summary: %s", summary_path)


def main() -> None:
    """
    Command-line entry point.
    """
    config = parse_args()
    logger = setup_logging(log_level=config.log_level, log_path=config.log_path)
    try:
        run(config=config, logger=logger)
    except Exception as exc:
        logger.exception("Fatal error: %s", str(exc))
        raise


if __name__ == "__main__":
    main()
