#!/usr/bin/env python3
"""
Run a transcriptome-wide GTEx isoform tissue-selectivity screen.

This script is designed for the transcript-level fertility-genomics workflow.
It loads the full GTEx RSEM transcript TPM matrix, maps samples to tissues,
computes transcript-level median TPM by tissue, computes within-gene isoform
usage by tissue, and ranks transcripts for target-tissue preferential usage.

The main biological question is:

    Which genes may be broadly expressed at gene level but have one transcript
    isoform that is preferentially expressed or preferentially used in testis?

The workflow intentionally writes tab-separated outputs only. Large all-gene
outputs are written as TSV.GZ files by default.

Expected core inputs
--------------------
1. GTEx transcript TPM matrix, for example:
   GTEx_Analysis_2025-08-22_v11_RSEMv1.3.3_transcripts_tpm.txt.gz

2. GTEx gene TPM GCT file, used only to map Ensembl gene IDs to gene symbols:
   GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_tpm.gct.gz

3. GTEx sample annotations:
   GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt

Primary outputs
---------------
<prefix>.transcript_tissue_median_tpm.tsv.gz
    Transcript-by-tissue median TPM matrix.

<prefix>.transcript_isoform_usage_median_by_tissue.tsv.gz
    Transcript-by-tissue median isoform-usage matrix. Isoform usage is
    calculated per sample as transcript TPM divided by total TPM for that gene.

<prefix>.gene_tissue_median_tpm.tsv.gz
    Gene-by-tissue median TPM matrix calculated from summed transcript TPM.

<prefix>.transcript_target_tissue_isoform_summary.tsv.gz
    Searchable transcriptome-wide summary ranked for target-tissue TPM and
    isoform usage.

<prefix>.candidate_target_tissue_isoforms.tsv
    Filtered candidate transcript isoforms with preferential target-tissue
    usage.

<prefix>.best_candidate_isoform_per_gene.tsv
    One best candidate transcript per gene.
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
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class Config:
    """
    Configuration for the transcriptome-wide isoform screen.

    Attributes
    ----------
    transcript_tpm_path
        Path to the GTEx transcript TPM matrix.
    gene_tpm_gct_path
        Path to the GTEx gene TPM GCT file, used to map Ensembl gene IDs to
        gene symbols.
    sample_attributes_path
        Path to GTEx SampleAttributesDS.txt.
    out_dir
        Output directory.
    out_prefix
        Prefix for output filenames.
    target_tissue
        Tissue of interest for target-tissue scoring.
    sample_id_col
        Sample ID column in the sample attributes file.
    tissue_col
        Tissue label column in the sample attributes file.
    float_dtype
        Floating-point dtype used for large expression arrays.
    min_tpm_present
        TPM threshold used to calculate target-tissue detection fraction.
    min_target_tpm_candidate
        Minimum target-tissue median TPM for candidate isoforms.
    min_target_usage_candidate
        Minimum target-tissue median isoform usage for candidate isoforms.
    min_log2_usage_ratio_candidate
        Minimum log2 target versus max non-target isoform-usage ratio for
        candidate isoforms.
    max_gene_log2_ratio_for_rescue
        Maximum gene-level log2 target versus max non-target ratio for flagging
        genes where isoform-level selectivity may rescue a non-selective gene.
    top_n_candidate_tissues
        Number of top TPM and usage tissues to list for candidate outputs.
    write_tissue_matrices
        Whether to write the transcript and gene tissue matrix outputs.
    log_path
        Optional path to a log file.
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
    target_tissue: str
    sample_id_col: str
    tissue_col: str
    float_dtype: str
    min_tpm_present: float
    min_target_tpm_candidate: float
    min_target_usage_candidate: float
    min_log2_usage_ratio_candidate: float
    max_gene_log2_ratio_for_rescue: float
    top_n_candidate_tissues: int
    write_tissue_matrices: bool
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
        Expression sample columns with matching sample annotations.
    ignored_expression_like_cols
        Columns that look like expression samples but were not present in the
        sample attributes file.
    """

    transcript_col: str
    gene_col: str
    metadata_cols: Tuple[str, ...]
    sample_cols: Tuple[str, ...]
    ignored_expression_like_cols: Tuple[str, ...]


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
    logger = logging.getLogger("gtex_transcriptome_isoform_screen")
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
        Parsed configuration.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Run a transcriptome-wide GTEx transcript isoform-usage screen."
        )
    )
    parser.add_argument(
        "--transcript_tpm",
        required=True,
        type=Path,
        help="GTEx RSEM transcript TPM matrix, optionally gzipped.",
    )
    parser.add_argument(
        "--gene_tpm_gct",
        required=True,
        type=Path,
        help="GTEx gene TPM GCT/GCT.GZ file used to map gene IDs to symbols.",
    )
    parser.add_argument(
        "--sample_attributes",
        required=True,
        type=Path,
        help="GTEx SampleAttributesDS.txt file.",
    )
    parser.add_argument(
        "--out_dir",
        required=True,
        type=Path,
        help="Output directory.",
    )
    parser.add_argument(
        "--out_prefix",
        default="gtex_v11_transcriptome_isoform_screen",
        help="Output filename prefix.",
    )
    parser.add_argument(
        "--target_tissue",
        default="Testis",
        help="Target tissue to score. Default: Testis.",
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
        "--float_dtype",
        default="float32",
        choices=("float32", "float64"),
        help="Floating dtype for expression arrays. Default: float32.",
    )
    parser.add_argument(
        "--min_tpm_present",
        default=1.0,
        type=float,
        help="TPM threshold used for target-tissue detection fraction.",
    )
    parser.add_argument(
        "--min_target_tpm_candidate",
        default=1.0,
        type=float,
        help="Minimum target-tissue median TPM for candidate isoforms.",
    )
    parser.add_argument(
        "--min_target_usage_candidate",
        default=0.25,
        type=float,
        help="Minimum target-tissue median isoform usage for candidates.",
    )
    parser.add_argument(
        "--min_log2_usage_ratio_candidate",
        default=1.0,
        type=float,
        help=(
            "Minimum log2 target versus max non-target isoform-usage ratio "
            "for candidates. Default 1 equals two-fold."
        ),
    )
    parser.add_argument(
        "--max_gene_log2_ratio_for_rescue",
        default=1.0,
        type=float,
        help=(
            "Maximum gene-level log2 target versus max non-target TPM ratio "
            "for flagging broad-gene isoform-rescue candidates."
        ),
    )
    parser.add_argument(
        "--top_n_candidate_tissues",
        default=5,
        type=int,
        help="Number of top tissues listed in candidate outputs.",
    )
    parser.add_argument(
        "--skip_tissue_matrices",
        action="store_true",
        help=(
            "Do not write the large transcript/gene tissue matrices. The "
            "searchable summary and candidate files are still written."
        ),
    )
    parser.add_argument(
        "--log_path",
        type=Path,
        default=None,
        help="Optional log-file path.",
    )
    parser.add_argument(
        "--log_level",
        default="INFO",
        help="Logging level. Default: INFO.",
    )

    args = parser.parse_args(argv)
    return Config(
        transcript_tpm_path=args.transcript_tpm,
        gene_tpm_gct_path=args.gene_tpm_gct,
        sample_attributes_path=args.sample_attributes,
        out_dir=args.out_dir,
        out_prefix=args.out_prefix,
        target_tissue=args.target_tissue,
        sample_id_col=args.sample_id_col,
        tissue_col=args.tissue_col,
        float_dtype=args.float_dtype,
        min_tpm_present=args.min_tpm_present,
        min_target_tpm_candidate=args.min_target_tpm_candidate,
        min_target_usage_candidate=args.min_target_usage_candidate,
        min_log2_usage_ratio_candidate=args.min_log2_usage_ratio_candidate,
        max_gene_log2_ratio_for_rescue=args.max_gene_log2_ratio_for_rescue,
        top_n_candidate_tissues=args.top_n_candidate_tissues,
        write_tissue_matrices=not args.skip_tissue_matrices,
        log_path=args.log_path,
        log_level=args.log_level,
    )


def open_text_maybe_gzip(*, path: Path):
    """
    Open a plain-text or gzip-compressed text file.

    Parameters
    ----------
    path
        Input path.

    Returns
    -------
    file object
        Text-mode file handle.
    """
    if str(path).endswith(".gz"):
        return gzip.open(filename=path, mode="rt")
    return path.open(mode="rt")


def strip_ensembl_version(*, identifier: object) -> str:
    """
    Remove an Ensembl version suffix from an identifier.

    Parameters
    ----------
    identifier
        Identifier value.

    Returns
    -------
    str
        Identifier with the trailing dot-version removed.
    """
    value = str(identifier).strip()
    if not value or value.lower() == "nan":
        return ""
    return value.split(".", maxsplit=1)[0]


def find_matrix_header_line(*, matrix_path: Path, logger: logging.Logger) -> int:
    """
    Find the zero-based header line for a transcript matrix.

    Parameters
    ----------
    matrix_path
        Matrix path.
    logger
        Logger instance.

    Returns
    -------
    int
        Zero-based header line index.
    """
    transcript_candidates = {"transcript_id", "target_id", "name", "transcript"}
    gene_candidates = {"gene_id", "gene", "description", "geneid"}

    with open_text_maybe_gzip(path=matrix_path) as handle:
        for line_index, line in enumerate(handle):
            fields = line.rstrip("\n").split("\t")
            lower_fields = {field.lower() for field in fields[:10]}
            has_transcript = bool(lower_fields.intersection(transcript_candidates))
            has_gene = bool(lower_fields.intersection(gene_candidates))
            if has_transcript and has_gene:
                logger.info("Detected matrix header at zero-based line %d", line_index)
                logger.info("First header columns: %s", fields[:8])
                return line_index

    raise ValueError(f"Could not detect transcript matrix header in {matrix_path}")


def detect_transcript_schema(
    *,
    columns: Sequence[str],
    sample_ids: Iterable[str],
) -> TranscriptSchema:
    """
    Detect transcript, gene, metadata, and annotated sample columns.

    Parameters
    ----------
    columns
        Matrix column names.
    sample_ids
        Sample IDs from the GTEx sample attributes file.

    Returns
    -------
    TranscriptSchema
        Detected matrix schema.
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

    sample_id_set = set(sample_ids)
    sample_cols = tuple(column for column in columns if column in sample_id_set)
    metadata_cols = tuple(column for column in columns if column not in sample_cols)

    ignored_expression_like_cols = tuple(
        column
        for column in columns
        if column.startswith("GTEX-") and column not in sample_id_set
    )

    if not sample_cols:
        raise ValueError(
            "No transcript matrix sample columns matched SampleAttributesDS.txt."
        )

    return TranscriptSchema(
        transcript_col=transcript_col,
        gene_col=gene_col,
        metadata_cols=metadata_cols,
        sample_cols=sample_cols,
        ignored_expression_like_cols=ignored_expression_like_cols,
    )


def load_sample_attributes(
    *,
    path: Path,
    sample_id_col: str,
    tissue_col: str,
    logger: logging.Logger,
) -> pd.DataFrame:
    """
    Load GTEx sample attributes and retain sample-to-tissue mapping.

    Parameters
    ----------
    path
        Sample attributes path.
    sample_id_col
        Sample ID column name.
    tissue_col
        Tissue column name.
    logger
        Logger instance.

    Returns
    -------
    pandas.DataFrame
        DataFrame containing unique sample IDs and tissue labels.
    """
    attributes = pd.read_csv(filepath_or_buffer=path, sep="\t", dtype=str)
    missing = [
        column for column in (sample_id_col, tissue_col) if column not in attributes.columns
    ]
    if missing:
        raise ValueError(
            f"Sample attributes file is missing required columns: {missing}"
        )

    attributes = attributes[[sample_id_col, tissue_col]].dropna()
    attributes = attributes.drop_duplicates(subset=[sample_id_col], keep="first")
    attributes[sample_id_col] = attributes[sample_id_col].astype(str)
    attributes[tissue_col] = attributes[tissue_col].astype(str)

    logger.info(
        "Loaded %d unique GTEx sample annotations across %d tissues",
        attributes.shape[0],
        attributes[tissue_col].nunique(),
    )
    return attributes


def load_gene_symbol_map(*, gct_path: Path, logger: logging.Logger) -> Dict[str, str]:
    """
    Load an Ensembl gene ID to gene-symbol map from a GTEx gene TPM GCT file.

    Parameters
    ----------
    gct_path
        GTEx gene TPM GCT path.
    logger
        Logger instance.

    Returns
    -------
    dict[str, str]
        Mapping from version-stripped Ensembl gene ID to gene symbol.
    """
    header_line = None
    with open_text_maybe_gzip(path=gct_path) as handle:
        for line_index, line in enumerate(handle):
            fields = line.rstrip("\n").split("\t")
            lower_fields = {field.lower() for field in fields[:5]}
            if {"name", "description"}.issubset(lower_fields):
                header_line = line_index
                break

    if header_line is None:
        raise ValueError(f"Could not detect a GCT header in {gct_path}")

    gene_table = pd.read_csv(
        filepath_or_buffer=gct_path,
        sep="\t",
        skiprows=header_line,
        usecols=["Name", "Description"],
        dtype=str,
    )
    gene_table = gene_table.dropna(subset=["Name"])
    gene_table["gene_id"] = gene_table["Name"].map(
        lambda value: strip_ensembl_version(identifier=value)
    )
    gene_table["gene_symbol"] = gene_table["Description"].fillna("").astype(str)
    gene_table = gene_table.drop_duplicates(subset=["gene_id"], keep="first")

    mapping = dict(zip(gene_table["gene_id"], gene_table["gene_symbol"]))
    logger.info("Loaded %d gene ID to symbol mappings", len(mapping))
    return mapping


def preview_transcript_header(
    *,
    path: Path,
    header_line: int,
    logger: logging.Logger,
) -> List[str]:
    """
    Read only the transcript matrix header.

    Parameters
    ----------
    path
        Transcript matrix path.
    header_line
        Zero-based header line index.
    logger
        Logger instance.

    Returns
    -------
    list[str]
        Header column names.
    """
    header = pd.read_csv(
        filepath_or_buffer=path,
        sep="\t",
        skiprows=header_line,
        nrows=0,
    )
    columns = list(header.columns)
    logger.info("Transcript matrix has %d columns", len(columns))
    return columns


def load_full_transcript_tpm(
    *,
    path: Path,
    schema: TranscriptSchema,
    header_line: int,
    float_dtype: str,
    logger: logging.Logger,
) -> pd.DataFrame:
    """
    Load the full GTEx transcript TPM matrix into memory.

    Parameters
    ----------
    path
        Transcript TPM path.
    schema
        Detected schema.
    header_line
        Zero-based header line index.
    float_dtype
        Floating dtype for sample-expression columns.
    logger
        Logger instance.

    Returns
    -------
    pandas.DataFrame
        Full matrix with metadata columns and annotated sample columns.
    """
    dtype_map = {column: str for column in schema.metadata_cols}
    dtype_map.update({column: float_dtype for column in schema.sample_cols})

    usecols = list(schema.metadata_cols) + list(schema.sample_cols)
    logger.info("Loading %d metadata columns", len(schema.metadata_cols))
    logger.info("Loading %d annotated sample expression columns", len(schema.sample_cols))

    dataframe = pd.read_csv(
        filepath_or_buffer=path,
        sep="\t",
        skiprows=header_line,
        usecols=usecols,
        dtype=dtype_map,
    )
    logger.info(
        "Loaded full transcript TPM matrix with %d transcripts and %d columns",
        dataframe.shape[0],
        dataframe.shape[1],
    )
    return dataframe


def normalise_transcript_dataframe(
    *,
    dataframe: pd.DataFrame,
    schema: TranscriptSchema,
    gene_id_to_symbol: Dict[str, str],
    sample_cols: Sequence[str],
    logger: logging.Logger,
) -> pd.DataFrame:
    """
    Add standard transcript and gene metadata columns to the matrix.

    Parameters
    ----------
    dataframe
        Loaded transcript TPM matrix.
    schema
        Detected schema.
    gene_id_to_symbol
        Mapping from version-stripped Ensembl gene IDs to gene symbols.
    sample_cols
        Annotated sample expression columns.
    logger
        Logger instance.

    Returns
    -------
    pandas.DataFrame
        Normalised matrix with standard metadata columns first.
    """
    transcript_raw = dataframe[schema.transcript_col].astype(str)
    gene_raw = dataframe[schema.gene_col].astype(str)
    gene_id = gene_raw.map(lambda value: strip_ensembl_version(identifier=value))
    transcript_id = transcript_raw.map(
        lambda value: strip_ensembl_version(identifier=value)
    )

    metadata = pd.DataFrame(
        {
            "transcript_id_with_version": transcript_raw.values,
            "transcript_id": transcript_id.values,
            "gene_id_with_version": gene_raw.values,
            "gene_id": gene_id.values,
            "gene_symbol": gene_id.map(gene_id_to_symbol).fillna("").values,
        }
    )

    extra_metadata_cols = [
        column
        for column in schema.metadata_cols
        if column not in {schema.transcript_col, schema.gene_col}
    ]
    extra_metadata = dataframe[extra_metadata_cols].reset_index(drop=True)
    expression = dataframe[list(sample_cols)].reset_index(drop=True)

    output = pd.concat(
        objs=[metadata, extra_metadata, expression],
        axis=1,
    )
    logger.info(
        "Normalised transcript matrix: %d rows, %d samples, %d extra metadata columns",
        output.shape[0],
        len(sample_cols),
        len(extra_metadata_cols),
    )
    return output


def make_expression_array(
    *,
    dataframe: pd.DataFrame,
    sample_cols: Sequence[str],
    float_dtype: str,
    logger: logging.Logger,
) -> np.ndarray:
    """
    Convert sample columns to a dense NumPy expression array.

    Parameters
    ----------
    dataframe
        Normalised transcript matrix.
    sample_cols
        Sample expression columns.
    float_dtype
        Floating dtype.
    logger
        Logger instance.

    Returns
    -------
    numpy.ndarray
        Expression array with shape transcripts by samples.
    """
    expression = dataframe[list(sample_cols)]
    array = expression.to_numpy(dtype=np.dtype(float_dtype), copy=False)
    if np.isnan(array).any():
        logger.warning("Expression array contains NaN values; replacing with 0.0")
        array = np.nan_to_num(array, nan=0.0)
    logger.info(
        "Expression array shape is %s and dtype is %s",
        array.shape,
        array.dtype,
    )
    return array


def make_tissue_sample_indices(
    *,
    sample_cols: Sequence[str],
    sample_attributes: pd.DataFrame,
    sample_id_col: str,
    tissue_col: str,
    logger: logging.Logger,
) -> Dict[str, np.ndarray]:
    """
    Build tissue to expression-column-index mappings.

    Parameters
    ----------
    sample_cols
        Expression sample columns in matrix order.
    sample_attributes
        Sample attributes with sample IDs and tissues.
    sample_id_col
        Sample ID column name.
    tissue_col
        Tissue column name.
    logger
        Logger instance.

    Returns
    -------
    dict[str, numpy.ndarray]
        Tissue labels mapped to integer column indices.
    """
    sample_to_tissue = dict(
        zip(sample_attributes[sample_id_col], sample_attributes[tissue_col])
    )
    tissue_to_indices: Dict[str, List[int]] = {}
    missing_samples = []

    for sample_index, sample_id in enumerate(sample_cols):
        tissue = sample_to_tissue.get(sample_id)
        if tissue is None:
            missing_samples.append(sample_id)
            continue
        tissue_to_indices.setdefault(tissue, []).append(sample_index)

    if missing_samples:
        logger.warning(
            "%d matrix samples were not present in the sample attributes file",
            len(missing_samples),
        )

    output = {
        tissue: np.array(indices, dtype=np.int64)
        for tissue, indices in sorted(tissue_to_indices.items())
    }
    logger.info("Mapped expression samples to %d tissues", len(output))
    return output


def compute_tissue_medians(
    *,
    values: np.ndarray,
    row_metadata: pd.DataFrame,
    tissue_to_indices: Dict[str, np.ndarray],
    logger: logging.Logger,
    value_label: str,
) -> pd.DataFrame:
    """
    Compute row-wise median values for each tissue.

    Parameters
    ----------
    values
        Row by sample array.
    row_metadata
        Metadata for rows in values.
    tissue_to_indices
        Tissue to sample-index mapping.
    logger
        Logger instance.
    value_label
        Label for logging.

    Returns
    -------
    pandas.DataFrame
        Metadata plus one median-value column per tissue.
    """
    medians: Dict[str, np.ndarray] = {}
    for tissue, indices in tissue_to_indices.items():
        if len(indices) == 0:
            continue
        logger.info(
            "Computing %s median for tissue %s using %d samples",
            value_label,
            tissue,
            len(indices),
        )
        medians[tissue] = np.median(values[:, indices], axis=1)

    median_frame = pd.DataFrame(medians)
    return pd.concat(
        objs=[row_metadata.reset_index(drop=True), median_frame.reset_index(drop=True)],
        axis=1,
    )


def compute_gene_expression_array(
    *,
    expression_array: np.ndarray,
    gene_ids: Sequence[str],
    float_dtype: str,
    logger: logging.Logger,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Sum transcript TPM values to gene-level expression per sample.

    Parameters
    ----------
    expression_array
        Transcript by sample expression array.
    gene_ids
        Gene ID for each transcript row.
    float_dtype
        Floating dtype for the output array.
    logger
        Logger instance.

    Returns
    -------
    tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]
        Unique gene IDs, transcript-to-gene inverse indices, and gene by sample
        expression array.
    """
    unique_gene_ids, inverse = np.unique(np.asarray(gene_ids), return_inverse=True)
    gene_array = np.zeros(
        shape=(len(unique_gene_ids), expression_array.shape[1]),
        dtype=np.dtype(float_dtype),
    )
    logger.info(
        "Summing transcript TPM to %d unique genes across %d samples",
        len(unique_gene_ids),
        expression_array.shape[1],
    )
    np.add.at(gene_array, inverse, expression_array)
    logger.info("Gene expression array shape is %s", gene_array.shape)
    return unique_gene_ids, inverse, gene_array


def make_gene_metadata(
    *,
    gene_ids: Sequence[str],
    gene_id_to_symbol: Dict[str, str],
) -> pd.DataFrame:
    """
    Make gene-level metadata rows.

    Parameters
    ----------
    gene_ids
        Version-stripped Ensembl gene IDs.
    gene_id_to_symbol
        Mapping from Ensembl gene IDs to symbols.

    Returns
    -------
    pandas.DataFrame
        Gene metadata table.
    """
    return pd.DataFrame(
        {
            "gene_id": list(gene_ids),
            "gene_symbol": [gene_id_to_symbol.get(gene_id, "") for gene_id in gene_ids],
        }
    )


def compute_isoform_usage_tissue_medians(
    *,
    expression_array: np.ndarray,
    gene_expression_array: np.ndarray,
    transcript_gene_inverse: np.ndarray,
    row_metadata: pd.DataFrame,
    tissue_to_indices: Dict[str, np.ndarray],
    epsilon: float,
    logger: logging.Logger,
) -> pd.DataFrame:
    """
    Compute transcript isoform-usage medians per tissue.

    Isoform usage is calculated per sample as transcript TPM divided by total
    gene TPM in that same sample.

    Parameters
    ----------
    expression_array
        Transcript by sample expression array.
    gene_expression_array
        Gene by sample expression array.
    transcript_gene_inverse
        Integer mapping from transcript rows to gene rows.
    row_metadata
        Transcript row metadata.
    tissue_to_indices
        Tissue to sample-index mapping.
    epsilon
        Small constant for safe division.
    logger
        Logger instance.

    Returns
    -------
    pandas.DataFrame
        Transcript metadata plus median isoform usage by tissue.
    """
    usage_medians: Dict[str, np.ndarray] = {}
    for tissue, indices in tissue_to_indices.items():
        logger.info(
            "Computing isoform-usage median for tissue %s using %d samples",
            tissue,
            len(indices),
        )
        expression_subset = expression_array[:, indices]
        gene_subset = gene_expression_array[transcript_gene_inverse][:, indices]
        usage_subset = np.divide(
            expression_subset,
            gene_subset + epsilon,
            out=np.zeros_like(expression_subset, dtype=expression_subset.dtype),
            where=gene_subset > 0,
        )
        usage_medians[tissue] = np.median(usage_subset, axis=1)

    usage_frame = pd.DataFrame(usage_medians)
    return pd.concat(
        objs=[row_metadata.reset_index(drop=True), usage_frame.reset_index(drop=True)],
        axis=1,
    )


def summarise_gene_target_tissue(
    *,
    gene_tissue_median: pd.DataFrame,
    target_tissue: str,
    epsilon: float,
) -> pd.DataFrame:
    """
    Summarise gene-level target-tissue specificity from tissue medians.

    Parameters
    ----------
    gene_tissue_median
        Gene metadata plus tissue median TPM values.
    target_tissue
        Target tissue name.
    epsilon
        Small constant for safe ratios.

    Returns
    -------
    pandas.DataFrame
        Gene-level target-tissue summary.
    """
    metadata_cols = ["gene_id", "gene_symbol"]
    tissue_cols = [
        column for column in gene_tissue_median.columns if column not in metadata_cols
    ]
    if target_tissue not in tissue_cols:
        raise ValueError(
            f"Target tissue {target_tissue!r} was not found in tissue columns."
        )
    non_target_cols = [column for column in tissue_cols if column != target_tissue]

    values = gene_tissue_median[tissue_cols]
    target = values[target_tissue].astype(float)
    max_non_target = values[non_target_cols].max(axis=1) if non_target_cols else 0.0
    median_non_target = values[non_target_cols].median(axis=1) if non_target_cols else 0.0
    max_tissue = values.idxmax(axis=1)

    summary = gene_tissue_median[metadata_cols].copy()
    summary["target_tissue"] = target_tissue
    summary["target_gene_median_tpm"] = target.values
    summary["max_non_target_gene_median_tpm"] = np.asarray(max_non_target)
    summary["median_non_target_gene_median_tpm"] = np.asarray(median_non_target)
    summary["log2_target_vs_max_non_target_gene_tpm"] = np.log2(
        (target + epsilon) / (np.asarray(max_non_target) + epsilon)
    )
    summary["gene_target_is_max_tissue"] = (max_tissue == target_tissue).astype(int)
    summary["gene_max_tpm_tissue"] = max_tissue.values
    return summary


def calculate_target_detection_fraction(
    *,
    expression_array: np.ndarray,
    target_indices: np.ndarray,
    min_tpm_present: float,
) -> np.ndarray:
    """
    Calculate the fraction of target-tissue samples with TPM above threshold.

    Parameters
    ----------
    expression_array
        Transcript by sample TPM array.
    target_indices
        Column indices for target-tissue samples.
    min_tpm_present
        TPM threshold for detection.

    Returns
    -------
    numpy.ndarray
        Fraction detected for each transcript row.
    """
    if len(target_indices) == 0:
        return np.zeros(expression_array.shape[0], dtype=float)
    return (expression_array[:, target_indices] >= min_tpm_present).mean(axis=1)


def summarise_transcript_target_tissue(
    *,
    transcript_tissue_median: pd.DataFrame,
    usage_tissue_median: pd.DataFrame,
    gene_summary: pd.DataFrame,
    expression_array: np.ndarray,
    tissue_to_indices: Dict[str, np.ndarray],
    target_tissue: str,
    min_tpm_present: float,
    min_target_tpm_candidate: float,
    min_target_usage_candidate: float,
    min_log2_usage_ratio_candidate: float,
    max_gene_log2_ratio_for_rescue: float,
    epsilon: float,
    logger: logging.Logger,
) -> pd.DataFrame:
    """
    Create a transcriptome-wide target-tissue isoform summary.

    Parameters
    ----------
    transcript_tissue_median
        Transcript metadata plus median TPM by tissue.
    usage_tissue_median
        Transcript metadata plus median isoform usage by tissue.
    gene_summary
        Gene-level target-tissue summary.
    expression_array
        Transcript by sample expression array.
    tissue_to_indices
        Tissue to sample-index mapping.
    target_tissue
        Target tissue name.
    min_tpm_present
        TPM threshold for target-tissue detection fraction.
    min_target_tpm_candidate
        Candidate threshold for target-tissue transcript TPM.
    min_target_usage_candidate
        Candidate threshold for target-tissue isoform usage.
    min_log2_usage_ratio_candidate
        Candidate threshold for target-tissue usage ratio.
    max_gene_log2_ratio_for_rescue
        Gene-level ratio ceiling for broad-gene rescue candidates.
    epsilon
        Small constant for safe ratios.
    logger
        Logger instance.

    Returns
    -------
    pandas.DataFrame
        Searchable target-tissue transcript summary.
    """
    metadata_cols = [
        "transcript_id_with_version",
        "transcript_id",
        "gene_id_with_version",
        "gene_id",
        "gene_symbol",
    ]
    transcript_tissue_cols = [
        column for column in transcript_tissue_median.columns if column not in metadata_cols
    ]
    usage_tissue_cols = [
        column for column in usage_tissue_median.columns if column not in metadata_cols
    ]
    if target_tissue not in transcript_tissue_cols:
        raise ValueError(f"Target tissue {target_tissue!r} was not found.")

    non_target_tpm_cols = [
        column for column in transcript_tissue_cols if column != target_tissue
    ]
    non_target_usage_cols = [
        column for column in usage_tissue_cols if column != target_tissue
    ]

    tpm_values = transcript_tissue_median[transcript_tissue_cols]
    usage_values = usage_tissue_median[usage_tissue_cols]

    target_tpm = tpm_values[target_tissue].astype(float)
    max_non_target_tpm = (
        tpm_values[non_target_tpm_cols].max(axis=1)
        if non_target_tpm_cols
        else pd.Series(np.zeros(tpm_values.shape[0]))
    )
    median_non_target_tpm = (
        tpm_values[non_target_tpm_cols].median(axis=1)
        if non_target_tpm_cols
        else pd.Series(np.zeros(tpm_values.shape[0]))
    )
    target_usage = usage_values[target_tissue].astype(float)
    max_non_target_usage = (
        usage_values[non_target_usage_cols].max(axis=1)
        if non_target_usage_cols
        else pd.Series(np.zeros(usage_values.shape[0]))
    )
    median_non_target_usage = (
        usage_values[non_target_usage_cols].median(axis=1)
        if non_target_usage_cols
        else pd.Series(np.zeros(usage_values.shape[0]))
    )

    summary = transcript_tissue_median[metadata_cols].copy()
    summary["target_tissue"] = target_tissue
    summary["target_median_tpm"] = target_tpm.values
    summary["max_non_target_median_tpm"] = np.asarray(max_non_target_tpm)
    summary["median_non_target_median_tpm"] = np.asarray(median_non_target_tpm)
    summary["log2_target_vs_max_non_target_tpm"] = np.log2(
        (target_tpm + epsilon) / (np.asarray(max_non_target_tpm) + epsilon)
    )
    summary["target_median_isoform_usage"] = target_usage.values
    summary["max_non_target_isoform_usage"] = np.asarray(max_non_target_usage)
    summary["median_non_target_isoform_usage"] = np.asarray(median_non_target_usage)
    summary["log2_target_vs_max_non_target_isoform_usage"] = np.log2(
        (target_usage + epsilon) / (np.asarray(max_non_target_usage) + epsilon)
    )
    summary["target_is_max_tpm_tissue"] = (
        tpm_values.idxmax(axis=1) == target_tissue
    ).astype(int)
    summary["target_is_max_usage_tissue"] = (
        usage_values.idxmax(axis=1) == target_tissue
    ).astype(int)

    target_indices = tissue_to_indices.get(target_tissue)
    if target_indices is None:
        raise ValueError(f"No samples were mapped to target tissue {target_tissue!r}.")
    summary["target_fraction_samples_tpm_ge_threshold"] = (
        calculate_target_detection_fraction(
            expression_array=expression_array,
            target_indices=target_indices,
            min_tpm_present=min_tpm_present,
        )
    )
    summary["target_present_tpm_threshold"] = min_tpm_present

    summary["target_tpm_rank_within_gene"] = summary.groupby(
        "gene_id", sort=False
    )["target_median_tpm"].rank(ascending=False, method="min")
    summary["target_usage_rank_within_gene"] = summary.groupby(
        "gene_id", sort=False
    )["target_median_isoform_usage"].rank(ascending=False, method="min")
    summary["n_transcripts_for_gene"] = summary.groupby("gene_id", sort=False)[
        "transcript_id"
    ].transform("count")

    summary = summary.merge(
        gene_summary,
        on=["gene_id", "gene_symbol", "target_tissue"],
        how="left",
        validate="many_to_one",
    )

    summary["is_target_tissue_isoform_candidate"] = (
        (summary["target_median_tpm"] >= min_target_tpm_candidate)
        & (summary["target_median_isoform_usage"] >= min_target_usage_candidate)
        & (
            summary["log2_target_vs_max_non_target_isoform_usage"]
            >= min_log2_usage_ratio_candidate
        )
        & (summary["target_usage_rank_within_gene"] == 1)
        & (summary["target_is_max_usage_tissue"] == 1)
    ).astype(int)

    summary["is_broad_gene_isoform_rescue_candidate"] = (
        (summary["is_target_tissue_isoform_candidate"] == 1)
        & (
            summary["log2_target_vs_max_non_target_gene_tpm"]
            <= max_gene_log2_ratio_for_rescue
        )
    ).astype(int)

    sort_cols = [
        "is_broad_gene_isoform_rescue_candidate",
        "is_target_tissue_isoform_candidate",
        "log2_target_vs_max_non_target_isoform_usage",
        "target_median_isoform_usage",
        "target_median_tpm",
    ]
    summary = summary.sort_values(by=sort_cols, ascending=[False] * len(sort_cols))
    logger.info(
        "Transcript summary contains %d target-tissue candidates and %d broad-gene rescue candidates",
        int(summary["is_target_tissue_isoform_candidate"].sum()),
        int(summary["is_broad_gene_isoform_rescue_candidate"].sum()),
    )
    return summary


def format_top_tissues(*, values: pd.Series, top_n: int) -> str:
    """
    Format top tissues and values as a compact semicolon-separated string.

    Parameters
    ----------
    values
        Tissue-value series.
    top_n
        Number of top tissues to include.

    Returns
    -------
    str
        Semicolon-separated tissue:value string.
    """
    sorted_values = values.sort_values(ascending=False).head(top_n)
    return ";".join(
        f"{tissue}:{value:.4g}" for tissue, value in sorted_values.items()
    )


def add_candidate_top_tissue_strings(
    *,
    candidates: pd.DataFrame,
    transcript_tissue_median: pd.DataFrame,
    usage_tissue_median: pd.DataFrame,
    top_n: int,
) -> pd.DataFrame:
    """
    Add top TPM and isoform-usage tissue strings to candidate rows.

    Parameters
    ----------
    candidates
        Candidate transcript summary rows.
    transcript_tissue_median
        Transcript median TPM matrix.
    usage_tissue_median
        Isoform-usage median matrix.
    top_n
        Number of tissues to report.

    Returns
    -------
    pandas.DataFrame
        Candidate rows with top tissue strings.
    """
    if candidates.empty:
        output = candidates.copy()
        output["top_tissues_by_median_tpm"] = ""
        output["top_tissues_by_median_isoform_usage"] = ""
        return output

    metadata_cols = [
        "transcript_id_with_version",
        "transcript_id",
        "gene_id_with_version",
        "gene_id",
        "gene_symbol",
    ]
    tpm_tissue_cols = [
        column for column in transcript_tissue_median.columns if column not in metadata_cols
    ]
    usage_tissue_cols = [
        column for column in usage_tissue_median.columns if column not in metadata_cols
    ]

    transcript_tpm_index = transcript_tissue_median.set_index("transcript_id")
    usage_index = usage_tissue_median.set_index("transcript_id")

    output = candidates.copy()
    output["top_tissues_by_median_tpm"] = [
        format_top_tissues(
            values=transcript_tpm_index.loc[transcript_id, tpm_tissue_cols],
            top_n=top_n,
        )
        for transcript_id in output["transcript_id"]
    ]
    output["top_tissues_by_median_isoform_usage"] = [
        format_top_tissues(
            values=usage_index.loc[transcript_id, usage_tissue_cols],
            top_n=top_n,
        )
        for transcript_id in output["transcript_id"]
    ]
    return output


def select_candidate_outputs(
    *,
    summary: pd.DataFrame,
    transcript_tissue_median: pd.DataFrame,
    usage_tissue_median: pd.DataFrame,
    top_n: int,
    logger: logging.Logger,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Select all candidate isoforms and the best candidate per gene.

    Parameters
    ----------
    summary
        Full transcript summary.
    transcript_tissue_median
        Transcript median TPM matrix.
    usage_tissue_median
        Isoform-usage median matrix.
    top_n
        Number of top tissues to report in candidate files.
    logger
        Logger instance.

    Returns
    -------
    tuple[pandas.DataFrame, pandas.DataFrame]
        Candidate isoforms and best candidate isoform per gene.
    """
    candidates = summary.loc[summary["is_target_tissue_isoform_candidate"] == 1].copy()
    candidates = add_candidate_top_tissue_strings(
        candidates=candidates,
        transcript_tissue_median=transcript_tissue_median,
        usage_tissue_median=usage_tissue_median,
        top_n=top_n,
    )

    if candidates.empty:
        best_per_gene = candidates.copy()
    else:
        best_per_gene = candidates.sort_values(
            by=[
                "is_broad_gene_isoform_rescue_candidate",
                "log2_target_vs_max_non_target_isoform_usage",
                "target_median_isoform_usage",
                "target_median_tpm",
            ],
            ascending=[False, False, False, False],
        ).drop_duplicates(subset=["gene_id"], keep="first")

    logger.info(
        "Selected %d candidate transcript rows across %d genes",
        candidates.shape[0],
        candidates["gene_id"].nunique() if not candidates.empty else 0,
    )
    logger.info("Selected %d best-candidate genes", best_per_gene.shape[0])
    return candidates, best_per_gene


def write_tsv(*, dataframe: pd.DataFrame, path: Path, logger: logging.Logger) -> None:
    """
    Write a DataFrame to a TSV or TSV.GZ file.

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
    logger.info("Wrote %d rows x %d columns to %s", dataframe.shape[0], dataframe.shape[1], path)


def run(*, config: Config, logger: logging.Logger) -> None:
    """
    Run the full transcriptome-wide isoform screen.

    Parameters
    ----------
    config
        Workflow configuration.
    logger
        Logger instance.
    """
    config.out_dir.mkdir(parents=True, exist_ok=True)

    with timed(logger=logger, label="Load sample annotations"):
        sample_attributes = load_sample_attributes(
            path=config.sample_attributes_path,
            sample_id_col=config.sample_id_col,
            tissue_col=config.tissue_col,
            logger=logger,
        )
        sample_ids = tuple(sample_attributes[config.sample_id_col].astype(str))

    with timed(logger=logger, label="Load gene ID to symbol map"):
        gene_id_to_symbol = load_gene_symbol_map(
            gct_path=config.gene_tpm_gct_path,
            logger=logger,
        )

    with timed(logger=logger, label="Detect transcript matrix schema"):
        header_line = find_matrix_header_line(
            matrix_path=config.transcript_tpm_path,
            logger=logger,
        )
        columns = preview_transcript_header(
            path=config.transcript_tpm_path,
            header_line=header_line,
            logger=logger,
        )
        schema = detect_transcript_schema(columns=columns, sample_ids=sample_ids)
        logger.info("Transcript ID column: %s", schema.transcript_col)
        logger.info("Gene ID column: %s", schema.gene_col)
        logger.info("Detected %d metadata columns", len(schema.metadata_cols))
        logger.info("Detected %d annotated sample columns", len(schema.sample_cols))
        if schema.ignored_expression_like_cols:
            logger.warning(
                "Ignoring %d GTEX-like columns absent from sample attributes",
                len(schema.ignored_expression_like_cols),
            )

    with timed(logger=logger, label="Load full transcript TPM matrix"):
        raw_transcript_df = load_full_transcript_tpm(
            path=config.transcript_tpm_path,
            schema=schema,
            header_line=header_line,
            float_dtype=config.float_dtype,
            logger=logger,
        )
        transcript_df = normalise_transcript_dataframe(
            dataframe=raw_transcript_df,
            schema=schema,
            gene_id_to_symbol=gene_id_to_symbol,
            sample_cols=schema.sample_cols,
            logger=logger,
        )
        del raw_transcript_df

    with timed(logger=logger, label="Build expression and tissue index arrays"):
        expression_array = make_expression_array(
            dataframe=transcript_df,
            sample_cols=schema.sample_cols,
            float_dtype=config.float_dtype,
            logger=logger,
        )
        tissue_to_indices = make_tissue_sample_indices(
            sample_cols=schema.sample_cols,
            sample_attributes=sample_attributes,
            sample_id_col=config.sample_id_col,
            tissue_col=config.tissue_col,
            logger=logger,
        )
        if config.target_tissue not in tissue_to_indices:
            available = ", ".join(sorted(tissue_to_indices)[:20])
            raise ValueError(
                f"Target tissue {config.target_tissue!r} was not found. "
                f"First available tissues are: {available}"
            )

    metadata_cols = [
        "transcript_id_with_version",
        "transcript_id",
        "gene_id_with_version",
        "gene_id",
        "gene_symbol",
    ]
    row_metadata = transcript_df[metadata_cols].copy()

    with timed(logger=logger, label="Compute transcript tissue median TPM"):
        transcript_tissue_median = compute_tissue_medians(
            values=expression_array,
            row_metadata=row_metadata,
            tissue_to_indices=tissue_to_indices,
            logger=logger,
            value_label="transcript TPM",
        )

    with timed(logger=logger, label="Compute gene expression and tissue median TPM"):
        unique_gene_ids, transcript_gene_inverse, gene_expression_array = (
            compute_gene_expression_array(
                expression_array=expression_array,
                gene_ids=row_metadata["gene_id"],
                float_dtype=config.float_dtype,
                logger=logger,
            )
        )
        gene_metadata = make_gene_metadata(
            gene_ids=unique_gene_ids,
            gene_id_to_symbol=gene_id_to_symbol,
        )
        gene_tissue_median = compute_tissue_medians(
            values=gene_expression_array,
            row_metadata=gene_metadata,
            tissue_to_indices=tissue_to_indices,
            logger=logger,
            value_label="gene TPM",
        )
        gene_summary = summarise_gene_target_tissue(
            gene_tissue_median=gene_tissue_median,
            target_tissue=config.target_tissue,
            epsilon=config.epsilon,
        )

    with timed(logger=logger, label="Compute isoform usage tissue medians"):
        usage_tissue_median = compute_isoform_usage_tissue_medians(
            expression_array=expression_array,
            gene_expression_array=gene_expression_array,
            transcript_gene_inverse=transcript_gene_inverse,
            row_metadata=row_metadata,
            tissue_to_indices=tissue_to_indices,
            epsilon=config.epsilon,
            logger=logger,
        )

    with timed(logger=logger, label="Build target-tissue transcript summary"):
        transcript_summary = summarise_transcript_target_tissue(
            transcript_tissue_median=transcript_tissue_median,
            usage_tissue_median=usage_tissue_median,
            gene_summary=gene_summary,
            expression_array=expression_array,
            tissue_to_indices=tissue_to_indices,
            target_tissue=config.target_tissue,
            min_tpm_present=config.min_tpm_present,
            min_target_tpm_candidate=config.min_target_tpm_candidate,
            min_target_usage_candidate=config.min_target_usage_candidate,
            min_log2_usage_ratio_candidate=config.min_log2_usage_ratio_candidate,
            max_gene_log2_ratio_for_rescue=config.max_gene_log2_ratio_for_rescue,
            epsilon=config.epsilon,
            logger=logger,
        )
        candidates, best_per_gene = select_candidate_outputs(
            summary=transcript_summary,
            transcript_tissue_median=transcript_tissue_median,
            usage_tissue_median=usage_tissue_median,
            top_n=config.top_n_candidate_tissues,
            logger=logger,
        )

    with timed(logger=logger, label="Write outputs"):
        prefix = config.out_dir / config.out_prefix
        if config.write_tissue_matrices:
            write_tsv(
                dataframe=transcript_tissue_median,
                path=Path(f"{prefix}.transcript_tissue_median_tpm.tsv.gz"),
                logger=logger,
            )
            write_tsv(
                dataframe=usage_tissue_median,
                path=Path(
                    f"{prefix}.transcript_isoform_usage_median_by_tissue.tsv.gz"
                ),
                logger=logger,
            )
            write_tsv(
                dataframe=gene_tissue_median,
                path=Path(f"{prefix}.gene_tissue_median_tpm.tsv.gz"),
                logger=logger,
            )
        write_tsv(
            dataframe=gene_summary,
            path=Path(f"{prefix}.gene_target_tissue_summary.tsv.gz"),
            logger=logger,
        )
        write_tsv(
            dataframe=transcript_summary,
            path=Path(f"{prefix}.transcript_target_tissue_isoform_summary.tsv.gz"),
            logger=logger,
        )
        write_tsv(
            dataframe=candidates,
            path=Path(f"{prefix}.candidate_target_tissue_isoforms.tsv"),
            logger=logger,
        )
        write_tsv(
            dataframe=best_per_gene,
            path=Path(f"{prefix}.best_candidate_isoform_per_gene.tsv"),
            logger=logger,
        )

    logger.info("GTEx transcriptome-wide isoform screen finished")


def main() -> None:
    """Run the command-line entry point."""
    config = parse_args()
    logger = setup_logging(log_level=config.log_level, log_path=config.log_path)
    logger.info("GTEx transcriptome-wide isoform screen starting")
    logger.info("Transcript TPM path: %s", config.transcript_tpm_path)
    logger.info("Gene TPM GCT path: %s", config.gene_tpm_gct_path)
    logger.info("Sample attributes path: %s", config.sample_attributes_path)
    logger.info("Output directory: %s", config.out_dir)
    logger.info("Target tissue: %s", config.target_tissue)
    logger.info("Float dtype: %s", config.float_dtype)
    run(config=config, logger=logger)


if __name__ == "__main__":
    main()
