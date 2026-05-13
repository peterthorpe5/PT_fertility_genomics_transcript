#!/usr/bin/env python3
"""
Annotate GTEx transcript isoform screen outputs with GENCODE metadata.

This script parses a GENCODE GTF file, creates transcript-level and
feature-level annotation tables, and joins transcript annotation onto the
GTEx transcriptome-wide isoform-screen outputs.

The intended use is the transcript-level fertility genomics workflow, where
GTEx v11 transcript TPM data have been screened for testis-preferential
isoform usage. For GTEx v11, use the matching GENCODE v47 primary assembly
comprehensive GTF where possible.

Primary outputs
---------------
<out_prefix>.gencode_transcript_annotation.tsv.gz
    Transcript-level annotation derived from the GTF.

<out_prefix>.gencode_transcript_features.tsv.gz
    Exon and CDS feature coordinates for downstream gene-model plotting.

<out_prefix>.transcript_target_tissue_isoform_summary.annotated.tsv.gz
    Annotated full transcript isoform screen summary, if supplied.

<out_prefix>.candidate_target_tissue_isoforms.annotated.tsv
    Annotated candidate isoform table, if supplied.

<out_prefix>.best_candidate_isoform_per_gene.annotated.tsv
    Annotated best-per-gene candidate table, if supplied.
"""

from __future__ import annotations

import argparse
import gzip
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


ATTRIBUTE_PATTERN = re.compile(r'([A-Za-z0-9_]+)\s+"([^"]*)"')


@dataclass(frozen=True)
class Config:
    """
    Configuration for GENCODE annotation of GTEx isoform-screen outputs.

    Attributes
    ----------
    gtf_path
        Path to the GENCODE GTF or GTF.GZ file.
    out_dir
        Output directory.
    out_prefix
        Prefix used for output filenames.
    screen_tsv
        Optional full transcript isoform-screen summary table.
    candidate_tsv
        Optional candidate isoform table.
    best_per_gene_tsv
        Optional best candidate per gene table.
    write_excel_outputs
        Whether to write formatted Excel copies of annotated browseable tables.
    excel_max_rows
        Maximum number of rows allowed for Excel output.
    log_path
        Optional log-file path.
    log_level
        Logging level.
    """

    gtf_path: Path
    out_dir: Path
    out_prefix: str
    screen_tsv: Optional[Path]
    candidate_tsv: Optional[Path]
    best_per_gene_tsv: Optional[Path]
    write_excel_outputs: bool
    excel_max_rows: int
    log_path: Optional[Path]
    log_level: str


@contextmanager
def timed(*, logger: logging.Logger, label: str) -> Iterator[None]:
    """
    Log elapsed time for a workflow block.

    Parameters
    ----------
    logger
        Logger instance.
    label
        Human-readable block label.

    Yields
    ------
    Iterator[None]
        Context manager body.
    """
    start_time = time.time()
    logger.info("%s: started", label)
    try:
        yield
    except Exception:
        elapsed = time.time() - start_time
        logger.exception("%s: failed after %.1f seconds", label, elapsed)
        raise
    elapsed = time.time() - start_time
    logger.info("%s: completed in %.1f seconds", label, elapsed)


def setup_logging(*, log_level: str, log_path: Optional[Path]) -> logging.Logger:
    """
    Configure console and optional file logging.

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
    logger = logging.getLogger("annotate_gtex_isoform_screen_with_gencode")
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
            "Annotate GTEx transcript isoform-screen outputs with GENCODE "
            "transcript, exon and CDS metadata."
        )
    )
    parser.add_argument(
        "--gtf",
        required=True,
        type=Path,
        help="GENCODE GTF or GTF.GZ file, ideally release-matched to GTEx.",
    )
    parser.add_argument(
        "--out_dir",
        required=True,
        type=Path,
        help="Output directory.",
    )
    parser.add_argument(
        "--out_prefix",
        default="gtex_v11_transcriptome_testis_isoform_screen",
        help="Output filename prefix.",
    )
    parser.add_argument(
        "--screen_tsv",
        type=Path,
        default=None,
        help="Optional full transcript isoform-screen summary TSV/TSV.GZ.",
    )
    parser.add_argument(
        "--candidate_tsv",
        type=Path,
        default=None,
        help="Optional candidate isoform TSV/TSV.GZ.",
    )
    parser.add_argument(
        "--best_per_gene_tsv",
        type=Path,
        default=None,
        help="Optional best candidate per gene TSV/TSV.GZ.",
    )
    parser.add_argument(
        "--write_excel_outputs",
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "Write formatted Excel copies of annotated candidate and "
            "best-per-gene tables. Default: True. Use "
            "--no-write_excel_outputs to disable."
        ),
    )
    parser.add_argument(
        "--excel_max_rows",
        type=int,
        default=200000,
        help=(
            "Maximum number of rows allowed for an Excel output. "
            "Default: 200000."
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
        gtf_path=args.gtf,
        out_dir=args.out_dir,
        out_prefix=args.out_prefix,
        screen_tsv=args.screen_tsv,
        candidate_tsv=args.candidate_tsv,
        best_per_gene_tsv=args.best_per_gene_tsv,
        write_excel_outputs=args.write_excel_outputs,
        excel_max_rows=args.excel_max_rows,
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
    Remove an Ensembl dot-version suffix from an identifier.

    Parameters
    ----------
    identifier
        Identifier value.

    Returns
    -------
    str
        Identifier without trailing dot-version suffix.
    """
    value = str(identifier).strip()
    if not value or value.lower() == "nan":
        return ""
    return value.split(".", maxsplit=1)[0]


def parse_gtf_attributes(*, attribute_text: str) -> Dict[str, List[str]]:
    """
    Parse the attribute column from a GTF record.

    Parameters
    ----------
    attribute_text
        Raw GTF attribute text.

    Returns
    -------
    Dict[str, List[str]]
        Attribute mapping. Repeated keys, such as tag, are stored as lists.
    """
    attributes: Dict[str, List[str]] = {}
    for key, value in ATTRIBUTE_PATTERN.findall(attribute_text):
        attributes.setdefault(key, []).append(value)
    return attributes


def get_first_attribute(
    *,
    attributes: Dict[str, List[str]],
    key: str,
    default: str = "",
) -> str:
    """
    Return the first value for an attribute key.

    Parameters
    ----------
    attributes
        Parsed GTF attributes.
    key
        Attribute key.
    default
        Value returned if the key is absent.

    Returns
    -------
    str
        First attribute value or default.
    """
    values = attributes.get(key)
    if not values:
        return default
    return values[0]


def join_unique_values(*, values: Iterable[str]) -> str:
    """
    Join unique non-empty values as a semicolon-separated string.

    Parameters
    ----------
    values
        Values to join.

    Returns
    -------
    str
        Semicolon-separated unique values in first-seen order.
    """
    seen = set()
    ordered = []
    for value in values:
        if value is None:
            continue
        clean_value = str(value).strip()
        if not clean_value or clean_value in seen:
            continue
        seen.add(clean_value)
        ordered.append(clean_value)
    return ";".join(ordered)


def update_first_non_empty(*, current: str, new_value: str) -> str:
    """
    Keep the current value unless it is empty.

    Parameters
    ----------
    current
        Existing value.
    new_value
        Candidate replacement value.

    Returns
    -------
    str
        Existing non-empty value or the new value.
    """
    if str(current).strip():
        return current
    return str(new_value).strip()


def make_transcript_record(
    *,
    seqname: str,
    source: str,
    start: int,
    end: int,
    strand: str,
    attributes: Dict[str, List[str]],
) -> Dict[str, object]:
    """
    Create a transcript-level metadata record from parsed GTF fields.

    Parameters
    ----------
    seqname
        Sequence or chromosome name.
    source
        GTF source field.
    start
        Transcript start coordinate.
    end
        Transcript end coordinate.
    strand
        Transcript strand.
    attributes
        Parsed GTF attributes.

    Returns
    -------
    Dict[str, object]
        Transcript metadata record.
    """
    gene_id_with_version = get_first_attribute(
        attributes=attributes,
        key="gene_id",
    )
    transcript_id_with_version = get_first_attribute(
        attributes=attributes,
        key="transcript_id",
    )
    tags = attributes.get("tag", [])

    return {
        "gencode_gene_id_with_version": gene_id_with_version,
        "gencode_gene_id": strip_ensembl_version(
            identifier=gene_id_with_version,
        ),
        "gencode_gene_name": get_first_attribute(
            attributes=attributes,
            key="gene_name",
        ),
        "gencode_gene_type": get_first_attribute(
            attributes=attributes,
            key="gene_type",
        ),
        "gencode_transcript_id_with_version": transcript_id_with_version,
        "gencode_transcript_id": strip_ensembl_version(
            identifier=transcript_id_with_version,
        ),
        "gencode_transcript_name": get_first_attribute(
            attributes=attributes,
            key="transcript_name",
        ),
        "gencode_transcript_type": get_first_attribute(
            attributes=attributes,
            key="transcript_type",
        ),
        "gencode_seqname": seqname,
        "gencode_source": source,
        "gencode_strand": strand,
        "gencode_transcript_start": start,
        "gencode_transcript_end": end,
        "gencode_transcript_support_level": get_first_attribute(
            attributes=attributes,
            key="transcript_support_level",
        ),
        "gencode_havana_transcript": get_first_attribute(
            attributes=attributes,
            key="havana_transcript",
        ),
        "gencode_ccdsid": get_first_attribute(
            attributes=attributes,
            key="ccdsid",
        ),
        "gencode_tags": join_unique_values(values=tags),
        "is_basic": int("basic" in tags),
        "is_mane_select": int("MANE_Select" in tags),
        "is_mane_plus_clinical": int("MANE_Plus_Clinical" in tags),
        "is_ensembl_canonical": int("Ensembl_canonical" in tags),
    }


def parse_gencode_gtf(
    *,
    gtf_path: Path,
    logger: logging.Logger,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Parse transcript, exon and CDS metadata from a GENCODE GTF file.

    Parameters
    ----------
    gtf_path
        GTF or GTF.GZ path.
    logger
        Logger instance.

    Returns
    -------
    Tuple[pandas.DataFrame, pandas.DataFrame]
        Transcript annotation table and feature-coordinate table.
    """
    transcript_records: Dict[str, Dict[str, object]] = {}
    feature_records: List[Dict[str, object]] = []
    exon_intervals: Dict[str, List[Tuple[int, int]]] = {}
    cds_intervals: Dict[str, List[Tuple[int, int]]] = {}
    protein_ids: Dict[str, List[str]] = {}
    feature_counts = {"gene": 0, "transcript": 0, "exon": 0, "CDS": 0}

    with open_text_maybe_gzip(path=gtf_path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                logger.warning(
                    "Skipping malformed GTF line %d with %d fields",
                    line_number,
                    len(fields),
                )
                continue

            seqname, source, feature_type, start_text, end_text = fields[:5]
            strand = fields[6]
            attributes = parse_gtf_attributes(attribute_text=fields[8])
            transcript_id_with_version = get_first_attribute(
                attributes=attributes,
                key="transcript_id",
            )
            gene_id_with_version = get_first_attribute(
                attributes=attributes,
                key="gene_id",
            )

            if feature_type in feature_counts:
                feature_counts[feature_type] += 1

            if feature_type not in {"transcript", "exon", "CDS"}:
                continue

            if not transcript_id_with_version:
                continue

            try:
                start = int(start_text)
                end = int(end_text)
            except ValueError:
                logger.warning(
                    "Skipping GTF line %d with non-integer coordinates",
                    line_number,
                )
                continue

            if transcript_id_with_version not in transcript_records:
                transcript_records[transcript_id_with_version] = make_transcript_record(
                    seqname=seqname,
                    source=source,
                    start=start,
                    end=end,
                    strand=strand,
                    attributes=attributes,
                )

            record = transcript_records[transcript_id_with_version]
            record["gencode_transcript_start"] = min(
                int(record["gencode_transcript_start"]),
                start,
            )
            record["gencode_transcript_end"] = max(
                int(record["gencode_transcript_end"]),
                end,
            )
            record["gencode_gene_id_with_version"] = update_first_non_empty(
                current=str(record.get("gencode_gene_id_with_version", "")),
                new_value=gene_id_with_version,
            )
            record["gencode_gene_id"] = update_first_non_empty(
                current=str(record.get("gencode_gene_id", "")),
                new_value=strip_ensembl_version(identifier=gene_id_with_version),
            )
            record["gencode_gene_name"] = update_first_non_empty(
                current=str(record.get("gencode_gene_name", "")),
                new_value=get_first_attribute(attributes=attributes, key="gene_name"),
            )
            record["gencode_gene_type"] = update_first_non_empty(
                current=str(record.get("gencode_gene_type", "")),
                new_value=get_first_attribute(attributes=attributes, key="gene_type"),
            )
            record["gencode_transcript_name"] = update_first_non_empty(
                current=str(record.get("gencode_transcript_name", "")),
                new_value=get_first_attribute(
                    attributes=attributes,
                    key="transcript_name",
                ),
            )
            record["gencode_transcript_type"] = update_first_non_empty(
                current=str(record.get("gencode_transcript_type", "")),
                new_value=get_first_attribute(
                    attributes=attributes,
                    key="transcript_type",
                ),
            )

            if feature_type in {"exon", "CDS"}:
                exon_number = get_first_attribute(
                    attributes=attributes,
                    key="exon_number",
                )
                feature_records.append(
                    {
                        "gencode_gene_id_with_version": gene_id_with_version,
                        "gencode_gene_id": strip_ensembl_version(
                            identifier=gene_id_with_version,
                        ),
                        "gencode_gene_name": get_first_attribute(
                            attributes=attributes,
                            key="gene_name",
                        ),
                        "gencode_transcript_id_with_version": (
                            transcript_id_with_version
                        ),
                        "gencode_transcript_id": strip_ensembl_version(
                            identifier=transcript_id_with_version,
                        ),
                        "gencode_transcript_name": get_first_attribute(
                            attributes=attributes,
                            key="transcript_name",
                        ),
                        "gencode_transcript_type": get_first_attribute(
                            attributes=attributes,
                            key="transcript_type",
                        ),
                        "gencode_seqname": seqname,
                        "gencode_source": source,
                        "gencode_feature_type": feature_type,
                        "gencode_start": start,
                        "gencode_end": end,
                        "gencode_strand": strand,
                        "gencode_exon_number": exon_number,
                        "gencode_exon_id": get_first_attribute(
                            attributes=attributes,
                            key="exon_id",
                        ),
                        "gencode_protein_id": get_first_attribute(
                            attributes=attributes,
                            key="protein_id",
                        ),
                    }
                )

            if feature_type == "exon":
                exon_intervals.setdefault(transcript_id_with_version, []).append(
                    (start, end)
                )
            elif feature_type == "CDS":
                cds_intervals.setdefault(transcript_id_with_version, []).append(
                    (start, end)
                )
                protein_id = get_first_attribute(
                    attributes=attributes,
                    key="protein_id",
                )
                if protein_id:
                    protein_ids.setdefault(transcript_id_with_version, []).append(
                        protein_id,
                    )

    transcript_table = pd.DataFrame.from_records(list(transcript_records.values()))
    if transcript_table.empty:
        raise ValueError(f"No transcript records were parsed from {gtf_path}")

    transcript_table["gencode_exon_count"] = transcript_table[
        "gencode_transcript_id_with_version"
    ].map(lambda value: len(exon_intervals.get(value, [])))
    transcript_table["gencode_transcript_length_bp"] = transcript_table[
        "gencode_transcript_id_with_version"
    ].map(
        lambda value: sum(
            end - start + 1
            for start, end in exon_intervals.get(value, [])
        )
    )
    transcript_table["gencode_cds_exon_count"] = transcript_table[
        "gencode_transcript_id_with_version"
    ].map(lambda value: len(cds_intervals.get(value, [])))
    transcript_table["gencode_cds_length_bp"] = transcript_table[
        "gencode_transcript_id_with_version"
    ].map(
        lambda value: sum(
            end - start + 1
            for start, end in cds_intervals.get(value, [])
        )
    )
    transcript_table["gencode_has_cds"] = (
        transcript_table["gencode_cds_length_bp"] > 0
    ).astype(int)
    transcript_table["gencode_protein_id"] = transcript_table[
        "gencode_transcript_id_with_version"
    ].map(lambda value: join_unique_values(values=protein_ids.get(value, [])))
    transcript_table["gencode_is_protein_coding_transcript"] = (
        transcript_table["gencode_transcript_type"].eq("protein_coding")
        & transcript_table["gencode_has_cds"].eq(1)
    ).astype(int)
    transcript_table["gencode_genomic_span_bp"] = (
        transcript_table["gencode_transcript_end"].astype(int)
        - transcript_table["gencode_transcript_start"].astype(int)
        + 1
    )

    feature_table = pd.DataFrame.from_records(feature_records)
    logger.info(
        "Parsed %d transcript records and %d exon/CDS feature records",
        transcript_table.shape[0],
        feature_table.shape[0],
    )
    logger.info("Observed GTF feature counts: %s", feature_counts)
    logger.info(
        "Protein-coding transcript annotations: %d",
        int(transcript_table["gencode_is_protein_coding_transcript"].sum()),
    )
    return transcript_table, feature_table


def read_tsv(*, path: Path, logger: logging.Logger) -> pd.DataFrame:
    """
    Read a TSV or TSV.GZ file.

    Parameters
    ----------
    path
        Input path.
    logger
        Logger instance.

    Returns
    -------
    pandas.DataFrame
        Loaded table.
    """
    dataframe = pd.read_csv(filepath_or_buffer=path, sep="\t", dtype=str)
    logger.info(
        "Read %d rows x %d columns from %s",
        dataframe.shape[0],
        dataframe.shape[1],
        path,
    )
    return dataframe


def write_tsv(*, dataframe: pd.DataFrame, path: Path, logger: logging.Logger) -> None:
    """
    Write a DataFrame as TSV or TSV.GZ.

    Parameters
    ----------
    dataframe
        Table to write.
    path
        Output path.
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
    invalid_characters = ["[", "]", ":", "*", "?", "/", "\\"]
    safe_name = sheet_name
    for invalid_character in invalid_characters:
        safe_name = safe_name.replace(invalid_character, "_")
    safe_name = safe_name.strip() or "results"
    return safe_name[:31]


def estimate_excel_column_widths(
    *,
    dataframe: pd.DataFrame,
    max_scan_rows: int = 2000,
    min_width: int = 8,
    max_text_width: int = 52,
    max_numeric_width: int = 18,
) -> Dict[str, int]:
    """
    Estimate readable Excel column widths.

    Parameters
    ----------
    dataframe
        Table to inspect.
    max_scan_rows
        Maximum number of rows to scan.
    min_width
        Minimum column width.
    max_text_width
        Maximum width for text columns.
    max_numeric_width
        Maximum width for numeric columns.

    Returns
    -------
    Dict[str, int]
        Column width mapping.
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
        numeric_series = pd.to_numeric(dataframe[column], errors="coerce")
        numeric_fraction = numeric_series.notna().mean() if len(dataframe) else 0
        if numeric_fraction > 0.8:
            width = min(max(raw_width, min_width), max_numeric_width)
        else:
            width = min(max(raw_width, min_width), max_text_width)
        widths[str(column)] = width
    return widths


def write_formatted_excel(
    *,
    dataframe: pd.DataFrame,
    path: Path,
    sheet_name: str,
    logger: logging.Logger,
    excel_max_rows: int,
) -> None:
    """
    Write a formatted Excel copy of a browseable result table.

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
    excel_max_rows
        Maximum allowed number of rows.
    """
    n_rows, n_cols = dataframe.shape
    if n_rows > excel_max_rows:
        logger.warning(
            "Skipping Excel output %s because %d rows exceeds limit %d",
            path,
            n_rows,
            excel_max_rows,
        )
        return
    if n_rows + 1 > 1_048_576 or n_cols > 16_384:
        logger.warning("Skipping Excel output %s because it exceeds Excel limits", path)
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
                properties={
                    "bold": True,
                    "text_wrap": True,
                    "valign": "top",
                    "border": 1,
                }
            )
            integer_format = workbook.add_format(properties={"num_format": "0"})
            float_format = workbook.add_format(properties={"num_format": "0.0000"})

            worksheet.freeze_panes(row=1, col=0)
            if n_cols > 0:
                worksheet.add_table(
                    first_row=0,
                    first_col=0,
                    last_row=n_rows,
                    last_col=n_cols - 1,
                    options={
                        "columns": [{"header": str(col)} for col in dataframe.columns],
                        "style": "Table Style Medium 2",
                        "autofilter": True,
                    },
                )
            worksheet.set_row(row=0, height=30, cell_format=header_format)
            widths = estimate_excel_column_widths(dataframe=dataframe)
            for column_index, column in enumerate(dataframe.columns):
                numeric_series = pd.to_numeric(dataframe[column], errors="coerce")
                numeric_fraction = (
                    numeric_series.notna().mean()
                    if len(dataframe)
                    else 0
                )
                if numeric_fraction > 0.8:
                    non_na = numeric_series.dropna()
                    if not non_na.empty and np.all(np.mod(non_na, 1) == 0):
                        cell_format = integer_format
                    else:
                        cell_format = float_format
                else:
                    cell_format = None
                worksheet.set_column(
                    first_col=column_index,
                    last_col=column_index,
                    width=widths[str(column)],
                    cell_format=cell_format,
                )
    except ImportError as error:
        logger.warning("Could not write Excel file %s: %s", path, error)
        return

    logger.info("Wrote formatted Excel file to %s", path)


def add_stripped_transcript_ids(*, dataframe: pd.DataFrame) -> pd.DataFrame:
    """
    Add helper transcript ID columns for robust annotation joins.

    Parameters
    ----------
    dataframe
        Input table.

    Returns
    -------
    pandas.DataFrame
        Table with helper join columns.
    """
    output = dataframe.copy()
    if "transcript_id_with_version" in output.columns:
        output["_join_transcript_id_with_version"] = output[
            "transcript_id_with_version"
        ].astype(str)
    else:
        output["_join_transcript_id_with_version"] = ""

    if "transcript_id" in output.columns:
        output["_join_transcript_id"] = output["transcript_id"].map(
            lambda value: strip_ensembl_version(identifier=value)
        )
    elif "transcript_id_with_version" in output.columns:
        output["_join_transcript_id"] = output["transcript_id_with_version"].map(
            lambda value: strip_ensembl_version(identifier=value)
        )
    else:
        raise ValueError(
            "Input table must contain transcript_id or transcript_id_with_version."
        )
    return output


def annotate_result_table(
    *,
    dataframe: pd.DataFrame,
    transcript_annotation: pd.DataFrame,
    logger: logging.Logger,
) -> pd.DataFrame:
    """
    Join GENCODE transcript metadata onto a GTEx result table.

    A versioned transcript ID join is attempted first. Remaining unmatched rows
    are filled using a version-stripped transcript ID join.

    Parameters
    ----------
    dataframe
        GTEx result table.
    transcript_annotation
        GENCODE transcript annotation table.
    logger
        Logger instance.

    Returns
    -------
    pandas.DataFrame
        Annotated result table.
    """
    left = add_stripped_transcript_ids(dataframe=dataframe)
    left["_row_id"] = np.arange(left.shape[0])

    annotation = transcript_annotation.drop_duplicates(
        subset=["gencode_transcript_id_with_version"],
        keep="first",
    ).copy()
    version_join = annotation.rename(
        columns={
            "gencode_transcript_id_with_version": "_join_transcript_id_with_version",
        }
    )
    version_merged = left.merge(
        version_join,
        on="_join_transcript_id_with_version",
        how="left",
        validate="many_to_one",
    )
    annotation_cols = [
        column
        for column in version_join.columns
        if column not in {"_join_transcript_id_with_version"}
    ]
    matched_by_version = version_merged["gencode_gene_id"].notna()

    fallback_annotation = transcript_annotation.drop_duplicates(
        subset=["gencode_transcript_id"],
        keep="first",
    ).rename(columns={"gencode_transcript_id": "_join_transcript_id"})
    fallback_merged = left.loc[~matched_by_version].merge(
        fallback_annotation,
        on="_join_transcript_id",
        how="left",
        validate="many_to_one",
    )

    output = version_merged.copy()
    if not fallback_merged.empty:
        fallback_by_row = fallback_merged.set_index("_row_id")
        for column in annotation_cols:
            if column not in fallback_by_row.columns:
                continue
            replacement = output["_row_id"].map(fallback_by_row[column])
            output[column] = output[column].where(output[column].notna(), replacement)

    output["gencode_annotation_match_type"] = "unmatched"
    output.loc[matched_by_version, "gencode_annotation_match_type"] = (
        "versioned_transcript_id"
    )
    fallback_match = (
        output["gencode_annotation_match_type"].eq("unmatched")
        & output["gencode_gene_id"].notna()
    )
    output.loc[fallback_match, "gencode_annotation_match_type"] = (
        "unversioned_transcript_id"
    )

    helper_cols = [
        "_join_transcript_id_with_version",
        "_join_transcript_id",
        "_row_id",
    ]
    output = output.drop(columns=[col for col in helper_cols if col in output.columns])
    matched_count = int(output["gencode_annotation_match_type"].ne("unmatched").sum())
    logger.info(
        "Annotated %d of %d result rows (%.2f%%)",
        matched_count,
        output.shape[0],
        100 * matched_count / max(output.shape[0], 1),
    )
    return output


def annotate_optional_table(
    *,
    input_path: Optional[Path],
    output_path: Path,
    transcript_annotation: pd.DataFrame,
    logger: logging.Logger,
) -> Optional[pd.DataFrame]:
    """
    Annotate an optional GTEx result table if a path was supplied.

    Parameters
    ----------
    input_path
        Optional input table path.
    output_path
        Output path.
    transcript_annotation
        GENCODE transcript annotation table.
    logger
        Logger instance.

    Returns
    -------
    Optional[pandas.DataFrame]
        Annotated table, or None if no input path was supplied.
    """
    if input_path is None:
        logger.info("No input supplied for %s; skipping", output_path.name)
        return None
    dataframe = read_tsv(path=input_path, logger=logger)
    annotated = annotate_result_table(
        dataframe=dataframe,
        transcript_annotation=transcript_annotation,
        logger=logger,
    )
    write_tsv(dataframe=annotated, path=output_path, logger=logger)
    return annotated


def run(*, config: Config, logger: logging.Logger) -> None:
    """
    Run GENCODE annotation workflow.

    Parameters
    ----------
    config
        Workflow configuration.
    logger
        Logger instance.
    """
    config.out_dir.mkdir(parents=True, exist_ok=True)
    prefix = config.out_dir / config.out_prefix

    with timed(logger=logger, label="Parse GENCODE GTF"):
        transcript_annotation, feature_table = parse_gencode_gtf(
            gtf_path=config.gtf_path,
            logger=logger,
        )

    with timed(logger=logger, label="Write parsed GENCODE annotation tables"):
        transcript_annotation_path = Path(
            f"{prefix}.gencode_transcript_annotation.tsv.gz"
        )
        feature_table_path = Path(f"{prefix}.gencode_transcript_features.tsv.gz")
        write_tsv(
            dataframe=transcript_annotation,
            path=transcript_annotation_path,
            logger=logger,
        )
        write_tsv(dataframe=feature_table, path=feature_table_path, logger=logger)

    with timed(logger=logger, label="Annotate GTEx result tables"):
        annotated_screen = annotate_optional_table(
            input_path=config.screen_tsv,
            output_path=Path(
                f"{prefix}.transcript_target_tissue_isoform_summary.annotated.tsv.gz"
            ),
            transcript_annotation=transcript_annotation,
            logger=logger,
        )
        annotated_candidates = annotate_optional_table(
            input_path=config.candidate_tsv,
            output_path=Path(
                f"{prefix}.candidate_target_tissue_isoforms.annotated.tsv"
            ),
            transcript_annotation=transcript_annotation,
            logger=logger,
        )
        annotated_best = annotate_optional_table(
            input_path=config.best_per_gene_tsv,
            output_path=Path(
                f"{prefix}.best_candidate_isoform_per_gene.annotated.tsv"
            ),
            transcript_annotation=transcript_annotation,
            logger=logger,
        )

    if config.write_excel_outputs:
        with timed(logger=logger, label="Write formatted Excel outputs"):
            if annotated_candidates is not None:
                write_formatted_excel(
                    dataframe=annotated_candidates,
                    path=Path(
                        f"{prefix}.candidate_target_tissue_isoforms.annotated.xlsx"
                    ),
                    sheet_name="candidate_isoforms_annotated",
                    logger=logger,
                    excel_max_rows=config.excel_max_rows,
                )
            if annotated_best is not None:
                write_formatted_excel(
                    dataframe=annotated_best,
                    path=Path(
                        f"{prefix}.best_candidate_isoform_per_gene.annotated.xlsx"
                    ),
                    sheet_name="best_isoforms_annotated",
                    logger=logger,
                    excel_max_rows=config.excel_max_rows,
                )
            if annotated_screen is not None and annotated_screen.shape[0] <= 50000:
                write_formatted_excel(
                    dataframe=annotated_screen,
                    path=Path(
                        f"{prefix}.transcript_target_tissue_isoform_summary"
                        ".annotated.xlsx"
                    ),
                    sheet_name="transcript_summary_annotated",
                    logger=logger,
                    excel_max_rows=config.excel_max_rows,
                )
            elif annotated_screen is not None:
                logger.info(
                    "Skipping Excel copy of full annotated screen because it "
                    "has %d rows; TSV.GZ remains the primary output",
                    annotated_screen.shape[0],
                )

    logger.info("GENCODE annotation workflow finished")


def main() -> None:
    """Run the command-line entry point."""
    config = parse_args()
    logger = setup_logging(log_level=config.log_level, log_path=config.log_path)
    logger.info("GENCODE annotation workflow starting")
    logger.info("GTF path: %s", config.gtf_path)
    logger.info("Output directory: %s", config.out_dir)
    logger.info("Excel outputs enabled: %s", config.write_excel_outputs)
    run(config=config, logger=logger)


if __name__ == "__main__":
    main()
