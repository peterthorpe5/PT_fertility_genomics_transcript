#!/usr/bin/env python3
"""
Audit transcript isoform prioritisation and gene-level evidence integration.

This script is intended to catch manuscript-critical inconsistencies before
publication. In particular, it compares the strict transcript-level tier-1
candidate table against the later integrated evidence table and reports row
inflation, duplicate keys, tier reclassification, and protein-coding flag
conflicts.

The script does not change the input data. It writes diagnostic TSV files and,
by default, formatted Excel copies for manual review.
"""

from __future__ import annotations

import argparse
import logging
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import pandas as pd


LOGGER_NAME = "audit_isoform_evidence_integration"


@dataclass(frozen=True)
class Config:
    """Runtime configuration for the integration audit."""

    pre_integration_tier1: Path
    post_integration_table: Path
    out_dir: Path
    out_prefix: str
    gene_col: str
    transcript_col: str
    write_excel_outputs: bool
    log_level: str
    log_path: Optional[Path]


def setup_logging(*, log_level: str, log_path: Optional[Path]) -> logging.Logger:
    """
    Configure console and optional file logging.

    Parameters
    ----------
    log_level
        Logging level name, for example ``INFO`` or ``DEBUG``.
    log_path
        Optional file path for log output.

    Returns
    -------
    logging.Logger
        Configured logger.
    """
    logger = logging.getLogger(LOGGER_NAME)
    logger.handlers.clear()
    logger.setLevel(getattr(logging, log_level.upper(), logging.INFO))
    formatter = logging.Formatter(
        fmt="%(asctime)s | %(levelname)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
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
    """
    Parse command-line arguments.

    Parameters
    ----------
    argv
        Optional command-line argument sequence for testing.

    Returns
    -------
    Config
        Parsed configuration.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Audit consistency between the transcript-level tier-1 candidate "
            "table and the later gene-level evidence integration output."
        )
    )
    parser.add_argument(
        "--pre_integration_tier1",
        required=True,
        type=Path,
        help="Strict transcript-level tier-1 candidate table TSV/TSV.GZ/XLSX.",
    )
    parser.add_argument(
        "--post_integration_table",
        required=True,
        type=Path,
        help="Integrated gene-level evidence table TSV/TSV.GZ/XLSX.",
    )
    parser.add_argument(
        "--out_dir",
        required=True,
        type=Path,
        help="Output directory.",
    )
    parser.add_argument(
        "--out_prefix",
        default="isoform_gene_evidence_integration_audit",
        help="Output filename prefix.",
    )
    parser.add_argument(
        "--gene_col",
        default="gene_symbol",
        help="Gene symbol column name.",
    )
    parser.add_argument(
        "--transcript_col",
        default="transcript_id_with_version",
        help="Versioned transcript ID column name.",
    )
    parser.add_argument(
        "--no_write_excel_outputs",
        action="store_true",
        help="Disable formatted Excel output files.",
    )
    parser.add_argument(
        "--log_level",
        default="INFO",
        help="Logging level.",
    )
    parser.add_argument(
        "--log_path",
        type=Path,
        default=None,
        help="Optional log file path.",
    )

    args = parser.parse_args(argv)
    return Config(
        pre_integration_tier1=args.pre_integration_tier1,
        post_integration_table=args.post_integration_table,
        out_dir=args.out_dir,
        out_prefix=args.out_prefix,
        gene_col=args.gene_col,
        transcript_col=args.transcript_col,
        write_excel_outputs=not args.no_write_excel_outputs,
        log_level=args.log_level,
        log_path=args.log_path,
    )


def read_table(*, path: Path, logger: logging.Logger) -> pd.DataFrame:
    """
    Read a TSV, TSV.GZ or XLSX table.

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
    logger.info("Reading %s", path)
    suffixes = "".join(path.suffixes).lower()
    if suffixes.endswith(".xlsx"):
        dataframe = pd.read_excel(path)
    else:
        dataframe = pd.read_csv(path, sep="\t", low_memory=False)
    logger.info("Loaded %d rows x %d columns from %s", *dataframe.shape, path)
    return dataframe


def validate_columns(
    *, dataframe: pd.DataFrame, columns: Iterable[str], table_name: str
) -> None:
    """
    Raise a clear error if required columns are missing.

    Parameters
    ----------
    dataframe
        Table to inspect.
    columns
        Required column names.
    table_name
        Human-readable table name for error messages.
    """
    missing = [column for column in columns if column not in dataframe.columns]
    if missing:
        raise ValueError(f"{table_name} is missing required columns: {missing}")


def normalise_key_columns(
    *, dataframe: pd.DataFrame, gene_col: str, transcript_col: str
) -> pd.DataFrame:
    """
    Add normalised key columns for robust comparison.

    Parameters
    ----------
    dataframe
        Input table.
    gene_col
        Gene symbol column name.
    transcript_col
        Transcript ID column name.

    Returns
    -------
    pandas.DataFrame
        Copy with ``_audit_gene_key``, ``_audit_transcript_key`` and
        ``_audit_gene_transcript_key`` columns.
    """
    output = dataframe.copy()
    output["_audit_gene_key"] = output[gene_col].astype(str).str.strip()
    output["_audit_transcript_key"] = output[transcript_col].astype(str).str.strip()
    output["_audit_gene_transcript_key"] = (
        output["_audit_gene_key"] + "||" + output["_audit_transcript_key"]
    )
    return output


def summarise_table(
    *, dataframe: pd.DataFrame, table_name: str, gene_col: str, transcript_col: str
) -> Dict[str, object]:
    """
    Summarise row counts and key duplication for one table.

    Parameters
    ----------
    dataframe
        Input table.
    table_name
        Name for the table.
    gene_col
        Gene symbol column name.
    transcript_col
        Transcript ID column name.

    Returns
    -------
    Dict[str, object]
        Summary metrics.
    """
    return {
        "table": table_name,
        "n_rows": len(dataframe),
        "n_unique_genes": dataframe[gene_col].nunique(dropna=True),
        "n_unique_transcripts": dataframe[transcript_col].nunique(dropna=True),
        "n_unique_gene_transcript_keys": dataframe[
            "_audit_gene_transcript_key"
        ].nunique(dropna=True),
        "n_duplicate_gene_transcript_rows": dataframe.duplicated(
            "_audit_gene_transcript_key", keep=False
        ).sum(),
        "n_duplicate_gene_rows": dataframe.duplicated(gene_col, keep=False).sum(),
    }


def compare_key_sets(
    *, pre: pd.DataFrame, post: pd.DataFrame
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Compare gene-transcript keys between pre- and post-integration tables.

    Parameters
    ----------
    pre
        Pre-integration table with normalised key columns.
    post
        Post-integration table with normalised key columns.

    Returns
    -------
    Tuple[pandas.DataFrame, pandas.DataFrame, pandas.DataFrame]
        Summary table, post-only rows, and pre-only rows.
    """
    pre_keys = set(pre["_audit_gene_transcript_key"])
    post_keys = set(post["_audit_gene_transcript_key"])

    post_only_keys = post_keys - pre_keys
    pre_only_keys = pre_keys - post_keys
    shared_keys = pre_keys & post_keys

    summary = pd.DataFrame(
        [
            {
                "comparison": "gene_transcript_key_overlap",
                "n_pre_keys": len(pre_keys),
                "n_post_keys": len(post_keys),
                "n_shared_keys": len(shared_keys),
                "n_post_only_keys": len(post_only_keys),
                "n_pre_only_keys": len(pre_only_keys),
            }
        ]
    )

    post_only = post.loc[
        post["_audit_gene_transcript_key"].isin(post_only_keys)
    ].copy()
    pre_only = pre.loc[pre["_audit_gene_transcript_key"].isin(pre_only_keys)].copy()

    return summary, post_only, pre_only


def find_duplicate_rows(*, dataframe: pd.DataFrame) -> pd.DataFrame:
    """
    Return rows with duplicated gene-transcript keys.

    Parameters
    ----------
    dataframe
        Input table with normalised key columns.

    Returns
    -------
    pandas.DataFrame
        Duplicate rows, if any.
    """
    duplicate_mask = dataframe.duplicated("_audit_gene_transcript_key", keep=False)
    return dataframe.loc[duplicate_mask].copy()


def value_counts_table(
    *, dataframe: pd.DataFrame, column: str, table_name: str
) -> pd.DataFrame:
    """
    Summarise values in a column if it exists.

    Parameters
    ----------
    dataframe
        Input table.
    column
        Column to count.
    table_name
        Name to include in output.

    Returns
    -------
    pandas.DataFrame
        Value counts table. Empty if the column is absent.
    """
    if column not in dataframe.columns:
        return pd.DataFrame()
    counts = dataframe[column].value_counts(dropna=False).reset_index()
    counts.columns = [column, "n_rows"]
    counts.insert(0, "table", table_name)
    return counts


def find_protein_coding_flag_conflicts(*, dataframe: pd.DataFrame) -> pd.DataFrame:
    """
    Identify rows where protein-coding classifications disagree.

    The current project has several related columns. This function highlights
    cases where a harmonised/integrated flag marks a transcript as coding while
    strict GENCODE-derived flags do not.

    Parameters
    ----------
    dataframe
        Integrated table.

    Returns
    -------
    pandas.DataFrame
        Rows with potentially inconsistent protein-coding status.
    """
    if "isoform_is_protein_coding" not in dataframe.columns:
        return pd.DataFrame()

    integrated_true = dataframe["isoform_is_protein_coding"].astype(str).str.lower().isin(
        {"true", "1", "yes"}
    )

    strict_columns = [
        column
        for column in (
            "gencode_is_protein_coding_transcript",
            "flag_protein_coding_transcript",
        )
        if column in dataframe.columns
    ]
    if not strict_columns:
        return pd.DataFrame()

    strict_true = pd.Series(False, index=dataframe.index)
    for column in strict_columns:
        strict_true = strict_true | dataframe[column].astype(str).str.lower().isin(
            {"true", "1", "yes"}
        )

    conflict_mask = integrated_true & ~strict_true
    return dataframe.loc[conflict_mask].copy()


def find_unexpected_tier_rows(*, post: pd.DataFrame) -> pd.DataFrame:
    """
    Identify rows in a protein-coding rescue table that are not tier 1.

    Parameters
    ----------
    post
        Integrated post-table.

    Returns
    -------
    pandas.DataFrame
        Rows with ``priority_tier`` not equal to the strict tier-1 label.
    """
    if "priority_tier" not in post.columns:
        return pd.DataFrame()
    expected = "tier_1_protein_coding_rescue_candidate"
    return post.loc[post["priority_tier"] != expected].copy()


def write_tsv(*, dataframe: pd.DataFrame, path: Path, logger: logging.Logger) -> None:
    """
    Write a DataFrame as a tab-separated file.

    Parameters
    ----------
    dataframe
        DataFrame to write.
    path
        Output path.
    logger
        Logger instance.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    dataframe.to_csv(path, sep="\t", index=False)
    logger.info("Wrote %d rows x %d columns to %s", *dataframe.shape, path)


def estimate_column_widths(*, dataframe: pd.DataFrame, max_rows: int = 1000) -> Dict[str, int]:
    """
    Estimate Excel column widths.

    Parameters
    ----------
    dataframe
        Input DataFrame.
    max_rows
        Maximum rows to scan.

    Returns
    -------
    Dict[str, int]
        Column-width mapping.
    """
    widths: Dict[str, int] = {}
    scan = dataframe.head(max_rows)
    for column in dataframe.columns:
        header_width = len(str(column)) + 2
        if scan.empty:
            value_width = 0
        else:
            value_width = scan[column].astype(str).str.len().max()
            if pd.isna(value_width):
                value_width = 0
        cap = 18 if pd.api.types.is_numeric_dtype(dataframe[column]) else 45
        widths[str(column)] = min(max(header_width, int(value_width) + 2, 8), cap)
    return widths


def write_formatted_excel(
    *, dataframe: pd.DataFrame, path: Path, sheet_name: str, logger: logging.Logger
) -> None:
    """
    Write a formatted Excel copy of a table.

    Parameters
    ----------
    dataframe
        DataFrame to write.
    path
        Output path.
    sheet_name
        Worksheet name.
    logger
        Logger instance.
    """
    if len(dataframe) + 1 > 1_048_576:
        logger.warning("Skipping Excel output for %s: too many rows", path)
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    safe_sheet_name = sheet_name[:31] or "results"
    with pd.ExcelWriter(path, engine="xlsxwriter") as writer:
        dataframe.to_excel(writer, sheet_name=safe_sheet_name, index=False)
        workbook = writer.book
        worksheet = writer.sheets[safe_sheet_name]
        header_format = workbook.add_format(
            {"bold": True, "text_wrap": True, "valign": "top", "border": 1}
        )
        float_format = workbook.add_format({"num_format": "0.0000"})
        integer_format = workbook.add_format({"num_format": "0"})
        worksheet.freeze_panes(1, 0)
        worksheet.set_row(0, 30, header_format)
        if dataframe.shape[1] > 0:
            columns = [{"header": str(column)} for column in dataframe.columns]
            worksheet.add_table(
                0,
                0,
                len(dataframe),
                dataframe.shape[1] - 1,
                {"columns": columns, "autofilter": True, "style": "Table Style Medium 2"},
            )
        widths = estimate_column_widths(dataframe=dataframe)
        for idx, column in enumerate(dataframe.columns):
            fmt = None
            if pd.api.types.is_integer_dtype(dataframe[column]):
                fmt = integer_format
            elif pd.api.types.is_float_dtype(dataframe[column]):
                fmt = float_format
            worksheet.set_column(idx, idx, widths[str(column)], fmt)
    logger.info("Wrote formatted Excel file to %s", path)


def write_outputs(
    *, outputs: Dict[str, pd.DataFrame], config: Config, logger: logging.Logger
) -> None:
    """
    Write all diagnostic outputs.

    Parameters
    ----------
    outputs
        Mapping of output suffix to DataFrame.
    config
        Runtime configuration.
    logger
        Logger instance.
    """
    for suffix, dataframe in outputs.items():
        tsv_path = config.out_dir / f"{config.out_prefix}.{suffix}.tsv"
        write_tsv(dataframe=dataframe, path=tsv_path, logger=logger)
        if config.write_excel_outputs:
            excel_path = config.out_dir / f"{config.out_prefix}.{suffix}.xlsx"
            write_formatted_excel(
                dataframe=dataframe,
                path=excel_path,
                sheet_name=suffix,
                logger=logger,
            )


def run(*, config: Config, logger: logging.Logger) -> Dict[str, pd.DataFrame]:
    """
    Run the integration consistency audit.

    Parameters
    ----------
    config
        Runtime configuration.
    logger
        Logger instance.

    Returns
    -------
    Dict[str, pandas.DataFrame]
        Output tables.
    """
    pre = read_table(path=config.pre_integration_tier1, logger=logger)
    post = read_table(path=config.post_integration_table, logger=logger)

    validate_columns(
        dataframe=pre,
        columns=[config.gene_col, config.transcript_col],
        table_name="pre-integration table",
    )
    validate_columns(
        dataframe=post,
        columns=[config.gene_col, config.transcript_col],
        table_name="post-integration table",
    )

    pre = normalise_key_columns(
        dataframe=pre, gene_col=config.gene_col, transcript_col=config.transcript_col
    )
    post = normalise_key_columns(
        dataframe=post, gene_col=config.gene_col, transcript_col=config.transcript_col
    )

    table_summary = pd.DataFrame(
        [
            summarise_table(
                dataframe=pre,
                table_name="pre_integration_tier1",
                gene_col=config.gene_col,
                transcript_col=config.transcript_col,
            ),
            summarise_table(
                dataframe=post,
                table_name="post_integration_table",
                gene_col=config.gene_col,
                transcript_col=config.transcript_col,
            ),
        ]
    )

    overlap_summary, post_only, pre_only = compare_key_sets(pre=pre, post=post)
    post_duplicates = find_duplicate_rows(dataframe=post)
    pre_duplicates = find_duplicate_rows(dataframe=pre)
    unexpected_tier_rows = find_unexpected_tier_rows(post=post)
    protein_coding_conflicts = find_protein_coding_flag_conflicts(dataframe=post)

    value_count_tables = []
    for table_name, dataframe in (
        ("pre_integration_tier1", pre),
        ("post_integration_table", post),
    ):
        for column in ("priority_tier", "integrated_evidence_class"):
            counts = value_counts_table(
                dataframe=dataframe, column=column, table_name=table_name
            )
            if not counts.empty:
                value_count_tables.append(counts)
    value_counts = (
        pd.concat(value_count_tables, ignore_index=True)
        if value_count_tables
        else pd.DataFrame()
    )

    interpretation_rows = [
        {
            "finding": "post_rows_minus_pre_rows",
            "value": int(len(post) - len(pre)),
            "interpretation": (
                "Positive values indicate that the compared post-integration "
                "table contains rows not present in the strict pre-integration "
                "tier-1 table. Inspect post_only_rows and unexpected_tier_rows."
            ),
        },
        {
            "finding": "post_only_gene_transcript_keys",
            "value": int(len(post_only)),
            "interpretation": (
                "Rows present after integration but absent from the strict tier-1 "
                "table. If these are tier 2/non-coding rows, the issue is likely "
                "definition drift rather than one-to-many merge duplication."
            ),
        },
        {
            "finding": "post_duplicate_gene_transcript_rows",
            "value": int(len(post_duplicates)),
            "interpretation": (
                "Duplicate gene-transcript keys in the post table. Non-zero values "
                "would support a one-to-many merge duplication problem."
            ),
        },
        {
            "finding": "unexpected_non_tier1_rows_in_post_table",
            "value": int(len(unexpected_tier_rows)),
            "interpretation": (
                "Rows in a protein-coding rescue output that are not strict tier 1. "
                "This is evidence of tier-definition drift and should be resolved "
                "before manuscript counts are finalised."
            ),
        },
        {
            "finding": "protein_coding_flag_conflict_rows",
            "value": int(len(protein_coding_conflicts)),
            "interpretation": (
                "Rows marked protein-coding by an integrated flag but not by strict "
                "GENCODE-derived protein-coding flags. These should be reviewed."
            ),
        },
    ]
    interpretation = pd.DataFrame(interpretation_rows)

    outputs = {
        "table_summary": table_summary,
        "key_overlap_summary": overlap_summary,
        "post_only_rows": post_only,
        "pre_only_rows": pre_only,
        "post_duplicate_gene_transcript_rows": post_duplicates,
        "pre_duplicate_gene_transcript_rows": pre_duplicates,
        "unexpected_tier_rows": unexpected_tier_rows,
        "protein_coding_flag_conflicts": protein_coding_conflicts,
        "value_counts": value_counts,
        "interpretation_summary": interpretation,
    }
    write_outputs(outputs=outputs, config=config, logger=logger)
    return outputs


def main(argv: Optional[Sequence[str]] = None) -> None:
    """Run the command-line program."""
    config = parse_args(argv)
    logger = setup_logging(log_level=config.log_level, log_path=config.log_path)
    logger.info("Starting isoform integration consistency audit")
    run(config=config, logger=logger)
    logger.info("Finished isoform integration consistency audit")


if __name__ == "__main__":
    main()
