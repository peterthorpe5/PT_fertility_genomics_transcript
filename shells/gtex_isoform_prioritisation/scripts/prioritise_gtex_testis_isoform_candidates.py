#!/usr/bin/env python3
"""
Prioritise GTEx testis-preferential isoform candidates.

This script ranks annotated GTEx transcript isoform candidates after the
transcriptome-wide isoform-usage screen and GENCODE annotation steps. It is
intended to help decide which genes are most worth manual review, gene-model
plotting, and later biological or structural follow-up.

The core biological use case is to identify genes where:

    1. the gene itself may not be strongly testis-specific at gene level; but
    2. one transcript isoform is preferentially used in testis; and
    3. that transcript is protein-coding or otherwise structurally interesting.

Primary inputs
--------------
An annotated best-per-gene table, for example:

    gtex_v11_transcriptome_testis_isoform_screen.
    best_candidate_isoform_per_gene.annotated.tsv

Primary outputs
---------------
<prefix>.prioritised_all.tsv
    All rows from the input table with prioritisation columns added.

<prefix>.summary_by_priority_tier.tsv
    Row counts and basic metric summaries per priority tier.

<prefix>.tier1_protein_coding_rescue_candidates.tsv
    Highest priority broad-gene rescue candidates that are protein-coding.

<prefix>.tier2_rescue_non_coding_or_uncertain_candidates.tsv
    Broad-gene rescue candidates where coding status is absent or uncertain.

<prefix>.tier3_protein_coding_high_confidence_candidates.tsv
    Strong protein-coding candidates that are not classified as broad-gene
    rescue candidates.

<prefix>.genes_to_plot_top_<N>.txt
    Ranked gene-symbol list for plotting.

The script writes tab-separated files for machine-readable outputs. Formatted
Excel copies of the main browseable outputs are written by default and can be
switched off using --no_write_excel_outputs.
"""

from __future__ import annotations

import argparse
import logging
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class Config:
    """
    Runtime configuration for candidate prioritisation.

    Attributes
    ----------
    input_tsv
        Annotated best-per-gene candidate TSV or TSV.GZ.
    out_dir
        Output directory.
    out_prefix
        Prefix used for output files.
    gene_symbol_col
        Column containing gene symbols.
    top_n_plot_genes
        Number of top-ranked genes to write to the primary plotting list.
    batch_size
        Number of genes per plotting batch list.
    selected_genes
        Optional genes of project interest to extract into a separate table.
    min_strong_tpm
        Minimum target-tissue median TPM for high-confidence expression.
    min_strong_usage
        Minimum target-tissue isoform usage for high-confidence expression.
    min_strong_log2_usage_ratio
        Minimum target versus max non-target log2 isoform-usage ratio.
    min_detection_fraction
        Minimum fraction of target-tissue samples with TPM above threshold.
    rescue_column
        Column flagging broad-gene isoform-rescue candidates.
    candidate_column
        Column flagging target-tissue isoform candidates.
    write_excel_outputs
        Whether to write formatted Excel copies of main browseable outputs.
    log_path
        Optional path to a log file.
    log_level
        Logging level.
    """

    input_tsv: Path
    out_dir: Path
    out_prefix: str
    gene_symbol_col: str
    top_n_plot_genes: int
    batch_size: int
    selected_genes: Tuple[str, ...]
    min_strong_tpm: float
    min_strong_usage: float
    min_strong_log2_usage_ratio: float
    min_detection_fraction: float
    rescue_column: str
    candidate_column: str
    write_excel_outputs: bool
    log_path: Optional[Path]
    log_level: str


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
    logger = logging.getLogger("prioritise_gtex_testis_isoform_candidates")
    logger.setLevel(getattr(logging, log_level.upper(), logging.INFO))
    logger.handlers = []
    logger.propagate = False

    formatter = logging.Formatter(
        fmt="%(asctime)s | %(levelname)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    stream_handler = logging.StreamHandler(stream=sys.stdout)
    stream_handler.setFormatter(fmt=formatter)
    stream_handler.setLevel(getattr(logging, log_level.upper(), logging.INFO))
    logger.addHandler(hdlr=stream_handler)

    if log_path is not None:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(filename=log_path)
        file_handler.setFormatter(fmt=formatter)
        file_handler.setLevel(getattr(logging, log_level.upper(), logging.INFO))
        logger.addHandler(hdlr=file_handler)

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
        description="Prioritise annotated GTEx testis isoform candidates."
    )
    parser.add_argument(
        "--input_tsv",
        required=True,
        type=Path,
        help="Annotated best-per-gene candidate TSV or TSV.GZ.",
    )
    parser.add_argument(
        "--out_dir",
        required=True,
        type=Path,
        help="Output directory.",
    )
    parser.add_argument(
        "--out_prefix",
        default="gtex_v11_testis_isoform_prioritised",
        help="Output filename prefix.",
    )
    parser.add_argument(
        "--gene_symbol_col",
        default="gene_symbol",
        help="Gene symbol column. Default: gene_symbol.",
    )
    parser.add_argument(
        "--top_n_plot_genes",
        default=100,
        type=int,
        help="Number of top-ranked genes to write for plotting. Default: 100.",
    )
    parser.add_argument(
        "--batch_size",
        default=25,
        type=int,
        help="Number of genes per plotting batch list. Default: 25.",
    )
    parser.add_argument(
        "--selected_genes",
        nargs="*",
        default=("AFG2B", "CFAP99", "SLC16A7", "ABCG4"),
        help=(
            "Optional gene symbols to extract into a separate check table. "
            "Default: AFG2B CFAP99 SLC16A7 ABCG4."
        ),
    )
    parser.add_argument(
        "--min_strong_tpm",
        default=1.0,
        type=float,
        help="Minimum target-tissue median TPM for high confidence.",
    )
    parser.add_argument(
        "--min_strong_usage",
        default=0.25,
        type=float,
        help="Minimum target-tissue isoform usage for high confidence.",
    )
    parser.add_argument(
        "--min_strong_log2_usage_ratio",
        default=1.0,
        type=float,
        help="Minimum log2 target versus max non-target usage ratio.",
    )
    parser.add_argument(
        "--min_detection_fraction",
        default=0.20,
        type=float,
        help=(
            "Minimum fraction of target-tissue samples above TPM threshold "
            "for high confidence. Default: 0.20."
        ),
    )
    parser.add_argument(
        "--rescue_column",
        default="is_broad_gene_isoform_rescue_candidate",
        help="Column flagging broad-gene isoform-rescue candidates.",
    )
    parser.add_argument(
        "--candidate_column",
        default="is_target_tissue_isoform_candidate",
        help="Column flagging target-tissue isoform candidates.",
    )
    parser.add_argument(
        "--no_write_excel_outputs",
        action="store_true",
        help="Do not write formatted Excel copies of browseable outputs.",
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

    args = parser.parse_args(args=argv)
    return Config(
        input_tsv=args.input_tsv,
        out_dir=args.out_dir,
        out_prefix=args.out_prefix,
        gene_symbol_col=args.gene_symbol_col,
        top_n_plot_genes=args.top_n_plot_genes,
        batch_size=args.batch_size,
        selected_genes=tuple(args.selected_genes),
        min_strong_tpm=args.min_strong_tpm,
        min_strong_usage=args.min_strong_usage,
        min_strong_log2_usage_ratio=args.min_strong_log2_usage_ratio,
        min_detection_fraction=args.min_detection_fraction,
        rescue_column=args.rescue_column,
        candidate_column=args.candidate_column,
        write_excel_outputs=not args.no_write_excel_outputs,
        log_path=args.log_path,
        log_level=args.log_level,
    )


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
    logger.info("Reading %s", path)
    dataframe = pd.read_csv(filepath_or_buffer=path, sep="\t", low_memory=False)
    logger.info(
        "Loaded %d rows x %d columns from %s",
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
        DataFrame to write.
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


def parse_boolean_series(*, series: pd.Series) -> pd.Series:
    """
    Convert a mixed boolean-like series to bool.

    Parameters
    ----------
    series
        Series containing values such as 1, 0, True, False, yes, or no.

    Returns
    -------
    pandas.Series
        Boolean series.
    """
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(value=False).astype(bool)

    if pd.api.types.is_numeric_dtype(series):
        return series.fillna(value=0).astype(float) != 0

    truthy = {"1", "true", "t", "yes", "y", "present"}
    falsey = {"0", "false", "f", "no", "n", "absent", "", "nan", "none"}

    values = series.fillna(value="").astype(str).str.strip().str.lower()
    output = values.isin(values=truthy)
    unresolved = ~(values.isin(values=truthy) | values.isin(values=falsey))
    if unresolved.any():
        output.loc[unresolved] = False
    return output.astype(bool)


def ensure_numeric_column(
    *, dataframe: pd.DataFrame, column: str, default: float = 0.0
) -> pd.Series:
    """
    Return a numeric column, using a default if the column is absent.

    Parameters
    ----------
    dataframe
        Input DataFrame.
    column
        Column name.
    default
        Default numeric value used when the column is absent.

    Returns
    -------
    pandas.Series
        Numeric series aligned to the DataFrame index.
    """
    if column not in dataframe.columns:
        return pd.Series(data=default, index=dataframe.index, dtype=float)
    return pd.to_numeric(arg=dataframe[column], errors="coerce").fillna(value=default)


def ensure_bool_column(*, dataframe: pd.DataFrame, column: str) -> pd.Series:
    """
    Return a boolean column, defaulting to False if absent.

    Parameters
    ----------
    dataframe
        Input DataFrame.
    column
        Column name.

    Returns
    -------
    pandas.Series
        Boolean series aligned to the DataFrame index.
    """
    if column not in dataframe.columns:
        return pd.Series(data=False, index=dataframe.index, dtype=bool)
    return parse_boolean_series(series=dataframe[column])


def cap_scale_series(
    *, series: pd.Series, lower: float, upper: float, default: float = 0.0
) -> pd.Series:
    """
    Scale a numeric series to the interval 0 to 1 after clipping.

    Parameters
    ----------
    series
        Numeric series.
    lower
        Lower clipping bound.
    upper
        Upper clipping bound.
    default
        Value used when upper equals lower.

    Returns
    -------
    pandas.Series
        Scaled series.
    """
    if math.isclose(a=upper, b=lower):
        return pd.Series(data=default, index=series.index, dtype=float)
    clipped = series.clip(lower=lower, upper=upper)
    return (clipped - lower) / (upper - lower)


def add_derived_flags(
    *, dataframe: pd.DataFrame,
    rescue_column: str,
    candidate_column: str,
    min_strong_tpm: float,
    min_strong_usage: float,
    min_strong_log2_usage_ratio: float,
    min_detection_fraction: float,
    logger: logging.Logger,
) -> pd.DataFrame:
    """
    Add derived biological and prioritisation flags.

    Parameters
    ----------
    dataframe
        Annotated candidate table.
    rescue_column
        Column flagging broad-gene isoform-rescue candidates.
    candidate_column
        Column flagging target-tissue isoform candidates.
    min_strong_tpm
        Minimum target median TPM for high-confidence expression.
    min_strong_usage
        Minimum target isoform usage for high-confidence expression.
    min_strong_log2_usage_ratio
        Minimum log2 usage ratio for high-confidence expression.
    min_detection_fraction
        Minimum target-tissue detection fraction.
    logger
        Logger instance.

    Returns
    -------
    pandas.DataFrame
        DataFrame with derived flag columns.
    """
    output = dataframe.copy()

    output["flag_target_isoform_candidate"] = ensure_bool_column(
        dataframe=output, column=candidate_column
    )
    output["flag_broad_gene_isoform_rescue"] = ensure_bool_column(
        dataframe=output, column=rescue_column
    )
    output["flag_has_cds"] = ensure_bool_column(
        dataframe=output, column="gencode_has_cds"
    )
    output["flag_protein_coding_transcript"] = ensure_bool_column(
        dataframe=output, column="gencode_is_protein_coding_transcript"
    )
    output["flag_basic_transcript"] = ensure_bool_column(
        dataframe=output, column="is_basic"
    )
    output["flag_mane_select"] = ensure_bool_column(
        dataframe=output, column="is_mane_select"
    )
    output["flag_ensembl_canonical"] = ensure_bool_column(
        dataframe=output, column="is_ensembl_canonical"
    )
    output["flag_primary_rank_1"] = (
        output.get("candidate_rank_tier", "") == "primary_rank_1"
    )
    output["flag_secondary_rank_2_to_3"] = (
        output.get("candidate_rank_tier", "") == "secondary_rank_2_to_3"
    )

    target_tpm = ensure_numeric_column(dataframe=output, column="target_median_tpm")
    target_usage = ensure_numeric_column(
        dataframe=output, column="target_median_isoform_usage"
    )
    log2_usage_ratio = ensure_numeric_column(
        dataframe=output,
        column="log2_target_vs_max_non_target_isoform_usage",
    )
    detection_fraction = ensure_numeric_column(
        dataframe=output,
        column="target_fraction_samples_tpm_ge_threshold",
    )

    output["flag_high_confidence_expression"] = (
        (target_tpm >= min_strong_tpm)
        & (target_usage >= min_strong_usage)
        & (log2_usage_ratio >= min_strong_log2_usage_ratio)
        & (detection_fraction >= min_detection_fraction)
    )
    output["flag_comp_chem_relevant_coding_isoform"] = (
        output["flag_protein_coding_transcript"] & output["flag_has_cds"]
    )
    output["flag_noncanonical_or_alternative_coding"] = (
        output["flag_comp_chem_relevant_coding_isoform"]
        & ~output["flag_mane_select"]
    )

    logger.info(
        "Derived flags: %d rescue, %d protein-coding, %d high-confidence rows",
        int(output["flag_broad_gene_isoform_rescue"].sum()),
        int(output["flag_comp_chem_relevant_coding_isoform"].sum()),
        int(output["flag_high_confidence_expression"].sum()),
    )
    return output


def assign_priority_tiers(*, dataframe: pd.DataFrame) -> pd.DataFrame:
    """
    Assign readable priority tiers.

    Parameters
    ----------
    dataframe
        Candidate DataFrame with derived flags.

    Returns
    -------
    pandas.DataFrame
        DataFrame with priority_tier and priority_tier_rank columns.
    """
    output = dataframe.copy()
    output["priority_tier_rank"] = 5
    output["priority_tier"] = "tier_5_other_candidate"

    tier4 = output["flag_target_isoform_candidate"]
    output.loc[tier4, "priority_tier_rank"] = 4
    output.loc[tier4, "priority_tier"] = "tier_4_other_testis_isoform_candidate"

    tier3 = (
        output["flag_comp_chem_relevant_coding_isoform"]
        & output["flag_high_confidence_expression"]
        & ~output["flag_broad_gene_isoform_rescue"]
    )
    output.loc[tier3, "priority_tier_rank"] = 3
    output.loc[tier3, "priority_tier"] = (
        "tier_3_protein_coding_high_confidence_candidate"
    )

    tier2 = (
        output["flag_broad_gene_isoform_rescue"]
        & output["flag_high_confidence_expression"]
        & ~output["flag_comp_chem_relevant_coding_isoform"]
    )
    output.loc[tier2, "priority_tier_rank"] = 2
    output.loc[tier2, "priority_tier"] = (
        "tier_2_rescue_non_coding_or_uncertain_candidate"
    )

    tier1 = (
        output["flag_broad_gene_isoform_rescue"]
        & output["flag_high_confidence_expression"]
        & output["flag_comp_chem_relevant_coding_isoform"]
    )
    output.loc[tier1, "priority_tier_rank"] = 1
    output.loc[tier1, "priority_tier"] = (
        "tier_1_protein_coding_rescue_candidate"
    )

    return output


def add_priority_score(*, dataframe: pd.DataFrame) -> pd.DataFrame:
    """
    Add a continuous priority score for within-tier ranking.

    Parameters
    ----------
    dataframe
        Candidate DataFrame with derived flags and tiers.

    Returns
    -------
    pandas.DataFrame
        DataFrame with priority_score added.
    """
    output = dataframe.copy()

    target_tpm = ensure_numeric_column(dataframe=output, column="target_median_tpm")
    target_usage = ensure_numeric_column(
        dataframe=output, column="target_median_isoform_usage"
    )
    log2_usage_ratio = ensure_numeric_column(
        dataframe=output,
        column="log2_target_vs_max_non_target_isoform_usage",
    )
    detection_fraction = ensure_numeric_column(
        dataframe=output,
        column="target_fraction_samples_tpm_ge_threshold",
    )
    target_tpm_rank = ensure_numeric_column(
        dataframe=output,
        column="target_tpm_rank_within_gene",
        default=99.0,
    )
    target_usage_rank = ensure_numeric_column(
        dataframe=output,
        column="target_usage_rank_within_gene",
        default=99.0,
    )

    scaled_log2_usage = cap_scale_series(
        series=log2_usage_ratio, lower=0.0, upper=8.0
    )
    scaled_tpm = cap_scale_series(
        series=np.log2(target_tpm + 1.0), lower=0.0, upper=8.0
    )
    scaled_usage = cap_scale_series(series=target_usage, lower=0.0, upper=1.0)
    scaled_detection = cap_scale_series(
        series=detection_fraction, lower=0.0, upper=1.0
    )
    scaled_rank_bonus = cap_scale_series(
        series=4.0 - target_usage_rank.clip(lower=1.0, upper=4.0),
        lower=0.0,
        upper=3.0,
    )

    score = pd.Series(data=0.0, index=output.index, dtype=float)
    score += output["flag_broad_gene_isoform_rescue"].astype(float) * 30.0
    score += output["flag_comp_chem_relevant_coding_isoform"].astype(float) * 20.0
    score += output["flag_high_confidence_expression"].astype(float) * 15.0
    score += output["flag_primary_rank_1"].astype(float) * 8.0
    score += output["flag_noncanonical_or_alternative_coding"].astype(float) * 4.0
    score += output["flag_basic_transcript"].astype(float) * 2.0
    score += output["flag_mane_select"].astype(float) * 1.0
    score += scaled_log2_usage * 12.0
    score += scaled_usage * 12.0
    score += scaled_tpm * 10.0
    score += scaled_detection * 7.0
    score += scaled_rank_bonus * 5.0
    score += (1.0 / target_tpm_rank.clip(lower=1.0, upper=99.0)) * 2.0

    output["priority_score"] = score.round(decimals=4)
    return output


def sort_prioritised_table(*, dataframe: pd.DataFrame) -> pd.DataFrame:
    """
    Sort prioritised candidates for manual review.

    Parameters
    ----------
    dataframe
        Prioritised DataFrame.

    Returns
    -------
    pandas.DataFrame
        Sorted DataFrame.
    """
    sort_columns = [
        "priority_tier_rank",
        "priority_score",
        "flag_broad_gene_isoform_rescue",
        "flag_comp_chem_relevant_coding_isoform",
        "log2_target_vs_max_non_target_isoform_usage",
        "target_median_isoform_usage",
        "target_median_tpm",
    ]
    available_sort_columns = [
        column for column in sort_columns if column in dataframe.columns
    ]
    ascending = [True, False, False, False, False, False, False][
        : len(available_sort_columns)
    ]
    return dataframe.sort_values(
        by=available_sort_columns,
        ascending=ascending,
        kind="mergesort",
    ).reset_index(drop=True)


def prioritise_candidates(
    *, dataframe: pd.DataFrame, config: Config, logger: logging.Logger
) -> pd.DataFrame:
    """
    Add derived flags, tiers, scores, and sorting to candidate table.

    Parameters
    ----------
    dataframe
        Annotated best-per-gene candidate table.
    config
        Runtime configuration.
    logger
        Logger instance.

    Returns
    -------
    pandas.DataFrame
        Prioritised candidate table.
    """
    if config.gene_symbol_col not in dataframe.columns:
        raise ValueError(
            f"Gene symbol column {config.gene_symbol_col!r} was not present."
        )

    logger.info("Adding derived flags")
    prioritised = add_derived_flags(
        dataframe=dataframe,
        rescue_column=config.rescue_column,
        candidate_column=config.candidate_column,
        min_strong_tpm=config.min_strong_tpm,
        min_strong_usage=config.min_strong_usage,
        min_strong_log2_usage_ratio=config.min_strong_log2_usage_ratio,
        min_detection_fraction=config.min_detection_fraction,
        logger=logger,
    )
    logger.info("Assigning priority tiers")
    prioritised = assign_priority_tiers(dataframe=prioritised)
    logger.info("Calculating priority score")
    prioritised = add_priority_score(dataframe=prioritised)
    prioritised = sort_prioritised_table(dataframe=prioritised)
    return prioritised


def summarise_by_tier(*, prioritised: pd.DataFrame) -> pd.DataFrame:
    """
    Summarise prioritised candidates by priority tier.

    Parameters
    ----------
    prioritised
        Prioritised candidate table.

    Returns
    -------
    pandas.DataFrame
        Summary table by priority tier.
    """
    metrics = [
        "priority_score",
        "target_median_tpm",
        "target_median_isoform_usage",
        "log2_target_vs_max_non_target_isoform_usage",
    ]
    available_metrics = [column for column in metrics if column in prioritised.columns]

    grouped = prioritised.groupby(
        by=["priority_tier_rank", "priority_tier"],
        dropna=False,
    )
    rows = []
    for (tier_rank, tier), group in grouped:
        row = {
            "priority_tier_rank": tier_rank,
            "priority_tier": tier,
            "n_rows": group.shape[0],
            "n_unique_genes": group["gene_symbol"].nunique()
            if "gene_symbol" in group.columns
            else group.shape[0],
            "n_rescue": int(group["flag_broad_gene_isoform_rescue"].sum())
            if "flag_broad_gene_isoform_rescue" in group.columns
            else 0,
            "n_protein_coding": int(
                group["flag_comp_chem_relevant_coding_isoform"].sum()
            )
            if "flag_comp_chem_relevant_coding_isoform" in group.columns
            else 0,
        }
        for metric in available_metrics:
            row[f"median_{metric}"] = group[metric].median()
            row[f"max_{metric}"] = group[metric].max()
        rows.append(row)

    return pd.DataFrame(data=rows).sort_values(
        by=["priority_tier_rank"],
        ascending=[True],
    )


def select_selected_genes(
    *, prioritised: pd.DataFrame, gene_symbol_col: str, selected_genes: Sequence[str]
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Extract project-selected genes and report missing selected genes.

    Parameters
    ----------
    prioritised
        Prioritised candidate table.
    gene_symbol_col
        Gene symbol column.
    selected_genes
        Gene symbols requested by the user.

    Returns
    -------
    tuple[pandas.DataFrame, pandas.DataFrame]
        Selected-gene rows and missing-gene report.
    """
    requested = pd.DataFrame(
        data={"requested_gene_symbol": [str(gene) for gene in selected_genes]}
    )
    if not selected_genes:
        return prioritised.iloc[0:0].copy(), requested.iloc[0:0].copy()

    gene_set = {str(gene).upper() for gene in selected_genes}
    observed = prioritised[gene_symbol_col].fillna(value="").astype(str)
    selected = prioritised.loc[observed.str.upper().isin(values=gene_set)].copy()
    found_set = {str(gene).upper() for gene in selected[gene_symbol_col].dropna()}
    missing = requested.loc[
        ~requested["requested_gene_symbol"].str.upper().isin(values=found_set)
    ].copy()
    return selected, missing


def write_gene_list(*, genes: Iterable[str], path: Path, logger: logging.Logger) -> None:
    """
    Write one gene symbol per line.

    Parameters
    ----------
    genes
        Iterable of gene symbols.
    path
        Output path.
    logger
        Logger instance.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    clean_genes = [str(gene).strip() for gene in genes if str(gene).strip()]
    with path.open(mode="wt") as handle:
        for gene in clean_genes:
            handle.write(f"{gene}\n")
    logger.info("Wrote %d gene symbols to %s", len(clean_genes), path)


def make_unique_gene_list(
    *, dataframe: pd.DataFrame, gene_symbol_col: str, limit: Optional[int] = None
) -> List[str]:
    """
    Make a unique gene-symbol list preserving row order.

    Parameters
    ----------
    dataframe
        Input DataFrame.
    gene_symbol_col
        Gene symbol column.
    limit
        Optional maximum number of genes.

    Returns
    -------
    list[str]
        Unique gene symbols.
    """
    genes: List[str] = []
    seen = set()
    for value in dataframe[gene_symbol_col].fillna(value="").astype(str):
        gene = value.strip()
        if not gene or gene in seen:
            continue
        genes.append(gene)
        seen.add(gene)
        if limit is not None and len(genes) >= limit:
            break
    return genes


def write_plot_batches(
    *, genes: Sequence[str], batch_size: int, out_dir: Path, logger: logging.Logger
) -> None:
    """
    Write batch gene lists for plotting jobs.

    Parameters
    ----------
    genes
        Ordered gene list.
    batch_size
        Number of genes per batch.
    out_dir
        Output directory for batch files.
    logger
        Logger instance.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    if batch_size <= 0:
        raise ValueError("batch_size must be positive")

    n_batches = int(math.ceil(len(genes) / batch_size)) if genes else 0
    for batch_index in range(n_batches):
        start = batch_index * batch_size
        end = start + batch_size
        batch_genes = genes[start:end]
        path = out_dir / f"genes_to_plot_batch_{batch_index + 1:03d}.txt"
        write_gene_list(genes=batch_genes, path=path, logger=logger)

    logger.info("Wrote %d plotting batch files to %s", n_batches, out_dir)


def make_safe_excel_sheet_name(*, sheet_name: str) -> str:
    """
    Make a safe Excel worksheet name.

    Parameters
    ----------
    sheet_name
        Proposed sheet name.

    Returns
    -------
    str
        Excel-safe sheet name no longer than 31 characters.
    """
    invalid_characters = ["[", "]", ":", "*", "?", "/", "\\"]
    safe_name = sheet_name
    for invalid_character in invalid_characters:
        safe_name = safe_name.replace(invalid_character, "_")
    safe_name = safe_name.strip()
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
    Estimate readable Excel column widths.

    Parameters
    ----------
    dataframe
        DataFrame to inspect.
    max_scan_rows
        Maximum number of rows to scan.
    min_width
        Minimum Excel column width.
    max_text_width
        Maximum width for text columns.
    max_numeric_width
        Maximum width for numeric columns.

    Returns
    -------
    dict[str, int]
        Column-width mapping.
    """
    widths: Dict[str, int] = {}
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
    Write a formatted Excel copy of a result table.

    Parameters
    ----------
    dataframe
        DataFrame to write.
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

    if n_rows + 1 > excel_max_rows:
        logger.warning(
            "Skipping Excel output for %s because it has %d rows, exceeding "
            "the Excel worksheet limit.",
            path,
            n_rows,
        )
        return
    if n_cols > excel_max_cols:
        logger.warning(
            "Skipping Excel output for %s because it has %d columns, exceeding "
            "the Excel worksheet limit.",
            path,
            n_cols,
        )
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
                table_columns = [{"header": str(column)} for column in dataframe.columns]
                worksheet.add_table(
                    first_row=0,
                    first_col=0,
                    last_row=n_rows,
                    last_col=n_cols - 1,
                    options={
                        "columns": table_columns,
                        "style": "Table Style Medium 2",
                        "autofilter": True,
                    },
                )
            worksheet.set_row(row=0, height=30, cell_format=header_format)
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
        logger.warning(
            "Could not write Excel output %s because xlsxwriter is unavailable: %s",
            path,
            error,
        )
        return

    logger.info(
        "Wrote formatted Excel file with %d rows x %d columns to %s",
        n_rows,
        n_cols,
        path,
    )


def write_outputs(
    *, prioritised: pd.DataFrame, config: Config, logger: logging.Logger
) -> None:
    """
    Write prioritised tables, summaries, and gene lists.

    Parameters
    ----------
    prioritised
        Prioritised candidate table.
    config
        Runtime configuration.
    logger
        Logger instance.
    """
    config.out_dir.mkdir(parents=True, exist_ok=True)
    prefix = config.out_dir / config.out_prefix

    summary = summarise_by_tier(prioritised=prioritised)
    tier1 = prioritised.loc[
        prioritised["priority_tier"] == "tier_1_protein_coding_rescue_candidate"
    ].copy()
    tier2 = prioritised.loc[
        prioritised["priority_tier"]
        == "tier_2_rescue_non_coding_or_uncertain_candidate"
    ].copy()
    tier3 = prioritised.loc[
        prioritised["priority_tier"]
        == "tier_3_protein_coding_high_confidence_candidate"
    ].copy()
    rescue = prioritised.loc[
        prioritised["flag_broad_gene_isoform_rescue"]
    ].copy()
    protein_coding = prioritised.loc[
        prioritised["flag_comp_chem_relevant_coding_isoform"]
    ].copy()
    selected, missing = select_selected_genes(
        prioritised=prioritised,
        gene_symbol_col=config.gene_symbol_col,
        selected_genes=config.selected_genes,
    )

    write_tsv(
        dataframe=prioritised,
        path=Path(f"{prefix}.prioritised_all.tsv"),
        logger=logger,
    )
    write_tsv(
        dataframe=summary,
        path=Path(f"{prefix}.summary_by_priority_tier.tsv"),
        logger=logger,
    )
    write_tsv(
        dataframe=tier1,
        path=Path(f"{prefix}.tier1_protein_coding_rescue_candidates.tsv"),
        logger=logger,
    )
    write_tsv(
        dataframe=tier2,
        path=Path(f"{prefix}.tier2_rescue_non_coding_or_uncertain_candidates.tsv"),
        logger=logger,
    )
    write_tsv(
        dataframe=tier3,
        path=Path(f"{prefix}.tier3_protein_coding_high_confidence_candidates.tsv"),
        logger=logger,
    )
    write_tsv(
        dataframe=rescue,
        path=Path(f"{prefix}.all_rescue_candidates.tsv"),
        logger=logger,
    )
    write_tsv(
        dataframe=protein_coding,
        path=Path(f"{prefix}.all_protein_coding_candidates.tsv"),
        logger=logger,
    )
    write_tsv(
        dataframe=selected,
        path=Path(f"{prefix}.selected_project_genes_priority_check.tsv"),
        logger=logger,
    )
    write_tsv(
        dataframe=missing,
        path=Path(f"{prefix}.selected_project_genes_missing_from_best.tsv"),
        logger=logger,
    )

    top_genes = make_unique_gene_list(
        dataframe=prioritised,
        gene_symbol_col=config.gene_symbol_col,
        limit=config.top_n_plot_genes,
    )
    tier1_genes = make_unique_gene_list(
        dataframe=tier1,
        gene_symbol_col=config.gene_symbol_col,
        limit=None,
    )
    rescue_genes = make_unique_gene_list(
        dataframe=rescue,
        gene_symbol_col=config.gene_symbol_col,
        limit=None,
    )

    write_gene_list(
        genes=top_genes,
        path=Path(f"{prefix}.genes_to_plot_top_{config.top_n_plot_genes}.txt"),
        logger=logger,
    )
    write_gene_list(
        genes=tier1_genes,
        path=Path(f"{prefix}.genes_to_plot_tier1_protein_coding_rescue.txt"),
        logger=logger,
    )
    write_gene_list(
        genes=rescue_genes,
        path=Path(f"{prefix}.genes_to_plot_all_rescue.txt"),
        logger=logger,
    )
    write_plot_batches(
        genes=top_genes,
        batch_size=config.batch_size,
        out_dir=config.out_dir / "plot_gene_batches_top_ranked",
        logger=logger,
    )

    if config.write_excel_outputs:
        write_formatted_excel(
            dataframe=prioritised,
            path=Path(f"{prefix}.prioritised_all.xlsx"),
            sheet_name="prioritised_all",
            logger=logger,
        )
        write_formatted_excel(
            dataframe=summary,
            path=Path(f"{prefix}.summary_by_priority_tier.xlsx"),
            sheet_name="summary_by_tier",
            logger=logger,
        )
        write_formatted_excel(
            dataframe=tier1,
            path=Path(f"{prefix}.tier1_protein_coding_rescue_candidates.xlsx"),
            sheet_name="tier1_coding_rescue",
            logger=logger,
        )
        write_formatted_excel(
            dataframe=selected,
            path=Path(f"{prefix}.selected_project_genes_priority_check.xlsx"),
            sheet_name="selected_project_genes",
            logger=logger,
        )

    logger.info(
        "Prioritisation output complete: %d total rows, %d tier 1 rows, %d top plot genes",
        prioritised.shape[0],
        tier1.shape[0],
        len(top_genes),
    )


def run(*, config: Config, logger: logging.Logger) -> None:
    """
    Run the prioritisation workflow.

    Parameters
    ----------
    config
        Runtime configuration.
    logger
        Logger instance.
    """
    dataframe = read_tsv(path=config.input_tsv, logger=logger)
    prioritised = prioritise_candidates(
        dataframe=dataframe,
        config=config,
        logger=logger,
    )
    write_outputs(prioritised=prioritised, config=config, logger=logger)


def main() -> None:
    """Run the command-line entry point."""
    config = parse_args()
    logger = setup_logging(log_level=config.log_level, log_path=config.log_path)
    logger.info("GTEx testis isoform candidate prioritisation starting")
    logger.info("Input table: %s", config.input_tsv)
    logger.info("Output directory: %s", config.out_dir)
    logger.info("Excel outputs enabled: %s", config.write_excel_outputs)
    run(config=config, logger=logger)
    logger.info("GTEx testis isoform candidate prioritisation finished")


if __name__ == "__main__":
    main()
