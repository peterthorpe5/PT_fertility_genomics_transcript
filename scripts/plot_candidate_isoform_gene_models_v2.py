#!/usr/bin/env python3
"""
Plot enhanced candidate isoform gene models.

This script adds a collapsed gene-level exon catalogue above the individual
transcript tracks. The collapsed track draws every annotated exon interval for
that gene as solid black boxes, so missing or alternative exons are easier to
see in each transcript model below.
"""

from __future__ import annotations

import argparse
import logging
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Rectangle
from matplotlib.transforms import blended_transform_factory


@dataclass(frozen=True)
class Config:
    """Configuration for enhanced isoform gene-model plotting."""

    features_tsv: Path
    candidate_tsv: Path
    out_dir: Path
    genes: Tuple[str, ...]
    gene_file: Optional[Path]
    output_formats: Tuple[str, ...]
    out_suffix: str
    max_transcripts: int
    label_union_exons: bool
    title_suffix: str
    log_path: Optional[Path]
    log_level: str
    dpi: int


def setup_logging(*, log_level: str, log_path: Optional[Path]) -> logging.Logger:
    """
    Configure stdout and optional file logging.

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
    logger = logging.getLogger("plot_candidate_isoform_gene_models_v2")
    logger.setLevel(getattr(logging, log_level.upper(), logging.INFO))
    logger.handlers = []
    logger.propagate = False

    formatter = logging.Formatter(
        fmt="%(asctime)s | %(levelname)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    stream_handler = logging.StreamHandler(stream=sys.stdout)
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)

    if log_path is not None:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(filename=log_path)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    return logger


def parse_args(argv: Optional[Sequence[str]] = None) -> Config:
    """
    Parse command-line arguments.

    Parameters
    ----------
    argv
        Optional command-line argument list.

    Returns
    -------
    Config
        Parsed configuration.
    """
    parser = argparse.ArgumentParser(
        description="Plot enhanced candidate isoform gene models."
    )
    parser.add_argument("--features_tsv", required=True, type=Path)
    parser.add_argument("--candidate_tsv", required=True, type=Path)
    parser.add_argument("--out_dir", required=True, type=Path)
    parser.add_argument("--genes", nargs="*", default=())
    parser.add_argument("--gene_file", type=Path, default=None)
    parser.add_argument("--output_formats", nargs="+", default=("pdf", "svg"))
    parser.add_argument(
        "--out_suffix",
        default="candidate_isoform_gene_model_v2",
    )
    parser.add_argument("--max_transcripts", type=int, default=40)
    parser.add_argument("--no_label_union_exons", action="store_true")
    parser.add_argument(
        "--title_suffix",
        default="GTEx v11 testis-preferential isoform usage",
    )
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--log_path", type=Path, default=None)
    parser.add_argument("--log_level", default="INFO")
    args = parser.parse_args(argv)

    return Config(
        features_tsv=args.features_tsv,
        candidate_tsv=args.candidate_tsv,
        out_dir=args.out_dir,
        genes=tuple(args.genes),
        gene_file=args.gene_file,
        output_formats=tuple(args.output_formats),
        out_suffix=args.out_suffix,
        max_transcripts=args.max_transcripts,
        label_union_exons=not args.no_label_union_exons,
        title_suffix=args.title_suffix,
        log_path=args.log_path,
        log_level=args.log_level,
        dpi=args.dpi,
    )


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
        Version-stripped identifier.
    """
    value = str(identifier).strip()
    if not value or value.lower() == "nan":
        return ""
    return value.split(".", maxsplit=1)[0]


def first_present_column(
    *, columns: Iterable[str], candidates: Sequence[str], required: bool = True
) -> Optional[str]:
    """
    Return the first matching column from a list of alternatives.

    Parameters
    ----------
    columns
        Available column names.
    candidates
        Candidate column names in priority order.
    required
        Whether absence should raise an error.

    Returns
    -------
    Optional[str]
        Matching column name or None.
    """
    available = set(columns)
    for candidate in candidates:
        if candidate in available:
            return candidate
    lower_to_original = {column.lower(): column for column in columns}
    for candidate in candidates:
        if candidate.lower() in lower_to_original:
            return lower_to_original[candidate.lower()]
    if required:
        raise ValueError(f"Missing required columns. Tried: {candidates}")
    return None


def read_tsv(*, path: Path, logger: logging.Logger) -> pd.DataFrame:
    """
    Read a tab-separated table.

    Parameters
    ----------
    path
        Input TSV or TSV.GZ file.
    logger
        Logger instance.

    Returns
    -------
    pandas.DataFrame
        Loaded table.
    """
    logger.info("Reading %s", path)
    dataframe = pd.read_csv(path, sep="\t", dtype=str)
    logger.info("Loaded %d rows x %d columns", dataframe.shape[0], dataframe.shape[1])
    return dataframe


def normalise_features(*, dataframe: pd.DataFrame) -> pd.DataFrame:
    """
    Normalise a GENCODE-derived feature table.

    The feature table may not always contain a gene-symbol column. This can
    happen when the feature output was written as a compact coordinate table.
    In that case, gene symbols are recovered later from the candidate table
    using gene IDs or candidate transcript IDs.

    Parameters
    ----------
    dataframe
        Raw feature table.

    Returns
    -------
    pandas.DataFrame
        Standardised feature table.
    """
    columns = dataframe.columns
    gene_symbol_col = first_present_column(
        columns=columns,
        candidates=(
            "gene_symbol",
            "gene_name",
            "gene",
            "gencode_gene_name",
        ),
        required=False,
    )
    gene_id_col = first_present_column(
        columns=columns, candidates=(
            "gene_id",
            "gene_id_with_version",
            "gencode_gene_id",
            "gencode_gene_id_with_version",
        )
    )
    transcript_id_col = first_present_column(
        columns=columns,
        candidates=(
            "transcript_id",
            "transcript_stable_id",
            "transcript",
            "gencode_transcript_id",
            "gencode_transcript_id_with_version",
        ),
    )
    transcript_id_version_col = first_present_column(
        columns=columns,
        candidates=(
            "transcript_id_with_version",
            "transcript_id_versioned",
            "transcript_id_full",
            "gencode_transcript_id_with_version",
            "gencode_transcript_id",
            transcript_id_col,
        ),
    )
    transcript_name_col = first_present_column(
        columns=columns,
        candidates=(
            "transcript_name",
            "transcript_label",
            "gencode_transcript_name",
        ),
        required=False,
    )
    feature_col = first_present_column(
        columns=columns, candidates=(
            "feature",
            "feature_type",
            "type",
            "gencode_feature_type",
        )
    )
    seqname_col = first_present_column(
        columns=columns, candidates=(
            "seqname",
            "chromosome",
            "chrom",
            "contig",
            "gencode_seqname",
        )
    )
    start_col = first_present_column(
        columns=columns, candidates=(
            "start",
            "feature_start",
            "genomic_start",
            "gencode_start",
        )
    )
    end_col = first_present_column(
        columns=columns, candidates=(
            "end",
            "feature_end",
            "genomic_end",
            "gencode_end",
        )
    )
    strand_col = first_present_column(
        columns=columns,
        candidates=("strand", "gencode_strand"),
    )

    if gene_symbol_col is None:
        gene_symbols = pd.Series([""] * dataframe.shape[0], index=dataframe.index)
    else:
        gene_symbols = dataframe[gene_symbol_col].fillna("").astype(str)

    output = pd.DataFrame(
        {
            "gene_symbol": gene_symbols,
            "gene_id": dataframe[gene_id_col].map(
                lambda value: strip_ensembl_version(identifier=value)
            ),
            "transcript_id_with_version": dataframe[
                transcript_id_version_col
            ].fillna("").astype(str),
            "transcript_id": dataframe[transcript_id_col].map(
                lambda value: strip_ensembl_version(identifier=value)
            ),
            "feature": dataframe[feature_col].fillna("").astype(str).str.lower(),
            "seqname": dataframe[seqname_col].fillna("").astype(str),
            "start": pd.to_numeric(dataframe[start_col], errors="coerce"),
            "end": pd.to_numeric(dataframe[end_col], errors="coerce"),
            "strand": dataframe[strand_col].fillna("").astype(str),
        }
    )
    if transcript_name_col is None:
        output["transcript_name"] = ""
    else:
        output["transcript_name"] = dataframe[transcript_name_col].fillna("").astype(str)
    output = output.dropna(subset=["start", "end"]).copy()
    output["start"] = output["start"].astype(int)
    output["end"] = output["end"].astype(int)
    return output


def normalise_candidates(*, dataframe: pd.DataFrame) -> pd.DataFrame:
    """
    Normalise an annotated candidate isoform table.

    Parameters
    ----------
    dataframe
        Raw candidate table.

    Returns
    -------
    pandas.DataFrame
        Standardised candidate table.
    """
    columns = dataframe.columns
    gene_symbol_col = first_present_column(
        columns=columns, candidates=("gene_symbol", "gene_name", "gene")
    )
    transcript_id_col = first_present_column(
        columns=columns,
        candidates=(
            "transcript_id",
            "transcript_stable_id",
            "transcript",
            "gencode_transcript_id",
            "gencode_transcript_id_with_version",
        ),
    )
    transcript_id_version_col = first_present_column(
        columns=columns,
        candidates=(
            "transcript_id_with_version",
            "transcript_id_versioned",
            "transcript_id_full",
            "gencode_transcript_id_with_version",
            "gencode_transcript_id",
            transcript_id_col,
        ),
    )
    gene_id_col = first_present_column(
        columns=columns,
        candidates=(
            "gene_id",
            "gene_id_with_version",
            "gencode_gene_id",
            "gencode_gene_id_with_version",
        ),
        required=False,
    )
    if gene_id_col is None:
        gene_ids = pd.Series([""] * dataframe.shape[0], index=dataframe.index)
    else:
        gene_ids = dataframe[gene_id_col].map(
            lambda value: strip_ensembl_version(identifier=value)
        )
    output = pd.DataFrame(
        {
            "gene_symbol": dataframe[gene_symbol_col].fillna("").astype(str),
            "gene_id": gene_ids,
            "transcript_id_with_version": dataframe[
                transcript_id_version_col
            ].fillna("").astype(str),
            "transcript_id": dataframe[transcript_id_col].map(
                lambda value: strip_ensembl_version(identifier=value)
            ),
        }
    )
    optional_map = {
        "candidate_rank_tier": ("candidate_rank_tier", "candidate_tier"),
        "target_median_tpm": ("target_median_tpm", "testis_median_tpm"),
        "target_median_isoform_usage": (
            "target_median_isoform_usage",
            "testis_median_isoform_usage",
        ),
        "log2_target_vs_max_non_target_isoform_usage": (
            "log2_target_vs_max_non_target_isoform_usage",
            "log2_testis_vs_max_non_testis_isoform_usage",
        ),
        "target_usage_rank_within_gene": ("target_usage_rank_within_gene",),
    }
    for output_col, possible_cols in optional_map.items():
        input_col = first_present_column(
            columns=columns, candidates=possible_cols, required=False
        )
        output[output_col] = "" if input_col is None else dataframe[input_col].fillna("").astype(str)
    return output


def merge_intervals(*, intervals: Iterable[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """
    Merge overlapping or adjacent genomic intervals.

    Parameters
    ----------
    intervals
        Start-end intervals.

    Returns
    -------
    list[tuple[int, int]]
        Merged intervals.
    """
    sorted_intervals = sorted((int(start), int(end)) for start, end in intervals)
    if not sorted_intervals:
        return []
    merged = [sorted_intervals[0]]
    for start, end in sorted_intervals[1:]:
        last_start, last_end = merged[-1]
        if start <= last_end + 1:
            merged[-1] = (last_start, max(last_end, end))
        else:
            merged.append((start, end))
    return merged


def make_candidate_lookup(*, candidates: pd.DataFrame) -> Dict[str, Dict[str, str]]:
    """
    Build transcript ID to candidate metadata lookup.

    Parameters
    ----------
    candidates
        Standardised candidate table.

    Returns
    -------
    dict[str, dict[str, str]]
        Candidate metadata keyed by transcript ID.
    """
    lookup: Dict[str, Dict[str, str]] = {}
    for _, row in candidates.iterrows():
        metadata = {
            "candidate_rank_tier": str(row.get("candidate_rank_tier", "")),
            "target_median_tpm": str(row.get("target_median_tpm", "")),
            "target_median_isoform_usage": str(row.get("target_median_isoform_usage", "")),
            "log2_target_vs_max_non_target_isoform_usage": str(row.get("log2_target_vs_max_non_target_isoform_usage", "")),
            "target_usage_rank_within_gene": str(row.get("target_usage_rank_within_gene", "")),
        }
        ids = {
            str(row.get("transcript_id", "")).strip(),
            str(row.get("transcript_id_with_version", "")).strip(),
            strip_ensembl_version(identifier=row.get("transcript_id_with_version", "")),
        }
        for transcript_id in ids:
            if transcript_id:
                lookup[transcript_id] = metadata
    return lookup


def value_to_float(*, value: object) -> Optional[float]:
    """
    Convert value to float where possible.

    Parameters
    ----------
    value
        Input value.

    Returns
    -------
    Optional[float]
        Float or None.
    """
    try:
        if value is None or str(value).strip() == "":
            return None
        numeric_value = float(value)
    except ValueError:
        return None
    if math.isnan(numeric_value):
        return None
    return numeric_value


def candidate_level(*, metadata: Optional[Mapping[str, str]]) -> str:
    """
    Classify candidate metadata for plotting.

    Parameters
    ----------
    metadata
        Candidate metadata or None.

    Returns
    -------
    str
        Candidate level.
    """
    if not metadata:
        return "other"
    tier = str(metadata.get("candidate_rank_tier", "")).lower()
    if "primary" in tier:
        return "primary"
    if "secondary" in tier or "2_to_3" in tier:
        return "secondary"
    return "candidate"


def candidate_metric_label(*, metadata: Mapping[str, str]) -> str:
    """
    Format candidate metrics for the right-hand labels.

    Parameters
    ----------
    metadata
        Candidate metadata.

    Returns
    -------
    str
        Compact metric label.
    """
    parts = []
    tpm = value_to_float(value=metadata.get("target_median_tpm"))
    usage = value_to_float(value=metadata.get("target_median_isoform_usage"))
    ratio = value_to_float(
        value=metadata.get("log2_target_vs_max_non_target_isoform_usage")
    )
    if tpm is not None:
        parts.append(f"TPM {tpm:.2f}")
    if usage is not None:
        parts.append(f"usage {usage:.3f}")
    if ratio is not None:
        parts.append(f"log2 ratio {ratio:.2f}")
    tier = str(metadata.get("candidate_rank_tier", "")).strip()
    if tier:
        parts.append(tier)
    return " | ".join(parts)


def read_gene_list(*, genes: Sequence[str], gene_file: Optional[Path]) -> Tuple[str, ...]:
    """
    Read and combine requested genes.

    Parameters
    ----------
    genes
        Genes from command line.
    gene_file
        Optional text file containing one gene per line.

    Returns
    -------
    tuple[str, ...]
        Unique ordered gene symbols.
    """
    output: List[str] = []
    seen = set()
    for gene in genes:
        clean = str(gene).strip()
        if clean and clean not in seen:
            output.append(clean)
            seen.add(clean)
    if gene_file is not None:
        with gene_file.open("rt") as handle:
            for line in handle:
                clean = line.strip()
                if clean and not clean.startswith("#") and clean not in seen:
                    output.append(clean)
                    seen.add(clean)
    return tuple(output)


def transcript_sort_key(
    *, transcript_id: str, rows: pd.DataFrame, lookup: Mapping[str, Mapping[str, str]]
) -> Tuple[int, float, int, str]:
    """
    Return a sort key that prioritises candidate transcripts.

    Parameters
    ----------
    transcript_id
        Transcript ID.
    rows
        Feature rows for the transcript.
    lookup
        Candidate metadata lookup.

    Returns
    -------
    tuple[int, float, int, str]
        Sort key.
    """
    metadata = lookup.get(transcript_id)
    level = candidate_level(metadata=metadata)
    priority = {"primary": 0, "secondary": 1, "candidate": 2}.get(level, 3)
    rank = float("inf")
    if metadata is not None:
        parsed_rank = value_to_float(value=metadata.get("target_usage_rank_within_gene"))
        if parsed_rank is not None:
            rank = parsed_rank
    return (priority, rank, int(rows["start"].min()), transcript_id)


def choose_transcripts(
    *, exon_features: pd.DataFrame, lookup: Mapping[str, Mapping[str, str]], max_transcripts: int
) -> List[str]:
    """
    Choose transcript tracks for plotting.

    Parameters
    ----------
    exon_features
        Exon rows for one gene.
    lookup
        Candidate metadata lookup.
    max_transcripts
        Maximum number of tracks.

    Returns
    -------
    list[str]
        Ordered transcript IDs.
    """
    transcript_ids = sorted(set(exon_features["transcript_id"]))
    sorted_ids = sorted(
        transcript_ids,
        key=lambda transcript_id: transcript_sort_key(
            transcript_id=transcript_id,
            rows=exon_features[exon_features["transcript_id"] == transcript_id],
            lookup=lookup,
        ),
    )
    return sorted_ids[:max_transcripts]


def draw_boxes(
    *, axis, intervals: Sequence[Tuple[int, int]], y: float, height: float,
    facecolour: str, edgecolour: str, linewidth: float, alpha: float = 1.0,
    zorder: int = 3,
) -> None:
    """
    Draw interval boxes.

    Parameters
    ----------
    axis
        Matplotlib axis.
    intervals
        Genomic intervals.
    y
        Track y-position.
    height
        Box height.
    facecolour
        Box fill colour.
    edgecolour
        Box edge colour.
    linewidth
        Edge linewidth.
    alpha
        Box opacity.
    zorder
        Matplotlib z-order.
    """
    for start, end in intervals:
        axis.add_patch(
            Rectangle(
                xy=(start, y - height / 2),
                width=max(end - start, 1),
                height=height,
                facecolor=facecolour,
                edgecolor=edgecolour,
                linewidth=linewidth,
                alpha=alpha,
                zorder=zorder,
            )
        )


def get_gene_features_for_plot(
    *, gene_symbol: str, features: pd.DataFrame, candidates: pd.DataFrame,
    logger: logging.Logger,
) -> pd.DataFrame:
    """
    Resolve feature rows for a requested gene symbol.

    The GENCODE feature table may not contain a gene-symbol column. When that
    happens, this function uses the candidate table to identify the gene ID or
    candidate transcript IDs for the requested gene, then returns all feature
    rows from the same gene locus.

    Parameters
    ----------
    gene_symbol
        Requested gene symbol.
    features
        Standardised feature table.
    candidates
        Standardised candidate table.
    logger
        Logger instance.

    Returns
    -------
    pandas.DataFrame
        Feature rows for the requested gene.
    """
    direct_matches = features[features["gene_symbol"] == gene_symbol].copy()
    if not direct_matches.empty:
        return direct_matches

    gene_candidates = candidates[candidates["gene_symbol"] == gene_symbol].copy()
    if gene_candidates.empty:
        logger.warning(
            "No candidate rows found for %s; cannot recover features by gene ID",
            gene_symbol,
        )
        return direct_matches

    gene_ids = {
        str(value).strip()
        for value in gene_candidates.get("gene_id", pd.Series(dtype=str))
        if str(value).strip()
    }

    candidate_transcript_ids = set()
    for column in ("transcript_id", "transcript_id_with_version"):
        if column in gene_candidates.columns:
            candidate_transcript_ids.update(
                str(value).strip()
                for value in gene_candidates[column]
                if str(value).strip()
            )
            candidate_transcript_ids.update(
                strip_ensembl_version(identifier=value)
                for value in gene_candidates[column]
                if str(value).strip()
            )

    if candidate_transcript_ids:
        matched_candidate_features = features[
            features["transcript_id"].isin(candidate_transcript_ids)
            | features["transcript_id_with_version"].isin(candidate_transcript_ids)
        ]
        gene_ids.update(
            str(value).strip()
            for value in matched_candidate_features["gene_id"]
            if str(value).strip()
        )

    if not gene_ids:
        logger.warning(
            "Could not infer a gene ID for %s from candidate transcripts",
            gene_symbol,
        )
        return direct_matches

    resolved = features[features["gene_id"].isin(gene_ids)].copy()
    if not resolved.empty:
        resolved["gene_symbol"] = gene_symbol
        logger.info(
            "Recovered %d feature rows for %s using gene IDs: %s",
            resolved.shape[0],
            gene_symbol,
            ", ".join(sorted(gene_ids)),
        )
    return resolved


def plot_gene_model(
    *, gene_symbol: str, features: pd.DataFrame, candidates: pd.DataFrame,
    out_dir: Path, output_formats: Sequence[str], out_suffix: str,
    max_transcripts: int, label_union_exons: bool, title_suffix: str,
    dpi: int, logger: logging.Logger,
) -> None:
    """
    Plot one enhanced gene model.

    Parameters
    ----------
    gene_symbol
        Gene symbol to plot.
    features
        Standardised feature table.
    candidates
        Standardised candidate table.
    out_dir
        Output directory.
    output_formats
        Output formats.
    out_suffix
        Filename suffix.
    max_transcripts
        Maximum number of transcript tracks.
    label_union_exons
        Whether to label merged union exons.
    title_suffix
        Secondary title line.
    dpi
        Output DPI for raster formats.
    logger
        Logger instance.
    """
    gene_features = get_gene_features_for_plot(
        gene_symbol=gene_symbol,
        features=features,
        candidates=candidates,
        logger=logger,
    )
    if gene_features.empty:
        logger.warning("No features found for %s", gene_symbol)
        return
    exon_features = gene_features[gene_features["feature"] == "exon"].copy()
    cds_features = gene_features[gene_features["feature"] == "cds"].copy()
    if exon_features.empty:
        logger.warning("No exon features found for %s", gene_symbol)
        return

    gene_candidates = candidates[candidates["gene_symbol"] == gene_symbol].copy()
    lookup = make_candidate_lookup(candidates=gene_candidates)
    transcript_ids = choose_transcripts(
        exon_features=exon_features, lookup=lookup, max_transcripts=max_transcripts
    )

    locus_start = int(exon_features["start"].min())
    locus_end = int(exon_features["end"].max())
    locus_span = max(locus_end - locus_start, 1)
    padding = locus_span * 0.02
    seqname = sorted(set(exon_features["seqname"]))[0]
    strand = sorted(set(exon_features["strand"]))[0]
    merged_exons = merge_intervals(
        intervals=zip(exon_features["start"], exon_features["end"])
    )

    n_tracks = len(transcript_ids)
    fig_height = max(6.5, 2.0 + 0.46 * n_tracks)
    fig, ax = plt.subplots(figsize=(17.0, fig_height))
    left_transform = blended_transform_factory(ax.transAxes, ax.transData)
    right_transform = blended_transform_factory(ax.transAxes, ax.transData)

    gene_y = n_tracks + 1.5
    ax.hlines(gene_y, locus_start, locus_end, color="black", linewidth=1.1, zorder=1)
    draw_boxes(
        axis=ax,
        intervals=merged_exons,
        y=gene_y,
        height=0.34,
        facecolour="black",
        edgecolour="black",
        linewidth=0.7,
        zorder=4,
    )
    ax.text(
        -0.012,
        gene_y,
        "Collapsed gene model\nall annotated exons",
        transform=left_transform,
        ha="right",
        va="center",
        fontsize=7.6,
        fontweight="bold",
    )
    if label_union_exons and len(merged_exons) <= 45:
        numbers = range(1, len(merged_exons) + 1)
        if strand == "-":
            numbers = range(len(merged_exons), 0, -1)
        for (start, end), number in zip(merged_exons, numbers):
            ax.text(
                start + (end - start) / 2,
                gene_y + 0.28,
                str(number),
                ha="center",
                va="bottom",
                fontsize=5.5,
                color="black",
                rotation=90 if len(merged_exons) > 25 else 0,
            )

    ax.hlines(
        n_tracks + 0.83,
        locus_start,
        locus_end,
        color="black",
        linewidth=0.5,
        alpha=0.35,
    )

    colours = {
        "primary": "#9e2f2f",
        "secondary": "#a66a1f",
        "candidate": "#4f6f99",
        "other": "#9aa3aa",
    }

    for index, transcript_id in enumerate(transcript_ids):
        y = n_tracks - index
        tx_exons = exon_features[exon_features["transcript_id"] == transcript_id].sort_values(["start", "end"])
        tx_cds = cds_features[cds_features["transcript_id"] == transcript_id].sort_values(["start", "end"])
        metadata = lookup.get(transcript_id)
        level = candidate_level(metadata=metadata)
        colour = colours[level]
        is_other = level == "other"
        alpha = 0.45 if is_other else 1.0
        intron_width = 0.65 if is_other else 1.2
        exon_height = 0.10 if is_other else 0.14
        cds_height = 0.22 if is_other else 0.32

        tx_start = int(tx_exons["start"].min())
        tx_end = int(tx_exons["end"].max())
        ax.hlines(y, tx_start, tx_end, color=colour, linewidth=intron_width, alpha=alpha)
        draw_boxes(
            axis=ax,
            intervals=list(zip(tx_exons["start"], tx_exons["end"])),
            y=y,
            height=exon_height,
            facecolour="white" if not is_other else "#eef1f2",
            edgecolour=colour,
            linewidth=0.55 if is_other else 0.9,
            alpha=alpha,
        )
        if not tx_cds.empty:
            draw_boxes(
                axis=ax,
                intervals=list(zip(tx_cds["start"], tx_cds["end"])),
                y=y,
                height=cds_height,
                facecolour=colour,
                edgecolour=colour,
                linewidth=0.5,
                alpha=alpha,
                zorder=5,
            )

        names = [x for x in tx_exons["transcript_name"].dropna().astype(str).unique() if x]
        prefix = "P: " if level == "primary" else "S: " if level == "secondary" else ""
        label = f"{prefix}{transcript_id}"
        if names:
            label += f"\n{names[0]}"
        ax.text(
            -0.012,
            y,
            label,
            transform=left_transform,
            ha="right",
            va="center",
            fontsize=6.7 if is_other else 7.4,
            color=colour,
            fontweight="bold" if not is_other else "normal",
        )
        if metadata:
            ax.text(
                1.012,
                y,
                candidate_metric_label(metadata=metadata),
                transform=right_transform,
                ha="left",
                va="center",
                fontsize=6.7,
                color=colour,
                bbox={
                    "boxstyle": "round,pad=0.2",
                    "facecolor": "white",
                    "edgecolor": colour,
                    "linewidth": 0.7,
                    "alpha": 0.95,
                },
            )

    ax.set_title(
        f"{gene_symbol}: collapsed exon catalogue and transcript isoforms\n{title_suffix}",
        loc="left",
        fontsize=14,
        fontweight="bold",
        pad=18,
    )
    ax.set_xlabel(f"Genomic coordinate on {seqname} ({strand} strand)", fontsize=9)
    ax.set_xlim(locus_start - padding, locus_end + padding)
    ax.set_ylim(0.25, n_tracks + 2.3)
    ax.set_yticks([])
    ax.spines[["left", "right", "top"]].set_visible(False)
    ax.tick_params(axis="x", labelsize=8)
    ax.grid(axis="x", linewidth=0.3, alpha=0.2)
    ax.text(
        0,
        0.02,
        "Top track: solid black boxes show the collapsed set of all annotated exon intervals. "
        "P: primary rank 1 candidate; S: secondary rank 2-3 candidate. "
        "Thick transcript boxes = CDS; thin boxes = exon/UTR; lines = introns.",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=7.2,
    )
    fig.subplots_adjust(left=0.25, right=0.76, top=0.88, bottom=0.12)

    out_dir.mkdir(parents=True, exist_ok=True)
    for output_format in output_formats:
        out_path = out_dir / f"{gene_symbol}.{out_suffix}.{output_format}"
        fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
        logger.info("Wrote %s", out_path)
    plt.close(fig)


def run(*, config: Config, logger: logging.Logger) -> None:
    """
    Run enhanced gene-model plotting.

    Parameters
    ----------
    config
        Workflow configuration.
    logger
        Logger instance.
    """
    features = normalise_features(dataframe=read_tsv(path=config.features_tsv, logger=logger))
    candidates = normalise_candidates(dataframe=read_tsv(path=config.candidate_tsv, logger=logger))
    genes = read_gene_list(genes=config.genes, gene_file=config.gene_file)
    if not genes:
        genes = tuple(sorted(candidates["gene_symbol"].dropna().unique()))
    logger.info("Plotting %d genes", len(genes))
    for gene_symbol in genes:
        plot_gene_model(
            gene_symbol=gene_symbol,
            features=features,
            candidates=candidates,
            out_dir=config.out_dir,
            output_formats=config.output_formats,
            out_suffix=config.out_suffix,
            max_transcripts=config.max_transcripts,
            label_union_exons=config.label_union_exons,
            title_suffix=config.title_suffix,
            dpi=config.dpi,
            logger=logger,
        )


def main() -> None:
    """Run command-line entry point."""
    config = parse_args()
    logger = setup_logging(log_level=config.log_level, log_path=config.log_path)
    logger.info("Enhanced candidate isoform gene-model plotting starting")
    run(config=config, logger=logger)


if __name__ == "__main__":
    main()
