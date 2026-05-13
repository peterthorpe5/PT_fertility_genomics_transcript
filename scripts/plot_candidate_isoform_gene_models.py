#!/usr/bin/env python3
"""
Plot publication-quality transcript model figures for candidate isoforms.

This script creates genome-browser-style gene model figures from the GENCODE
feature table produced by annotate_gtex_isoform_screen_with_gencode.py and the
annotated GTEx candidate isoform table.

The figure highlights candidate testis-preferential isoforms while still
showing the other annotated transcripts for the gene. Exons are drawn as boxes,
CDS segments as thicker boxes, and introns as connecting lines.
"""

from __future__ import annotations

import argparse
import logging
import math
import sys
import time
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Sequence, Tuple

import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


@dataclass(frozen=True)
class Config:
    """
    Configuration for candidate isoform gene-model plotting.

    Attributes
    ----------
    features_tsv
        GENCODE transcript feature table from the annotation script.
    candidates_tsv
        Annotated candidate isoform table.
    out_dir
        Output directory for figures.
    genes
        Optional list of gene symbols to plot. If omitted, top-ranked candidate
        genes are selected from the candidate table.
    max_genes
        Maximum number of genes to plot if genes are not supplied.
    max_transcripts_per_gene
        Maximum number of transcript models to show per gene.
    formats
        Output formats, for example pdf and svg.
    title_suffix
        Optional text appended to figure titles.
    log_path
        Optional log-file path.
    log_level
        Logging level.
    """

    features_tsv: Path
    candidates_tsv: Path
    out_dir: Path
    genes: Optional[Tuple[str, ...]]
    max_genes: int
    max_transcripts_per_gene: int
    formats: Tuple[str, ...]
    title_suffix: str
    log_path: Optional[Path]
    log_level: str


@contextmanager
def timed(*, logger: logging.Logger, label: str) -> Iterator[None]:
    """
    Log elapsed time for a block.

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
        Optional log-file path.

    Returns
    -------
    logging.Logger
        Configured logger.
    """
    logger = logging.getLogger("plot_candidate_isoform_gene_models")
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
        description="Plot candidate GTEx isoform gene models from GENCODE features."
    )
    parser.add_argument(
        "--features_tsv",
        required=True,
        type=Path,
        help="GENCODE transcript feature TSV/TSV.GZ from annotation script.",
    )
    parser.add_argument(
        "--candidates_tsv",
        required=True,
        type=Path,
        help="Annotated candidate isoform TSV/TSV.GZ.",
    )
    parser.add_argument(
        "--out_dir",
        required=True,
        type=Path,
        help="Output directory for figures.",
    )
    parser.add_argument(
        "--genes",
        nargs="*",
        default=None,
        help=(
            "Optional gene symbols to plot. If omitted, the highest-ranked "
            "candidate genes are selected automatically."
        ),
    )
    parser.add_argument(
        "--max_genes",
        type=int,
        default=25,
        help="Maximum number of genes to plot when --genes is omitted.",
    )
    parser.add_argument(
        "--max_transcripts_per_gene",
        type=int,
        default=20,
        help="Maximum transcript rows shown per gene. Default: 20.",
    )
    parser.add_argument(
        "--formats",
        nargs="+",
        default=("pdf", "svg"),
        help="Output formats. Default: pdf svg.",
    )
    parser.add_argument(
        "--title_suffix",
        default="GTEx testis-preferential isoform screen",
        help="Optional title suffix.",
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
    genes = tuple(args.genes) if args.genes else None
    return Config(
        features_tsv=args.features_tsv,
        candidates_tsv=args.candidates_tsv,
        out_dir=args.out_dir,
        genes=genes,
        max_genes=args.max_genes,
        max_transcripts_per_gene=args.max_transcripts_per_gene,
        formats=tuple(args.formats),
        title_suffix=args.title_suffix,
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
    dataframe = pd.read_csv(filepath_or_buffer=path, sep="\t", dtype=str)
    logger.info(
        "Read %d rows x %d columns from %s",
        dataframe.shape[0],
        dataframe.shape[1],
        path,
    )
    return dataframe


def to_float(*, value: object, default: float = 0.0) -> float:
    """
    Convert a value to float with fallback.

    Parameters
    ----------
    value
        Value to convert.
    default
        Returned if conversion fails.

    Returns
    -------
    float
        Converted value.
    """
    try:
        if pd.isna(value):
            return default
        return float(value)
    except (TypeError, ValueError):
        return default


def normalise_gene_name_columns(*, candidates: pd.DataFrame) -> pd.DataFrame:
    """
    Add stable gene and transcript display columns to the candidate table.

    Parameters
    ----------
    candidates
        Annotated candidate table.

    Returns
    -------
    pandas.DataFrame
        Candidate table with helper display columns.
    """
    output = candidates.copy()
    if "gene_symbol" in output.columns:
        output["plot_gene_symbol"] = output["gene_symbol"].fillna("")
    elif "gencode_gene_name" in output.columns:
        output["plot_gene_symbol"] = output["gencode_gene_name"].fillna("")
    else:
        raise ValueError("Candidate table needs gene_symbol or gencode_gene_name")

    if "transcript_id_with_version" in output.columns:
        output["plot_transcript_id"] = output["transcript_id_with_version"].fillna("")
    elif "gencode_transcript_id_with_version" in output.columns:
        output["plot_transcript_id"] = output[
            "gencode_transcript_id_with_version"
        ].fillna("")
    elif "transcript_id" in output.columns:
        output["plot_transcript_id"] = output["transcript_id"].fillna("")
    else:
        raise ValueError("Candidate table needs a transcript ID column")
    return output


def select_genes_to_plot(*, candidates: pd.DataFrame, config: Config) -> List[str]:
    """
    Select genes for plotting.

    Parameters
    ----------
    candidates
        Annotated candidate table with helper columns.
    config
        Plotting configuration.

    Returns
    -------
    List[str]
        Gene symbols to plot.
    """
    if config.genes:
        return list(config.genes)

    ranking_cols = [
        "is_broad_gene_isoform_rescue_candidate",
        "is_target_tissue_isoform_candidate",
        "target_median_isoform_usage",
        "log2_target_vs_max_non_target_isoform_usage",
        "target_median_tpm",
    ]
    available = [col for col in ranking_cols if col in candidates.columns]
    ranked = candidates.copy()
    for column in available:
        ranked[column] = ranked[column].map(lambda value: to_float(value=value))
    if available:
        ranked = ranked.sort_values(by=available, ascending=[False] * len(available))
    genes = ranked["plot_gene_symbol"].dropna().astype(str)
    genes = [gene for gene in genes if gene]
    return list(dict.fromkeys(genes))[: config.max_genes]


def prepare_gene_features(*, features: pd.DataFrame, gene_symbol: str) -> pd.DataFrame:
    """
    Subset and type-convert GENCODE features for one gene.

    Parameters
    ----------
    features
        GENCODE feature table.
    gene_symbol
        Gene symbol.

    Returns
    -------
    pandas.DataFrame
        Feature subset for the gene.
    """
    subset = features.loc[features["gencode_gene_name"].eq(gene_symbol)].copy()
    if subset.empty:
        return subset
    subset["gencode_start"] = pd.to_numeric(subset["gencode_start"], errors="coerce")
    subset["gencode_end"] = pd.to_numeric(subset["gencode_end"], errors="coerce")
    subset = subset.dropna(subset=["gencode_start", "gencode_end"])
    subset["gencode_start"] = subset["gencode_start"].astype(int)
    subset["gencode_end"] = subset["gencode_end"].astype(int)
    return subset


def rank_transcripts_for_gene(
    *,
    gene_features: pd.DataFrame,
    gene_candidates: pd.DataFrame,
    max_transcripts: int,
) -> List[str]:
    """
    Rank transcripts for display in a gene-model plot.

    Candidate transcripts are shown first, followed by non-candidate
    transcripts with more exon evidence.

    Parameters
    ----------
    gene_features
        GENCODE features for one gene.
    gene_candidates
        Candidate rows for one gene.
    max_transcripts
        Maximum transcript IDs to return.

    Returns
    -------
    List[str]
        Ordered transcript IDs with versions.
    """
    candidate_ids = list(
        dict.fromkeys(
            gene_candidates["plot_transcript_id"].dropna().astype(str).tolist()
        )
    )
    feature_counts = (
        gene_features.groupby("gencode_transcript_id_with_version")[
            "gencode_feature_type"
        ]
        .count()
        .sort_values(ascending=False)
    )
    non_candidate_ids = [
        transcript_id
        for transcript_id in feature_counts.index.astype(str).tolist()
        if transcript_id not in candidate_ids
    ]
    ordered = candidate_ids + non_candidate_ids
    return ordered[:max_transcripts]


def make_candidate_lookup(
    *,
    gene_candidates: pd.DataFrame,
) -> Dict[str, Dict[str, str]]:
    """
    Make a transcript ID to candidate metadata lookup.

    Parameters
    ----------
    gene_candidates
        Candidate rows for one gene.

    Returns
    -------
    Dict[str, Dict[str, str]]
        Candidate metadata keyed by transcript ID with version.
    """
    lookup: Dict[str, Dict[str, str]] = {}
    for _, row in gene_candidates.iterrows():
        transcript_id = str(row.get("plot_transcript_id", ""))
        if not transcript_id:
            continue
        lookup[transcript_id] = {
            column: str(row.get(column, ""))
            for column in row.index
        }
    return lookup


def get_track_style(*, transcript_id: str, candidate_lookup: Dict[str, Dict[str, str]]):
    """
    Return drawing style for a transcript track.

    Parameters
    ----------
    transcript_id
        Transcript ID with version.
    candidate_lookup
        Candidate metadata lookup.

    Returns
    -------
    tuple
        Face colour, edge colour, line width, alpha and label prefix.
    """
    if transcript_id not in candidate_lookup:
        return "#D8DDE6", "#7A8391", 0.8, 0.65, ""

    tier = candidate_lookup[transcript_id].get("candidate_rank_tier", "")
    if tier == "primary_rank_1":
        return "#D84A4A", "#8E1F1F", 1.8, 0.95, "P: "
    if tier == "secondary_rank_2_to_3":
        return "#F2A14A", "#A95C00", 1.5, 0.9, "S: "
    return "#7B68EE", "#483D8B", 1.4, 0.9, "C: "


def format_metric(*, value: str, digits: int = 3) -> str:
    """
    Format a numeric candidate metric for display.

    Parameters
    ----------
    value
        Metric value.
    digits
        Number of significant digits.

    Returns
    -------
    str
        Formatted metric.
    """
    try:
        return f"{float(value):.{digits}g}"
    except (TypeError, ValueError):
        return ""


def draw_gene_model(
    *,
    gene_symbol: str,
    gene_features: pd.DataFrame,
    gene_candidates: pd.DataFrame,
    transcript_order: List[str],
    output_paths: Sequence[Path],
    title_suffix: str,
    logger: logging.Logger,
) -> None:
    """
    Draw and save a transcript model figure for one gene.

    Parameters
    ----------
    gene_symbol
        Gene symbol.
    gene_features
        GENCODE features for the gene.
    gene_candidates
        Candidate rows for the gene.
    transcript_order
        Ordered transcript IDs to display.
    output_paths
        Output figure paths.
    title_suffix
        Text appended to the figure title.
    logger
        Logger instance.
    """
    if not transcript_order:
        logger.warning("No transcripts to plot for %s", gene_symbol)
        return

    gene_features = gene_features.loc[
        gene_features["gencode_transcript_id_with_version"].isin(transcript_order)
    ].copy()
    candidate_lookup = make_candidate_lookup(gene_candidates=gene_candidates)

    min_pos = int(gene_features["gencode_start"].min())
    max_pos = int(gene_features["gencode_end"].max())
    seqname = str(gene_features["gencode_seqname"].dropna().iloc[0])
    strand = str(gene_features["gencode_strand"].dropna().iloc[0])
    span = max_pos - min_pos + 1

    n_tracks = len(transcript_order)
    figure_height = max(4.5, 1.0 + 0.42 * n_tracks)
    figure_width = 15
    fig, ax = plt.subplots(figsize=(figure_width, figure_height))

    ax.set_xlim(min_pos - span * 0.03, max_pos + span * 0.28)
    ax.set_ylim(-0.8, n_tracks + 1.2)
    ax.set_yticks([])
    ax.spines[["left", "right", "top"]].set_visible(False)
    ax.ticklabel_format(style="plain", axis="x")
    ax.set_xlabel(f"Genomic coordinate on {seqname} ({strand} strand)")

    title = f"{gene_symbol}: transcript isoform structure"
    if title_suffix:
        title = f"{title}\n{title_suffix}"
    ax.set_title(title, loc="left", fontsize=14, fontweight="bold")

    row_lookup = {
        transcript_id: n_tracks - index
        for index, transcript_id in enumerate(transcript_order)
    }

    for transcript_id in transcript_order:
        y_pos = row_lookup[transcript_id]
        transcript_features = gene_features.loc[
            gene_features["gencode_transcript_id_with_version"].eq(transcript_id)
        ]
        exon_features = transcript_features.loc[
            transcript_features["gencode_feature_type"].eq("exon")
        ].sort_values("gencode_start")
        cds_features = transcript_features.loc[
            transcript_features["gencode_feature_type"].eq("CDS")
        ].sort_values("gencode_start")
        if exon_features.empty:
            continue

        face_colour, edge_colour, line_width, alpha, label_prefix = get_track_style(
            transcript_id=transcript_id,
            candidate_lookup=candidate_lookup,
        )
        tx_start = int(exon_features["gencode_start"].min())
        tx_end = int(exon_features["gencode_end"].max())
        ax.plot(
            [tx_start, tx_end],
            [y_pos, y_pos],
            color=edge_colour,
            linewidth=1.0,
            alpha=0.75,
            zorder=1,
        )

        for _, exon in exon_features.iterrows():
            start = int(exon["gencode_start"])
            end = int(exon["gencode_end"])
            width = max(end - start + 1, span * 0.0008)
            rectangle = patches.Rectangle(
                (start, y_pos - 0.10),
                width,
                0.20,
                facecolor=face_colour,
                edgecolor=edge_colour,
                linewidth=line_width,
                alpha=alpha,
                zorder=2,
            )
            ax.add_patch(rectangle)

        for _, cds in cds_features.iterrows():
            start = int(cds["gencode_start"])
            end = int(cds["gencode_end"])
            width = max(end - start + 1, span * 0.0008)
            rectangle = patches.Rectangle(
                (start, y_pos - 0.18),
                width,
                0.36,
                facecolor=face_colour,
                edgecolor=edge_colour,
                linewidth=line_width + 0.4,
                alpha=min(alpha + 0.05, 1.0),
                zorder=3,
            )
            ax.add_patch(rectangle)

        transcript_name = str(
            exon_features["gencode_transcript_name"].replace("nan", "").iloc[0]
        )
        label = f"{label_prefix}{transcript_id}"
        if transcript_name and transcript_name != "nan":
            label = f"{label}\n{transcript_name}"
        ax.text(
            min_pos - span * 0.02,
            y_pos,
            label,
            ha="right",
            va="center",
            fontsize=7.5,
            color=edge_colour,
        )

        if transcript_id in candidate_lookup:
            row = candidate_lookup[transcript_id]
            tpm_text = format_metric(
                value=row.get("target_median_tpm", ""),
            )
            usage_text = format_metric(
                value=row.get("target_median_isoform_usage", ""),
            )
            ratio_text = format_metric(
                value=row.get(
                    "log2_target_vs_max_non_target_isoform_usage",
                    "",
                ),
            )
            tier_text = row.get("candidate_rank_tier", "")
            metric_text = (
                f"TPM {tpm_text} | usage {usage_text} | "
                f"log2 ratio {ratio_text} | {tier_text}"
            )
            ax.text(
                max_pos + span * 0.03,
                y_pos,
                metric_text,
                ha="left",
                va="center",
                fontsize=7.5,
                color=edge_colour,
                bbox={
                    "boxstyle": "round,pad=0.25",
                    "facecolor": "white",
                    "edgecolor": edge_colour,
                    "alpha": 0.9,
                },
            )

    legend_y = 0.2
    ax.text(
        min_pos,
        legend_y,
        "P: primary rank 1 candidate    S: secondary rank 2-3 candidate    "
        "pale tracks: other annotated transcripts\n"
        "Thick boxes = CDS; thin boxes = exon/UTR; connecting lines = introns",
        fontsize=8,
        va="bottom",
        ha="left",
        color="#333333",
    )

    fig.tight_layout()
    for output_path in output_paths:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, bbox_inches="tight")
        logger.info("Wrote figure to %s", output_path)
    plt.close(fig)


def safe_filename(*, value: str) -> str:
    """
    Convert a string to a safe filename component.

    Parameters
    ----------
    value
        Input string.

    Returns
    -------
    str
        Safe filename component.
    """
    safe = "".join(
        char if char.isalnum() or char in {"_", "-"} else "_"
        for char in value
    )
    return safe.strip("_") or "gene"


def run(*, config: Config, logger: logging.Logger) -> None:
    """
    Run candidate isoform gene-model plotting.

    Parameters
    ----------
    config
        Plotting configuration.
    logger
        Logger instance.
    """
    config.out_dir.mkdir(parents=True, exist_ok=True)
    with timed(logger=logger, label="Read feature and candidate tables"):
        features = read_tsv(path=config.features_tsv, logger=logger)
        candidates = read_tsv(path=config.candidates_tsv, logger=logger)
        candidates = normalise_gene_name_columns(candidates=candidates)

    genes = select_genes_to_plot(candidates=candidates, config=config)
    logger.info("Selected %d genes for plotting: %s", len(genes), ", ".join(genes))

    with timed(logger=logger, label="Plot gene models"):
        for gene_symbol in genes:
            gene_features = prepare_gene_features(
                features=features,
                gene_symbol=gene_symbol,
            )
            if gene_features.empty:
                logger.warning("No GENCODE features found for gene %s", gene_symbol)
                continue
            gene_candidates = candidates.loc[
                candidates["plot_gene_symbol"].eq(gene_symbol)
            ].copy()
            transcript_order = rank_transcripts_for_gene(
                gene_features=gene_features,
                gene_candidates=gene_candidates,
                max_transcripts=config.max_transcripts_per_gene,
            )
            output_paths = [
                config.out_dir
                / (
                    f"{safe_filename(value=gene_symbol)}"
                    f".candidate_isoform_gene_model.{fmt}"
                )
                for fmt in config.formats
            ]
            draw_gene_model(
                gene_symbol=gene_symbol,
                gene_features=gene_features,
                gene_candidates=gene_candidates,
                transcript_order=transcript_order,
                output_paths=output_paths,
                title_suffix=config.title_suffix,
                logger=logger,
            )

    logger.info("Candidate isoform gene-model plotting finished")


def main() -> None:
    """Run the command-line entry point."""
    config = parse_args()
    logger = setup_logging(log_level=config.log_level, log_path=config.log_path)
    logger.info("Candidate isoform gene-model plotting starting")
    logger.info("Feature table: %s", config.features_tsv)
    logger.info("Candidate table: %s", config.candidates_tsv)
    logger.info("Output directory: %s", config.out_dir)
    run(config=config, logger=logger)


if __name__ == "__main__":
    main()
