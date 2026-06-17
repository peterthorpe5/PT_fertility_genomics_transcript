#!/usr/bin/env python3
"""
Classify coding-sequence consequences of candidate testis isoforms.

The transcript-level prioritisation identifies candidate testis-preferential
isoforms. For manuscript interpretation, the important next question is whether
those isoforms change coding sequence or only non-coding/UTR structure. This
script compares each candidate transcript to a reference transcript from the
same gene using GENCODE exon and CDS features.

Reference transcript selection is deliberately explicit:
1. MANE Select protein-coding transcript with CDS, if available.
2. Ensembl canonical protein-coding transcript with CDS, if available.
3. Basic protein-coding transcript with the longest CDS, if available.
4. Any transcript with the longest CDS.

If the selected reference is the candidate itself, a secondary non-candidate
reference is also selected where possible and used for the final classification.
"""

from __future__ import annotations

import argparse
import logging
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import pandas as pd


LOGGER_NAME = "classify_candidate_isoform_cds_consequences"


@dataclass(frozen=True)
class Config:
    """Runtime configuration."""

    candidate_table: Path
    transcript_annotation: Path
    transcript_features: Path
    out_dir: Path
    out_prefix: str
    write_excel_outputs: bool
    log_level: str
    log_path: Optional[Path]


def setup_logging(*, log_level: str, log_path: Optional[Path]) -> logging.Logger:
    """Configure and return a logger."""
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
        description="Classify candidate isoforms by CDS/exon structural change."
    )
    parser.add_argument("--candidate_table", required=True, type=Path)
    parser.add_argument("--transcript_annotation", required=True, type=Path)
    parser.add_argument("--transcript_features", required=True, type=Path)
    parser.add_argument("--out_dir", required=True, type=Path)
    parser.add_argument(
        "--out_prefix",
        default="gtex_v11_candidate_isoform_cds_consequence",
    )
    parser.add_argument("--no_write_excel_outputs", action="store_true")
    parser.add_argument("--log_level", default="INFO")
    parser.add_argument("--log_path", type=Path, default=None)
    args = parser.parse_args(argv)
    return Config(
        candidate_table=args.candidate_table,
        transcript_annotation=args.transcript_annotation,
        transcript_features=args.transcript_features,
        out_dir=args.out_dir,
        out_prefix=args.out_prefix,
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
    logger.info("Loaded %d rows x %d columns from %s", *dataframe.shape, path)
    return dataframe


def as_bool(series: pd.Series) -> pd.Series:
    """Convert common bool-like values to Boolean."""
    return series.astype(str).str.strip().str.lower().isin({"true", "1", "yes", "y"})


def get_first_column(*, columns: Iterable[str], candidates: Sequence[str]) -> str:
    """Return the first available column from a list of candidates."""
    available = set(columns)
    for column in candidates:
        if column in available:
            return column
    raise ValueError(f"Missing required columns. Tried: {candidates}")


def normalise_candidate_columns(*, dataframe: pd.DataFrame) -> pd.DataFrame:
    """Normalise candidate table identifiers to common column names."""
    output = dataframe.copy()
    transcript_col = get_first_column(
        columns=output.columns,
        candidates=("transcript_id_with_version", "gencode_transcript_id_with_version"),
    )
    gene_col = get_first_column(
        columns=output.columns,
        candidates=("gene_symbol", "gencode_gene_name", "gene_name"),
    )
    output["candidate_transcript_id_with_version"] = output[transcript_col].astype(str)
    output["candidate_gene_symbol"] = output[gene_col].astype(str)
    if "gencode_gene_id" in output.columns:
        output["candidate_gene_id"] = output["gencode_gene_id"].astype(str)
    elif "gene_id" in output.columns:
        output["candidate_gene_id"] = output["gene_id"].astype(str)
    else:
        output["candidate_gene_id"] = ""
    return output


def strip_version(identifier: object) -> str:
    """Strip an Ensembl dot-version suffix from an identifier."""
    text = str(identifier)
    return text.split(".")[0]


def interval_length(intervals: Iterable[Tuple[int, int]]) -> int:
    """Return total inclusive length of intervals."""
    return int(sum((end - start + 1) for start, end in intervals))


def merge_intervals(intervals: Iterable[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """Merge overlapping or adjacent intervals."""
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


def interval_intersection_length(
    *, first: Sequence[Tuple[int, int]], second: Sequence[Tuple[int, int]]
) -> int:
    """Calculate inclusive intersection length between merged intervals."""
    i = 0
    j = 0
    total = 0
    first = list(first)
    second = list(second)
    while i < len(first) and j < len(second):
        start = max(first[i][0], second[j][0])
        end = min(first[i][1], second[j][1])
        if start <= end:
            total += end - start + 1
        if first[i][1] < second[j][1]:
            i += 1
        else:
            j += 1
    return int(total)


def normalise_annotation(*, annotation: pd.DataFrame) -> pd.DataFrame:
    """Normalise transcript annotation columns and flags."""
    output = annotation.copy()
    output["transcript_id_with_version"] = output[
        get_first_column(
            columns=output.columns,
            candidates=("gencode_transcript_id_with_version", "transcript_id_with_version"),
        )
    ].astype(str)
    output["gene_symbol"] = output[
        get_first_column(columns=output.columns, candidates=("gencode_gene_name", "gene_symbol"))
    ].astype(str)
    output["gene_id"] = output[
        get_first_column(columns=output.columns, candidates=("gencode_gene_id", "gene_id"))
    ].astype(str)
    for column in ("is_mane_select", "is_ensembl_canonical", "is_basic", "gencode_has_cds"):
        if column not in output.columns:
            output[column] = 0
    if "gencode_cds_length_bp" not in output.columns:
        output["gencode_cds_length_bp"] = 0
    if "gencode_is_protein_coding_transcript" not in output.columns:
        output["gencode_is_protein_coding_transcript"] = 0
    output["gencode_cds_length_bp"] = pd.to_numeric(
        output["gencode_cds_length_bp"], errors="coerce"
    ).fillna(0)
    return output


def normalise_features(*, features: pd.DataFrame) -> pd.DataFrame:
    """Normalise feature columns."""
    output = features.copy()
    output["transcript_id_with_version"] = output[
        get_first_column(
            columns=output.columns,
            candidates=("gencode_transcript_id_with_version", "transcript_id_with_version"),
        )
    ].astype(str)
    output["gene_symbol"] = output[
        get_first_column(columns=output.columns, candidates=("gencode_gene_name", "gene_symbol"))
    ].astype(str)
    output["gene_id"] = output[
        get_first_column(columns=output.columns, candidates=("gencode_gene_id", "gene_id"))
    ].astype(str)
    output["feature_type"] = output[
        get_first_column(columns=output.columns, candidates=("gencode_feature_type", "feature_type"))
    ].astype(str)
    output["start"] = pd.to_numeric(
        output[get_first_column(columns=output.columns, candidates=("gencode_start", "start"))],
        errors="coerce",
    )
    output["end"] = pd.to_numeric(
        output[get_first_column(columns=output.columns, candidates=("gencode_end", "end"))],
        errors="coerce",
    )
    output = output.dropna(subset=["start", "end"])
    output["start"] = output["start"].astype(int)
    output["end"] = output["end"].astype(int)
    return output


def build_interval_lookup(*, features: pd.DataFrame, feature_type: str) -> Dict[str, List[Tuple[int, int]]]:
    """Build transcript-to-interval lookup for one feature type."""
    subset = features.loc[features["feature_type"] == feature_type]
    lookup: Dict[str, List[Tuple[int, int]]] = {}
    for transcript_id, group in subset.groupby("transcript_id_with_version"):
        intervals = merge_intervals(zip(group["start"], group["end"]))
        lookup[str(transcript_id)] = intervals
    return lookup


def select_reference_transcript(
    *, gene_annotation: pd.DataFrame, candidate_transcript_id: str
) -> Tuple[str, str]:
    """Select the best available reference transcript for a candidate gene."""
    annotation = gene_annotation.copy()
    annotation["_has_cds"] = as_bool(annotation["gencode_has_cds"]) | (
        pd.to_numeric(annotation["gencode_cds_length_bp"], errors="coerce").fillna(0) > 0
    )
    coding = annotation.loc[annotation["_has_cds"]].copy()
    if coding.empty:
        return "", "no_cds_reference_available"

    mane = coding.loc[as_bool(coding["is_mane_select"])]
    if not mane.empty:
        return str(mane.iloc[0]["transcript_id_with_version"]), "mane_select"

    canonical = coding.loc[as_bool(coding["is_ensembl_canonical"])]
    if not canonical.empty:
        return str(canonical.iloc[0]["transcript_id_with_version"]), "ensembl_canonical"

    basic = coding.loc[as_bool(coding["is_basic"])]
    if not basic.empty:
        basic = basic.sort_values("gencode_cds_length_bp", ascending=False)
        return str(basic.iloc[0]["transcript_id_with_version"]), "basic_longest_cds"

    coding = coding.sort_values("gencode_cds_length_bp", ascending=False)
    return str(coding.iloc[0]["transcript_id_with_version"]), "longest_cds"


def select_secondary_reference(
    *, gene_annotation: pd.DataFrame, candidate_transcript_id: str
) -> Tuple[str, str]:
    """Select a non-candidate reference transcript where possible."""
    annotation = gene_annotation.copy()
    annotation["_has_cds"] = as_bool(annotation["gencode_has_cds"]) | (
        pd.to_numeric(annotation["gencode_cds_length_bp"], errors="coerce").fillna(0) > 0
    )
    coding = annotation.loc[
        annotation["_has_cds"]
        & (annotation["transcript_id_with_version"] != candidate_transcript_id)
    ].copy()
    if coding.empty:
        return "", "no_non_candidate_reference_available"
    coding = coding.sort_values("gencode_cds_length_bp", ascending=False)
    return str(coding.iloc[0]["transcript_id_with_version"]), "non_candidate_longest_cds"


def compare_cds_intervals(
    *, candidate_cds: Sequence[Tuple[int, int]], reference_cds: Sequence[Tuple[int, int]]
) -> Dict[str, object]:
    """Compare candidate and reference CDS interval sets."""
    candidate_cds = merge_intervals(candidate_cds)
    reference_cds = merge_intervals(reference_cds)
    candidate_len = interval_length(candidate_cds)
    reference_len = interval_length(reference_cds)
    shared_len = interval_intersection_length(first=candidate_cds, second=reference_cds)
    union_len = candidate_len + reference_len - shared_len
    candidate_specific = candidate_len - shared_len
    reference_specific = reference_len - shared_len
    jaccard = shared_len / union_len if union_len else 0.0

    if candidate_len == 0:
        cds_change_class = "candidate_without_cds"
    elif reference_len == 0:
        cds_change_class = "no_reference_cds"
    elif candidate_specific == 0 and reference_specific == 0:
        cds_change_class = "cds_identical_to_reference"
    elif candidate_specific > 0 and reference_specific == 0:
        cds_change_class = "candidate_cds_superset_of_reference"
    elif candidate_specific == 0 and reference_specific > 0:
        cds_change_class = "candidate_cds_subset_of_reference"
    else:
        cds_change_class = "cds_partial_overlap_change"

    return {
        "candidate_cds_bp": candidate_len,
        "reference_cds_bp": reference_len,
        "shared_cds_bp": shared_len,
        "candidate_specific_cds_bp": candidate_specific,
        "reference_specific_cds_bp": reference_specific,
        "cds_jaccard": round(jaccard, 6),
        "cds_change_class": cds_change_class,
        "cds_changed_relative_to_reference": cds_change_class
        not in {"cds_identical_to_reference", "candidate_without_cds", "no_reference_cds"},
    }


def classify_candidates(
    *, candidates: pd.DataFrame, annotation: pd.DataFrame, features: pd.DataFrame, logger: logging.Logger
) -> pd.DataFrame:
    """Classify CDS consequences for all candidate rows."""
    candidates = normalise_candidate_columns(dataframe=candidates)
    annotation = normalise_annotation(annotation=annotation)
    features = normalise_features(features=features)

    candidate_genes = set(candidates["candidate_gene_symbol"].astype(str))
    annotation_subset = annotation.loc[annotation["gene_symbol"].isin(candidate_genes)].copy()
    feature_subset = features.loc[features["gene_symbol"].isin(candidate_genes)].copy()
    logger.info(
        "Subset annotation to %d rows and features to %d rows for %d candidate genes",
        len(annotation_subset),
        len(feature_subset),
        len(candidate_genes),
    )

    cds_lookup = build_interval_lookup(features=feature_subset, feature_type="CDS")
    exon_lookup = build_interval_lookup(features=feature_subset, feature_type="exon")
    annotation_by_gene = dict(tuple(annotation_subset.groupby("gene_symbol")))

    rows: List[Dict[str, object]] = []
    for _, candidate in candidates.iterrows():
        gene_symbol = str(candidate["candidate_gene_symbol"])
        candidate_tx = str(candidate["candidate_transcript_id_with_version"])
        gene_annotation = annotation_by_gene.get(gene_symbol, pd.DataFrame())
        if gene_annotation.empty:
            rows.append(
                {
                    "gene_symbol": gene_symbol,
                    "transcript_id_with_version": candidate_tx,
                    "reference_transcript_id": "",
                    "reference_selection_reason": "gene_not_found_in_annotation",
                    "final_reference_transcript_id": "",
                    "final_reference_reason": "gene_not_found_in_annotation",
                    "cds_change_class": "annotation_missing",
                    "cds_changed_relative_to_reference": False,
                }
            )
            continue

        ref_tx, ref_reason = select_reference_transcript(
            gene_annotation=gene_annotation, candidate_transcript_id=candidate_tx
        )
        final_ref_tx = ref_tx
        final_ref_reason = ref_reason
        candidate_is_primary_reference = ref_tx == candidate_tx
        if candidate_is_primary_reference:
            secondary_tx, secondary_reason = select_secondary_reference(
                gene_annotation=gene_annotation, candidate_transcript_id=candidate_tx
            )
            if secondary_tx:
                final_ref_tx = secondary_tx
                final_ref_reason = secondary_reason

        comparison = compare_cds_intervals(
            candidate_cds=cds_lookup.get(candidate_tx, []),
            reference_cds=cds_lookup.get(final_ref_tx, []),
        )
        candidate_exons = exon_lookup.get(candidate_tx, [])
        reference_exons = exon_lookup.get(final_ref_tx, [])
        exon_comparison = compare_cds_intervals(
            candidate_cds=candidate_exons,
            reference_cds=reference_exons,
        )
        row = candidate.to_dict()
        row.update(
            {
                "gene_symbol": gene_symbol,
                "transcript_id_with_version": candidate_tx,
                "reference_transcript_id": ref_tx,
                "reference_selection_reason": ref_reason,
                "candidate_is_primary_reference": candidate_is_primary_reference,
                "final_reference_transcript_id": final_ref_tx,
                "final_reference_reason": final_ref_reason,
                "candidate_exon_bp": exon_comparison["candidate_cds_bp"],
                "reference_exon_bp": exon_comparison["reference_cds_bp"],
                "candidate_specific_exon_bp": exon_comparison[
                    "candidate_specific_cds_bp"
                ],
                "reference_specific_exon_bp": exon_comparison[
                    "reference_specific_cds_bp"
                ],
                "exon_jaccard": exon_comparison["cds_jaccard"],
                "exon_structure_class": exon_comparison["cds_change_class"].replace(
                    "cds", "exon"
                ),
            }
        )
        row.update(comparison)
        if row["cds_change_class"] == "cds_identical_to_reference" and row[
            "exon_structure_class"
        ] != "exon_identical_to_reference":
            row["biological_review_class"] = "utr_or_non_coding_exon_change_only"
        elif row["cds_change_class"] in {
            "candidate_cds_superset_of_reference",
            "candidate_cds_subset_of_reference",
            "cds_partial_overlap_change",
        }:
            row["biological_review_class"] = "cds_changing_candidate"
        elif row["cds_change_class"] in {"candidate_without_cds", "no_reference_cds"}:
            row["biological_review_class"] = "ambiguous_or_no_cds_reference"
        else:
            row["biological_review_class"] = "cds_identical_or_not_structurally_compelling"
        rows.append(row)

    return pd.DataFrame(rows)


def summarise_classification(*, dataframe: pd.DataFrame) -> pd.DataFrame:
    """Summarise CDS classification results."""
    group_cols = ["biological_review_class", "cds_change_class"]
    rows = []
    for keys, group in dataframe.groupby(group_cols, dropna=False):
        if not isinstance(keys, tuple):
            keys = (keys,)
        rows.append(
            {
                "biological_review_class": keys[0],
                "cds_change_class": keys[1],
                "n_rows": len(group),
                "n_unique_genes": group.get("gene_symbol", pd.Series(dtype=str)).nunique(),
                "median_candidate_cds_bp": pd.to_numeric(
                    group.get("candidate_cds_bp", pd.Series(dtype=float)), errors="coerce"
                ).median(),
                "median_candidate_specific_cds_bp": pd.to_numeric(
                    group.get("candidate_specific_cds_bp", pd.Series(dtype=float)), errors="coerce"
                ).median(),
            }
        )
    return pd.DataFrame(rows).sort_values("n_rows", ascending=False)


def write_tsv(*, dataframe: pd.DataFrame, path: Path, logger: logging.Logger) -> None:
    """Write a TSV output."""
    path.parent.mkdir(parents=True, exist_ok=True)
    dataframe.to_csv(path, sep="\t", index=False)
    logger.info("Wrote %d rows x %d columns to %s", *dataframe.shape, path)


def write_excel(*, dataframe: pd.DataFrame, path: Path, sheet_name: str, logger: logging.Logger) -> None:
    """Write a formatted Excel file."""
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
        for idx, column in enumerate(dataframe.columns):
            value_width = dataframe[column].head(1000).astype(str).str.len().max()
            if pd.isna(value_width):
                value_width = 0
            width = min(max(len(str(column)) + 2, int(value_width) + 2, 8), 45)
            fmt = float_format if pd.api.types.is_float_dtype(dataframe[column]) else None
            worksheet.set_column(idx, idx, width, fmt)
    logger.info("Wrote formatted Excel output to %s", path)


def write_outputs(*, outputs: Dict[str, pd.DataFrame], config: Config, logger: logging.Logger) -> None:
    """Write all outputs as TSV and optionally XLSX."""
    for suffix, dataframe in outputs.items():
        path = config.out_dir / f"{config.out_prefix}.{suffix}.tsv"
        write_tsv(dataframe=dataframe, path=path, logger=logger)
        if config.write_excel_outputs:
            write_excel(
                dataframe=dataframe,
                path=config.out_dir / f"{config.out_prefix}.{suffix}.xlsx",
                sheet_name=suffix,
                logger=logger,
            )


def run(*, config: Config, logger: logging.Logger) -> Dict[str, pd.DataFrame]:
    """Run candidate CDS consequence classification."""
    candidates = read_table(path=config.candidate_table, logger=logger)
    annotation = read_table(path=config.transcript_annotation, logger=logger)
    features = read_table(path=config.transcript_features, logger=logger)
    classified = classify_candidates(
        candidates=candidates,
        annotation=annotation,
        features=features,
        logger=logger,
    )
    summary = summarise_classification(dataframe=classified)
    cds_changed = classified.loc[
        classified["biological_review_class"] == "cds_changing_candidate"
    ].copy()
    outputs = {
        "candidate_cds_consequence_classification": classified,
        "summary_by_cds_consequence": summary,
        "cds_changing_candidates": cds_changed,
    }
    write_outputs(outputs=outputs, config=config, logger=logger)
    return outputs


def main(argv: Optional[Sequence[str]] = None) -> None:
    """Run the command-line program."""
    config = parse_args(argv)
    logger = setup_logging(log_level=config.log_level, log_path=config.log_path)
    logger.info("Starting candidate isoform CDS consequence classification")
    run(config=config, logger=logger)
    logger.info("Finished candidate isoform CDS consequence classification")


if __name__ == "__main__":
    main()
