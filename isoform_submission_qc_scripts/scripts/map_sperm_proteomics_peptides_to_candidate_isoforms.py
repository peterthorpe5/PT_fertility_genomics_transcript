#!/usr/bin/env python3
"""
Map observed sperm proteomics peptides to candidate transcript isoforms.

Gene-level sperm proteomics supports translation of a locus, but it does not
prove that a specific testis-preferential isoform is translated. This script is
an optional validation layer for MaxQuant peptide outputs. It asks whether any
observed peptide sequence maps to the candidate protein product and whether that
peptide is unique to the candidate isoform among protein products of the same
gene.

Inputs are deliberately flexible because MaxQuant tables can differ slightly.
The peptide table can be ``peptides.txt`` or ``evidence.txt`` and must contain a
peptide sequence column such as ``Sequence`` or ``Modified sequence``.
"""

from __future__ import annotations

import argparse
import gzip
import logging
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

import pandas as pd


LOGGER_NAME = "map_sperm_proteomics_peptides_to_candidate_isoforms"
AMINO_ACID_PATTERN = re.compile(r"[^A-Z]")


@dataclass(frozen=True)
class Config:
    """Runtime configuration."""

    candidate_table: Path
    gencode_translations_fasta: Path
    peptide_table: Path
    out_dir: Path
    out_prefix: str
    sequence_col: Optional[str]
    gene_col: str
    transcript_col: str
    protein_col: str
    min_peptide_length: int
    write_excel_outputs: bool
    log_level: str
    log_path: Optional[Path]


@dataclass(frozen=True)
class ProteinRecord:
    """A parsed protein FASTA record."""

    protein_id: str
    transcript_id: str
    gene_id: str
    transcript_name: str
    gene_symbol: str
    sequence: str


def setup_logging(*, log_level: str, log_path: Optional[Path]) -> logging.Logger:
    """Configure logging."""
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
        description="Map MaxQuant sperm peptides to candidate isoform protein sequences."
    )
    parser.add_argument("--candidate_table", required=True, type=Path)
    parser.add_argument("--gencode_translations_fasta", required=True, type=Path)
    parser.add_argument("--peptide_table", required=True, type=Path)
    parser.add_argument("--out_dir", required=True, type=Path)
    parser.add_argument(
        "--out_prefix",
        default="candidate_isoform_sperm_peptide_support",
    )
    parser.add_argument(
        "--sequence_col",
        default=None,
        help="Peptide sequence column. If omitted, common MaxQuant names are detected.",
    )
    parser.add_argument("--gene_col", default="gene_symbol")
    parser.add_argument("--transcript_col", default="transcript_id_with_version")
    parser.add_argument("--protein_col", default="gencode_protein_id")
    parser.add_argument("--min_peptide_length", type=int, default=7)
    parser.add_argument("--no_write_excel_outputs", action="store_true")
    parser.add_argument("--log_level", default="INFO")
    parser.add_argument("--log_path", type=Path, default=None)
    args = parser.parse_args(argv)
    return Config(
        candidate_table=args.candidate_table,
        gencode_translations_fasta=args.gencode_translations_fasta,
        peptide_table=args.peptide_table,
        out_dir=args.out_dir,
        out_prefix=args.out_prefix,
        sequence_col=args.sequence_col,
        gene_col=args.gene_col,
        transcript_col=args.transcript_col,
        protein_col=args.protein_col,
        min_peptide_length=args.min_peptide_length,
        write_excel_outputs=not args.no_write_excel_outputs,
        log_level=args.log_level,
        log_path=args.log_path,
    )


def open_text(path: Path):
    """Open plain or gzipped text."""
    if "".join(path.suffixes).lower().endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt", encoding="utf-8")


def strip_version(identifier: object) -> str:
    """Strip dot-version suffix."""
    return str(identifier).split(".")[0]


def clean_peptide_sequence(sequence: object) -> str:
    """
    Clean a MaxQuant peptide sequence.

    Modified sequences may include punctuation or modification labels. This
    function keeps only uppercase letters and removes underscores.
    """
    text = str(sequence).replace("_", "")
    text = re.sub(r"\([^)]*\)", "", text)
    text = re.sub(r"\[[^]]*\]", "", text)
    text = text.upper()
    return AMINO_ACID_PATTERN.sub("", text)


def parse_gencode_header(header: str) -> Tuple[str, str, str, str, str]:
    """
    Parse a GENCODE protein FASTA header.

    Parameters
    ----------
    header
        FASTA header without the leading ``>``.

    Returns
    -------
    Tuple[str, str, str, str, str]
        Protein ID, transcript ID, gene ID, transcript name and gene symbol.
    """
    fields = header.split("|")
    protein_id = fields[0].split()[0] if fields else ""
    transcript_id = fields[1] if len(fields) > 1 else ""
    gene_id = fields[2] if len(fields) > 2 else ""
    transcript_name = fields[5] if len(fields) > 5 else ""
    gene_symbol = fields[6] if len(fields) > 6 else ""
    return protein_id, transcript_id, gene_id, transcript_name, gene_symbol


def read_gencode_translations(*, path: Path, logger: logging.Logger) -> List[ProteinRecord]:
    """Read GENCODE protein translations FASTA."""
    logger.info("Reading GENCODE translations from %s", path)
    records: List[ProteinRecord] = []
    header: Optional[str] = None
    sequence_parts: List[str] = []
    with open_text(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    protein_id, transcript_id, gene_id, transcript_name, gene_symbol = parse_gencode_header(header)
                    records.append(
                        ProteinRecord(
                            protein_id=protein_id,
                            transcript_id=transcript_id,
                            gene_id=gene_id,
                            transcript_name=transcript_name,
                            gene_symbol=gene_symbol,
                            sequence="".join(sequence_parts),
                        )
                    )
                header = line[1:]
                sequence_parts = []
            else:
                sequence_parts.append(line.strip())
        if header is not None:
            protein_id, transcript_id, gene_id, transcript_name, gene_symbol = parse_gencode_header(header)
            records.append(
                ProteinRecord(
                    protein_id=protein_id,
                    transcript_id=transcript_id,
                    gene_id=gene_id,
                    transcript_name=transcript_name,
                    gene_symbol=gene_symbol,
                    sequence="".join(sequence_parts),
                )
            )
    logger.info("Loaded %d protein translation records", len(records))
    return records


def read_table(*, path: Path, logger: logging.Logger) -> pd.DataFrame:
    """Read TSV/TSV.GZ/XLSX table."""
    logger.info("Reading %s", path)
    if "".join(path.suffixes).lower().endswith(".xlsx"):
        dataframe = pd.read_excel(path)
    else:
        dataframe = pd.read_csv(path, sep="\t", low_memory=False)
    logger.info("Loaded %d rows x %d columns from %s", *dataframe.shape, path)
    return dataframe


def detect_sequence_column(*, dataframe: pd.DataFrame, requested: Optional[str]) -> str:
    """Detect peptide sequence column."""
    if requested:
        if requested not in dataframe.columns:
            raise ValueError(f"Requested sequence column not found: {requested}")
        return requested
    for candidate in ("Sequence", "Modified sequence", "Peptide sequence", "sequence"):
        if candidate in dataframe.columns:
            return candidate
    raise ValueError(
        "Could not find peptide sequence column. Use --sequence_col to specify it."
    )


def prepare_observed_peptides(
    *, peptide_table: pd.DataFrame, sequence_col: str, min_peptide_length: int
) -> pd.DataFrame:
    """Clean and collapse observed peptide sequences."""
    output = peptide_table.copy()
    output["clean_peptide_sequence"] = output[sequence_col].map(clean_peptide_sequence)
    output = output.loc[output["clean_peptide_sequence"].str.len() >= min_peptide_length]
    summary_cols = ["clean_peptide_sequence"]
    optional_cols = [
        column
        for column in ("Gene names", "Proteins", "Leading razor protein", "Protein group IDs")
        if column in output.columns
    ]
    collapsed = (
        output.groupby("clean_peptide_sequence", dropna=False)
        .agg(
            n_evidence_rows=("clean_peptide_sequence", "size"),
            **{
                f"source_{column.replace(' ', '_').lower()}": (column, lambda x: ";".join(sorted(set(map(str, x.dropna()))))[:3000])
                for column in optional_cols
            },
        )
        .reset_index()
    )
    return collapsed.loc[:, summary_cols + [c for c in collapsed.columns if c not in summary_cols]]


def prepare_candidate_table(*, dataframe: pd.DataFrame, config: Config) -> pd.DataFrame:
    """Normalise candidate identifiers."""
    output = dataframe.copy()
    for column in (config.gene_col, config.transcript_col):
        if column not in output.columns:
            raise ValueError(f"Candidate table is missing required column: {column}")
    if config.protein_col not in output.columns:
        output[config.protein_col] = ""
    output["candidate_gene_symbol"] = output[config.gene_col].astype(str)
    output["candidate_transcript_id_with_version"] = output[config.transcript_col].astype(str)
    output["candidate_transcript_id"] = output[
        "candidate_transcript_id_with_version"
    ].map(strip_version)
    output["candidate_protein_id_with_version"] = output[config.protein_col].astype(str)
    output["candidate_protein_id"] = output["candidate_protein_id_with_version"].map(strip_version)
    return output


def records_to_dataframe(*, records: Sequence[ProteinRecord]) -> pd.DataFrame:
    """Convert protein records to DataFrame."""
    return pd.DataFrame([record.__dict__ for record in records])


def map_peptides_to_candidates(
    *, candidates: pd.DataFrame, proteins: pd.DataFrame, peptides: pd.DataFrame, logger: logging.Logger
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Map observed peptides to candidate proteins and same-gene proteins."""
    proteins = proteins.copy()
    proteins["protein_id_no_version"] = proteins["protein_id"].map(strip_version)
    proteins["transcript_id_no_version"] = proteins["transcript_id"].map(strip_version)

    candidate_genes = set(candidates["candidate_gene_symbol"].astype(str))
    gene_proteins = proteins.loc[proteins["gene_symbol"].isin(candidate_genes)].copy()
    protein_by_transcript = gene_proteins.set_index("transcript_id_no_version")
    protein_by_protein = gene_proteins.set_index("protein_id_no_version")

    evidence_rows: List[Dict[str, object]] = []
    summary_rows: List[Dict[str, object]] = []
    peptide_sequences = peptides["clean_peptide_sequence"].dropna().astype(str).unique().tolist()
    logger.info("Mapping %d unique observed peptides", len(peptide_sequences))

    for _, candidate in candidates.iterrows():
        gene_symbol = str(candidate["candidate_gene_symbol"])
        candidate_tx = str(candidate["candidate_transcript_id"])
        candidate_protein = str(candidate["candidate_protein_id"])

        candidate_record = None
        if candidate_protein and candidate_protein != "nan" and candidate_protein in protein_by_protein.index:
            candidate_record = protein_by_protein.loc[candidate_protein]
            if isinstance(candidate_record, pd.DataFrame):
                candidate_record = candidate_record.iloc[0]
        elif candidate_tx in protein_by_transcript.index:
            candidate_record = protein_by_transcript.loc[candidate_tx]
            if isinstance(candidate_record, pd.DataFrame):
                candidate_record = candidate_record.iloc[0]

        if candidate_record is None:
            summary_rows.append(
                {
                    "gene_symbol": gene_symbol,
                    "transcript_id_with_version": candidate["candidate_transcript_id_with_version"],
                    "candidate_protein_id": candidate_protein,
                    "candidate_protein_found_in_fasta": False,
                    "n_observed_peptides_mapping_candidate": 0,
                    "n_gene_unique_observed_peptides": 0,
                    "has_candidate_peptide_support": False,
                    "has_gene_unique_candidate_peptide_support": False,
                }
            )
            continue

        candidate_sequence = str(candidate_record["sequence"])
        same_gene = gene_proteins.loc[gene_proteins["gene_symbol"] == gene_symbol]
        candidate_hits = 0
        unique_hits = 0
        for peptide in peptide_sequences:
            if peptide not in candidate_sequence:
                continue
            candidate_hits += 1
            containing = same_gene.loc[same_gene["sequence"].str.contains(peptide, regex=False)]
            containing_transcripts = sorted(set(containing["transcript_id"].astype(str)))
            unique_within_gene = len(containing_transcripts) == 1
            if unique_within_gene:
                unique_hits += 1
            peptide_meta = peptides.loc[peptides["clean_peptide_sequence"] == peptide].iloc[0].to_dict()
            evidence_rows.append(
                {
                    "gene_symbol": gene_symbol,
                    "candidate_transcript_id_with_version": candidate[
                        "candidate_transcript_id_with_version"
                    ],
                    "candidate_protein_id": candidate_record["protein_id"],
                    "peptide_sequence": peptide,
                    "peptide_unique_within_gene": unique_within_gene,
                    "n_same_gene_proteins_containing_peptide": len(containing_transcripts),
                    "same_gene_transcripts_containing_peptide": ";".join(containing_transcripts),
                    **peptide_meta,
                }
            )
        summary_rows.append(
            {
                "gene_symbol": gene_symbol,
                "transcript_id_with_version": candidate["candidate_transcript_id_with_version"],
                "candidate_protein_id": candidate_record["protein_id"],
                "candidate_protein_found_in_fasta": True,
                "n_observed_peptides_mapping_candidate": candidate_hits,
                "n_gene_unique_observed_peptides": unique_hits,
                "has_candidate_peptide_support": candidate_hits > 0,
                "has_gene_unique_candidate_peptide_support": unique_hits > 0,
            }
        )

    return pd.DataFrame(summary_rows), pd.DataFrame(evidence_rows)


def write_tsv(*, dataframe: pd.DataFrame, path: Path, logger: logging.Logger) -> None:
    """Write a TSV table."""
    path.parent.mkdir(parents=True, exist_ok=True)
    dataframe.to_csv(path, sep="\t", index=False)
    logger.info("Wrote %d rows x %d columns to %s", *dataframe.shape, path)


def write_excel(*, dataframe: pd.DataFrame, path: Path, sheet_name: str, logger: logging.Logger) -> None:
    """Write a formatted Excel table."""
    if len(dataframe) + 1 > 1_048_576:
        logger.warning("Skipping Excel output for %s because it exceeds Excel limits", path)
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(path, engine="xlsxwriter") as writer:
        dataframe.to_excel(writer, sheet_name=sheet_name[:31], index=False)
        workbook = writer.book
        worksheet = writer.sheets[sheet_name[:31]]
        header_format = workbook.add_format({"bold": True, "text_wrap": True, "border": 1})
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
            worksheet.set_column(idx, idx, width)
    logger.info("Wrote formatted Excel output to %s", path)


def write_outputs(*, outputs: Dict[str, pd.DataFrame], config: Config, logger: logging.Logger) -> None:
    """Write all outputs."""
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
    """Run peptide-to-isoform mapping."""
    candidates = prepare_candidate_table(
        dataframe=read_table(path=config.candidate_table, logger=logger), config=config
    )
    peptide_raw = read_table(path=config.peptide_table, logger=logger)
    sequence_col = detect_sequence_column(dataframe=peptide_raw, requested=config.sequence_col)
    peptides = prepare_observed_peptides(
        peptide_table=peptide_raw,
        sequence_col=sequence_col,
        min_peptide_length=config.min_peptide_length,
    )
    proteins = records_to_dataframe(
        records=read_gencode_translations(path=config.gencode_translations_fasta, logger=logger)
    )
    summary, evidence = map_peptides_to_candidates(
        candidates=candidates, proteins=proteins, peptides=peptides, logger=logger
    )
    evidence_positive = evidence.loc[evidence.get("peptide_unique_within_gene", pd.Series(dtype=bool)) == True].copy() if not evidence.empty else evidence
    outputs = {
        "candidate_peptide_support_summary": summary,
        "candidate_observed_peptide_evidence": evidence,
        "candidate_gene_unique_peptide_evidence": evidence_positive,
    }
    write_outputs(outputs=outputs, config=config, logger=logger)
    return outputs


def main(argv: Optional[Sequence[str]] = None) -> None:
    """Run the command-line program."""
    config = parse_args(argv)
    logger = setup_logging(log_level=config.log_level, log_path=config.log_path)
    logger.info("Starting sperm peptide to candidate isoform mapping")
    run(config=config, logger=logger)
    logger.info("Finished sperm peptide to candidate isoform mapping")


if __name__ == "__main__":
    main()
