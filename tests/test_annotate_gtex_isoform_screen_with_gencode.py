"""Unit tests for annotate_gtex_isoform_screen_with_gencode.py."""

from __future__ import annotations

import gzip
import tempfile
import unittest
from pathlib import Path

import pandas as pd

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import annotate_gtex_isoform_screen_with_gencode as annotate


class TestAnnotateGtexIsoformScreenWithGencode(unittest.TestCase):
    """Tests for GENCODE annotation helpers."""

    def test_parse_gtf_attributes_preserves_repeated_tags(self) -> None:
        """Repeated GTF tag attributes should be retained as a list."""
        attributes = annotate.parse_gtf_attributes(
            attribute_text=(
                'gene_id "ENSG1.1"; transcript_id "ENST1.1"; '
                'gene_name "GENE1"; tag "basic"; tag "MANE_Select";'
            )
        )
        self.assertEqual(attributes["gene_id"], ["ENSG1.1"])
        self.assertEqual(attributes["tag"], ["basic", "MANE_Select"])

    def test_strip_ensembl_version(self) -> None:
        """Ensembl version suffixes should be stripped safely."""
        self.assertEqual(
            annotate.strip_ensembl_version(identifier="ENST000001.12"),
            "ENST000001",
        )
        self.assertEqual(annotate.strip_ensembl_version(identifier=""), "")

    def test_parse_gencode_gtf_builds_annotation_and_features(self) -> None:
        """A small GTF should produce transcript and feature tables."""
        base_attrs = (
            'gene_id "ENSG1.1"; gene_name "GENE1"; '
            'transcript_id "ENST1.1"; transcript_name "GENE1-201"; '
            'transcript_type "protein_coding";'
        )
        transcript_attrs = (
            'gene_id "ENSG1.1"; gene_name "GENE1"; '
            'gene_type "protein_coding"; transcript_id "ENST1.1"; '
            'transcript_name "GENE1-201"; '
            'transcript_type "protein_coding"; tag "basic"; '
            'tag "MANE_Select";'
        )
        gtf_text = "\n".join(
            [
                f"chr1\tsrc\ttranscript\t100\t500\t.\t+\t.\t{transcript_attrs}",
                f"chr1\tsrc\texon\t100\t200\t.\t+\t.\t{base_attrs} exon_number \"1\";",
                (
                    f"chr1\tsrc\tCDS\t130\t200\t.\t+\t0\t{base_attrs} "
                    "protein_id \"ENSP1.1\";"
                ),
                f"chr1\tsrc\texon\t300\t500\t.\t+\t.\t{base_attrs} exon_number \"2\";",
            ]
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "mini.gtf.gz"
            with gzip.open(path, "wt") as handle:
                handle.write(gtf_text)
            logger = annotate.setup_logging(log_level="ERROR", log_path=None)
            transcripts, features = annotate.parse_gencode_gtf(
                gtf_path=path,
                logger=logger,
            )
        self.assertEqual(transcripts.shape[0], 1)
        row = transcripts.iloc[0]
        self.assertEqual(row["gencode_exon_count"], 2)
        self.assertEqual(row["gencode_cds_exon_count"], 1)
        self.assertEqual(row["gencode_has_cds"], 1)
        self.assertEqual(row["is_basic"], 1)
        self.assertEqual(row["is_mane_select"], 1)
        self.assertEqual(features.shape[0], 3)

    def test_annotate_result_table_matches_versioned_transcript(self) -> None:
        """Result rows should be annotated by versioned transcript ID."""
        result = pd.DataFrame(
            {
                "transcript_id_with_version": ["ENST1.1", "ENST2.1"],
                "transcript_id": ["ENST1", "ENST2"],
                "gene_symbol": ["GENE1", "GENE2"],
            }
        )
        annotation = pd.DataFrame(
            {
                "gencode_transcript_id_with_version": ["ENST1.1"],
                "gencode_transcript_id": ["ENST1"],
                "gencode_gene_id": ["ENSG1"],
                "gencode_gene_name": ["GENE1"],
            }
        )
        logger = annotate.setup_logging(log_level="ERROR", log_path=None)
        output = annotate.annotate_result_table(
            dataframe=result,
            transcript_annotation=annotation,
            logger=logger,
        )
        self.assertEqual(
            output.loc[0, "gencode_annotation_match_type"],
            "versioned_transcript_id",
        )
        self.assertEqual(output.loc[1, "gencode_annotation_match_type"], "unmatched")

    def test_make_safe_excel_sheet_name(self) -> None:
        """Excel sheet names should be shortened and cleaned."""
        safe = annotate.make_safe_excel_sheet_name(
            sheet_name="candidate/isoforms:annotated*long_name_excessive"
        )
        self.assertLessEqual(len(safe), 31)
        self.assertNotIn("/", safe)
        self.assertNotIn(":", safe)


if __name__ == "__main__":
    unittest.main()
