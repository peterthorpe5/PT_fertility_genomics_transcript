"""Unit tests for the GTEx transcript isoform first-pass workflow."""

from __future__ import annotations

import logging
import tempfile
import unittest
from pathlib import Path

import pandas as pd

import sys

SCRIPT_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

import gtex_transcript_isoform_first_pass as isoform  # noqa: E402


class GtexTranscriptIsoformFirstPassTests(unittest.TestCase):
    """Tests for target resolution, extraction, and isoform usage."""

    def setUp(self) -> None:
        """Create a quiet logger for tests."""
        self.logger = logging.getLogger("test_gtex_transcript_isoform")
        self.logger.handlers = []
        self.logger.addHandler(logging.NullHandler())
        self.logger.setLevel(logging.CRITICAL)

    def test_strip_ensembl_version(self) -> None:
        """Version suffixes should be removed from Ensembl IDs."""
        self.assertEqual(
            isoform.strip_ensembl_version(identifier="ENSG000001.12"),
            "ENSG000001",
        )
        self.assertEqual(
            isoform.strip_ensembl_version(identifier="ENST000001.3"),
            "ENST000001",
        )

    def test_load_gene_symbol_map_from_gct(self) -> None:
        """Gene-symbol maps should be loaded from a minimal GCT file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            gct_path = Path(tmpdir) / "mini.gct"
            gct_path.write_text(
                "#1.2\n"
                "2\t2\n"
                "Name\tDescription\tS1\tS2\n"
                "ENSG000001.1\tGENEA\t1\t2\n"
                "ENSG000002.5\tGENEB\t3\t4\n"
            )
            gene_id_to_symbol, symbol_to_gene_ids = isoform.load_gene_symbol_map(
                gene_tpm_gct_path=gct_path,
                logger=self.logger,
            )
            self.assertEqual(gene_id_to_symbol["ENSG000001"], "GENEA")
            self.assertEqual(symbol_to_gene_ids["GENEB"], {"ENSG000002"})

    def test_resolve_targets_accepts_symbols_and_ensembl_ids(self) -> None:
        """Target resolver should accept symbols and Ensembl gene IDs."""
        resolved = isoform.resolve_targets(
            requested_terms=("GENEA", "ENSG000002.5", "MISSING"),
            symbol_to_gene_ids={"GENEA": {"ENSG000001"}},
            logger=self.logger,
        )
        self.assertEqual(set(resolved.resolved_gene_ids), {"ENSG000001", "ENSG000002"})
        self.assertEqual(resolved.unresolved_terms, ("MISSING",))

    def test_extract_transcript_matrix_subset(self) -> None:
        """The extractor should retain all transcripts from requested genes."""
        with tempfile.TemporaryDirectory() as tmpdir:
            matrix_path = Path(tmpdir) / "transcripts_tpm.txt"
            matrix_path.write_text(
                "transcript_id\tgene_id\tS1\tS2\n"
                "ENST000001.1\tENSG000001.1\t10\t20\n"
                "ENST000002.1\tENSG000001.1\t30\t40\n"
                "ENST000003.1\tENSG000009.1\t50\t60\n"
            )
            subset = isoform.extract_transcript_matrix_subset(
                matrix_path=matrix_path,
                target_gene_ids={"ENSG000001"},
                gene_id_to_symbol={"ENSG000001": "GENEA"},
                chunk_size=2,
                logger=self.logger,
                matrix_label="TPM",
            )
            self.assertEqual(subset.shape[0], 2)
            self.assertEqual(set(subset["transcript_id"]), {"ENST000001", "ENST000002"})
            self.assertEqual(set(subset["gene_symbol"]), {"GENEA"})

    def test_compute_isoform_usage_and_summary(self) -> None:
        """Isoform usage should identify a testis-preferential transcript."""
        metadata = pd.DataFrame(
            {
                "transcript_id_with_version": ["ENST1.1", "ENST2.1", "ENST3.1"],
                "transcript_id": ["ENST1", "ENST2", "ENST3"],
                "gene_id_with_version": ["ENSG1.1", "ENSG1.1", "ENSG2.1"],
                "gene_id": ["ENSG1", "ENSG1", "ENSG2"],
                "gene_symbol": ["GENEA", "GENEA", "GENEB"],
            }
        )
        expression = pd.DataFrame(
            {
                "T1": [90.0, 10.0, 5.0],
                "T2": [80.0, 20.0, 5.0],
                "L1": [10.0, 90.0, 5.0],
                "L2": [20.0, 80.0, 5.0],
            }
        )
        sample_to_tissue = pd.Series(
            {"T1": "Testis", "T2": "Testis", "L1": "Liver", "L2": "Liver"}
        )
        usage = isoform.compute_isoform_usage(
            metadata=metadata,
            expression=expression,
            epsilon=1e-8,
            logger=self.logger,
        )
        self.assertGreater(float(usage.loc[0, "T1"]), 0.89)
        self.assertLess(float(usage.loc[0, "L1"]), 0.11)

        tissue_medians = isoform.compute_tissue_medians(
            metadata=metadata,
            expression=expression,
            sample_to_tissue=sample_to_tissue,
            logger=self.logger,
        )
        usage_medians = isoform.compute_usage_tissue_medians(
            metadata=metadata,
            usage=usage,
            sample_to_tissue=sample_to_tissue,
            logger=self.logger,
        )
        summary = isoform.compute_summary(
            metadata=metadata,
            expression=expression,
            tissue_medians=tissue_medians,
            usage_tissue_medians=usage_medians,
            sample_to_tissue=sample_to_tissue,
            target_tissue="Testis",
            min_tpm_present=1.0,
            epsilon=1e-8,
            logger=self.logger,
        )
        row = summary.loc[summary["transcript_id"] == "ENST1", :].iloc[0]
        self.assertEqual(int(row["target_is_max_usage_tissue"]), 1)
        self.assertGreater(float(row["log2_target_vs_max_non_target_isoform_usage"]), 2.0)


if __name__ == "__main__":
    unittest.main()
