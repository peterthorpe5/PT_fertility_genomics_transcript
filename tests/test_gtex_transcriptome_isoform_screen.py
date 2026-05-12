"""Unit tests for gtex_transcriptome_isoform_screen.py."""

from __future__ import annotations

import importlib.util
import logging
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

SCRIPT_PATH = (
    Path(__file__).resolve().parents[1]
    / "scripts"
    / "gtex_transcriptome_isoform_screen.py"
)
SPEC = importlib.util.spec_from_file_location("gtex_screen", SCRIPT_PATH)
gtex_screen = importlib.util.module_from_spec(SPEC)
sys.modules["gtex_screen"] = gtex_screen
assert SPEC.loader is not None
SPEC.loader.exec_module(gtex_screen)


class TestGtexTranscriptomeIsoformScreen(unittest.TestCase):
    """Test core functions in the transcriptome-wide isoform screen."""

    def setUp(self) -> None:
        """Create a quiet logger for tests."""
        self.logger = logging.getLogger("test_gtex_screen")
        self.logger.handlers = []
        self.logger.addHandler(logging.NullHandler())
        self.logger.setLevel(logging.CRITICAL)

    def test_strip_ensembl_version(self) -> None:
        """Ensembl version suffixes should be removed safely."""
        self.assertEqual(
            gtex_screen.strip_ensembl_version(identifier="ENSG000001.12"),
            "ENSG000001",
        )
        self.assertEqual(
            gtex_screen.strip_ensembl_version(identifier="ENST000001"),
            "ENST000001",
        )
        self.assertEqual(gtex_screen.strip_ensembl_version(identifier=None), "None")

    def test_detect_transcript_schema(self) -> None:
        """Transcript, gene, metadata, and sample columns should be detected."""
        columns = ["transcript_id", "gene_id", "length", "S1", "S2", "GTEX-X"]
        schema = gtex_screen.detect_transcript_schema(
            columns=columns,
            sample_ids=["S1", "S2"],
        )
        self.assertEqual(schema.transcript_col, "transcript_id")
        self.assertEqual(schema.gene_col, "gene_id")
        self.assertEqual(schema.sample_cols, ("S1", "S2"))
        self.assertIn("length", schema.metadata_cols)
        self.assertEqual(schema.ignored_expression_like_cols, ("GTEX-X",))

    def test_compute_gene_expression_array(self) -> None:
        """Transcript expression should be summed correctly to gene level."""
        expression = np.array(
            [
                [10.0, 20.0],
                [30.0, 40.0],
                [5.0, 6.0],
            ],
            dtype=np.float32,
        )
        genes, inverse, gene_array = gtex_screen.compute_gene_expression_array(
            expression_array=expression,
            gene_ids=["G1", "G1", "G2"],
            float_dtype="float32",
            logger=self.logger,
        )
        self.assertEqual(list(genes), ["G1", "G2"])
        self.assertEqual(list(inverse), [0, 0, 1])
        np.testing.assert_allclose(gene_array[0], [40.0, 60.0])
        np.testing.assert_allclose(gene_array[1], [5.0, 6.0])

    def test_isoform_usage_tissue_medians(self) -> None:
        """Median isoform usage should be calculated within each gene."""
        expression = np.array(
            [
                [80.0, 60.0, 10.0],
                [20.0, 40.0, 90.0],
                [5.0, 5.0, 5.0],
            ],
            dtype=np.float32,
        )
        genes, inverse, gene_array = gtex_screen.compute_gene_expression_array(
            expression_array=expression,
            gene_ids=["G1", "G1", "G2"],
            float_dtype="float32",
            logger=self.logger,
        )
        self.assertEqual(list(genes), ["G1", "G2"])
        row_metadata = pd.DataFrame(
            {
                "transcript_id_with_version": ["T1.1", "T2.1", "T3.1"],
                "transcript_id": ["T1", "T2", "T3"],
                "gene_id_with_version": ["G1.1", "G1.1", "G2.1"],
                "gene_id": ["G1", "G1", "G2"],
                "gene_symbol": ["GENE1", "GENE1", "GENE2"],
            }
        )
        usage = gtex_screen.compute_isoform_usage_tissue_medians(
            expression_array=expression,
            gene_expression_array=gene_array,
            transcript_gene_inverse=inverse,
            row_metadata=row_metadata,
            tissue_to_indices={"Testis": np.array([0, 1]), "Brain": np.array([2])},
            epsilon=1e-8,
            logger=self.logger,
        )
        t1_testis = usage.loc[usage["transcript_id"] == "T1", "Testis"].iloc[0]
        t2_brain = usage.loc[usage["transcript_id"] == "T2", "Brain"].iloc[0]
        t3_testis = usage.loc[usage["transcript_id"] == "T3", "Testis"].iloc[0]
        self.assertAlmostEqual(t1_testis, 0.7, places=6)
        self.assertAlmostEqual(t2_brain, 0.9, places=6)
        self.assertAlmostEqual(t3_testis, 1.0, places=6)

    def test_end_to_end_small_run(self) -> None:
        """A tiny GTEx-like input should produce candidate outputs."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            transcript_path = tmp_path / "transcripts_tpm.tsv"
            gene_gct_path = tmp_path / "gene_tpm.gct"
            sample_attributes_path = tmp_path / "sample_attributes.tsv"
            out_dir = tmp_path / "out"

            transcript_path.write_text(
                "transcript_id\tgene_id\tS_testis_1\tS_testis_2\tS_brain_1\tS_heart_1\n"
                "T1.1\tG1.1\t80\t60\t10\t20\n"
                "T2.1\tG1.1\t20\t40\t90\t80\n"
                "T3.1\tG2.1\t5\t5\t5\t5\n"
            )
            gene_gct_path.write_text(
                "#1.2\n"
                "3\t4\n"
                "Name\tDescription\tS_testis_1\tS_testis_2\tS_brain_1\tS_heart_1\n"
                "G1.1\tGENE1\t100\t100\t100\t100\n"
                "G2.1\tGENE2\t5\t5\t5\t5\n"
            )
            sample_attributes_path.write_text(
                "SAMPID\tSMTSD\n"
                "S_testis_1\tTestis\n"
                "S_testis_2\tTestis\n"
                "S_brain_1\tBrain\n"
                "S_heart_1\tHeart\n"
            )

            config = gtex_screen.Config(
                transcript_tpm_path=transcript_path,
                gene_tpm_gct_path=gene_gct_path,
                sample_attributes_path=sample_attributes_path,
                out_dir=out_dir,
                out_prefix="test",
                target_tissue="Testis",
                sample_id_col="SAMPID",
                tissue_col="SMTSD",
                float_dtype="float32",
                min_tpm_present=1.0,
                min_target_tpm_candidate=1.0,
                min_target_usage_candidate=0.25,
                min_log2_usage_ratio_candidate=0.5,
                max_gene_log2_ratio_for_rescue=1.0,
                top_n_candidate_tissues=3,
                write_tissue_matrices=True,
                log_path=None,
                log_level="CRITICAL",
            )
            gtex_screen.run(config=config, logger=self.logger)

            candidate_path = out_dir / "test.candidate_target_tissue_isoforms.tsv"
            best_path = out_dir / "test.best_candidate_isoform_per_gene.tsv"
            summary_path = out_dir / "test.transcript_target_tissue_isoform_summary.tsv.gz"

            self.assertTrue(candidate_path.exists())
            self.assertTrue(best_path.exists())
            self.assertTrue(summary_path.exists())

            candidates = pd.read_csv(candidate_path, sep="\t")
            self.assertIn("T1", set(candidates["transcript_id"]))
            t1 = candidates.loc[candidates["transcript_id"] == "T1"].iloc[0]
            self.assertEqual(t1["is_broad_gene_isoform_rescue_candidate"], 1)
            self.assertGreater(t1["target_median_isoform_usage"], 0.5)


if __name__ == "__main__":
    unittest.main()
