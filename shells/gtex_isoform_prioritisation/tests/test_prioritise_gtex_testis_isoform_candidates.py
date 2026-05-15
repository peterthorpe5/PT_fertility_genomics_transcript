"""Tests for prioritise_gtex_testis_isoform_candidates.py."""

from __future__ import annotations

import importlib.util
import logging
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd


SCRIPT_PATH = (
    Path(__file__).resolve().parents[1]
    / "scripts"
    / "prioritise_gtex_testis_isoform_candidates.py"
)

spec = importlib.util.spec_from_file_location("priority", SCRIPT_PATH)
priority = importlib.util.module_from_spec(spec)
sys.modules["priority"] = priority
spec.loader.exec_module(priority)


class TestPrioritiseGtexTestisIsoformCandidates(unittest.TestCase):
    """Unit tests for GTEx isoform candidate prioritisation."""

    def setUp(self):
        """Create a quiet logger for tests."""
        self.logger = logging.getLogger("test_priority")
        self.logger.handlers = []
        self.logger.addHandler(logging.NullHandler())
        self.logger.setLevel(logging.CRITICAL)

    def make_input_dataframe(self):
        """Make a small annotated candidate-like input table."""
        return pd.DataFrame(
            data={
                "gene_symbol": ["GENE1", "GENE2", "GENE3", "GENE4"],
                "transcript_id": ["ENST1", "ENST2", "ENST3", "ENST4"],
                "is_target_tissue_isoform_candidate": [1, 1, 1, 1],
                "is_broad_gene_isoform_rescue_candidate": [1, 1, 0, 0],
                "candidate_rank_tier": [
                    "primary_rank_1",
                    "primary_rank_1",
                    "primary_rank_1",
                    "secondary_rank_2_to_3",
                ],
                "target_median_tpm": [5.0, 4.0, 3.0, 1.5],
                "target_median_isoform_usage": [0.60, 0.55, 0.70, 0.30],
                "log2_target_vs_max_non_target_isoform_usage": [3.0, 2.0, 4.0, 1.2],
                "target_fraction_samples_tpm_ge_threshold": [0.90, 0.80, 0.70, 0.30],
                "target_tpm_rank_within_gene": [1, 1, 1, 2],
                "target_usage_rank_within_gene": [1, 1, 1, 2],
                "gencode_has_cds": [1, 0, 1, 1],
                "gencode_is_protein_coding_transcript": [1, 0, 1, 1],
                "is_basic": [1, 1, 1, 0],
                "is_mane_select": [0, 0, 1, 0],
                "is_ensembl_canonical": [0, 0, 1, 0],
            }
        )

    def make_config(self, out_dir):
        """Make a default test configuration."""
        return priority.Config(
            input_tsv=out_dir / "input.tsv",
            out_dir=out_dir,
            out_prefix="test_priority",
            gene_symbol_col="gene_symbol",
            top_n_plot_genes=3,
            batch_size=2,
            selected_genes=("GENE1", "MISSING"),
            min_strong_tpm=1.0,
            min_strong_usage=0.25,
            min_strong_log2_usage_ratio=1.0,
            min_detection_fraction=0.20,
            rescue_column="is_broad_gene_isoform_rescue_candidate",
            candidate_column="is_target_tissue_isoform_candidate",
            write_excel_outputs=False,
            log_path=None,
            log_level="CRITICAL",
        )

    def test_parse_boolean_series(self):
        """Mixed boolean-like values should be parsed correctly."""
        values = pd.Series(data=[1, 0, "true", "False", "yes", "", None])
        parsed = priority.parse_boolean_series(series=values)
        self.assertEqual(parsed.tolist(), [True, False, True, False, True, False, False])

    def test_assign_priority_tiers(self):
        """Rescue coding candidates should be tier 1."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_dir = Path(tmp_dir)
            config = self.make_config(out_dir=out_dir)
            prioritised = priority.prioritise_candidates(
                dataframe=self.make_input_dataframe(),
                config=config,
                logger=self.logger,
            )

        row_gene1 = prioritised.loc[prioritised["gene_symbol"] == "GENE1"].iloc[0]
        row_gene2 = prioritised.loc[prioritised["gene_symbol"] == "GENE2"].iloc[0]
        row_gene3 = prioritised.loc[prioritised["gene_symbol"] == "GENE3"].iloc[0]

        self.assertEqual(
            row_gene1["priority_tier"],
            "tier_1_protein_coding_rescue_candidate",
        )
        self.assertEqual(
            row_gene2["priority_tier"],
            "tier_2_rescue_non_coding_or_uncertain_candidate",
        )
        self.assertEqual(
            row_gene3["priority_tier"],
            "tier_3_protein_coding_high_confidence_candidate",
        )

    def test_priority_sorting_puts_tier1_first(self):
        """Tier 1 rows should sort before lower-priority tiers."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_dir = Path(tmp_dir)
            config = self.make_config(out_dir=out_dir)
            prioritised = priority.prioritise_candidates(
                dataframe=self.make_input_dataframe(),
                config=config,
                logger=self.logger,
            )

        self.assertEqual(prioritised.iloc[0]["gene_symbol"], "GENE1")
        self.assertLessEqual(
            prioritised.iloc[0]["priority_tier_rank"],
            prioritised.iloc[-1]["priority_tier_rank"],
        )

    def test_gene_list_preserves_order_and_uniqueness(self):
        """Gene lists should preserve ranked order and remove duplicates."""
        data = pd.DataFrame(data={"gene_symbol": ["A", "B", "A", "C"]})
        genes = priority.make_unique_gene_list(
            dataframe=data,
            gene_symbol_col="gene_symbol",
            limit=3,
        )
        self.assertEqual(genes, ["A", "B", "C"])

    def test_end_to_end_writes_expected_outputs(self):
        """A small run should write prioritised tables and gene lists."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_dir = Path(tmp_dir)
            config = self.make_config(out_dir=out_dir)
            input_df = self.make_input_dataframe()
            input_df.to_csv(path_or_buf=config.input_tsv, sep="\t", index=False)

            priority.run(config=config, logger=self.logger)

            expected_files = [
                out_dir / "test_priority.prioritised_all.tsv",
                out_dir / "test_priority.summary_by_priority_tier.tsv",
                out_dir / "test_priority.tier1_protein_coding_rescue_candidates.tsv",
                out_dir / "test_priority.genes_to_plot_top_3.txt",
                out_dir / "plot_gene_batches_top_ranked" / "genes_to_plot_batch_001.txt",
                out_dir / "plot_gene_batches_top_ranked" / "genes_to_plot_batch_002.txt",
            ]
            for path in expected_files:
                self.assertTrue(path.exists(), f"Expected file missing: {path}")

            prioritised = pd.read_csv(
                filepath_or_buffer=out_dir / "test_priority.prioritised_all.tsv",
                sep="\t",
            )
            self.assertIn("priority_score", prioritised.columns)
            self.assertIn("priority_tier", prioritised.columns)


if __name__ == "__main__":
    unittest.main()
