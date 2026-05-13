"""Unit tests for plot_candidate_isoform_gene_models.py."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import pandas as pd

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import plot_candidate_isoform_gene_models as plotter


class TestPlotCandidateIsoformGeneModels(unittest.TestCase):
    """Tests for candidate isoform plotting helpers."""

    def make_features(self) -> pd.DataFrame:
        """Return a small feature table for plotting tests."""
        return pd.DataFrame(
            {
                "gencode_gene_name": ["GENE1"] * 6,
                "gencode_transcript_id_with_version": [
                    "ENST1.1",
                    "ENST1.1",
                    "ENST1.1",
                    "ENST2.1",
                    "ENST2.1",
                    "ENST2.1",
                ],
                "gencode_transcript_name": ["GENE1-201"] * 3 + ["GENE1-202"] * 3,
                "gencode_transcript_type": ["protein_coding"] * 6,
                "gencode_feature_type": ["exon", "CDS", "exon", "exon", "CDS", "exon"],
                "gencode_seqname": ["chr1"] * 6,
                "gencode_start": [100, 130, 300, 120, 150, 420],
                "gencode_end": [200, 200, 500, 240, 220, 550],
                "gencode_strand": ["+"] * 6,
            }
        )

    def test_select_genes_to_plot_uses_user_genes(self) -> None:
        """Explicit genes should be used before automatic ranking."""
        candidates = pd.DataFrame(
            {"gene_symbol": ["A", "B"], "transcript_id": ["T1", "T2"]}
        )
        candidates = plotter.normalise_gene_name_columns(candidates=candidates)
        config = plotter.Config(
            features_tsv=Path("features.tsv"),
            candidates_tsv=Path("candidates.tsv"),
            out_dir=Path("out"),
            genes=("CFAP99",),
            max_genes=10,
            max_transcripts_per_gene=20,
            formats=("pdf",),
            title_suffix="",
            log_path=None,
            log_level="ERROR",
        )
        self.assertEqual(
            plotter.select_genes_to_plot(candidates=candidates, config=config),
            ["CFAP99"],
        )

    def test_prepare_gene_features_subsets_and_converts_coordinates(self) -> None:
        """Gene feature preparation should subset and convert coordinates."""
        features = self.make_features().astype(str)
        subset = plotter.prepare_gene_features(features=features, gene_symbol="GENE1")
        self.assertEqual(subset.shape[0], 6)
        self.assertTrue(pd.api.types.is_integer_dtype(subset["gencode_start"]))

    def test_draw_gene_model_writes_pdf(self) -> None:
        """The drawing function should create a figure file."""
        features = self.make_features()
        candidates = pd.DataFrame(
            {
                "plot_gene_symbol": ["GENE1"],
                "plot_transcript_id": ["ENST1.1"],
                "candidate_rank_tier": ["primary_rank_1"],
                "target_median_tpm": ["12.4"],
                "target_median_isoform_usage": ["0.72"],
                "log2_target_vs_max_non_target_isoform_usage": ["3.1"],
            }
        )
        transcript_order = ["ENST1.1", "ENST2.1"]
        logger = plotter.setup_logging(log_level="ERROR", log_path=None)
        with tempfile.TemporaryDirectory() as tmpdir:
            output = Path(tmpdir) / "gene_model.pdf"
            plotter.draw_gene_model(
                gene_symbol="GENE1",
                gene_features=features,
                gene_candidates=candidates,
                transcript_order=transcript_order,
                output_paths=[output],
                title_suffix="test",
                logger=logger,
            )
            self.assertTrue(output.exists())
            self.assertGreater(output.stat().st_size, 0)


if __name__ == "__main__":
    unittest.main()
