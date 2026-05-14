"""Unit tests for enhanced isoform gene-model plotting."""

from __future__ import annotations

import logging
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from scripts.plot_candidate_isoform_gene_models_v2 import (
    choose_transcripts,
    make_candidate_lookup,
    merge_intervals,
    normalise_candidates,
    normalise_features,
    plot_gene_model,
    strip_ensembl_version,
)


class TestEnhancedGeneModelPlots(unittest.TestCase):
    """Test enhanced candidate isoform plotting helpers."""

    def setUp(self) -> None:
        """Create small synthetic inputs."""
        self.features_raw = pd.DataFrame(
            {
                "gene_name": ["GENE1"] * 9,
                "gene_id": ["ENSG1.1"] * 9,
                "transcript_id": [
                    "ENST1.1",
                    "ENST1.1",
                    "ENST1.1",
                    "ENST2.1",
                    "ENST2.1",
                    "ENST2.1",
                    "ENST2.1",
                    "ENST3.1",
                    "ENST3.1",
                ],
                "transcript_name": [
                    "GENE1-201",
                    "GENE1-201",
                    "GENE1-201",
                    "GENE1-202",
                    "GENE1-202",
                    "GENE1-202",
                    "GENE1-202",
                    "GENE1-203",
                    "GENE1-203",
                ],
                "feature": [
                    "exon",
                    "exon",
                    "CDS",
                    "exon",
                    "exon",
                    "exon",
                    "CDS",
                    "exon",
                    "exon",
                ],
                "seqname": ["chr1"] * 9,
                "start": [100, 300, 120, 100, 200, 300, 210, 500, 700],
                "end": [150, 360, 140, 150, 240, 360, 230, 550, 750],
                "strand": ["+"] * 9,
            }
        )
        self.candidates_raw = pd.DataFrame(
            {
                "gene_symbol": ["GENE1", "GENE1"],
                "transcript_id_with_version": ["ENST1.1", "ENST2.1"],
                "transcript_id": ["ENST1", "ENST2"],
                "candidate_rank_tier": [
                    "primary_rank_1",
                    "secondary_rank_2_to_3",
                ],
                "target_median_tpm": ["4.1", "2.2"],
                "target_median_isoform_usage": ["0.55", "0.31"],
                "log2_target_vs_max_non_target_isoform_usage": ["3.2", "2.1"],
                "target_usage_rank_within_gene": ["1", "2"],
            }
        )

    def test_strip_ensembl_version(self) -> None:
        """Ensembl dot-version suffixes are stripped."""
        self.assertEqual(
            strip_ensembl_version(identifier="ENST000001.5"),
            "ENST000001",
        )

    def test_merge_intervals(self) -> None:
        """Overlapping and adjacent intervals are merged."""
        merged = merge_intervals(
            intervals=[(1, 10), (8, 20), (21, 25), (40, 50)]
        )
        self.assertEqual(merged, [(1, 25), (40, 50)])

    def test_normalise_features(self) -> None:
        """Feature tables are normalised to expected columns."""
        features = normalise_features(dataframe=self.features_raw)
        self.assertEqual(features.loc[0, "gene_id"], "ENSG1")
        self.assertEqual(features.loc[0, "transcript_id"], "ENST1")
        self.assertEqual(set(features["feature"]), {"exon", "cds"})

    def test_candidate_priority(self) -> None:
        """Primary and secondary candidates are prioritised in track order."""
        features = normalise_features(dataframe=self.features_raw)
        candidates = normalise_candidates(dataframe=self.candidates_raw)
        lookup = make_candidate_lookup(candidates=candidates)
        chosen = choose_transcripts(
            exon_features=features[features["feature"] == "exon"],
            lookup=lookup,
            max_transcripts=2,
        )
        self.assertEqual(chosen, ["ENST1", "ENST2"])

    def test_plot_gene_model_writes_files(self) -> None:
        """A synthetic gene model plot is written."""
        logger = logging.getLogger("test_plot_gene_model")
        logger.handlers = []
        logger.addHandler(logging.NullHandler())
        features = normalise_features(dataframe=self.features_raw)
        candidates = normalise_candidates(dataframe=self.candidates_raw)
        with tempfile.TemporaryDirectory() as tmpdir:
            out_dir = Path(tmpdir)
            plot_gene_model(
                gene_symbol="GENE1",
                features=features,
                candidates=candidates,
                out_dir=out_dir,
                output_formats=("pdf", "svg"),
                out_suffix="test_model",
                max_transcripts=20,
                label_union_exons=True,
                title_suffix="unit test",
                dpi=100,
                logger=logger,
            )
            self.assertTrue((out_dir / "GENE1.test_model.pdf").exists())
            self.assertTrue((out_dir / "GENE1.test_model.svg").exists())


if __name__ == "__main__":
    unittest.main()
