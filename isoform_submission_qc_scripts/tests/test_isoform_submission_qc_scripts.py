"""Unit tests for isoform submission QC helper scripts."""

from __future__ import annotations

import gzip
import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def load_module(script_name: str):
    """Load a script module by filename."""
    path = ROOT / "scripts" / script_name
    spec = importlib.util.spec_from_file_location(script_name.replace(".py", ""), path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class TestAuditIntegrationConsistency(unittest.TestCase):
    """Tests for integration consistency audit functions."""

    @classmethod
    def setUpClass(cls):
        cls.mod = load_module("audit_isoform_evidence_integration_consistency.py")

    def test_compare_key_sets_detects_post_only_rows(self):
        """Post-only gene-transcript keys should be detected."""
        pre = pd.DataFrame(
            {
                "gene_symbol": ["A", "B"],
                "transcript_id_with_version": ["T1.1", "T2.1"],
            }
        )
        post = pd.DataFrame(
            {
                "gene_symbol": ["A", "B", "C"],
                "transcript_id_with_version": ["T1.1", "T2.1", "T3.1"],
                "priority_tier": [
                    "tier_1_protein_coding_rescue_candidate",
                    "tier_1_protein_coding_rescue_candidate",
                    "tier_2_rescue_non_coding_or_uncertain_candidate",
                ],
            }
        )
        pre = self.mod.normalise_key_columns(
            dataframe=pre,
            gene_col="gene_symbol",
            transcript_col="transcript_id_with_version",
        )
        post = self.mod.normalise_key_columns(
            dataframe=post,
            gene_col="gene_symbol",
            transcript_col="transcript_id_with_version",
        )
        summary, post_only, pre_only = self.mod.compare_key_sets(pre=pre, post=post)
        self.assertEqual(int(summary.loc[0, "n_post_only_keys"]), 1)
        self.assertEqual(len(post_only), 1)
        self.assertEqual(len(pre_only), 0)

    def test_protein_coding_conflict_detected(self):
        """Integrated protein-coding true with strict false should be flagged."""
        dataframe = pd.DataFrame(
            {
                "isoform_is_protein_coding": [True, True],
                "gencode_is_protein_coding_transcript": [0, 1],
                "flag_protein_coding_transcript": [False, True],
            }
        )
        conflicts = self.mod.find_protein_coding_flag_conflicts(dataframe=dataframe)
        self.assertEqual(len(conflicts), 1)


class TestStrictRebuild(unittest.TestCase):
    """Tests for strict integrated table rebuilding."""

    @classmethod
    def setUpClass(cls):
        cls.mod = load_module("rebuild_strict_integrated_isoform_tables.py")

    def test_strict_mask_excludes_tier2_and_non_cds_rows(self):
        """Only strict tier-1 protein-coding CDS rescue rows should pass."""
        dataframe = pd.DataFrame(
            {
                "priority_tier": [
                    "tier_1_protein_coding_rescue_candidate",
                    "tier_2_rescue_non_coding_or_uncertain_candidate",
                    "tier_1_protein_coding_rescue_candidate",
                ],
                "gencode_transcript_type": [
                    "protein_coding",
                    "protein_coding_CDS_not_defined",
                    "protein_coding",
                ],
                "gencode_has_cds": [1, 0, 0],
                "gencode_cds_length_bp": [300, 0, 0],
                "isoform_is_rescue_candidate": [True, True, True],
            }
        )
        config = self.mod.Config(
            integrated_table=Path("dummy.tsv"),
            out_dir=Path("dummy"),
            out_prefix="dummy",
            selected_genes=[],
            require_priority_tier1=True,
            require_strict_protein_coding=True,
            require_cds=True,
            require_rescue=True,
            write_excel_outputs=False,
            log_level="CRITICAL",
            log_path=None,
        )
        mask = self.mod.build_strict_mask(dataframe=dataframe, config=config)
        self.assertEqual(mask.tolist(), [True, False, False])


class TestCdsClassification(unittest.TestCase):
    """Tests for CDS consequence classification."""

    @classmethod
    def setUpClass(cls):
        cls.mod = load_module("classify_candidate_isoform_cds_consequences.py")

    def test_interval_comparison_identifies_partial_change(self):
        """CDS intervals with candidate- and reference-specific bases should differ."""
        result = self.mod.compare_cds_intervals(
            candidate_cds=[(1, 100), (201, 300)],
            reference_cds=[(1, 100), (401, 500)],
        )
        self.assertEqual(result["cds_change_class"], "cds_partial_overlap_change")
        self.assertTrue(result["cds_changed_relative_to_reference"])

    def test_classify_candidates_detects_cds_change(self):
        """A candidate with different CDS intervals should be classified as CDS-changing."""
        candidates = pd.DataFrame(
            {
                "gene_symbol": ["GENE1"],
                "transcript_id_with_version": ["TX2.1"],
            }
        )
        annotation = pd.DataFrame(
            {
                "gencode_gene_id": ["G1", "G1"],
                "gencode_gene_name": ["GENE1", "GENE1"],
                "gencode_transcript_id_with_version": ["TX1.1", "TX2.1"],
                "is_mane_select": [1, 0],
                "is_ensembl_canonical": [1, 0],
                "is_basic": [1, 1],
                "gencode_has_cds": [1, 1],
                "gencode_cds_length_bp": [200, 200],
            }
        )
        features = pd.DataFrame(
            {
                "gencode_gene_id": ["G1", "G1", "G1", "G1"],
                "gencode_gene_name": ["GENE1", "GENE1", "GENE1", "GENE1"],
                "gencode_transcript_id_with_version": ["TX1.1", "TX1.1", "TX2.1", "TX2.1"],
                "gencode_feature_type": ["CDS", "CDS", "CDS", "CDS"],
                "gencode_start": [1, 201, 1, 401],
                "gencode_end": [100, 300, 100, 500],
            }
        )
        logger = self.mod.setup_logging(log_level="CRITICAL", log_path=None)
        result = self.mod.classify_candidates(
            candidates=candidates,
            annotation=annotation,
            features=features,
            logger=logger,
        )
        self.assertEqual(result.loc[0, "biological_review_class"], "cds_changing_candidate")


class TestPeptideMapping(unittest.TestCase):
    """Tests for peptide-to-isoform mapping."""

    @classmethod
    def setUpClass(cls):
        cls.mod = load_module("map_sperm_proteomics_peptides_to_candidate_isoforms.py")

    def test_clean_peptide_sequence_removes_modification_text(self):
        """Modified MaxQuant sequences should be converted to bare amino acids."""
        self.assertEqual(self.mod.clean_peptide_sequence("_ACD(Oxidation)EF_"), "ACDEF")

    def test_peptide_mapping_finds_candidate_unique_peptide(self):
        """A peptide present only in the candidate protein should be reported."""
        candidates = pd.DataFrame(
            {
                "gene_symbol": ["GENE1"],
                "transcript_id_with_version": ["TX2.1"],
                "gencode_protein_id": ["P2.1"],
            }
        )
        proteins = pd.DataFrame(
            {
                "protein_id": ["P1.1", "P2.1"],
                "transcript_id": ["TX1.1", "TX2.1"],
                "gene_id": ["G1.1", "G1.1"],
                "transcript_name": ["GENE1-201", "GENE1-202"],
                "gene_symbol": ["GENE1", "GENE1"],
                "sequence": ["AAAAAAAKKKKK", "AAAAAAAPPPPP"],
            }
        )
        peptides = pd.DataFrame(
            {
                "clean_peptide_sequence": ["AAAAAAA", "PPPPP"],
                "n_evidence_rows": [1, 1],
            }
        )
        prepared = self.mod.prepare_candidate_table(
            dataframe=candidates,
            config=self.mod.Config(
                candidate_table=Path("dummy"),
                gencode_translations_fasta=Path("dummy"),
                peptide_table=Path("dummy"),
                out_dir=Path("dummy"),
                out_prefix="dummy",
                sequence_col=None,
                gene_col="gene_symbol",
                transcript_col="transcript_id_with_version",
                protein_col="gencode_protein_id",
                min_peptide_length=5,
                write_excel_outputs=False,
                log_level="CRITICAL",
                log_path=None,
            ),
        )
        logger = self.mod.setup_logging(log_level="CRITICAL", log_path=None)
        summary, evidence = self.mod.map_peptides_to_candidates(
            candidates=prepared,
            proteins=proteins,
            peptides=peptides,
            logger=logger,
        )
        self.assertEqual(int(summary.loc[0, "n_gene_unique_observed_peptides"]), 1)
        self.assertEqual(len(evidence), 2)


if __name__ == "__main__":
    unittest.main()
