"""Tests for isoform candidate gene-level evidence integration."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import pandas as pd

from scripts import integrate_isoform_candidates_with_gene_level_evidence as integ


class TestIntegrateIsoformCandidatesWithGeneEvidence(unittest.TestCase):
    """Unit tests for the evidence integration helpers."""

    def test_normalise_gene_symbol(self) -> None:
        """Gene symbols are stripped and upper-cased."""
        self.assertEqual(integ.normalise_gene_symbol(value=" slc16a7 "), "SLC16A7")
        self.assertEqual(integ.normalise_gene_symbol(value=None), "")

    def test_strip_ensembl_version(self) -> None:
        """Ensembl dot-version suffixes are removed."""
        self.assertEqual(
            integ.strip_ensembl_version(identifier="ENSG00000123456.7"),
            "ENSG00000123456",
        )

    def test_merge_by_gene_keys_uses_gene_id_fallback(self) -> None:
        """Gene ID fallback is used when symbols do not match."""
        logger = integ.setup_logging(log_level="CRITICAL", log_path=None)
        base = pd.DataFrame(
            {
                "gene_symbol": ["OLDNAME"],
                "gene_id": ["ENSG000001.4"],
                "join_gene_symbol": ["OLDNAME"],
                "join_gene_id": ["ENSG000001"],
            }
        )
        evidence = pd.DataFrame(
            {
                "join_gene_symbol": ["NEWNAME"],
                "join_gene_id": ["ENSG000001"],
                "evidence_flag": [True],
            }
        )
        merged, diagnostics = integ.merge_by_gene_keys(
            base=base,
            evidence=evidence,
            source_name="test_source",
            logger=logger,
        )
        self.assertTrue(bool(merged.loc[0, "evidence_flag"]))
        self.assertEqual(diagnostics["matched_by_gene_id_fallback"], 1)

    def test_derive_integrated_evidence_class(self) -> None:
        """Strong rescue candidates are classified with sperm/druggable support."""
        dataframe = pd.DataFrame(
            {
                "gene_symbol": ["GENE1"],
                "priority_tier": ["tier_1_protein_coding_rescue_candidate"],
                "gencode_is_protein_coding_transcript": [1],
                "is_broad_gene_isoform_rescue_candidate": [1],
                "biochem_sperm_present_any": [True],
                "biochem_public_prot_present_any": [False],
                "geneprio_ot_any_small_molecule_tractable": [True],
                "priority_score": [100.0],
            }
        )
        result = integ.derive_integrated_evidence(dataframe=dataframe)
        self.assertEqual(
            result.loc[0, "integrated_evidence_class"],
            "protein_coding_isoform_rescue_sperm_supported_druggable",
        )
        self.assertGreater(result.loc[0, "integrated_review_score"], 100.0)

    def test_end_to_end_small_run(self) -> None:
        """A tiny input set should produce integrated output tables."""
        logger = integ.setup_logging(log_level="CRITICAL", log_path=None)
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            isoform_path = tmp_path / "isoforms.tsv"
            biochem_path = tmp_path / "biochem.xlsx"
            priority_path = tmp_path / "priority.xlsx"
            tract_path = tmp_path / "tract.tsv"
            out_dir = tmp_path / "out"

            isoforms = pd.DataFrame(
                {
                    "gene_symbol": ["GENE1", "GENE2"],
                    "gene_id": ["ENSG000001", "ENSG000002"],
                    "transcript_id_with_version": ["ENST1.1", "ENST2.1"],
                    "target_median_tpm": [10.0, 5.0],
                    "target_median_isoform_usage": [0.6, 0.7],
                    "priority_tier": [
                        "tier_1_protein_coding_rescue_candidate",
                        "tier_4_other_testis_isoform_candidate",
                    ],
                    "priority_score": [100.0, 50.0],
                    "is_broad_gene_isoform_rescue_candidate": [1, 0],
                    "gencode_is_protein_coding_transcript": [1, 1],
                }
            )
            isoforms.to_csv(path_or_buf=isoform_path, sep="\t", index=False)

            biochem = pd.DataFrame(
                {
                    "gene_key": ["GENE1", "GENE2"],
                    "ensembl_gene_id": ["ENSG000001.1", "ENSG000002.1"],
                    "sperm_present_any": [True, False],
                    "public_prot_present_any": [True, False],
                    "public_proteomics_evidence_level": ["Strong", "None"],
                    "is_cell_surface_candidate": [False, False],
                    "is_membrane": [True, False],
                    "predicted_target_class": ["Enzyme", "Unknown"],
                }
            )
            with pd.ExcelWriter(path=biochem_path) as writer:
                biochem.to_excel(
                    excel_writer=writer,
                    sheet_name="Genes_Master",
                    index=False,
                )

            priority = pd.DataFrame(
                {
                    "gene_key": ["GENE1"],
                    "gene_id": ["ENSG000001.1"],
                    "candidate_druggable_sperm_protein": [True],
                    "ot_any_small_molecule_tractable": [True],
                    "ot_any_antibody_tractable": [False],
                    "ot_any_protac_tractable": [False],
                    "ot_any_tractable": [True],
                }
            )
            with pd.ExcelWriter(path=priority_path) as writer:
                priority.to_excel(
                    excel_writer=writer,
                    sheet_name="ranked",
                    index=False,
                )

            tract = pd.DataFrame(
                {
                    "gene_key": ["GENE1"],
                    "gene_id": ["ENSG000001.1"],
                    "gene_name": ["GENE1"],
                    "ot_any_small_molecule_tractable": [True],
                    "ot_any_antibody_tractable": [False],
                    "ot_any_protac_tractable": [False],
                    "ot_tractability_summary": ["SM:Druggable Family:True"],
                }
            )
            tract.to_csv(path_or_buf=tract_path, sep="\t", index=False)

            config = integ.Config(
                isoform_candidates_tsv=isoform_path,
                biochem_xlsx=biochem_path,
                biochem_sheet="Genes_Master",
                gene_priority_xlsx=priority_path,
                gene_priority_sheet="ranked",
                tractability_tsv=tract_path,
                out_dir=out_dir,
                out_prefix="test",
                selected_genes=("GENE1", "MISSING"),
                top_n_review=10,
                write_excel_outputs=False,
                log_path=None,
                log_level="CRITICAL",
            )
            integ.run(config=config, logger=logger)

            full_path = out_dir / "test.isoform_candidates_with_gene_level_evidence.tsv"
            top_path = out_dir / (
                "test.top_sperm_supported_druggable_isoform_rescue_candidates.tsv"
            )
            missing_path = out_dir / "test.selected_project_genes_missing.tsv"

            self.assertTrue(full_path.exists())
            self.assertTrue(top_path.exists())
            self.assertTrue(missing_path.exists())

            full = pd.read_csv(filepath_or_buffer=full_path, sep="\t")
            top = pd.read_csv(filepath_or_buffer=top_path, sep="\t")
            missing = pd.read_csv(filepath_or_buffer=missing_path, sep="\t")

            self.assertEqual(full.shape[0], 2)
            self.assertEqual(top.shape[0], 1)
            self.assertIn("MISSING", set(missing["gene_symbol"]))


if __name__ == "__main__":
    unittest.main()
