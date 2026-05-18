#!/usr/bin/env bash
#$ -cwd
#$ -j y
#$ -jc long
#$ -mods l_hard mfree 32G
#$ -adds l_hard h_vmem 32G
#$ -N iso_gene_evidence
set -euo pipefail

WORKDIR="/home/pthorpe001/data/2026_sperm_Gates_transcript_level"
SCRIPT_DIR="${WORKDIR}/PT_fertility_genomics_transcript/scripts"
EVIDENCE_DIR="${WORKDIR}/gene_level_evidence"
OUT_DIR="${WORKDIR}/results/05_gtex_isoform_gene_level_evidence_top3"
LOG_DIR="${OUT_DIR}/logs"

mkdir -p "${OUT_DIR}" "${LOG_DIR}"

ISOFORM_TSV="${WORKDIR}/results/04_gtex_testis_isoform_prioritisation_top3/gtex_v11_testis_isoform_prioritised.prioritised_all.tsv"
BIOCHEM_XLSX="${EVIDENCE_DIR}/SUMMARY_fertility_evidence_reduced.biochem.xlsx"
GENE_PRIORITY_XLSX="${EVIDENCE_DIR}/sperm_target_priorities_from_master.xlsx"
TRACTABILITY_TSV="${EVIDENCE_DIR}/gene_context_features_universe_plus_tractability.tsv.zip"

python "${SCRIPT_DIR}/integrate_isoform_candidates_with_gene_level_evidence.py" \
  --isoform_candidates_tsv "${ISOFORM_TSV}" \
  --biochem_xlsx "${BIOCHEM_XLSX}" \
  --biochem_sheet "Genes_Master" \
  --gene_priority_xlsx "${GENE_PRIORITY_XLSX}" \
  --gene_priority_sheet "ranked" \
  --tractability_tsv "${TRACTABILITY_TSV}" \
  --out_dir "${OUT_DIR}" \
  --out_prefix "gtex_v11_isoform_candidates_gene_level_evidence" \
  --selected_genes AFG2B CFAP99 SLC16A7 ABCG4 \
  --top_n_review 50 \
  --log_path "${LOG_DIR}/gtex_v11_isoform_candidates_gene_level_evidence.log" \
  --log_level "INFO"
