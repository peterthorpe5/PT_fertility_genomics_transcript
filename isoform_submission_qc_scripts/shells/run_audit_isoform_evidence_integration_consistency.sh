#!/usr/bin/env bash
#$ -cwd
#$ -V
#$ -N isoform_integration_audit
#$ -pe smp 1
#$ -l h_vmem=8G
#$ -l h_rt=04:00:00
#$ -j y
#$ -o logs/isoform_integration_audit.$JOB_ID.log

set -euo pipefail

WORKDIR="/home/${USER}/data/2026_sperm_Gates_transcript_level"
SCRIPT_DIR="${WORKDIR}/PT_fertility_genomics_transcript/scripts"
OUT_DIR="${WORKDIR}/results/06_isoform_submission_qc"
LOG_DIR="${WORKDIR}/logs"

PRE_TIER1="${WORKDIR}/results/04_gtex_testis_isoform_prioritisation_top3/gtex_v11_testis_isoform_prioritised.tier1_protein_coding_rescue_candidates.tsv"
POST_TABLE="${WORKDIR}/results/05_gtex_isoform_gene_level_evidence_top3/gtex_v11_isoform_candidates_gene_level_evidence.tier1_protein_coding_rescue_with_gene_evidence.tsv"

mkdir -p "${OUT_DIR}" "${LOG_DIR}" logs

python "${SCRIPT_DIR}/audit_isoform_evidence_integration_consistency.py" \
    --pre_integration_tier1 "${PRE_TIER1}" \
    --post_integration_table "${POST_TABLE}" \
    --out_dir "${OUT_DIR}" \
    --out_prefix "gtex_v11_isoform_gene_evidence_integration_audit" \
    --log_path "${LOG_DIR}/gtex_v11_isoform_gene_evidence_integration_audit.log" \
    --log_level "INFO"
