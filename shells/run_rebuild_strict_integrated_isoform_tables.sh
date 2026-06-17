#!/usr/bin/env bash
#$ -cwd
#$ -V
#$ -N isoform_strict_rebuild
#$ -pe smp 2
#$ -j y
#$ -o logs/isoform_strict_rebuild.$JOB_ID.log

set -euo pipefail

WORKDIR="/home/${USER}/data/2026_sperm_Gates_transcript_level"
SCRIPT_DIR="${WORKDIR}/PT_fertility_genomics_transcript/scripts"
OUT_DIR="${WORKDIR}/results/06_isoform_submission_qc"
LOG_DIR="${WORKDIR}/logs"

INTEGRATED_TABLE="${WORKDIR}/results/05_gtex_isoform_gene_level_evidence_top3/gtex_v11_isoform_candidates_gene_level_evidence.tier1_protein_coding_rescue_with_gene_evidence.tsv"

mkdir -p "${OUT_DIR}" "${LOG_DIR}" logs

python "${SCRIPT_DIR}/rebuild_strict_integrated_isoform_tables.py" \
    --integrated_table "${INTEGRATED_TABLE}" \
    --out_dir "${OUT_DIR}" \
    --out_prefix "gtex_v11_isoform_candidates_gene_level_evidence.strict_rebuilt" \
    --selected_genes AFG2B CFAP99 SLC16A7 ABCG4 \
    --log_path "${LOG_DIR}/gtex_v11_isoform_candidates_gene_level_evidence.strict_rebuilt.log" \
    --log_level "INFO"
