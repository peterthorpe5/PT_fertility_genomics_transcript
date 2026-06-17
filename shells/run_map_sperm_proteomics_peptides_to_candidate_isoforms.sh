#!/usr/bin/env bash
#$ -cwd
#$ -V
#$ -N isoform_peptides
#$ -pe smp 2
#$ -j y
#$ -o logs/isoform_peptides.$JOB_ID.log

set -euo pipefail

WORKDIR="/home/${USER}/data/2026_sperm_Gates_transcript_level"
SCRIPT_DIR="${WORKDIR}/PT_fertility_genomics_transcript/scripts"
OUT_DIR="${WORKDIR}/results/07_isoform_sperm_peptide_support"
LOG_DIR="${WORKDIR}/logs"

CANDIDATE_TABLE="${WORKDIR}/results/06_isoform_submission_qc/gtex_v11_isoform_candidates_gene_level_evidence.strict_rebuilt.strict_top_supported_druggable_isoform_rescue_candidates.tsv"
TRANSLATIONS_FASTA="${WORKDIR}/genome_resources/gencode.v47.pc_translations.fa.gz"

# Use peptides.txt or evidence.txt after the MaxQuant output has finished downloading.
PEPTIDE_TABLE="/cluster/gjb_lab/pthorpe001/2026_sperm_Gates/PXD037531/peptides.txt"

mkdir -p "${OUT_DIR}" "${LOG_DIR}" logs

python "${SCRIPT_DIR}/map_sperm_proteomics_peptides_to_candidate_isoforms.py" \
    --candidate_table "${CANDIDATE_TABLE}" \
    --gencode_translations_fasta "${TRANSLATIONS_FASTA}" \
    --peptide_table "${PEPTIDE_TABLE}" \
    --out_dir "${OUT_DIR}" \
    --out_prefix "gtex_v11_strict_supported_druggable_isoform_rescue.sperm_peptide_support" \
    --log_path "${LOG_DIR}/gtex_v11_strict_supported_druggable_isoform_rescue.sperm_peptide_support.log" \
    --log_level "INFO"
