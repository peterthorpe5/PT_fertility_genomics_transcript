#!/usr/bin/env bash
#$ -cwd
#$ -V
#$ -N isoform_cds_class
#$ -pe smp 1
#$ -l h_vmem=16G
#$ -l h_rt=08:00:00
#$ -j y
#$ -o logs/isoform_cds_class.$JOB_ID.log

set -euo pipefail

WORKDIR="/home/${USER}/data/2026_sperm_Gates_transcript_level"
SCRIPT_DIR="${WORKDIR}/PT_fertility_genomics_transcript/scripts"
OUT_DIR="${WORKDIR}/results/06_isoform_submission_qc"
LOG_DIR="${WORKDIR}/logs"

CANDIDATE_TABLE="${OUT_DIR}/gtex_v11_isoform_candidates_gene_level_evidence.strict_rebuilt.strict_top_supported_druggable_isoform_rescue_candidates.tsv"
ANNOTATION_TSV="${WORKDIR}/results/02_gtex_transcript_isoform_annotation_top3/gtex_v11_transcriptome_testis_isoform_screen.gencode_transcript_annotation.tsv.gz"
FEATURES_TSV="${WORKDIR}/results/02_gtex_transcript_isoform_annotation_top3/gtex_v11_transcriptome_testis_isoform_screen.gencode_transcript_features.tsv.gz"

mkdir -p "${OUT_DIR}" "${LOG_DIR}" logs

python "${SCRIPT_DIR}/classify_candidate_isoform_cds_consequences.py" \
    --candidate_table "${CANDIDATE_TABLE}" \
    --transcript_annotation "${ANNOTATION_TSV}" \
    --transcript_features "${FEATURES_TSV}" \
    --out_dir "${OUT_DIR}" \
    --out_prefix "gtex_v11_strict_supported_druggable_isoform_rescue.cds_consequence" \
    --log_path "${LOG_DIR}/gtex_v11_strict_supported_druggable_isoform_rescue.cds_consequence.log" \
    --log_level "INFO"
