#!/usr/bin/env bash
#$ -cwd
#$ -j y
#$ -jc long
#$ -mods l_hard mfree 16G
#$ -adds l_hard h_vmem 16G
#$ -N isoform_priority
set -euo pipefail

WORKDIR="/home/pthorpe001/data/2026_sperm_Gates_transcript_level"
SCRIPT_DIR="${WORKDIR}/PT_fertility_genomics_transcript/scripts"
IN_DIR="${WORKDIR}/results/02_gtex_transcript_isoform_annotation_top3"
OUT_DIR="${WORKDIR}/results/04_gtex_testis_isoform_prioritisation_top3"
LOG_DIR="${OUT_DIR}/logs"

INPUT_TSV="${IN_DIR}/gtex_v11_transcriptome_testis_isoform_screen.best_candidate_isoform_per_gene.annotated.tsv"

mkdir -p "${OUT_DIR}" "${LOG_DIR}"

python "${SCRIPT_DIR}/prioritise_gtex_testis_isoform_candidates.py" \
  --input_tsv "${INPUT_TSV}" \
  --out_dir "${OUT_DIR}" \
  --out_prefix "gtex_v11_testis_isoform_prioritised" \
  --top_n_plot_genes 100 \
  --batch_size 25 \
  --selected_genes AFG2B CFAP99 SLC16A7 ABCG4 \
  --min_strong_tpm 1 \
  --min_strong_usage 0.25 \
  --min_strong_log2_usage_ratio 1 \
  --min_detection_fraction 0.20 \
  --log_path "${LOG_DIR}/gtex_v11_testis_isoform_prioritisation.log" \
  --log_level "INFO"
