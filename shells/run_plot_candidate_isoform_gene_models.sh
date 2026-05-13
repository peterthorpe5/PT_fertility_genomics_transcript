#!/usr/bin/env bash
#$ -cwd
#$ -j y
#$ -jc short
#$ -mods l_hard mfree 40G
#$ -adds l_hard h_vmem 40G
#$ -N plot_isoforms
set -euo pipefail

WORKDIR="/home/pthorpe001/data/2026_sperm_Gates_transcript_level"
SCRIPT_DIR="${WORKDIR}/PT_fertility_genomics_transcript/scripts"
ANNOT_DIR="${WORKDIR}/results/02_gtex_transcript_isoform_annotation_top3"
OUT_DIR="${WORKDIR}/results/03_gtex_isoform_gene_model_figures_top3"
LOG_DIR="${OUT_DIR}/logs"
PREFIX="gtex_v11_transcriptome_testis_isoform_screen"

mkdir -p "${OUT_DIR}" "${LOG_DIR}"

python "${SCRIPT_DIR}/plot_candidate_isoform_gene_models.py" \
  --features_tsv "${ANNOT_DIR}/${PREFIX}.gencode_transcript_features.tsv.gz" \
  --candidates_tsv "${ANNOT_DIR}/${PREFIX}.candidate_target_tissue_isoforms.annotated.tsv" \
  --out_dir "${OUT_DIR}" \
  --genes AFG2B CFAP99 SLC16A7 ABCG4 \
  --max_transcripts_per_gene 24 \
  --formats pdf svg \
  --title_suffix "GTEx v11 testis-preferential isoform usage" \
  --log_path "${LOG_DIR}/${PREFIX}.candidate_isoform_gene_models.log" \
  --log_level "INFO"
