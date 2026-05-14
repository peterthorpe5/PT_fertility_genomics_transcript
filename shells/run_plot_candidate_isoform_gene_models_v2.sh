#!/usr/bin/env bash
#$ -cwd
#$ -j y
#$ -jc long
#$ -mods l_hard mfree 16G
#$ -adds l_hard h_vmem 16G
#$ -N isoform_fig_v2
set -euo pipefail

WORKDIR="/home/pthorpe001/data/2026_sperm_Gates_transcript_level"
SCRIPT_DIR="${WORKDIR}/PT_fertility_genomics_transcript/scripts"

ANNOT_DIR="${WORKDIR}/results/02_gtex_transcript_isoform_annotation_top3"
FIGURE_DIR="${WORKDIR}/results/03_gtex_isoform_gene_model_figures_top3_v2"
LOG_DIR="${FIGURE_DIR}/logs"

mkdir -p "${FIGURE_DIR}" "${LOG_DIR}"

python "${SCRIPT_DIR}/plot_candidate_isoform_gene_models_v2.py" \
  --features_tsv "${ANNOT_DIR}/gtex_v11_transcriptome_testis_isoform_screen.gencode_transcript_features.tsv.gz" \
  --candidate_tsv "${ANNOT_DIR}/gtex_v11_transcriptome_testis_isoform_screen.candidate_target_tissue_isoforms.annotated.tsv" \
  --out_dir "${FIGURE_DIR}" \
  --genes AFG2B CFAP99 SLC16A7 ABCG4 \
  --output_formats pdf svg \
  --out_suffix "candidate_isoform_gene_model_v2" \
  --max_transcripts 40 \
  --title_suffix "GTEx v11 testis-preferential isoform usage" \
  --log_path "${LOG_DIR}/plot_candidate_isoform_gene_models_v2.log" \
  --log_level "INFO"
