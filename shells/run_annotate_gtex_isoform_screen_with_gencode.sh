#!/usr/bin/env bash
#$ -cwd
#$ -j y
#$ -jc long
#$ -mods l_hard mfree 80G
#$ -adds l_hard h_vmem 80G
#$ -N annotate_isoforms
set -euo pipefail

WORKDIR="/home/pthorpe001/data/2026_sperm_Gates_transcript_level"
SCRIPT_DIR="${WORKDIR}/PT_fertility_genomics_transcript/scripts"
GENOME_DIR="${WORKDIR}/genome_resources"
SCREEN_DIR="${WORKDIR}/results/01_gtex_transcript_isoform_first_pass_top3"
OUT_DIR="${WORKDIR}/results/02_gtex_transcript_isoform_annotation_top3"
LOG_DIR="${OUT_DIR}/logs"

GTF="${GENOME_DIR}/gencode.v47.primary_assembly.annotation.gtf.gz"
PREFIX="gtex_v11_transcriptome_testis_isoform_screen"

mkdir -p "${OUT_DIR}" "${LOG_DIR}"

python "${SCRIPT_DIR}/annotate_gtex_isoform_screen_with_gencode.py" \
  --gtf "${GTF}" \
  --screen_tsv "${SCREEN_DIR}/${PREFIX}.transcript_target_tissue_isoform_summary.tsv.gz" \
  --candidate_tsv "${SCREEN_DIR}/${PREFIX}.candidate_target_tissue_isoforms.tsv" \
  --best_per_gene_tsv "${SCREEN_DIR}/${PREFIX}.best_candidate_isoform_per_gene.tsv" \
  --out_dir "${OUT_DIR}" \
  --out_prefix "${PREFIX}" \
  --write_excel_outputs \
  --log_path "${LOG_DIR}/${PREFIX}.gencode_annotation.log" \
  --log_level "INFO"
