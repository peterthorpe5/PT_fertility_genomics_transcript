#!/usr/bin/env bash
set -euo pipefail

WORKDIR="/home/pthorpe001/data/2026_sperm_Gates_transcript_level"
GTEX_DIR="${WORKDIR}/GTEx"
SCRIPT_DIR="${WORKDIR}/PT_fertility_genomics_transcript/scripts"
OUT_DIR="${WORKDIR}/results/01_gtex_transcript_isoform_first_pass"
LOG_DIR="${OUT_DIR}/logs"

mkdir -p "${OUT_DIR}" "${LOG_DIR}"

python "${SCRIPT_DIR}/gtex_transcriptome_isoform_screen.py" \
  --transcript_tpm "${GTEX_DIR}/GTEx_Analysis_2025-08-22_v11_RSEMv1.3.3_transcripts_tpm.txt.gz" \
  --gene_tpm_gct "${GTEX_DIR}/GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_tpm.gct.gz" \
  --sample_attributes "${GTEX_DIR}/GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt" \
  --out_dir "${OUT_DIR}" \
  --out_prefix "gtex_v11_transcriptome_testis_isoform_screen" \
  --target_tissue "Testis" \
  --float_dtype "float32" \
  --min_tpm_present 1 \
  --min_target_tpm_candidate 1 \
  --min_target_usage_candidate 0.25 \
  --min_log2_usage_ratio_candidate 1 \
  --max_gene_log2_ratio_for_rescue 1 \
  --log_path "${LOG_DIR}/gtex_v11_transcriptome_testis_isoform_screen.log" \
  --log_level "INFO"
