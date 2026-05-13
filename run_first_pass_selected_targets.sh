#!/usr/bin/env bash
#$ -cwd
#$ -j y
#$ -jc long
#$ -mods l_hard mfree 100G
#$ -adds l_hard h_vmem 100G
#$ -N focused_first_pass

set -euo pipefail

WORKDIR="/home/pthorpe001/data/2026_sperm_Gates_transcript_level"
GTEX_DIR="${WORKDIR}/GTEx"
OUT_DIR="${WORKDIR}/results/01_gtex_transcript_isoform_first_pass_top3"
SCRIPT_DIR="${WORKDIR}/PT_fertility_genomics_transcript/scripts"


mkdir -p "${OUT_DIR}"

TPM_FILE="${GTEX_DIR}/GTEx_Analysis_2025-08-22_v11_RSEMv1.3.3_transcripts_tpm.txt.gz"
COUNT_FILE="${GTEX_DIR}/GTEx_Analysis_2025-08-22_v11_RSEMv1.3.3_transcripts_expected_count.txt.gz"
GENE_GCT="${GTEX_DIR}/GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_tpm.gct.gz"
SAMPLE_ATTR="${GTEX_DIR}/GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt"

COUNT_ARGS=()
if [[ -s "${COUNT_FILE}" ]] && gzip -t "${COUNT_FILE}" 2>/dev/null; then
    COUNT_ARGS=(--expected_count_path "${COUNT_FILE}")
else
    echo "Expected-count file not found, incomplete, or still downloading; running TPM-only first pass" >&2
fi

python "${SCRIPT_DIR}/gtex_transcript_isoform_first_pass.py" \
    --transcript_tpm_path "${TPM_FILE}" \
    "${COUNT_ARGS[@]}" \
    --gene_tpm_gct_path "${GENE_GCT}" \
    --sample_attributes_path "${SAMPLE_ATTR}" \
    --target_genes AFG2B CFAP99 SLC16A7 ABCG4 \
    --target_tissue Testis \
    --min_tpm_present 1.0 \
    --chunk_size 250 \
    --out_dir "${OUT_DIR}" \
    --out_prefix selected_targets_AFG2B_CFAP99_SLC16A7_ABCG4_gtex_v11 \
    --log_path "${OUT_DIR}/selected_targets_AFG2B_CFAP99_SLC16A7_ABCG4_gtex_v11.log" \
    --log_level INFO
