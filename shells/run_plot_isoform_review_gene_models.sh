#!/usr/bin/env bash
#$ -cwd
#$ -V
#$ -N plot_isoform_gene_models
#$ -pe smp 2
#$ -l h_rt=04:00:00
#$ -l h_vmem=8G
#$ -o logs/plot_isoform_gene_models.o$JOB_ID
#$ -e logs/plot_isoform_gene_models.e$JOB_ID

set -euo pipefail

# -----------------------------------------------------------------------------
# Plot candidate isoform gene models for manuscript/supplementary review.
#
# This wrapper calls the existing Python plotting script. It does not contain
# inline Python. It writes one PDF, SVG and PNG per gene and, if ImageMagick is
# available, a PNG contact sheet for rapid visual review.
#
# Recommended use:
#   qsub /home/${USER}/data/2026_sperm_Gates_transcript_level/\
#     PT_fertility_genomics_transcript/shells/run_plot_isoform_review_gene_models.sh
#
# Optional overrides:
#   WORKDIR=/path/to/project qsub shells/run_plot_isoform_review_gene_models.sh
#   GENE_FILE=/path/to/gene_list.txt qsub shells/run_plot_isoform_review_gene_models.sh
#   CANDIDATE_TSV=/path/to/candidates.tsv qsub shells/run_plot_isoform_review_gene_models.sh
#   PLOT_SCRIPT=/path/to/plot_candidate_isoform_gene_models_v2.py qsub shells/run_plot_isoform_review_gene_models.sh
# -----------------------------------------------------------------------------

WORKDIR="${WORKDIR:-/home/${USER}/data/2026_sperm_Gates_transcript_level}"
SCRIPT_DIR="${SCRIPT_DIR:-${WORKDIR}/PT_fertility_genomics_transcript/scripts}"
SHELL_DIR="${SHELL_DIR:-${WORKDIR}/PT_fertility_genomics_transcript/shells}"
RESULTS_DIR="${RESULTS_DIR:-${WORKDIR}/results}"
LOG_DIR="${LOG_DIR:-${WORKDIR}/logs}"

mkdir -p "${LOG_DIR}"

# Input files from the existing transcript-level analysis.
FEATURES_TSV="${FEATURES_TSV:-${RESULTS_DIR}/02_gtex_transcript_isoform_annotation_top3/gtex_v11_transcriptome_testis_isoform_screen.gencode_transcript_features.tsv.gz}"

# Use the full annotated candidate table for plotting. This highlights candidate
# isoforms and also supports selected-project genes such as AFG2B/SLC16A7.
# The strict/QC files should remain the source of manuscript counts; this table
# is used only for gene-model visualisation labels.
CANDIDATE_TSV="${CANDIDATE_TSV:-${RESULTS_DIR}/02_gtex_transcript_isoform_annotation_top3/gtex_v11_transcriptome_testis_isoform_screen.candidate_target_tissue_isoforms.annotated.tsv}"

OUT_DIR="${OUT_DIR:-${RESULTS_DIR}/08_candidate_gene_model_review_figures}"
GENE_FILE="${GENE_FILE:-${OUT_DIR}/isoform_gene_model_review_genes.txt}"

mkdir -p "${OUT_DIR}"

# Prefer the current repository script. Fall back to the patched script name used
# during development if that is what exists in the scripts directory.
PLOT_SCRIPT="${PLOT_SCRIPT:-${SCRIPT_DIR}/plot_candidate_isoform_gene_models_v2.py}"
if [[ ! -f "${PLOT_SCRIPT}" ]] && [[ -f "${SCRIPT_DIR}/plot_candidate_isoform_gene_models_v2_patched2.py" ]]; then
    PLOT_SCRIPT="${SCRIPT_DIR}/plot_candidate_isoform_gene_models_v2_patched2.py"
fi

if [[ ! -f "${PLOT_SCRIPT}" ]]; then
    echo "ERROR: plotting script not found: ${PLOT_SCRIPT}" >&2
    echo "Copy plot_candidate_isoform_gene_models_v2.py into ${SCRIPT_DIR}, or set PLOT_SCRIPT=/path/to/script.py" >&2
    exit 1
fi

if [[ ! -f "${FEATURES_TSV}" ]]; then
    echo "ERROR: feature table missing: ${FEATURES_TSV}" >&2
    exit 1
fi

if [[ ! -f "${CANDIDATE_TSV}" ]]; then
    echo "ERROR: candidate table missing: ${CANDIDATE_TSV}" >&2
    exit 1
fi

# Create the default gene list unless the user supplied a GENE_FILE path that
# already exists.
if [[ ! -f "${GENE_FILE}" ]]; then
    cat > "${GENE_FILE}" <<'GENES'
C2orf74
PPP1CC
NDRG3
BRD8
PRDX4
NXT2
PFKP
CAPZB
IQGAP2
SPATA20
ACSL1
PSMA6
RIOK3
ANO1
HK1
CPT1B
AFG2B
SLC16A7
GENES
fi

printf "Project directory: %s\n" "${WORKDIR}"
printf "Plotting script: %s\n" "${PLOT_SCRIPT}"
printf "Features table: %s\n" "${FEATURES_TSV}"
printf "Candidate table: %s\n" "${CANDIDATE_TSV}"
printf "Gene list: %s\n" "${GENE_FILE}"
printf "Output directory: %s\n" "${OUT_DIR}"
printf "Gene count: "
grep -cv '^\s*$' "${GENE_FILE}" || true

python "${PLOT_SCRIPT}" \
    --features_tsv "${FEATURES_TSV}" \
    --candidate_tsv "${CANDIDATE_TSV}" \
    --gene_file "${GENE_FILE}" \
    --out_dir "${OUT_DIR}" \
    --output_formats pdf svg png \
    --out_suffix "candidate_isoform_gene_model_v2_review" \
    --max_transcripts 18 \
    --title_suffix "GTEx v11 testis-preferential isoform usage" \
    --dpi 300 \
    --log_path "${LOG_DIR}/plot_isoform_review_gene_models.log" \
    --log_level INFO

# Optional contact sheet for quick visual triage.
if command -v montage >/dev/null 2>&1; then
    montage \
        "${OUT_DIR}"/*.candidate_isoform_gene_model_v2_review.png \
        -tile 3x \
        -geometry +24+24 \
        -background white \
        "${OUT_DIR}/isoform_gene_model_review_contact_sheet.png"
    printf "Contact sheet written to: %s\n" "${OUT_DIR}/isoform_gene_model_review_contact_sheet.png"
else
    echo "ImageMagick montage was not found; skipping contact-sheet creation."
fi

printf "Finished plotting candidate isoform gene models.\n"
