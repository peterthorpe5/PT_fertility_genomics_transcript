# Candidate isoform gene-model plotting wrapper

This bundle contains a shell wrapper and gene lists for generating manuscript and supplementary gene-model diagrams for the transcript-level sperm/fertility isoform project.

## Files

- `shells/run_plot_isoform_review_gene_models.sh`
  - SGE/qsub wrapper to plot all selected genes.
- `gene_lists/isoform_gene_model_review_genes.txt`
  - Default 18-gene review list.
- `gene_lists/isoform_unique_peptide_genes.txt`
  - The eight candidates with within-gene candidate-isoform-specific peptide evidence.
- `gene_lists/high_ranked_cds_review_genes.txt`
  - High-ranked CDS-changing review genes.
- `gene_lists/selected_project_example_genes.txt`
  - AFG2B and SLC16A7 selected-project examples.
- `scripts/plot_candidate_isoform_gene_models_v2.py`
  - Convenience copy of the patched plotting script. Use your repository copy if newer.

## Recommended target location

Copy into the project repository:

```bash
WORKDIR="/home/${USER}/data/2026_sperm_Gates_transcript_level"
SCRIPT_DIR="${WORKDIR}/PT_fertility_genomics_transcript/scripts"
SHELL_DIR="${WORKDIR}/PT_fertility_genomics_transcript/shells"

cp scripts/plot_candidate_isoform_gene_models_v2.py "${SCRIPT_DIR}/"
cp shells/run_plot_isoform_review_gene_models.sh "${SHELL_DIR}/"
```

## Run all suggested genes

```bash
WORKDIR="/home/${USER}/data/2026_sperm_Gates_transcript_level"
SHELL_DIR="${WORKDIR}/PT_fertility_genomics_transcript/shells"

qsub "${SHELL_DIR}/run_plot_isoform_review_gene_models.sh"
```

Outputs go to:

```bash
/home/${USER}/data/2026_sperm_Gates_transcript_level/results/08_candidate_gene_model_review_figures
```

Each gene should get PDF, SVG and PNG outputs.

## Run only a subset

For the eight isoform-specific peptide-supported candidates:

```bash
WORKDIR="/home/${USER}/data/2026_sperm_Gates_transcript_level"
GENE_FILE="${WORKDIR}/PT_fertility_genomics_transcript/gene_lists/isoform_unique_peptide_genes.txt" \
qsub "${WORKDIR}/PT_fertility_genomics_transcript/shells/run_plot_isoform_review_gene_models.sh"
```

For selected project examples only:

```bash
WORKDIR="/home/${USER}/data/2026_sperm_Gates_transcript_level"
GENE_FILE="${WORKDIR}/PT_fertility_genomics_transcript/gene_lists/selected_project_example_genes.txt" \
OUT_DIR="${WORKDIR}/results/08_candidate_gene_model_review_figures_selected_project" \
qsub "${WORKDIR}/PT_fertility_genomics_transcript/shells/run_plot_isoform_review_gene_models.sh"
```

## Notes

The wrapper uses the full annotated candidate table for plotting labels:

```text
results/02_gtex_transcript_isoform_annotation_top3/gtex_v11_transcriptome_testis_isoform_screen.candidate_target_tissue_isoforms.annotated.tsv
```

The strict/QC tables remain the source of manuscript counts. The full annotated candidate table is used here because it labels candidate isoforms and also supports selected-project examples such as AFG2B and SLC16A7.

