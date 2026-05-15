# GTEx testis isoform candidate prioritisation

This bundle adds a prioritisation step after the GTEx transcriptome-wide isoform screen and the GENCODE annotation step.

The script ranks annotated best-per-gene isoform candidates so that only the most biologically interesting genes need to be plotted or reviewed manually.

## Main script

```text
scripts/prioritise_gtex_testis_isoform_candidates.py
```

## Cluster wrapper

```text
shells/run_prioritise_gtex_testis_isoform_candidates.sh
```

The wrapper assumes the project layout:

```bash
WORKDIR="/home/pthorpe001/data/2026_sperm_Gates_transcript_level"
SCRIPT_DIR="${WORKDIR}/PT_fertility_genomics_transcript/scripts"
IN_DIR="${WORKDIR}/results/02_gtex_transcript_isoform_annotation_top3"
OUT_DIR="${WORKDIR}/results/04_gtex_testis_isoform_prioritisation_top3"
```

## Input

The expected input is the annotated best-per-gene candidate table:

```text
gtex_v11_transcriptome_testis_isoform_screen.best_candidate_isoform_per_gene.annotated.tsv
```

## Outputs

The main outputs are:

```text
gtex_v11_testis_isoform_prioritised.prioritised_all.tsv
gtex_v11_testis_isoform_prioritised.summary_by_priority_tier.tsv
gtex_v11_testis_isoform_prioritised.tier1_protein_coding_rescue_candidates.tsv
gtex_v11_testis_isoform_prioritised.tier2_rescue_non_coding_or_uncertain_candidates.tsv
gtex_v11_testis_isoform_prioritised.tier3_protein_coding_high_confidence_candidates.tsv
gtex_v11_testis_isoform_prioritised.all_rescue_candidates.tsv
gtex_v11_testis_isoform_prioritised.all_protein_coding_candidates.tsv
gtex_v11_testis_isoform_prioritised.selected_project_genes_priority_check.tsv
gtex_v11_testis_isoform_prioritised.selected_project_genes_missing_from_best.tsv
gtex_v11_testis_isoform_prioritised.genes_to_plot_top_100.txt
gtex_v11_testis_isoform_prioritised.genes_to_plot_tier1_protein_coding_rescue.txt
gtex_v11_testis_isoform_prioritised.genes_to_plot_all_rescue.txt
plot_gene_batches_top_ranked/genes_to_plot_batch_001.txt
```

Formatted Excel copies of the main browseable outputs are written by default. Disable these with:

```bash
--no_write_excel_outputs
```

## Priority tiers

The script uses these broad tiers:

```text
tier_1_protein_coding_rescue_candidate
    Broad-gene isoform-rescue candidate, high-confidence expression, and coding/CDS support.

tier_2_rescue_non_coding_or_uncertain_candidate
    Broad-gene isoform-rescue candidate with high-confidence expression, but not clearly protein-coding.

tier_3_protein_coding_high_confidence_candidate
    Strong protein-coding testis isoform candidate, but not classified as a broad-gene rescue case.

tier_4_other_testis_isoform_candidate
    Other candidates from the upstream isoform screen.

tier_5_other_candidate
    Fallback category for rows that are present but are not marked as candidates.
```

## Recommended run

Copy the files into the repo:

```bash
WORKDIR="/home/pthorpe001/data/2026_sperm_Gates_transcript_level"
SCRIPT_DIR="${WORKDIR}/PT_fertility_genomics_transcript/scripts"
SHELL_DIR="${WORKDIR}/PT_fertility_genomics_transcript/shells"

mkdir -p "${SCRIPT_DIR}" "${SHELL_DIR}"

cp scripts/prioritise_gtex_testis_isoform_candidates.py "${SCRIPT_DIR}/"
cp shells/run_prioritise_gtex_testis_isoform_candidates.sh "${SHELL_DIR}/"
```

Run tests from the unpacked bundle or from the repo root if the tests are copied there:

```bash
python -m unittest discover -s tests -v
```

Submit or run the wrapper:

```bash
qsub "${SHELL_DIR}/run_prioritise_gtex_testis_isoform_candidates.sh"
```

or:

```bash
bash "${SHELL_DIR}/run_prioritise_gtex_testis_isoform_candidates.sh"
```

## Suggested next use

Use the top gene list to drive the v2 gene-model plotting script:

```text
gtex_v11_testis_isoform_prioritised.genes_to_plot_top_100.txt
```

The batch files can be used later to split plotting into smaller jobs.
