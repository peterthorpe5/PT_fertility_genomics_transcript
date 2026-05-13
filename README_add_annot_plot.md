# GTEx GENCODE annotation and candidate isoform gene-model plotting

This bundle adds the second stage of the transcript-level fertility genomics workflow:

1. annotate GTEx transcript isoform-screen outputs using the matching GENCODE GTF
2. make publication-style transcript/gene model figures for candidate isoforms

The scripts are designed for this project layout:

```bash
WORKDIR="/home/pthorpe001/data/2026_sperm_Gates_transcript_level"
SCRIPT_DIR="${WORKDIR}/PT_fertility_genomics_transcript/scripts"
```

## Recommended GTF

For GTEx v11 transcript quantification, use the release-matched GENCODE v47 primary assembly comprehensive GTF:

```bash
cd /home/pthorpe001/data/2026_sperm_Gates_transcript_level/genome_resources

wget -c \
  https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.primary_assembly.annotation.gtf.gz
```

Use `annotation.gtf.gz`, not `basic.annotation.gtf.gz`, because the basic file is a subset and may drop transcripts that are useful for isoform matching.

## Files in this bundle

```text
scripts/annotate_gtex_isoform_screen_with_gencode.py
scripts/plot_candidate_isoform_gene_models.py
tests/test_annotate_gtex_isoform_screen_with_gencode.py
tests/test_plot_candidate_isoform_gene_models.py
run_annotate_gtex_isoform_screen_with_gencode.sh
run_plot_candidate_isoform_gene_models.sh
```

## Install/copy into the project repo

```bash
WORKDIR="/home/pthorpe001/data/2026_sperm_Gates_transcript_level"
SCRIPT_DIR="${WORKDIR}/PT_fertility_genomics_transcript/scripts"

mkdir -p "${SCRIPT_DIR}"

cp scripts/annotate_gtex_isoform_screen_with_gencode.py "${SCRIPT_DIR}/"
cp scripts/plot_candidate_isoform_gene_models.py "${SCRIPT_DIR}/"
cp run_annotate_gtex_isoform_screen_with_gencode.sh "${SCRIPT_DIR}/"
cp run_plot_candidate_isoform_gene_models.sh "${SCRIPT_DIR}/"
```

## Run the unit tests

From the extracted bundle directory:

```bash
python -m unittest discover -s tests -v
```

Expected result:

```text
Ran 8 tests
OK
```

## Stage 2A: annotate the GTEx isoform screen

The wrapper assumes the top-3 rerun output folder:

```bash
/home/pthorpe001/data/2026_sperm_Gates_transcript_level/results/01_gtex_transcript_isoform_first_pass_top3
```

Run:

```bash
qsub "${SCRIPT_DIR}/run_annotate_gtex_isoform_screen_with_gencode.sh"
```

or run directly:

```bash
bash "${SCRIPT_DIR}/run_annotate_gtex_isoform_screen_with_gencode.sh"
```

Main output directory:

```text
/home/pthorpe001/data/2026_sperm_Gates_transcript_level/results/02_gtex_transcript_isoform_annotation_top3
```

Key outputs:

```text
gtex_v11_transcriptome_testis_isoform_screen.gencode_transcript_annotation.tsv.gz
gtex_v11_transcriptome_testis_isoform_screen.gencode_transcript_features.tsv.gz
gtex_v11_transcriptome_testis_isoform_screen.transcript_target_tissue_isoform_summary.annotated.tsv.gz
gtex_v11_transcriptome_testis_isoform_screen.candidate_target_tissue_isoforms.annotated.tsv
gtex_v11_transcriptome_testis_isoform_screen.best_candidate_isoform_per_gene.annotated.tsv
gtex_v11_transcriptome_testis_isoform_screen.candidate_target_tissue_isoforms.annotated.xlsx
gtex_v11_transcriptome_testis_isoform_screen.best_candidate_isoform_per_gene.annotated.xlsx
```

Excel output is enabled by default. Disable it with:

```bash
--no-write_excel_outputs
```

## Stage 2B: plot candidate isoform gene models

Run:

```bash
qsub "${SCRIPT_DIR}/run_plot_candidate_isoform_gene_models.sh"
```

or run directly:

```bash
bash "${SCRIPT_DIR}/run_plot_candidate_isoform_gene_models.sh"
```

Main output directory:

```text
/home/pthorpe001/data/2026_sperm_Gates_transcript_level/results/03_gtex_isoform_gene_model_figures_top3
```

By default, the wrapper plots:

```text
AFG2B
CFAP99
SLC16A7
ABCG4
```

Outputs are written as PDF and SVG files.

## Figure design

The gene model figure shows:

- one row per transcript isoform
- exons as boxes
- CDS regions as thicker boxes
- introns as connecting lines
- primary rank 1 candidate isoforms highlighted
- secondary rank 2-3 candidate isoforms highlighted
- non-candidate transcripts shown in pale background tracks
- a right-hand metric label for candidate isoforms

Candidate metrics shown on the right include:

```text
target_median_tpm
target_median_isoform_usage
log2_target_vs_max_non_target_isoform_usage
candidate_rank_tier
```

## Notes

The annotation script writes primary outputs as tab-separated files. Excel files are only convenience copies for manual browsing.

The plotting script uses the parsed GENCODE feature table rather than reparsing the GTF. This keeps the plotting stage faster and makes the gene-model figures reproducible from the same annotation table used for the candidate outputs.
