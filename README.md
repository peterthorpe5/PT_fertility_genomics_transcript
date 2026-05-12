# GTEx transcript-level isoform first-pass workflow

This mini-workflow is the transcript-level companion to the existing gene-level GTEx tissue-specificity analysis.

It is designed for the new working directory:

```bash
/home/pthorpe001/data/2026_sperm_Gates_transcript_level
```

The first-pass question is:

```text
Do putative target genes that are not sperm/testis-specific at gene level have individual transcript isoforms that are preferentially expressed or preferentially used in testis?
```

## Files expected in the GTEx directory

```bash
/home/pthorpe001/data/2026_sperm_Gates_transcript_level/GTEx/GTEx_Analysis_2025-08-22_v11_RSEMv1.3.3_transcripts_tpm.txt.gz
/home/pthorpe001/data/2026_sperm_Gates_transcript_level/GTEx/GTEx_Analysis_2025-08-22_v11_RSEMv1.3.3_transcripts_expected_count.txt.gz
/home/pthorpe001/data/2026_sperm_Gates_transcript_level/GTEx/GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_tpm.gct.gz
/home/pthorpe001/data/2026_sperm_Gates_transcript_level/GTEx/GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt
```

The expected-count file is optional for this first pass. If it is not present yet, the shell wrapper runs the TPM-based analysis and skips count extraction.

## Install/copy

Copy the folder into the project scripts area:

```bash
mkdir -p /home/pthorpe001/data/2026_sperm_Gates_transcript_level/scripts
cp -r gtex_transcript_isoform_first_pass \
  /home/pthorpe001/data/2026_sperm_Gates_transcript_level/scripts/
```

## Run unit tests

```bash
cd /home/pthorpe001/data/2026_sperm_Gates_transcript_level/scripts/gtex_transcript_isoform_first_pass
python -m unittest discover -s tests -v
```

## Run the selected-target first pass

```bash
cd /home/pthorpe001/data/2026_sperm_Gates_transcript_level/scripts/gtex_transcript_isoform_first_pass
bash run_first_pass_selected_targets.sh
```

## Direct Python command

```bash
python scripts/gtex_transcript_isoform_first_pass.py \
  --transcript_tpm_path /home/pthorpe001/data/2026_sperm_Gates_transcript_level/GTEx/GTEx_Analysis_2025-08-22_v11_RSEMv1.3.3_transcripts_tpm.txt.gz \
  --expected_count_path /home/pthorpe001/data/2026_sperm_Gates_transcript_level/GTEx/GTEx_Analysis_2025-08-22_v11_RSEMv1.3.3_transcripts_expected_count.txt.gz \
  --gene_tpm_gct_path /home/pthorpe001/data/2026_sperm_Gates_transcript_level/GTEx/GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_tpm.gct.gz \
  --sample_attributes_path /home/pthorpe001/data/2026_sperm_Gates_transcript_level/GTEx/GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt \
  --target_genes AFG2B CFAP99 SLC16A7 ABCG4 \
  --target_tissue Testis \
  --min_tpm_present 1.0 \
  --chunk_size 250 \
  --out_dir /home/pthorpe001/data/2026_sperm_Gates_transcript_level/results/01_gtex_transcript_isoform_first_pass \
  --out_prefix selected_targets_AFG2B_CFAP99_SLC16A7_ABCG4_gtex_v11 \
  --log_path /home/pthorpe001/data/2026_sperm_Gates_transcript_level/results/01_gtex_transcript_isoform_first_pass/selected_targets_AFG2B_CFAP99_SLC16A7_ABCG4_gtex_v11.log
```

## Outputs

The selected-target run writes only TSV or TSV.GZ files:

```text
selected_targets_AFG2B_CFAP99_SLC16A7_ABCG4_gtex_v11.transcript_subset_tpm_matrix.tsv.gz
selected_targets_AFG2B_CFAP99_SLC16A7_ABCG4_gtex_v11.transcript_sample_tpm_long.tsv.gz
selected_targets_AFG2B_CFAP99_SLC16A7_ABCG4_gtex_v11.transcript_tissue_median_tpm.tsv
selected_targets_AFG2B_CFAP99_SLC16A7_ABCG4_gtex_v11.transcript_isoform_usage_median_by_tissue.tsv
selected_targets_AFG2B_CFAP99_SLC16A7_ABCG4_gtex_v11.transcript_testis_isoform_summary.tsv
selected_targets_AFG2B_CFAP99_SLC16A7_ABCG4_gtex_v11.transcript_subset_expected_count_matrix.tsv.gz
```

The main file to inspect first is:

```text
selected_targets_AFG2B_CFAP99_SLC16A7_ABCG4_gtex_v11.transcript_testis_isoform_summary.tsv
```

Useful columns include:

```text
gene_symbol
transcript_id
target_median_tpm
max_non_target_tpm
log2_target_vs_max_non_target_tpm
target_median_isoform_usage
max_non_target_isoform_usage
log2_target_vs_max_non_target_isoform_usage
target_is_max_usage_tissue
target_tpm_rank_within_gene
target_usage_rank_within_gene
top_tissues_by_median_tpm
top_tissues_by_median_isoform_usage
```

## Interpretation

A potentially interesting transcript-level rescue signal would look like:

```text
Gene is not strongly testis-specific at gene level,
but one transcript has high target_median_tpm,
high target_median_isoform_usage,
target_is_max_usage_tissue = 1,
and positive log2_target_vs_max_non_target_isoform_usage.
```

The expected-count subset is saved for later formal transcript-usage modelling. The first-pass summaries are deliberately based on TPM and isoform fraction because they are easier to inspect and are closer to the immediate biological prioritisation question.
