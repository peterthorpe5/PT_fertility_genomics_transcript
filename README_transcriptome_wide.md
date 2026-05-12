# GTEx transcriptome-wide isoform screen

This folder contains a transcriptome-wide extension of the selected-target GTEx transcript isoform workflow.

The workflow loads the full GTEx v11 RSEM transcript TPM matrix, computes transcript median TPM by tissue, computes within-gene isoform usage by tissue, and ranks transcripts for target-tissue preferential usage.

## Files

```text
gtex_transcriptome_isoform_screen/
├── README.md
├── run_gtex_transcriptome_isoform_screen.sh
├── scripts/
│   └── gtex_transcriptome_isoform_screen.py
└── tests/
    └── test_gtex_transcriptome_isoform_screen.py
```

## Suggested project locations

For your current project, copy the Python script into the GitHub repo scripts folder:

```bash
WORKDIR="/home/pthorpe001/data/2026_sperm_Gates_transcript_level"
SCRIPT_DIR="${WORKDIR}/PT_fertility_genomics_transcript/scripts"

mkdir -p "${SCRIPT_DIR}"
cp scripts/gtex_transcriptome_isoform_screen.py "${SCRIPT_DIR}/"
```

The outputs are written to:

```bash
OUT_DIR="${WORKDIR}/results/01_gtex_transcript_isoform_first_pass"
```

## Run unit tests

From this folder:

```bash
python -m unittest discover -s tests -v
```

## Run the transcriptome-wide screen

After copying the Python script into the GitHub repo scripts folder, run:

```bash
bash run_gtex_transcriptome_isoform_screen.sh
```

The wrapper assumes this layout:

```bash
WORKDIR="/home/pthorpe001/data/2026_sperm_Gates_transcript_level"
GTEX_DIR="${WORKDIR}/GTEx"
SCRIPT_DIR="${WORKDIR}/PT_fertility_genomics_transcript/scripts"
OUT_DIR="${WORKDIR}/results/01_gtex_transcript_isoform_first_pass"
```

## Main outputs

```text
gtex_v11_transcriptome_testis_isoform_screen.transcript_tissue_median_tpm.tsv.gz
gtex_v11_transcriptome_testis_isoform_screen.transcript_isoform_usage_median_by_tissue.tsv.gz
gtex_v11_transcriptome_testis_isoform_screen.gene_tissue_median_tpm.tsv.gz
gtex_v11_transcriptome_testis_isoform_screen.gene_target_tissue_summary.tsv.gz
gtex_v11_transcriptome_testis_isoform_screen.transcript_target_tissue_isoform_summary.tsv.gz
gtex_v11_transcriptome_testis_isoform_screen.candidate_target_tissue_isoforms.tsv
gtex_v11_transcriptome_testis_isoform_screen.best_candidate_isoform_per_gene.tsv
```

## Candidate definition

By default, a candidate transcript must satisfy:

```text
target_median_tpm >= 1
target_median_isoform_usage >= 0.25
log2_target_vs_max_non_target_isoform_usage >= 1
target_usage_rank_within_gene == 1
target_is_max_usage_tissue == 1
```

A broader-gene isoform-rescue candidate additionally has:

```text
log2_target_vs_max_non_target_gene_tpm <= 1
```

This is intended to flag genes where gene-level expression is not strongly testis-selective but one transcript isoform may be preferentially used in testis.

## Notes

The script deliberately loads the full transcript TPM matrix into memory. It uses `float32` by default for expression arrays to reduce RAM use while keeping the code simple and transparent.

All tabular outputs are TSV or TSV.GZ files, not comma-separated files.
