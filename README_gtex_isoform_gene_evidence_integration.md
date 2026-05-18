# GTEx isoform candidate integration with gene-level evidence

This bundle adds the next major evidence layer to the GTEx transcript-level
fertility genomics workflow. It joins transcript-level testis isoform candidates
to gene-level sperm RNA, sperm proteomics, biochemical accessibility, Open
Targets tractability, HPO, ClinVar and literature-derived evidence.

The main aim is to identify candidate isoforms that are not only
transcript-level rescue signals, but are also supported by sperm/gene-level and
therapeutic evidence.

## Files

```text
scripts/
  integrate_isoform_candidates_with_gene_level_evidence.py
shells/
  run_integrate_isoform_candidates_with_gene_level_evidence.sh
tests/
  test_integrate_isoform_candidates_with_gene_level_evidence.py
```

## Expected project layout

```bash
WORKDIR="/home/pthorpe001/data/2026_sperm_Gates_transcript_level"
SCRIPT_DIR="${WORKDIR}/PT_fertility_genomics_transcript/scripts"
SHELL_DIR="${WORKDIR}/PT_fertility_genomics_transcript/shells"
EVIDENCE_DIR="${WORKDIR}/gene_level_evidence"
```

Copy the uploaded gene-level evidence files into:

```bash
mkdir -p "${EVIDENCE_DIR}"

cp SUMMARY_fertility_evidence_reduced.biochem.xlsx \
  "${EVIDENCE_DIR}/"
cp sperm_target_priorities_from_master.xlsx \
  "${EVIDENCE_DIR}/"
cp gene_context_features_universe_plus_tractability.tsv.zip \
  "${EVIDENCE_DIR}/"
```

## Install into the repository

```bash
WORKDIR="/home/pthorpe001/data/2026_sperm_Gates_transcript_level"
SCRIPT_DIR="${WORKDIR}/PT_fertility_genomics_transcript/scripts"
SHELL_DIR="${WORKDIR}/PT_fertility_genomics_transcript/shells"

mkdir -p "${SCRIPT_DIR}" "${SHELL_DIR}"

cp scripts/integrate_isoform_candidates_with_gene_level_evidence.py \
  "${SCRIPT_DIR}/"
cp shells/run_integrate_isoform_candidates_with_gene_level_evidence.sh \
  "${SHELL_DIR}/"
```

## Run tests

From the unpacked bundle directory:

```bash
python -m unittest discover -s tests -v
```

Expected result:

```text
Ran 5 tests
OK
```

## Run the integration

```bash
qsub "${SHELL_DIR}/run_integrate_isoform_candidates_with_gene_level_evidence.sh"
```

or interactively:

```bash
bash "${SHELL_DIR}/run_integrate_isoform_candidates_with_gene_level_evidence.sh"
```

## Main outputs

The wrapper writes to:

```text
/home/pthorpe001/data/2026_sperm_Gates_transcript_level/results/05_gtex_isoform_gene_level_evidence_top3
```

Outputs:

```text
gtex_v11_isoform_candidates_gene_level_evidence.isoform_candidates_with_gene_level_evidence.tsv
gtex_v11_isoform_candidates_gene_level_evidence.tier1_protein_coding_rescue_with_gene_evidence.tsv
gtex_v11_isoform_candidates_gene_level_evidence.top_sperm_supported_druggable_isoform_rescue_candidates.tsv
gtex_v11_isoform_candidates_gene_level_evidence.selected_project_genes_isoform_gene_evidence.tsv
gtex_v11_isoform_candidates_gene_level_evidence.selected_project_genes_missing.tsv
gtex_v11_isoform_candidates_gene_level_evidence.summary_by_integrated_evidence_class.tsv
gtex_v11_isoform_candidates_gene_level_evidence.merge_diagnostics.tsv
gtex_v11_isoform_candidates_gene_level_evidence.genes_to_review_top_50.txt
gtex_v11_isoform_candidates_gene_level_evidence.tier1_genes_to_review_top_50.txt
```

Formatted Excel files are written by default for the main review tables. Disable
this with:

```bash
--no_write_excel_outputs
```

## Integrated evidence classes

The script creates a derived column called `integrated_evidence_class`.
Important classes include:

```text
protein_coding_isoform_rescue_sperm_supported_druggable
protein_coding_isoform_rescue_sperm_supported
protein_coding_isoform_rescue_druggable_no_sperm_support
protein_coding_isoform_rescue_limited_gene_support
non_coding_or_uncertain_isoform_rescue
other_isoform_candidate
```

The strongest immediate review group is:

```text
protein_coding_isoform_rescue_sperm_supported_druggable
```

This means:

```text
protein-coding testis isoform rescue candidate
+
sperm RNA or proteomics support
+
Open Targets tractability or biochemical accessibility signal
```

## Notes

The integrated review score is intended for triage, not as a formal statistical
model. It adds a transparent gene-level evidence score to the transcript-level
priority score where available. The score components are recorded in
`integrated_score_components`.
