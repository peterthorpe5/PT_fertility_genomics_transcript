# GTEx transcript isoform submission QC scripts

This bundle adds the next review-critical scripts for the GTEx transcript-level sperm/testis isoform project.

The scripts are designed to resolve the main reviewer-facing weaknesses identified in the current manuscript draft:

1. audit the apparent 371 versus 612 candidate-row discrepancy after gene-level evidence integration;
2. rebuild strict integrated candidate tables using stable definitions;
3. classify whether candidate isoforms alter CDS structure or are mainly UTR/non-coding exon changes;
4. optionally map observed MaxQuant sperm peptides to candidate isoform protein sequences when peptide-level output is available.

All scripts use tab-separated output files. Excel review files are written by default unless `--no_write_excel_outputs` is supplied.

## Scripts

### `audit_isoform_evidence_integration_consistency.py`

Compares the strict pre-integration tier-1 table against the later integrated evidence output. It reports:

- row counts;
- unique gene/transcript counts;
- duplicate gene-transcript keys;
- post-integration rows not present in the strict tier-1 table;
- unexpected non-tier-1 rows in a purported protein-coding rescue output;
- conflicts between harmonised protein-coding flags and strict GENCODE-derived protein-coding flags.

Initial sandbox run on the uploaded results showed that the 612-row post-integration table did **not** contain duplicated gene-transcript keys. Instead, it contained 371 strict tier-1 rows plus 241 tier-2 rows that were reclassified by later harmonised protein-coding logic. This is definition drift, not a simple one-to-many merge duplication.

### `rebuild_strict_integrated_isoform_tables.py`

Rebuilds manuscript-safe integrated tables using strict definitions:

- original `priority_tier == tier_1_protein_coding_rescue_candidate`;
- strict GENCODE protein-coding transcript annotation;
- CDS support;
- broad-gene isoform-rescue status.

Initial sandbox run on the uploaded results recovered:

- 371 strict tier-1 protein-coding rescue rows;
- 222 strict tier-1 rows with sperm/protein support plus druggability/accessibility evidence.

These strict counts are safer for manuscript use than the previous 612-row harmonised table.

### `classify_candidate_isoform_cds_consequences.py`

Compares each candidate isoform to a reference transcript from the same gene and classifies whether the candidate changes CDS structure. Reference choice is:

1. MANE Select protein-coding transcript with CDS;
2. Ensembl canonical protein-coding transcript with CDS;
3. basic protein-coding transcript with the longest CDS;
4. any longest-CDS transcript.

If the candidate itself is the primary reference, the script uses the longest non-candidate CDS transcript as the final comparator where possible.

Initial sandbox run on the 222 strict supported/druggable isoform-rescue candidates found:

- 179 CDS-changing candidates;
- 42 UTR/non-coding exon-only changes;
- 1 CDS-identical/not structurally compelling candidate.

This is the biological climax that the manuscript currently needs.

### `map_sperm_proteomics_peptides_to_candidate_isoforms.py`

Optional peptide-level validation script for when MaxQuant peptide output is available. It accepts `peptides.txt` or `evidence.txt`, candidate isoform tables, and `gencode.v47.pc_translations.fa.gz`. It reports:

- observed peptides mapping to each candidate protein;
- peptides unique to the candidate isoform among proteins from the same gene;
- candidate-level peptide-support summaries.

This does **not** replace manual proteomics curation, but it provides a defensible first screen for isoform-specific sperm proteomic support.

## Suggested run order

```bash
WORKDIR="/home/${USER}/data/2026_sperm_Gates_transcript_level"
SCRIPT_DIR="${WORKDIR}/PT_fertility_genomics_transcript/scripts"
SHELL_DIR="${WORKDIR}/PT_fertility_genomics_transcript/shells"

mkdir -p "${SCRIPT_DIR}" "${SHELL_DIR}"

cp scripts/*.py "${SCRIPT_DIR}/"
cp shells/*.sh "${SHELL_DIR}/"

python -m unittest discover -s tests -v

qsub "${SHELL_DIR}/run_audit_isoform_evidence_integration_consistency.sh"
qsub "${SHELL_DIR}/run_rebuild_strict_integrated_isoform_tables.sh"
qsub "${SHELL_DIR}/run_classify_candidate_isoform_cds_consequences.sh"
```

After MaxQuant peptide outputs are available:

```bash
qsub "${SHELL_DIR}/run_map_sperm_proteomics_peptides_to_candidate_isoforms.sh"
```

## Manuscript implications

The current manuscript should be updated to remove the unresolved 371 versus 612 caveat once the audit/rebuild outputs are generated on the cluster. Based on the uploaded results, the safest current story is:

- 371 strict tier-1 protein-coding isoform-rescue candidates;
- 222 strict candidates retain gene-level sperm/protein support plus tractability/accessibility evidence;
- among those 222, an initial CDS-structure comparison classified 179 as CDS-changing relative to the selected reference transcript.

The phrase `candidate isoform has sperm proteomic support` should be avoided unless peptide-level mapping supports the exact isoform. Safer wording is: `candidate isoform occurs in a gene with sperm proteomic support`.
