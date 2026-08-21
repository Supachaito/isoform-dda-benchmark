# Isoform-aware DDA proteomics benchmark — reproducibility release

This repository is the public-code release structure for a four-workflow DDA benchmark comparing
AlphaPept (AP), FragPipe (FP), MetaMorpheus (MM), and MaxQuant (MQ) under a common
isoform-aware reference framework.

## Scientific principle

The benchmark separates three claims that should not be conflated:

1. **Structural specificity** — what protein entries a measured peptide sequence is compatible with.
2. **Replicate support** — how consistently that evidence is observed across biological replicates.
3. **Quantitative reproducibility** — how consistently peptide/SEPEP abundance is measured.

Native protein or isoform reporting is therefore not treated as equivalent to exact isoform attribution.

## What this release contains

The release is organized around the manuscript analysis path:

- `scripts/01_theoretical/` — theoretical isoform resolvability and accessibility ceiling
- `scripts/02_observed_mapping/` — extraction, I/L-equivalent common-reference remapping, and classification
- `scripts/03_ambiguity_replicates/` — ambiguity decomposition and replicate-threshold robustness
- `scripts/04_quantitative_sepep/` — quantitative benchmark, SEPEP analyses, reciprocal-event analyses
- `scripts/05_entrapment_fdr/` — protein-level shuffled-entrapment QC and audit
- `scripts/06_native_isoform_audit/` — native suffix-bearing accession and sibling-family audit
- `scripts/supplementary/` — analyses not required for the main result
- `config/` — frozen analysis parameters and validation anchors
- `environment/` — Python/R package records and software-version template
- `docs/` — run order, release checklist, and canonical-script manifest
- `tools/` — deterministic collector and release validator

## Important release rule

The public repository must contain the **exact canonical scripts that produced the final manuscript
results**. Older development versions are not substituted. The collector in `tools/build_public_repo.ps1`
copies only the requested canonical filenames and refuses ambiguous matches.

## Core verified analysis settings

- Primary analyses: MBR OFF
- Cell lines: C33A, HeLa, SiHa
- Replicates: 3 per cell line (9 DDA runs total)
- Common isoform-inclusive UniProt FASTA: 169,637 entries
- Explicitly suffixed isoform entries: 22,131
- I/L-equivalent common-reference remapping: enabled
- Theoretical digestion:
  - missed cleavages: 2
  - peptide length: 7–50 aa
  - proteases: trypsin, LysC, GluC-E, chymotrypsin F/Y/W
  - conservative sequence remapping: enabled
- Entrapment database: 1:1 protein-level shuffled entrapment
- FDRBench seed: 20260812
- Nominal protein FDR: 0.01
- Strong reciprocal SEPEP event: one unit with log2FC >= +1 and another <= -1

See `config/analysis_parameters.yaml` and `config/expected_validation_anchors.yaml`.

## Reproducibility philosophy

Validation anchors are treated as hard checks. A downstream script should not be interpreted if it fails
to reproduce the frozen upstream counts expected by the final manuscript. The release therefore includes
known expected peptide counts and key summary values as explicit audit targets.

## Data availability

Large raw and processed proteomics data should **not** be committed to GitHub. Deposit them in an
appropriate proteomics repository and provide the accession in `CODE_AVAILABILITY.md`.

## License

No license is selected automatically. Choose a license before public release; see
`LICENSE_CHOOSE_BEFORE_RELEASE.txt`.
