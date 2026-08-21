# Reproduction run order

This document describes the **logical** run order. Some search-engine runs are external applications rather
than Python scripts, so exact parameter files/version records must accompany the public release.

## Stage 0 — external workflow outputs

Generate/obtain the nine-run MBR-OFF outputs from AP, FP, MM, and MQ using the same
isoform-inclusive UniProt FASTA and the manuscript's native workflow-specific filtering procedures.

Do not impose an artificial common scoring model across the four software packages.

## Stage 1 — theoretical reference space

1. `Step01_isoform_resolvability_v1.0.3.py`
2. `Step01B_extract_and_compare_observed_v1.1.1.py`
3. `Step01C_accessible_theoretical_ceiling_v1.0.0.py`
4. `Step01C2_consensus_accessibility_v1.0.0.py`
5. `Step01C3_resolved_vs_unresolved_v1.0.0.py`

Canonical Step01 command used for the frozen theoretical analysis:

```powershell
python .\Step01_isoform_resolvability_v1.0.3.py `
  --fasta .\uniprotkb_proteome_UP000005640_2026_08_04.fasta `
  --outdir .\Step01_results_v103 `
  --missed-cleavages 2 `
  --min-length 7 `
  --max-length 50 `
  --proteases trypsin lysc gluc_e chymo_fyw `
  --sequence-remap
```

## Stage 2 — common-reference mapping and observed evidence

Run the frozen common-FASTA mapping/normalization chain and validate the final workflow totals against
`config/expected_validation_anchors.yaml`.

## Stage 3 — ambiguity and replicate support

1. `Step03_ambiguity_decomposition_v1.0.1.py`
2. `Step04_replicate_threshold_robustness_v1.0.0.py`

Interpret downstream results only if the frozen common-reference and primary-peptide anchors match.

## Stage 4 — quantitative / SEPEP / reciprocal analyses

1. `STEP4A_QUANTITATIVE_BENCHMARK_V5_FINAL.py`
2. canonical `STEP4B_SEPEPQUANT_INSPIRED...` script selected from the real project tree
3. canonical reciprocal-event analysis script selected from the real project tree

The collector refuses to pick among multiple wildcard matches automatically.

## Stage 5 — shuffled-entrapment audit

1. `Step02A_REBUILD_ONECLICK_v1.0.5.ps1`
2. `Step02B_HOMEWORK_CHECK_v1.0.0.py`
3. `Step02C_AUDIT_UNMATCHED_v1.0.0.py`
4. `Step02D_PROTEIN_ENTRAPMENT_FDR_v1.0.3.py`
5. `Step02E_SALVAGE_MM_MQ_v1.0.0.py`
6. `Step02F_MQ_FINAL_MM_SEQUENCE_SALVAGE_v1.0.0.py`
7. canonical Step02G independent-validation script

Scientific distinction:
- MQ direct evaluation uses native protein-group FASTA-header information.
- MM sequence reconstruction is a diagnostic/post-hoc reconstruction because the relevant native MM
  protein table did not retain the entrapment marker.
- AP/MM bounds straddling 1% are classified as inconclusive, not failed.

## Stage 6 — native isoform/sibling audit

Run `AUDIT_NATIVE_SIBLING_ISOFORMS_v1.0.0.py`.

The headline comparison to validate is:
- 1,660 native base-accession families with >=2 explicit suffix-bearing entries
- 2 families with >=2 exact sibling isoforms after strict peptide-level resolution
