# Reproducibility / release checklist

## Canonical code
- [ ] Run `tools/build_public_repo.ps1`.
- [ ] Every required manifest row is `COPIED`.
- [ ] No required row is `MISSING` or `AMBIGUOUS`.
- [ ] Step02D is **v1.0.3**, not v1.0.0–v1.0.2.
- [ ] Step01 observed extraction is **v1.1.1**.
- [ ] Quantitative benchmark is **STEP4A_QUANTITATIVE_BENCHMARK_V5_FINAL.py**.
- [ ] Step02G source code has been recovered and included.

## Scientific validation
- [ ] AP/FP/MM/MQ common-reference peptide totals match 23,767 / 22,674 / 24,543 / 18,962.
- [ ] Primary discriminative totals match 180 / 101 / 149 / 70.
- [ ] Primary union is 297.
- [ ] Theoretical tryptic primary fraction reproduces 4.10%.
- [ ] Strict pooled accessibility sequence reproduces 245 -> 195 -> 21.
- [ ] Entrapment FDP bounds reproduce manuscript values.
- [ ] Reciprocal event sequence reproduces 342 -> 110 -> 44 -> 15.
- [ ] Native sibling audit reproduces 1,660 -> 2.

## Portability / privacy
- [ ] Run `python tools/validate_release.py --release-root <release>`.
- [ ] Remove or parameterize personal absolute paths.
- [ ] No local usernames remain in public scripts/configs.
- [ ] No raw MS data are committed.
- [ ] No proprietary binaries/JARs are redistributed without permission.
- [ ] External software download/version instructions are documented.

## Environments
- [ ] Python exact version recorded.
- [ ] Python package versions pinned.
- [ ] R version and `sessionInfo()` recorded.
- [ ] AP, FP/MSFragger/Philosopher, MM, MQ versions recorded.
- [ ] FDRBench version/source recorded.

## GitHub / Zenodo
- [ ] Choose a software license.
- [ ] Add final GitHub URL to `CODE_AVAILABILITY.md`.
- [ ] Create tag `v1.0.0` for the submitted manuscript.
- [ ] Archive the exact tag in Zenodo.
- [ ] Add Zenodo DOI.
- [ ] Link proteomics repository accession.
- [ ] Do not modify the archived release after submission; create a new version instead.
