# Code Availability

## Manuscript-ready statement â€” fill placeholders before submission

Custom analysis scripts, exact analysis parameters, and software-version records supporting this study
are available at **https://github.com/Supachaito/isoform-dda-benchmark**. A versioned archival release is available through
**<ZENODO_DOI_OR_ARCHIVE_URL>**. The repository contains the scripts used for theoretical
isoform-resolvability analysis, common-reference peptide remapping, replicate-support analyses,
quantitative/SEPEP analyses, shuffled-entrapment evaluation, and native isoform-output auditing.

Raw and processed mass-spectrometry data are available from **<PROTEOMICS_REPOSITORY>** under
accession **<ACCESSION>**.

## Before replacing the placeholders

- Tag the exact manuscript release, e.g. `v1.0.0`.
- Archive that tag in Zenodo and obtain the DOI.
- Run `tools/validate_release.py`.
- Confirm `COLLECTION_REPORT.tsv` contains no unresolved required scripts.
- Confirm no personal absolute paths or local-only identifiers remain in the public code.
- Add search-engine/software version numbers to `environment/software_versions.tsv`.
