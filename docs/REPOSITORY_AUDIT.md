# Repository release audit — 2026-10-04

This audit was performed before public release of the UCEI reproducibility repository.

## Passed checks

- Repository remains private during audit.
- MIT license is present.
- Raw third-party datasets are not intentionally redistributed as primary raw resources.
- Stable dataset releases/accessions are documented in `docs/DATASETS.md`.
- Checksum-covered historical scoring files are present.
- Standalone scorer test input, expected output, session information, and MD5 manifest are present.
- Standalone test reports maximum absolute error 7.7715611723761e-15 at tolerance 1e-10.
- Full-precision UCEI loadings are present.
- Historical burden-residualization and score parameters are present.
- Reference PCA variance and self-reproduction summaries are present.
- Figure 1–6 source-data bundle, manifest, index, and 19-check critical audit are present.
- Publication scripts 01–20 contain no hard-coded workstation drive paths.
- No publication script 01–20 calls `rm(list=ls())`, `setwd()`, or `install.packages()`.
- The historical workstation-specific Figure 1–6 harvester was moved out of `scripts/publication/` into `scripts/provenance/`.
- A portable Figure 1–6 verifier was added to `scripts/publication/`.
- Reader-facing scoring wrapper `scoring/score_ucei.R` was added without changing the checksum-covered historical scorer.
- `.gitattributes` now prevents line-ending conversion of checksum-covered historical scoring files.
- A portable `01_verify_ucei_scoring.R` script now checks both numerical reproduction and the historical MD5 manifest.
- During code audit, grouped `data.table::.SD` scoping errors in the reconstructed DepMap and CN-aware transcription scripts were identified and corrected before release.
- Reader-facing documentation now uses **UCEI** while preserving internal `v2`/`frozen` labels only where required for traceability.

## Documented provenance limitations

Some publication-facing scripts are source-supported reconstructions because the exact standalone historical function bodies were not preserved. Corresponding files under `scripts/provenance/` state this explicitly and retain numerical verification targets. These scripts must not be described as byte-for-byte recovery of every historical workstation script.

The exact historical `glmnet` package version used for the nested Cox analysis was not found in the retained checksum-covered scorer session. The model specification is preserved, but the package version is not inferred.

Historical source-data manifests retain workstation paths and internal object names for provenance. Portable users should use repository-relative files and the public verifier instead.

## Audit status

**PASS WITH DOCUMENTED PROVENANCE CAVEATS.** No remaining blocker was identified in the intended publication-facing paths. The repository should remain private until the owner performs the final visual file-tree/privacy check below.

## Public-release gate

The repository is suitable for public release once the owner confirms that no additional private files, credentials, unpublished third-party raw data, or author-sensitive material have been added outside the audited paths.

Recommended final manual checks immediately before changing repository visibility:

1. Open the GitHub file tree and confirm no accidental uploads outside the documented directories.
2. Confirm that `source_data/` contains only intended processed/derived source data.
3. Confirm that manuscript wording matches the repository terminology and Data Availability statement.
4. Create a tagged release or archive DOI after the final submission commit, if desired.

Do not change repository visibility until those owner-level checks are complete.