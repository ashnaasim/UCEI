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
- GitHub Actions reproducibility workflow completed successfully on R 4.5.2: standalone scoring/checksum verification PASS and Figure 1–6 source-data verification PASS.
- During code audit, grouped `data.table::.SD` scoping errors in the reconstructed DepMap and CN-aware transcription scripts were identified and corrected before release.
- Reader-facing documentation now uses **UCEI** while preserving internal `v2`/`frozen` labels only where required for traceability.

## Documented provenance limitations

Some publication-facing scripts are source-supported reconstructions because the exact standalone historical function bodies were not preserved. Corresponding files under `scripts/provenance/` state this explicitly and retain numerical verification targets. These scripts must not be described as byte-for-byte recovery of every historical workstation script.

The exact historical `glmnet` package version used for the nested Cox analysis was not found in the retained checksum-covered scorer session. The model specification is preserved, but the package version is not inferred.

Historical source-data manifests retain workstation paths and internal object names for provenance. Portable users should use repository-relative files and the public verifier instead.

## Audit status

**PASS — CLEARED FOR PUBLIC RELEASE (2026-10-04).**

The previously identified commit-metadata privacy blocker has been resolved. All rewritten commits on `main` use the GitHub noreply author/committer identity `aa246 <204020622+ashnaasim@users.noreply.github.com>`. The rewritten history was force-pushed to remote `main`, and the post-rewrite GitHub Actions reproducibility workflow completed successfully.

The repository file tree remains clean, with no accidental personal documents, credentials, unpublished third-party raw data, or unrelated private files identified in the audited release state.

## Public-release gate

The repository has passed the final pre-publication privacy, file-tree, and reproducibility checks.

Final privacy/file-tree check completed on 2026-10-04:

- Root contains only the intended release structures: `.github/`, `docs/`, `parameters/`, `reproducibility/`, `results/`, `scoring/`, `scripts/`, and `source_data/`, plus repository control files.
- No manuscript DOCX/PDF files, EndNote libraries, CVs, personal images, local caches, environment files, API credentials, passwords, tokens, or unrelated project files were found in the tracked tree.
- Historical `E:/cnv/...` paths remain only in provenance/source-manifest material and do not contain a Windows username or credential.
- TCGA identifiers in scoring/reference/source-data files are public de-identified study identifiers, not personal names.
- Git commit author and committer metadata were rewritten to use the GitHub noreply identity `aa246 <204020622+ashnaasim@users.noreply.github.com>`.
- The rewritten `main` history was force-pushed successfully.
- The post-rewrite GitHub Actions reproducibility workflow completed successfully.

Recommended final release steps:

1. Confirm that manuscript terminology and the Data Availability statement match the repository.
2. Change repository visibility to Public.
3. Verify anonymous public access to the repository and key reproducibility files.
4. Create a tagged release corresponding to the submitted manuscript version.
5. Archive that release in a persistent repository such as Zenodo and record the DOI.
6. Cite the public GitHub repository and archival DOI in the manuscript and submission system.

**Release status: CLEARED FOR PUBLIC RELEASE.**
