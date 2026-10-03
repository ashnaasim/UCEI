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

**HOLD FOR ONE PRIVACY FIX.** The repository file tree itself is clean and no accidental personal documents or raw private files were found. However, the Git commit history currently records a non-noreply personal author email in the commit metadata. Because that metadata becomes inspectable when a repository is public, repository visibility should remain private until the commit history is rewritten with a GitHub noreply email (or the owner explicitly accepts that exposure).

## Public-release gate

The repository is suitable for public release once the owner confirms that no additional private files, credentials, unpublished third-party raw data, or author-sensitive material have been added outside the audited paths.

Final privacy/file-tree check completed on 2026-10-04:

- Root contains only the intended release structures: `.github/`, `docs/`, `parameters/`, `reproducibility/`, `results/`, `scoring/`, `scripts/`, and `source_data/`, plus repository control files.
- No manuscript DOCX/PDF files, EndNote libraries, CVs, personal images, local caches, environment files, API credentials, passwords, tokens, or unrelated project files were found in the tracked tree.
- Historical `E:/cnv/...` paths remain only in provenance/source-manifest material and do not contain a Windows username or credential.
- TCGA identifiers in scoring/reference/source-data files are public de-identified study identifiers, not personal names.
- **Privacy blocker:** commit metadata uses a personal email rather than a GitHub noreply address.

Recommended final checks immediately before changing repository visibility:

1. Rewrite commit author/committer email metadata to a GitHub noreply address and force-push the rewritten `main` history.
2. Re-run the automated reproducibility workflow after that history rewrite.
3. Confirm that manuscript wording matches the repository terminology and Data Availability statement.
4. Only then change repository visibility to Public.
5. Create a tagged release or archive DOI after the final submission commit, if desired.

Do not change repository visibility until those owner-level checks are complete.