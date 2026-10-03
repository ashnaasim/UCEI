# Recovery and archive manifest

This manifest records final publication assets recovered from the canonical UCEI project archive.

## Verified scoring assets on the original analysis workstation

Original directory:

`E:/cnv/UCEI_BENCHMARK_RESCUE/B2_CORE_UCEI_V2_BENCHMARK/B7_RUNTIME_REPRODUCIBILITY/UCEI_v2_STANDALONE_SCORER/`

Verified files:

- `score_ucei_v2.R`
- `UCEI_v2_PCA_loadings.csv`
- `UCEI_v2_residualization_parameters.csv`
- `UCEI_v2_score_parameters.csv`
- `README.txt`
- `data_dictionary.csv`
- `sessionInfo.txt`
- `MD5SUMS.txt`
- `standalone_reproduction_test.csv`
- `test_input_20_samples.csv`
- `test_expected_output.csv`

The canonical transcript contains the beginning of the standalone scorer and independently records its required inputs and scoring sequence. The exact local files should be copied into this repository before public release rather than silently reconstructed where the transcript is incomplete.

## Verified derivation/reference package

Original directory:

`E:/cnv/UCEI_BENCHMARK_RESCUE/UCEI_V2_FROZEN_REFERENCE/`

Verified files include:

- `UCEI_v2_frozen_burden_residualization_parameters.csv`
- `UCEI_v2_frozen_PCA_loadings.csv`
- `UCEI_v2_frozen_score_parameters.csv`
- `UCEI_v2_PCA_variance_explained.csv`
- `UCEI_v2_reference_scores.csv`
- `UCEI_v2_self_reproduction_by_cancer.csv`
- `UCEI_v2_self_reproduction_summary.csv`
- `UCEI_v2_FROZEN_PARAMETERS.rds`
- `report_UCEI_v2_freeze.rds`
- `sessionInfo.txt`

Reader-facing documentation refers to these as the retained UCEI derivation/scoring parameters; historical filenames are preserved only for provenance.

## Verified final figure-source package

Original directory:

`E:/cnv/UCEI_BIODATAMINING_SUBMISSION/FINAL_FIGURE_DATA/`

Verified outputs:

- `UCEI_FIGURES_1_TO_6_SOURCE_DATA_BUNDLE.rds`
- `UCEI_FIGURES_1_TO_6_SOURCE_MANIFEST.csv`
- `UCEI_FIGURES_1_TO_6_CRITICAL_AUDIT.csv`
- `UCEI_FIGURE_DATA_INDEX.rds`

The final figure-data harvest reported all critical source checks passing and did not refit analyses.

## Final-analysis script families to archive

The manuscript-facing script set will cover:

1. UCEI scoring and derivation parameters
2. construction stability / bootstrap / leave-one-cancer-out analyses
3. conventional CNA comparison
4. published CN-signature benchmarking and repeated nested validation
5. complete training-fold UCEI reconstruction
6. ablation, censoring-aware metrics, and held-out-cancer transportability
7. UCEI–ArchSig state analysis
8. latent architecture / transcriptional-manifestation analysis
9. reduced-feature external validation
10. DepMap genome-wide dependency analysis
11. CN-aware transcriptional decomposition
12. threshold/discordance and external CRISPR sensitivity analyses
13. proteomic integration and pathway analysis
14. single-cell and spatial sensitivity analyses
15. PRISM pharmacologic sensitivity analysis
16. figure-source-data assembly and audit

Only final supported branches should be promoted into `scripts/`; failed, exploratory, superseded, or manuscript-excluded branches remain in the historical archive.
