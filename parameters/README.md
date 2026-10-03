# UCEI parameters

This directory separates reader-facing parameter summaries from exact historical derivation/reference files.

## Reader-facing file

- `UCEI_PCA_loadings.csv` — full-precision six-feature PC1 loadings with descriptive feature labels. Historical source-column labels are retained in the `feature_code` field for traceability.

## Exact reference files

`reference/` preserves the retained derivation outputs used for reconstruction, including:

- cancer-specific burden-residualization parameters;
- full-precision PC1 loadings;
- score center/scale and orientation parameters;
- PCA variance explained;
- reference scores;
- numerical self-reproduction summary.

Historical filenames contain internal labels such as `v2` and `frozen`. These labels are retained solely for provenance. The manuscript and reader-facing method are referred to as **UCEI**.

The standalone scorer uses checksum-covered copies of the necessary parameter files under `../scoring/`.
