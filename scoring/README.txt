UCEI-v2 Standalone Frozen Scorer
================================

Version: UCEI-v2.0

Purpose
-------
Computes the frozen UCEI-v2 score from six precomputed copy-number architecture features and burden_total.

Required input columns
----------------------
cancer
burden_total
jsd_neutral
ent_state_size
ent_chr_sd
ent_chr_mean
ent_ampclass
ent_transition

Scoring procedure
-----------------
1. Each raw feature is residualized against burden_total using frozen cancer-specific intercept and slope parameters.
2. Residuals are centered by the frozen residual mean and divided by the frozen residual SD.
3. The six standardized residual features are combined using the frozen PCA PC1 loadings.
4. The resulting PC1 score is transformed using the frozen score center and scale.

Important scope limitation
--------------------------
The scorer supports cancer labels represented in the frozen reference parameter table.
It does not silently extrapolate cancer-specific residualization parameters to unseen histologies.

Minimal R usage
---------------
source('score_ucei_v2.R')
x <- read.csv('test_input_20_samples.csv', check.names=FALSE)
x$UCEI_v2_z <- score_ucei_v2(x, '.')

Reproducibility test
--------------------
Compare calculated UCEI_v2_z values with test_expected_output.csv.
Expected maximum absolute difference is <= 1e-10.

Files
-----
score_ucei_v2.R
UCEI_v2_residualization_parameters.csv
UCEI_v2_PCA_loadings.csv
UCEI_v2_score_parameters.csv
test_input_20_samples.csv
test_expected_output.csv
data_dictionary.csv
sessionInfo.txt
MD5SUMS.txt
