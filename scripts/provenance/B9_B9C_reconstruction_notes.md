# B9 / B9C reconstruction provenance

The publication script combines two final Stage-3 comparator analyses.

## B9 — direct Steele21 versus UCEI comparison

The preserved evidence map reports bidirectional incremental information on a common patient set. Adding UCEI to the Steele signature model increased concordance by approximately 0.0201 and improved fit with likelihood-ratio chi-square approximately 38.88 (P approximately 4.50e-10). Conversely, adding the Steele signature block to UCEI increased concordance by approximately 0.0229 (P approximately 0.0044). The interpretation is complementarity, not replacement.

## B9C — repeated nested benchmark

The final comparator panel contained 39 published predictors: Steele CN1-CN21, Drews CX1-CX17 and cDITHER. The principal benchmark used 10 repeats of five-fold outer validation with identical outer partitions for all competing models. Penalized Cox fitting used glmnet::cv.glmnet with family='cox', alpha=0.5, five-fold inner CV, Cox deviance, lambda.min, training-derived scaling, standardize=FALSE, unpenalized baseline covariates and penalized candidate predictors.

Locked mean C-indices were:

- AllSOTA39 + UCEI: 0.72610
- UCEI: 0.72458
- AllSOTA39: 0.71985
- cDITHER: 0.71809
- Drews17: 0.71565
- Steele21: 0.71255
- clinical baseline: 0.70468

Locked mean paired deltas were approximately:

- combined versus AllSOTA39: +0.006252, positive in 10/10 repeats
- combined versus UCEI: +0.001519, positive in 7/10
- UCEI versus AllSOTA39: +0.004733, positive in 10/10
- UCEI versus cDITHER: +0.006494, positive in 10/10
- UCEI versus Drews17: +0.008933, positive in 10/10
- UCEI versus Steele21: +0.01203, positive in 10/10

Selection stability in the combined model included CX13, UCEI and cDITHER in 50/50 folds; CX17 in 47/50; CN2 and CN9 in 44/50; CX12 in 34/50; CN5 in 32/50; and CX10 in 26/50.

## Recovery boundary

The exact historical B9C function names survive in the code-run archive (including B9C_make_folds, B9C_scale_fold, B9C_baseline_predict, B9C_penalized_predict and B9C_make_x), but their complete function bodies and the exact historical outer-fold seed/assignment generator are not available in the searchable transcript.

For that reason, the publication script requires an explicit historical fold-assignment table rather than silently inventing the outer partitions. The inner-fold generator in the reconstructed script is deterministic but should be checked against the historical B9C output tables. Exact public release should include the recovered historical outer-fold assignments, or the original B9C report object if it contains them.

The numerical audit in the script is intended to detect mismatch before public release.
