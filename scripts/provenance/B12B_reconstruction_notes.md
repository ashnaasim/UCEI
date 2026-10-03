# B12B reconstruction provenance

The file `05_training_fold_ucei_reconstruction.R` is a publication-facing reconstruction of the final executed B12B Stage-3 analysis. No standalone historical B12B `.R` file survived on disk, so the script was reconstructed from the preserved code-run record, manuscript methods, final output schemas, and the locked evidence map.

## Source-supported final settings

- Analysis cohort: 1,574 samples, 439 PFI events, 9 cancers.
- Six raw architecture features: `jsd_neutral`, `ent_state_size`, `ent_chr_sd`, `ent_chr_mean`, `ent_ampclass`, `ent_transition`.
- Comparator panel: 21 Steele CN signatures (CN1-CN21), 17 Drews signatures (CX1-CX17), and cDITHER = 39 predictors.
- Complete UCEI derivation was re-estimated inside every outer training fold: cancer-specific burden regression, residual mean/SD estimation, standardization, PCA, training-derived PC1 orientation, and score centering/scaling.
- Outer validation: 10 repeats x 5 folds.
- Fold construction: cancer-by-event stratified.
- Outer seeds: `20260813 + repeat`.
- Inner seeds: `20260813 + repeat*100 + fold`.
- Inner CV: 5 folds.
- Penalized model: `glmnet::cv.glmnet`, Cox family, alpha 0.5, Cox deviance, `lambda.min`, `standardize=FALSE`, `nlambda=100`.
- Baseline covariates were unpenalized; candidate predictors were penalized.
- Candidate scaling was learned from the outer training partition.
- Nonzero coefficient threshold for selection: absolute coefficient > 1e-10.

## Locked historical result targets

Mean repeat-level C-index:

- baseline: 0.7082636
- training-fold UCEI: 0.7248864
- AllSOTA39: 0.7215963
- AllSOTA39 + training-fold UCEI: 0.7272679

All four recorded paired comparisons were positive in 10/10 repeats. UCEI was selected in 50/50 folds both alone and in the combined model.

## Important denominator note

The preserved B12B fold audit contains approximately 1,574 observations per repeat (for example 1,253 training + 321 test in fold 1) and 439 total events, not 3,243/781. The 3,243/781 denominator belongs to the broader clinical benchmark/held-out-cancer analyses. Repository documentation therefore retains the B12B-specific denominator rather than silently substituting the broader benchmark denominator.

This script should be validated against the retained historical B12B CSV outputs before the repository is made public.
