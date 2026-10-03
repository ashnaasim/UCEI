# B12C/B12C2 reconstruction provenance

`06_construction_ablation_and_higher_dimensional.R` reconstructs the final Stage-3 construction-sensitivity analyses from the preserved manuscript methods, final object/output schemas, and the audited evidence map.

## Directly supported by the preserved sources

The final B12C analysis compared:

- the complete one-dimensional UCEI;
- six leave-one-feature-out reconstructions;
- PCA without prior CNA-burden residualization;
- an equal-weight burden-residualized score;
- a six-residualized-component joint model.

The same repeated outer-validation framework was used rather than choosing a preferred variant on the full dataset.

The locked historical B12C mean repeat-level C-indices were:

- six residualized components jointly: 0.7303709
- leave out amplification-class entropy: 0.7268158
- leave out state-size entropy: 0.7261437
- complete UCEI: 0.7251728
- leave out transition entropy: 0.7248910
- leave out chromosome-mean architecture: 0.7245005
- leave out chromosome-SD architecture: 0.7238749
- leave out neutral-state JSD: 0.7222620
- PCA without burden residualization: 0.7201334
- equal-weight residualized score: 0.7173920

Historical mean differences versus complete UCEI were approximately +0.005198 for the six-component joint model, +0.001643 after omitting amplification-class entropy, +0.000971 after omitting state-size entropy, -0.005039 without burden residualization, and -0.007781 with equal weighting.

The B12C2 higher-dimensional comparison retained:

- AllSOTA39: 0.7215963
- AllSOTA39 + UCEI: 0.7272679
- AllSOTA39 + ARCH6: 0.7318942
- AllSOTA39 + ARCH5 without state-size: 0.7318520

These results support UCEI as a parsimonious one-dimensional representation, not as the maximally predictive architecture model.

## Reconstruction boundary

The standalone historical B12C/B12C2 R script was not preserved. Names of the original in-memory functions (`B12C_prepare_score`, `B12C_residual_matrices`, `B12C_cox_score`, `B12C_cox_block`) survive, but their complete function bodies are not present in the searchable transcript. Accordingly, the publication script reconstructs the documented training-only design and must be numerically checked against the retained historical B12C CSVs before public release.

In particular, implementation details for the equal-weight orientation and the exact Cox helper internals are reconstruction choices until historical-output verification is completed. Do not describe this file as a byte-for-byte historical script.
