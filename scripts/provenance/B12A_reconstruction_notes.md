# B12A reconstruction provenance

`02_ucei_construction_stability.R` was reconstructed from the final executed Stage-3 B12A record because no standalone B12A `.R` file survived on disk.

Recovered final settings include: six architecture features; cancer-specific residualization against `burden_total`; residual centering/scaling within each refit; PCA with `prcomp(..., center = FALSE, scale. = FALSE)`; 200 cancer-stratified bootstrap reconstructions; leave-one-cancer-out reconstruction; and seed `20260813`.

The historical B12A outputs retained on disk are used as the verification target: `B12A_bootstrap_all_loadings.csv`, `B12A_bootstrap_global_summary.csv`, `B12A_bootstrap_loading_summary.csv`, `B12A_leave_one_cancer_out_feature_summary.csv`, `B12A_leave_one_cancer_out_loadings.csv`, and `B12A_leave_one_cancer_out_summary.csv`.

Final archive evidence reports approximately mean bootstrap loading congruence 0.99997 and minimum LOCO congruence 0.99993. This is a reconstruction of the final executed logic, not a claimed byte-for-byte recovery of a missing historical script.
