# B2 reconstruction provenance

The B2 conventional-CNA benchmark is reconstructed from the final manuscript methods, the retained root-level B2 output tables, and the audited evidence map.

Source-supported final design:

- pairwise Spearman comparison of UCEI with conventional CNA and individual architecture descriptors;
- univariable dependence of each metric on total CNA burden;
- multivariable redundancy model using burden_total, n_segments, sd_cn, mean_cn, frac_amp and frac_del;
- nested Cox comparison of burden-based PFI modeling with and without UCEI;
- cancer-stratified pooled Cox models, with age included when available in the prespecified benchmark.

Locked headline values include R2=0.6739 for the six conventional CNA metrics jointly explaining UCEI variance, and Spearman correlations of approximately 0.7766 with amplification-class entropy, 0.6670 with chromosome-mean architecture, 0.6543 with chromosome-SD architecture, 0.5849 with state-size entropy and 0.1459 with total CNA burden.

The original monolithic B2 script was not recovered as a standalone R file. The historical B2 directory does retain its exported tables, including 03_UCEI_v2_vs_comparator_correlations.csv, 04_UCEI_explained_by_conventional_CN.csv, 05_metric_dependence_on_copy_number_burden.csv, fixed-burden summaries and PFI benchmark tables. This publication script therefore reconstructs the documented analysis from a prepared patient-level table and contains a numerical headline audit.
