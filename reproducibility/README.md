# Reproducibility

The repository contains several independent reproducibility layers.

## Standalone UCEI scoring

`../scoring/` contains the checksum-verified historical scorer, the full parameter files, a 20-sample test input, expected output, session information, and `MD5SUMS.txt`.

Retained standalone test:

- n = 20
- maximum absolute error = 7.7715611723761e-15
- tolerance = 1e-10
- PASS

## Derivation reconstruction

`../parameters/reference/` contains the retained reference outputs. The stored self-reproduction summary reports 500 tested scores with a worst maximum absolute difference of approximately 7.99e-15.

## Construction stability

The publication scripts and source-data package retain:

- 200 cancer-stratified bootstrap reconstructions;
- leave-one-cancer-out reconstruction across 32 cancers;
- feature/construction ablation analyses.

These evaluate computational construction stability, not prospective clinical validation.

## Figure source data

`../source_data/UCEI_FIGURES_1_TO_6_CRITICAL_AUDIT.csv` records 19/19 successful historical source-data checks. Run:

```r
Rscript scripts/publication/20_verify_figure_source_data.R source_data
```

to verify that the public bundle, manifest, index, and audit are present and internally readable. This verification does not refit any analysis.

## Analysis-code provenance

Where exact historical scripts survived, they are preserved or wrapped without modifying checksum-covered files. Where exact standalone function bodies were not recoverable, `../scripts/provenance/` explicitly records the reconstruction boundary and retained numerical targets.
