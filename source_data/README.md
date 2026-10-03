# Figure 1–6 source data

This directory contains the final Figure 1–6 source-data package generated from already-computed analysis outputs.

Files:

- `UCEI_FIGURES_1_TO_6_SOURCE_DATA_BUNDLE.rds`
- `UCEI_FIGURES_1_TO_6_SOURCE_MANIFEST.csv`
- `UCEI_FIGURES_1_TO_6_CRITICAL_AUDIT.csv`
- `UCEI_FIGURE_DATA_INDEX.rds`

The historical harvest explicitly performed **no model refitting or re-analysis**. Its critical audit passed 19/19 checks.

The manifest preserves original workstation paths and internal historical object labels as provenance. Those paths are not required to use the public repository and should not be interpreted as portable input locations.

For a portable integrity check, run:

```r
Rscript scripts/publication/20_verify_figure_source_data.R source_data
```

The verifier checks file availability, the historical audit, and the Figure1–Figure6 bundle structure; it does not rerun the analyses.
