# Unified Copy-number Entropy Index (UCEI)

Reproducibility repository for the manuscript **“Unified Copy-number Entropy Index for Pan-Cancer Analysis of Copy-number Architecture, Transcriptional Manifestation and Functional Dependencies.”**

UCEI is an outcome-independent, burden-adjusted one-dimensional representation of higher-order copy-number (CN) architecture. The manuscript derivation used six cancer-specific burden-residualized architecture features in 10,495 tumors across 32 cancer types; PC1 explained 61.7377% of residual feature variance.

## Repository status

This repository is currently a **private pre-submission staging archive**. Files are being recovered and audited against the canonical Stage 3 / UCEI Deep / benchmark records before public release. Only source-supported, final-analysis material will be included. Historical exploratory or superseded branches are not part of the publication pipeline.

## Planned contents

- `scoring/` — UCEI scoring implementation, parameter files, test input/output, and scorer documentation.
- `parameters/` — derivation/reference parameters used to reproduce the retained UCEI score.
- `scripts/` — final analysis scripts corresponding to the submitted manuscript.
- `source_data/` — source data and manifests underlying manuscript Figures 1–6.
- `results/` — processed non-identifying result tables required for reproducibility.
- `reproducibility/` — session information, checksums, reconstruction tests, and audit notes.
- `docs/` — dataset identifiers, provenance, and repository documentation.

## Core scoring inputs

The standalone scorer requires:

`cancer`, `burden_total`, `jsd_neutral`, `ent_state_size`, `ent_chr_sd`, `ent_chr_mean`, `ent_ampclass`, and `ent_transition`.

The scoring procedure is: cancer-specific residualization against total CNA burden; standardization using retained residual means/SDs; projection onto the retained six-feature PC1 loading vector; and final score scaling.

## Public datasets

Raw third-party datasets are **not redistributed** here. Analyses used public resources including TCGA/GDC, TCGA-CDR, DepMap Public 26Q1, PRISM Repurposing 19Q4 v4, GEO GSE264573, GEO GSE278936, and external cBioPortal prostate cancer studies. Stable identifiers are documented in `docs/DATASETS.md`.

## Reproducibility policy

The repository will preserve the exact scoring parameters and final analysis logic used for the manuscript. Internal historical filenames containing labels such as `v2` or `frozen` may be retained where necessary for byte-level provenance, but the reader-facing method is referred to simply as **UCEI**.

## License

Code in this repository is released under the MIT License. Third-party datasets remain subject to the terms of their originating repositories.
