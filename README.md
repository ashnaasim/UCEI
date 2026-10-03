# Unified Copy-number Entropy Index (UCEI)

Reproducibility repository for the manuscript **“Unified Copy-number Entropy Index for Pan-Cancer Analysis of Copy-number Architecture, Transcriptional Manifestation and Functional Dependencies.”**

UCEI is an outcome-independent, burden-adjusted one-dimensional representation of higher-order copy-number (CN) architecture. The derivation used six cancer-specific burden-residualized architecture features in 10,495 tumors across 32 cancer types; PC1 explained 61.7377% of residual feature variance.

## Repository contents

- `scoring/` — standalone UCEI scorer, exact historical scoring files, parameters, and numerical reproduction test.
- `parameters/` — reader-facing loading table plus retained derivation/reference parameters.
- `scripts/publication/` — publication-facing analysis scripts.
- `scripts/provenance/` — reconstruction notes and workstation-specific historical code retained only for provenance.
- `source_data/` — Figure 1–6 source-data bundle, manifest, index, and critical audit.
- `reproducibility/` — reproducibility notes and audit guidance.
- `docs/` — dataset identifiers, software requirements, and repository audit.
- `results/` — documentation for generated result tables; raw third-party datasets are not redistributed.

## Quick-start scoring

The public entry point is:

```r
source("scoring/score_ucei.R")
x <- read.csv("scoring/test_input_20_samples.csv", check.names = FALSE)
x$UCEI <- score_ucei(x)
```

The expected 20-sample test output is in `scoring/test_expected_output.csv`. The retained standalone test reports a maximum absolute error of approximately `7.77e-15`, well below the prespecified `1e-10` tolerance.

The checksum-verified historical scorer is retained unchanged as `scoring/score_ucei_v2.R`. Historical filenames and internal column labels containing `v2` or `frozen` are preserved only where required for provenance and checksum consistency; the reader-facing method is **UCEI**.

## Scoring inputs

The standalone scorer requires:

`cancer`, `burden_total`, `jsd_neutral`, `ent_state_size`, `ent_chr_sd`, `ent_chr_mean`, `ent_ampclass`, and `ent_transition`.

The scoring sequence is cancer-specific residualization against total CNA burden, residual standardization using retained derivation parameters, projection onto the retained six-feature PC1 loading vector, and final score scaling.

## Public datasets

Raw third-party datasets are not redistributed. The analyses used TCGA/GDC and TCGA-CDR, DepMap Public 26Q1, PRISM Repurposing 19Q4 v4, GEO GSE264573, GEO GSE278936, CCLE quantitative proteomics, and independent cBioPortal prostate cancer cohorts. Stable accessions/releases are listed in `docs/DATASETS.md`.

## Reproducibility and provenance

The repository separates three types of material:

1. **Checksum-verified historical files** — retained unchanged where available.
2. **Publication-facing scripts** — portable scripts corresponding to the final manuscript analyses.
3. **Source-supported reconstructions** — used only where exact standalone historical function bodies were not preserved; these are explicitly documented in `scripts/provenance/` and, where possible, include numerical checks against retained outputs.

The Figure 1–6 source-data package was assembled from already-computed outputs; its historical critical audit passed 19/19 checks. The portable public verifier is `scripts/publication/20_verify_figure_source_data.R`.

## Interpretation boundaries

UCEI is a statistical representation of copy-number architecture, not a direct instantaneous measurement of chromosomal-instability rate. Information-theoretic entropy is not thermodynamic entropy. External reduced-feature scores are proxy analyses rather than exact six-feature UCEI reconstruction, and functional-dependency associations are not equivalent to validated therapeutic targets.

## License

Code is released under the MIT License. Third-party datasets remain subject to the terms of their originating repositories.
