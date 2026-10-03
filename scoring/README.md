# UCEI scoring

This directory contains the standalone scoring implementation used to reproduce the retained UCEI values.

## Public entry point

Use `score_ucei.R`:

```r
source("scoring/score_ucei.R")
x <- read.csv("scoring/test_input_20_samples.csv", check.names = FALSE)
x$UCEI <- score_ucei(x)
```

By default, the wrapper reads the retained parameter files from the same `scoring/` directory.

## Required input columns

- `cancer`
- `burden_total`
- `jsd_neutral`
- `ent_state_size`
- `ent_chr_sd`
- `ent_chr_mean`
- `ent_ampclass`
- `ent_transition`

## Scoring sequence

1. Residualize each architecture feature against `burden_total` using the retained cancer-specific intercept and slope.
2. Center and scale the residual using the retained cancer/feature-specific residual mean and SD.
3. Project the six standardized residual features onto the retained PC1 loading vector.
4. Center and scale the resulting PC1 score using the retained score parameters.

The scorer does not extrapolate cancer-specific residualization parameters to unsupported histologies.

## Historical checksum-verified scorer

`score_ucei_v2.R`, `README.txt`, `data_dictionary.csv`, the three `UCEI_v2_*.csv` parameter files, and the associated test/session files are retained under their original names because they are covered by `MD5SUMS.txt`.

These historical labels are provenance identifiers, not the reader-facing method name. Do not rename or edit those checksum-covered files without regenerating the checksum manifest.

## Numerical reproduction test

The retained test is:

- samples: 20
- maximum absolute error: 7.7715611723761e-15
- tolerance: 1e-10
- status: PASS

See `standalone_reproduction_test.csv`.

## Environment

The retained scorer session used R 4.5.2 on Windows x86_64. The scorer itself uses base R utilities; broader analysis-script dependencies are documented in `../docs/SOFTWARE.md`.
