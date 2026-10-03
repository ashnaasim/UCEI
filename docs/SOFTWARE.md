# Software environment

## Retained scorer environment

The checksum-covered standalone scorer session records:

- R 4.5.2
- Windows x86_64
- data.table 1.18.0
- survival 3.8-6
- readxl 1.4.5

The standalone scoring function itself uses base R utilities.

## Publication-script dependencies

The final publication-facing scripts use:

- `data.table`
- `survival`
- `glmnet`

The nested Cox benchmark uses `glmnet::cv.glmnet` with `family="cox"`, `alpha=0.5`, Cox deviance for inner cross-validation, `lambda.min`, `standardize=FALSE`, and unpenalized prespecified baseline covariates.

The exact historical `glmnet` package version was not preserved in the checksum-covered scorer session and is therefore not asserted here. This is a documented provenance limitation rather than an inferred version.

## Reproducibility note

For exact numerical reproduction of the standalone UCEI score, use the checksum-covered files in `scoring/`. For manuscript analyses, consult the corresponding provenance note under `scripts/provenance/` when a publication-facing script is a source-supported reconstruction rather than a byte-for-byte historical file.
