# UCEI scoring package

The final standalone scorer was originally validated in the benchmark/reproducibility workflow.

## Required columns

- `cancer`
- `burden_total`
- `jsd_neutral`
- `ent_state_size`
- `ent_chr_sd`
- `ent_chr_mean`
- `ent_ampclass`
- `ent_transition`

## Scoring sequence

1. Residualize each raw architecture feature against `burden_total` using cancer-specific retained intercept and slope parameters.
2. Center and scale residuals using retained cancer/feature-specific residual means and SDs.
3. Combine the six standardized residual features using the retained PC1 loading vector.
4. Transform the resulting PC1 score using the retained score center and scale.

The retained scorer does not silently extrapolate cancer-specific residualization parameters to unsupported histologies.

## Important provenance note

The canonical recovery transcript contains only the first portion of `score_ucei_v2.R`, while the original workstation archive records the complete tested file and its checksum. Therefore the runtime scorer itself must be copied from the verified local file before public release rather than completed by inference.

The parameter tables recovered independently from the locked manuscript tables are being staged under `parameters/`.
