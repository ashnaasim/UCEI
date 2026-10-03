# B10 held-out-cancer transportability: reconstruction provenance

The B10 analysis evaluated cross-histology transportability by internal-external validation: each of the nine cancers in the common clinical benchmark was withheld in turn, model construction used the remaining cancers, and the omitted cancer was used only for evaluation.

## Preserved final analysis objects

The Stage-3 workspace records the objects `B10`, `B10_age_col`, `B10_burden_col`, `B10_cancers`, `B10_cindex`, `B10_event_col`, `B10_meta`, `B10_pick`, `B10_results`, `B10_rows`, `B10_summary`, `B10_train_scale`, `B10_ucei_col` and `B10_z`.

The final figure-source package also explicitly harvested `B10_cindex`, `B10_meta`, `B10_results` and `B10_summary` for Figure 4.

## Locked final results

The retained B10 result table contains 3,243 tumors and 781 PFI events across nine cancers:

| Cancer | n | Events | Baseline C | UCEI C | Delta C | UCEI HR |
|---|---:|---:|---:|---:|---:|---:|
| ACC | 75 | 39 | 0.5252374 | 0.6876562 | 0.1624188 | 2.565250 |
| BRCA | 776 | 91 | 0.5587691 | 0.5616034 | 0.0028343 | 1.107721 |
| KIRC | 506 | 152 | 0.5495625 | 0.5983089 | 0.0487464 | 1.402346 |
| KIRP | 226 | 39 | 0.5595044 | 0.6648190 | 0.1053146 | 1.864863 |
| LGG | 502 | 190 | 0.6178772 | 0.6611143 | 0.0432371 | 1.265126 |
| LIHC | 254 | 109 | 0.5000969 | 0.5440445 | 0.0439475 | 1.116356 |
| PRAD | 443 | 84 | 0.6128100 | 0.6099216 | -0.0028885 | 1.319906 |
| READ | 98 | 21 | 0.4619003 | 0.5851364 | 0.1232361 | 1.657108 |
| UCEC | 363 | 56 | 0.5684187 | 0.5998327 | 0.0314140 | 1.267160 |

Summary: positive Delta C in 8/9 cancers, UCEI HR >1 in 9/9, mean Delta C 0.06202895, median 0.04394752, range -0.002888451 to 0.1624188.

These results support broad directional transportability and substantial between-cancer heterogeneity; they do not support a claim of universal improvement.

## Reconstruction boundary

The complete historical B10 function bodies were not preserved in the searchable transcript. However, the surviving object names establish that training-derived scaling and held-out C-index evaluation were part of the implementation. The final manuscript states that model construction was performed using the remaining cancer types and the held-out cancer was used only for evaluation.

The publication script therefore reconstructs the documented analysis using training-derived scaling, a burden + age baseline model, a burden + age + UCEI model, prediction in the omitted cancer, and a cancer-specific Cox model for held-out UCEI directionality.

Because the exact historical function bodies are unavailable, the script contains a full row-level numerical audit against the retained B10 result table and fails if the reconstruction does not reproduce the historical values within tolerance. Until that audit is run successfully on the canonical 3,243-sample benchmark table, this should be described as a source-supported reconstruction rather than a byte-for-byte recovery of the historical script.
