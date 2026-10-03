# B12D censoring-aware discrimination: reconstruction provenance

The final B12D analysis is preserved as a supportive survival-discrimination sensitivity analysis, not a clinical-utility claim.

## Directly preserved evidence

The Stage-3 record retains B12D_FOLD, B12D_REPEAT, B12D_SUMMARY, B12D_PRIMARY, B12D_DELTA, B12D_DELTA_SUMMARY, B12D_SUPPORT and B12D_SUPPORT_SUMMARY, together with the historical CSV outputs on the analysis workstation.

The manuscript-facing targets are Uno C-index at 5 years plus time-dependent AUC at 1, 3 and 5 years. These analyses test whether the benchmark conclusion depends on Harrell concordance alone.

Locked mean values are: Uno C at 5 y — UCEI 0.6941794, combined 0.6933424, comparator 0.6870363, baseline 0.6733994; AUC 1 y — 0.7526530, 0.7575508, 0.7517836, 0.7377908; AUC 3 y — 0.7320194, 0.7321776, 0.7314566, 0.7256819; AUC 5 y — 0.7484906, 0.7494020, 0.7431265, 0.7334740, respectively.

The final support audit reported 50 folds, with minimum evaluable event/control counts of 32/214 at 1 year, 61/85 at 3 years and 75/33 at 5 years.

## Recovery boundary

The exact historical R function body used to calculate fold-level IPCW/Uno C and time-dependent AUC is not recoverable from the searchable Stage-3 transcript currently available. The record does not establish which specific R package/function implementation produced those fold-level metrics.

For that reason, 07_censoring_aware_discrimination.R deliberately does not substitute an inferred survAUC, timeROC, riskRegression or other implementation and call it historical code. Instead it reads the retained fold/repeat-level B12D metrics, regenerates summary and paired-delta tables, regenerates the horizon-support summary, and verifies all manuscript-facing values against locked numerical targets.

A full estimator-level recomputation script should only be added if the original B12D estimator code or a definitive package/function record is recovered.
