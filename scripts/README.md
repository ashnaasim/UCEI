# Analysis scripts

Publication-facing analysis scripts are stored under `scripts/publication/`. Historical workstation-specific recovery code and provenance notes are stored under `scripts/provenance/`.

## Publication script map

1. `02_ucei_construction_stability.R` — B12A bootstrap and leave-one-cancer-out construction stability.
2. `03_conventional_cna_benchmark.R` — B2 conventional CNA overlap, burden dependence and redundancy benchmark.
3. `04_published_signature_nested_benchmark.R` — B9/B9C Steele/Drews/cDITHER/39-feature nested benchmark.
4. `05_training_fold_ucei_reconstruction.R` — B12B complete training-fold UCEI reconstruction.
5. `06_construction_ablation_and_higher_dimensional.R` — B12C/B12C2 construction ablation and higher-dimensional alternatives.
6. `07_censoring_aware_discrimination.R` — B12D censoring-aware and time-dependent discrimination analysis.
7. `08_heldout_cancer_transportability.R` — B10 internal-external held-out-cancer transportability.
8. `09_archsig_state_analysis.R` — UCEI–ArchSig architecture–transcriptome states.
9. `10_latent_architecture_manifestation.R` — latent architecture and manifestation factors.
10. `11_depmap_genomewide_dependency.R` — DepMap genome-wide dependency modeling.
11. `12_cn_aware_transcription.R` — CN-aware transcriptional decomposition.
12. `13_proteomic_pathway_analysis.R` — proteomic integration and pathway over-representation analysis.
13. `14_singlecell_spatial_sensitivity.R` — single-cell and spatial sensitivity analyses.
14. `15_prism_sensitivity.R` — PRISM pharmacologic sensitivity analysis.
15. `16_wgcna_string_network_context.R` — WGCNA/STRING network context.
16. `17_discordance_hinge_models.R` — architecture–manifestation discordance and linear-versus-hinge modeling.
17. `18_external_ranked_crispr_validation.R` — external ranked CRISPR validation.
18. `19_external_reduced_ucei_validation.R` — independent reduced-feature UCEI proxy analyses.
19. `20_verify_figure_source_data.R` — portable verification of the published Figure 1–6 source-data bundle and historical 19-check audit.

## Provenance

Where an exact historical standalone script survived, the retained implementation is preserved. Where only executed code-run records and outputs survived, the public script is explicitly documented as a source-supported reconstruction and retained numerical results are used as verification targets.

The executed workstation-specific Figure 1–6 harvester is retained separately as `scripts/provenance/20_build_figure_source_data_historical.R`; it is not the portable public entry point.

Failed, exploratory, superseded or manuscript-excluded branches are not silently merged into the publication workflow.
