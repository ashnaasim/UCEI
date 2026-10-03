# Analysis scripts

Publication-facing scripts reconstructed from the final supported Stage 3 / UCEI Deep / benchmark analyses are stored under `scripts/publication/`.

## Current publication script map

1. `02_ucei_construction_stability.R` — B12A bootstrap and leave-one-cancer-out construction stability.
2. `03_conventional_cna_benchmark.R` — B2 conventional CNA overlap, burden dependence and redundancy benchmark.
3. `04_published_signature_nested_benchmark.R` — B9/B9C Steele/Drews/cDITHER/AllSOTA39 nested benchmark.
4. `05_training_fold_ucei_reconstruction.R` — B12B complete training-fold UCEI reconstruction.
5. `06_construction_ablation_and_higher_dimensional.R` — B12C/B12C2 construction ablation and higher-dimensional alternatives.
6. `07_censoring_aware_discrimination.R` — B12D censoring-aware/time-dependent discrimination recovery and audit.
7. `08_heldout_cancer_transportability.R` — B10 internal-external held-out-cancer transportability.
8. `09_archsig_state_analysis.R` — UCEI–ArchSig architecture–transcriptome states.
9. `10_latent_architecture_manifestation.R` — latent architecture/manifestation factors.
10. `11_depmap_genomewide_dependency.R` — DepMap genome-wide dependency modeling.
11. `12_cn_aware_transcription.R` — CN-aware transcriptional decomposition.
12. `13_proteomic_pathway_analysis.R` — proteomic integration and MSigDB pathway ORA.
13. `14_singlecell_spatial_sensitivity.R` — single-cell CopyKAT and spatial sensitivity analyses.
14. `15_prism_sensitivity.R` — PRISM pharmacologic sensitivity.
15. `16_wgcna_string_network_context.R` — recovered WGCNA/STRING network context.
16. `17_discordance_hinge_models.R` — architecture–manifestation discordance and hinge-response modeling.
17. `18_external_ranked_crispr_validation.R` — external ranked paired-CRISPR validation.
18. `19_external_reduced_ucei_validation.R` — external reduced-feature UCEI/copy-number entropy proxy validation.
19. `20_build_figure_source_data.R` — recovered final Figure 1–6 source-data harvester.

## Provenance

Reconstruction notes are stored in `scripts/provenance/`. Where an exact historical standalone script survived, the retained implementation is used. Where only the executed code-run record and outputs survived, the repository explicitly labels the publication script as a source-supported reconstruction and uses retained numerical outputs as verification targets.

Failed, exploratory, superseded or manuscript-excluded branches are not silently merged into the publication workflow.
