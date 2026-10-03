# Public datasets and stable identifiers

This file records public resources used in the UCEI manuscript. Raw third-party data should be obtained from the originating repositories and are not redistributed here.

## Primary pan-cancer data

- **The Cancer Genome Atlas (TCGA) / NCI Genomic Data Commons (GDC)** — segment-level copy-number and matched transcriptomic resources used for pan-cancer derivation and downstream analyses.
- **TCGA Pan-Cancer Clinical Data Resource (TCGA-CDR)** — harmonized progression-free interval and clinical outcome information used in survival analyses.

## Functional genomics

- **DepMap Public 26Q1** — CRISPR gene-effect, gene-expression, WGS-derived gene-level copy-number, and model metadata used in functional-genomic analyses.
- **PRISM Repurposing 19Q4 Dataset, version 4** — pharmacologic sensitivity resource; Figshare DOI: **10.6084/m9.figshare.9393293**. The UCEI workflow used the secondary-screen dose-response, treatment-information, and cell-line-information files.

## Single-cell and spatial transcriptomics

- **GEO GSE264573** — single-cell prostate cancer dataset used for RNA-inferred CNA/CopyKAT sensitivity analysis. BioProject: **PRJNA1103171**.
- **GEO GSE278936** — prostate cancer spatial transcriptomic dataset used for spatial CNA-proxy sensitivity analysis. BioProject: **PRJNA1169790**.

## Independent external prostate cancer cohorts

External reduced-feature UCEI analyses used cBioPortal-format study packs with the following study identifiers:

- `prad_mskcc`
- `prad_p1000`
- `prad_su2c_2019`
- `prostate_msk_2024`

These analyses are reduced-feature/proxy support and are not presented as exact six-feature UCEI reconstruction.

## Proteomics

Quantitative CCLE proteomic data were taken from the Nusinow et al. mass-spectrometry resource, **“Quantitative Proteomics of the Cancer Cell Line Encyclopedia”** (Cell, 2020; PMID **31978347**; DOI **10.1016/j.cell.2019.12.023**). The retained analysis used the normalized quantitative proteomics file `protein_quant_current_normalized.csv.gz`. The publication workflow aligned this resource to the functional-genomic analysis and used 7,408 genes as the measurable proteomic background.

## Redistribution

Repository files contain derived parameters, processed analysis outputs, and figure source data. Users should retrieve raw third-party datasets directly from the original repositories and comply with their respective terms of use.
