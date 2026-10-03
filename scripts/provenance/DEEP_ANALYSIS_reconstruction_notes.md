# Deep-analysis reconstruction boundaries

The publication-facing scripts in this block were reconstructed from the canonical Stage-3/UCEI Deep code-run documents, final manuscript methods, Supplementary Results, final evidence map and retained output summaries.

## 09 UCEI–ArchSig states
Source-supported elements include cancer-specific median state definition; four-state event counts/rates; burden- and cancer-adjusted Cox comparisons; separate joint-high-versus-all-others models; cancer-specific UCEI–ArchSig Spearman correlations; and the non-significant continuous interaction. Locked headline values are 3,265 tumors / 783 events; event rates 14.2%, 21.1%, 25.8%, 34.3%; four-state HRs 1.51, 2.01, 2.90; binary joint-high HRs 2.13 and 2.09; clinical complete-case HR 1.37.

## 10 Latent architecture–manifestation
The canonical record supports a six-feature architecture PC1 and a seven-module manifestation PC1, explaining 62.55% and 88.24% of variance. The adjusted architecture-to-manifestation coefficient was beta=0.2286, SE=0.0118, P=1.06e-79, n=3,983, full-model R2=0.578 and Spearman rho=0.309. Exact historical column labels varied across runs, so the publication script resolves prepared final columns rather than claiming byte-for-byte recovery.

## 11 DepMap genome-wide dependency
The final method used DepMap Public 26Q1 and 852 aligned models. For each of 17,566 genes, CRISPR gene effect was modeled against architecture–manifestation activity with lineage, target expression, target CN and background/global dependency adjustment. Candidate dependencies required beta_AM<0 and BH FDR<0.10. Locked totals: 638 FDR-supported and 753 additional nominal directional genes. The script expects the harmonized long-format table and does not redistribute raw DepMap data.

## 12 CN-aware transcription
The final 11A code-run record explicitly states that each gene's expression was modeled as a function of matched local gene-level CN, reduced UCEI architecture, global CN burden and lineage. Locked totals are 18,393 evaluable genes, 15,561 positive local-CN effects at FDR<0.10 and zero individual positive residual UCEI effects at FDR<0.10.

## 13 Proteomics/pathway analysis
The proteomic source was CCLE/Nusinow quantitative mass spectrometry, with the Gygi normalized file. The aligned proteomic background was 7,408 genes. The exact four-class directional/significance decision table must accompany source data; the publication script therefore requires the already classified gene table instead of re-inventing class assignments. Locked counts are 68 hidden dependencies, 64 protein-over-RNA decoupling genes, 28 protein-linked dependencies and 5 compensation candidates (165 total). Final pathway ORA used MSigDB 7.5.1, a 7,408-gene background, 5,076 eligible pathways, pathway sizes 10–500, overlap >=2, one-sided hypergeometric tests and global BH correction per query class.

## 14 Single-cell/spatial sensitivity
The exact raw-data workflows are large and were explicitly supplementary. The final evidence is mixed/weak: mapped CopyKAT on eight GSE264573 biopsies was technically complete but biologically weak; GSE278936 processed 32/32 spatial tumors but paired within-sample high/low CNV-like contrasts were not significant. The publication script therefore operates on final post-inference tables and reproduces the reported statistical layers rather than redistributing raw GEO matrices.

## 15 PRISM
The final pharmacologic analysis used PRISM Repurposing Secondary Screen AUC, sign-reversed after within-compound standardization. The aligned resource contained 342 cell lines and 1,502 compounds; compound models required at least 80 complete observations. Models tested architecture–manifestation discordance (and complementary UCEI models) with ArchSig, global dependency and lineage adjustment followed by BH correction. No FDR-supported architecture-associated drug-sensitivity class was identified.

These scripts must be checked against the retained final outputs before the repository is made public. Where exact historical function bodies are not recoverable, the repository labels the code as a source-supported publication reconstruction rather than a byte-for-byte historical script.
