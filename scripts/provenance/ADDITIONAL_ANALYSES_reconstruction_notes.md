# Additional analysis-script provenance

## 16 WGCNA/STRING
The final network layer was recovered from exported WGCNA CSV tables plus a local STRING v12 physical-network file. The historical recovery reported 256 candidate STRING associations and 59 ranked candidate genes. WGCNA is treated as co-expression context; STRING as database-supported network context. No transcription-factor activity inference or experimental physical-interaction claim is made.

## 17 Discordance/hinge
The manuscript defines architecture–manifestation discordance from standardized architecture and manifestation scores and compares linear versus one-breakpoint hinge models. The final checkpoint/DNA-damage result was tau about 0.966, Delta AIC 17.44, hinge beta about 0.366 and BH FDR about 1.08e-4. Other major systems were interpreted as graded. The repository script reconstructs this documented model class from prepared system-level data; it does not claim a physical phase transition.

## 18 External ranked CRISPR
The historical 18D code is partially recoverable verbatim and confirms one-sided Wilcoxon testing of lower external rank percentiles, BH correction, whole-panel and module-level testing. Final WGL results: 28/28 panel genes present, median rank percentile 0.31382284 versus background 0.5005553, P=0.004423512, FDR=0.008847024. DG full-panel enrichment was not significant after correction; complementary top-N overlap was strongest for checkpoint/replication/chromosome-inheritance genes.

## 19 External reduced-feature UCEI
The executed 09A workflow is preserved in the Stage-3 deep code-run record. It computed a reduced CNA-entropy proxy from cBioPortal discrete gene-level CNA states:
state entropy + 0.50*altered fraction + 0.25*amp/del balance + 0.25*(high-amplification fraction + deep-deletion fraction), followed by cohort z-scaling. This is explicitly a reduced proxy and not exact six-feature UCEI. The publication script uses the final prostate study identifiers and the same proxy logic on extracted cBioPortal study packs.

## 20 Figure-source-data harvester
20_build_figure_source_data.R was recovered directly from the final executed figure-data harvester in the canonical September code-run record. Its purpose is source-data collection only: no analysis rerun, no refitting and no manual result fabrication. The historical run generated the four files already archived under source_data/ and passed 19/19 critical checks.
