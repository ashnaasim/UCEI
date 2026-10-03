# Reproducibility

The final repository should contain:

- R session information for the retained scoring implementation and major analysis pipeline;
- checksums for scorer/package files;
- the 20-sample standalone scorer test input and expected output;
- the standalone reproduction result;
- construction-stability outputs (200 cancer-stratified bootstraps and leave-one-cancer-out reconstruction);
- figure-source-data audit outputs.

The canonical project records report numerical score reconstruction to approximately machine precision and highly stable construction under bootstrap and cancer omission. These checks establish computational reproducibility/construction stability, not prospective clinical validation.
