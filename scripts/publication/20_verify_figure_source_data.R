############################################################
## 20: VERIFY PUBLISHED FIGURE SOURCE DATA
##
## Portable public verification entry point.
## No analyses are rerun or refitted.
############################################################

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly=TRUE)
source_dir <- if(length(args)>=1L) args[[1L]] else "source_data"

manifest_file <- file.path(source_dir,"UCEI_FIGURES_1_TO_6_SOURCE_MANIFEST.csv")
audit_file <- file.path(source_dir,"UCEI_FIGURES_1_TO_6_CRITICAL_AUDIT.csv")
bundle_file <- file.path(source_dir,"UCEI_FIGURES_1_TO_6_SOURCE_DATA_BUNDLE.rds")
index_file <- file.path(source_dir,"UCEI_FIGURE_DATA_INDEX.rds")

required <- c(manifest_file,audit_file,bundle_file,index_file)
missing <- required[!file.exists(required)]
if(length(missing)) stop("Missing source-data file(s): ",paste(missing,collapse=", "))

MANIFEST <- fread(manifest_file)
AUDIT <- fread(audit_file)
FIGDATA <- readRDS(bundle_file)
INDEX <- readRDS(index_file)

if(!all(c("figure","check","observed","pass") %in% names(AUDIT))) {
    stop("Critical-audit file does not have the expected schema.")
}
if(any(is.na(AUDIT$pass)) || any(!as.logical(AUDIT$pass))) {
    print(AUDIT[is.na(pass) | pass==FALSE])
    stop("One or more historical figure-data audit checks failed.")
}

expected_figures <- paste0("Figure",1:6)
if(!all(expected_figures %in% names(FIGDATA))) {
    stop("Figure-data bundle is missing one or more Figure1-Figure6 entries.")
}

COUNTS <- data.table(
    figure=expected_figures,
    n_items=vapply(FIGDATA[expected_figures],length,integer(1))
)

cat("\nUCEI figure source-data verification\n")
cat("-----------------------------------\n")
cat("Manifest rows:",nrow(MANIFEST),"\n")
cat("Critical checks passed:",sum(as.logical(AUDIT$pass)),"/",nrow(AUDIT),"\n")
print(COUNTS)
cat("\nPASS: source-data bundle and historical audit are internally available.\n")
cat("No model fitting or re-analysis was performed.\n")
