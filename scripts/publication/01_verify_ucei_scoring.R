############################################################
## 01: VERIFY STANDALONE UCEI SCORING
##
## Uses only repository scoring files.
## No analysis model is fitted.
############################################################

args <- commandArgs(trailingOnly=TRUE)
repo_root <- if(length(args)>=1L) args[[1L]] else "."
scoring_dir <- file.path(repo_root,"scoring")

required <- c(
    "score_ucei.R",
    "score_ucei_v2.R",
    "test_input_20_samples.csv",
    "test_expected_output.csv",
    "standalone_reproduction_test.csv",
    "MD5SUMS.txt",
    "UCEI_v2_residualization_parameters.csv",
    "UCEI_v2_PCA_loadings.csv",
    "UCEI_v2_score_parameters.csv"
)

missing <- required[!file.exists(file.path(scoring_dir,required))]
if(length(missing)) {
    stop("Missing scoring file(s): ",paste(missing,collapse=", "))
}

source(file.path(scoring_dir,"score_ucei.R"))

x <- utils::read.csv(
    file.path(scoring_dir,"test_input_20_samples.csv"),
    stringsAsFactors=FALSE,
    check.names=FALSE
)
expected <- utils::read.csv(
    file.path(scoring_dir,"test_expected_output.csv"),
    stringsAsFactors=FALSE,
    check.names=FALSE
)

if(!"sample_id"%in%names(x) || !"sample_id"%in%names(expected)) {
    stop("Test files must contain sample_id.")
}

observed <- score_ucei(x,parameter_dir=scoring_dir)
m <- match(x$sample_id,expected$sample_id)
if(anyNA(m)) stop("Expected-output file does not contain every test sample.")

err <- abs(observed-as.numeric(expected$expected_UCEI_v2_z[m]))
max_err <- max(err,na.rm=TRUE)
tolerance <- 1e-10

cat("UCEI standalone scoring verification\n")
cat("------------------------------------\n")
cat("Samples:",length(observed),"\n")
cat("Maximum absolute error:",format(max_err,digits=16),"\n")
cat("Tolerance:",format(tolerance,scientific=TRUE),"\n")

if(!is.finite(max_err) || max_err>tolerance) {
    stop("Standalone UCEI scoring verification FAILED.")
}

## Check the historical checksum manifest.
manifest <- utils::read.delim(
    file.path(scoring_dir,"MD5SUMS.txt"),
    stringsAsFactors=FALSE,
    check.names=FALSE
)
if(!all(c("file","md5")%in%names(manifest))) {
    stop("MD5SUMS.txt does not have the expected file/md5 columns.")
}
paths <- file.path(scoring_dir,manifest$file)
exists <- file.exists(paths)
actual <- rep(NA_character_,length(paths))
actual[exists] <- as.character(tools::md5sum(paths[exists]))
checksum_ok <- exists & tolower(actual)==tolower(manifest$md5)

cat("Checksum files present:",sum(exists),"/",length(exists),"\n")
cat("Checksum matches:",sum(checksum_ok),"/",length(checksum_ok),"\n")

if(!all(checksum_ok)) {
    print(data.frame(
        file=manifest$file,
        expected=manifest$md5,
        observed=actual,
        pass=checksum_ok
    )[!checksum_ok,,drop=FALSE])
    stop("Historical scoring-file checksum verification FAILED.")
}

cat("PASS: scoring values and historical checksums verified.\n")
