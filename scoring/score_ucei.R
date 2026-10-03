# Reader-facing UCEI scoring wrapper.
#
# The checksum-verified historical implementation is retained unchanged in
# score_ucei_v2.R. This file exposes reader-facing function names without
# modifying the historical scorer or its checksum-covered parameter files.

.ucei_scoring_dir <- local({
    ofile <- tryCatch(sys.frame(1)$ofile, error=function(e) NULL)
    if (!is.null(ofile) && nzchar(ofile)) {
        dirname(normalizePath(ofile, winslash="/", mustWork=FALSE))
    } else {
        getwd()
    }
})

source(file.path(.ucei_scoring_dir, "score_ucei_v2.R"), local=environment())

score_ucei <- function(input, parameter_dir=.ucei_scoring_dir) {
    score_ucei_v2(input=input, parameter_dir=parameter_dir)
}

score_ucei_file <- function(input_csv, output_csv, parameter_dir=.ucei_scoring_dir) {
    d <- utils::read.csv(input_csv, stringsAsFactors=FALSE, check.names=FALSE)
    d$UCEI <- score_ucei(d, parameter_dir=parameter_dir)
    utils::write.csv(d, output_csv, row.names=FALSE)
    invisible(d)
}
