############################################################
## B12D: CENSORING-AWARE AND TIME-DEPENDENT DISCRIMINATION
##
## Publication-facing recovery/audit script.
##
## IMPORTANT PROVENANCE BOUNDARY
## -----------------------------
## The complete historical function body used to estimate
## fold-level Uno C and time-dependent AUC was not preserved
## in the searchable code-run record. The historical B12D
## output CSVs did survive on disk.
##
## Therefore this script DOES NOT silently substitute a new
## IPCW/AUC implementation. It reads the retained historical
## fold/repeat-level metrics and regenerates manuscript
## summaries, paired deltas, support checks, and numerical QC.
##
## Usage:
## Rscript 07_censoring_aware_discrimination.R historical_B12D_dir output_dir
############################################################

suppressPackageStartupMessages({
    library(data.table)
})

args <- commandArgs(trailingOnly=TRUE)
input_dir <- if (length(args)>=1L) args[[1L]] else Sys.getenv("UCEI_B12D_HISTORICAL_DIR",unset="")
output_dir <- if (length(args)>=2L) args[[2L]] else file.path("results","benchmark","B12D_censoring_sensitivity")

if (!nzchar(input_dir) || !dir.exists(input_dir)) {
    stop("Supply the historical B12D result directory as argument 1 or set UCEI_B12D_HISTORICAL_DIR.")
}
dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)

find_one <- function(pattern) {
    z <- list.files(input_dir,pattern=pattern,full.names=TRUE,ignore.case=TRUE)
    if (length(z)!=1L) stop("Expected exactly one file matching /",pattern,"/; found ",length(z),".")
    z
}

fold_file <- find_one("^B12D_fold_censoring_adjusted_metrics\\.csv$")
repeat_file <- find_one("^B12D_repeat_censoring_adjusted_metrics.*\\.csv$")
support_file <- find_one("^B12D_fold_horizon_support\\.csv$")

B12D_FOLD <- fread(fold_file)
B12D_REPEAT <- fread(repeat_file)
B12D_SUPPORT <- fread(support_file)

normalize_metric_table <- function(D) {
    nms <- names(D)
    metric_col <- intersect(c("metric","Metric"),nms)[1L]
    model_col <- intersect(c("model","Model"),nms)[1L]
    value_col <- intersect(c("value","Value","estimate","Estimate","metric_value"),nms)[1L]
    rep_col <- intersect(c("rep_id","repeat","Repeat","rep"),nms)[1L]
    fold_col <- intersect(c("fold","outer_fold","Fold"),nms)[1L]

    if (any(is.na(c(metric_col,model_col,value_col)))) {
        stop("Could not identify metric/model/value columns in: ",paste(nms,collapse=", "))
    }

    out <- data.table(
        metric=as.character(D[[metric_col]]),
        model=as.character(D[[model_col]]),
        value=as.numeric(D[[value_col]])
    )
    if (!is.na(rep_col)) out[,rep_id:=D[[rep_col]]]
    if (!is.na(fold_col)) out[,fold:=D[[fold_col]]]
    out
}

RPT <- normalize_metric_table(B12D_REPEAT)
FLD <- normalize_metric_table(B12D_FOLD)

metric_order <- c("UnoC_full","UnoC_5y","AUC_1y","AUC_3y","AUC_5y")

B12D_SUMMARY <- RPT[
    is.finite(value),
    .(
        mean_value=mean(value),
        sd_value=sd(value),
        min_value=min(value),
        max_value=max(value),
        valid_repeats=.N
    ),
    by=.(metric,model)
]

if ("fold" %in% names(FLD) && "rep_id" %in% names(FLD)) {
    fs <- FLD[
        is.finite(value),
        .(valid_folds=uniqueN(fold)),
        by=.(metric,model,rep_id)
    ][
        ,
        .(min_valid_folds=min(valid_folds)),
        by=.(metric,model)
    ]
    B12D_SUMMARY <- merge(B12D_SUMMARY,fs,by=c("metric","model"),all.x=TRUE)
} else {
    B12D_SUMMARY[,min_valid_folds:=NA_integer_]
}

B12D_SUMMARY[,metric:=factor(metric,levels=metric_order)]
setorder(B12D_SUMMARY,metric,-mean_value)

B12D_PRIMARY <- B12D_SUMMARY[
    metric %in% c("UnoC_5y","AUC_1y","AUC_3y","AUC_5y")
]

if (!"rep_id" %in% names(RPT)) stop("Repeat-level file lacks repeat identifier.")

W <- dcast(RPT,rep_id + metric ~ model,value.var="value")

required_models <- c("baseline","strict_UCEI","AllSOTA39","AllSOTA39_plus_strict_UCEI")
miss_models <- setdiff(required_models,names(W))
if (length(miss_models)) stop("Missing B12D model columns: ",paste(miss_models,collapse=", "))

B12D_DELTA <- rbindlist(list(
    W[,.(rep_id,metric,comparison="combo - AllSOTA39",
          delta=AllSOTA39_plus_strict_UCEI-AllSOTA39)],
    W[,.(rep_id,metric,comparison="combo - strict_UCEI",
          delta=AllSOTA39_plus_strict_UCEI-strict_UCEI)],
    W[,.(rep_id,metric,comparison="strict_UCEI - AllSOTA39",
          delta=strict_UCEI-AllSOTA39)],
    W[,.(rep_id,metric,comparison="strict_UCEI - baseline",
          delta=strict_UCEI-baseline)]
))

B12D_DELTA_SUMMARY <- B12D_DELTA[
    is.finite(delta),
    .(
        mean_delta=mean(delta),
        sd_delta=sd(delta),
        min_delta=min(delta),
        max_delta=max(delta),
        positive_repeats=sum(delta>0),
        negative_repeats=sum(delta<0),
        total_repeats=.N
    ),
    by=.(metric,comparison)
]

support_names <- names(B12D_SUPPORT)
first_present <- function(x) {
    y <- intersect(x,support_names)
    if (length(y)) y[1L] else NA_character_
}

fold_col <- first_present(c("fold","outer_fold","Fold"))
e1 <- first_present(c("events_1y","n_events_1y"))
c1 <- first_present(c("controls_1y","n_controls_1y"))
e3 <- first_present(c("events_3y","n_events_3y"))
c3 <- first_present(c("controls_3y","n_controls_3y"))
e5 <- first_present(c("events_5y","n_events_5y"))
c5 <- first_present(c("controls_5y","n_controls_5y"))

if (any(is.na(c(e1,c1,e3,c3,e5,c5)))) {
    stop("Could not identify 1/3/5-y event/control support columns.")
}

B12D_SUPPORT_SUMMARY <- data.table(
    folds=if (!is.na(fold_col)) uniqueN(B12D_SUPPORT[[fold_col]]) else nrow(B12D_SUPPORT),
    min_events_1y=min(B12D_SUPPORT[[e1]],na.rm=TRUE),
    min_controls_1y=min(B12D_SUPPORT[[c1]],na.rm=TRUE),
    min_events_3y=min(B12D_SUPPORT[[e3]],na.rm=TRUE),
    min_controls_3y=min(B12D_SUPPORT[[c3]],na.rm=TRUE),
    min_events_5y=min(B12D_SUPPORT[[e5]],na.rm=TRUE),
    min_controls_5y=min(B12D_SUPPORT[[c5]],na.rm=TRUE)
)

LOCKED <- data.table(
    metric=c(rep("UnoC_5y",4),rep("AUC_1y",4),rep("AUC_3y",4),rep("AUC_5y",4)),
    model=c(
        "strict_UCEI","AllSOTA39_plus_strict_UCEI","AllSOTA39","baseline",
        "AllSOTA39_plus_strict_UCEI","strict_UCEI","AllSOTA39","baseline",
        "AllSOTA39_plus_strict_UCEI","strict_UCEI","AllSOTA39","baseline",
        "AllSOTA39_plus_strict_UCEI","strict_UCEI","AllSOTA39","baseline"
    ),
    expected_mean=c(
        0.6941794,0.6933424,0.6870363,0.6733994,
        0.7575508,0.7526530,0.7517836,0.7377908,
        0.7321776,0.7320194,0.7314566,0.7256819,
        0.7494020,0.7484906,0.7431265,0.7334740
    ),
    expected_sd=c(
        0.005032872,0.004359947,0.002541810,0.004470388,
        0.003018891,0.002466247,0.003153611,0.003192277,
        0.002941699,0.001736475,0.004061018,0.002015790,
        0.003815747,0.003003693,0.003501308,0.003580478
    )
)

AUDIT <- merge(
    LOCKED,
    B12D_PRIMARY[,.(metric=as.character(metric),model,observed_mean=mean_value,observed_sd=sd_value)],
    by=c("metric","model"),
    all.x=TRUE
)

AUDIT[,abs_mean_diff:=abs(observed_mean-expected_mean)]
AUDIT[,abs_sd_diff:=abs(observed_sd-expected_sd)]
AUDIT[,pass:=is.finite(abs_mean_diff)&is.finite(abs_sd_diff)&abs_mean_diff<5e-7&abs_sd_diff<5e-7]

fwrite(FLD,file.path(output_dir,"B12D_fold_censoring_adjusted_metrics.csv"))
fwrite(B12D_SUPPORT,file.path(output_dir,"B12D_fold_horizon_support.csv"))
fwrite(B12D_SUPPORT_SUMMARY,file.path(output_dir,"B12D_horizon_support_summary.csv"))
fwrite(B12D_SUMMARY,file.path(output_dir,"B12D_metric_summary.csv"))
fwrite(B12D_DELTA_SUMMARY,file.path(output_dir,"B12D_paired_delta_summary.csv"))
fwrite(B12D_DELTA,file.path(output_dir,"B12D_paired_deltas_by_repeat.csv"))
fwrite(B12D_PRIMARY,file.path(output_dir,"B12D_PRIMARY_censoring_adjusted_results.csv"))
fwrite(RPT,file.path(output_dir,"B12D_repeat_censoring_adjusted_metrics.csv"))
fwrite(AUDIT,file.path(output_dir,"B12D_locked_numerical_audit.csv"))

cat("\nB12D historical-result reconstruction complete.\n")
cat("Locked numerical checks passed: ",sum(AUDIT$pass,na.rm=TRUE),"/",nrow(AUDIT),"\n",sep="")
print(B12D_PRIMARY)
print(B12D_SUPPORT_SUMMARY)

if (!all(AUDIT$pass)) {
    stop("B12D locked numerical audit failed; reconcile historical inputs before public release.")
}
