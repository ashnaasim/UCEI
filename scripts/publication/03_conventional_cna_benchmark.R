############################################################
## B2: CONVENTIONAL CNA BENCHMARK AND NON-REDUNDANCY
##
## Publication-facing reconstruction of the final B2 analysis.
##
## This script operates on prepared patient-level tables rather
## than redistributing raw TCGA data.
##
## Required reference-table columns:
##   cancer, burden_total, UCEI
## plus the conventional/architecture metrics used below.
##
## Optional clinical benchmark columns:
##   patient_id, cancer, pfi_time, pfi_event, burden_total,
##   UCEI and age (when available).
##
## Usage:
## Rscript 03_conventional_cna_benchmark.R \
##   reference_table.csv clinical_benchmark.csv output_dir
############################################################

suppressPackageStartupMessages({
    library(data.table)
    library(survival)
})

args <- commandArgs(trailingOnly=TRUE)
reference_file <- if (length(args)>=1L) args[[1L]] else Sys.getenv("UCEI_B2_REFERENCE",unset="")
clinical_file  <- if (length(args)>=2L) args[[2L]] else Sys.getenv("UCEI_B2_CLINICAL",unset="")
output_dir     <- if (length(args)>=3L) args[[3L]] else file.path("results","benchmark","B2_conventional_cna")

if (!nzchar(reference_file) || !file.exists(reference_file)) {
    stop("Supply the prepared UCEI reference table as argument 1 or set UCEI_B2_REFERENCE.")
}
dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)

read_table <- function(f) {
    if (grepl("\\.rds$",f,ignore.case=TRUE)) as.data.table(readRDS(f)) else fread(f)
}
pick_ucei <- function(D) {
    z <- intersect(c("UCEI_v2_z","UCEI_z","UCEI","ucei"),names(D))
    if (!length(z)) stop("Could not identify UCEI score column.")
    z[1L]
}
zsafe <- function(x) {
    x <- as.numeric(x)
    s <- sd(x,na.rm=TRUE)
    m <- mean(x,na.rm=TRUE)
    if (!is.finite(s)||s<=0) return(rep(NA_real_,length(x)))
    (x-m)/s
}
cindex <- function(time,event,risk) {
    as.numeric(survival::concordance(Surv(time,event)~risk,reverse=TRUE)$concordance)
}

D <- read_table(reference_file)
ucei_col <- pick_ucei(D)
D[, UCEI_for_B2:=as.numeric(get(ucei_col))]

############################################################
## 1. Pairwise UCEI correlations
############################################################

B2_comparators <- c(
    "ent_ampclass",
    "ent_chr_mean",
    "ent_chr_sd",
    "ent_state_size",
    "ent_state",
    "n_segments",
    "sd_cn",
    "frac_del",
    "ent_transition",
    "burden_total",
    "jsd_neutral",
    "frac_amp",
    "mean_cn"
)

available <- B2_comparators[B2_comparators %in% names(D)]

B2_COR <- rbindlist(lapply(available,function(v) {
    ok <- is.finite(D$UCEI_for_B2) & is.finite(as.numeric(D[[v]]))
    data.table(
        metric=v,
        n=sum(ok),
        spearman_rho=cor(D$UCEI_for_B2[ok],as.numeric(D[[v]][ok]),method="spearman"),
        pearson_r=cor(D$UCEI_for_B2[ok],as.numeric(D[[v]][ok]),method="pearson")
    )
}))
setorder(B2_COR,-abs(spearman_rho))

############################################################
## 2. Dependence of each metric on total CNA burden
############################################################

B2_BURDEN <- rbindlist(lapply(available,function(v) {
    if (v=="burden_total") {
        return(data.table(metric=v,n=sum(is.finite(D$burden_total)),R2=1,pearson_r=1,spearman_rho=1))
    }
    x <- as.numeric(D$burden_total)
    y <- as.numeric(D[[v]])
    ok <- is.finite(x)&is.finite(y)
    fit <- lm(y[ok]~x[ok])
    data.table(
        metric=v,
        n=sum(ok),
        R2=summary(fit)$r.squared,
        pearson_r=cor(x[ok],y[ok],method="pearson"),
        spearman_rho=cor(x[ok],y[ok],method="spearman")
    )
}))

############################################################
## 3. Multivariable redundancy
## Final conventional six-metric panel:
## burden_total + n_segments + sd_cn + mean_cn + frac_amp +
## frac_del
############################################################

conv6 <- c("burden_total","n_segments","sd_cn","mean_cn","frac_amp","frac_del")
miss6 <- setdiff(conv6,names(D))
if (length(miss6)) stop("Missing conventional CNA metrics: ",paste(miss6,collapse=", "))

R <- D[,c("UCEI_for_B2",conv6),with=FALSE]
R <- R[complete.cases(R)]
fit_red <- lm(
    UCEI_for_B2 ~ burden_total + n_segments + sd_cn + mean_cn + frac_amp + frac_del,
    data=R
)
sm_red <- summary(fit_red)

B2_REDUNDANCY <- data.table(
    n=nrow(R),
    R2=sm_red$r.squared,
    adjusted_R2=sm_red$adj.r.squared,
    residual_fraction=1-sm_red$r.squared
)

############################################################
## 4. Fixed-burden architecture spread
############################################################

F <- D[
    is.finite(UCEI_for_B2) & is.finite(as.numeric(burden_total)),
    .(cancer,burden_total=as.numeric(burden_total),UCEI=UCEI_for_B2)
]

F[, burden_decile:=cut(
    burden_total,
    breaks=quantile(burden_total,probs=seq(0,1,0.1),na.rm=TRUE),
    include.lowest=TRUE,
    duplicates="drop"
)]

B2_FIXED_BURDEN <- F[
    ,
    .(
        n=.N,
        burden_median=median(burden_total),
        UCEI_mean=mean(UCEI),
        UCEI_sd=sd(UCEI),
        UCEI_IQR=IQR(UCEI),
        UCEI_min=min(UCEI),
        UCEI_max=max(UCEI)
    ),
    by=burden_decile
]

############################################################
## 5. Incremental PFI information beyond CNA burden
############################################################

B2_INCREMENTAL <- NULL
if (nzchar(clinical_file) && file.exists(clinical_file)) {
    C <- read_table(clinical_file)
    uc <- pick_ucei(C)
    required <- c("cancer","pfi_time","pfi_event","burden_total",uc)
    miss <- setdiff(required,names(C))
    if (length(miss)) stop("Clinical benchmark missing: ",paste(miss,collapse=", "))

    C <- C[complete.cases(C[,..required])]
    C[, cancer:=as.character(cancer)]
    C[, pfi_time:=as.numeric(pfi_time)]
    C[, pfi_event:=as.integer(pfi_event)]
    C[, burden_z:=zsafe(burden_total)]
    C[, ucei_z:=zsafe(get(uc))]

    use_age <- "age" %in% names(C) && sum(is.finite(as.numeric(C$age)))==nrow(C)
    if (use_age) C[, age_z:=zsafe(age)]

    base_formula <- if (use_age) {
        Surv(pfi_time,pfi_event) ~ burden_z + age_z + strata(cancer)
    } else {
        Surv(pfi_time,pfi_event) ~ burden_z + strata(cancer)
    }
    full_formula <- if (use_age) {
        Surv(pfi_time,pfi_event) ~ burden_z + age_z + ucei_z + strata(cancer)
    } else {
        Surv(pfi_time,pfi_event) ~ burden_z + ucei_z + strata(cancer)
    }

    m0 <- coxph(base_formula,data=C,ties="efron",x=TRUE)
    m1 <- coxph(full_formula,data=C,ties="efron",x=TRUE)

    lr <- 2*(as.numeric(logLik(m1))-as.numeric(logLik(m0)))
    ddf <- length(coef(m1))-length(coef(m0))
    p <- pchisq(lr,df=ddf,lower.tail=FALSE)

    C0 <- as.numeric(summary(m0)$concordance[1])
    C1 <- as.numeric(summary(m1)$concordance[1])
    b <- coef(m1)["ucei_z"]
    se <- sqrt(vcov(m1)["ucei_z","ucei_z"])

    B2_INCREMENTAL <- data.table(
        n=m1$n,
        events=m1$nevent,
        beta_UCEI=as.numeric(b),
        HR_UCEI=exp(as.numeric(b)),
        CI_low=exp(as.numeric(b)-1.96*se),
        CI_high=exp(as.numeric(b)+1.96*se),
        p_UCEI=summary(m1)$coefficients["ucei_z","Pr(>|z|)"],
        LRT_chisq=lr,
        LRT_df=ddf,
        LRT_p=p,
        base_C=C0,
        full_C=C1,
        delta_C=C1-C0
    )
}

############################################################
## Export
############################################################

fwrite(B2_COR,file.path(output_dir,"B2_UCEI_vs_comparator_correlations.csv"))
fwrite(B2_BURDEN,file.path(output_dir,"B2_metric_dependence_on_burden.csv"))
fwrite(B2_REDUNDANCY,file.path(output_dir,"B2_UCEI_explained_by_conventional_CN.csv"))
fwrite(B2_FIXED_BURDEN,file.path(output_dir,"B2_fixed_burden_architecture_spread.csv"))
if (!is.null(B2_INCREMENTAL)) fwrite(B2_INCREMENTAL,file.path(output_dir,"B2_incremental_PFI_beyond_burden.csv"))

############################################################
## Locked headline audit
############################################################

AUDIT <- data.table(
    quantity=c(
        "R2_conventional_six",
        "rho_ent_ampclass",
        "rho_ent_chr_mean",
        "rho_ent_chr_sd",
        "rho_ent_state_size",
        "rho_burden_total"
    ),
    expected=c(
        0.6739,
        0.7766,
        0.6670,
        0.6543,
        0.5849,
        0.1459
    )
)

lookup_rho <- function(v) B2_COR[metric==v,spearman_rho][1L]
AUDIT[, observed:=c(
    B2_REDUNDANCY$R2,
    lookup_rho("ent_ampclass"),
    lookup_rho("ent_chr_mean"),
    lookup_rho("ent_chr_sd"),
    lookup_rho("ent_state_size"),
    lookup_rho("burden_total")
)]
AUDIT[, abs_diff:=abs(observed-expected)]
AUDIT[, pass:=is.finite(abs_diff)&abs_diff<5e-4]

fwrite(AUDIT,file.path(output_dir,"B2_locked_headline_audit.csv"))

cat("\nB2 complete.\n")
print(B2_REDUNDANCY)
print(B2_COR)
print(AUDIT)
