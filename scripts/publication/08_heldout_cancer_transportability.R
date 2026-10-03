############################################################
## B10: HELD-OUT-CANCER TRANSPORTABILITY
##
## Publication-facing reconstruction of the final Stage-3
## internal-external validation by cancer type.
##
## Design supported by the preserved manuscript/analysis record:
##   - common clinical benchmark: 3,243 tumors / 781 PFI events
##   - nine cancers
##   - one cancer withheld at a time
##   - model coefficients learned from the remaining cancers
##   - held-out cancer used only for evaluation
##   - baseline: CNA burden + age
##   - UCEI model: CNA burden + age + UCEI
##   - predictor scaling learned from the training cancers
##   - held-out C-index used for discrimination
##   - within-held-out-cancer Cox model used to summarize
##     UCEI directionality (HR, CI, P)
##
## IMPORTANT
## ---------
## The historical B10 object names survive (B10_train_scale,
## B10_cindex, B10_results, B10_summary), but the complete
## original standalone function bodies were not preserved.
## Therefore this reconstruction contains a locked row-level
## numerical audit and stops if the documented reconstruction
## does not reproduce the retained historical B10 results.
##
## Usage:
## Rscript 08_heldout_cancer_transportability.R \
##   clinical_benchmark.csv output_dir
############################################################

suppressPackageStartupMessages({
    library(data.table)
    library(survival)
})

args <- commandArgs(trailingOnly=TRUE)

input_file <- if (length(args)>=1L) args[[1L]] else
    Sys.getenv("UCEI_B10_INPUT", unset="")

output_dir <- if (length(args)>=2L) args[[2L]] else
    file.path("results","benchmark","B10_heldout_cancer_transportability")

if (!nzchar(input_file) || !file.exists(input_file)) {
    stop("Supply the nine-cancer clinical benchmark as argument 1 or set UCEI_B10_INPUT.")
}

dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)

read_table <- function(f) {
    if (grepl("\\.rds$",f,ignore.case=TRUE)) {
        as.data.table(readRDS(f))
    } else {
        fread(f)
    }
}

B10_pick <- function(D,candidates,label) {
    hit <- intersect(candidates,names(D))
    if (!length(hit)) stop("Could not identify ",label," column.")
    hit[1L]
}

B10_cindex <- function(time,event,risk) {
    as.numeric(
        survival::concordance(
            Surv(time,event) ~ risk,
            reverse=TRUE
        )$concordance
    )
}

B10_train_scale <- function(train,test,column) {
    xtr <- as.numeric(train[[column]])
    xte <- as.numeric(test[[column]])
    mu <- mean(xtr,na.rm=TRUE)
    ss <- sd(xtr,na.rm=TRUE)
    if (!is.finite(ss) || ss<=0) stop("Invalid training SD for ",column)
    list(
        train=(xtr-mu)/ss,
        test=(xte-mu)/ss,
        center=mu,
        scale=ss
    )
}

B10 <- read_table(input_file)

B10_ucei_col <- B10_pick(
    B10,
    c("UCEI_v2_z","UCEI_z","UCEI","ucei"),
    "UCEI"
)
B10_burden_col <- B10_pick(
    B10,
    c("burden_total","burden","CNA_burden"),
    "CNA burden"
)
B10_age_col <- B10_pick(
    B10,
    c("age","age_at_diagnosis"),
    "age"
)
B10_event_col <- B10_pick(
    B10,
    c("pfi_event","PFI_event"),
    "PFI event"
)
B10_time_col <- B10_pick(
    B10,
    c("pfi_time","PFI_time"),
    "PFI time"
)

required <- c(
    "cancer",
    B10_time_col,
    B10_event_col,
    B10_burden_col,
    B10_age_col,
    B10_ucei_col
)

B10 <- B10[complete.cases(B10[,..required])]
B10[, cancer:=tolower(as.character(cancer))]
B10[, B10_time:=as.numeric(get(B10_time_col))]
B10[, B10_event:=as.integer(get(B10_event_col))]
B10[, B10_burden:=as.numeric(get(B10_burden_col))]
B10[, B10_age:=as.numeric(get(B10_age_col))]
B10[, B10_ucei:=as.numeric(get(B10_ucei_col))]

B10_cancers <- sort(unique(B10$cancer))

if (length(B10_cancers)!=9L) {
    stop("B10 expects nine cancers; found ",length(B10_cancers),": ",
         paste(B10_cancers,collapse=", "))
}

if (nrow(B10)!=3243L || sum(B10$B10_event)!=781L) {
    warning(
        "B10 historical denominator was 3,243 samples / 781 events; current input has ",
        nrow(B10)," / ",sum(B10$B10_event),"."
    )
}

B10_rows <- vector("list",length(B10_cancers))

for (i in seq_along(B10_cancers)) {

    cc <- B10_cancers[[i]]

    tr <- B10[cancer!=cc]
    te <- B10[cancer==cc]

    ## All preprocessing constants are learned from the training
    ## cancer set and applied unchanged to the held-out cancer.
    sc_b <- B10_train_scale(tr,te,"B10_burden")
    sc_a <- B10_train_scale(tr,te,"B10_age")
    sc_u <- B10_train_scale(tr,te,"B10_ucei")

    tr2 <- copy(tr)
    te2 <- copy(te)

    tr2[, burden_z:=sc_b$train]
    te2[, burden_z:=sc_b$test]
    tr2[, age_z:=sc_a$train]
    te2[, age_z:=sc_a$test]
    tr2[, UCEI_z:=sc_u$train]
    te2[, UCEI_z:=sc_u$test]

    ## A held-out cancer is an unseen stratum/level. The
    ## transport model therefore uses only predictors that can
    ## be applied to an unseen histology.
    m0 <- coxph(
        Surv(B10_time,B10_event) ~ burden_z + age_z,
        data=tr2,
        ties="efron",
        x=TRUE
    )

    m1 <- coxph(
        Surv(B10_time,B10_event) ~ burden_z + age_z + UCEI_z,
        data=tr2,
        ties="efron",
        x=TRUE
    )

    lp0 <- as.numeric(
        predict(m0,newdata=te2,type="lp",reference="zero")
    )
    lp1 <- as.numeric(
        predict(m1,newdata=te2,type="lp",reference="zero")
    )

    C0 <- B10_cindex(
        te2$B10_time,
        te2$B10_event,
        lp0
    )

    C1 <- B10_cindex(
        te2$B10_time,
        te2$B10_event,
        lp1
    )

    ## Directionality in the held-out cancer itself.
    ## Predictors remain on the training-derived scale.
    heldout_fit <- coxph(
        Surv(B10_time,B10_event) ~ burden_z + age_z + UCEI_z,
        data=te2,
        ties="efron",
        x=TRUE
    )

    sm <- summary(heldout_fit)
    beta <- unname(coef(heldout_fit)["UCEI_z"])
    se <- sqrt(vcov(heldout_fit)["UCEI_z","UCEI_z"])

    B10_rows[[i]] <- data.table(
        held_out_cancer=cc,
        test_n=nrow(te2),
        test_events=sum(te2$B10_event),
        train_n=nrow(tr2),
        train_events=sum(tr2$B10_event),
        C_baseline=C0,
        C_UCEI=C1,
        delta_C=C1-C0,
        heldout_UCEI_beta=beta,
        heldout_UCEI_SE=se,
        heldout_UCEI_HR=exp(beta),
        CI_low=exp(beta-1.96*se),
        CI_high=exp(beta+1.96*se),
        heldout_p=sm$coefficients["UCEI_z","Pr(>|z|)"]
    )
}

B10_results <- rbindlist(B10_rows)

B10_summary <- B10_results[
    ,
    .(
        n=nrow(B10),
        events=sum(B10$B10_event),
        cancers=.N,
        positive_delta_C=sum(delta_C>0),
        total_cancers=.N,
        mean_delta_C=mean(delta_C),
        median_delta_C=median(delta_C),
        min_delta_C=min(delta_C),
        max_delta_C=max(delta_C),
        positive_HR=sum(heldout_UCEI_HR>1)
    )
]

############################################################
## Locked historical row-level audit
############################################################

LOCKED <- data.table(
    held_out_cancer=c(
        "acc","brca","kirc","kirp","lgg",
        "lihc","prad","read","ucec"
    ),
    test_n=c(75,776,506,226,502,254,443,98,363),
    test_events=c(39,91,152,39,190,109,84,21,56),
    C_baseline=c(
        0.5252374,0.5587691,0.5495625,0.5595044,0.6178772,
        0.5000969,0.6128100,0.4619003,0.5684187
    ),
    C_UCEI=c(
        0.6876562,0.5616034,0.5983089,0.6648190,0.6611143,
        0.5440445,0.6099216,0.5851364,0.5998327
    ),
    delta_C=c(
        0.162418791,0.002834336,0.048746436,0.105314640,
        0.043237141,0.043947521,-0.002888451,0.123236124,
        0.031414011
    ),
    heldout_UCEI_beta=c(
        0.9420559,0.1023045,0.3381466,0.6231875,0.2351715,
        0.1100697,0.2775603,0.5050737,0.2367780
    ),
    heldout_UCEI_SE=c(
        0.21595721,0.10813617,0.07808597,0.14448835,0.06066027,
        0.09836339,0.10894490,0.25880094,0.14083988
    ),
    heldout_UCEI_HR=c(
        2.565250,1.107721,1.402346,1.864863,1.265126,
        1.116356,1.319906,1.657108,1.267160
    ),
    CI_low=c(
        1.6799764,0.8961546,1.2033365,1.4049340,1.1233074,
        0.9206065,1.0661226,0.9978263,0.9614928
    ),
    CI_high=c(
        3.917023,1.369234,1.634268,2.475357,1.424849,
        1.353728,1.634100,2.751988,1.670001
    ),
    heldout_p=c(
        1.287411e-05,3.441123e-01,1.488113e-05,1.610074e-05,
        1.058123e-04,2.631356e-01,1.084317e-02,5.098672e-02,
        9.272686e-02
    )
)

AUDIT <- merge(
    LOCKED,
    B10_results,
    by="held_out_cancer",
    suffixes=c("_expected","_observed"),
    all.x=TRUE
)

for (v in c(
    "C_baseline","C_UCEI","delta_C",
    "heldout_UCEI_beta","heldout_UCEI_SE",
    "heldout_UCEI_HR","CI_low","CI_high","heldout_p"
)) {
    AUDIT[[paste0(v,"_abs_diff")]] <-
        abs(
            AUDIT[[paste0(v,"_observed")]] -
            AUDIT[[paste0(v,"_expected")]]
        )
}

AUDIT[, pass_counts:=
    test_n_expected==test_n_observed &
    test_events_expected==test_events_observed
]

AUDIT[, pass_numeric:=
    C_baseline_abs_diff<5e-6 &
    C_UCEI_abs_diff<5e-6 &
    delta_C_abs_diff<5e-6 &
    heldout_UCEI_beta_abs_diff<5e-5 &
    heldout_UCEI_SE_abs_diff<5e-5 &
    heldout_UCEI_HR_abs_diff<5e-5 &
    CI_low_abs_diff<5e-5 &
    CI_high_abs_diff<5e-5
]

AUDIT[, pass:=pass_counts & pass_numeric]

SUMMARY_LOCKED <- data.table(
    n=3243L,
    events=781L,
    cancers=9L,
    positive_delta_C=8L,
    total_cancers=9L,
    mean_delta_C=0.06202895,
    median_delta_C=0.04394752,
    min_delta_C=-0.002888451,
    max_delta_C=0.1624188,
    positive_HR=9L
)

SUMMARY_AUDIT <- cbind(
    expected=SUMMARY_LOCKED,
    observed=B10_summary
)

############################################################
## Export
############################################################

fwrite(
    B10_results,
    file.path(output_dir,"B10_heldout_cancer_results.csv")
)

fwrite(
    B10_summary,
    file.path(output_dir,"B10_heldout_cancer_summary.csv")
)

fwrite(
    AUDIT,
    file.path(output_dir,"B10_locked_row_level_audit.csv")
)

cat("\nB10 complete.\n")
print(B10_results)
print(B10_summary)
cat("Row-level locked checks passed: ",
    sum(AUDIT$pass,na.rm=TRUE),"/",nrow(AUDIT),"\n",sep="")

if (!all(AUDIT$pass)) {
    stop(
        "The reconstructed B10 implementation does not reproduce every ",
        "locked historical row within tolerance. Do not treat this file ",
        "as exact historical code until the remaining implementation ",
        "detail is reconciled."
    )
}
