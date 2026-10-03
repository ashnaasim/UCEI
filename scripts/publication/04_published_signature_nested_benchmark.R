############################################################
## B9 / B9C: PUBLISHED CN SIGNATURE COMPARISON AND
## REPEATED NESTED BENCHMARK
##
## B9:
##   direct same-data Steele21 versus UCEI nested Cox models
##
## B9C:
##   repeated nested benchmark of:
##     clinical baseline
##     UCEI
##     cDITHER
##     Steele21
##     Drews17
##     AllSOTA39
##     AllSOTA39 + UCEI
##
## IMPORTANT:
## Exact historical outer-fold assignments/seeds are not
## recoverable from the searchable transcript. To reproduce
## the submitted numbers exactly, provide the retained fold
## assignment table as argument 2.
##
## Usage:
## Rscript 04_published_signature_nested_benchmark.R \
##   locked_benchmark.csv fold_assignments.csv output_dir
############################################################

suppressPackageStartupMessages({
    library(data.table)
    library(survival)
    library(glmnet)
})

args <- commandArgs(trailingOnly=TRUE)
input_file <- if (length(args)>=1L) args[[1L]] else Sys.getenv("UCEI_B9C_INPUT",unset="")
fold_file <- if (length(args)>=2L) args[[2L]] else Sys.getenv("UCEI_B9C_FOLDS",unset="")
output_dir <- if (length(args)>=3L) args[[3L]] else file.path("results","benchmark","B9_B9C_published_comparators")

if (!nzchar(input_file)||!file.exists(input_file)) stop("Supply the locked B9/B9C benchmark table.")
if (!nzchar(fold_file)||!file.exists(fold_file)) {
    stop(
        "Exact B9C fold assignments are required for historical reproduction. ",
        "Supply fold_assignments.csv as argument 2 or set UCEI_B9C_FOLDS."
    )
}
dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)

read_table <- function(f) {
    if (grepl("\\.rds$",f,ignore.case=TRUE)) as.data.table(readRDS(f)) else fread(f)
}
zsafe <- function(x) {
    x <- as.numeric(x)
    s <- sd(x,na.rm=TRUE); m <- mean(x,na.rm=TRUE)
    if (!is.finite(s)||s<=0) rep(NA_real_,length(x)) else (x-m)/s
}
cindex <- function(time,event,risk) {
    as.numeric(concordance(Surv(time,event)~risk,reverse=TRUE)$concordance)
}

D <- read_table(input_file)
FOLDS <- read_table(fold_file)

ucei_col <- intersect(c("UCEI_v2_z","UCEI_z","UCEI","ucei"),names(D))[1L]
if (is.na(ucei_col)) stop("UCEI column not found.")

steele <- paste0("CN",1:21)
drews <- paste0("CX",1:17)
sota <- c(steele,drews,"cDITHER")

required <- c(
    "patient_id","cancer","pfi_time","pfi_event",
    "burden_total","age",ucei_col,sota
)
miss <- setdiff(required,names(D))
if (length(miss)) stop("Benchmark input missing: ",paste(miss,collapse=", "))

D <- D[complete.cases(D[,..required])]
D[, cancer:=as.character(cancer)]
D[, pfi_time:=as.numeric(pfi_time)]
D[, pfi_event:=as.integer(pfi_event)]
D[, burden_total:=as.numeric(burden_total)]
D[, age:=as.numeric(age)]
D[, UCEI_B9C:=as.numeric(get(ucei_col))]
for (nm in sota) set(D,j=nm,value=as.numeric(D[[nm]]))

if (anyDuplicated(D$patient_id)) stop("Duplicate patient_id in benchmark input.")

############################################################
## B9: same-data bidirectional Steele comparison
############################################################

B9 <- copy(D)
B9[, burden_z:=zsafe(burden_total)]
B9[, age_z:=zsafe(age)]
B9[, UCEI_z2:=zsafe(UCEI_B9C)]
for (nm in steele) B9[[paste0(nm,"_z")]] <- zsafe(B9[[nm]])

steele_z <- paste0(steele,"_z")
base_terms <- c("burden_z","age_z","strata(cancer)")

make_formula <- function(extra) {
    as.formula(paste("Surv(pfi_time,pfi_event) ~",paste(c(base_terms,extra),collapse=" + ")))
}

m_base <- coxph(make_formula(character()),data=B9,ties="efron",x=TRUE)
m_ucei <- coxph(make_formula("UCEI_z2"),data=B9,ties="efron",x=TRUE)
m_steele <- coxph(make_formula(steele_z),data=B9,ties="efron",x=TRUE)
m_both <- coxph(make_formula(c("UCEI_z2",steele_z)),data=B9,ties="efron",x=TRUE)

lrt <- function(full,reduced,label) {
    lr <- 2*(as.numeric(logLik(full))-as.numeric(logLik(reduced)))
    df <- length(coef(full))-length(coef(reduced))
    data.table(
        comparison=label,
        chisq=lr,
        df=df,
        p=pchisq(lr,df=df,lower.tail=FALSE),
        reduced_C=as.numeric(summary(reduced)$concordance[1]),
        full_C=as.numeric(summary(full)$concordance[1]),
        delta_C=as.numeric(summary(full)$concordance[1])-as.numeric(summary(reduced)$concordance[1])
    )
}

B9_LRT <- rbindlist(list(
    lrt(m_both,m_steele,"UCEI added to Steele21"),
    lrt(m_both,m_ucei,"Steele21 added to UCEI")
))

############################################################
## B9C: exact historical outer folds supplied externally
############################################################

fold_required <- c("patient_id","repeat","fold")
fm <- setdiff(fold_required,names(FOLDS))
if (length(fm)) stop("Fold file must contain patient_id, repeat, fold.")

FOLDS[, patient_id:=as.character(patient_id)]
D[, patient_id:=as.character(patient_id)]

if (anyDuplicated(FOLDS[,.(patient_id,repeat)])) {
    stop("Fold assignment table has duplicate patient_id/repeat pairs.")
}

repeats <- sort(unique(FOLDS$repeat))
if (length(repeats)!=10L) warning("Expected 10 validation repeats; found ",length(repeats),".")
if (any(vapply(split(FOLDS$fold,FOLDS$repeat),function(x) uniqueN(x),integer(1))!=5L)) {
    warning("At least one repeat does not contain exactly five folds.")
}

scale_pair <- function(train,test) {
    train <- as.matrix(train); test <- as.matrix(test)
    mu <- colMeans(train); ss <- apply(train,2,sd)
    ss[!is.finite(ss)|ss<=0] <- 1
    list(
        train=sweep(sweep(train,2,mu,"-"),2,ss,"/"),
        test=sweep(sweep(test,2,mu,"-"),2,ss,"/")
    )
}

make_inner_folds <- function(cancer,event,k=5L,seed=1L) {
    set.seed(seed)
    strata <- interaction(cancer,event,drop=TRUE,lex.order=TRUE)
    out <- integer(length(strata))
    for (s in levels(strata)) {
        ii <- which(strata==s)
        labels <- rep(seq_len(k),length.out=length(ii))
        out[ii] <- sample(labels,length(labels),replace=FALSE)
    }
    out
}

baseline_predict <- function(tr,te) {
    tr <- copy(tr); te <- copy(te)
    tr[, burden_z:=zsafe(burden_total)]
    tr[, age_z:=zsafe(age)]
    bmu <- mean(tr$burden_total); bsd <- sd(tr$burden_total)
    amu <- mean(tr$age); asd <- sd(tr$age)
    te[, burden_z:=(burden_total-bmu)/bsd]
    te[, age_z:=(age-amu)/asd]
    fit <- coxph(
        Surv(pfi_time,pfi_event)~burden_z+age_z+strata(cancer),
        data=tr,ties="efron",x=TRUE
    )
    as.numeric(predict(fit,newdata=te,type="lp"))
}

penalized_predict <- function(tr,te,candidates,inner_seed) {
    base_num <- scale_pair(
        tr[,.(burden_total,age)],
        te[,.(burden_total,age)]
    )

    lev <- sort(unique(tr$cancer))
    ctr <- model.matrix(~factor(tr$cancer,levels=lev))[,-1L,drop=FALSE]
    cte <- model.matrix(~factor(te$cancer,levels=lev))[,-1L,drop=FALSE]
    colnames(ctr) <- paste0("cancer_",lev[-1L])
    colnames(cte) <- colnames(ctr)

    C <- scale_pair(
        as.matrix(tr[,..candidates]),
        as.matrix(te[,..candidates])
    )

    Xtr <- cbind(
        burden_total=base_num$train[,"burden_total"],
        age=base_num$train[,"age"],
        ctr,
        C$train
    )
    Xte <- cbind(
        burden_total=base_num$test[,"burden_total"],
        age=base_num$test[,"age"],
        cte,
        C$test
    )

    nbase <- 2L+ncol(ctr)
    penalty <- c(rep(0,nbase),rep(1,length(candidates)))
    inner <- make_inner_folds(tr$cancer,tr$pfi_event,5L,inner_seed)

    fit <- cv.glmnet(
        Xtr,
        Surv(tr$pfi_time,tr$pfi_event),
        family="cox",
        alpha=0.5,
        foldid=inner,
        penalty.factor=penalty,
        standardize=FALSE,
        type.measure="deviance"
    )

    cf <- as.matrix(coef(fit,s="lambda.min"))
    list(
        pred=as.numeric(predict(fit,newx=Xte,s="lambda.min",type="link")),
        coef=cf,
        selected=rownames(cf)[abs(cf[,1])>1e-10],
        lambda=fit$lambda.min
    )
}

models <- list(
    UCEI="UCEI_B9C",
    cDITHER="cDITHER",
    Steele21=steele,
    Drews17=drews,
    AllSOTA39=sota,
    AllSOTA39_plus_UCEI=c(sota,"UCEI_B9C")
)

FOLD_RES <- list()
PRED_RES <- list()
SEL_RES <- list()
counter <- 0L

for (rr in repeats) {
    fr <- FOLDS[repeat==rr]
    X <- merge(D,fr[,.(patient_id,fold)],by="patient_id",all=FALSE)

    for (ff in sort(unique(X$fold))) {
        tr <- X[fold!=ff]
        te <- X[fold==ff]

        counter <- counter+1L
        pb <- baseline_predict(tr,te)

        FOLD_RES[[length(FOLD_RES)+1L]] <- data.table(
            repeat=rr,fold=ff,model="Clinical_baseline",
            n=nrow(te),events=sum(te$pfi_event),
            C=cindex(te$pfi_time,te$pfi_event,pb)
        )
        PRED_RES[[length(PRED_RES)+1L]] <- data.table(
            patient_id=te$patient_id,repeat=rr,fold=ff,
            model="Clinical_baseline",pfi_time=te$pfi_time,
            pfi_event=te$pfi_event,lp=pb
        )

        for (mn in names(models)) {
            cand <- models[[mn]]
            fit <- penalized_predict(
                tr,te,cand,
                inner_seed=100000L + as.integer(rr)*100L + as.integer(ff)
            )

            FOLD_RES[[length(FOLD_RES)+1L]] <- data.table(
                repeat=rr,fold=ff,model=mn,
                n=nrow(te),events=sum(te$pfi_event),
                C=cindex(te$pfi_time,te$pfi_event,fit$pred)
            )
            PRED_RES[[length(PRED_RES)+1L]] <- data.table(
                patient_id=te$patient_id,repeat=rr,fold=ff,
                model=mn,pfi_time=te$pfi_time,
                pfi_event=te$pfi_event,lp=fit$pred
            )

            if (mn=="AllSOTA39_plus_UCEI") {
                candidate_selected <- intersect(
                    fit$selected,
                    c(sota,"UCEI_B9C")
                )
                SEL_RES[[length(SEL_RES)+1L]] <- data.table(
                    repeat=rr,fold=ff,
                    predictor=c(sota,"UCEI"),
                    selected=as.integer(c(sota,"UCEI_B9C") %in% candidate_selected)
                )
            }
        }
    }
}

B9C_fold_results <- rbindlist(FOLD_RES)
B9C_predictions <- rbindlist(PRED_RES)
B9C_selection_results <- rbindlist(SEL_RES)

B9C_repeat_results <- B9C_predictions[
    ,
    .(C=cindex(pfi_time,pfi_event,lp)),
    by=.(repeat,model)
]

B9C_summary <- B9C_repeat_results[
    ,
    .(mean_C=mean(C),SD=sd(C),min_C=min(C),max_C=max(C)),
    by=model
][order(-mean_C)]

W <- dcast(B9C_repeat_results,repeat~model,value.var="C")

B9C_delta <- rbindlist(list(
    W[,.(repeat,comparison="combo vs AllSOTA39",
         delta_C=AllSOTA39_plus_UCEI-AllSOTA39)],
    W[,.(repeat,comparison="combo vs UCEI",
         delta_C=AllSOTA39_plus_UCEI-UCEI)],
    W[,.(repeat,comparison="UCEI vs AllSOTA39",
         delta_C=UCEI-AllSOTA39)],
    W[,.(repeat,comparison="UCEI vs cDITHER",
         delta_C=UCEI-cDITHER)],
    W[,.(repeat,comparison="UCEI vs Drews17",
         delta_C=UCEI-Drews17)],
    W[,.(repeat,comparison="UCEI vs Steele21",
         delta_C=UCEI-Steele21)]
))

B9C_delta_summary <- B9C_delta[
    ,
    .(
        mean_delta_C=mean(delta_C),
        SD=sd(delta_C),
        positive_repeats=sum(delta_C>0),
        total_repeats=.N
    ),
    by=comparison
]

B9C_selection_frequency <- B9C_selection_results[
    ,
    .(outer_folds_selected=sum(selected),total_outer_folds=.N),
    by=predictor
][order(-outer_folds_selected,predictor)]

############################################################
## Locked numerical audit
############################################################

LOCKED <- data.table(
    model=c(
        "AllSOTA39_plus_UCEI","UCEI","AllSOTA39",
        "cDITHER","Drews17","Steele21","Clinical_baseline"
    ),
    expected_mean=c(
        0.72610,0.72458,0.71985,
        0.71809,0.71565,0.71255,0.70468
    )
)

AUDIT <- merge(
    LOCKED,
    B9C_summary[,.(model,observed_mean=mean_C)],
    by="model",
    all.x=TRUE
)
AUDIT[,abs_diff:=abs(observed_mean-expected_mean)]
AUDIT[,pass:=is.finite(abs_diff)&abs_diff<5e-4]

############################################################
## Export
############################################################

fwrite(B9_LRT,file.path(output_dir,"B9_Steele_UCEI_bidirectional_LRT.csv"))
fwrite(B9C_fold_results,file.path(output_dir,"B9C_fold_results.csv"))
fwrite(B9C_repeat_results,file.path(output_dir,"B9C_repeat_results.csv"))
fwrite(B9C_summary,file.path(output_dir,"B9C_summary.csv"))
fwrite(B9C_delta,file.path(output_dir,"B9C_delta_by_repeat.csv"))
fwrite(B9C_delta_summary,file.path(output_dir,"B9C_delta_summary.csv"))
fwrite(B9C_selection_results,file.path(output_dir,"B9C_selection_results.csv"))
fwrite(B9C_selection_frequency,file.path(output_dir,"B9C_selection_frequency.csv"))
fwrite(AUDIT,file.path(output_dir,"B9C_locked_headline_audit.csv"))

cat("\nB9/B9C complete.\n")
print(B9_LRT)
print(B9C_summary)
print(B9C_delta_summary)
print(B9C_selection_frequency[1:min(.N,15)])
print(AUDIT)
