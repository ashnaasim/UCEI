############################################################
## B12C / B12C2: UCEI CONSTRUCTION ABLATION AND
## HIGHER-DIMENSIONAL ARCHITECTURE ALTERNATIVES
##
## Publication-facing reconstruction from the final Stage-3
## B12C/B12C2 records.
##
## IMPORTANT
## ---------
## The final result tables and model definitions are preserved
## in the project archive, but the original standalone B12C
## script did not survive on disk. This file reconstructs the
## analysis described in the manuscript and retained outputs.
## It MUST be checked against the historical CSVs before the
## repository is made public.
##
## Primary B12C variants:
##   1. full_UCEI
##   2. six leave-one-feature-out UCEI reconstructions
##   3. PCA_without_burden_residualization
##   4. equal_weight_residualized
##   5. six_residualized_components_joint
##
## B12C2 higher-dimensional comparison:
##   - AllSOTA39
##   - AllSOTA39 + training-fold UCEI
##   - AllSOTA39 + six residualized architecture components
##   - AllSOTA39 + five components excluding state-size
##
## Validation:
##   10 repeats x 5 outer folds
##   folds stratified by cancer x event
##   training-only score construction in every outer fold
############################################################

suppressPackageStartupMessages({
    library(data.table)
    library(survival)
    library(glmnet)
})

args <- commandArgs(trailingOnly=TRUE)

input_file <- if (length(args)>=1L) args[[1L]] else
    Sys.getenv("UCEI_B12C_INPUT", unset="")

output_dir <- if (length(args)>=2L) args[[2L]] else
    file.path("results","benchmark","B12C_construction_ablation")

if (!nzchar(input_file) || !file.exists(input_file)) {
    stop("Supply the merged B12C benchmark table as argument 1 or set UCEI_B12C_INPUT.")
}

dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)

B12C_features <- c(
    "jsd_neutral",
    "ent_state_size",
    "ent_chr_sd",
    "ent_chr_mean",
    "ent_ampclass",
    "ent_transition"
)

B12C_SOTA <- c(
    paste0("CN",1:21),
    paste0("CX",1:17),
    "cDITHER"
)

B12C_required <- c(
    "patient_id","cancer","pfi_time","pfi_event",
    "age","burden_total",
    B12C_features,B12C_SOTA
)

B12C_D <- if (grepl("\\.rds$",input_file,ignore.case=TRUE)) {
    as.data.table(readRDS(input_file))
} else {
    fread(input_file)
}

miss <- setdiff(B12C_required,names(B12C_D))
if (length(miss)) stop("Missing required columns: ",paste(miss,collapse=", "))

B12C_D <- B12C_D[complete.cases(B12C_D[,..B12C_required])]
B12C_D[, cancer:=as.character(cancer)]
for (nm in setdiff(B12C_required,c("patient_id","cancer"))) {
    set(B12C_D,j=nm,value=as.numeric(B12C_D[[nm]]))
}
B12C_D[, pfi_event:=as.integer(pfi_event)]

if (anyDuplicated(B12C_D$patient_id)) stop("Duplicate patient_id values in B12C input.")

############################################################
## Helpers
############################################################

B12C_Cindex <- function(time,event,risk) {
    unname(survival::concordance(
        Surv(time,event)~risk,
        reverse=TRUE
    )$concordance)
}

B12C_make_folds <- function(cancer,event,k=5L,seed) {
    set.seed(seed)
    strata <- interaction(cancer,event,drop=TRUE,lex.order=TRUE)
    fold <- integer(length(strata))
    for (s in levels(strata)) {
        ii <- which(strata==s)
        lab <- rep(seq_len(k),length.out=length(ii))
        fold[ii] <- sample(lab,length(lab),replace=FALSE)
    }
    fold
}

B12C_scale_pair <- function(train,test) {
    train <- as.matrix(train)
    test <- as.matrix(test)
    mu <- colMeans(train,na.rm=TRUE)
    ss <- apply(train,2L,sd,na.rm=TRUE)
    ss[!is.finite(ss)|ss<=0] <- 1
    list(
        train=sweep(sweep(train,2L,mu,"-"),2L,ss,"/"),
        test=sweep(sweep(test,2L,mu,"-"),2L,ss,"/"),
        center=mu,scale=ss
    )
}

B12C_residual_matrices <- function(train,test,features=B12C_features) {
    Ztr <- matrix(NA_real_,nrow(train),length(features),
                  dimnames=list(NULL,features))
    Zte <- matrix(NA_real_,nrow(test),length(features),
                  dimnames=list(NULL,features))

    for (cc in sort(unique(train$cancer))) {
        itr <- which(train$cancer==cc)
        ite <- which(test$cancer==cc)
        if (!length(itr)) next

        Xtr <- cbind(1,train$burden_total[itr])
        Xte <- cbind(1,test$burden_total[ite])

        for (f in features) {
            fit <- lm.fit(Xtr,train[[f]][itr])
            cf <- as.numeric(fit$coefficients)
            rtr <- train[[f]][itr]-as.numeric(Xtr%*%cf)
            mu <- mean(rtr,na.rm=TRUE)
            ss <- sd(rtr,na.rm=TRUE)
            if (!is.finite(ss)||ss<=0) stop("Invalid residual SD: ",cc," / ",f)
            Ztr[itr,f] <- (rtr-mu)/ss
            if (length(ite)) {
                rte <- test[[f]][ite]-as.numeric(Xte%*%cf)
                Zte[ite,f] <- (rte-mu)/ss
            }
        }
    }

    if (any(!is.finite(Ztr))||any(!is.finite(Zte))) {
        stop("Non-finite residual matrices.")
    }
    list(train=Ztr,test=Zte)
}

B12C_raw_standardized_matrices <- function(train,test,features=B12C_features) {
    ## No burden regression; feature centering/scaling is learned
    ## exclusively from the outer training partition.
    B12C_scale_pair(
        as.matrix(train[,..features]),
        as.matrix(test[,..features])
    )
}

B12C_pca_score <- function(Ztr,Zte,orientation_target=NULL) {
    pc <- prcomp(Ztr,center=FALSE,scale.=FALSE)
    load <- pc$rotation[,1L]
    tr <- as.numeric(Ztr%*%load)
    te <- as.numeric(Zte%*%load)

    if (!is.null(orientation_target)) {
        cc <- suppressWarnings(cor(tr,orientation_target,use="complete.obs"))
        if (is.finite(cc)&&cc<0) {
            load <- -load
            tr <- -tr
            te <- -te
        }
    }

    mu <- mean(tr)
    ss <- sd(tr)
    list(
        train=(tr-mu)/ss,
        test=(te-mu)/ss,
        raw_train=tr,
        raw_test=te,
        loadings=load,
        pc1_variance=pc$sdev[1L]^2/sum(pc$sdev^2)
    )
}

B12C_cox_score <- function(train,test,score_train,score_test) {
    tr <- copy(train)
    te <- copy(test)
    tr[, score:=score_train]
    te[, score:=score_test]

    fit <- coxph(
        Surv(pfi_time,pfi_event)~
            score + burden_total + age + strata(cancer),
        data=tr,
        ties="efron",
        x=TRUE
    )

    pred <- as.numeric(predict(fit,newdata=te,type="lp"))
    b <- unname(coef(fit)["score"])

    list(
        C=B12C_Cindex(te$pfi_time,te$pfi_event,pred),
        beta=b,
        HR=exp(b),
        pred=pred
    )
}

B12C_cox_block <- function(train,test,Ztr,Zte) {
    tr <- copy(train)
    te <- copy(test)

    zz <- colnames(Ztr)
    for (j in seq_along(zz)) {
        tr[[paste0("arch_",j)]] <- Ztr[,j]
        te[[paste0("arch_",j)]] <- Zte[,j]
    }

    arch_terms <- paste0("arch_",seq_along(zz))

    f <- as.formula(
        paste0(
            "Surv(pfi_time,pfi_event) ~ ",
            paste(arch_terms,collapse=" + "),
            " + burden_total + age + strata(cancer)"
        )
    )

    fit <- coxph(f,data=tr,ties="efron",x=TRUE)
    pred <- as.numeric(predict(fit,newdata=te,type="lp"))

    list(
        C=B12C_Cindex(te$pfi_time,te$pfi_event,pred),
        pred=pred
    )
}

B12C_fit_penalized <- function(Xtr,Xte,ytr,penalty,foldid) {
    cv <- cv.glmnet(
        Xtr,ytr,
        family="cox",
        alpha=0.5,
        foldid=foldid,
        penalty.factor=penalty,
        type.measure="deviance",
        standardize=FALSE,
        nlambda=100
    )
    cf <- as.matrix(coef(cv,s="lambda.min"))
    list(
        pred=as.numeric(predict(cv,newx=Xte,s="lambda.min",type="link")),
        coef=cf,
        lambda=cv$lambda.min
    )
}

############################################################
## Repeated outer validation
############################################################

B12C_fold_rows <- list()
B12C_beta_rows <- list()
B12C_corr_rows <- list()
B12C_pred_rows <- list()

B12C2_pred_rows <- list()
B12C2_select_rows <- list()

repeats <- 10L
outer_k <- 5L

for (rr in seq_len(repeats)) {

    outer <- B12C_make_folds(
        B12C_D$cancer,
        B12C_D$pfi_event,
        k=outer_k,
        seed=20260813+rr
    )

    for (ff in seq_len(outer_k)) {

        tr <- B12C_D[outer!=ff]
        te <- B12C_D[outer==ff]

        R <- B12C_residual_matrices(tr,te,B12C_features)

        orient_target <- rowMeans(
            cbind(tr$jsd_neutral,tr$ent_state_size),
            na.rm=TRUE
        )

        full <- B12C_pca_score(
            R$train,R$test,
            orientation_target=orient_target
        )

        variants <- list(
            full_UCEI=full
        )

        for (omit in B12C_features) {
            keep <- setdiff(B12C_features,omit)
            variants[[paste0("LOFO_minus_",omit)]] <- B12C_pca_score(
                R$train[,keep,drop=FALSE],
                R$test[,keep,drop=FALSE],
                orientation_target=orient_target
            )
        }

        RAW <- B12C_raw_standardized_matrices(
            tr,te,B12C_features
        )
        variants[["PCA_without_burden_residualization"]] <-
            B12C_pca_score(
                RAW$train,RAW$test,
                orientation_target=orient_target
            )

        ## Equal-magnitude weighting with signs retained from the
        ## training-fold full-PC1 orientation.
        sgn <- sign(full$loadings[B12C_features])
        sgn[sgn==0] <- 1
        eq_tr_raw <- as.numeric(R$train%*%(sgn/sqrt(length(sgn))))
        eq_te_raw <- as.numeric(R$test%*%(sgn/sqrt(length(sgn))))
        eq_mu <- mean(eq_tr_raw)
        eq_sd <- sd(eq_tr_raw)
        variants[["equal_weight_residualized"]] <- list(
            train=(eq_tr_raw-eq_mu)/eq_sd,
            test=(eq_te_raw-eq_mu)/eq_sd
        )

        ########################################################
        ## One-dimensional variants
        ########################################################

        for (mn in names(variants)) {
            v <- variants[[mn]]
            z <- B12C_cox_score(
                tr,te,
                v$train,
                v$test
            )

            B12C_fold_rows[[length(B12C_fold_rows)+1L]] <- data.table(
                rep_id=rr,fold=ff,model=mn,
                n=nrow(te),events=sum(te$pfi_event),C=z$C
            )

            B12C_beta_rows[[length(B12C_beta_rows)+1L]] <- data.table(
                rep_id=rr,fold=ff,model=mn,
                beta=z$beta,HR=z$HR
            )

            B12C_corr_rows[[length(B12C_corr_rows)+1L]] <- data.table(
                rep_id=rr,fold=ff,model=mn,
                train_cor_to_full=suppressWarnings(
                    cor(v$train,full$train,use="complete.obs")
                ),
                test_spearman_to_full=suppressWarnings(
                    cor(v$test,full$test,method="spearman",use="complete.obs")
                )
            )

            B12C_pred_rows[[length(B12C_pred_rows)+1L]] <- data.table(
                patient_id=te$patient_id,
                pfi_time=te$pfi_time,
                pfi_event=te$pfi_event,
                rep_id=rr,fold=ff,model=mn,
                linear_predictor=z$pred
            )
        }

        ########################################################
        ## Six residualized architecture components jointly
        ########################################################

        z6 <- B12C_cox_block(
            tr,te,
            R$train,R$test
        )

        B12C_fold_rows[[length(B12C_fold_rows)+1L]] <- data.table(
            rep_id=rr,fold=ff,
            model="six_residualized_components_joint",
            n=nrow(te),events=sum(te$pfi_event),C=z6$C
        )

        B12C_pred_rows[[length(B12C_pred_rows)+1L]] <- data.table(
            patient_id=te$patient_id,
            pfi_time=te$pfi_time,
            pfi_event=te$pfi_event,
            rep_id=rr,fold=ff,
            model="six_residualized_components_joint",
            linear_predictor=z6$pred
        )

        ########################################################
        ## B12C2: higher-dimensional architecture + AllSOTA39
        ########################################################

        base_scaled <- B12C_scale_pair(
            as.matrix(tr[,.(burden_total,age)]),
            as.matrix(te[,.(burden_total,age)])
        )

        lev <- sort(unique(tr$cancer))
        ctr <- model.matrix(~factor(tr$cancer,levels=lev))[,-1L,drop=FALSE]
        cte <- model.matrix(~factor(te$cancer,levels=lev))[,-1L,drop=FALSE]
        colnames(ctr) <- paste0("cancer_",lev[-1L])
        colnames(cte) <- colnames(ctr)

        Xbase_tr <- cbind(
            burden_total=base_scaled$train[,"burden_total"],
            age=base_scaled$train[,"age"],
            ctr
        )
        Xbase_te <- cbind(
            burden_total=base_scaled$test[,"burden_total"],
            age=base_scaled$test[,"age"],
            cte
        )

        S <- B12C_scale_pair(
            as.matrix(tr[,..B12C_SOTA]),
            as.matrix(te[,..B12C_SOTA])
        )

        inner <- B12C_make_folds(
            tr$cancer,tr$pfi_event,
            k=5L,
            seed=20260813+rr*100L+ff
        )

        ytr <- Surv(tr$pfi_time,tr$pfi_event)

        arch5 <- setdiff(B12C_features,"ent_state_size")

        model_mats <- list(
            AllSOTA39=list(
                tr=cbind(Xbase_tr,S$train),
                te=cbind(Xbase_te,S$test),
                penalty=c(rep(0,ncol(Xbase_tr)),rep(1,length(B12C_SOTA)))
            ),
            AllSOTA39_plus_strict_UCEI=list(
                tr=cbind(Xbase_tr,S$train,UCEI=full$train),
                te=cbind(Xbase_te,S$test,UCEI=full$test),
                penalty=c(rep(0,ncol(Xbase_tr)),rep(1,length(B12C_SOTA)),1)
            ),
            AllSOTA39_plus_ARCH6=list(
                tr=cbind(Xbase_tr,S$train,R$train),
                te=cbind(Xbase_te,S$test,R$test),
                penalty=c(rep(0,ncol(Xbase_tr)),rep(1,length(B12C_SOTA)+6L))
            ),
            AllSOTA39_plus_ARCH5_no_state_size=list(
                tr=cbind(Xbase_tr,S$train,R$train[,arch5,drop=FALSE]),
                te=cbind(Xbase_te,S$test,R$test[,arch5,drop=FALSE]),
                penalty=c(rep(0,ncol(Xbase_tr)),rep(1,length(B12C_SOTA)+5L))
            )
        )

        for (mn in names(model_mats)) {
            mm <- model_mats[[mn]]
            fit <- B12C_fit_penalized(
                mm$tr,mm$te,ytr,mm$penalty,inner
            )

            B12C2_pred_rows[[length(B12C2_pred_rows)+1L]] <- data.table(
                patient_id=te$patient_id,
                pfi_time=te$pfi_time,
                pfi_event=te$pfi_event,
                rep_id=rr,fold=ff,model=mn,
                linear_predictor=fit$pred
            )

            nz <- rownames(fit$coef)[abs(fit$coef[,1L])>1e-10]

            if (mn %in% c(
                "AllSOTA39_plus_ARCH6",
                "AllSOTA39_plus_ARCH5_no_state_size"
            )) {
                B12C2_select_rows[[length(B12C2_select_rows)+1L]] <- data.table(
                    rep_id=rr,fold=ff,model=mn,
                    selected=paste(nz,collapse=";")
                )
            }
        }
    }
}

############################################################
## B12C summaries
############################################################

B12C_FOLD <- rbindlist(B12C_fold_rows)
B12C_BETA <- rbindlist(B12C_beta_rows)
B12C_CORR <- rbindlist(B12C_corr_rows)
B12C_PRED <- rbindlist(B12C_pred_rows)

B12C_REPEAT <- B12C_PRED[
    ,
    .(C=B12C_Cindex(pfi_time,pfi_event,linear_predictor)),
    by=.(rep_id,model)
]

B12C_SUMMARY <- B12C_REPEAT[
    ,
    .(
        mean_C=mean(C),
        sd_C=sd(C),
        min_C=min(C),
        max_C=max(C)
    ),
    by=model
][order(-mean_C)]

B12C_WIDE <- dcast(
    B12C_REPEAT,
    rep_id~model,
    value.var="C"
)

variant_names <- setdiff(
    names(B12C_WIDE),
    c("rep_id","full_UCEI")
)

B12C_DELTA <- rbindlist(lapply(
    variant_names,
    function(mn) {
        data.table(
            rep_id=B12C_WIDE$rep_id,
            model=mn,
            C=B12C_WIDE[[mn]],
            full_C=B12C_WIDE$full_UCEI,
            delta_vs_full=B12C_WIDE[[mn]]-B12C_WIDE$full_UCEI
        )
    }
))

B12C_DELTA_SUMMARY <- B12C_DELTA[
    ,
    .(
        mean_delta_vs_full=mean(delta_vs_full),
        sd_delta=sd(delta_vs_full),
        min_delta=min(delta_vs_full),
        max_delta=max(delta_vs_full),
        better_than_full_repeats=sum(delta_vs_full>0),
        worse_than_full_repeats=sum(delta_vs_full<0),
        total_repeats=.N
    ),
    by=model
][order(-mean_delta_vs_full)]

B12C_BETA_SUMMARY <- B12C_BETA[
    ,
    .(
        folds=.N,
        positive_beta_folds=sum(beta>0),
        median_HR=median(HR),
        min_HR=min(HR),
        max_HR=max(HR)
    ),
    by=model
]

B12C_CORR_SUMMARY <- B12C_CORR[
    ,
    .(
        mean_train_cor=mean(train_cor_to_full),
        min_train_cor=min(train_cor_to_full),
        mean_test_spearman=mean(test_spearman_to_full),
        min_test_spearman=min(test_spearman_to_full)
    ),
    by=model
]

############################################################
## B12C2 summaries
############################################################

B12C2_PRED <- rbindlist(B12C2_pred_rows)

B12C2_REPEAT <- B12C2_PRED[
    ,
    .(C=B12C_Cindex(pfi_time,pfi_event,linear_predictor)),
    by=.(rep_id,model)
]

B12C2_SUMMARY <- B12C2_REPEAT[
    ,
    .(
        mean_C=mean(C),
        sd_C=sd(C),
        min_C=min(C),
        max_C=max(C)
    ),
    by=model
][order(-mean_C)]

B12C2_WIDE <- dcast(
    B12C2_REPEAT,
    rep_id~model,
    value.var="C"
)

B12C2_DELTA <- rbindlist(list(
    B12C2_WIDE[,.(rep_id,comparison="UCEI_combo - AllSOTA39",
                    delta_C=AllSOTA39_plus_strict_UCEI-AllSOTA39)],
    B12C2_WIDE[,.(rep_id,comparison="ARCH6_combo - AllSOTA39",
                    delta_C=AllSOTA39_plus_ARCH6-AllSOTA39)],
    B12C2_WIDE[,.(rep_id,comparison="ARCH5_combo - AllSOTA39",
                    delta_C=AllSOTA39_plus_ARCH5_no_state_size-AllSOTA39)],
    B12C2_WIDE[,.(rep_id,comparison="ARCH6_combo - UCEI_combo",
                    delta_C=AllSOTA39_plus_ARCH6-AllSOTA39_plus_strict_UCEI)],
    B12C2_WIDE[,.(rep_id,comparison="ARCH5_combo - UCEI_combo",
                    delta_C=AllSOTA39_plus_ARCH5_no_state_size-AllSOTA39_plus_strict_UCEI)],
    B12C2_WIDE[,.(rep_id,comparison="ARCH6_combo - ARCH5_combo",
                    delta_C=AllSOTA39_plus_ARCH6-AllSOTA39_plus_ARCH5_no_state_size)]
))

B12C2_DELTA_SUMMARY <- B12C2_DELTA[
    ,
    .(
        mean_delta_C=mean(delta_C),
        sd_delta_C=sd(delta_C),
        min_delta_C=min(delta_C),
        max_delta_C=max(delta_C),
        positive_repeats=sum(delta_C>0),
        negative_repeats=sum(delta_C<0),
        total_repeats=.N
    ),
    by=comparison
]

############################################################
## Export
############################################################

fwrite(B12C_BETA,
       file.path(output_dir,"B12C_coefficient_stability.csv"))
fwrite(B12C_DELTA,
       file.path(output_dir,"B12C_delta_vs_full_UCEI.csv"))
fwrite(B12C_FOLD,
       file.path(output_dir,"B12C_fold_performance.csv"))
fwrite(B12C_DELTA_SUMMARY,
       file.path(output_dir,"B12C_KEY_design_ablations.csv"))
fwrite(B12C_SUMMARY,
       file.path(output_dir,"B12C_performance_summary.csv"))
fwrite(B12C_REPEAT,
       file.path(output_dir,"B12C_repeat_performance.csv"))
fwrite(B12C_CORR_SUMMARY,
       file.path(output_dir,"B12C_score_similarity_summary.csv"))

fwrite(B12C2_SUMMARY,
       file.path(output_dir,"B12C2_higher_dimensional_summary.csv"))
fwrite(B12C2_DELTA,
       file.path(output_dir,"B12C2_delta_by_repeat.csv"))
fwrite(B12C2_DELTA_SUMMARY,
       file.path(output_dir,"B12C2_delta_summary.csv"))

cat("\nB12C/B12C2 complete.\n")
print(B12C_SUMMARY)
print(B12C_DELTA_SUMMARY)
print(B12C2_SUMMARY)
print(B12C2_DELTA_SUMMARY)
