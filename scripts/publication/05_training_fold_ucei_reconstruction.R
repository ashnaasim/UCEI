############################################################
## B12B: COMPLETE TRAINING-FOLD UCEI RECONSTRUCTION
##
## Publication-facing reconstruction of the final executed
## Stage-3 B12B analysis.
##
## Historical analysis target:
##   n = 1,574 complete cases
##   events = 439
##   cancers = 9
##
## Models:
##   baseline
##   training-fold UCEI
##   AllSOTA39
##   AllSOTA39 + training-fold UCEI
##
## Recovered final settings:
##   10 repeats x 5 outer folds
##   outer folds stratified by cancer x event
##   5-fold inner CV
##   glmnet Cox elastic net, alpha = 0.5
##   type.measure = "deviance"
##   lambda.min
##   nlambda = 100
##   standardize = FALSE
##   baseline covariates unpenalized
##   candidate predictors penalized
##   outer seed = 20260813 + repeat
##   inner seed = 20260813 + repeat*100 + fold
##
## Required input:
## A merged patient-level table containing:
##   patient_id, cancer, pfi_time, pfi_event, age,
##   burden_total,
##   six raw UCEI architecture features,
##   CN1-CN21, CX1-CX17, cDITHER,
##   and optionally a reference/frozen UCEI column for audit.
##
## Usage:
## Rscript 05_training_fold_ucei_reconstruction.R input.csv output_dir
############################################################

suppressPackageStartupMessages({
    library(data.table)
    library(survival)
    library(glmnet)
})

args <- commandArgs(trailingOnly = TRUE)

input_file <- if (length(args) >= 1L) args[[1L]] else
    Sys.getenv("UCEI_B12B_INPUT", unset = "")

output_dir <- if (length(args) >= 2L) args[[2L]] else
    file.path("results", "benchmark", "B12B_train_fold_reconstruction")

if (!nzchar(input_file) || !file.exists(input_file)) {
    stop("Supply the merged B12B benchmark table as argument 1 or set UCEI_B12B_INPUT.")
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

B12B_features <- c(
    "jsd_neutral",
    "ent_state_size",
    "ent_chr_sd",
    "ent_chr_mean",
    "ent_ampclass",
    "ent_transition"
)

B12B_SOTA <- c(
    paste0("CN", 1:21),
    paste0("CX", 1:17),
    "cDITHER"
)

B12B_required <- c(
    "patient_id",
    "cancer",
    "pfi_time",
    "pfi_event",
    "age",
    "burden_total",
    B12B_features,
    B12B_SOTA
)

B12B_D <- if (grepl("\\.rds$", input_file, ignore.case = TRUE)) {
    as.data.table(readRDS(input_file))
} else {
    fread(input_file)
}

missing_cols <- setdiff(B12B_required, names(B12B_D))
if (length(missing_cols)) {
    stop("Missing required B12B columns: ", paste(missing_cols, collapse = ", "))
}

B12B_D <- B12B_D[complete.cases(B12B_D[, ..B12B_required])]
B12B_D[, cancer := as.character(cancer)]
B12B_D[, pfi_time := as.numeric(pfi_time)]
B12B_D[, pfi_event := as.integer(pfi_event)]
B12B_D[, age := as.numeric(age)]
B12B_D[, burden_total := as.numeric(burden_total)]

for (nm in c(B12B_features, B12B_SOTA)) {
    set(B12B_D, j = nm, value = as.numeric(B12B_D[[nm]]))
}

if (anyDuplicated(B12B_D$patient_id)) {
    stop("B12B input contains duplicate patient_id values.")
}

############################################################
## Helpers
############################################################

B12B_Cindex <- function(time, event, risk) {
    unname(
        survival::concordance(
            Surv(time, event) ~ risk,
            reverse = TRUE
        )$concordance
    )
}

B12B_make_folds <- function(cancer, event, k = 5L, seed) {
    set.seed(seed)

    strata <- interaction(
        cancer,
        event,
        drop = TRUE,
        lex.order = TRUE
    )

    fold <- integer(length(strata))

    for (s in levels(strata)) {
        ii <- which(strata == s)

        lab <- rep(
            seq_len(k),
            length.out = length(ii)
        )

        fold[ii] <- sample(
            lab,
            length(lab),
            replace = FALSE
        )
    }

    fold
}

B12B_train_scale <- function(train, test) {
    train <- as.matrix(train)
    test <- as.matrix(test)

    mu <- colMeans(train, na.rm = TRUE)
    ss <- apply(train, 2L, sd, na.rm = TRUE)

    bad <- !is.finite(ss) | ss <= 0
    ss[bad] <- 1

    train_z <- sweep(
        sweep(train, 2L, mu, "-"),
        2L, ss, "/"
    )

    test_z <- sweep(
        sweep(test, 2L, mu, "-"),
        2L, ss, "/"
    )

    list(
        train = train_z,
        test = test_z,
        center = mu,
        scale = ss
    )
}

B12B_construct_ucei <- function(train, test) {
    train <- as.data.table(train)
    test <- as.data.table(test)

    Ztr <- matrix(
        NA_real_,
        nrow = nrow(train),
        ncol = length(B12B_features),
        dimnames = list(NULL, B12B_features)
    )

    Zte <- matrix(
        NA_real_,
        nrow = nrow(test),
        ncol = length(B12B_features),
        dimnames = list(NULL, B12B_features)
    )

    cancer_levels <- sort(unique(train$cancer))

    if (any(!test$cancer %in% cancer_levels)) {
        stop("Outer test fold contains a cancer absent from its training fold.")
    }

    parameter_rows <- vector(
        "list",
        length(cancer_levels) * length(B12B_features)
    )
    pr <- 0L

    for (cc in cancer_levels) {
        itr <- which(train$cancer == cc)
        ite <- which(test$cancer == cc)

        Xtr <- cbind(
            intercept = 1,
            burden_total = train$burden_total[itr]
        )

        Xte <- cbind(
            intercept = 1,
            burden_total = test$burden_total[ite]
        )

        for (f in B12B_features) {
            ytr <- train[[f]][itr]

            fit <- lm.fit(
                x = Xtr,
                y = ytr
            )

            coef <- as.numeric(fit$coefficients)
            if (length(coef) != 2L || any(!is.finite(coef))) {
                stop("Invalid burden-regression coefficients for ", cc, " / ", f)
            }

            rtr <- ytr - as.numeric(Xtr %*% coef)

            rmu <- mean(rtr, na.rm = TRUE)
            rsd <- sd(rtr, na.rm = TRUE)

            if (!is.finite(rsd) || rsd <= 0) {
                stop("Invalid residual SD for ", cc, " / ", f)
            }

            Ztr[itr, f] <- (rtr - rmu) / rsd

            if (length(ite)) {
                rte <- test[[f]][ite] - as.numeric(Xte %*% coef)
                Zte[ite, f] <- (rte - rmu) / rsd
            }

            pr <- pr + 1L
            parameter_rows[[pr]] <- data.table(
                cancer = cc,
                feature = f,
                intercept = coef[1L],
                slope = coef[2L],
                residual_mean = rmu,
                residual_sd = rsd,
                n_train = length(itr)
            )
        }
    }

    if (any(!is.finite(Ztr)) || any(!is.finite(Zte))) {
        stop("Non-finite residualized architecture values encountered.")
    }

    pca <- prcomp(
        Ztr,
        center = FALSE,
        scale. = FALSE
    )

    load <- pca$rotation[, 1L]
    raw_tr <- as.numeric(Ztr %*% load)
    raw_te <- as.numeric(Zte %*% load)

    ## Final executed B12B used a training-only architecture
    ## orientation convention, not the global reference loadings.
    orient_target <- rowMeans(
        cbind(
            train$jsd_neutral,
            train$ent_state_size
        ),
        na.rm = TRUE
    )

    orient_cor <- suppressWarnings(
        cor(raw_tr, orient_target, use = "complete.obs")
    )

    if (is.finite(orient_cor) && orient_cor < 0) {
        load <- -load
        raw_tr <- -raw_tr
        raw_te <- -raw_te
    }

    score_center <- mean(raw_tr, na.rm = TRUE)
    score_scale <- sd(raw_tr, na.rm = TRUE)

    if (!is.finite(score_scale) || score_scale <= 0) {
        stop("Invalid training UCEI score scale.")
    }

    list(
        train_score = (raw_tr - score_center) / score_scale,
        test_score = (raw_te - score_center) / score_scale,
        loadings = load,
        pc1_variance = pca$sdev[1L]^2 / sum(pca$sdev^2),
        score_center = score_center,
        score_scale = score_scale,
        parameters = rbindlist(parameter_rows[seq_len(pr)])
    )
}

B12B_fit_predict <- function(
    Xtr,
    Xte,
    ytr,
    penalty,
    inner_fold
) {
    fit <- cv.glmnet(
        x = Xtr,
        y = ytr,
        family = "cox",
        alpha = 0.5,
        foldid = inner_fold,
        penalty.factor = penalty,
        standardize = FALSE,
        type.measure = "deviance",
        nlambda = 100
    )

    beta <- as.matrix(
        coef(
            fit,
            s = "lambda.min"
        )
    )

    pred <- as.numeric(
        predict(
            fit,
            newx = Xte,
            s = "lambda.min",
            type = "link"
        )
    )

    list(
        pred = pred,
        lambda = fit$lambda.min,
        coef = beta
    )
}

############################################################
## Reference loadings for audit only
############################################################

ref_loading_file <- file.path(
    "parameters",
    "reference",
    "UCEI_v2_frozen_PCA_loadings.csv"
)

B12B_frozen_loading <- NULL

if (file.exists(ref_loading_file)) {
    ref <- fread(ref_loading_file)
    B12B_frozen_loading <- setNames(
        ref$loading,
        sub("_v2_resid_z$", "", ref$feature)
    )
}

frozen_score_candidates <- c(
    "UCEI_v2_z",
    "UCEI_z",
    "UCEI"
)

frozen_score_col <- frozen_score_candidates[
    frozen_score_candidates %in% names(B12B_D)
][1L]

if (!length(frozen_score_col)) {
    frozen_score_col <- NA_character_
}

############################################################
## Repeated nested validation
############################################################

B12B_repeats <- 10L
B12B_outer_k <- 5L

B12B_fold_rows <- list()
B12B_pred_rows <- list()
B12B_audit_rows <- list()
B12B_select_rows <- list()

row_id <- 0L

for (rr in seq_len(B12B_repeats)) {

    outer_fold <- B12B_make_folds(
        cancer = B12B_D$cancer,
        event = B12B_D$pfi_event,
        k = B12B_outer_k,
        seed = 20260813 + rr
    )

    for (ff in seq_len(B12B_outer_k)) {

        row_id <- row_id + 1L

        tr_idx <- which(outer_fold != ff)
        te_idx <- which(outer_fold == ff)

        tr <- B12B_D[tr_idx]
        te <- B12B_D[te_idx]

        construction <- B12B_construct_ucei(
            train = tr,
            test = te
        )

        ######################################################
        ## Baseline matrix:
        ## training-scaled burden and age + cancer indicators
        ######################################################

        base_numeric <- B12B_train_scale(
            train = as.matrix(
                tr[, .(burden_total, age)]
            ),
            test = as.matrix(
                te[, .(burden_total, age)]
            )
        )

        cancer_levels <- sort(unique(tr$cancer))

        ftr <- factor(
            tr$cancer,
            levels = cancer_levels
        )

        fte <- factor(
            te$cancer,
            levels = cancer_levels
        )

        cancer_tr <- model.matrix(
            ~ ftr
        )[, -1L, drop = FALSE]

        cancer_te <- model.matrix(
            ~ fte
        )[, -1L, drop = FALSE]

        colnames(cancer_tr) <- paste0(
            "cancer_",
            cancer_levels[-1L]
        )
        colnames(cancer_te) <- colnames(cancer_tr)

        Xbase_tr <- cbind(
            burden_total = base_numeric$train[, "burden_total"],
            age = base_numeric$train[, "age"],
            cancer_tr
        )

        Xbase_te <- cbind(
            burden_total = base_numeric$test[, "burden_total"],
            age = base_numeric$test[, "age"],
            cancer_te
        )

        ######################################################
        ## SOTA predictor scaling learned in outer training
        ######################################################

        sota_scaled <- B12B_train_scale(
            train = as.matrix(tr[, ..B12B_SOTA]),
            test = as.matrix(te[, ..B12B_SOTA])
        )

        Xsota_tr <- sota_scaled$train
        Xsota_te <- sota_scaled$test

        Xu_tr <- cbind(
            Xbase_tr,
            UCEI = construction$train_score
        )
        Xu_te <- cbind(
            Xbase_te,
            UCEI = construction$test_score
        )

        Xs_tr <- cbind(
            Xbase_tr,
            Xsota_tr
        )
        Xs_te <- cbind(
            Xbase_te,
            Xsota_te
        )

        Xc_tr <- cbind(
            Xbase_tr,
            Xsota_tr,
            UCEI = construction$train_score
        )
        Xc_te <- cbind(
            Xbase_te,
            Xsota_te,
            UCEI = construction$test_score
        )

        ytr <- Surv(
            tr$pfi_time,
            tr$pfi_event
        )

        inner_fold <- B12B_make_folds(
            cancer = tr$cancer,
            event = tr$pfi_event,
            k = 5L,
            seed = 20260813 + rr * 100L + ff
        )

        ######################################################
        ## Baseline: unpenalized Cox through glmnet framework
        ######################################################

        fit_b <- B12B_fit_predict(
            Xtr = Xbase_tr,
            Xte = Xbase_te,
            ytr = ytr,
            penalty = rep(0, ncol(Xbase_tr)),
            inner_fold = inner_fold
        )

        fit_u <- B12B_fit_predict(
            Xtr = Xu_tr,
            Xte = Xu_te,
            ytr = ytr,
            penalty = c(
                rep(0, ncol(Xbase_tr)),
                1
            ),
            inner_fold = inner_fold
        )

        fit_s <- B12B_fit_predict(
            Xtr = Xs_tr,
            Xte = Xs_te,
            ytr = ytr,
            penalty = c(
                rep(0, ncol(Xbase_tr)),
                rep(1, length(B12B_SOTA))
            ),
            inner_fold = inner_fold
        )

        fit_c <- B12B_fit_predict(
            Xtr = Xc_tr,
            Xte = Xc_te,
            ytr = ytr,
            penalty = c(
                rep(0, ncol(Xbase_tr)),
                rep(1, length(B12B_SOTA)),
                1
            ),
            inner_fold = inner_fold
        )

        pred_list <- list(
            baseline = fit_b$pred,
            strict_UCEI = fit_u$pred,
            AllSOTA39 = fit_s$pred,
            AllSOTA39_plus_strict_UCEI = fit_c$pred
        )

        for (mn in names(pred_list)) {
            pp <- pred_list[[mn]]

            B12B_fold_rows[[length(B12B_fold_rows) + 1L]] <- data.table(
                rep_id = rr,
                fold = ff,
                model = mn,
                n = length(te_idx),
                events = sum(te$pfi_event),
                C = B12B_Cindex(
                    te$pfi_time,
                    te$pfi_event,
                    pp
                )
            )

            B12B_pred_rows[[length(B12B_pred_rows) + 1L]] <- data.table(
                patient_id = te$patient_id,
                rep_id = rr,
                fold = ff,
                model = mn,
                pfi_time = te$pfi_time,
                pfi_event = te$pfi_event,
                linear_predictor = pp
            )
        }

        load_congruence <- NA_real_
        if (!is.null(B12B_frozen_loading)) {
            common <- intersect(
                names(B12B_frozen_loading),
                names(construction$loadings)
            )

            if (length(common) == 6L) {
                x <- construction$loadings[common]
                y <- B12B_frozen_loading[common]
                load_congruence <- sum(x * y) /
                    sqrt(sum(x^2) * sum(y^2))
                load_congruence <- abs(load_congruence)
            }
        }

        strict_frozen_rho <- NA_real_
        if (!is.na(frozen_score_col)) {
            strict_frozen_rho <- suppressWarnings(
                cor(
                    construction$test_score,
                    te[[frozen_score_col]],
                    method = "spearman",
                    use = "complete.obs"
                )
            )
        }

        B12B_audit_rows[[row_id]] <- data.table(
            rep_id = rr,
            fold = ff,
            n_train = nrow(tr),
            n_test = nrow(te),
            events_train = sum(tr$pfi_event),
            events_test = sum(te$pfi_event),
            pc1_variance = construction$pc1_variance,
            loading_congruence = load_congruence,
            strict_vs_frozen_spearman = strict_frozen_rho,
            lambda_UCEI = fit_u$lambda,
            lambda_SOTA = fit_s$lambda,
            lambda_combo = fit_c$lambda
        )

        coef_u <- fit_u$coef
        coef_s <- fit_s$coef
        coef_c <- fit_c$coef

        nz_u <- rownames(coef_u)[abs(coef_u[, 1L]) > 1e-10]
        nz_s <- rownames(coef_s)[abs(coef_s[, 1L]) > 1e-10]
        nz_c <- rownames(coef_c)[abs(coef_c[, 1L]) > 1e-10]

        B12B_select_rows[[row_id]] <- data.table(
            rep_id = rr,
            fold = ff,
            UCEI_selected_alone = as.integer("UCEI" %in% nz_u),
            UCEI_selected_combo = as.integer("UCEI" %in% nz_c),
            n_SOTA_selected = sum(B12B_SOTA %in% nz_s),
            n_SOTA_selected_combo = sum(B12B_SOTA %in% nz_c)
        )
    }
}

B12B_FOLD_C <- rbindlist(B12B_fold_rows)
B12B_PRED <- rbindlist(B12B_pred_rows)
B12B_AUDIT <- rbindlist(B12B_audit_rows)
B12B_SELECT <- rbindlist(B12B_select_rows)

############################################################
## Aggregate held-out predictions within each repeat
############################################################

B12B_REPEAT_C <- B12B_PRED[
    ,
    .(
        C = B12B_Cindex(
            pfi_time,
            pfi_event,
            linear_predictor
        ),
        mean_fold_C = mean(
            B12B_FOLD_C[
                rep_id == .BY$rep_id &
                model == .BY$model,
                C
            ]
        ),
        min_fold_C = min(
            B12B_FOLD_C[
                rep_id == .BY$rep_id &
                model == .BY$model,
                C
            ]
        ),
        max_fold_C = max(
            B12B_FOLD_C[
                rep_id == .BY$rep_id &
                model == .BY$model,
                C
            ]
        )
    ),
    by = .(rep_id, model)
]

B12B_C_SUMMARY <- B12B_REPEAT_C[
    ,
    .(
        mean_C = mean(C),
        sd_C = sd(C),
        min_C = min(C),
        max_C = max(C)
    ),
    by = model
][order(-mean_C)]

B12B_WIDE <- dcast(
    B12B_REPEAT_C,
    rep_id ~ model,
    value.var = "C"
)

B12B_DELTA <- rbindlist(list(
    B12B_WIDE[, .(
        rep_id,
        comparison = "strict_UCEI - baseline",
        delta_C = strict_UCEI - baseline
    )],
    B12B_WIDE[, .(
        rep_id,
        comparison = "strict_UCEI - AllSOTA39",
        delta_C = strict_UCEI - AllSOTA39
    )],
    B12B_WIDE[, .(
        rep_id,
        comparison = "combo - AllSOTA39",
        delta_C = AllSOTA39_plus_strict_UCEI - AllSOTA39
    )],
    B12B_WIDE[, .(
        rep_id,
        comparison = "combo - strict_UCEI",
        delta_C = AllSOTA39_plus_strict_UCEI - strict_UCEI
    )]
))

B12B_DELTA_SUMMARY <- B12B_DELTA[
    ,
    .(
        mean_delta_C = mean(delta_C),
        sd_delta_C = sd(delta_C),
        min_delta_C = min(delta_C),
        max_delta_C = max(delta_C),
        positive_repeats = sum(delta_C > 0),
        total_repeats = .N
    ),
    by = comparison
]

B12B_CONSTRUCTION <- B12B_AUDIT[
    ,
    .(
        folds = .N,
        pc1_variance_mean = mean(pc1_variance),
        pc1_variance_min = min(pc1_variance),
        pc1_variance_max = max(pc1_variance),
        loading_congruence_mean = mean(
            loading_congruence,
            na.rm = TRUE
        ),
        loading_congruence_min = min(
            loading_congruence,
            na.rm = TRUE
        ),
        strict_vs_frozen_spearman_mean = mean(
            strict_vs_frozen_spearman,
            na.rm = TRUE
        ),
        strict_vs_frozen_spearman_min = min(
            strict_vs_frozen_spearman,
            na.rm = TRUE
        )
    )
]

B12B_SELECTION_SUMMARY <- B12B_SELECT[
    ,
    .(
        folds = .N,
        UCEI_selected_alone = sum(UCEI_selected_alone),
        UCEI_selected_combo = sum(UCEI_selected_combo),
        median_SOTA_selected = median(n_SOTA_selected),
        range_SOTA_selected = paste(
            range(n_SOTA_selected),
            collapse = "-"
        ),
        median_SOTA_selected_combo = median(n_SOTA_selected_combo),
        range_SOTA_selected_combo = paste(
            range(n_SOTA_selected_combo),
            collapse = "-"
        )
    )
]

############################################################
## Write historical B12B outputs
############################################################

fwrite(
    B12B_C_SUMMARY,
    file.path(output_dir, "B12B_C_summary.csv")
)

fwrite(
    B12B_AUDIT,
    file.path(output_dir, "B12B_construction_audit_50_folds.csv")
)

fwrite(
    B12B_CONSTRUCTION,
    file.path(output_dir, "B12B_construction_summary.csv")
)

fwrite(
    B12B_DELTA,
    file.path(output_dir, "B12B_delta_C_by_repeat.csv")
)

fwrite(
    B12B_DELTA_SUMMARY,
    file.path(output_dir, "B12B_delta_C_summary.csv")
)

fwrite(
    B12B_FOLD_C,
    file.path(output_dir, "B12B_fold_specific_C.csv")
)

fwrite(
    B12B_PRED,
    file.path(output_dir, "B12B_outer_test_predictions.csv")
)

fwrite(
    B12B_REPEAT_C,
    file.path(output_dir, "B12B_repeat_level_C.csv")
)

fwrite(
    B12B_SELECT,
    file.path(output_dir, "B12B_selection_50_folds.csv")
)

fwrite(
    B12B_SELECTION_SUMMARY,
    file.path(output_dir, "B12B_selection_summary.csv")
)

cat("\nB12B complete.\n")
print(B12B_C_SUMMARY)
print(B12B_DELTA_SUMMARY)
print(B12B_SELECTION_SUMMARY)
print(B12B_CONSTRUCTION)
