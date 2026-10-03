############################################################
## B12A: UCEI CONSTRUCTION STABILITY
##
## Publication-facing reconstruction from the final executed
## Stage-3/B12A analysis record.
##
## Final historical settings retained:
##   - six raw architecture features
##   - cancer-specific burden residualization re-estimated
##     within every refit
##   - residual centering/scaling re-estimated within refit
##   - PCA: prcomp(Z, center = FALSE, scale. = FALSE)
##   - 200 cancer-stratified bootstrap reconstructions
##   - leave-one-cancer-out reconstruction across cancers
##   - seed = 20260813
##
## Required input columns:
##   cancer
##   burden_total
##   jsd_neutral
##   ent_state_size
##   ent_chr_sd
##   ent_chr_mean
##   ent_ampclass
##   ent_transition
##
## Usage:
##   Rscript 02_ucei_construction_stability.R input.csv output_dir
##
## The input table is the derivation feature table before
## burden residualization. Third-party raw TCGA data are not
## redistributed in the repository.
############################################################

suppressPackageStartupMessages({
    library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

input_file <- if (length(args) >= 1L) {
    args[[1L]]
} else {
    Sys.getenv("UCEI_DERIVATION_FEATURE_TABLE", unset = "")
}

output_dir <- if (length(args) >= 2L) {
    args[[2L]]
} else {
    file.path("results", "benchmark", "B12A_construction_stability")
}

if (!nzchar(input_file) || !file.exists(input_file)) {
    stop(
        "Derivation feature table not found. Supply it as the first command-line ",
        "argument or set UCEI_DERIVATION_FEATURE_TABLE."
    )
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

B12A_features <- c(
    "jsd_neutral",
    "ent_state_size",
    "ent_chr_sd",
    "ent_chr_mean",
    "ent_ampclass",
    "ent_transition"
)

B12A_required <- c(
    "cancer",
    "burden_total",
    B12A_features
)

B12A_D <- if (grepl("\\.rds$", input_file, ignore.case = TRUE)) {
    as.data.table(readRDS(input_file))
} else {
    fread(input_file)
}

B12A_missing <- setdiff(B12A_required, names(B12A_D))
if (length(B12A_missing)) {
    stop(
        "Missing required columns: ",
        paste(B12A_missing, collapse = ", ")
    )
}

B12A_D <- B12A_D[
    complete.cases(B12A_D[, ..B12A_required])
]

B12A_D[, cancer := as.character(cancer)]
B12A_D[, burden_total := as.numeric(burden_total)]

for (f in B12A_features) {
    set(B12A_D, j = f, value = as.numeric(B12A_D[[f]]))
}

############################################################
## Core UCEI refit
############################################################

B12A_fit_ucei <- function(train) {

    train <- as.data.table(train)

    Z <- matrix(
        NA_real_,
        nrow = nrow(train),
        ncol = length(B12A_features),
        dimnames = list(NULL, B12A_features)
    )

    cancers <- unique(train$cancer)

    for (cc in cancers) {

        ii <- which(train$cancer == cc)
        x <- train$burden_total[ii]

        X <- cbind(
            intercept = 1,
            burden_total = x
        )

        for (f in B12A_features) {

            y <- train[[f]][ii]

            fit <- lm.fit(
                x = X,
                y = y
            )

            resid <- as.numeric(fit$residuals)

            resid_mean <- mean(
                resid,
                na.rm = TRUE
            )

            resid_sd <- sd(
                resid,
                na.rm = TRUE
            )

            if (!is.finite(resid_sd) || resid_sd <= 0) {
                stop(
                    "Non-finite/zero residual SD for cancer ",
                    cc,
                    ", feature ",
                    f
                )
            }

            Z[ii, f] <- (
                resid - resid_mean
            ) / resid_sd
        }
    }

    if (any(!is.finite(Z))) {
        stop("Non-finite standardized residuals encountered.")
    }

    pca <- prcomp(
        Z,
        center = FALSE,
        scale. = FALSE
    )

    loading <- pca$rotation[, 1L]
    score <- pca$x[, 1L]

    pc1_variance <- (
        pca$sdev[1L]^2 /
            sum(pca$sdev^2)
    )

    ## Historical full-reference orientation rule:
    ## orient PC1 to the direction of the mean of
    ## jsd_neutral and ent_state_size.
    orientation_target <- rowMeans(
        cbind(
            train$jsd_neutral,
            train$ent_state_size
        ),
        na.rm = TRUE
    )

    orientation_cor <- suppressWarnings(
        cor(
            score,
            orientation_target,
            use = "complete.obs"
        )
    )

    if (is.finite(orientation_cor) && orientation_cor < 0) {
        loading <- -loading
        score <- -score
    }

    list(
        loading = loading,
        score = score,
        pc1_variance = pc1_variance
    )
}

B12A_loading_congruence <- function(x, ref) {

    x <- as.numeric(x)
    ref <- as.numeric(ref)

    sum(x * ref) /
        sqrt(
            sum(x^2) *
                sum(ref^2)
        )
}

############################################################
## Full derivation fit
############################################################

B12A_full <- B12A_fit_ucei(
    B12A_D
)

B12A_frozen_loading <- B12A_full$loading

############################################################
## Cancer-stratified bootstrap
############################################################

set.seed(20260813)

B12A_nboot <- 200L
B12A_cancers <- unique(B12A_D$cancer)

B12A_BOOT_list <- vector(
    "list",
    B12A_nboot
)

for (b in seq_len(B12A_nboot)) {

    boot_idx <- unlist(
        lapply(
            B12A_cancers,
            function(cc) {
                ii <- which(
                    B12A_D$cancer == cc
                )

                sample(
                    ii,
                    length(ii),
                    replace = TRUE
                )
            }
        ),
        use.names = FALSE
    )

    fit_b <- B12A_fit_ucei(
        B12A_D[boot_idx]
    )

    ld <- fit_b$loading

    ## Align bootstrap component sign to the full-reference
    ## loading vector.
    sign_cor <- suppressWarnings(
        cor(
            ld,
            B12A_frozen_loading,
            use = "complete.obs"
        )
    )

    if (is.finite(sign_cor) && sign_cor < 0) {
        ld <- -ld
    }

    congruence <- B12A_loading_congruence(
        ld,
        B12A_frozen_loading
    )

    B12A_BOOT_list[[b]] <- data.table(
        replicate = b,
        feature = B12A_features,
        loading = as.numeric(
            ld[B12A_features]
        ),
        pc1_variance = fit_b$pc1_variance,
        loading_congruence = congruence
    )
}

B12A_BOOT <- rbindlist(
    B12A_BOOT_list,
    use.names = TRUE
)

B12A_BOOT_SUMMARY <- B12A_BOOT[
    ,
    .(
        frozen_loading = B12A_frozen_loading[feature[1L]],
        bootstrap_mean = mean(loading),
        bootstrap_sd = sd(loading),
        CI025 = as.numeric(
            quantile(
                loading,
                0.025,
                names = FALSE
            )
        ),
        CI975 = as.numeric(
            quantile(
                loading,
                0.975,
                names = FALSE
            )
        ),
        sign_consistency = mean(
            sign(loading) ==
                sign(
                    B12A_frozen_loading[
                        feature[1L]
                    ]
                )
        )
    ),
    by = feature
]

B12A_BOOT_GLOBAL <- B12A_BOOT[
    ,
    .(
        n_replicates = uniqueN(replicate),
        pc1_variance_mean = mean(
            pc1_variance
        ),
        pc1_variance_CI025 = as.numeric(
            quantile(
                pc1_variance,
                0.025,
                names = FALSE
            )
        ),
        pc1_variance_CI975 = as.numeric(
            quantile(
                pc1_variance,
                0.975,
                names = FALSE
            )
        ),
        congruence_mean = mean(
            loading_congruence
        ),
        congruence_min = min(
            loading_congruence
        ),
        congruence_CI025 = as.numeric(
            quantile(
                loading_congruence,
                0.025,
                names = FALSE
            )
        )
    )
]

############################################################
## Leave-one-cancer-out reconstruction
############################################################

B12A_LOCO_list <- vector(
    "list",
    length(B12A_cancers)
)

for (i in seq_along(B12A_cancers)) {

    omitted <- B12A_cancers[[i]]

    fit_l <- B12A_fit_ucei(
        B12A_D[
            cancer != omitted
        ]
    )

    ld <- fit_l$loading

    sign_cor <- suppressWarnings(
        cor(
            ld,
            B12A_frozen_loading,
            use = "complete.obs"
        )
    )

    if (is.finite(sign_cor) && sign_cor < 0) {
        ld <- -ld
    }

    congruence <- B12A_loading_congruence(
        ld,
        B12A_frozen_loading
    )

    B12A_LOCO_list[[i]] <- data.table(
        omitted_cancer = omitted,
        feature = B12A_features,
        loading = as.numeric(
            ld[B12A_features]
        ),
        frozen_loading = as.numeric(
            B12A_frozen_loading[
                B12A_features
            ]
        ),
        abs_loading_difference = abs(
            as.numeric(
                ld[B12A_features]
            ) -
                as.numeric(
                    B12A_frozen_loading[
                        B12A_features
                    ]
                )
        ),
        loading_congruence = congruence,
        pc1_variance = fit_l$pc1_variance
    )
}

B12A_LOCO <- rbindlist(
    B12A_LOCO_list,
    use.names = TRUE
)

B12A_LOCO_SUMMARY <- unique(
    B12A_LOCO[
        ,
        .(
            omitted_cancer,
            loading_congruence,
            pc1_variance
        )
    ]
)

B12A_LOCO_FEATURE <- B12A_LOCO[
    ,
    .(
        max_abs_loading_difference = max(
            abs_loading_difference
        ),
        mean_abs_loading_difference = mean(
            abs_loading_difference
        )
    ),
    by = feature
]

############################################################
## Export — historical B12A filenames retained
############################################################

fwrite(
    B12A_BOOT,
    file.path(
        output_dir,
        "B12A_bootstrap_all_loadings.csv"
    )
)

fwrite(
    B12A_BOOT_GLOBAL,
    file.path(
        output_dir,
        "B12A_bootstrap_global_summary.csv"
    )
)

fwrite(
    B12A_BOOT_SUMMARY,
    file.path(
        output_dir,
        "B12A_bootstrap_loading_summary.csv"
    )
)

fwrite(
    B12A_LOCO_FEATURE,
    file.path(
        output_dir,
        "B12A_leave_one_cancer_out_feature_summary.csv"
    )
)

fwrite(
    B12A_LOCO,
    file.path(
        output_dir,
        "B12A_leave_one_cancer_out_loadings.csv"
    )
)

fwrite(
    B12A_LOCO_SUMMARY,
    file.path(
        output_dir,
        "B12A_leave_one_cancer_out_summary.csv"
    )
)

cat(
    "\nB12A complete.\n",
    "Bootstrap reconstructions: ",
    B12A_nboot,
    "\nMean loading congruence: ",
    sprintf(
        "%.7f",
        B12A_BOOT_GLOBAL$congruence_mean
    ),
    "\nMinimum bootstrap congruence: ",
    sprintf(
        "%.7f",
        B12A_BOOT_GLOBAL$congruence_min
    ),
    "\nMinimum LOCO congruence: ",
    sprintf(
        "%.7f",
        min(
            B12A_LOCO_SUMMARY$loading_congruence
        )
    ),
    "\nOutput directory: ",
    normalizePath(
        output_dir,
        winslash = "/",
        mustWork = FALSE
    ),
    "\n",
    sep = ""
)
