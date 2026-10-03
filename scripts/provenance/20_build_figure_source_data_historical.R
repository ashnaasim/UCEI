############################################################
## ARCHIVAL HISTORICAL FIGURE-DATA HARVESTER
##
## This file is retained for provenance. It contains workstation-
## specific paths from the executed September 2026 harvest and is
## NOT the portable public entry point.
############################################################

## UCEI-v2 FINAL FIGURE DATA HARVESTER
##
## FIGURES 1–6
##
## PURPOSE:
## Collect ALL source data needed for final figures BEFORE
## any plotting is attempted.
##
## IMPORTANT:
## - NO rm()
## - NO gc()
## - NO analysis rerun
## - NO model refitting
## - NO manual result fabrication
## - NO [[ syntax
## - current benchmark workspace remains untouched
##
############################################################
############################################################

suppressPackageStartupMessages({
    library(data.table)
})

############################################################
## 0. OUTPUT DIRECTORY
############################################################

FD_OUT <- paste0(
    "E:/cnv/UCEI_BIODATAMINING_SUBMISSION/",
    "FINAL_FIGURE_DATA"
)

dir.create(
    FD_OUT,
    recursive = TRUE,
    showWarnings = FALSE
)

############################################################
## 1. MANIFEST
############################################################

FD_MANIFEST <- data.table(
    figure = character(),
    piece = character(),
    source_type = character(),
    source = character(),
    status = character(),
    class = character(),
    rows = integer(),
    cols = integer(),
    size_MB = numeric(),
    columns = character()
)

############################################################
## 2. HELPERS
############################################################

fd_describe <- function(x) {
    
    if (is.null(x)) {
        
        return(
            list(
                class = NA_character_,
                rows = NA_integer_,
                cols = NA_integer_,
                size_MB = NA_real_,
                columns = ""
            )
        )
    }
    
    cls <- paste(
        class(x),
        collapse = "/"
    )
    
    nr <- tryCatch(
        nrow(x),
        error = function(e) NA_integer_
    )
    
    nc <- tryCatch(
        ncol(x),
        error = function(e) NA_integer_
    )
    
    if (is.null(nr)) {
        nr <- NA_integer_
    }
    
    if (is.null(nc)) {
        nc <- NA_integer_
    }
    
    cn <- tryCatch(
        colnames(x),
        error = function(e) NULL
    )
    
    if (is.null(cn) && is.list(x)) {
        cn <- names(x)
    }
    
    if (is.null(cn)) {
        cn <- character()
    }
    
    list(
        class = cls,
        rows = as.integer(nr),
        cols = as.integer(nc),
        size_MB = round(
            as.numeric(
                object.size(x)
            ) / 1024^2,
            3
        ),
        columns = paste(
            head(
                cn,
                100
            ),
            collapse = " | "
        )
    )
}

fd_record <- function(
        figure,
        piece,
        source_type,
        source,
        value,
        status = NULL
) {
    
    info <- fd_describe(
        value
    )
    
    if (is.null(status)) {
        
        status <- if (
            is.null(value)
        ) {
            "MISSING"
        } else {
            "PASS"
        }
    }
    
    FD_MANIFEST <<- rbind(
        FD_MANIFEST,
        data.table(
            figure = figure,
            piece = piece,
            source_type = source_type,
            source = source,
            status = status,
            class = info$class,
            rows = info$rows,
            cols = info$cols,
            size_MB = info$size_MB,
            columns = info$columns
        ),
        fill = TRUE
    )
    
    invisible(value)
}

############################################################
## Get a current-session object.
############################################################

fd_object <- function(
        object_name,
        figure,
        required = FALSE
) {
    
    if (
        exists(
            object_name,
            envir = .GlobalEnv,
            inherits = FALSE
        )
    ) {
        
        x <- get(
            object_name,
            envir = .GlobalEnv,
            inherits = FALSE
        )
        
        fd_record(
            figure = figure,
            piece = object_name,
            source_type = "current_R_object",
            source = object_name,
            value = x,
            status = "PASS"
        )
        
        return(x)
    }
    
    fd_record(
        figure = figure,
        piece = object_name,
        source_type = "current_R_object",
        source = object_name,
        value = NULL,
        status = if (
            required
        ) {
            "MISSING_REQUIRED"
        } else {
            "NOT_PRESENT_OPTIONAL"
        }
    )
    
    NULL
}

############################################################
## Read one file.
############################################################

fd_file <- function(
        path,
        figure,
        piece,
        required = FALSE
) {
    
    if (!file.exists(path)) {
        
        fd_record(
            figure = figure,
            piece = piece,
            source_type = "disk_file",
            source = path,
            value = NULL,
            status = if (
                required
            ) {
                "MISSING_REQUIRED"
            } else {
                "NOT_PRESENT_OPTIONAL"
            }
        )
        
        return(NULL)
    }
    
    ext <- tolower(
        tools::file_ext(path)
    )
    
    x <- tryCatch(
        
        {
            
            if (ext == "rds") {
                
                readRDS(path)
                
            } else if (
                ext %in%
                c(
                    "csv",
                    "tsv",
                    "txt"
                )
            ) {
                
                suppressWarnings(
                    fread(
                        path,
                        showProgress = FALSE
                    )
                )
                
            } else {
                
                NULL
            }
        },
        
        error = function(e) {
            NULL
        }
    )
    
    fd_record(
        figure = figure,
        piece = piece,
        source_type = "disk_file",
        source = path,
        value = x,
        status = if (
            is.null(x)
        ) {
            "READ_ERROR"
        } else {
            "PASS"
        }
    )
    
    x
}

############################################################
## Add named item to a list without [[.
############################################################

fd_append_named <- function(
        x,
        name,
        value
) {
    
    x[
        length(x) + 1L
    ] <- list(
        value
    )
    
    names(x)[
        length(x)
    ] <- name
    
    x
}

############################################################
## Read all CSV tables from a known result directory.
############################################################

fd_directory_csv <- function(
        directory,
        figure,
        prefix,
        recursive = FALSE,
        max_file_MB = 75
) {
    
    out <- list()
    
    if (!dir.exists(directory)) {
        
        fd_record(
            figure = figure,
            piece = paste0(
                prefix,
                "_directory"
            ),
            source_type = "disk_directory",
            source = directory,
            value = NULL,
            status = "DIRECTORY_MISSING"
        )
        
        return(out)
    }
    
    files <- list.files(
        directory,
        pattern = "\\.csv$",
        recursive = recursive,
        full.names = TRUE,
        ignore.case = TRUE
    )
    
    if (length(files) == 0L) {
        
        fd_record(
            figure = figure,
            piece = paste0(
                prefix,
                "_directory"
            ),
            source_type = "disk_directory",
            source = directory,
            value = NULL,
            status = "NO_CSV_FILES"
        )
        
        return(out)
    }
    
    for (i in seq_along(files)) {
        
        f <- files[i]
        
        mb <- as.numeric(
            file.info(f)$size
        ) / 1024^2
        
        if (
            !is.finite(mb) ||
            mb > max_file_MB
        ) {
            
            fd_record(
                figure = figure,
                piece = basename(f),
                source_type = "disk_file",
                source = f,
                value = NULL,
                status = "SKIPPED_TOO_LARGE"
            )
            
            next
        }
        
        z <- tryCatch(
            suppressWarnings(
                fread(
                    f,
                    showProgress = FALSE
                )
            ),
            error = function(e) NULL
        )
        
        nm <- basename(f)
        
        if (!is.null(z)) {
            
            out <- fd_append_named(
                out,
                nm,
                z
            )
        }
        
        fd_record(
            figure = figure,
            piece = paste0(
                prefix,
                "::",
                nm
            ),
            source_type = "disk_file",
            source = f,
            value = z,
            status = if (
                is.null(z)
            ) {
                "READ_ERROR"
            } else {
                "PASS"
            }
        )
    }
    
    out
}

############################################################
## Recursively extract tables from report/list objects.
############################################################

fd_flatten_tables <- function(
        x,
        prefix = "report",
        depth = 0L,
        max_depth = 8L
) {
    
    out <- list()
    
    if (
        is.data.frame(x) ||
        inherits(
            x,
            "data.table"
        )
    ) {
        
        out <- fd_append_named(
            out,
            prefix,
            as.data.table(
                copy(x)
            )
        )
        
        return(out)
    }
    
    if (
        is.list(x) &&
        depth < max_depth
    ) {
        
        nms <- names(x)
        
        if (is.null(nms)) {
            
            nms <- paste0(
                "element_",
                seq_along(x)
            )
        }
        
        for (i in seq_along(x)) {
            
            nm <- nms[i]
            
            if (
                is.na(nm) ||
                nm == ""
            ) {
                
                nm <- paste0(
                    "element_",
                    i
                )
            }
            
            child <- getElement(
                x,
                i
            )
            
            tmp <- fd_flatten_tables(
                child,
                prefix = paste0(
                    prefix,
                    "$",
                    nm
                ),
                depth = depth + 1L,
                max_depth = max_depth
            )
            
            if (length(tmp) > 0L) {
                out <- c(
                    out,
                    tmp
                )
            }
        }
    }
    
    out
}

############################################################
## Recursively extract short scalar/vector metadata.
############################################################

fd_flatten_metadata <- function(
        x,
        prefix = "report",
        depth = 0L,
        max_depth = 8L
) {
    
    out <- list()
    
    if (
        is.atomic(x) &&
        length(x) > 0L &&
        length(x) <= 100L
    ) {
        
        out <- fd_append_named(
            out,
            prefix,
            x
        )
        
        return(out)
    }
    
    if (
        is.list(x) &&
        !is.data.frame(x) &&
        depth < max_depth
    ) {
        
        nms <- names(x)
        
        if (is.null(nms)) {
            
            nms <- paste0(
                "element_",
                seq_along(x)
            )
        }
        
        for (i in seq_along(x)) {
            
            nm <- nms[i]
            
            if (
                is.na(nm) ||
                nm == ""
            ) {
                
                nm <- paste0(
                    "element_",
                    i
                )
            }
            
            child <- getElement(
                x,
                i
            )
            
            tmp <- fd_flatten_metadata(
                child,
                prefix = paste0(
                    prefix,
                    "$",
                    nm
                ),
                depth = depth + 1L,
                max_depth = max_depth
            )
            
            if (length(tmp) > 0L) {
                out <- c(
                    out,
                    tmp
                )
            }
        }
    }
    
    out
}

############################################################
## Collect multiple current objects.
############################################################

fd_objects <- function(
        object_names,
        figure
) {
    
    out <- list()
    
    for (nm in object_names) {
        
        z <- fd_object(
            nm,
            figure = figure,
            required = FALSE
        )
        
        if (!is.null(z)) {
            
            out <- fd_append_named(
                out,
                nm,
                z
            )
        }
    }
    
    out
}

############################################################
## Resolve exact basename anywhere below a root directory.
## Used only for one final 15B file whose timestamped folder
## can vary.
############################################################

fd_find_basename <- function(
        root,
        target
) {
    
    if (!dir.exists(root)) {
        return(
            character()
        )
    }
    
    hits <- list.files(
        root,
        recursive = TRUE,
        full.names = TRUE
    )
    
    hits[
        basename(hits) ==
            target
    ]
}

############################################################
## 3. CREATE MASTER FIGURE DATA LIST
############################################################

FIGDATA <- list(
    Figure1 = list(),
    Figure2 = list(),
    Figure3 = list(),
    Figure4 = list(),
    Figure5 = list(),
    Figure6 = list()
)

############################################################
############################################################
##
## FIGURE 1
## UCEI-v2 CONSTRUCTION + STABILITY
##
############################################################
############################################################

FD_REF <- paste0(
    "E:/cnv/UCEI_BENCHMARK_RESCUE/",
    "UCEI_V2_FROZEN_REFERENCE"
)

FD_B12A <- paste0(
    "E:/cnv/UCEI_BENCHMARK_RESCUE/",
    "B12A_UCEI_CONSTRUCTION_STABILITY"
)

FD_B7 <- paste0(
    "E:/cnv/UCEI_BENCHMARK_RESCUE/",
    "B2_CORE_UCEI_V2_BENCHMARK/",
    "B7_RUNTIME_REPRODUCIBILITY"
)

FIGDATA$Figure1$reference_loadings <- fd_file(
    file.path(
        FD_REF,
        "UCEI_v2_frozen_PCA_loadings.csv"
    ),
    "Figure1",
    "reference_loadings",
    TRUE
)

FIGDATA$Figure1$pca_variance <- fd_file(
    file.path(
        FD_REF,
        "UCEI_v2_PCA_variance_explained.csv"
    ),
    "Figure1",
    "pca_variance",
    TRUE
)

FIGDATA$Figure1$reference_scores <- fd_file(
    file.path(
        FD_REF,
        "UCEI_v2_reference_scores.csv"
    ),
    "Figure1",
    "reference_scores",
    TRUE
)

FIGDATA$Figure1$score_parameters <- fd_file(
    file.path(
        FD_REF,
        "UCEI_v2_frozen_score_parameters.csv"
    ),
    "Figure1",
    "score_parameters",
    TRUE
)

FIGDATA$Figure1$residualization_parameters <- fd_file(
    file.path(
        FD_REF,
        "UCEI_v2_frozen_burden_residualization_parameters.csv"
    ),
    "Figure1",
    "burden_residualization_parameters",
    TRUE
)

FIGDATA$Figure1$self_reproduction <- fd_file(
    file.path(
        FD_REF,
        "UCEI_v2_self_reproduction_summary.csv"
    ),
    "Figure1",
    "self_reproduction_summary",
    FALSE
)

FIGDATA$Figure1$bootstrap_all <- fd_file(
    file.path(
        FD_B12A,
        "B12A_bootstrap_all_loadings.csv"
    ),
    "Figure1",
    "bootstrap_all_loadings",
    TRUE
)

FIGDATA$Figure1$bootstrap_global <- fd_file(
    file.path(
        FD_B12A,
        "B12A_bootstrap_global_summary.csv"
    ),
    "Figure1",
    "bootstrap_global_summary",
    TRUE
)

FIGDATA$Figure1$bootstrap_summary <- fd_file(
    file.path(
        FD_B12A,
        "B12A_bootstrap_loading_summary.csv"
    ),
    "Figure1",
    "bootstrap_loading_summary",
    TRUE
)

FIGDATA$Figure1$loco_summary <- fd_file(
    file.path(
        FD_B12A,
        "B12A_leave_one_cancer_out_summary.csv"
    ),
    "Figure1",
    "construction_LOCO_summary",
    TRUE
)

FIGDATA$Figure1$loco_loadings <- fd_file(
    file.path(
        FD_B12A,
        "B12A_leave_one_cancer_out_loadings.csv"
    ),
    "Figure1",
    "construction_LOCO_loadings",
    FALSE
)

FIGDATA$Figure1$loco_feature_summary <- fd_file(
    file.path(
        FD_B12A,
        "B12A_leave_one_cancer_out_feature_summary.csv"
    ),
    "Figure1",
    "construction_LOCO_feature_summary",
    FALSE
)

FIGDATA$Figure1$standalone_reproduction <- fd_file(
    file.path(
        FD_B7,
        "47_UCEI_v2_test_reproduction_summary.csv"
    ),
    "Figure1",
    "standalone_scorer_reproduction",
    FALSE
)

############################################################
############################################################
##
## FIGURE 2
## INFORMATION BEYOND CONVENTIONAL CNA SUMMARIES
##
############################################################
############################################################

FD_B2 <- paste0(
    "E:/cnv/UCEI_BENCHMARK_RESCUE/",
    "B2_CORE_UCEI_V2_BENCHMARK"
)

FIGDATA$Figure2$clinical_common <- fd_object(
    "B2_COMMON",
    "Figure2",
    TRUE
)

FIGDATA$Figure2$size_audit <- fd_object(
    "B2_size_audit",
    "Figure2",
    FALSE
)

FD_B2_REPORT <- fd_object(
    "B2_REPORT",
    "Figure2",
    TRUE
)

FIGDATA$Figure2$report_tables <- fd_flatten_tables(
    FD_B2_REPORT,
    prefix = "B2_REPORT"
)

FIGDATA$Figure2$report_metadata <- fd_flatten_metadata(
    FD_B2_REPORT,
    prefix = "B2_REPORT"
)

############################################################
## Audit every table recovered from B2_REPORT.
############################################################

if (
    length(
        FIGDATA$Figure2$report_tables
    ) > 0L
) {
    
    for (
        nm in names(
            FIGDATA$Figure2$report_tables
        )
    ) {
        
        fd_record(
            figure = "Figure2",
            piece = nm,
            source_type = "B2_REPORT_nested_table",
            source = nm,
            value = getElement(
                FIGDATA$Figure2$report_tables,
                nm
            ),
            status = "PASS"
        )
    }
}

############################################################
## Original root-level B2 exported tables.
############################################################

FIGDATA$Figure2$disk_tables <- fd_directory_csv(
    FD_B2,
    figure = "Figure2",
    prefix = "B2_disk",
    recursive = FALSE,
    max_file_MB = 25
)

############################################################
## Reference cohort provenance comes from exact v2 scores.
############################################################

FIGDATA$Figure2$reference_scores <- FIGDATA$Figure1$reference_scores

############################################################
############################################################
##
## FIGURE 3
## SOTA BENCHMARKING + NESTED VALIDATION
##
############################################################
############################################################

FIGDATA$Figure3$current_objects <- fd_objects(
    c(
        "B3_COMMON_DITHER",
        
        "B4_block_lrt",
        "B4_individual",
        "B4_joint",
        "B4_sig_cor",
        "B4_steele_redundancy",
        
        "B5_COMMON_DREWS",
        
        "B9_audit",
        "B9_all_sota",
        "B9_coef",
        "B9_full_coef",
        "B9_key_coef",
        "B9_LRT",
        "B9_red_summary",
        
        "B9C_cindex",
        "B9C_delta",
        "B9C_delta_summary",
        "B9C_fold_results",
        "B9C_repeat_results",
        "B9C_selection",
        "B9C_selection_frequency",
        "B9C_selection_results",
        "B9C_summary",
        "B9C_totals",
"B9C_ucei_selection"
    ),
    "Figure3"
)

############################################################
## Capture B9/B9C report tables too.
############################################################

FD_B9_REPORT <- fd_object(
    "B9_report",
    "Figure3",
    FALSE
)

FD_B9C_REPORT <- fd_object(
    "B9C_report",
    "Figure3",
    FALSE
)

FIGDATA$Figure3$B9_report_tables <- fd_flatten_tables(
    FD_B9_REPORT,
    "B9_report"
)

FIGDATA$Figure3$B9C_report_tables <- fd_flatten_tables(
    FD_B9C_REPORT,
    "B9C_report"
)

############################################################
############################################################
##
## FIGURE 4
## STRICT TRAINING-ONLY VALIDATION + ABLATION +
## TRANSPORTABILITY + EXTERNAL SUPPORT
##
############################################################
############################################################

FIGDATA$Figure4$current_objects <- fd_objects(
    c(
        "B10_cindex",
        "B10_meta",
        "B10_results",
        "B10_summary",
        
        "B12B_AUDIT",
        "B12B_C_SUMMARY",
        "B12B_CONSTRUCTION",
        "B12B_DELTA",
        "B12B_DELTA_SUMMARY",
        "B12B_FOLD_C",
        "B12B_REPEAT_C",
        "B12B_SELECT",
        "B12B_SELECTION_SUMMARY",
        
        "B12C_BETA",
        "B12C_BETA_SUMMARY",
        "B12C_CORR",
        "B12C_CORR_SUMMARY",
        "B12C_DELTA",
        "B12C_DELTA_SUMMARY",
        "B12C_FOLD",
        "B12C_SUMMARY",
        
        "B12C2_DELTA",
        "B12C2_DELTA_SUMMARY",
        "B12C2_SELECT",
        "B12C2_SELECT_SUMMARY",
        "B12C2_SUMMARY",
        
        "B12D_DELTA",
        "B12D_DELTA_SUMMARY",
        "B12D_LONG",
        "B12D_PRIMARY",
        "B12D_SUMMARY",
        "B12D_SUPPORT",
        "B12D_SUPPORT_SUMMARY"
    ),
    "Figure4"
)

############################################################
## External reduced-UCEI datasets/results
############################################################

FD_09A <- paste0(
    "E:/cnv/UCEI_DEEP_RESULTS/",
    "09A_LOCAL_EXTERNAL_VALIDATION_reduced_UCEI"
)

FD_09B <- paste0(
    "E:/cnv/UCEI_DEEP_RESULTS/",
    "09B_STRICT_EXTERNAL_VALIDATION_SUMMARY"
)

FIGDATA$Figure4$external_09A <- fd_directory_csv(
    FD_09A,
    figure = "Figure4",
    prefix = "09A_external",
    recursive = FALSE,
    max_file_MB = 50
)

FIGDATA$Figure4$external_09B <- fd_directory_csv(
    FD_09B,
    figure = "Figure4",
    prefix = "09B_strict_external",
    recursive = FALSE,
    max_file_MB = 50
)

############################################################
############################################################
##
## FIGURE 5
## UCEI–ArchSig STATE + PROGRESSION
##
############################################################
############################################################

FIGDATA$Figure5$current_objects <- fd_objects(
    c(
        "B8_archsig_by_cancer",
        "B8_archsig_overall",
        "B8_cancer_summary",
        "B8_COMMON_V2_BIOLOGY",
        "B8_cont_lrt",
        "B8_joint",
        "B8_joint_ci",
        "B8_joint_coef",
        "B8_module_coupling",
        "B8_module_coverage",
        
        "B8K_clin0",
        "B8K_clin_a",
        "B8K_clin_u",
        "B8K_clin_joint",
        "B8K_clin_state",
        "B8K_coef",
        "B8K_cohort",
        "B8K_joint_compare",
        "B8K_LRT",
        "B8K_pre0",
        "B8K_pre_a",
        "B8K_pre_u",
        "B8K_pre_joint",
        "B8K_state_cox"
    ),
    "Figure5"
)

FD_B8K_REPORT <- fd_object(
    "B8K_report",
    "Figure5",
    FALSE
)

FIGDATA$Figure5$report_tables <- fd_flatten_tables(
    FD_B8K_REPORT,
    "B8K_report"
)

FIGDATA$Figure5$report_metadata <- fd_flatten_metadata(
    FD_B8K_REPORT,
    "B8K_report"
)

############################################################
############################################################
##
## FIGURE 6
## SET 3 BIOLOGY
##
############################################################
############################################################

FD_08A <- paste0(
    "E:/cnv/UCEI_DEEP_RESULTS/",
    "08A_UCEI_survival_mediation"
)

FD_08B <- paste0(
    "E:/cnv/UCEI_DEEP_RESULTS/",
    "08B_UCEI_incremental_coupling"
)

FD_08D <- paste0(
    "E:/cnv/UCEI_DEEP_RESULTS/",
    "08D_STRICT_latent_order_manifest_drivers"
)

FD_08E <- paste0(
    "E:/cnv/UCEI_DEEP_RESULTS/",
    "08E_FORMAL_latent_order_manifest_model"
)

FD_08F2 <- paste0(
    "E:/cnv/UCEI_DEEP_RESULTS/",
    "08F2_MODULE_driver_axis_bridge_model"
)

FD_08G <- paste0(
    "E:/cnv/UCEI_DEEP_RESULTS/",
    "08G_ANTI_CIRCULARITY_leave_axis_out_specificity"
)

FD_08H <- paste0(
    "E:/cnv/UCEI_DEEP_RESULTS/",
    "08H_FINAL_manuscript_evidence_package"
)

FD_10B <- paste0(
    "E:/cnv/UCEI_DEEP_RESULTS/",
    "10B_DEPMAP_GENOMEWIDE_DEPENDENCY_SIGNATURE"
)

FD_22A <- paste0(
    "E:/cnv/UCEI_DEEP_RESULTS/",
    "22A_GENOMEWIDE_PROTEOMIC_ENTROPY_BUFFER_MAP_20260731_051236"
)

FD_22B <- paste0(
    "E:/cnv/UCEI_DEEP_RESULTS/",
    "22B_EXTERNAL_PATHWAY_ENRICHMENT_FINAL_20260902"
)

############################################################
## 6A — mediation
############################################################

FIGDATA$Figure6$mediation_08A <- fd_directory_csv(
    FD_08A,
    "Figure6",
    "08A_mediation",
    FALSE,
    75
)

############################################################
## 6B — incremental architecture/manifestation
############################################################

FIGDATA$Figure6$incremental_08B <- fd_directory_csv(
    FD_08B,
    "Figure6",
    "08B_incremental",
    FALSE,
    75
)

############################################################
## 6C — strict driver atlas
############################################################

FIGDATA$Figure6$drivers_08D <- fd_directory_csv(
    FD_08D,
    "Figure6",
    "08D_drivers",
    FALSE,
    75
)

############################################################
## 6D — formal latent-order / manifest-mechanism model
############################################################

FIGDATA$Figure6$latent_manifest_08E <- fd_directory_csv(
    FD_08E,
    "Figure6",
    "08E_latent_manifest",
    FALSE,
    75
)

############################################################
## 6E — driver-axis bridge
############################################################

FIGDATA$Figure6$axis_bridge_08F2 <- fd_directory_csv(
    FD_08F2,
    "Figure6",
    "08F2_axis_bridge",
    FALSE,
    75
)

############################################################
## 6F — anti-circularity
############################################################

FIGDATA$Figure6$anti_circularity_08G <- fd_directory_csv(
    FD_08G,
    "Figure6",
    "08G_anti_circularity",
    FALSE,
    75
)

############################################################
## 6G — existing final Set 3 evidence map
############################################################

FIGDATA$Figure6$evidence_package_08H <- fd_directory_csv(
    FD_08H,
    "Figure6",
    "08H_evidence_package",
    FALSE,
    75
)

############################################################
## 6H — DepMap genome-wide dependency signature
############################################################

FIGDATA$Figure6$depmap_report <- fd_file(
    file.path(
        FD_10B,
        "report10B_DepMap_genomewide_dependency_signature.rds"
    ),
    "Figure6",
    "10B_DepMap_report",
    TRUE
)

FIGDATA$Figure6$depmap_tables <- fd_directory_csv(
    FD_10B,
    "Figure6",
    "10B_DepMap",
    FALSE,
    75
)

############################################################
## 6I — 15B triadic entropy mechanics
##
## Folder name varied during repair, so find exact FINAL
## filename rather than guessing timestamp.
############################################################

FD_15B_HITS <- fd_find_basename(
    "E:/cnv/UCEI_DEEP_RESULTS",
    "15B_triadic_entropy_mechanics_pathway_results_FIXED.csv"
)

if (length(FD_15B_HITS) > 0L) {
    
    ########################################################
    ## Prefer latest modified copy if duplicates exist.
    ########################################################
    
    FD_15B_INFO <- file.info(
        FD_15B_HITS
    )
    
    FD_15B_ORDER <- order(
        FD_15B_INFO$mtime,
        decreasing = TRUE
    )
    
    FD_15B_FILE <- FD_15B_HITS[
        FD_15B_ORDER[1]
    ]
    
    FIGDATA$Figure6$triadic_15B <- fd_file(
        FD_15B_FILE,
        "Figure6",
        "15B_triadic_entropy_mechanics",
        FALSE
    )
    
} else {
    
    FIGDATA$Figure6$triadic_15B <- NULL
    
    fd_record(
        "Figure6",
        "15B_triadic_entropy_mechanics",
        "disk_file_search",
        "15B_triadic_entropy_mechanics_pathway_results_FIXED.csv",
        NULL,
        "NOT_FOUND_OPTIONAL"
    )
}

############################################################
## 6J — 22A genome-wide proteomic entropy map
############################################################

FIGDATA$Figure6$proteomic_22A <- fd_object(
    "G22A_gene_protein_entropy_map",
    "Figure6",
    TRUE
)

FIGDATA$Figure6$proteomic_22A_counts <- fd_objects(
    c(
        "G22A_class_audit",
        "G22A_expected_counts",
        "G22A_observed_counts",
        "G22A_n_genes",
        "G22A_n_rows",
        "G22A_prioritized_n"
    ),
    "Figure6"
)

############################################################
## Also collect exported 22A tables.
############################################################

FIGDATA$Figure6$proteomic_22A_disk <- fd_directory_csv(
    FD_22A,
    "Figure6",
    "22A_proteomic",
    FALSE,
    75
)

############################################################
## 6K — final 22B pathway enrichment
############################################################

FD_G22B_OBJECT <- fd_object(
    "G22B_final_bundle",
    "Figure6",
    FALSE
)

if (!is.null(FD_G22B_OBJECT)) {
    
    FIGDATA$Figure6$pathways_22B <- FD_G22B_OBJECT
    
} else {
    
    FIGDATA$Figure6$pathways_22B <- fd_file(
        file.path(
            FD_22B,
            "22B_FINAL_complete_bundle.rds"
        ),
        "Figure6",
        "22B_final_bundle",
        TRUE
    )
}

FIGDATA$Figure6$pathway_final_objects <- fd_objects(
    c(
        "G22B_class_table",
        "G22B_cross_class",
        "G22B_decision",
        "G22B_pathway_enrichment_final",
        "G22B_results",
        "G22B_sig_counts",
        "G22B_significant",
        "G22B_theme_summary",
        "G22B_top_combined",
        "G22B_top_nonredundant",
        "G22B_top_subclasses"
    ),
    "Figure6"
)

############################################################
## Exported 22B CSVs as additional provenance.
############################################################

FIGDATA$Figure6$pathway_22B_disk <- fd_directory_csv(
    FD_22B,
    "Figure6",
    "22B_pathway",
    FALSE,
    75
)

############################################################
############################################################
##
## 4. CRITICAL DATA AUDIT
##
############################################################
############################################################

fd_nrow <- function(x) {
    
    if (is.null(x)) {
        return(
            NA_integer_
        )
    }
    
    z <- tryCatch(
        nrow(x),
        error = function(e) NA_integer_
    )
    
    if (is.null(z)) {
        z <- NA_integer_
    }
    
    as.integer(z)
}

############################################################
## Figure 1 exact reference audit
############################################################

FD_F1_N <- fd_nrow(
    FIGDATA$Figure1$reference_scores
)

FD_F1_BOOT_N <- fd_nrow(
    FIGDATA$Figure1$bootstrap_all
)

FD_F1_BOOT_SUM_N <- fd_nrow(
    FIGDATA$Figure1$bootstrap_summary
)

FD_F1_LOCO_N <- fd_nrow(
    FIGDATA$Figure1$loco_summary
)

############################################################
## Figure 2
############################################################

FD_F2_CLIN_N <- fd_nrow(
    FIGDATA$Figure2$clinical_common
)

FD_F2_REPORT_TABLES_N <- length(
    FIGDATA$Figure2$report_tables
)

############################################################
## Figure 3
############################################################

FD_F3_NAMES <- names(
    FIGDATA$Figure3$current_objects
)

FD_F3_HAS_B9C <- "B9C_summary" %in%
    FD_F3_NAMES

FD_F3_HAS_B9 <- "B9_LRT" %in%
    FD_F3_NAMES

############################################################
## Figure 4
############################################################

FD_F4_NAMES <- names(
    FIGDATA$Figure4$current_objects
)

FD_F4_HAS_B10 <- "B10_results" %in%
    FD_F4_NAMES

FD_F4_HAS_B12B <- "B12B_C_SUMMARY" %in%
    FD_F4_NAMES

FD_F4_HAS_B12C <- "B12C_SUMMARY" %in%
    FD_F4_NAMES

FD_F4_HAS_B12D <- "B12D_SUMMARY" %in%
    FD_F4_NAMES

############################################################
## Figure 5
############################################################

FD_F5_NAMES <- names(
    FIGDATA$Figure5$current_objects
)

FD_F5_HAS_COHORT <- "B8K_cohort" %in%
    FD_F5_NAMES

FD_F5_HAS_LRT <- "B8K_LRT" %in%
    FD_F5_NAMES

############################################################
## Figure 6
############################################################

FD_F6_22A_N <- fd_nrow(
    FIGDATA$Figure6$proteomic_22A
)

FD_F6_HAS_DEPMAP <- !is.null(
    FIGDATA$Figure6$depmap_report
)

FD_F6_HAS_22B <- !is.null(
    FIGDATA$Figure6$pathways_22B
)

FD_F6_HAS_08E <- length(
    FIGDATA$Figure6$latent_manifest_08E
) > 0L

FD_F6_HAS_08G <- length(
    FIGDATA$Figure6$anti_circularity_08G
) > 0L

############################################################
## Final check table
############################################################

FD_CHECKS <- data.table(
    figure = c(
        "Figure1",
        "Figure1",
        "Figure1",
        "Figure1",
        
        "Figure2",
        "Figure2",
        
        "Figure3",
        "Figure3",
        
        "Figure4",
        "Figure4",
        "Figure4",
        "Figure4",
        
        "Figure5",
        "Figure5",
        
        "Figure6",
        "Figure6",
        "Figure6",
        "Figure6",
        "Figure6"
    ),
    
    check = c(
        "Reference cohort has 10,495 rows",
        "Bootstrap long table has 1,200 rows",
        "Bootstrap summary has 6 components",
        "Construction LOCO has 32 cancers",
        
        "Clinical benchmark cohort has 3,243 rows",
        "B2_REPORT yielded result tables",
        
        "B9 combined benchmark available",
        "B9C nested-validation summary available",
        
        "B10 cross-histology results available",
        "B12B strict-training results available",
        "B12C ablation results available",
        "B12D censoring-aware results available",
        
        "B8K state cohort available",
        "B8K likelihood-ratio results available",
        
        "08E formal latent-manifest files available",
        "08G anti-circularity files available",
        "10B DepMap report available",
        "22A gene map has 7,408 genes",
        "22B final pathway bundle available"
    ),
    
    observed = c(
        as.character(
            FD_F1_N
        ),
        as.character(
            FD_F1_BOOT_N
        ),
        as.character(
            FD_F1_BOOT_SUM_N
        ),
        as.character(
            FD_F1_LOCO_N
        ),
        
        as.character(
            FD_F2_CLIN_N
        ),
        as.character(
            FD_F2_REPORT_TABLES_N
        ),
        
        as.character(
            FD_F3_HAS_B9
        ),
        as.character(
            FD_F3_HAS_B9C
        ),
        
        as.character(
            FD_F4_HAS_B10
        ),
        as.character(
            FD_F4_HAS_B12B
        ),
        as.character(
            FD_F4_HAS_B12C
        ),
        as.character(
            FD_F4_HAS_B12D
        ),
        
        as.character(
            FD_F5_HAS_COHORT
        ),
        as.character(
            FD_F5_HAS_LRT
        ),
        
        as.character(
            FD_F6_HAS_08E
        ),
        as.character(
            FD_F6_HAS_08G
        ),
        as.character(
            FD_F6_HAS_DEPMAP
        ),
        as.character(
            FD_F6_22A_N
        ),
        as.character(
            FD_F6_HAS_22B
        )
    ),
    
    pass = c(
        identical(
            FD_F1_N,
            10495L
        ),
        identical(
            FD_F1_BOOT_N,
            1200L
        ),
        identical(
            FD_F1_BOOT_SUM_N,
            6L
        ),
        identical(
            FD_F1_LOCO_N,
            32L
        ),
        
        identical(
            FD_F2_CLIN_N,
            3243L
        ),
        FD_F2_REPORT_TABLES_N > 0L,
        
        FD_F3_HAS_B9,
        FD_F3_HAS_B9C,
        
        FD_F4_HAS_B10,
        FD_F4_HAS_B12B,
        FD_F4_HAS_B12C,
        FD_F4_HAS_B12D,
        
        FD_F5_HAS_COHORT,
        FD_F5_HAS_LRT,
        
        FD_F6_HAS_08E,
        FD_F6_HAS_08G,
        FD_F6_HAS_DEPMAP,
        identical(
            FD_F6_22A_N,
            7408L
        ),
        FD_F6_HAS_22B
    )
)

############################################################
############################################################
##
## 5. SAVE THE AUTHORITATIVE FIGURE-DATA BUNDLE
##
############################################################
############################################################

FD_BUNDLE_FILE <- file.path(
    FD_OUT,
    "UCEI_FIGURES_1_TO_6_SOURCE_DATA_BUNDLE.rds"
)

FD_MANIFEST_FILE <- file.path(
    FD_OUT,
    "UCEI_FIGURES_1_TO_6_SOURCE_MANIFEST.csv"
)

FD_CHECK_FILE <- file.path(
    FD_OUT,
    "UCEI_FIGURES_1_TO_6_CRITICAL_AUDIT.csv"
)

saveRDS(
    FIGDATA,
    FD_BUNDLE_FILE,
    compress = TRUE
)

fwrite(
    FD_MANIFEST,
    FD_MANIFEST_FILE
)

fwrite(
    FD_CHECKS,
    FD_CHECK_FILE
)

############################################################
## Also save one compact navigation object so future plotting
## never searches the filesystem again.
############################################################

FD_INDEX <- list(
    bundle_file = FD_BUNDLE_FILE,
    manifest_file = FD_MANIFEST_FILE,
    audit_file = FD_CHECK_FILE,
    figure1 = names(
        FIGDATA$Figure1
    ),
    figure2 = names(
        FIGDATA$Figure2
    ),
    figure3 = names(
        FIGDATA$Figure3
    ),
    figure4 = names(
        FIGDATA$Figure4
    ),
    figure5 = names(
        FIGDATA$Figure5
    ),
    figure6 = names(
        FIGDATA$Figure6
    )
)

saveRDS(
    FD_INDEX,
    file.path(
        FD_OUT,
        "UCEI_FIGURE_DATA_INDEX.rds"
    ),
    compress = TRUE
)

############################################################
############################################################
##
## 6. FINAL REPORT
##
############################################################
############################################################

cat(
    "\n\n============================================================\n"
)

cat(
    "UCEI FIGURE-DATA HARVEST COMPLETE\n"
)

cat(
    "============================================================\n"
)

print(
    FD_CHECKS
)

cat(
    "\n------------------------------------------------------------\n"
)

cat(
    "PASS:",
    sum(
        FD_CHECKS$pass
    ),
    "/",
    nrow(
        FD_CHECKS
    ),
    "\n"
)

cat(
    "FAILED CHECKS:",
    sum(
        !FD_CHECKS$pass
    ),
    "\n"
)

if (
    any(
        !FD_CHECKS$pass
    )
) {
    
    cat(
        "\nCHECKS REQUIRING ATTENTION:\n"
    )
    
    print(
        FD_CHECKS[
            pass == FALSE
        ]
    )
}

cat(
    "\n------------------------------------------------------------\n"
)

cat(
    "Bundle:\n",
    FD_BUNDLE_FILE,
    "\n",
    sep = ""
)

cat(
    "\nManifest:\n",
    FD_MANIFEST_FILE,
    "\n",
    sep = ""
)

cat(
    "\nCritical audit:\n",
    FD_CHECK_FILE,
    "\n",
    sep = ""
)

cat(
    "\n============================================================\n"
)

cat(
    "NO ANALYSES WERE REFITTED.\n"
)

cat(
    "NO WORKSPACE OBJECTS WERE DELETED.\n"
)