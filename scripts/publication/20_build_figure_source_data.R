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
