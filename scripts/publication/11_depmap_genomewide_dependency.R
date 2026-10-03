############################################################
## 11: DEPMAP GENOME-WIDE DEPENDENCY MODEL
##
## Prepared long-format input:
## ModelID, gene, gene_effect, AM_score, lineage,
## target_expression, target_CN, global_dependency
##
## One row per model x gene.
############################################################
suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly=TRUE)
input_file <- if(length(args)>=1) args[[1]] else Sys.getenv("UCEI_DEPMAP_LONG",unset="")
output_dir <- if(length(args)>=2) args[[2]] else file.path("results","functional","depmap_genomewide")
if(!nzchar(input_file)||!file.exists(input_file)) stop("Supply prepared DepMap long table.")
dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)

D <- if(grepl("\\.rds$",input_file,ignore.case=TRUE)) as.data.table(readRDS(input_file)) else fread(input_file)
req <- c("ModelID","gene","gene_effect","AM_score","lineage","target_expression","target_CN","global_dependency")
miss <- setdiff(req,names(D))
if(length(miss)) stop("Missing columns: ",paste(miss,collapse=", "))

fit_one <- function(z){
    z <- z[complete.cases(z[,..req])]
    if(nrow(z)<80L) return(NULL)
    f <- try(lm(gene_effect ~ AM_score + target_expression + target_CN + global_dependency + factor(lineage),data=z),silent=TRUE)
    if(inherits(f,"try-error")) return(NULL)
    s <- summary(f)$coefficients
    if(!"AM_score"%in%rownames(s)) return(NULL)
    data.table(
        n=nrow(z),
        beta_AM=s["AM_score","Estimate"],
        SE=s["AM_score","Std. Error"],
        p=s["AM_score","Pr(>|t|)"]
    )
}

RES <- D[,fit_one(.SD),by=gene]
RES[,FDR:=p.adjust(p,"BH")]
RES[,supported:=FDR<0.10 & beta_AM<0]
setorder(RES,FDR,beta_AM)

SUMMARY <- RES[,.(genes_tested=.N,FDR_supported=sum(supported),nominal_directional=sum(p<0.05 & beta_AM<0))]
fwrite(RES,file.path(output_dir,"depmap_genomewide_dependency_results.csv"))
fwrite(SUMMARY,file.path(output_dir,"depmap_genomewide_dependency_summary.csv"))

LOCKED <- data.table(expected_models=852L,expected_genes=17566L,expected_FDR_supported=638L,expected_nominal_directional=753L)
fwrite(LOCKED,file.path(output_dir,"depmap_locked_targets.csv"))
cat("\nDepMap genome-wide dependency model complete.\n")
print(SUMMARY)
