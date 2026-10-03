############################################################
## 12: CN-AWARE TRANSCRIPTIONAL DECOMPOSITION
##
## Prepared long-format input:
## ModelID, gene, expression, local_CN, reduced_UCEI,
## global_burden, lineage
############################################################
suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly=TRUE)
input_file <- if(length(args)>=1) args[[1]] else Sys.getenv("UCEI_CN_AWARE_LONG",unset="")
output_dir <- if(length(args)>=2) args[[2]] else file.path("results","functional","cn_aware_transcription")
if(!nzchar(input_file)||!file.exists(input_file)) stop("Supply prepared CN-aware long table.")
dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)

D <- if(grepl("\\.rds$",input_file,ignore.case=TRUE)) as.data.table(readRDS(input_file)) else fread(input_file)
req <- c("ModelID","gene","expression","local_CN","reduced_UCEI","global_burden","lineage")
miss <- setdiff(req,names(D))
if(length(miss)) stop("Missing columns: ",paste(miss,collapse=", "))

fit_vars <- setdiff(req,"gene")

fit_one <- function(z){
    z <- z[complete.cases(z[,..fit_vars])]
    if(nrow(z)<80L) return(NULL)
    f <- try(lm(expression ~ local_CN + reduced_UCEI + global_burden + factor(lineage),data=z),silent=TRUE)
    if(inherits(f,"try-error")) return(NULL)
    s <- summary(f)$coefficients
    if(!all(c("local_CN","reduced_UCEI")%in%rownames(s))) return(NULL)
    data.table(
        n=nrow(z),
        beta_local_CN=s["local_CN","Estimate"],
        p_local_CN=s["local_CN","Pr(>|t|)"],
        beta_UCEI=s["reduced_UCEI","Estimate"],
        p_UCEI=s["reduced_UCEI","Pr(>|t|)"]
    )
}

RES <- D[,fit_one(.SD),by=gene]
RES[,FDR_local_CN:=p.adjust(p_local_CN,"BH")]
RES[,FDR_UCEI:=p.adjust(p_UCEI,"BH")]
RES[,positive_local_CN:=beta_local_CN>0 & FDR_local_CN<0.10]
RES[,positive_residual_UCEI:=beta_UCEI>0 & FDR_UCEI<0.10]

SUMMARY <- RES[,.(genes_tested=.N,
                   positive_local_CN=sum(positive_local_CN),
                   positive_residual_UCEI=sum(positive_residual_UCEI))]
fwrite(RES,file.path(output_dir,"cn_aware_gene_results.csv"))
fwrite(SUMMARY,file.path(output_dir,"cn_aware_summary.csv"))

LOCKED <- data.table(expected_genes=18393L,expected_positive_local_CN=15561L,expected_positive_residual_UCEI=0L)
fwrite(LOCKED,file.path(output_dir,"cn_aware_locked_targets.csv"))
cat("\nCN-aware transcriptional decomposition complete.\n")
print(SUMMARY)