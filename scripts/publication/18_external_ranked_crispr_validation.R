############################################################
## 18: EXTERNAL RANKED CRISPR VALIDATION
##
## Inputs:
## 1. target_panel.csv: gene, mechanism_module
## 2. ranked_hits.csv: external_dataset, gene, external_rank
##
## Smaller rank = stronger external hit.
############################################################

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly=TRUE)
if(length(args)<2L) stop("Usage: 18_external_ranked_crispr_validation.R target_panel.csv ranked_hits.csv [output_dir]")
TARGET <- fread(args[[1L]])
EXT <- fread(args[[2L]])
output_dir <- if(length(args)>=3L) args[[3L]] else
    file.path("results","functional","external_ranked_crispr")
dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)

stopifnot(all(c("gene","mechanism_module")%in%names(TARGET)))
stopifnot(all(c("external_dataset","gene","external_rank")%in%names(EXT)))

clean_gene <- function(x) toupper(trimws(as.character(x)))
TARGET[,gene:=clean_gene(gene)]
EXT[,gene:=clean_gene(gene)]
EXT[,external_rank:=as.numeric(external_rank)]
EXT <- EXT[is.finite(external_rank) & nzchar(gene)]
EXT <- EXT[order(external_dataset,external_rank)]
EXT <- EXT[, .SD[1L], by=.(external_dataset,gene)]
EXT[, rank_percentile:=frank(external_rank,ties.method="average")/.N,by=external_dataset]

target_genes <- unique(TARGET$gene)
module_sets <- split(TARGET$gene,TARGET$mechanism_module)
datasets_use <- unique(EXT$external_dataset)

rank_enrich_one <- function(ds,geneset_name,geneset) {
    x <- EXT[external_dataset==ds]
    geneset <- intersect(unique(geneset),x$gene)
    if(length(geneset)==0L) {
        return(data.table(
            external_dataset=ds,geneset_name=geneset_name,
            geneset_n_in_external=0L,
            median_rank_percentile=NA_real_,
            background_median_rank_percentile=median(x$rank_percentile),
            wilcox_p=NA_real_,overlap_genes=NA_character_
        ))
    }
    x[,in_set:=gene%in%geneset]
    wt <- wilcox.test(
        x[in_set==TRUE,rank_percentile],
        x[in_set==FALSE,rank_percentile],
        alternative="less",
        exact=FALSE
    )
    data.table(
        external_dataset=ds,
        geneset_name=geneset_name,
        geneset_n_in_external=length(geneset),
        median_rank_percentile=median(x[in_set==TRUE,rank_percentile],na.rm=TRUE),
        background_median_rank_percentile=median(x[in_set==FALSE,rank_percentile],na.rm=TRUE),
        rank_enrichment_direction="lower rank percentile = stronger external hit",
        wilcox_p=wt$p.value,
        overlap_genes=paste(x[in_set==TRUE][order(external_rank),gene],collapse="; ")
    )
}

PANEL <- rbindlist(lapply(datasets_use,function(ds)
    rank_enrich_one(ds,"Full entropy-bottleneck panel",target_genes)),fill=TRUE)

MODULE <- rbindlist(lapply(datasets_use,function(ds)
    rbindlist(lapply(names(module_sets),function(m)
        rank_enrich_one(ds,m,module_sets[[m]])),fill=TRUE)),fill=TRUE)

PANEL[,FDR:=p.adjust(wilcox_p,"BH")]
MODULE[,FDR:=p.adjust(wilcox_p,"BH")]
setorder(PANEL,FDR,median_rank_percentile)
setorder(MODULE,FDR,median_rank_percentile)

## Direct top-N overlap used as the complementary 18C layer.
topNs <- c(100L,500L)
TOPN <- rbindlist(lapply(datasets_use,function(ds){
    x <- EXT[external_dataset==ds]
    rbindlist(lapply(topNs,function(N){
        hits <- head(x[order(external_rank)],N)$gene
        ov <- intersect(target_genes,hits)
        M <- uniqueN(x$gene); K <- min(N,M); n <- length(intersect(target_genes,x$gene)); k <- length(ov)
        p <- if(n>0) phyper(k-1,K,M-K,n,lower.tail=FALSE) else NA_real_
        data.table(external_dataset=ds,topN=N,panel_n=n,overlap_n=k,
                   overlap_genes=paste(ov,collapse="; "),p=p)
    }),fill=TRUE)
}),fill=TRUE)
TOPN[,FDR:=p.adjust(p,"BH")]

fwrite(PANEL,file.path(output_dir,"18D_full_panel_external_rank_enrichment.csv"))
fwrite(MODULE,file.path(output_dir,"18D_module_external_rank_enrichment.csv"))
fwrite(TOPN,file.path(output_dir,"18C_topN_overlap_results.csv"))

LOCKED <- data.table(
    external_dataset=c("WGL","DG"),
    expected_panel_n=c(28L,8L),
    expected_median=c(0.31382284,0.08234759),
    expected_background_median=c(0.5005553,0.5004550),
    expected_p=c(0.004423512,0.106188103),
    expected_FDR=c(0.008847024,0.106188103)
)
fwrite(LOCKED,file.path(output_dir,"external_CRISPR_locked_panel_targets.csv"))

cat("\nExternal ranked CRISPR validation complete.\n")
print(PANEL)
