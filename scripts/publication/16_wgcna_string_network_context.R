############################################################
## 16: WGCNA / STRING NETWORK CONTEXT
##
## Publication-facing reconstruction of the final recovered
## network-context layer.
##
## IMPORTANT:
## - This does NOT rerun WGCNA.
## - It consumes the surviving WGCNA CSV exports and the
##   retained STRING v12 physical-network edge table.
## - STRING associations are database-supported network
##   associations; they are not treated as universal proof
##   of direct physical binding.
############################################################

suppressPackageStartupMessages({
    library(data.table)
})

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 3L) {
    stop(
        "Usage: 16_wgcna_string_network_context.R ",
        "WGCNA_hub_gene_table.csv WGCNA_target_mechanism_gene_hub_table.csv ",
        "local_PPI_interactions.tsv [output_dir]"
    )
}

hub_file <- args[[1L]]
target_hub_file <- args[[2L]]
ppi_file <- args[[3L]]
output_dir <- if (length(args)>=4L) args[[4L]] else
    file.path("results","network","wgcna_string")

dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)

HUB <- fread(hub_file)
TARGET <- fread(target_hub_file)
PPI <- fread(ppi_file)

clean_gene <- function(x) {
    x <- toupper(trimws(as.character(x)))
    x <- sub("\\s*\\(.*\\)$","",x)
    x <- sub("\\|.*$","",x)
    x
}

pick_col <- function(nms,candidates,required=TRUE) {
    exact <- candidates[candidates %in% nms]
    if(length(exact)) return(exact[1L])
    for(p in candidates) {
        z <- grep(p,nms,ignore.case=TRUE,value=TRUE)
        if(length(z)) return(z[1L])
    }
    if(required) stop("Could not resolve required column from: ",paste(candidates,collapse=", "))
    NA_character_
}

target_gene_col <- pick_col(names(TARGET),c("gene","symbol","gene_symbol","Hugo_Symbol"))
TARGET[, gene:=clean_gene(get(target_gene_col))]
candidate_genes <- unique(TARGET$gene)

p1 <- pick_col(names(PPI),c("preferredName_A","protein1","gene1","from","node1","A"))
p2 <- pick_col(names(PPI),c("preferredName_B","protein2","gene2","to","node2","B"))

PPI[, gene_A:=clean_gene(get(p1))]
PPI[, gene_B:=clean_gene(get(p2))]

PPI_SUB <- PPI[
    gene_A %in% candidate_genes &
    gene_B %in% candidate_genes &
    gene_A != gene_B
]

PPI_SUB[, edge_key:=ifelse(gene_A<gene_B,
                           paste(gene_A,gene_B,sep="::"),
                           paste(gene_B,gene_A,sep="::"))]
PPI_SUB <- unique(PPI_SUB,by="edge_key")

DEGREE <- rbindlist(list(
    PPI_SUB[,.(gene=gene_A)],
    PPI_SUB[,.(gene=gene_B)]
))[,.(STRING_degree=.N),by=gene]

TARGET2 <- merge(TARGET,DEGREE,by="gene",all.x=TRUE)
TARGET2[is.na(STRING_degree),STRING_degree:=0L]

## Recover an exported WGCNA hub score if present.
hub_score_col <- pick_col(
    names(TARGET2),
    c("strict_coexpression_hub_score","coexpression_hub_score",
      "kME","module_membership","hub_score"),
    required=FALSE
)

if (!is.na(hub_score_col)) {
    TARGET2[, WGCNA_score:=as.numeric(get(hub_score_col))]
    TARGET2[, WGCNA_score_z:=as.numeric(scale(WGCNA_score))]
} else {
    TARGET2[, WGCNA_score:=NA_real_]
    TARGET2[, WGCNA_score_z:=0]
}

TARGET2[, STRING_degree_z:=if(sd(STRING_degree)>0)
            as.numeric(scale(STRING_degree)) else 0]

TARGET2[, STRING_plus_WGCNA_score:=
            STRING_degree_z + fifelse(is.finite(WGCNA_score_z),WGCNA_score_z,0)]

setorder(TARGET2,-STRING_plus_WGCNA_score,-STRING_degree)

AUDIT <- data.table(
    item=c(
        "WGCNA CSV outputs loaded",
        "STRING physical-network file loaded",
        "candidate STRING network edges",
        "candidate genes ranked",
        "TF activity inference"
    ),
    observed=c(
        nrow(HUB),
        nrow(PPI),
        nrow(PPI_SUB),
        nrow(TARGET2),
        0
    ),
    expected=c(NA,NA,256,59,0),
    interpretation=c(
        "WGCNA co-expression evidence recovered from exported tables",
        "STRING v12 network evidence supplied locally",
        "Historical final recovery reported 256 candidate edges",
        "Historical final recovery reported 59 ranked candidate genes",
        "No regulon/motif TF-activity claim"
    )
)

fwrite(PPI_SUB,file.path(output_dir,"STRING_candidate_network_edges.csv"))
fwrite(TARGET2,file.path(output_dir,"STRING_WGCNA_candidate_ranking.csv"))
fwrite(AUDIT,file.path(output_dir,"WGCNA_STRING_network_audit.csv"))

cat("\nWGCNA/STRING network-context recovery complete.\n")
print(AUDIT)
print(TARGET2[1:min(.N,30)])
