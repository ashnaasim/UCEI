############################################################
## 13: PROTEOMIC INTEGRATION AND PATHWAY ORA
##
## The exact historical four-class decision table must be
## supplied rather than re-inferred.
##
## Inputs:
## 1) classified_gene_file: gene, proteomic_class
## 2) pathway_membership_file: pathway, gene, source
## 3) background_file: one column named gene
############################################################
suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly=TRUE)
if(length(args)<3L) stop("Usage: script classified_genes.csv pathway_membership.csv background.csv [output_dir]")
class_file <- args[[1]]
path_file <- args[[2]]
bg_file <- args[[3]]
output_dir <- if(length(args)>=4L) args[[4]] else file.path("results","multiomic","proteomics_pathways")
dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)

C <- fread(class_file)
P <- fread(path_file)
B <- fread(bg_file)
if(!all(c("gene","proteomic_class")%in%names(C))) stop("classified file needs gene, proteomic_class")
if(!all(c("pathway","gene")%in%names(P))) stop("pathway file needs pathway, gene")
if(!"gene"%in%names(B)) stop("background file needs gene")

C[,gene:=toupper(gene)]
P[,gene:=toupper(gene)]
B[,gene:=toupper(gene)]

bg <- unique(B$gene)
C <- unique(C[gene%in%bg,.(gene,proteomic_class)])
P <- unique(P[gene%in%bg,.(pathway,gene,source=if("source"%in%names(P)) source else NA_character_)])

CLASS_COUNTS <- C[,.N,by=proteomic_class][order(-N)]
combined <- unique(C$gene)

test_set <- function(query,label){
    query <- intersect(unique(query),bg)
    N <- length(bg); n <- length(query)
    z <- P[,{
        members <- unique(gene)
        K <- length(members)
        k <- sum(query%in%members)
        .(K=K,k=k)
    },by=.(pathway,source)]
    z <- z[K>=10 & K<=500 & k>=2]
    z[,p:=phyper(k-1,K,N-K,n,lower.tail=FALSE)]
    z[,fold_enrichment:=(k/n)/(K/N)]
    z[,FDR:=p.adjust(p,"BH")]
    z[,query_class:=label]
    z
}

ORA <- rbindlist(c(
    list(test_set(combined,"combined_165")),
    lapply(split(C$gene,C$proteomic_class),function(g){
        test_set(g,unique(C[gene%in%g,proteomic_class])[1])
    })
),fill=TRUE)
setorder(ORA,query_class,FDR,-fold_enrichment)

fwrite(CLASS_COUNTS,file.path(output_dir,"proteomic_class_counts.csv"))
fwrite(ORA,file.path(output_dir,"proteomic_pathway_ORA.csv"))
fwrite(ORA[FDR<0.05],file.path(output_dir,"proteomic_pathway_ORA_significant.csv"))

LOCKED_COUNTS <- data.table(
    proteomic_class=c(
        "Protein/RNA-independent hidden dependencies",
        "Protein-over-RNA decoupling",
        "Protein-linked dependencies",
        "Protein-compensation candidates"
    ),
    expected_n=c(68L,64L,28L,5L)
)
fwrite(LOCKED_COUNTS,file.path(output_dir,"proteomic_locked_class_counts.csv"))
cat("\nProteomic/pathway analysis complete.\n")
print(CLASS_COUNTS)
cat("Significant combined pathways:",ORA[query_class=="combined_165" & FDR<0.05,.N],"\n")
