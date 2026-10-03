############################################################
## 14: SINGLE-CELL / SPATIAL SENSITIVITY ANALYSES
##
## This publication-facing script operates on the final
## prepared post-inference tables, not on raw GEO matrices.
##
## Single-cell input:
##   sample, clone_entropy, clone_cnv_amplitude,
##   ArchSig, CentromereKinetochore, Regulator, Proliferation
##
## Spatial input:
##   sample, CNV_like, Manifestation, ArchSig,
##   CentromereKinetochore, Regulator, Proliferation
############################################################
suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly=TRUE)
if(length(args)<2L) stop("Usage: script singlecell.csv spatial.csv [output_dir]")
SC <- fread(args[[1]])
SP <- fread(args[[2]])
output_dir <- if(length(args)>=3L) args[[3]] else file.path("results","sensitivity","singlecell_spatial")
dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)

modules <- c("ArchSig","CentromereKinetochore","Regulator","Proliferation")

## Sample-level clone-entropy associations
SC_SAMPLE <- rbindlist(lapply(modules,function(m){
    if(!all(c("sample","clone_entropy",m)%in%names(SC))) return(NULL)
    d <- unique(SC[,c("sample","clone_entropy",m),with=FALSE])
    z <- suppressWarnings(cor.test(d$clone_entropy,d[[m]],method="spearman",exact=FALSE))
    data.table(level="sample_entropy",module=m,n=nrow(d),rho=unname(z$estimate),p=z$p.value)
}),fill=TRUE)

## Clone-level amplitude associations
SC_CLONE <- rbindlist(lapply(modules,function(m){
    if(!all(c("clone_cnv_amplitude",m)%in%names(SC))) return(NULL)
    d <- SC[is.finite(clone_cnv_amplitude)&is.finite(get(m))]
    z <- suppressWarnings(cor.test(d$clone_cnv_amplitude,d[[m]],method="spearman",exact=FALSE))
    data.table(level="clone_amplitude",module=m,n=nrow(d),rho=unname(z$estimate),p=z$p.value)
}),fill=TRUE)

## Clone-level sample fixed-effect models
SC_FE <- rbindlist(lapply(modules,function(m){
    if(!all(c("sample","clone_cnv_amplitude",m)%in%names(SC))) return(NULL)
    d <- SC[complete.cases(SC[,c("sample","clone_cnv_amplitude",m),with=FALSE])]
    f <- lm(d[[m]] ~ d$clone_cnv_amplitude + factor(d$sample))
    s <- summary(f)$coefficients
    data.table(level="clone_fixed_effect",module=m,n=nrow(d),
               beta=s[2,"Estimate"],p=s[2,"Pr(>|t|)"])
}),fill=TRUE)

## Spatial sample-fixed-effect models
spmods <- c("Manifestation",modules)
SP_FE <- rbindlist(lapply(spmods,function(m){
    if(!all(c("sample","CNV_like",m)%in%names(SP))) return(NULL)
    d <- SP[complete.cases(SP[,c("sample","CNV_like",m),with=FALSE])]
    f <- lm(d[[m]] ~ d$CNV_like + factor(d$sample))
    s <- summary(f)$coefficients
    data.table(module=m,n=nrow(d),beta=s[2,"Estimate"],p=s[2,"Pr(>|t|)"])
}),fill=TRUE)
SP_FE[,FDR:=p.adjust(p,"BH")]

## Within-sample high/low contrasts
SP[,CNV_group:=ifelse(CNV_like > median(CNV_like,na.rm=TRUE),"high","low"),by=sample]
SP_PAIR <- rbindlist(lapply(spmods,function(m){
    if(!m%in%names(SP)) return(NULL)
    q <- SP[,.(diff=mean(get(m)[CNV_group=="high"],na.rm=TRUE)-
                    mean(get(m)[CNV_group=="low"],na.rm=TRUE)),by=sample]
    q <- q[is.finite(diff)]
    wt <- suppressWarnings(wilcox.test(q$diff,mu=0,exact=FALSE))
    bt <- binom.test(sum(q$diff>0),length(q$diff),p=0.5)
    data.table(module=m,samples=nrow(q),median_diff=median(q$diff),
               positive_samples=sum(q$diff>0),
               sign_p=bt$p.value,wilcox_p=wt$p.value)
}),fill=TRUE)

fwrite(SC_SAMPLE,file.path(output_dir,"singlecell_sample_entropy_associations.csv"))
fwrite(SC_CLONE,file.path(output_dir,"singlecell_clone_amplitude_associations.csv"))
fwrite(SC_FE,file.path(output_dir,"singlecell_clone_fixed_effect_models.csv"))
fwrite(SP_FE,file.path(output_dir,"spatial_fixed_effect_models.csv"))
fwrite(SP_PAIR,file.path(output_dir,"spatial_within_sample_contrasts.csv"))

cat("\nSingle-cell/spatial sensitivity analysis complete.\n")
