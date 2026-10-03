############################################################
## 10: LATENT ARCHITECTURE-MANIFESTATION ANALYSIS
##
## Prepared input must contain:
## cancer, burden_total,
## six burden-adjusted architecture features,
## seven transcriptomic module scores,
## and optionally PFI variables.
############################################################
suppressPackageStartupMessages({
    library(data.table)
    library(survival)
})

args <- commandArgs(trailingOnly=TRUE)
input_file <- if(length(args)>=1) args[[1]] else Sys.getenv("UCEI_LATENT_INPUT",unset="")
output_dir <- if(length(args)>=2) args[[2]] else file.path("results","biology","latent_architecture_manifestation")
if(!nzchar(input_file)||!file.exists(input_file)) stop("Supply prepared latent-factor table.")
dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)

D <- if(grepl("\\.rds$",input_file,ignore.case=TRUE)) as.data.table(readRDS(input_file)) else fread(input_file)

arch <- c(
    "jsd_neutral_resid_z","ent_state_size_resid_z","ent_chr_sd_resid_z",
    "ent_chr_mean_resid_z","ent_ampclass_resid_z","ent_transition_resid_z"
)
mods <- c(
    "ArchSig","Core_Mitotic","Spindle","Kinetochore",
    "Centromere","E2F_MYC_FOXM1","Replication_Stress_DNARepair"
)

resolve <- function(targets,nms){
    out <- character(length(targets))
    for(i in seq_along(targets)){
        t <- targets[i]
        cands <- unique(c(
            t,
            sub("_resid_z$","_v2_resid_z",t),
            gsub("_","",tolower(t))
        ))
        hit <- nms[tolower(nms)%in%tolower(cands)]
        if(!length(hit)) {
            ## permissive normalized-name fallback
            nn <- gsub("[^a-z0-9]","",tolower(nms))
            tt <- gsub("[^a-z0-9]","",tolower(t))
            hit <- nms[nn==tt]
        }
        if(!length(hit)) stop("Could not resolve column: ",t)
        out[i] <- hit[1]
    }
    out
}
arch_cols <- resolve(arch,names(D))
mod_cols <- resolve(mods,names(D))

Xarch <- scale(as.matrix(D[,..arch_cols]))
Xmod <- scale(as.matrix(D[,..mod_cols]))
ok <- complete.cases(Xarch)&complete.cases(Xmod)&complete.cases(D[,.(cancer,burden_total)])
D2 <- D[ok]
Xarch <- Xarch[ok,,drop=FALSE]
Xmod <- Xmod[ok,,drop=FALSE]

pa <- prcomp(Xarch,center=FALSE,scale.=FALSE)
pm <- prcomp(Xmod,center=FALSE,scale.=FALSE)

latent_arch <- as.numeric(pa$x[,1])
manifest <- as.numeric(pm$x[,1])

## Orient both factors to positive relation with intuitive summaries.
if(cor(latent_arch,rowMeans(Xarch),use="complete.obs")<0) latent_arch <- -latent_arch
if(cor(manifest,rowMeans(Xmod),use="complete.obs")<0) manifest <- -manifest

D2[,latent_architecture:=as.numeric(scale(latent_arch))]
D2[,manifestation:=as.numeric(scale(manifest))]
D2[,burden_z:=as.numeric(scale(burden_total))]

fit <- lm(manifestation ~ latent_architecture + burden_z + factor(cancer),data=D2)
sf <- summary(fit)
LATENT_ASSOC <- data.table(
    n=nrow(D2),
    beta=coef(fit)["latent_architecture"],
    SE=sf$coefficients["latent_architecture","Std. Error"],
    p=sf$coefficients["latent_architecture","Pr(>|t|)"],
    R2=sf$r.squared,
    spearman_rho=cor(D2$latent_architecture,D2$manifestation,method="spearman")
)

VARIANCE <- data.table(
    factor=c("latent_architecture","manifestation"),
    variance_explained=c(
        pa$sdev[1]^2/sum(pa$sdev^2),
        pm$sdev[1]^2/sum(pm$sdev^2)
    )
)

fwrite(VARIANCE,file.path(output_dir,"latent_factor_variance.csv"))
fwrite(LATENT_ASSOC,file.path(output_dir,"latent_architecture_manifestation_association.csv"))
fwrite(data.table(feature=arch,loading=pa$rotation[,1]),file.path(output_dir,"latent_architecture_loadings.csv"))
fwrite(data.table(module=mods,loading=pm$rotation[,1]),file.path(output_dir,"manifestation_loadings.csv"))

if(all(c("pfi_time","pfi_event")%in%names(D2))){
    S <- D2[complete.cases(pfi_time,pfi_event)]
    m0 <- coxph(Surv(pfi_time,pfi_event)~burden_z+strata(cancer),data=S)
    m1 <- coxph(Surv(pfi_time,pfi_event)~burden_z+latent_architecture+manifestation+strata(cancer),data=S)
    SURV <- data.table(
        model=c("baseline","architecture_plus_manifestation"),
        Cindex=c(summary(m0)$concordance[1],summary(m1)$concordance[1])
    )
    fwrite(SURV,file.path(output_dir,"latent_survival_comparison.csv"))
}
cat("\nLatent architecture-manifestation analysis complete.\n")
print(VARIANCE)
print(LATENT_ASSOC)
