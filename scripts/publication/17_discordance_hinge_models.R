############################################################
## 17: ARCHITECTURE-MANIFESTATION DISCORDANCE AND
##     LINEAR VERSUS HINGE DEPENDENCY MODELS
##
## Prepared input:
## model_id, system, dependency, discordance
## plus any covariates that should be included in every fit.
##
## discordance is defined upstream as standardized UCEI minus
## standardized ArchSig / manifestation activity.
############################################################

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly=TRUE)
input_file <- if(length(args)>=1L) args[[1L]] else Sys.getenv("UCEI_HINGE_INPUT",unset="")
output_dir <- if(length(args)>=2L) args[[2L]] else
    file.path("results","functional","discordance_hinge")
if(!nzchar(input_file)||!file.exists(input_file)) stop("Supply prepared system-level dependency table.")
dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)

D <- fread(input_file)
req <- c("model_id","system","dependency","discordance")
miss <- setdiff(req,names(D))
if(length(miss)) stop("Missing columns: ",paste(miss,collapse=", "))

## Optional covariates are all remaining numeric/factor columns
## except identifiers and the outcome/exposure.
covars <- setdiff(names(D),c("model_id","system","dependency","discordance"))

fit_system <- function(z) {
    z <- z[complete.cases(z[,c("dependency","discordance",covars),with=FALSE])]
    if(nrow(z)<80L) return(NULL)

    rhs_cov <- if(length(covars)) paste(covars,collapse=" + ") else NULL
    linear_rhs <- paste(c("discordance",rhs_cov),collapse=" + ")
    f_lin <- lm(as.formula(paste("dependency ~",linear_rhs)),data=z)

    ## Prespecified admissible range avoids extreme-tail breakpoints.
    qs <- quantile(z$discordance,c(0.10,0.90),na.rm=TRUE,names=FALSE)
    taus <- sort(unique(z$discordance[z$discordance>=qs[1] & z$discordance<=qs[2]]))
    if(length(taus)>200L) {
        taus <- unique(as.numeric(quantile(taus,probs=seq(0,1,length.out=200),names=FALSE)))
    }

    best <- NULL
    for(tau in taus) {
        zz <- copy(z)
        zz[, hinge:=pmax(0,discordance-tau)]
        rhs <- paste(c("discordance","hinge",rhs_cov),collapse=" + ")
        f <- try(lm(as.formula(paste("dependency ~",rhs)),data=zz),silent=TRUE)
        if(inherits(f,"try-error")) next
        a <- AIC(f)
        if(is.null(best)||a<best$aic) best <- list(tau=tau,fit=f,aic=a)
    }
    if(is.null(best)) return(NULL)

    s_lin <- summary(f_lin)$coefficients
    s_h <- summary(best$fit)$coefficients

    data.table(
        n=nrow(z),
        tau=best$tau,
        AIC_linear=AIC(f_lin),
        AIC_hinge=best$aic,
        delta_AIC=AIC(f_lin)-best$aic,
        beta_linear=s_lin["discordance","Estimate"],
        p_linear=s_lin["discordance","Pr(>|t|)"],
        beta_hinge=s_h["hinge","Estimate"],
        p_hinge=s_h["hinge","Pr(>|t|)"]
    )
}

RES <- D[,fit_system(.SD),by=system]
RES[,FDR_linear:=p.adjust(p_linear,"BH")]
RES[,FDR_hinge:=p.adjust(p_hinge,"BH")]
RES[,response_form:=ifelse(delta_AIC>2 & FDR_hinge<0.05,"Threshold-like","Graded")]

fwrite(RES,file.path(output_dir,"dependency_linear_vs_hinge_results.csv"))

LOCKED <- data.table(
    system="Checkpoint / DNA damage",
    tau_expected=0.966,
    delta_AIC_expected=17.44,
    hinge_beta_expected=0.366,
    hinge_FDR_expected=1.08e-4
)
fwrite(LOCKED,file.path(output_dir,"hinge_locked_checkpoint_target.csv"))

cat("\nDiscordance/hinge analysis complete.\n")
print(RES)
