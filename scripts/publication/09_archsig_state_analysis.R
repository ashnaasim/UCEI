############################################################
## 09: UCEI-ARCHSIG ARCHITECTURE-TRANSCRIPTOME STATES
##
## Publication-facing reconstruction from the final Stage-3
## biology analysis.
##
## Required columns:
## patient_id, cancer, pfi_time, pfi_event,
## UCEI, ArchSig, burden_total
##
## Optional clinical-complete-case columns may be supplied in
## a second file with a precomputed clinical_covariate_score
## or explicit cancer-specific covariates for separate models.
############################################################
suppressPackageStartupMessages({
    library(data.table)
    library(survival)
})

args <- commandArgs(trailingOnly=TRUE)
input_file <- if(length(args)>=1) args[[1]] else Sys.getenv("UCEI_ARCHSIG_INPUT",unset="")
output_dir <- if(length(args)>=2) args[[2]] else file.path("results","biology","archsig_states")
if(!nzchar(input_file)||!file.exists(input_file)) stop("Supply prepared UCEI-ArchSig table.")
dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)

D <- if(grepl("\\.rds$",input_file,ignore.case=TRUE)) as.data.table(readRDS(input_file)) else fread(input_file)

pick <- function(x,cands,label){
    z <- intersect(cands,names(x))
    if(!length(z)) stop("Missing ",label," column.")
    z[1]
}
uc <- pick(D,c("UCEI_z","UCEI_v2_z","UCEI","ucei"),"UCEI")
ac <- pick(D,c("ArchSig","ArchSig_Core","archsig","archsig_score"),"ArchSig")
bc <- pick(D,c("burden_total_z","burden_total","CNA_burden"),"burden")
req <- c("cancer","pfi_time","pfi_event",uc,ac,bc)
D <- D[complete.cases(D[,..req])]
D[,cancer:=tolower(as.character(cancer))]
D[,UCEI:=as.numeric(get(uc))]
D[,ArchSig_score:=as.numeric(get(ac))]
D[,burden:=as.numeric(get(bc))]
D[,pfi_time:=as.numeric(pfi_time)]
D[,pfi_event:=as.integer(pfi_event)]

## Cancer-specific medians, exactly as used for state definition.
D[, ucei_high := UCEI > median(UCEI,na.rm=TRUE), by=cancer]
D[, archsig_high := ArchSig_score > median(ArchSig_score,na.rm=TRUE), by=cancer]

D[, state4 := fifelse(!ucei_high & !archsig_high,"Low UCEI / Low ArchSig",
              fifelse(ucei_high & !archsig_high,"High UCEI / Low ArchSig",
              fifelse(!ucei_high & archsig_high,"Low UCEI / High ArchSig",
                      "High UCEI / High ArchSig")))]
D[, state4:=factor(state4,levels=c(
    "Low UCEI / Low ArchSig",
    "High UCEI / Low ArchSig",
    "Low UCEI / High ArchSig",
    "High UCEI / High ArchSig"
))]

STATE_COUNTS <- D[,.(n=.N,events=sum(pfi_event),event_rate=mean(pfi_event)),by=state4]

fit4 <- coxph(
    Surv(pfi_time,pfi_event) ~ state4 + burden + strata(cancer),
    data=D,ties="efron",x=TRUE
)
s4 <- summary(fit4)
STATE_HR <- data.table(
    term=rownames(s4$coefficients),
    beta=s4$coefficients[,"coef"],
    HR=s4$conf.int[,"exp(coef)"],
    CI_low=s4$conf.int[,"lower .95"],
    CI_high=s4$conf.int[,"upper .95"],
    p=s4$coefficients[,"Pr(>|z|)"]
)[grepl("^state4",term)]

D[, joint_high := as.integer(ucei_high & archsig_high)]

fit_jh_cancer <- coxph(
    Surv(pfi_time,pfi_event) ~ joint_high + strata(cancer),
    data=D,ties="efron",x=TRUE
)
fit_jh_burden <- coxph(
    Surv(pfi_time,pfi_event) ~ joint_high + burden + strata(cancer),
    data=D,ties="efron",x=TRUE
)

extract_jh <- function(fit,label){
    s <- summary(fit)
    data.table(
        model=label,
        HR=s$conf.int["joint_high","exp(coef)"],
        CI_low=s$conf.int["joint_high","lower .95"],
        CI_high=s$conf.int["joint_high","upper .95"],
        p=s$coefficients["joint_high","Pr(>|z|)"],
        Cindex=as.numeric(s$concordance[1]),
        n=fit$n,
        events=fit$nevent
    )
}
JOINT_HIGH <- rbindlist(list(
    extract_jh(fit_jh_cancer,"cancer adjusted"),
    extract_jh(fit_jh_burden,"cancer + burden adjusted")
))

## Continuous interaction test
m_add <- coxph(
    Surv(pfi_time,pfi_event) ~ scale(UCEI)+scale(ArchSig_score)+burden+strata(cancer),
    data=D,ties="efron"
)
m_int <- coxph(
    Surv(pfi_time,pfi_event) ~ scale(UCEI)*scale(ArchSig_score)+burden+strata(cancer),
    data=D,ties="efron"
)
INTERACTION <- data.table(
    LRT=2*(as.numeric(logLik(m_int))-as.numeric(logLik(m_add))),
    df=length(coef(m_int))-length(coef(m_add))
)
INTERACTION[,p:=pchisq(LRT,df=df,lower.tail=FALSE)]

## Cancer-specific rank coupling
CANCER_COR <- D[,{
    z <- suppressWarnings(cor.test(UCEI,ArchSig_score,method="spearman",exact=FALSE))
    .(n=.N,rho=unname(z$estimate),p=z$p.value)
},by=cancer]
CANCER_COR[,FDR:=p.adjust(p,"BH")]

fwrite(STATE_COUNTS,file.path(output_dir,"archsig_state_counts.csv"))
fwrite(STATE_HR,file.path(output_dir,"archsig_state_HR.csv"))
fwrite(JOINT_HIGH,file.path(output_dir,"archsig_joint_high_models.csv"))
fwrite(INTERACTION,file.path(output_dir,"archsig_continuous_interaction.csv"))
fwrite(CANCER_COR,file.path(output_dir,"archsig_by_cancer_correlations.csv"))

LOCKED <- data.table(
    quantity=c("n","events","LL_rate","HL_rate","LH_rate","HH_rate"),
    expected=c(3265,783,0.142,0.211,0.258,0.343)
)
cat("\nArchSig state analysis complete.\n")
print(STATE_COUNTS)
print(JOINT_HIGH)
print(INTERACTION)
