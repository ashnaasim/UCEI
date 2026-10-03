############################################################
## 15: PRISM PHARMACOLOGIC SENSITIVITY ANALYSIS
##
## Prepared long-format input:
## ModelID, compound, AUC, discordance, ArchSig,
## global_dependency, lineage, UCEI
##
## Optional: drug_class
############################################################
suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly=TRUE)
input_file <- if(length(args)>=1L) args[[1]] else Sys.getenv("UCEI_PRISM_LONG",unset="")
output_dir <- if(length(args)>=2L) args[[2]] else file.path("results","sensitivity","prism")
if(!nzchar(input_file)||!file.exists(input_file)) stop("Supply prepared PRISM long table.")
dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)

D <- fread(input_file)
req <- c("ModelID","compound","AUC","discordance","ArchSig","global_dependency","lineage","UCEI")
miss <- setdiff(req,names(D))
if(length(miss)) stop("Missing columns: ",paste(miss,collapse=", "))

D[,AUC:=as.numeric(AUC)]
## Higher score = greater sensitivity, as in final workflow.
D[,sensitivity:=-as.numeric(scale(AUC)),by=compound]

fit_compound <- function(z,predictor){
    vars <- c("sensitivity",predictor,"ArchSig","global_dependency","lineage")
    z <- z[complete.cases(z[,..vars])]
    if(nrow(z)<80L || sd(z$sensitivity)==0) return(NULL)
    f <- lm(
        as.formula(paste("sensitivity ~",predictor,"+ ArchSig + global_dependency + factor(lineage)")),
        data=z
    )
    s <- summary(f)$coefficients
    if(!predictor%in%rownames(s)) return(NULL)
    data.table(
        n=nrow(z),
        beta=s[predictor,"Estimate"],
        SE=s[predictor,"Std. Error"],
        p=s[predictor,"Pr(>|t|)"]
    )
}

R1 <- D[,fit_compound(.SD,"discordance"),by=compound]
R1[,predictor:="discordance"]
R2 <- D[,fit_compound(.SD,"UCEI"),by=compound]
R2[,predictor:="UCEI"]
RES <- rbindlist(list(R1,R2),fill=TRUE)
RES[,FDR:=p.adjust(p,"BH"),by=predictor]
RES[,supported:=FDR<0.10]

fwrite(RES,file.path(output_dir,"prism_compound_models.csv"))

if("drug_class"%in%names(D)){
    MAP <- unique(D[,.(compound,drug_class)])
    C <- merge(RES,MAP,by="compound",all.x=TRUE)
    CLASS <- C[,.(compounds=.N,
                  median_beta=median(beta,na.rm=TRUE),
                  FDR_supported=sum(supported,na.rm=TRUE)),by=.(predictor,drug_class)]
    fwrite(CLASS,file.path(output_dir,"prism_drug_class_summary.csv"))
}

LOCKED <- data.table(expected_aligned_models=342L,expected_unique_compounds=1502L,
                     expected_FDR_supported_classes=0L)
fwrite(LOCKED,file.path(output_dir,"prism_locked_targets.csv"))
cat("\nPRISM sensitivity analysis complete.\n")
print(RES[,.(tested=.N,FDR_supported=sum(supported)),by=predictor])
