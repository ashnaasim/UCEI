############################################################
## 19: EXTERNAL REDUCED-FEATURE UCEI VALIDATION
##
## Adapted from the executed SCRIPT 09A cBioPortal workflow.
## This is a reduced CNA-entropy proxy and MUST NOT be described
## as exact reconstruction of the six-feature UCEI.
############################################################

suppressPackageStartupMessages({
    library(data.table)
    library(survival)
})

set.seed(20260710)

args <- commandArgs(trailingOnly=TRUE)
study_root <- if(length(args)>=1L) args[[1L]] else Sys.getenv("UCEI_EXTERNAL_STUDY_ROOT",unset="")
output_dir <- if(length(args)>=2L) args[[2L]] else file.path("results","external","reduced_ucei")
if(!nzchar(study_root)||!dir.exists(study_root)) {
    stop("Supply a directory containing extracted cBioPortal study packs.")
}
dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)

study_ids <- c(
    "prad_mskcc",
    "prad_p1000",
    "prad_su2c_2019",
    "prostate_msk_2024"
)

safe_num <- function(x) suppressWarnings(as.numeric(x))
zscale <- function(x) {
    x <- safe_num(x); s <- sd(x,na.rm=TRUE); m <- mean(x,na.rm=TRUE)
    if(!is.finite(s)||s==0) return(rep(NA_real_,length(x)))
    as.numeric((x-m)/s)
}
shannon <- function(p) {
    p <- p[p>0 & is.finite(p)]
    if(!length(p)) return(NA_real_)
    -sum(p*log(p))
}
norm_id <- function(x) gsub("\\.","-",as.character(x))

find_file <- function(root,patterns) {
    ff <- list.files(root,recursive=TRUE,full.names=TRUE)
    hits <- ff[Reduce(\`|\`,lapply(patterns,function(p) grepl(p,basename(ff),ignore.case=TRUE)))]
    if(!length(hits)) return(NA_character_)
    hits[1L]
}
read_cbio <- function(path) {
    if(is.na(path)||!file.exists(path)) return(data.table())
    lines <- readLines(path,warn=FALSE)
    tmp <- tempfile(fileext=".txt")
    writeLines(lines[!grepl("^#",lines)],tmp)
    tryCatch(fread(tmp,sep="\t",header=TRUE,fill=TRUE),error=function(e)data.table())
}

compute_gene_cna_entropy <- function(cna_dt) {
    if(!nrow(cna_dt)) return(data.table())
    sample_cols <- setdiff(names(cna_dt),c(
        "Hugo_Symbol","Gene","gene","Hugo Symbol",
        "Entrez_Gene_Id","Entrez Gene Id","Cytoband"
    ))
    if(length(sample_cols)<20L) return(data.table())
    mat <- as.matrix(cna_dt[,..sample_cols])
    storage.mode(mat) <- "numeric"

    rbindlist(lapply(seq_along(sample_cols),function(j){
        x <- mat[,j]; x <- x[is.finite(x)]
        if(length(x)<100L) return(NULL)
        states <- c(-2,-1,0,1,2)
        counts <- table(factor(x,levels=states))
        p <- as.numeric(counts)/sum(counts)
        amp_frac <- mean(x>0)
        del_frac <- mean(x<0)
        high_amp_frac <- mean(x>=2)
        deep_del_frac <- mean(x<=-2)
        altered_frac <- mean(x!=0)
        state_entropy <- shannon(p)
        ampdel_balance <- 1-abs(amp_frac-del_frac)
        proxy <- state_entropy +
            0.50*altered_frac +
            0.25*ampdel_balance +
            0.25*(high_amp_frac+deep_del_frac)

        data.table(
            SAMPLE_ID=sample_cols[j],
            n_genes_cna=length(x),
            state_entropy=state_entropy,
            amp_frac=amp_frac,
            del_frac=del_frac,
            high_amp_frac=high_amp_frac,
            deep_del_frac=deep_del_frac,
            altered_frac=altered_frac,
            ampdel_balance=ampdel_balance,
            reduced_UCEI_proxy_raw=proxy
        )
    }),fill=TRUE)[,reduced_UCEI_proxy_z:=zscale(reduced_UCEI_proxy_raw)]
}

detect_binary <- function(clin) {
    if(!nrow(clin)) return(data.table())
    id_cols <- intersect(c("SAMPLE_ID","PATIENT_ID","Sample Identifier","Patient Identifier"),names(clin))
    if(!length(id_cols)) return(data.table())
    id <- if("SAMPLE_ID"%in%names(clin)) "SAMPLE_ID" else id_cols[1L]

    z <- copy(clin)
    z[,endpoint_text:=apply(.SD,1,function(r)paste(tolower(as.character(r)),collapse=" | "))]
    z[,metastatic_binary:=fifelse(
        grepl("metastatic|metastasis|metastases|bone|liver|lung|lymph node",endpoint_text) &
        !grepl("non.metastatic|no metast",endpoint_text),1L,
        fifelse(grepl("primary|localized|locoregional|prostate gland|radical prostatectomy",endpoint_text),0L,NA_integer_)
    )]
    z[,aggressive_binary:=metastatic_binary]

    gleason_cols <- grep("gleason|grade",names(z),value=TRUE,ignore.case=TRUE)
    if(length(gleason_cols)) {
        g <- apply(z[,..gleason_cols],1,function(r){
            q <- suppressWarnings(as.numeric(gsub("[^0-9.]","",as.character(r))))
            q <- q[is.finite(q)]
            if(!length(q)) NA_real_ else max(q)
        })
        z[is.na(aggressive_binary)&is.finite(g)&g>=8,aggressive_binary:=1L]
        z[is.na(aggressive_binary)&is.finite(g)&g<8,aggressive_binary:=0L]
    }

    z[,.(SAMPLE_ID=norm_id(get(id)),metastatic_binary,aggressive_binary,endpoint_text)]
}

detect_survival <- function(clin) {
    if(!nrow(clin)) return(data.table())
    id_cols <- intersect(c("SAMPLE_ID","PATIENT_ID","Sample Identifier","Patient Identifier"),names(clin))
    if(!length(id_cols)) return(data.table())
    id <- if("SAMPLE_ID"%in%names(clin)) "SAMPLE_ID" else id_cols[1L]
    nms <- names(clin)
    tc <- grep("OS_MONTHS|OVERALL_SURVIVAL_MONTHS|SURVIVAL_MONTHS|PFS_MONTHS|DFS_MONTHS|RFS_MONTHS",
               nms,value=TRUE,ignore.case=TRUE)[1L]
    sc <- grep("OS_STATUS|OVERALL_SURVIVAL_STATUS|PFS_STATUS|DFS_STATUS|RFS_STATUS|STATUS",
               nms,value=TRUE,ignore.case=TRUE)[1L]
    if(is.na(tc)||is.na(sc)) return(data.table())

    out <- clin[,.(SAMPLE_ID=norm_id(get(id)),
                   survival_time=safe_num(get(tc)),
                   survival_status_raw=as.character(get(sc)))]
    out[,survival_event:=fifelse(
        grepl("1|deceased|dead|progressed|recurred|relapse|event",tolower(survival_status_raw)),1L,
        fifelse(grepl("0|living|alive|diseasefree|disease free|censored",tolower(survival_status_raw)),0L,NA_integer_)
    )]
    out[is.finite(survival_time)&survival_time>0&survival_event%in%c(0,1)]
}

analyze_one <- function(study_id) {
    root <- file.path(study_root,study_id)
    if(!dir.exists(root)) return(list(
        audit=data.table(study_id=study_id,status="MISSING_STUDY_PACK"),
        binary=data.table(),survival=data.table(),scores=data.table()
    ))

    cna_file <- find_file(root,c("^data_cna\\.txt$","data_discrete_cna","gistic"))
    cs <- read_cbio(find_file(root,c("data_clinical_sample")))
    cp <- read_cbio(find_file(root,c("data_clinical_patient")))
    clin <- rbindlist(list(cs,cp),fill=TRUE)
    cna <- read_cbio(cna_file)
    scores <- compute_gene_cna_entropy(cna)
    if(!nrow(scores)) return(list(
        audit=data.table(study_id=study_id,status="FAIL_NO_USABLE_CNA"),
        binary=data.table(),survival=data.table(),scores=data.table()
    ))
    scores[,SAMPLE_ID:=norm_id(SAMPLE_ID)]
    b <- detect_binary(clin)
    s <- detect_survival(clin)
    if(nrow(b)) scores <- merge(scores,b,by="SAMPLE_ID",all.x=TRUE)
    if(nrow(s)) scores <- merge(scores,s,by="SAMPLE_ID",all.x=TRUE)

    BR <- list()
    for(ep in c("metastatic_binary","aggressive_binary")) {
        if(!ep%in%names(scores)) next
        d <- scores[get(ep)%in%c(0,1)&is.finite(reduced_UCEI_proxy_z)]
        if(nrow(d)>=30L && length(unique(d[[ep]]))==2L && min(table(d[[ep]]))>=5L) {
            wt <- wilcox.test(reduced_UCEI_proxy_z~get(ep),data=d)
            gl <- glm(as.formula(paste(ep,"~ reduced_UCEI_proxy_z")),data=d,family=binomial())
            sm <- summary(gl)$coefficients
            BR[[length(BR)+1L]] <- data.table(
                study_id=study_id,endpoint=ep,n=nrow(d),
                n_event=sum(d[[ep]]==1),n_control=sum(d[[ep]]==0),
                delta_event_minus_control=
                    mean(d[get(ep)==1,reduced_UCEI_proxy_z])-
                    mean(d[get(ep)==0,reduced_UCEI_proxy_z]),
                wilcox_p=wt$p.value,
                logistic_beta=sm["reduced_UCEI_proxy_z","Estimate"],
                logistic_OR=exp(sm["reduced_UCEI_proxy_z","Estimate"]),
                logistic_p=sm["reduced_UCEI_proxy_z","Pr(>|z|)"]
            )
        }
    }
    BR <- rbindlist(BR,fill=TRUE)

    SR <- data.table()
    if(all(c("survival_time","survival_event")%in%names(scores))) {
        d <- scores[is.finite(reduced_UCEI_proxy_z)&is.finite(survival_time)&survival_time>0&survival_event%in%c(0,1)]
        if(nrow(d)>=50L && sum(d$survival_event)>=10L) {
            fit <- coxph(Surv(survival_time,survival_event)~reduced_UCEI_proxy_z,data=d,ties="efron")
            sm <- summary(fit)
            SR <- data.table(
                study_id=study_id,n=fit$n,events=fit$nevent,
                logHR=sm$coefficients["reduced_UCEI_proxy_z","coef"],
                HR=sm$conf.int["reduced_UCEI_proxy_z","exp(coef)"],
                lower95=sm$conf.int["reduced_UCEI_proxy_z","lower .95"],
                upper95=sm$conf.int["reduced_UCEI_proxy_z","upper .95"],
                p=sm$coefficients["reduced_UCEI_proxy_z","Pr(>|z|)"]
            )
        }
    }

    list(
        audit=data.table(study_id=study_id,status="PASS",n_cna_samples=nrow(scores),
                         n_binary_results=nrow(BR),n_survival_results=nrow(SR)),
        binary=BR,survival=SR,scores=scores
    )
}

R <- setNames(lapply(study_ids,analyze_one),study_ids)
AUDIT <- rbindlist(lapply(R,\`[[\`,"audit"),fill=TRUE)
BINARY <- rbindlist(lapply(R,\`[[\`,"binary"),fill=TRUE)
SURVIVAL <- rbindlist(lapply(R,\`[[\`,"survival"),fill=TRUE)

for(id in names(R)) {
    if(nrow(R[[id]]$scores)) fwrite(R[[id]]$scores,file.path(output_dir,paste0("sample_scores_",id,".csv")))
}

if(nrow(BINARY)) {
    BINARY[,direction_support:=delta_event_minus_control>0]
    BINARY[,nominal_support:=wilcox_p<0.05|logistic_p<0.05]
}
if(nrow(SURVIVAL)) {
    SURVIVAL[,direction_support:=HR>1]
    SURVIVAL[,nominal_support:=p<0.05]
}

SUMMARY <- data.table(
    validation_layer=c("External binary phenotype","External survival phenotype"),
    n_results=c(nrow(BINARY),nrow(SURVIVAL)),
    n_direction_supported=c(
        if(nrow(BINARY))sum(BINARY$direction_support)else 0,
        if(nrow(SURVIVAL))sum(SURVIVAL$direction_support)else 0
    ),
    n_nominal_supported=c(
        if(nrow(BINARY))sum(BINARY$nominal_support)else 0,
        if(nrow(SURVIVAL))sum(SURVIVAL$nominal_support)else 0
    )
)

fwrite(AUDIT,file.path(output_dir,"external_reduced_ucei_audit.csv"))
fwrite(BINARY,file.path(output_dir,"external_reduced_ucei_binary_results.csv"))
fwrite(SURVIVAL,file.path(output_dir,"external_reduced_ucei_survival_results.csv"))
fwrite(SUMMARY,file.path(output_dir,"external_reduced_ucei_summary.csv"))

cat("\nExternal reduced-feature UCEI analysis complete.\n")
print(AUDIT)
print(SUMMARY)