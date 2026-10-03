score_ucei_v2 <- function(input, parameter_dir='.') {
  features <- c('jsd_neutral','ent_state_size','ent_chr_sd','ent_chr_mean','ent_ampclass','ent_transition')
  required <- c('cancer','burden_total',features)
  missing <- setdiff(required,names(input))
  if (length(missing)) stop(paste('Missing required columns:',paste(missing,collapse=', ')))
  rp <- utils::read.csv(file.path(parameter_dir,'UCEI_v2_residualization_parameters.csv'),stringsAsFactors=FALSE,check.names=FALSE)
  ld <- utils::read.csv(file.path(parameter_dir,'UCEI_v2_PCA_loadings.csv'),stringsAsFactors=FALSE,check.names=FALSE)
  sp <- utils::read.csv(file.path(parameter_dir,'UCEI_v2_score_parameters.csv'),stringsAsFactors=FALSE,check.names=FALSE)
  cancer <- tolower(as.character(input$cancer))
  burden <- suppressWarnings(as.numeric(input$burden_total))
  if (any(!is.finite(burden))) stop('burden_total contains non-finite values')
  raw <- numeric(nrow(input))
  for (f in features) {
    p <- rp[rp$feature==f,,drop=FALSE]
    ii <- match(cancer,tolower(as.character(p$cancer)))
    if (anyNA(ii)) stop(paste('Unsupported cancer label(s):',paste(unique(cancer[is.na(ii)]),collapse=', ')))
    x <- suppressWarnings(as.numeric(input[[f]]))
    if (any(!is.finite(x))) stop(paste(f,'contains non-finite values'))
    rsd <- p$residual_sd[ii]
    if (any(!is.finite(rsd) | rsd<=0)) stop(paste('Invalid frozen residual SD for',f))
    pred <- p$intercept[ii] + p$slope[ii]*burden
    rz <- (x-pred-p$residual_mean[ii])/rsd
    wi <- match(f,ld$raw_feature)
    if (is.na(wi)) stop(paste('Missing PCA loading for',f))
    raw <- raw + rz*as.numeric(ld$loading[wi])
  }
  (raw-as.numeric(sp$score_center[1]))/as.numeric(sp$score_scale[1])
}

score_ucei_v2_file <- function(input_csv, output_csv, parameter_dir='.') {
  d <- utils::read.csv(input_csv,stringsAsFactors=FALSE,check.names=FALSE)
  d$UCEI_v2_z <- score_ucei_v2(d,parameter_dir)
  utils::write.csv(d,output_csv,row.names=FALSE)
  invisible(d)
}
