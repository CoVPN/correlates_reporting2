# COR="D31toM12_nextgen_mock_sera";
# COR="D31toM12_nextgen_mock_tcell";
Sys.setenv(TRIAL = "nextgen_mock")
Sys.setenv(VERBOSE = 1)
source(here::here("..", "_common.R")) 

marker_sets = unique(assay_metadata$panel)
# hack, for now only work on one marker set
marker_sets=c('pseudoneutid50_sera')

{
  quiet_library("SuperLearner")
  quiet_library("tidyverse")
  quiet_library("dplyr")
  quiet_library("glue")
  source("code/format_utils.R")
  get.short.name=function(assays){
    all.markers.names.short = sub("Pseudovirus-", "", assay_metadata$assay_label_short[match(assays,assay_metadata$assay)])
    all.markers.names.short = sub(" \\(AU/ml\\)", "", sub("Anti Spike ", "", all.markers.names.short))
    all.markers.names.short
  }
  
  time.start = Sys.time()
  myprint(verbose)
  
  cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ running cop_mediation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
  
  # path to save output
  save.results.to = here::here("output")
  if (!dir.exists(save.results.to)) dir.create(save.results.to)
  save.results.to = paste0(save.results.to, "/", attr(config, "config"))
  if (!dir.exists(save.results.to)) dir.create(save.results.to)
  save.results.to = paste0(save.results.to, "/", COR, "/")
  if (!dir.exists(save.results.to)) dir.create(save.results.to)
  myprint(save.results.to)
    
  dat = subset(dat_proc, ph1==1)

  # remove .+ from right hand side
  f.adj=gsub("~\\s*\\.\\s*\\+\\s*", "~", config$covariates)
  # -1 to remove intercept. better to do it after calling model.matrix b/c if done on the formula with a factor, it does not work right
  cov.df = as.data.frame(model.matrix(as.formula(f.adj), dat)[,-1,drop=F])
  # to ensure that the code works when config$covariates contains factor, we need to rename the columns of cov.df
  names(cov.df) = "cov_adj_"%.%1:ncol(cov.df)
  cov.df.f = concatList(names(cov.df),"+")
  
  # set Super Learner libraries
  sl_library <- c(
    "SL.mean",
    "SL.glm",
    # "SL.glmnet",
    # "SL.xgboost",
    # "SL.ranger",
    "SL.gam",
    "SL.earth"
  )
  
}

for (marker_set in marker_sets) {
  
  fname.suffix = marker_set
  
  assays = subset(assay_metadata, panel==marker_set, assay, drop=T)
  if (len(assays)==0) {
    warning("no assays found")
    next
  }
  
  all.markers=c(paste0("Day", tpeak, assays))
  
  tmp = get.short.name(assays); all.markers.names.short = c(glue("D{tpeak} {tmp}"))
  names(all.markers.names.short) = all.markers; 
  all.markers.names.long = c(as.matrix(labels.title)[DayPrefix%.%tpeak, assays])
  names(all.markers.names.long) = all.markers
  
  names(all.markers) = all.markers.names.short
  
  results_list = lapply (all.markers, function(a) {
    if (verbose) myprint(a)
    
    # set marker values outside ph2 to NA
    dat[[a]][dat$ph2==0]=NA
  
    set.seed(612); 
    if (verbose) cat("fitting (1,1) ...")
    fit_trt1_medtrt1 <- survtmle::hazard_tmle(
      ftime = dat$EventTimePrimary,
      ftype = dat$EventIndPrimary,
      trt = dat$Trt,
      # adjust for covariates of interest - I think this needs to be a dataframe
      adjustVars = cov.df,
      # mediator must be a df
      mediator = data.frame(dat[a]),
      # mediator trt value = M(a)
      mediatorTrtVal = 1,
      # trt of interest = T(a, M(1)),
      trtOfInterest = 1,
      # weights depend on whether marker is a bAb or nAb and stage
      mediatorSampProb = 1 / dat$wt,
      # keep these
      mediatorInCensMod = FALSE,
      mediatorStratify.ftime = TRUE,
      # depends on stage
      t0 = tfinal.tpeak,
      # get values from SL library
      # SL.ctime = sl_library,
      # SL.ftime = sl_library,
      # SL.mediator = sl_library,
      # SL.trtMediator = sl_library,
      # SL.eif = sl_library,
      glm.ctime = paste0("splines::ns(t, 3) + ", cov.df.f),
      glm.ftime = paste0("splines::ns(t, 3) + ", cov.df.f, " + ", a),
      glm.mediator = cov.df.f,
      glm.trtMediator = paste0(cov.df.f, " + ", a),
      glm.eif = ".",
      verbose = TRUE,
      # changing these might help with model fit
      maxIter = 1,
      gtol = 0.05,
      gtolCens = 0.05,
      truncateH = 0.9,
      returnModels = TRUE
    )
        
    # P[T(0, M(0)) < tau]
    set.seed(828)
    if (verbose) cat("fitting (0,0) ...")
    fit_trt0_medtrt0 <- survtmle::hazard_tmle(
      ftime = dat$EventTimePrimary,
      ftype = dat$EventIndPrimary,
      trt = dat$Trt,
      adjustVars = cov.df,
      mediator = data.frame(dat[a]),
      mediatorTrtVal = 0,
      trtOfInterest = 0,
      mediatorSampProb = 1 / dat$wt,
      mediatorInCensMod = FALSE,
      mediatorStratify.ftime = TRUE,
      t0 = tfinal.tpeak,
      # SL.ctime = sl_library,
      # SL.ftime = sl_library,
      # SL.mediator = sl_library,
      # SL.trtMediator = sl_library,
      # SL.eif = sl_library,
      glm.ctime = "splines::ns(t, 3) + "%.%cov.df.f,
      glm.ftime = paste0("splines::ns(t, 3) + ", cov.df.f, " + ", a),
      glm.mediator = cov.df.f,
      glm.trtMediator = paste0(cov.df.f, " + ", a),
      glm.eif = ".",
      verbose = TRUE,
      maxIter = 1,
      gtol = 0.05,
      gtolCens = 0.05,
      truncateH = 0.9,
      returnModels = TRUE
    )
        
    # print(fit_trt1_medtrt1); print(fit_trt0_medtrt0)
        
    # bind first two fits together
    fit_trt_equal_medtrt <- list(
      est = rbind(fit_trt0_medtrt0$est, fit_trt1_medtrt1$est),
      ic = cbind(fit_trt0_medtrt0$ic, fit_trt1_medtrt1$ic)
    )
        
    # P[T(1, M(0)) < tau]
    # set.seed(901)
    set.seed(777)
    if (verbose) cat("fitting (0,1) ...")
    fit_trt1_medtrt0 <- survtmle::hazard_tmle(
      ftime = dat$EventTimePrimary,
      ftype = dat$EventIndPrimary,
      trt = dat$Trt,
      adjustVars = cov.df,
      mediator = data.frame(dat[a]),
      mediatorTrtVal = 0,
      trtOfInterest = 1,
      mediatorSampProb = 1 / dat$wt,
      mediatorInCensMod = FALSE,
      mediatorStratify.ftime = TRUE,
      t0 = tfinal.tpeak,
      # SL.ctime = sl_library,
      # SL.ftime = sl_library,
      # SL.mediator = sl_library,
      # SL.trtMediator = sl_library,
      # SL.eif = sl_library,
      glm.ctime = "splines::ns(t, 3) + "%.%cov.df.f,
      glm.ftime = paste0("splines::ns(t, 3) + ", cov.df.f, " + ", a),
      glm.mediator = cov.df.f,
      glm.trtMediator = paste0(cov.df.f, " + ", a),
      glm.eif = ".",
      verbose = TRUE,
      maxIter = 1,
      gtol = 0.05,
      gtolCens = 0.05,
      truncateH = 0.9,
      returnModels = TRUE
    )
        
    # print(fit_trt1_medtrt0)
        
    fit <- compute_mediation_params(fit_trt_equal_medtrt, fit_trt1_medtrt0)
        
    fit$eff
    
  }) # returns a list of fit$eff for a marker_set
  
  saveRDS(results_list, file = glue('{save.results.to}eff_{marker_set}.RDS'))
  
  
  
} # end marker_set


