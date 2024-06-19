# COR="D35prevent19_stage2_severe";
# COR="D35prevent19_stage2_delta";
renv::activate(project = here::here(".."))
Sys.setenv(TRIAL = "prevent19_stage2")
Sys.setenv(VERBOSE = 1)
source(here::here("..", "_common.R")) # dat_proc is made


{
  library(kyotil) # p.adj.perm, getFormattedSummary
  library(marginalizedRisk)
  library(tools) # toTitleCase
  library(survey)
  library(plotrix) # weighted.hist
  library(parallel)
  library(forestplot)
  library(Hmisc) # wtd.quantile, cut2
  library(xtable) # this is a dependency of kyotil
  source(here::here("code", "params.R"))
  time.start = Sys.time()
  myprint(study_name)
  myprint(verbose)
  
  
  # path for figures and tables etc
  save.results.to = here::here("output")
  if (!dir.exists(save.results.to))
    dir.create(save.results.to)
  save.results.to = paste0(save.results.to, "/", attr(config, "config"))
  if (!dir.exists(save.results.to))
    dir.create(save.results.to)
  save.results.to = paste0(save.results.to, "/", COR, "/")
  if (!dir.exists(save.results.to))
    dir.create(save.results.to)
  print(paste0("save.results.to equals ", save.results.to))
  
  # append to file names for figures and tables
  fname.suffix = ""
  
  myprint(tfinal.tpeak)
  write(tfinal.tpeak, file = paste0(save.results.to, "timepoints_cum_risk_" %.% study_name))
  
  
  if (COR == "D35prevent19_stage2_delta") {
    # # if doing competing risk estimation
    # comp.risk = T
    # dat_proc$EventIndOfInterest = dat_proc$DeltaCOVIDIndD35_108to21Dec10
    # dat_proc$EventIndCompeting  = dat_proc$NonDeltaCOVIDIndD35_108to21Dec10
    # form.0 = update(
    #   Surv(COVIDTimeD35_to21Dec10, EventIndOfInterest) ~ 1,
    #   as.formula(config$covariates)
    # )
    # dat_proc$yy = dat_proc$EventIndOfInterest
    
    form.0 = update(
      Surv(COVIDTimeD35to21Dec10, KnownOrImputedDeltaCOVIDIndD35_108to21Dec10) ~ 1,
      as.formula(config$covariates)
    )
    dat_proc$yy = dat_proc$KnownOrImputedDeltaCOVIDIndD35_108to21Dec10
    
  } else if (COR == "D35prevent19_stage2_severe") {
    form.0 = update(
      Surv(SevereCOVIDTimeD35to21Dec10, SevereCOVIDIndD35_108to21Dec10) ~ 1,
      as.formula(config$covariates)
    )
    dat_proc$yy = dat_proc$SevereCOVIDIndD35_108to21Dec10
  }
  
  for (a in c("Day35"%.%assays)) {
    dat_proc[[a%.%"centered"]] = scale(dat_proc[[a]], scale=F)
  }
  
  dat.vac.seroneg = subset(dat_proc, Trt == 1 & ph1)
  dat.pla.seroneg = subset(dat_proc, Trt == 0 & ph1)
  dat.vac.seroneg.ph2 = subset(dat.vac.seroneg, ph2)
  
  
  # define trichotomized markers
  if (is.null(attr(dat_proc, "marker.cutpoints"))) {
    dat.vac.seroneg = add.trichotomized.markers (dat.vac.seroneg, all.markers, wt.col.name =
                                                   "wt")
    marker.cutpoints = attr(dat.vac.seroneg, "marker.cutpoints")
  } else {
    marker.cutpoints = attr(dat_proc, "marker.cutpoints")
  }
  # save cutpoints to files
  for (a in all.markers) {
    q.a = marker.cutpoints[[a]]
    if (startsWith(a, "Day")) {
      # not fold change
      write(
        paste0(labels.axis[1, marker.name.to.assay(a)], " [", concatList(round(q.a, 2), ", "), ")%"),
        file = paste0(save.results.to, "cutpoints_", a)
      )
    } else {
      # fold change
      write(
        paste0(escape(a), " [", concatList(round(q.a, 2), ", "), ")%"),
        file = paste0(save.results.to, "cutpoints_", a)
      )
    }
  }
  
  #create twophase design object
  design.vacc.seroneg <-
    twophase(
      id = list( ~ 1,  ~ 1),
      strata = list(NULL,  ~ Wstratum),
      subset =  ~ ph2,
      data = dat.vac.seroneg
    )
  with(dat.vac.seroneg, table(Wstratum, ph2))
  
  # table of ph1 and ph2 cases
  tab = with(dat.vac.seroneg, table(ph2, EventIndPrimary))
  names(dimnames(tab))[2] = "Event Indicator"
  print(tab)
  mytex(
    tab,
    file.name = "tab1_" %.% fname.suffix,
    save2input.only = T,
    input.foldername = save.results.to
  )
  
  begin = Sys.time()
}


###################################################################################################
# estimate overall VE in the placebo and vaccine arms

cor_coxph_risk_no_marker (
  form.0,
  dat=dat.vac.seroneg,
  fname.suffix, 
  save.results.to,
  config,
  config.cor,
  tfinal.tpeak,
  
  dat.plac = NULL,
  verbose=FALSE
) 



###################################################################################################
# Univariate models

cor_coxph_coef_1(
  form.0,
  design_or_dat = design.vacc.seroneg,
  fname.suffix="",
  save.results.to,
  config,
  config.cor,
  all.markers,
  all.markers.names.short,
  
  dat.pla.seroneg = NULL,
  show.q = F,
  verbose = T
)



###################################################################################################
# Quadratic models

if (COR == "D35prevent19_stage2_severe") {
  
  # quadratic models
  cor_coxph_coef_n (
    form.0,
    design_or_dat = design.vacc.seroneg,
    fname.suffix="Quadratic",
    save.results.to,
    config,
    config.cor,
    all.markers = sapply(assays, function (a) paste0("Day35",a, "centered + I(Day35",a, "centered^2)")),
    all.markers.names.short,
    
    nCoef=2,
    col.headers=c("D35","D35**2"),
    verbose = T)
  
}


###################################################################################################
# marginalized risk and controlled VE
###################################################################################################

# # if competing risk
# form.0 = list(form.0, as.formula(
#   sub(
#     "EventIndOfInterest",
#     "EventIndCompeting",
#     paste0(deparse(form.0, width.cutoff = 500))
#   )
# )),

markers = "Day35" %.% c("pseudoneutid50_D614G", "pseudoneutid50_Delta", "bindSpike_D614", "bindSpike_Delta1") # save time, only need these

cor_coxph_risk_bootstrap(
  form.0,
  dat = dat.vac.seroneg,
  fname.suffix,
  save.results.to,
  config,
  config.cor,
  tfinal.tpeak,
  
  markers = markers,

  run.Sgts = F # whether to get risk conditional on continuous S>=s
)


cor_coxph_risk_plotting (
  form.0,
  dat = dat.vac.seroneg,
  fname.suffix,
  save.results.to,
  config,
  config.cor,
  tfinal.tpeak,
  
  markers = markers,
  markers.names.short = all.markers.names.short[markers],
  markers.names.long = all.markers.names.long[markers],
  marker.cutpoints,
  assay_metadata,
  
  dat.plac = NULL,
  res.plac.cont = NULL,
  prev.plac = NULL,
  overall.ve=NULL,
  
  show.ve.curves = F,
  plot.geq = F,
  plot.w.plac = F,
  for.title = ""
)



################################################################################
print(date())
print("cor_coxph run time: " %.% format(Sys.time() - time.start, digits = 1))
