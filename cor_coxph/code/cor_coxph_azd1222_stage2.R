# COR="D57azd1222_stage2_severe_nAb";
# COR="D57azd1222_stage2_severe_bAb";
# COR="D57azd1222_stage2_delta_nAb";
# COR="D57azd1222_stage2_delta_bAb";
# COR="D57azd1222_stage2_delta_nAb_sens";
# COR="D57azd1222_stage2_delta_bAb_sens";
renv::activate(project = here::here(".."))
Sys.setenv(TRIAL = "azd1222_stage2")
Sys.setenv(VERBOSE = 1)
source(here::here("..", "_common.R")) 


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
  
  form.0 = update(
    as.formula(paste0("Surv(",config.cor$EventTimePrimary,", ",config.cor$EventIndPrimary,") ~ 1")),
    as.formula(config.cor$covariates)
  )
  
  dat_proc$yy = dat_proc[[config.cor$EventIndPrimary]]
  
  for (a in c("Day57"%.%assays)) {
    dat_proc[[a%.%"centered"]] = scale(dat_proc[[a]], scale=F)
  }
  
  dat.vac.seroneg = subset(dat_proc, Trt == 1 & ph1)
  dat.pla.seroneg = subset(dat_proc, Trt == 0 & ph1)
  dat.vac.seroneg.ph2 = subset(dat.vac.seroneg, ph2)
  
  
  marker.cutpoints = attr(dat_proc, "marker.cutpoints")
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
  # there are only 4 ph2 ptids in non-US, senior population, none of who are cases
  with(dat.vac.seroneg, table(Wstratum, ph2))
  # subset(dat.vac.seroneg, Country!=2 & Senior & ph2, EventIndPrimary)
  
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
# marginalized risk and controlled VE

if (endsWith(COR, "nAb")) {
  markers = all.markers 
} else {
  markers = "Day57" %.% c("bindSpike_D614", "bindSpike_Delta1") # save time, only need these
}


if (contain(COR, "severe")) {

  cor_coxph_risk_bootstrap(
    form.0,
    dat = dat.vac.seroneg,
    fname.suffix,
    save.results.to,
    config,
    config.cor,
    tfinal.tpeak,
    
    markers=all.markers,

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
    
    markers = all.markers,
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
  
  cor_coxph_risk_tertile_incidence_curves (
    form.0,
    dat = dat.vac.seroneg,
    fname.suffix,
    save.results.to,
    config,
    config.cor,
    tfinal.tpeak,
    
    markers = all.markers,
    markers.names.short = all.markers.names.short[markers],
    markers.names.long = all.markers.names.long[markers],
    marker.cutpoints,
    assay_metadata,
    
    dat.plac = NULL,
    for.title = ""
  )
  
  
} else if (contain(COR, "delta")) {
  # competing risk
  
  cor_coxph_risk_bootstrap(
    form.0,
    dat = dat.vac.seroneg,
    fname.suffix,
    save.results.to,
    config,
    config.cor,
    tfinal.tpeak,
    
    markers = all.markers,

    run.Sgts = F # whether to get risk conditional on continuous S>=s
  )
  
  cor_coxph_risk_plotting(
    form.0,
    dat = dat.vac.seroneg,
    fname.suffix,
    save.results.to,
    config,
    config.cor,
    tfinal.tpeak,
    
    markers = all.markers,
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
  
  cor_coxph_risk_tertile_incidence_curves(
    form.0,
    dat = dat.vac.seroneg,
    fname.suffix,
    save.results.to,
    config,
    config.cor,
    tfinal.tpeak,
    
    markers = all.markers,
    markers.names.short = all.markers.names.short[markers],
    markers.names.long = all.markers.names.long[markers],
    marker.cutpoints,
    assay_metadata,
    
    dat.plac = NULL,
    for.title = ""
  )
  
}


print(date())
print("cor_coxph run time: " %.% format(Sys.time() - time.start, digits =
                                          1))
