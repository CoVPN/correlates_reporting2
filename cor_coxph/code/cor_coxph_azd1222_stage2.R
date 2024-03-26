# COR="D57azd1222_stage2_severe_nAb";
# COR="D57azd1222_stage2_severe_bAb";
renv::activate(project = here::here(".."))
Sys.setenv(TRIAL = "azd1222_stage2")
Sys.setenv(VERBOSE = 1)
source(here::here("..", "_common.R")) # dat.mock is made


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
  write(tfinal.tpeak,
        file = paste0(save.results.to, "timepoints_cum_risk_" %.% study_name))
  
  
  if (COR == "D57azd1222_stage2_delta_nAb" | COR == "D57azd1222_stage2_delta_bAb") {
    # for use in competing risk estimation
    comp.risk = T
    dat.mock$EventIndOfInterest = dat.mock$DeltaEventIndD57_120to21Dec10
    dat.mock$EventIndCompeting  = dat.mock$NonDeltaEventIndD57_120to21Dec10
    
    form.0 = update(
      Surv(EventTimeD57_to21Dec10, EventIndOfInterest) ~ 1,
      as.formula(config$covariates)
    )
    dat.mock$yy = dat.mock$EventIndOfInterest
    
  } else if (COR == "D57azd1222_stage2_severe_nAb" | COR == "D57azd1222_stage2_severe_bAb") {
    form.0 = update(
      Surv(
        SevereEventTimeD57_to21Dec10,
        SevereEventIndD57_120to21Dec10
      ) ~ 1,
      as.formula(config$covariates)
    )
    dat.mock$yy = dat.mock$SevereEventIndD57_120to21Dec10
  }
  
  for (a in c("Day57"%.%assays)) {
    dat.mock[[a%.%"centered"]] = scale(dat.mock[[a]], scale=F)
  }
  
  dat.vac.seroneg = subset(dat.mock, Trt == 1 & ph1)
  dat.pla.seroneg = subset(dat.mock, Trt == 0 & ph1)
  dat.vac.seroneg.ph2 = subset(dat.vac.seroneg, ph2)
  
  
  # define trichotomized markers
  if (is.null(attr(dat.mock, "marker.cutpoints"))) {
    dat.vac.seroneg = add.trichotomized.markers (dat.vac.seroneg, all.markers, wt.col.name =
                                                   "wt")
    marker.cutpoints = attr(dat.vac.seroneg, "marker.cutpoints")
  } else {
    marker.cutpoints = attr(dat.mock, "marker.cutpoints")
  }
  
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

source(here::here("code", "cor_coxph_risk_no_marker.R"))

if (Sys.getenv("COR_COXPH_NO_MARKER_ONLY") == 1) q("no")



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
###################################################################################################

if (endsWith(COR, "nAb")) {
  markers = all.markers 
} else {
  markers = "Day57" %.% c("bindSpike_D614", "bindSpike_Delta1") # save time, only need these
}


if (contain(COR, "severe")) {
  # non-competing risk
  
  cor_coxph_risk_bootstrap(
    form.0,
    dat = dat.vac.seroneg,
    fname.suffix,
    save.results.to,
    config,
    config.cor,

    tfinal.tpeak,
    all.markers = markers,

    comp.risk = F,
    run.Sgts = F # whether to get risk conditional on continuous S>=s
  )

  
  cor_coxph_risk_plotting (
    form.0,
    dat = dat.vac.seroneg,
    fname.suffix,
    save.results.to,
    config,
    config.cor,
    
    assay_metadata,
    
    tfinal.tpeak,
    all.markers = markers,
    all.markers.names.short[markers],
    all.markers.names.long[markers],
    marker.cutpoints,

    multi.imp = F,
    comp.risk = F,

    dat.pla.seroneg = NULL,
    res.plac.cont = NULL,
    prev.plac = NULL,

    variant = NULL,

    show.ve.curves = F,
    plot.geq = F,
    plot.w.plac = F,
    for.title = ""
  )


} else if (contain(COR, "delta")) {
  # competing risk
  
  cor_coxph_risk_bootstrap(
    form.0 = list(form.0, as.formula(
      sub(
        "EventIndOfInterest",
        "EventIndCompeting",
        paste0(deparse(form.0, width.cutoff = 500))
      )
    )),
    dat,
    fname.suffix,
    save.results.to,

    tfinal.tpeak,
    all.markers = markers,

    numCores,
    B,

    comp.risk = T,
    run.Sgts = F # whether to get risk conditional on continuous S>=s
  )

  cor_coxph_risk_plotting(
    form.0 = list(form.0, as.formula(
      sub(
        "EventIndOfInterest",
        "EventIndCompeting",
        paste0(deparse(form.0, width.cutoff = 500))
      )
    )),
    dat,
    fname.suffix,
    save.results.to,
    config,
    config.cor,
    
    assay_metadata,

    tfinal.tpeak,
    all.markers = markers,
    all.markers.names.short[markers],
    all.markers.names.long[markers],
    marker.cutpoints,

    multi.imp = F,
    comp.risk = T,

    dat.pla.seroneg = NULL,
    res.plac.cont = NULL,
    prev.plac = NULL,

    variant = NULL,

    show.ve.curves = F,
    plot.geq = F,
    plot.w.plac = F,
    for.title = ""
  )

}


print(date())
print("cor_coxph run time: " %.% format(Sys.time() - time.start, digits =
                                          1))
