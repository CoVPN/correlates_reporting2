# COR="D35prevent19nvx";
renv::activate(project = here::here(".."))
Sys.setenv(TRIAL = "prevent19nvx")
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
  # defined differently in cor_coxph_xx.R
  fname.suffix = ""
  
  myprint(tfinal.tpeak)
  write(tfinal.tpeak,
        file = paste0(save.results.to, "timepoints_cum_risk_" %.% study_name))
  
  dat_proc$yy = dat_proc$EventIndPrimary
  
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
  
  
  # some exploratory code
  if (config$is_ows_trial) source(here::here("code", "cor_coxph_misc.R"))
  
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

# don't run this, instead 
# source(here::here("code", "cor_coxph_risk_no_marker.R"))

cat("read from prevent19 estimates b/c we have to do case control sampling, which does not work for placebo\n")

load("output/prevent19/D35/marginalized.risk.no.marker..Rdata")

print(cbind(prev.plac, prev.vacc, overall.ve))


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
  
  dat.pla.seroneg = dat.pla.seroneg,
  show.q = F,
  verbose = T
)




###################################################################################################
# marginalized risk and controlled VE
###################################################################################################


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

cor_coxph_risk_plotting (
  form.0,
  dat = dat.vac.seroneg,
  fname.suffix,
  save.results.to,
  config,
  config.cor,
  tfinal.tpeak,
  
  markers = all.markers,
  markers.names.short = all.markers.names.short,
  markers.names.long = all.markers.names.long,
  marker.cutpoints,
  assay_metadata,
  
  dat.plac = dat.pla.seroneg,
  res.plac.cont,
  prev.plac,
  overall.ve,
  
  show.ve.curves = T,
  plot.geq = F,
  plot.w.plac = T,
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
  markers.names.short = all.markers.names.short,
  markers.names.long = all.markers.names.long,
  marker.cutpoints,
  assay_metadata,
  
  dat.plac = dat.pla.seroneg,
  for.title = ""
)


print(date())
print("cor_coxph run time: " %.% format(Sys.time() - time.start, digits=1))
