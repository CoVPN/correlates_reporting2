# COR="D31nextgen_mock";
# COR="D31nextgen_mock_tcell";
renv::activate(project = here::here(".."))
Sys.setenv(TRIAL = "nextgen_mock")
Sys.setenv(VERBOSE = 1)
source(here::here("..", "_common.R")) 


# hack, also hack in _common.R L143 and in Rmd L26
all.markers=setdiff(all.markers, c("Day31bindSpike_IgG_N","Day31bindSpike_IgA_N"))

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
  fname.suffix = "Obj1"
  
  myprint(tfinal.tpeak)
  write(tfinal.tpeak, file = paste0(save.results.to, "timepoints_cum_risk_" %.% fname.suffix))
  
  
  dat_proc$yy = dat_proc$COVIDIndD31_7toM12

  for (a in c("Day31"%.%assays)) {
    dat_proc[[a%.%"centered"]] = scale(dat_proc[[a]], scale=F)
  }
  
  dat.vac = subset(dat_proc, Trt == 1 & ph1)
  dat.plac = subset(dat_proc, Trt == 0 & ph1)
  dat.vac.ph2 = subset(dat.vac, ph2==1)
  
  
  # define trichotomized markers
  if (is.null(attr(dat_proc, "marker.cutpoints"))) {
    # do it in _common.R
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
        file = paste0(save.results.to, "cutpoints_", a,".txt")
      )
    } else {
      # fold change
      write(
        paste0(escape(a), " [", concatList(round(q.a, 2), ", "), ")%"),
        file = paste0(save.results.to, "cutpoints_", a,".txt")
      )
    }
  }
  
  #create twophase design object
  design.vacc <-
    twophase(
      id = list( ~ 1,  ~ 1),
      strata = list(NULL,  ~ Wstratum),
      subset =  ~ ph2,
      data = dat.vac
    )
  with(dat.vac, table(Wstratum, ph2))
  
  # table of ph1 and ph2 cases
  tab = with(dat.vac, table(ph2, EventIndPrimary))
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
  dat=dat.vac,
  fname.suffix, 
  save.results.to,
  config,
  config.cor,
  tfinal.tpeak,
  
  dat.plac = dat.plac,
  verbose=FALSE
) 



###################################################################################################
# Univariate models

cor_coxph_coef_1(
  form.0,
  design_or_dat = design.vacc,
  fname.suffix,
  save.results.to,
  config,
  config.cor,
  all.markers,
  all.markers.names.short,
  
  dat.plac,
  show.q = F,
  verbose = T
)


cor_coxph_risk_tertile_incidence_curves (
  form.0,
  dat = dat.vac,
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

  dat.plac,
  for.title = ""
)


# ###################################################################################################
# # marginalized risk and controlled VE
# ###################################################################################################
# 
# # # if competing risk
# # form.0 = list(form.0, as.formula(
# #   sub(
# #     "EventIndOfInterest",
# #     "EventIndCompeting",
# #     paste0(deparse(form.0, width.cutoff = 500))
# #   )
# # )),
# 
# markers = "Day31" %.% c("pseudoneutid50_D614G", "pseudoneutid50_Delta", "bindSpike_D614", "bindSpike_Delta1") # save time, only need these
# 
# cor_coxph_risk_bootstrap(
#   form.0,
#   dat = dat.vac,
#   fname.suffix,
#   save.results.to,
#   config,
#   config.cor,
#   tfinal.tpeak,
#   
#   markers = markers,
# 
#   run.Sgts = F # whether to get risk conditional on continuous S>=s
# )
# 
# cor_coxph_risk_plotting (
#   form.0,
#   dat = dat.vac,
#   fname.suffix,
#   save.results.to,
#   config,
#   config.cor,
#   tfinal.tpeak,
#   
#   markers = markers,
#   markers.names.short = all.markers.names.short[markers],
#   markers.names.long = all.markers.names.long[markers],
#   marker.cutpoints,
#   assay_metadata,
#   
#   dat.plac = NULL,
#   res.plac.cont = NULL,
#   prev.plac = NULL,
#   overall.ve=NULL,
#   
#   show.ve.curves = F,
#   plot.geq = F,
#   plot.w.plac = F,
#   for.title = ""
# )




################################################################################
print(date())
print("cor_coxph run time: " %.% format(Sys.time() - time.start, digits = 1))
