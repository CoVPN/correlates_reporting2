# COR="D31toM6_nextgen_mock";
# COR="D31toM6_nextgen_mock_tcell";
renv::activate(project = here::here(".."))
Sys.setenv(TRIAL = "nextgen_mock")
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
  
  # hack
  # source("~/copcor/R/cor_coxph_coef_1.R")
  
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
  
  for (a in c("Day31"%.%assays)) {
    dat_proc[[a%.%"centered"]] = scale(dat_proc[[a]], scale=F)
  }
  
  dat.vacc = subset(dat_proc, Trt == 1 & ph1)
  dat.plac = subset(dat_proc, Trt == 0 & ph1)
  # dat.vacc.ph2 = subset(dat.vacc, ph2==1)
  
  design.vacc <-
    twophase(
      id = list( ~ 1,  ~ 1),
      strata = list(NULL,  ~ Wstratum),
      subset =  ~ ph2,
      data = dat.vacc
    )

  design.plac <-
    twophase(
      id = list( ~ 1,  ~ 1),
      strata = list(NULL,  ~ Wstratum),
      subset =  ~ ph2,
      data = dat.plac
    )
  
  # define trichotomized markers
  if (is.null(attr(dat_proc, "marker.cutpoints"))) {
    # better create marker.cutpoints attr in _common.R
  } else {
    marker.cutpoints = attr(dat_proc, "marker.cutpoints")
  }
  # save cutpoints to files
  for (a in c("Day31"%.%assays, "B"%.%assays, "Delta31overB"%.%assays)) 
    write(paste0(labels.axis[1, marker.name.to.assay(a)], " [", concatList(round(marker.cutpoints[[a]], 2), ", "), ")%"),
      file = paste0(save.results.to, "cutpoints_", a,".txt"))

  begin = Sys.time()
}


###################################################################################################
# estimate overall VE in the placebo and vaccine arms

# append to file names for figures and tables
fname.suffix = "ExpVacc"

cor_coxph_risk_no_marker (
  form.0,
  dat=dat.vacc,
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

trts=c(1,0); marker_sets = 1:3 
# trts=1; marker_sets = 1 

for (trt in trts) {
  if (trt==1) {
    dat.1=dat.vacc; design.1 = design.vacc
    dat.0=dat.plac
    fname.suffix.0 = "ExpVacc"
    trt.label="ExpVacc"
    cmp.label="CtlVacc"

  } else {
    dat.1=dat.plac; design.1 = design.plac
    dat.0=dat.vacc
    fname.suffix.0 = "CtlVacc"
    trt.label="CtlVacc"
    cmp.label="ExpVacc"
  }
  
  
  # table of ph1 and ph2 cases
  tab1 = with(dat.1, table(ph2, EventIndPrimary))
  names(dimnames(tab1))[2] = "Event Indicator"; print(tab1)
  
  for (marker_set in marker_sets) {
    
    if (marker_set==1) {
      fname.suffix = fname.suffix.0%.%"_D31"
      
      all.markers=c(paste0("Day", tpeak, assays))
      all.markers.names.short = sub("Pseudovirus-", "", assay_metadata$assay_label_short[match(assays,assay_metadata$assay)])
      all.markers.names.short = sub(" \\(AU/ml\\)", "", sub("Anti Spike ", "", all.markers.names.short))
      all.markers.names.short = c("D"%.%tpeak%.%" "%.%all.markers.names.short)
      
    } else if (marker_set==2) {
      fname.suffix = fname.suffix.0%.%"_B"
      
      all.markers=c(paste0("B", assays) )
      all.markers.names.short = sub("Pseudovirus-", "", assay_metadata$assay_label_short[match(assays,assay_metadata$assay)])
      all.markers.names.short = sub(" \\(AU/ml\\)", "", sub("Anti Spike ", "", all.markers.names.short))
      all.markers.names.short = c("B "%.%all.markers.names.short)
      
    } else if (marker_set==3) {
      fname.suffix = fname.suffix.0%.%"_D31overB"
      
      all.markers=c(paste0("Delta", tpeak, "overB", assays))
      all.markers.names.short = sub("Pseudovirus-", "", assay_metadata$assay_label_short[match(assays,assay_metadata$assay)])
      all.markers.names.short = sub(" \\(AU/ml\\)", "", sub("Anti Spike ", "", all.markers.names.short))
      all.markers.names.short = c("D"%.%tpeak%.%"/B "%.%all.markers.names.short    )
      
    }
    names(all.markers.names.short) = all.markers
    
    # need to save tab1 for each distinct fname.suffix
    mytex(tab1, file.name = "tab1_" %.% fname.suffix, save2input.only = T, input.foldername = save.results.to)
  
    cat("\n\n"); myprint(fname.suffix)
    
    cor_coxph_coef_1(
      form.0,
      design_or_dat = design.vacc,
      fname.suffix,
      save.results.to,
      config,
      config.cor,
      markers = all.markers,
      markers.names.short = all.markers.names.short,
  
      dat.plac = dat.0,
      show.q = F,
      
      forestplot.markers=NULL, 
      for.title="",
      run.trichtom=TRUE,
      cmp.label = cmp.label,
      verbose = T
    )
  
    cor_coxph_risk_tertile_incidence_curves (
      form.0,
      dat = dat.1,
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
  
      dat.plac = dat.0,
      for.title = "",
      
      trt.label = trt.label,
      cmp.label = cmp.label
    )
  
  
  }

}

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
#   dat = dat.vacc,
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
#   dat = dat.vacc,
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
