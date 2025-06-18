# COR="D31toM12_nextgen_mock_sera";
# COR="D31toM12_nextgen_mock_tcell";
Sys.setenv(TRIAL = "nextgen_mock")
Sys.setenv(VERBOSE = 1)
source(here::here("..", "_common.R")) 


trts=c(1,0)
marker_sets = unique(assay_metadata$panel)
marker_sets='pseudoneutid50_sera'
# marker_set='pseudoneutid50_sera'


{
  library("kyotil") # p.adj.perm, getFormattedSummary
  quiet_library("marginalizedRisk")
  quiet_library("tools") # toTitleCase
  quiet_library("survey")
  quiet_library("glue")
  quiet_library("plotrix") # weighted.hist
  quiet_library("parallel")
  quiet_library("forestplot")
  quiet_library("Hmisc") # wtd.quantile, cut2
  quiet_library("xtable") # this is a dependency of kyotil
  
  source(here::here("code", "params.R"))
  time.start = Sys.time()
  myprint(verbose)
  
  # hack
  # source("~/copcor/R/cor_coxph_coef_1.R")
  
  cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ running cor_coxph ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
  
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
  myprint(save.results.to)
  
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

  
  get.short.name=function(assays){
    all.markers.names.short = sub("Pseudovirus-", "", assay_metadata$assay_label_short[match(assays,assay_metadata$assay)])
    all.markers.names.short = sub(" \\(AU/ml\\)", "", sub("Anti Spike ", "", all.markers.names.short))
    all.markers.names.short
  }
  
  begin = Sys.time()
}


###################################################################################################
# estimate overall VE in the placebo and vaccine arms

# append to file names for figures and tables
fname.suffix = "InvVacc"

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

for (trt in trts) {
  
  if (trt==1) {
    dat.1=dat.vacc; design.1 = design.vacc
    dat.0=dat.plac
    fname.suffix.0 = "InvVacc"
    trt.label="InvVacc"
    cmp.label="CtlVacc"

  } else {
    dat.1=dat.plac; design.1 = design.plac
    dat.0=dat.vacc
    fname.suffix.0 = "CtlVacc"
    trt.label="CtlVacc"
    cmp.label="InvVacc"
  }
  
  # table of ph1 and ph2 cases
  tab1 = with(dat.1, table(ph2, EventIndPrimary))
  names(dimnames(tab1))[2] = "Event Indicator"; print(tab1)
  
  for (marker_set in marker_sets) {
    
    fname.suffix = fname.suffix.0%.%"_"%.%marker_set
    assays = subset(assay_metadata, panel==marker_set, assay, drop=T)
    if (len(assays)==0) {
      warning("no assays found")
      next
    }
    
    all.markers=c(paste0("Day", tpeak, assays), paste0("B", assays), paste0("Delta", tpeak, "overB", assays))
    tmp = get.short.name(assays); all.markers.names.short = c(glue("D{tpeak} {tmp}"), glue("D01 {tmp}"), glue("D{tpeak}/D01 {tmp}"))
    names(all.markers.names.short) = all.markers
    all.markers.names.long = 
      c(as.matrix(labels.title)[DayPrefix%.%tpeak, assays], 
        as.matrix(labels.title)["B", assays], 
        as.matrix(labels.title)["Delta"%.%tpeak%.%"overB", assays])
    names(all.markers.names.long) = all.markers
    
    # need to save tab1 for each distinct fname.suffix
    mytex(tab1, file.name = "tab1_" %.% fname.suffix, save2input.only = T, input.foldername = save.results.to)
  
    cat("\n\n"); myprint(fname.suffix)
    
    cor_coxph_coef_1(
      form.0,
      design_or_dat = design.1,
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
      cmp.label,
      verbose = F
    )
    
    # put six curves in the same plot
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
  
  } # end marker_set loop

} # end trt loop


# putting Cox tables from two arms in the same table
for (marker_set in marker_sets) {

  load(glue("{save.results.to}/coxph_slopes_InvVacc_{marker_set}.Rdata"))
  tab.cont.1 = tab.cont
  load(glue("{save.results.to}/coxph_slopes_CtlVacc_{marker_set}.Rdata"))
  tab.cont.2 = tab.cont

  header=paste0("\\hline\n 
         \\multicolumn{1}{l}{} & \\multicolumn{4}{c}{Investigational} & \\multicolumn{4}{c}{Comparator} \\\\
         \\multicolumn{1}{l}{} & \\multicolumn{1}{c}{No. cases /} & \\multicolumn{2}{c}{HR} & \\multicolumn{1}{c}{P-value}  &
                                 \\multicolumn{1}{c}{No. cases /} & \\multicolumn{2}{c}{HR} & \\multicolumn{1}{c}{P-value}  \\\\ 
         \\multicolumn{1}{l}{Immunologic Marker} & \\multicolumn{1}{c}{No. at-risk**} & \\multicolumn{1}{c}{Pt. Est.} & \\multicolumn{1}{c}{95\\% CI} & \\multicolumn{1}{c}{} &
                                                   \\multicolumn{1}{c}{No. at-risk**} & \\multicolumn{1}{c}{Pt. Est.} & \\multicolumn{1}{c}{95\\% CI} & \\multicolumn{1}{c}{} \\\\ 
         \\hline\n 
    ")
  
  fname.suffix = "2arms_"%.%marker_set
  
  mytex(cbind(tab.cont.1, tab.cont.2), 
        file.name="CoR_univariable_svycoxph_pretty_"%.%fname.suffix, align="c", include.colnames = F, save2input.only=T, input.foldername=save.results.to,
        col.headers=header,
        longtable=T, 
        label=paste0("tab:CoR_univariable_svycoxph_pretty_2arms"), 
        caption.placement = "top", 
        caption = escape_latex(glue("
Inferences for D01 and D31 antibody marker covariate-adjusted correlates of risk of COVID-19 through {tfinal.tpeak} days post D31 for the Investigational and Comparator Vaccine arms. \\
Non-cases have no evidence of SARS-CoV-2 infection (i.e., never tested nucleic acid amplification/PCR positive) after D01 up to the date by which the last enrolled participant reached {tfinal.tpeak} days post D31. \\
HRs were estimated using inverse probability sampling weighted Cox regression models; 95% confidence intervals (CIs) and Wald-based p-values are shown. \\
Analyses adjust for the randomization strata and baseline risk score via an inverse probability sampling weighted Cox model (cor_coxph module at CoVPN GitHub) (ccIAS-sera). \\
No. cases: number of COVID-19 endpoints eligible for correlates analysis regardless of whether antibody markers are available (phase 1/PPI set cases). \\
No. at-risk: number of participants eligible for correlates analysis regardless of whether antibody markers are available (Phase 1/PPI Set cases). \\
N: Nucleocapsid protein. ID50: 50% inhibitory serum dilution neutralizing antibody titer. nAb: neutralizing antibody.
\n"))
        # caption=paste0("Inference for Day ", config.cor$tpeak, " antibody marker covariate-adjusted correlates of risk of ", config.cor$txt.endpoint, 
                       # " in the ", escape(marker_set), " group: Hazard ratios. Baseline covariates adjusted for: ", escape(paste(deparse(form.0[[3]]), collapse = " ")), 
                       # ", endpoint variable: ", escape(config.cor$EventIndPrimary), ".")
  )
  
}


# putting trichotomized incidence curves from two arms on the same plot
for (marker_set in marker_sets) {

  assays = subset(assay_metadata, panel==marker_set, assay, drop=T)
  all.markers=c(paste0("Day", tpeak, assays), paste0("B", assays), paste0("Delta", tpeak, "overB", assays))
  tmp = get.short.name(assays); all.markers.names.short = c(glue("D{tpeak} {tmp}"), glue("D01 {tmp}"), glue("D{tpeak}/D01 {tmp}"))
  names(all.markers.names.short) = all.markers
  all.markers.names.long =
    c(as.matrix(labels.title)[DayPrefix%.%tpeak, assays],
      as.matrix(labels.title)["B", assays],
      as.matrix(labels.title)["Delta"%.%tpeak%.%"overB", assays])
  names(all.markers.names.long) = all.markers

  cor_coxph_risk_tertile_incidence_curves_2arms (
    form.0,
    dat = dat.vacc,
    fname.suffix = "InvVacc_"%.%marker_set,
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
    for.title = "",

    trt.label = trt.label,
    cmp.label = cmp.label,

    fname.suffix.2 = "CtlVacc_"%.%marker_set
  )

}



################################################################################
print(date())
print("cor_coxph run time: " %.% format(Sys.time() - time.start, digits = 1))
