COR="D43vat08_combined_M5_bAb"; Sys.setenv(stage = 2) 
# COR="D43vat08_combined_M6_st1.nAb.batch0and1"; Sys.setenv(stage = 1) 

Sys.setenv(TRIAL = "vat08_combined")
Sys.setenv(VERBOSE = 1) 
renv::activate(project = here::here(".."))     
source(here::here("..", "_common.R")) 
source(here::here("code", "params.R"))
source(here::here("code", "cor_coxph_coef_1_mi.R"))
source(here::here("code", "cor_coxph_coef_n_mi.R"))


{
library(kyotil) # p.adj.perm, getFormattedSummary
library(xtable) # this is a dependency of kyotil
library(marginalizedRisk)
library(tools) # toTitleCase
library(survey)
library(plotrix) # weighted.hist
library(parallel)
library(Hmisc) # wtd.quantile, cut2
library(mitools)
library(glue)

time.start=Sys.time()
print(date())
for (i in 1:5) cat(COR, " "); cat("\n")

tp = substr(COR,2,3) 
    
# path for figures and tables etc
save.results.to.0 = here::here("output"); if (!dir.exists(save.results.to.0))  dir.create(save.results.to.0)
save.results.to.0 = paste0(save.results.to.0, "/", attr(config,"config")); if (!dir.exists(save.results.to.0))  dir.create(save.results.to.0)
save.results.to.0 = paste0(save.results.to.0, "/", COR, "/"); if (!dir.exists(save.results.to.0))  dir.create(save.results.to.0)

# There are 3 analyses
# Stage 1 trial Non-naives: Omicron COVID-19 14 days post D22 (or post D43) through to 180 days post dose 2;
# Stage 2 trial Non-naives: Same as above.
# save results in separate folders
save.results.to = paste0(save.results.to.0, "/stage1nnaive/"); if (!dir.exists(save.results.to))  dir.create(save.results.to)
save.results.to = paste0(save.results.to.0, "/stage2nnaive/"); if (!dir.exists(save.results.to))  dir.create(save.results.to)

# B=1e3 and numPerm=1e4 take 10 min to run with 30 CPUS for one analysis
B <-       config$num_boot_replicates 
numPerm <- config$num_perm_replicates # number permutation replicates 1e4
myprint(B, numPerm)


# toggle between nnaive and prev_inf
# dat_proc$kp = dat_proc$prev_inf
dat_proc$kp = dat_proc$Bserostatus # nnaive

dat.vac.seropos.st1 = subset(dat_proc, Trt==1 & kp==1 & Trialstage==1 & ph1)
dat.pla.seropos.st1 = subset(dat_proc, Trt==0 & kp==1 & Trialstage==1 & ph1)
dat.vac.seropos.st2 = subset(dat_proc, Trt==1 & kp==1 & Trialstage==2 & ph1)
dat.pla.seropos.st2 = subset(dat_proc, Trt==0 & kp==1 & Trialstage==2 & ph1)

for (a in c("Day"%.%tpeak%.%assays, "B"%.%assays, "Delta"%.%tpeak%.%"overB"%.%assays)) {
  dat.vac.seropos.st1[[a%.%"centered"]] = scale(dat.vac.seropos.st1[[a]], scale=F)
  dat.pla.seropos.st1[[a%.%"centered"]] = scale(dat.pla.seropos.st1[[a]], scale=F)
  dat.vac.seropos.st2[[a%.%"centered"]] = scale(dat.vac.seropos.st2[[a]], scale=F)
  dat.pla.seropos.st2[[a%.%"centered"]] = scale(dat.pla.seropos.st2[[a]], scale=F)
  if(endsWith(COR, "st1.nAb.batch0and1")) {
    for (imp in 1:10){
      dat.vac.seropos.st1[[a%.%"centered_"%.%imp]] = scale(dat.vac.seropos.st1[[a%.%"_"%.%imp]], scale=F)
      dat.pla.seropos.st1[[a%.%"centered_"%.%imp]] = scale(dat.pla.seropos.st1[[a%.%"_"%.%imp]], scale=F)
    }
  }
}

}

# loop through stage 1 and 2 non-naive
# for st2 sensitivity analysis, only do stage 2
# forgo the naive populations from mono- and bi-valent trials

if (Sys.getenv("stage")=="") {
  stages=1:2
  if(endsWith(COR, "st1.nAb.batch0and1") | contain(COR, "M5")) {
    stages=1
  } else if(endsWith(COR, "st2.nAb.sen")) {
    stages=2
  }
} else {
  stages=as.numeric(Sys.getenv("stage"))
}
myprint(stages)




################################################################################

for (iSt in stages) {
  # iSt=2
  
  cat("\n\n\n\n")
  myprint(iSt)
  
  if (iSt==1) {dat.vacc=dat.vac.seropos.st1; dat.plac=dat.pla.seropos.st1; save.results.to=save.results.to.0%.%"stage1nnaive/"}
  if (iSt==2) {dat.vacc=dat.vac.seropos.st2; dat.plac=dat.pla.seropos.st2; save.results.to=save.results.to.0%.%"stage2nnaive/"}
  # if (iSt==3) {dat.vacc=dat.vac.seroneg.st2; dat.plac=dat.pla.seroneg.st2; save.results.to=save.results.to.0%.%"stage2naive/"}
  
  form.0 = update(Surv(EventTimeOfInterest, EventIndOfInterest) ~ 1, as.formula(config$covariates))
  
  # get cutpoints and dichcutpoints, and turn discrete markers into factors
  {
    
  # need to do it within iSt loop because marker.cutpoints are needed in later function calls
    
  # vaccine, trichotomized
  marker.cutpoints = list()
  for (a in c(paste0("Day", tpeak, assays),
              paste0("B", assays),
              paste0("Delta", tpeak, "overB", assays))) {
    # get cut points
    tmpname = names(table(dat.vacc[[a%.%"cat"]]))[2]; tmpname
    stopifnot(!is.na(tmpname))
    tmpname = substr(tmpname, 2, nchar(tmpname)-1)
    tmpname = as.numeric(strsplit(tmpname, ",")[[1]])
    tmpname = setdiff(tmpname,Inf) # if there are two categories, remove the second cut point, which is Inf
    marker.cutpoints[[a]] <- tmpname
    
    dat.vacc[[a%.%"cat"]] = as.factor(dat.vacc[[a%.%"cat"]])

    write(paste0(escape(a),     " [", concatList(round(marker.cutpoints[[a]], 2), ", "), ")%"), 
          file=paste0(save.results.to, "cutpoints_", a, ".txt"))
  }
  

  # vaccine, dichotomized
  marker.cutpoints.dich = list()
  for (a in "B"%.%assays) {
    # get cut points
    levels = names(table(dat.vacc[[a%.%"dich"]]))
    tmpname = levels[2]
    tmpname = substr(tmpname, 2, nchar(tmpname)-1)
    tmpname = as.numeric(strsplit(tmpname, ",")[[1]])
    tmpname = setdiff(tmpname,Inf) # the second cut point  is Inf
    marker.cutpoints.dich[[a]] <- tmpname
    
    # define two factor variables with different reference levels
    dat.vacc[[a%.%"dich1"]] = factor(dat.vacc[[a%.%"dich"]], levels=rev(levels))
    dat.vacc[[a%.%"dich"]]  = factor(dat.vacc[[a%.%"dich"]], levels=levels)
      
    write(paste0(escape(a), " ", concatList(round(marker.cutpoints.dich[[a]], 2), ", "), "%"), 
          file=paste0(save.results.to, "dichcutpoints_", a, ".txt"))
  }
  
    
  # placebo, trichotomized
  marker.cutpoints.plac = list()
  for (a in c(paste0("Day", tpeak, assays))) {
    # get cut points
    tmpname = names(table(dat.plac[[a%.%"cat"]]))[2]; tmpname
    stopifnot(!is.na(tmpname))
    tmpname = substr(tmpname, 2, nchar(tmpname)-1)
    tmpname = as.numeric(strsplit(tmpname, ",")[[1]])
    tmpname = setdiff(tmpname,Inf) # if there are two categories, remove the second cut point, which is Inf
    marker.cutpoints.plac[[a]] <- tmpname
    
    dat.plac[[a%.%"cat"]] = as.factor(dat.plac[[a%.%"cat"]])
    
    write(paste0(escape(a),     " [", concatList(round(marker.cutpoints.plac[[a]], 2), ", "), ")%"), 
          file=paste0(save.results.to, "cutpoints_", a, "_plac.txt"))
  }
  
  } # end marker cut points
  
  
  
  # univariate
  if(T){
  ###################################
  # Univariate: Dxx, B, Dxx/B, vaccine arm
  
  all.markers=c(paste0("Day", tpeak, assays),
                paste0("B", assays),
                paste0("Delta", tpeak, "overB", assays))
  
  all.markers.names.short = sub("Pseudovirus-", "", assay_metadata$assay_label_short[match(assays,assay_metadata$assay)])
  all.markers.names.short = sub(" \\(AU/ml\\)", "", sub("Anti Spike ", "", all.markers.names.short))
  all.markers.names.short = c("D"%.%tpeak%.%" "%.%all.markers.names.short,
                              "B "%.%all.markers.names.short,
                              "D"%.%tpeak%.%"/B "%.%all.markers.names.short)
  names(all.markers.names.short) = all.markers
  
  multivariate_assays = config$multivariate_assays
  
  
  cor_coxph_coef_1_mi (
    form.0,
    dat=dat.vacc,
    fname.suffix="D"%.%tpeak,
    save.results.to,
    config,
    config.cor,
    
    markers=all.markers,
    markers.names.short=all.markers.names.short,
    
    dat.pla.seroneg = dat.plac,
    show.q=FALSE,
    
    forestplot.markers=list(1:length(assays), 1:length(assays)+length(assays), 1:length(assays)+2*length(assays)),
    for.title=paste0("Stage ",iSt," NN, Vaccine"),
    
    verbose=T
  )
  
  
  
  ###################################
  # Univariate trichotomized curves: D, D/B
  
  all.markers=c(paste0("Day", tpeak, assays),
                paste0("B", assays),
                paste0("Delta", tpeak, "overB", assays))
  
  all.markers.names.short = sub("Pseudovirus-", "", assay_metadata$assay_label_short[match(assays,assay_metadata$assay)])
  all.markers.names.short = sub(" \\(AU/ml\\)", "", sub("Anti Spike ", "", all.markers.names.short))
  all.markers.names.short = c("D"%.%tpeak%.%" "%.%all.markers.names.short,
                              "B"%.%" "%.%all.markers.names.short,
                              "D"%.%tpeak%.%"/B "%.%all.markers.names.short)
  names(all.markers.names.short) = all.markers
  
  all.markers.names.long  = sub("Pseudovirus-", "", assay_metadata$assay_label[match(assays,assay_metadata$assay)])
  all.markers.names.long  = c("D"%.%tpeak%.%" "%.%all.markers.names.long,
                              "B"%.%" "%.%all.markers.names.long,
                              "D"%.%tpeak%.%"/B "%.%all.markers.names.long)
  names(all.markers.names.long) = all.markers
  
  cor_coxph_risk_tertile_incidence_curves(
    form.0,
    dat=dat.vacc,
    fname.suffix="D"%.%tpeak,
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
    for.title=paste0("Stage ",iSt," NN, Vaccine")
  )
  
  
  
  ###################################
  # Univariate Placebo: Dxx
  
  all.markers=c(paste0("Day", tpeak, assays))
  
  all.markers.names.short = sub("Pseudovirus-", "", assay_metadata$assay_label_short[match(assays,assay_metadata$assay)])
  all.markers.names.short = sub(" \\(AU/ml\\)", "", sub("Anti Spike ", "", all.markers.names.short))
  all.markers.names.short = c("D"%.%tpeak%.%" "%.%all.markers.names.short)
  names(all.markers.names.short) = all.markers
  multivariate_assays = config$multivariate_assays
  
  cor_coxph_coef_1_mi (
    form.0,
    dat=dat.plac,
    fname.suffix="D"%.%tpeak%.%"_plac",
    save.results.to,
    config,
    config.cor,
    markers=all.markers,
    markers.names.short=all.markers.names.short,
    
    dat.pla.seroneg = NULL,
    show.q=FALSE,
    
    forestplot.markers=1:length(assays),
    for.title=paste0("Stage ",iSt," NN, Placebo"),
    
    verbose=T
  )
  
  ## comment out since fails to run with prev_inf due to zero sampling fraction
  # # repeat stage 1 placebo, keeping only the countries that also appear in stage 2
  # if(iSt==1) {
  #   cor_coxph_coef_1_mi (
  #     form.0,
  #     dat=subset(dat.plac, cc %in% c("Columbia", "Ghana", "Kenya", "Nepal", "India")),
  #     fname.suffix="D"%.%tpeak%.%"_plac_alt2",
  #     save.results.to,
  #     config,
  #     config.cor,
  #     markers=all.markers,
  #     markers.names.short=all.markers.names.short,
  #     
  #     dat.pla.seroneg = NULL,
  #     show.q=FALSE,
  #     
  #     forestplot.markers=1:length(assays),
  #     for.title=paste0("Stage ",iSt," NN, Placebo"),
  #     
  #     verbose=T
  #   )
  #   
  # }
  
  
  ###################################
  # Univariate trichotomized curves, placebo: D
  
  all.markers=c(paste0("Day", tpeak, assays))
  
  all.markers.names.short = sub("Pseudovirus-", "", assay_metadata$assay_label_short[match(assays,assay_metadata$assay)])
  all.markers.names.short = sub(" \\(AU/ml\\)", "", sub("Anti Spike ", "", all.markers.names.short))
  all.markers.names.short = c("D"%.%tpeak%.%" "%.%all.markers.names.short)
  names(all.markers.names.short) = all.markers
  
  all.markers.names.long  = sub("Pseudovirus-", "", assay_metadata$assay_label[match(assays,assay_metadata$assay)])
  all.markers.names.short = sub(" \\(AU/ml\\)", "", sub("Anti Spike ", "", all.markers.names.short))
  all.markers.names.long  = c("D"%.%tpeak%.%" "%.%all.markers.names.long)
  names(all.markers.names.long) = all.markers
  
  cor_coxph_risk_tertile_incidence_curves(
    form.0,
    dat=dat.plac,
    fname.suffix="D"%.%tpeak%.%"_plac",
    save.results.to,
    config,
    config.cor,
    tfinal.tpeak,
    
    markers = all.markers,
    markers.names.short = all.markers.names.short,
    markers.names.long = all.markers.names.long,
    marker.cutpoints.plac,
    assay_metadata,
    
    dat.vacc,
    for.title=paste0("Stage ",iSt," NN, Placebo"),
    plac.actually=T
  )
  }
  
  
  # baseline detectable * Dxx
  if(F) {
    
  ###################################
  # Bhigh x D15 (D15/B), Vacc
  
  all.markers = c(sapply(assays, function (a) paste0("Bhigh * Day"%.%tpeak%.%"",a,"centered")),
                  sapply(assays, function (a) paste0("Bhigh * Delta"%.%tpeak%.%"overB",a,"centered"))
  )
  all.markers.names.short = sub("Pseudovirus-", "", assay_metadata$assay_label_short[match(assays,assay_metadata$assay)])
  all.markers.names.short = sub(" \\(AU/ml\\)", "", sub("Anti Spike ", "", all.markers.names.short))
  all.markers.names.short = c("D"%.%tpeak%.%" "%.%all.markers.names.short, "D"%.%tpeak%.%"/B "%.%all.markers.names.short)
  names(all.markers.names.short) = all.markers
  
  # parameters for R script
  nCoef=3
  col.headers=c("Bhigh", "center(D"%.%tpeak%.%" or fold)", "Bhigh:center(D"%.%tpeak%.%" or fold)")
  
  cor_coxph_coef_n_mi (
    form.0,
    dat=dat.vacc,
    fname.suffix="Bhigh*D"%.%tpeak,
    save.results.to,
    config,
    config.cor,
    all.markers,
    all.markers.names.short,
    
    nCoef,
    col.headers,
    verbose=verbose
  )
  
  
  ###################################
  # (1 - Bhigh) x D15 (D15/B), Vacc
  
  all.markers = c(sapply(assays, function (a) paste0("I(1-Bhigh) * Day"%.%tpeak%.%"",a,"centered")),
                  sapply(assays, function (a) paste0("I(1-Bhigh) * Delta"%.%tpeak%.%"overB",a,"centered"))
  )
  all.markers.names.short = sub("Pseudovirus-", "", assay_metadata$assay_label_short[match(assays,assay_metadata$assay)])
  all.markers.names.short = sub(" \\(AU/ml\\)", "", sub("Anti Spike ", "", all.markers.names.short))
  all.markers.names.short = c("D"%.%tpeak%.%" "%.%all.markers.names.short, "D"%.%tpeak%.%"/B "%.%all.markers.names.short)
  names(all.markers.names.short) = all.markers
  
  # parameters for R script
  nCoef=3
  col.headers=c("Blow", "center(D"%.%tpeak%.%" or fold)", "Blow:center(D"%.%tpeak%.%" or fold)")
  
  cor_coxph_coef_n_mi (
    form.0,
    dat=dat.vacc,
    fname.suffix="(1-Bhigh)*D"%.%tpeak,
    save.results.to,
    config,
    config.cor,
    all.markers,
    all.markers.names.short,
    
    nCoef,
    col.headers,
    verbose=verbose
  )
  }
  
  
  # Bxx + Dxx, Bxx * Dxx
  if(T) { # for st1.nAb.batch0and1, cor_coxph_coef_n_mi needs to be modified with respect to f
    
  ###################################
  # B + Dxx (Dxx/B), vaccine arm
  
  all.markers = c(sapply(assays, function (a) paste0("B",a, "centered + Day"%.%tpeak%.%"",a, "centered")),
                  sapply(assays, function (a) paste0("B",a, "centered + Delta"%.%tpeak%.%"overB",a, "centered"))
  )
  
  all.markers.names.short = sub("Pseudovirus-", "", assay_metadata$assay_label_short[match(assays,assay_metadata$assay)])
  all.markers.names.short = sub(" \\(AU/ml\\)", "", sub("Anti Spike ", "", all.markers.names.short))
  all.markers.names.short = c("D"%.%tpeak%.%" "%.%all.markers.names.short, "D"%.%tpeak%.%"/B "%.%all.markers.names.short)
  names(all.markers.names.short) = all.markers
  nCoef=2
  col.headers=c("center(B)", "center(D"%.%tpeak%.%" or fold)")
  
  cor_coxph_coef_n_mi (
    form.0,
    dat=dat.vacc,
    fname.suffix="B+D"%.%tpeak,
    save.results.to,
    config,
    config.cor,
    all.markers,
    all.markers.names.short,

    nCoef,
    col.headers,
    verbose=verbose
  )


  ###################################
  # B x D15 (D15/B)
  
  all.markers = c(sapply(assays, function (a) paste0("B",a, "centered * Day"%.%tpeak%.%"",a,"centered")),
                  sapply(assays, function (a) paste0("B",a, "centered * Delta"%.%tpeak%.%"overB",a,"centered"))
  )
    
  all.markers.names.short = sub("Pseudovirus-", "", assay_metadata$assay_label_short[match(assays,assay_metadata$assay)])
  all.markers.names.short = sub(" \\(AU/ml\\)", "", sub("Anti Spike ", "", all.markers.names.short))
  all.markers.names.short = c("D"%.%tpeak%.%" "%.%all.markers.names.short, "D"%.%tpeak%.%"/B "%.%all.markers.names.short)
  names(all.markers.names.short) = all.markers
  
  # parameters for R script
  nCoef=3
  col.headers=c("center(B)", "center(D"%.%tpeak%.%" or fold)", "center(B):center(D"%.%tpeak%.%" or fold)")
  
  # vaccine arm
  
  cor_coxph_coef_n_mi (
    form.0,
    dat=dat.vacc,
    fname.suffix="B*D"%.%tpeak,
    save.results.to,
    config,
    config.cor,
    all.markers,
    all.markers.names.short,

    nCoef,
    col.headers,
    verbose=verbose
  )
  
  }
  
  
  # dich_B x 
  if(F) {

  ###################################
  # dich_B x D15 (D15/B)
  
  all.markers = c(sapply(assays, function (a) paste0("B",a, "dich * Day"%.%tpeak%.%"",a,"centered")),
                  sapply(assays, function (a) paste0("B",a, "dich * Delta"%.%tpeak%.%"overB",a,"centered"))
  )
  
  all.markers.names.short = sub("Pseudovirus-", "", assay_metadata$assay_label_short[match(assays,assay_metadata$assay)])
  all.markers.names.short = sub(" \\(AU/ml\\)", "", sub("Anti Spike ", "", all.markers.names.short))
  all.markers.names.short = c("D"%.%tpeak%.%" "%.%all.markers.names.short, "D"%.%tpeak%.%"/B "%.%all.markers.names.short)
  names(all.markers.names.short) = all.markers
  
  # parameters for R script
  nCoef=3
  col.headers=c("dich(B)", "center(D"%.%tpeak%.%" or fold)", "dich(B):center(D"%.%tpeak%.%" or fold)")
  
  # vaccine arm
  
  cor_coxph_coef_n_mi (
    form.0,
    dat=dat.vacc,
    fname.suffix="dichB*D"%.%tpeak,
    save.results.to,
    config,
    config.cor,
    all.markers,
    all.markers.names.short,
    
    nCoef,
    col.headers,
    verbose=verbose
  )
  
  
  
  ###################################
  # (1 - dich_B) x D15 (D15/B)
  
  all.markers = c(sapply(assays, function (a) paste0("B",a, "dich1 * Day"%.%tpeak%.%"",a,"centered")),
                  sapply(assays, function (a) paste0("B",a, "dich1 * Delta"%.%tpeak%.%"overB",a,"centered"))
  )
  
  all.markers.names.short = sub("Pseudovirus-", "", assay_metadata$assay_label_short[match(assays,assay_metadata$assay)])
  all.markers.names.short = sub(" \\(AU/ml\\)", "", sub("Anti Spike ", "", all.markers.names.short))
  all.markers.names.short = c("D"%.%tpeak%.%" "%.%all.markers.names.short, "D"%.%tpeak%.%"/B "%.%all.markers.names.short)
  names(all.markers.names.short) = all.markers
  
  # parameters for R script
  nCoef=3
  col.headers=c("1-dich(B)", "center(D"%.%tpeak%.%" or fold)", "(1-dich(B)):center(D"%.%tpeak%.%" or fold)")
  
  # vaccine arm
  
  cor_coxph_coef_n_mi (
    form.0,
    dat=dat.vacc,
    fname.suffix="(1-dichB)*D"%.%tpeak,
    save.results.to,
    config,
    config.cor,
    all.markers,
    all.markers.names.short,
    
    nCoef,
    col.headers,
    verbose=verbose
  )
  
  }
  

  
  ############################
  # count ph1 and ph2 cases
  
  # vaccine arm
  fname.suffix = "D"%.%tpeak
  
  # imputed events of interest
  tabs=sapply(1:10, simplify="array", function (imp) {
    dat.vacc$EventIndOfInterest = dat.vacc[[config.cor$EventIndPrimary%.%imp]]
    with(dat.vacc, table(ph2, EventIndOfInterest))
  }) 
  tab =apply(tabs, c(1,2), mean)
  names(dimnames(tab))[2]="Event Indicator"
  tab
  mytex(tab,     
        # for some unknown reason, kyoti::mytex line 160 runs into file not found error if file.name is not full path.
        file.name = save.results.to %.% "tab1_" %.% fname.suffix, 
        save2input.only=T, 
        input.foldername=save.results.to, 
        digits=1)
  
  # placebo arm
  fname.suffix = "D"%.%tpeak%.%"_plac"
  
  # imputed events of interest
  tabs=sapply(1:10, simplify="array", function (imp) {
    dat.plac$EventIndOfInterest = dat.plac[[config.cor$EventIndPrimary%.%imp]]
    with(dat.plac, table(ph2, EventIndOfInterest))
  }) 
  tab =apply(tabs, c(1,2), mean)
  names(dimnames(tab))[2]="Event Indicator"
  tab
  mytex(tab,     
        # for some unknown reason, kyoti::mytex line 160 runs into file not found error if file.name is not full path.
        file.name = save.results.to %.% "tab1_" %.% fname.suffix, 
        save2input.only=T, 
        input.foldername=save.results.to, 
        digits=1)
  
  
  if(iSt==1) {
    # placebo arm stage 1 with Columbia, Ghana, Kenya, Nepal, India
    fname.suffix = "D"%.%tpeak%.%"_plac_alt2"
    
    # imputed events of interest
    tabs=sapply(1:10, simplify="array", function (imp) {
      dat.plac$EventIndOfInterest = dat.plac[[config.cor$EventIndPrimary%.%imp]]
      with(dat.plac[dat.plac$cc %in% c("Columbia", "Ghana", "Kenya", "Nepal", "India"),], table(ph2, EventIndOfInterest))
    }) 
    tab =apply(tabs, c(1,2), mean)
    names(dimnames(tab))[2]="Event Indicator"
    tab
    mytex(tab,     
          file.name = save.results.to %.% "tab1_" %.% fname.suffix, 
          save2input.only=T, 
          input.foldername=save.results.to, 
          digits=1)
    

  }  
    
  
  
} # iSt loop


print(date())

print("cor_coxph run time: " %.% format(Sys.time()-time.start, digits=1) )
