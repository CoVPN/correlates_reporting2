# COR="D22M6ominAb"
# COR="D43M12ominAb"
Sys.setenv(TRIAL = "vat08_combined")

Sys.setenv(VERBOSE = 1) 
renv::activate(project = here::here(".."))     
source(here::here("..", "_common.R")) 
source(here::here("code", "params.R"))


library(kyotil) # p.adj.perm, getFormattedSummary
library(xtable) # this is a dependency of kyotil
library(marginalizedRisk)
library(tools) # toTitleCase
library(survey)
library(plotrix) # weighted.hist
library(parallel)
library(forestplot)
library(Hmisc) # wtd.quantile, cut2
library(mitools)
library(glue)

time.start=Sys.time()
print(date())


# path for figures and tables etc
save.results.to.0 = here::here("output"); if (!dir.exists(save.results.to.0))  dir.create(save.results.to.0)
save.results.to.0 = paste0(save.results.to.0, "/", attr(config,"config")); if (!dir.exists(save.results.to.0))  dir.create(save.results.to.0)
save.results.to.0 = paste0(save.results.to.0, "/", COR, "/"); if (!dir.exists(save.results.to.0))  dir.create(save.results.to.0)

# There are 3 analyses
# Stage 1 trial Non-naives: Omicron COVID-19 14 days post D22 (or post D43) through to 180 days post dose 2;
# Stage 2 trial Naives: Same as above;
# Stage 2 trial Non-naives: Same as above.
# save results in separate folders
save.results.to = paste0(save.results.to.0, "/stage1nnaive/"); if (!dir.exists(save.results.to))  dir.create(save.results.to)
save.results.to = paste0(save.results.to.0, "/stage2naive/");  if (!dir.exists(save.results.to))  dir.create(save.results.to)
save.results.to = paste0(save.results.to.0, "/stage2nnaive/"); if (!dir.exists(save.results.to))  dir.create(save.results.to)


# B=1e3 and numPerm=1e4 take 10 min to run with 30 CPUS for one analysis
B <-       config$num_boot_replicates 
numPerm <- config$num_perm_replicates # number permutation replicates 1e4
myprint(B, numPerm)


## add trichotomized markers 

summary(subset(dat.mock, Trt==1 & Trialstage==1 & Bserostatus==0, Day22pseudoneutid50_BA.1, drop=T))
summary(subset(dat.mock, Trt==1 & Trialstage==1 & Bserostatus==1, Day22pseudoneutid50_BA.1, drop=T))
summary(subset(dat.mock, Trt==1 & Trialstage==2 & Bserostatus==0, Day22pseudoneutid50_BA.1, drop=T))
summary(subset(dat.mock, Trt==1 & Trialstage==2 & Bserostatus==1, Day22pseudoneutid50_BA.1, drop=T))

summary(subset(dat.mock, Trt==1 & Trialstage==1 & Bserostatus==0, Day43pseudoneutid50_BA.1, drop=T))
summary(subset(dat.mock, Trt==1 & Trialstage==1 & Bserostatus==1, Day43pseudoneutid50_BA.1, drop=T))
summary(subset(dat.mock, Trt==1 & Trialstage==2 & Bserostatus==0, Day43pseudoneutid50_BA.1, drop=T))
summary(subset(dat.mock, Trt==1 & Trialstage==2 & Bserostatus==1, Day43pseudoneutid50_BA.1, drop=T))

# based on the above descriptives, will use the same cut points for nnaive stage 1 and stage 2
# and a separate set of cut points for naive stage 2

dat.vac.seropos = subset(dat.mock, Trt==1 & ph1 & Bserostatus==1)
dat.vac.seropos = add.trichotomized.markers (dat.vac.seropos, all.markers, wt.col.name="wt")

dat.vac.seroneg.2 = subset(dat.mock, Trt==1 & ph1 & Trialstage==2 & Bserostatus==0)
dat.vac.seroneg.2 = add.trichotomized.markers (dat.vac.seroneg.2, all.markers, wt.col.name="wt")

if (COR=="D22M6omi" | COR=="D22M6ominAb") {
  # stage 2, naive, Omicron ID50 markers are 90% below limit, change to cut into to low and high
  for (a in all.markers[startsWith(all.markers, 'Day22pseudoneutid50')]) {        
    cutpoint=min(dat.vac.seroneg.2[[a]], na.rm = T)
    dat.vac.seroneg.2[[a%.%'cat']] = cut(dat.vac.seroneg.2[[a]], breaks = c(-Inf, cutpoint, Inf))
    attr(dat.vac.seroneg.2, "marker.cutpoints")[[a]] = cutpoint
  }
}

# save cut points to files
for (a in all.markers) {        
  marker.cutpoints=attr(dat.vac.seropos, "marker.cutpoints")[[a]]
  write(paste0(gsub("_", "\\_", a, fixed = TRUE),     " [", concatList(round(marker.cutpoints, 2), ", "), ")%"), 
        file=paste0(save.results.to.0%.%"/stage1nnaive/", "cutpoints_", a, "_"%.%study_name))
  write(paste0(gsub("_", "\\_", a, fixed = TRUE),     " [", concatList(round(marker.cutpoints, 2), ", "), ")%"), 
        file=paste0(save.results.to.0%.%"/stage2nnaive/", "cutpoints_", a, "_"%.%study_name))
  
  marker.cutpoints=attr(dat.vac.seroneg.2, "marker.cutpoints")[[a]]
  write(paste0(gsub("_", "\\_", a, fixed = TRUE),     " [", concatList(round(marker.cutpoints, 2), ", "), ")%"), 
        file=paste0(save.results.to.0%.%"/stage2naive/", "cutpoints_", a, "_"%.%study_name))
}

# separate nnaive by trial stage
dat.vac.seropos.1 = subset(dat.vac.seropos, Trialstage==1)
dat.vac.seropos.2 = subset(dat.vac.seropos, Trialstage==2)

# add placebo counterpart
dat.pla.seropos.1=subset(dat.mock, Trt==0 & ph1 & Trialstage==1 & Bserostatus==1)
dat.pla.seroneg.2=subset(dat.mock, Trt==0 & ph1 & Trialstage==2 & Bserostatus==0)
dat.pla.seropos.2=subset(dat.mock, Trt==0 & ph1 & Trialstage==2 & Bserostatus==1)

# for validation use
rv=list() 



################################################################################
# loop through 3 analyses

for (iAna in 1:3) {
  # iAna=3
  cat("\n\n\n\n")
  myprint(iAna)
  
  if (iAna==1) {dat.vac=dat.vac.seropos.1; dat.pla=dat.pla.seropos.1; save.results.to=save.results.to.0%.%"stage1nnaive/"}
  if (iAna==2) {dat.vac=dat.vac.seroneg.2; dat.pla=dat.pla.seroneg.2; save.results.to=save.results.to.0%.%"stage2naive/"}
  if (iAna==3) {dat.vac=dat.vac.seropos.2; dat.pla=dat.pla.seropos.2; save.results.to=save.results.to.0%.%"stage2nnaive/"}

  fname.suffix = study_name
  
  ############################
  # count ph1 and ph2 cases
  
  # imputed events of interest
  tabs=sapply(1:10, simplify="array", function (imp) {
    dat.vac$EventIndOfInterest = dat.vac[[config.cor$EventIndPrimary%.%imp]]
    with(dat.vac, table(ph2, EventIndOfInterest))
  })
  tab =apply(tabs, c(1,2), mean)
  names(dimnames(tab))[2]="Event Indicator"
  tab
  mytex(tab, file.name="tab1", save2input.only=T, input.foldername=save.results.to, digits=1)
  
  
  ############################
  # formula for coxph

  form.0 = update(Surv(EventTimeOfInterest, EventIndOfInterest) ~ 1, as.formula(config$covariates_riskscore))

  source(here::here("code", "cor_coxph_ph_MI.R"))
  
  
  
  # #####################################
  # # formula for competing risk analysis
  # 
  # # if there are very few competing events, the coxph for competing event may throw warnings
  # 
  # form.0=list(
  #   update(Surv(EventTimePrimaryD29, EventIndOfInterest) ~ 1, as.formula(config$covariates_riskscore)),
  #   update(Surv(EventTimePrimaryD29, EventIndCompeting)  ~ 1, as.formula(config$covariates_riskscore))
  # )
  # 
  # tfinal.tpeak = tfinal.tpeak.ls[[region]][[variant]]
  # write(tfinal.tpeak, file=paste0(save.results.to, "timepoints_cum_risk_", fname.suffix))
  # 
  # # run analyses
  # source(here::here("code", "cor_coxph_risk_no_marker.R"))
  # source(here::here("code", "cor_coxph_risk_bootstrap.R"))
  # 
  # # make tables and figures
  # source(here::here("code", "cor_coxph_risk_plotting.R"))
  # 

} # for iRegion





print(date())
print("cor_coxph run time: "%.%format(Sys.time()-time.start, digits=1))
