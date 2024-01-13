# COR="D15to181"

Sys.setenv(TRIAL = "covail")
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
save.results.to = here::here("output"); if (!dir.exists(save.results.to))  dir.create(save.results.to)
save.results.to = paste0(save.results.to, "/", attr(config,"config")); if (!dir.exists(save.results.to))  dir.create(save.results.to)
save.results.to = paste0(save.results.to, "/", COR, "/"); if (!dir.exists(save.results.to))  dir.create(save.results.to)


# B=1e3 and numPerm=1e4 take 10 min to run with 30 CPUS for one analysis
B <-       config$num_boot_replicates 
numPerm <- config$num_perm_replicates # number permutation replicates 1e4
myprint(B, numPerm)


# redefine assays to focus on the 12 markers
assays = c("pseudoneutid50_D614G", "pseudoneutid50_Delta", "pseudoneutid50_Beta", "pseudoneutid50_BA.1", "pseudoneutid50_BA.4.BA.5", "pseudoneutid50_MDW")
all.markers = c("B"%.%assays, "Day15"%.%assays)

# add trichotomized markers 
dat = subset(dat.mock, TrtmRNA==1 & arm!=3 & ph1.D15) # TrtonedosemRNA := TrtmRNA==1 & arm!=3
dat$ph2.D15=T
dat = add.trichotomized.markers (dat, all.markers, ph2.col.name="ph2.D15", wt.col.name="wt.D15")

# save cut points to files
for (a in all.markers) {        
  write(paste0(gsub("_", "\\_", a, fixed = TRUE),     " [", concatList(round(marker.cutpoints[[a]], 2), ", "), ")%"), 
        file=paste0(save.results.to, "cutpoints_", a, "_"%.%study_name))
}


################################################################################
# loop through 3 analyses

for (iAna in 3:3) {
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
  
  if (iAna==3) dat.vac=subset(dat.vac, !(Country %in% c(10 ) & Trialstage==2 & Bserostatus==1 & Age<60))

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
