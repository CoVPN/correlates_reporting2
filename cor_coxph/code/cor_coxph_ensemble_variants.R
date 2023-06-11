#Sys.setenv(TRIAL = "janssen_partA_VL"); COR="D29VL"; Sys.setenv(VERBOSE = 1) 
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

begin=Sys.time()
print(date())

# with(subset(dat.mock, EventIndPrimaryIncludeNotMolecConfirmedD29 & Trt==1 & ph1.D29), table(is.na(seq1.spike.weighted.hamming.hotdeck1), is.na(seq1.log10vl), EventIndPrimaryMolecConfirmedD29))
# with(subset(dat.mock, EventIndPrimaryIncludeNotMolecConfirmedD29 & Trt==1 & ph1.D29 & is.na(seq1.spike.weighted.hamming.hotdeck1) & EventIndPrimaryMolecConfirmedD29==0), summary(EventTimePrimary))
# with(subset(dat.mock, EventIndPrimaryIncludeNotMolecConfirmedD29 & Trt==1 & ph1.D29 & is.na(seq1.spike.weighted.hamming.hotdeck1) & EventIndPrimaryMolecConfirmedD29==1), summary(EventTimePrimary))
# with(subset(dat.mock, EventIndPrimaryIncludeNotMolecConfirmedD29 & Trt==1 & ph1.D29 & is.na(seq1.spike.weighted.hamming.hotdeck1) & EventIndPrimaryMolecConfirmedD29==0), summary(EventTimePrimary+CalendarDateEnrollment))
# with(subset(dat.mock, EventIndPrimaryIncludeNotMolecConfirmedD29 & Trt==1 & ph1.D29 & is.na(seq1.spike.weighted.hamming.hotdeck1) & EventIndPrimaryMolecConfirmedD29==1), summary(EventTimePrimary+CalendarDateEnrollment))
# subset(dat.mock, EventIndPrimaryIncludeNotMolecConfirmedD29 & Trt==1 & ph1.D29 & is.na(seq1.spike.weighted.hamming.hotdeck1) & EventIndPrimaryMolecConfirmedD29==0, select=Ptid, drop=T)


# path for figures and tables etc
save.results.to = here::here("output"); if (!dir.exists(save.results.to))  dir.create(save.results.to)
save.results.to = paste0(save.results.to, "/", attr(config,"config")); if (!dir.exists(save.results.to))  dir.create(save.results.to)
save.results.to = paste0(save.results.to, "/", COR, "/"); if (!dir.exists(save.results.to))  dir.create(save.results.to)

# B=1e3 and numPerm=1e4 take 10 min to run with 30 CPUS for one analysis
B <-       config$num_boot_replicates 
numPerm <- config$num_perm_replicates # number permutation replicates 1e4
myprint(B, numPerm)


# uloq censoring, done here b/c should not be done for immunogenicity reports
# note that if delta are used, delta needs to be recomputed
for (a in assays) {
  uloq=assay_metadata$uloq[assay_metadata$assay==a]
  for (t in c("Day29")  ) {
    dat.mock[[t %.% a]] <- ifelse(dat.mock[[t %.% a]] > log10(uloq), log10(uloq), dat.mock[[t %.% a]])
  }
}    

# "Ancestral.Lineage", "Alpha", "Beta", "Delta", "Epsilon", "Gamma", "Lambda", "Mu", "Zeta", "Iota"
variants=lapply(tfinal.tpeak.ls, function(x) names(x))

# add trichotomized markers using the same cutoffs for all regions
dat.vac.seroneg.allregions=subset(dat.mock, Trt==1 & ph1)
dat.vac.seroneg.allregions = add.trichotomized.markers (dat.vac.seroneg.allregions, all.markers, wt.col.name="wt")
marker.cutpoints=attr(dat.vac.seroneg.allregions, "marker.cutpoints")
for (a in all.markers) {        
  q.a=marker.cutpoints[[a]]
  if (startsWith(a, "Day")) {
    write(paste0(gsub("_", "\\\\_", a, fixed = TRUE),     " [", concatList(round(q.a, 2), ", "), ")%"), file=paste0(save.results.to, "cutpoints_", a, "_"%.%study_name))
  }
}

# add placebo counterpart
dat.pla.seroneg.allregions=subset(dat.mock, Trt==0 & ph1)

# for validation use
rv=list() 

regions=c("US","LatAm","RSA")



################################################################################
# loop through regions. 1: US, 2: LatAm, 3: RSA

for (iRegion in 1:3) {
# iRegion=1; variant=variants$US[1]
# iRegion=3; variant=variants$RSA[1]
  region=regions[iRegion]
  
  # subset dataset to region
  dat.vac.seroneg=subset(dat.vac.seroneg.allregions, Region==iRegion-1)
  dat.pla.seroneg=subset(dat.pla.seroneg.allregions, Region==iRegion-1)
  
  # loop through variants within this region
  for (variant in variants[[iRegion]]) {
    print("==========================================")
    myprint(region, variant)
    
    ############################
    # count ph1 and ph2 cases
    
    # imputed
    tabs=sapply(1:10, simplify="array", function (imp) {
      dat.vac.seroneg$EventIndOfInterest = ifelse(dat.vac.seroneg$EventIndPrimary==1 & dat.vac.seroneg[["seq1.variant.hotdeck"%.%imp]]==variant, 1, 0)
      with(dat.vac.seroneg, table(ph2, EventIndOfInterest))
    })
    tab =apply(tabs, c(1,2), mean)
    names(dimnames(tab))[2]="Event Indicator"
    mytex(tab, file.name=paste0("tab1_",region,"_",variant), save2input.only=T, input.foldername=save.results.to, digits=1)

    # non-imputed
    dat.vac.seroneg$EventIndOfInterest = ifelse(dat.vac.seroneg$EventIndPrimary==1 & dat.vac.seroneg[["seq1.variant"]]==variant, 1, 0)
    tab = with(dat.vac.seroneg, table(ph2, EventIndOfInterest)) # NA not counted, which is what we want
    names(dimnames(tab))[2]="Event Indicator"
    mytex(tab, file.name=paste0("tab1_nonimputed_",region,"_",variant), save2input.only=T, input.foldername=save.results.to, digits=0)


    ############################
    # formula for coxph

    form.0 = update(Surv(EventTimePrimaryD29, EventIndOfInterest) ~ 1, as.formula(config$covariates_riskscore))

    # source(here::here("code", "cor_coxph_ph_MI.R"))


    #####################################
    # formula for competing risk analysis
    
    # if there are very few competing events, the coxph for competing event may throw warnings
    
    form.0=list(
      update(Surv(EventTimePrimaryD29, EventIndOfInterest) ~ 1, as.formula(config$covariates_riskscore)),
      update(Surv(EventTimePrimaryD29, EventIndCompeting)  ~ 1, as.formula(config$covariates_riskscore))
    )
  
    tfinal.tpeak = tfinal.tpeak.ls[[region]][[variant]]
    write(tfinal.tpeak, file=paste0(save.results.to, "timepoints_cum_risk_", region, "_", variant))

    source(here::here("code", "cor_coxph_risk_no_marker.R"))

    source(here::here("code", "cor_coxph_risk_bootstrap.R"))
    source(here::here("code", "cor_coxph_risk_plotting.R"))

  } # for variant
  
} # for iRegion

warnings() # print out warnings