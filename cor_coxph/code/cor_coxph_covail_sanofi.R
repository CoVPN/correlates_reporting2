# COR="D15to181"
# COR="D15to91"
# COR="D92to181"

# modified from cor_coxph_covail for the COVAIL Sanofi manuscript

renv::activate(project = here::here(".."))
Sys.setenv(TRIAL = "covail_sanofi")
source(here::here("..", "_common.R")) 
source(here::here("code", "params.R"))


{
Sys.setenv(VERBOSE = 1) 
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

multi.imp=FALSE

# path for figures and tables etc
save.results.to = here::here("output"); if (!dir.exists(save.results.to))  dir.create(save.results.to)
save.results.to = paste0(save.results.to, "/", attr(config,"config")); if (!dir.exists(save.results.to))  dir.create(save.results.to)
save.results.to = paste0(save.results.to, "/", COR, "/"); if (!dir.exists(save.results.to))  dir.create(save.results.to)


# B=1e3 and numPerm=1e4 take 10 min to run with 30 CPUS for one analysis
B <-       config$num_boot_replicates 
numPerm <- config$num_perm_replicates # number permutation replicates 1e4
myprint(B, numPerm)


marker.cutpoints = attr(dat_proc, "marker.cutpoints")
# save cut points to files
for (a in names(marker.cutpoints)) {        
  write(paste0(escape(a),     " [", concatList(round(marker.cutpoints[[a]], 2), ", "), ")%"), 
        file=paste0(save.results.to, "cutpoints_", a))
}

# create centered version of markers for later use, which is necessary because we do not want to do scaling within naive and non-naive separately
for (a in c("Day15"%.%assays, "B"%.%assays, "Delta15overB"%.%assays)) {
  dat_proc[[a%.%"centered"]] = scale(dat_proc[[a]], scale=F)
}

assays = c("pseudoneutid50_D614G", "pseudoneutid50_Delta", "pseudoneutid50_Beta", "pseudoneutid50_BA.1", "pseudoneutid50_BA.4.BA.5", "pseudoneutid50_MDW")
dat.sanofi = subset(dat_proc, ph1.D15 & TrtSanofi==1)
dat.sanofi$ph2=1
dat.onedosemRNA = subset(dat_proc, ph1.D15 & TrtonedosemRNA==1) 
dat.onedosemRNA$ph2=1

dat.onedosemRNA$Day15_MDW_L = dat.onedosemRNA$Day15pseudoneutid50_MDWcat=="(-Inf,3.62]"
dat.onedosemRNA$Day15_MDW_M = dat.onedosemRNA$Day15pseudoneutid50_MDWcat=="(3.62,4.07]"
dat.onedosemRNA$Day15_MDW_H = dat.onedosemRNA$Day15pseudoneutid50_MDWcat=="(4.07, Inf]"
dat.onedosemRNA$Day15_MDW_M_H = dat.onedosemRNA$Day15_MDW_M | dat.onedosemRNA$Day15_MDW_H
dat.onedosemRNA$B_MDW_L = dat.onedosemRNA$Bpseudoneutid50_MDWcat=="(-Inf,2.56]"
dat.onedosemRNA$B_MDW_M = dat.onedosemRNA$Bpseudoneutid50_MDWcat=="(2.56,3.27]"
dat.onedosemRNA$B_MDW_H = dat.onedosemRNA$Bpseudoneutid50_MDWcat=="(3.27, Inf]"
dat.onedosemRNA$B_MDW_M_H = dat.onedosemRNA$B_MDW_M | dat.onedosemRNA$B_MDW_H


# all cases have covid lineage observed

if (COR=="D15to181") {
  form.0 = update(Surv(COVIDtimeD22toD181, COVIDIndD22toD181) ~ 1, as.formula(config$covariates))
  dat.onedosemRNA$yy=dat.onedosemRNA$COVIDIndD22toD181
  dat.sanofi$yy     =dat.sanofi$COVIDIndD22toD181
  
} else if (COR=="D15to91") {
  form.0 = update(Surv(COVIDtimeD22toD91, COVIDIndD22toD91) ~ 1, as.formula(config$covariates))
  dat.onedosemRNA$yy=dat.onedosemRNA$COVIDIndD22toD91
  dat.sanofi$yy     =dat.sanofi$COVIDIndD22toD91
  
} else if (COR=="D92to181") {
  form.0 = update(Surv(COVIDtimeD92toD181, COVIDIndD92toD181) ~ 1, as.formula(config$covariates))
  dat.onedosemRNA$yy=dat.onedosemRNA$COVIDIndD92toD181
  dat.sanofi$yy     =dat.sanofi$COVIDIndD92toD181

} else if (COR=="D15to181BA45") {
  form.0 = update(Surv(COVIDtimeD22toD181, EventIndOfInterest) ~ 1, as.formula(config$covariates))
  dat.onedosemRNA$EventIndOfInterest = ifelse(dat.onedosemRNA$COVIDIndD22toD181==1 &  dat.onedosemRNA$COVIDlineage %in% c("BA.4","BA.5"), 1, 0)
  dat.onedosemRNA$EventIndCompeting  = ifelse(dat.onedosemRNA$COVIDIndD22toD181==1 & !dat.onedosemRNA$COVIDlineage %in% c("BA.4","BA.5"), 1, 0)
  dat.onedosemRNA$yy=dat.onedosemRNA$EventIndOfInterest
  
} else if (COR=="D15to91BA45") {
  form.0 = update(Surv(COVIDtimeD22toD91, EventIndOfInterest) ~ 1, as.formula(config$covariates))
  dat.onedosemRNA$EventIndOfInterest = ifelse(dat.onedosemRNA$COVIDIndD22toD91==1 &  dat.onedosemRNA$COVIDlineage %in% c("BA.4","BA.5"), 1, 0)
  dat.onedosemRNA$EventIndCompeting  = ifelse(dat.onedosemRNA$COVIDIndD22toD91==1 & !dat.onedosemRNA$COVIDlineage %in% c("BA.4","BA.5"), 1, 0)
  dat.onedosemRNA$yy=dat.onedosemRNA$EventIndOfInterest
  
} else if (COR=="D92to181BA45") {
  form.0 = update(Surv(COVIDtimeD92toD181, EventIndOfInterest) ~ 1, as.formula(config$covariates))
  dat.onedosemRNA$EventIndOfInterest = ifelse(dat.onedosemRNA$COVIDIndD92toD181==1 &  dat.onedosemRNA$COVIDlineage %in% c("BA.4","BA.5"), 1, 0)
  dat.onedosemRNA$EventIndCompeting  = ifelse(dat.onedosemRNA$COVIDIndD92toD181==1 & !dat.onedosemRNA$COVIDlineage %in% c("BA.4","BA.5"), 1, 0)
  dat.onedosemRNA$yy=dat.onedosemRNA$EventIndOfInterest
  
} else stop("Wrong COR")


form.1 = update(form.0, ~.-naive)

dat.n = subset(dat.onedosemRNA, naive==1)
dat.nn = subset(dat.onedosemRNA, naive==0)

assays = c("pseudoneutid50_D614G", "pseudoneutid50_Delta", "pseudoneutid50_Beta", "pseudoneutid50_BA.1", "pseudoneutid50_BA.4.BA.5", "pseudoneutid50_MDW")

# parameters for ph R scripts
show.q=F # no multiplicity adjustment
use.svy = F 
has.plac = F
}



################################################################################
# Peak Obj 1 for Sanofi vaccines

for (iObj in c(1,11)) {
  
  # define all.markers
  if(iObj==1) {
    all.markers = c("B"%.%assays, "Day15"%.%assays, "Delta15overB"%.%assays)
    all.markers.names.short = assay_metadata$assay_label_short[match(assays,assay_metadata$assay)]
    all.markers.names.short = c("B "%.%all.markers.names.short, "D15 "%.%all.markers.names.short, "D15/B "%.%all.markers.names.short)
    
  } else if(iObj==11){
    # B marker + D15/B
    all.markers = sapply(assays, function (a) paste0("B",a, "centered + Delta15overB",a, "centered")
    )
    all.markers.names.short = sub("Pseudovirus-", "", assay_metadata$assay_label_short[match(assays,assay_metadata$assay)])
    all.markers.names.short = all.markers.names.short
    # parameters for R script
    nCoef=2
    col.headers=c("center(B)", "center(D15/B)")
  }
  
  dat=dat.sanofi
  fname.suffix = 'sanofi'
  
  if(iObj==1) {
    source(here::here("code", "cor_coxph_ph.R"))
    
    # forest plot
    fits = lapply ("Day15"%.%assays, function (a) coxph(update(form.0, as.formula(paste0("~.+", a))), dat) )
    forest.covail (fits, names=assays, fname.suffix, save.results.to)
    
  } else if(iObj==11) {
    fname.suffix = paste0(fname.suffix, "_B+D15overB")
    source(here::here("code", "cor_coxph_ph_coef.R"))
  }
  
}





################################################################################
print(date())
print("cor_coxph run time: "%.%format(Sys.time()-time.start, digits=1))
