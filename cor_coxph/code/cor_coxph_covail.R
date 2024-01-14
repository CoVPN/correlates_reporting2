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


# add trichotomized markers to TrtonedosemRNA
dat.onedosemRNA = subset(dat.mock, ph1.D15 & TrtonedosemRNA==1) # TrtonedosemRNA := TrtmRNA==1 & arm!=3
dat.onedosemRNA$ph2=T

assays = c("pseudoneutid50_D614G", "pseudoneutid50_Delta", "pseudoneutid50_Beta", "pseudoneutid50_BA.1", "pseudoneutid50_BA.4.BA.5", "pseudoneutid50_MDW")
all.markers = c("B"%.%assays, "Day15"%.%assays)
dat.onedosemRNA = add.trichotomized.markers (dat.onedosemRNA, all.markers, ph2.col.name="ph2", wt.col.name="wt.D15")
marker.cutpoints = attr(dat.onedosemRNA, "marker.cutpoints")

# save cut points to files
for (a in all.markers) {        
  write(paste0(gsub("_", "\\_", a, fixed = TRUE),     " [", concatList(round(marker.cutpoints[[a]], 2), ", "), ")%"), 
        file=paste0(save.results.to, "cutpoints_", a))
}

show.q=F # no multiplicity adjustment

if (COR=="D15to181") {
  form.0 = update(Surv(COVIDtimeD22toD181, COVIDIndD22toD181) ~ 1, as.formula(config$covariates_riskscore))
  dat.onedosemRNA$yy=dat.onedosemRNA$COVIDIndD22toD181
  
} else if (COR=="D15to91") {
  form.0 = update(Surv(COVIDtimeD22toD91, COVIDIndD22toD91) ~ 1, as.formula(config$covariates_riskscore))
  dat.onedosemRNA$yy=dat.onedosemRNA$COVIDIndD22toD91
  
} else if (COR=="D92to181") {
  form.0 = update(Surv(COVIDtimeD92toD181, COVIDIndD92toD181) ~ 1, as.formula(config$covariates_riskscore))
  dat.onedosemRNA$yy=dat.onedosemRNA$COVIDIndD92toD181
} 

assays = c("pseudoneutid50_D614G", "pseudoneutid50_Delta", "pseudoneutid50_Beta", "pseudoneutid50_BA.1", "pseudoneutid50_BA.4.BA.5", "pseudoneutid50_MDW")


################################################################################
# Peak Obj 1-3

for (iObj in c(1,2,3)) {
  
  # define all.markers
  if(iObj==1) {
    all.markers = c("B"%.%assays, "Day15"%.%assays)
    all.markers.names.short = assay_metadata$assay_label_short[match(assays,assay_metadata$assay)]
    all.markers.names.short = c("B "%.%all.markers.names.short, "D15 "%.%all.markers.names.short)
    
  } else if(iObj==2){
    all.markers = sapply(assays, function (a) paste0("naive * scale(Day15",a, ",scale=F)"))
    all.markers.names.short = sub("Pseudovirus-", "", assay_metadata$assay_label_short[match(assays,assay_metadata$assay)])
    # table col names
    col1="naive"
    col2="center(D15)"
    col3="naive:center(D15)"
    
  } else if(iObj==3){
    all.markers = sapply(assays, function (a) paste0("scale(B",a, ") * scale(Day15",a, ",scale=F)"))
    all.markers.names.short = sub("Pseudovirus-", "", assay_metadata$assay_label_short[match(assays,assay_metadata$assay)])
    # table col names
    col1="center(B)"
    col2="center(D15)"
    col3="center(B):center(D15)"
  }
  
  # same formula, repeated 7 times
  for (iTask in 1:7) {
    if (iTask==1) {
      dat=dat.onedosemRNA
      fname.suffix = 'mRNA_onedose'
      
    } else if (iTask==2) {
      dat=subset(dat.onedosemRNA, TrtA==1)
      fname.suffix = 'mRNA_Moderna'
    } else if (iTask==3) {
      dat=subset(dat.onedosemRNA, TrtA==0)
      fname.suffix = 'mRNA_Pfizer'
      
    } else if (iTask==4) {
      dat=subset(dat.onedosemRNA, TrtB==1)
      fname.suffix = 'mRNA_Prototype'
    } else if (iTask==5) {
      dat=subset(dat.onedosemRNA, TrtB==0)
      fname.suffix = 'mRNA_Omicron-Containing'
      
    } else if (iTask==6) {
      dat=subset(dat.onedosemRNA, TrtC==1)
      fname.suffix = 'mRNA_Bivalent'
    } else if (iTask==7) {
      dat=subset(dat.onedosemRNA, TrtC==0)
      fname.suffix = 'mRNA_Monovalent'
      
    } 
    
    if(iObj==1) {
      source(here::here("code", "cor_coxph_ph_cohort.R"))
      
    } else if(iObj==2) {
      fname.suffix = paste0(fname.suffix, "_NxD15")
      source(here::here("code", "cor_coxph_ph_cohort_itxn.R"))
      
    } else if(iObj==3) {
      fname.suffix = paste0(fname.suffix, "_BxD15")
      source(here::here("code", "cor_coxph_ph_cohort_itxn.R"))
      
    }
    
  }
}


################################################################################
# Peak Obj 4

# same formula, repeated 3 times
for (iTask in 1:3) {
  if (iTask==1) {
    dat=subset(dat.onedosemRNA, TrtA %in% c(1,0))
    fname.suffix = 'mRNA_Mod_Pfi'
    
    all.markers = sapply(assays, function (a) paste0("TrtA * scale(Day15",a, ",scale=F)"))
    all.markers.names.short = sub("Pseudovirus-", "", assay_metadata$assay_label_short[match(assays,assay_metadata$assay)])
    # table col names
    col1="TrtA Moderna~Pfizer"
    col2="center(D15)"
    col3="TrtA:center(D15)"
    
  } else if (iTask==2) {
    dat=subset(dat.onedosemRNA, TrtB %in% c(1,0))
    fname.suffix = 'mRNA_Pro_Omi'
    
    all.markers = sapply(assays, function (a) paste0("TrtB * scale(Day15",a, ",scale=F)"))
    all.markers.names.short = sub("Pseudovirus-", "", assay_metadata$assay_label_short[match(assays,assay_metadata$assay)])
    # table col names
    col1="TrtB Prot~Omicron"
    col2="center(D15)"
    col3="TrtB:center(D15)"
    
  } else if (iTask==3) {
    dat=subset(dat.onedosemRNA, TrtC %in% c(1,0))
    fname.suffix = 'mRNA_Bi_Mono'
    
    all.markers = sapply(assays, function (a) paste0("TrtC * scale(Day15",a, ",scale=F)"))
    all.markers.names.short = sub("Pseudovirus-", "", assay_metadata$assay_label_short[match(assays,assay_metadata$assay)])
    # table col names
    col1="TrtC Biv~Monovalent"
    col2="center(D15)"
    col3="TrtC:center(D15)"
  } 
  
  source(here::here("code", "cor_coxph_ph_cohort_itxn.R"))

}
  
print(date())
print("cor_coxph run time: "%.%format(Sys.time()-time.start, digits=1))
