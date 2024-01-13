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
# all.markers.names.short=all.markers.names.short[assays]

# add trichotomized markers 
dat.onedosemRNA = subset(dat.mock, ph1.D15 & TrtmRNA==1) # TrtonedosemRNA := TrtmRNA==1 & arm!=3
dat.onedosemRNA$ph2=T
dat.onedosemRNA = add.trichotomized.markers (dat.onedosemRNA, all.markers, ph2.col.name="ph2", wt.col.name="wt.D15")
marker.cutpoints = attr(dat.onedosemRNA, "marker.cutpoints")

# save cut points to files
for (a in all.markers) {        
  write(paste0(gsub("_", "\\_", a, fixed = TRUE),     " [", concatList(round(marker.cutpoints[[a]], 2), ", "), ")%"), 
        file=paste0(save.results.to, "cutpoints_", a, "_"%.%study_name))
}


################################################################################
# Obj 1


fname.suffix = study_name

form.0 = update(Surv(COVIDtimeD22toD181, COVIDIndD22toD181) ~ 1, as.formula(config$covariates_riskscore))

dat=dat.onedosemRNA
dat$yy=dat$COVIDIndD22toD181

tab=with(dat, table(ph2, COVIDIndD22toD181))
names(dimnames(tab))[2]="Event Indicator"
print(tab)
mytex(tab, file.name="tab1", save2input.only=T, input.foldername=save.results.to)


source(here::here("code", "cor_coxph_ph_cohort.R"))
  



print(date())
print("cor_coxph run time: "%.%format(Sys.time()-time.start, digits=1))
