# COR="D43vat08_combined_M12_nAb"
# COR="D22vat08_combined_M12_nAb"
# COR="D43vat08_combined_M6_nAb.st2.sen"

Sys.setenv(TRIAL = "vat08_combined")
Sys.setenv(VERBOSE = 1) 
renv::activate(project = here::here(".."))     
source(here::here("..", "_common.R")) 
source(here::here("code", "params.R"))


{
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
}


{
tp = substr(COR,2,3) 
    
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


dat.vac.seropos.st1 = subset(dat_proc, Trt==1 & Bserostatus==1 & Trialstage==1 & ph1)
dat.vac.seropos.st2 = subset(dat_proc, Trt==1 & Bserostatus==1 & Trialstage==2 & ph1)
dat.vac.seroneg.st2 = subset(dat_proc, Trt==1 & Bserostatus==0 & Trialstage==2 & ph1)

## get cutpoints and turn trichotomized markers into factors
if(!contain(COR, "st2.sen")) {
  marker.cutpoints = list()
  for (a in all.markers) {
    # get cut points
    tmpname = names(table(dat.vac.seropos.st1[[a%.%"cat"]]))[2]
    tmpname = substr(tmpname, 2, nchar(tmpname)-1)
    tmpname = as.numeric(strsplit(tmpname, ",")[[1]])
    tmpname = setdiff(tmpname,Inf) # if there are two categories, remove the second cut point, which is Inf
    marker.cutpoints[[a]] <- tmpname
    
    dat.vac.seropos.st1[[a%.%"cat"]] = as.factor(dat.vac.seropos.st1[[a%.%"cat"]])
    
    write(paste0(escape(a),     " [", concatList(round(marker.cutpoints[[a]], 2), ", "), ")%"), 
          file=paste0(save.results.to.0%.%"/stage1nnaive/", "cutpoints_", a))
  }
}

marker.cutpoints = list()
for (a in all.markers) {
  # get cut points
  tmpname = names(table(dat.vac.seropos.st2[[a%.%"cat"]]))[2]
  tmpname = substr(tmpname, 2, nchar(tmpname)-1)
  tmpname = as.numeric(strsplit(tmpname, ",")[[1]])
  tmpname = setdiff(tmpname,Inf) # if there are two categories, remove the second cut point, which is Inf
  marker.cutpoints[[a]] <- tmpname
  
  dat.vac.seropos.st2[[a%.%"cat"]] = as.factor(dat.vac.seropos.st2[[a%.%"cat"]])
  
  write(paste0(escape(a),     " [", concatList(round(marker.cutpoints[[a]], 2), ", "), ")%"), 
        file=paste0(save.results.to.0%.%"/stage2nnaive/", "cutpoints_", a))
}

marker.cutpoints = list()
for (a in all.markers) {
  # get cut points
  tmpname = names(table(dat.vac.seroneg.st2[[a%.%"cat"]]))[2]
  tmpname = substr(tmpname, 2, nchar(tmpname)-1)
  tmpname = as.numeric(strsplit(tmpname, ",")[[1]])
  tmpname = setdiff(tmpname,Inf) # if there are two categories, remove the second cut point, which is Inf
  marker.cutpoints[[a]] <- tmpname
  
  dat.vac.seroneg.st2[[a%.%"cat"]] = as.factor(dat.vac.seroneg.st2[[a%.%"cat"]])
  
  write(paste0(escape(a),     " [", concatList(round(marker.cutpoints[[a]], 2), ", "), ")%"), 
        file=paste0(save.results.to.0%.%"/stage2naive/", "cutpoints_", a))
}



# add placebo counterpart
dat.pla.seropos.st1=subset(dat_proc, Trt==0 & ph1 & Trialstage==1 & Bserostatus==1)
dat.pla.seropos.st2=subset(dat_proc, Trt==0 & ph1 & Trialstage==2 & Bserostatus==1)
dat.pla.seroneg.st2=subset(dat_proc, Trt==0 & ph1 & Trialstage==2 & Bserostatus==0)

# for validation use
rv=list() 

}

################################################################################
# loop through 3 analyses

for (iAna in 1:3) {
  # iAna=3
  cat("\n\n\n\n")
  myprint(iAna)
  
  # st2.sen 
  if (iAna==1 & contain(COR,"st2.sen")) next
    
  if (iAna==1) {dat.vacc=dat.vac.seropos.st1; dat.plac=dat.pla.seropos.st1; save.results.to=save.results.to.0%.%"stage1nnaive/"}
  if (iAna==2) {dat.vacc=dat.vac.seropos.st2; dat.plac=dat.pla.seropos.st2; save.results.to=save.results.to.0%.%"stage2nnaive/"}
  if (iAna==3) {dat.vacc=dat.vac.seroneg.st2; dat.plac=dat.pla.seroneg.st2; save.results.to=save.results.to.0%.%"stage2naive/"}

  fname.suffix = study_name
  
  ############################
  # count ph1 and ph2 cases
  
  # imputed events of interest
  tabs=sapply(1:10, simplify="array", function (imp) {
    dat.vacc$EventIndOfInterest = dat.vacc[[config.cor$EventIndPrimary%.%imp]]
    with(dat.vacc, table(ph2, EventIndOfInterest))
  })
  tab =apply(tabs, c(1,2), mean)
  names(dimnames(tab))[2]="Event Indicator"
  tab
  mytex(tab,     
        file.name = "tab1_" %.% fname.suffix, 
        save2input.only=T, 
        input.foldername=save.results.to, 
        digits=1)
  
  
  ############################

  form.0 = update(Surv(EventTimeOfInterest, EventIndOfInterest) ~ 1, as.formula(config$covariates))

  multivariate_assays = config$multivariate_assays
  
  all.markers.names.short = assay_metadata$assay_label_short[match(sub("Day"%.%tp,"",all.markers),assays)]
  all.markers.names.long  = assay_metadata$assay_label[match(sub("Day"%.%tp,"",all.markers),assays)]
  names(all.markers.names.short) = all.markers
  names(all.markers.names.long) = all.markers
  
  source(here::here("code", "cor_coxph_ph_MI.R"))
  


} # iAna loop





print(date())
print("cor_coxph run time: "%.%format(Sys.time()-time.start, digits=1))
