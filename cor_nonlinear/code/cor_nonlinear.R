#Sys.setenv(TRIAL = "moderna_mock"); Args=c(COR="D57"); Sys.setenv(VERBOSE = 1) # TRIAL: moderna_mock  moderna_real  janssen_pooled_mock  janssen_pooled_real  janssen_na_mock  hvtn705
renv::activate(project = here::here(".."))    
    # There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
    if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))
    
source(here::here("..", "_common.R"))
#-----------------------------------------------
time.start=Sys.time()
myprint(study_name)
    
library(kyotil) # p.adj.perm, getFormattedSummary
library(parallel)
library(chngpt)
library(marginalizedRisk)
library(xtable) # this is a dependency of kyotil
    
source(here::here("code", "params.R"))

## COR defines the analysis to be done, e.g. D29, D57, D29start1
#Args <- commandArgs(trailingOnly=TRUE)
#if (length(Args)==0) Args=c(COR="D29") 
#COR=Args[1]; myprint(COR)
## COR defines the analysis to be done, e.g. D14
#Args <- commandArgs(trailingOnly=TRUE)
#if (length(Args)==0) Args=c(COR="D57") 
#COR=Args[1]; myprint(COR)
## COR has a set of analysis-specific parameters defined in the config file
#config.cor <- config::get(config = COR)
#tpeak=as.integer(paste0(config.cor$tpeak))
#tpeaklag=as.integer(paste0(config.cor$tpeaklag))
#tfinal.tpeak=as.integer(paste0(config.cor$tfinal.tpeak))
#myprint(tpeak, tpeaklag, tfinal.tpeak)
#if (length(tpeak)==0 | length(tpeaklag)==0 | length(tfinal.tpeak)==0) stop("config "%.%COR%.%" misses some fields")


# save tables and figures to analysis-specific folders
save.results.to = here::here("output")
if (!dir.exists(save.results.to))  dir.create(save.results.to)
save.results.to = paste0(here::here("output"), "/", attr(config,"config"));
if (!dir.exists(save.results.to))  dir.create(save.results.to)
save.results.to = paste0(save.results.to, "/", COR,"/");
if (!dir.exists(save.results.to))  dir.create(save.results.to)
print(paste0("save.results.to equals ", save.results.to))
    
# B=1e3 and numPerm=1e4 take 10 min to run with 30 CPUS for one analysis
B <-       config$num_boot_replicates 
numPerm <- config$num_perm_replicates # number permutation replicates 1e4
myprint(B)
myprint(numPerm)
    
    
###################################################################################################
# uloq censoring
# note that if delta are used, delta needs to be recomputed
    
for (a in assays) {
  for (t in "Day"%.%tpeak ) {
    dat.mock[[t %.% a]] <- ifelse(dat.mock[[t %.% a]] > log10(uloqs[a]), log10(uloqs[a]), dat.mock[[t %.% a]])
  }
}
    


# the following data frame define the phase 1 ptids
# do this after uloq censoring
dat.vac.seroneg=subset(dat.mock, Trt==1 & ph1==1)
dat.pla.seroneg=subset(dat.mock, Trt==0 & ph1==1)

    
# define an alias for EventIndPrimaryDxx
dat.vac.seroneg$yy=dat.vac.seroneg[["EventIndPrimary"]]
dat.pla.seroneg$yy=dat.pla.seroneg[["EventIndPrimary"]]

#hist(dat.vac.seroneg$EventTimePrimaryD29)
#hist(dat.vac.seroneg$EventTimePrimaryD29[dat.vac.seroneg$EventIndPrimaryD29==1])
    
form.0.logistic = update (EventIndPrimary ~ 1, as.formula(config$covariates_riskscore))
print(form.0.logistic)




####################################################################################################

dat.vacc.pop.ph2 = subset(dat.vac.seroneg, ph2==1)

# there are two dependencies on cor_coxph

# load prev.plac, prev.vacc, with bootstrap CI
tmp=paste0(here::here(".."), "/cor_coxph/output/",attr(config,"config"),"/", COR,"/", "marginalized.risk.no.marker.Rdata")
if (file.exists(tmp)) load(tmp) else stop("cannot load prev.plac, prev.vacc")
# if this does not exist, the code will throw error

# load ylims.cor, which is a list of two: 1 with placebo lines, 2 without placebo lines.
tmp=paste0(here::here(".."), "/cor_coxph/output/",attr(config,"config"),"/", COR,"/", "ylims.cor.Rdata")
if (file.exists(tmp)) load(tmp)
# if this does not exist, the code will find alternative ylim



####################################################################################################
# GAM
source(here::here("code", "cor_nonlinear_gam.R"))

print("cor_nonlinear run time: "%.%format(Sys.time()-time.start, digits=1))
