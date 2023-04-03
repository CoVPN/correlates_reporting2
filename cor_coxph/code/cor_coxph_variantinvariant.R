library(kyotil) # p.adj.perm, getFormattedSummary
library(marginalizedRisk)
library(tools) # toTitleCase
library(survey)
library(plotrix) # weighted.hist
library(parallel)
library(forestplot)
library(Hmisc) # wtd.quantile, cut2
library(xtable) # this is a dependency of kyotil

begin=Sys.time()
print(date())

COR="D29xxx"; Sys.setenv(VERBOSE = 1) 

for (TRIAL in c("janssen_pooled_partA")) {


renv::activate(project = here::here(".."))     
source(here::here("..", "_common.R")) # dat.mock is made
source(here::here("code", "params.R"))

# path for figures and tables etc
save.results.to = here::here("output");                                if (!dir.exists(save.results.to))  dir.create(save.results.to)
save.results.to = paste0(save.results.to, "/", attr(config,"config")); if (!dir.exists(save.results.to))  dir.create(save.results.to)
save.results.to = paste0(save.results.to, "/", COR,"/");               if (!dir.exists(save.results.to))  dir.create(save.results.to)
print(paste0("save.results.to equals ", save.results.to))

# B=1e3 and numPerm=1e4 take 10 min to run with 30 CPUS for one analysis
B <-       config$num_boot_replicates 
numPerm <- config$num_perm_replicates # number permutation replicates 1e4
myprint(B)
myprint(numPerm)

# uloq censoring, done here b/c should not be done for immunogenicity reports
# note that if delta are used, delta needs to be recomputed
for (a in assays) {
  for (t in "Day"%.%tpeak ) {
    dat.mock[[t %.% a]] <- ifelse(dat.mock[[t %.% a]] > log10(uloqs[a]), log10(uloqs[a]), dat.mock[[t %.% a]])
  }
}    

# define an alias for EventIndPrimaryDxx
dat.mock$yy=dat.mock[[config.cor$EventIndPrimary]]
dat.mock$yy=dat.mock[[config.cor$EventIndPrimary]]

myprint(tfinal.tpeak)
write(tfinal.tpeak, file=paste0(save.results.to, "timepoints_cum_risk_"%.%study_name))
    
dat.vac.seroneg=subset(dat.mock, Trt==1 & ph1)
dat.pla.seroneg=subset(dat.mock, Trt==0 & ph1)


# define trichotomized markers
dat.vac.seroneg = add.trichotomized.markers (dat.vac.seroneg, all.markers, wt.col.name="wt")
marker.cutpoints=attr(dat.vac.seroneg, "marker.cutpoints")
for (a in all.markers) {        
    q.a=marker.cutpoints[[a]]
    if (startsWith(a, "Day")) {
        # not fold change
        write(paste0(labels.axis[1,get.assay.from.name(a)], " [", concatList(round(q.a, 2), ", "), ")%"), file=paste0(save.results.to, "cutpoints_", a, "_"%.%study_name))
    } else {
        # fold change
        write(paste0(a, " [", concatList(round(q.a, 2), ", "), ")%"), file=paste0(save.results.to, "cutpoints_", a, "_"%.%study_name))
    }
}

    
# table of ph1 and ph2 cases
tab=with(dat.vac.seroneg, table(ph2, EventIndPrimary))
names(dimnames(tab))[2]="Event Indicator"
print(tab)
mytex(tab, file.name="tab1", save2input.only=T, input.foldername=save.results.to)

# for use in competing risk estimation
dat.vac.seroneg.ph2=subset(dat.vac.seroneg, ph2)

f= Surv(EventTimePrimaryIncludeNotMolecConfirmedD29, EventIndPrimaryXXX) ~ as.factor(Region) + risk_score + Day29pseudoneutid50 + Day29ADCP
design.vacc.seroneg<-twophase(id=list(~1,~1), weights=list(NULL,~wt.XXX), subset=~ph2, data=dat.vac.seroneg, method="approx")

fit=svycoxph(f, design=design.vacc.seroneg) 

fits=list(fit)
est=getFormattedSummary(fits, exp=T, robust=T, type=1)
ci= getFormattedSummary(fits, exp=T, robust=T, type=13)
est = paste0(est, " ", ci)
p=  getFormattedSummary(fits, exp=T, robust=T, type=10)

tab=cbind(est, p)
colnames(tab)=c("HR", "P value")
tab

mytex(tab, file.name=paste0("CoR_lineage1", align="c", include.colnames = T, save2input.only=T, input.foldername=save.results.to)

    




} # end for TRIAL
