#Sys.setenv(TRIAL = "vat08m_naive"); COR="D43"; Sys.setenv(VERBOSE = 1)
#Sys.setenv(TRIAL = "moderna_mock"); COR="D29"; Sys.setenv(VERBOSE = 1) 
#Sys.setenv(TRIAL = "janssen_pooled_partA"); COR="D29IncludeNotMolecConfirmedstart1"; Sys.setenv(VERBOSE = 1) 
#Sys.setenv(TRIAL = "azd1222"); COR="D57"; Sys.setenv(VERBOSE = 1) 
#Sys.setenv(TRIAL = "moderna_real"); COR="D57over29"; Sys.setenv(VERBOSE = 1) 
#Sys.setenv(TRIAL = "moderna_real"); COR="D57"; Sys.setenv(VERBOSE = 1) 
#Sys.setenv(TRIAL = "azd1222_bAb"); COR="D57"; Sys.setenv(VERBOSE = 1) 
#Sys.setenv(TRIAL = "janssen_pooled_EUA"); COR="D29IncludeNotMolecConfirmedstart1"; Sys.setenv(VERBOSE = 1) 
#Sys.setenv(TRIAL = "profiscov"); COR="D91"; Sys.setenv(VERBOSE = 1) 
#Sys.setenv(TRIAL = "profiscov_lvmn"); COR="D43start48"; Sys.setenv(VERBOSE = 1) 
#Sys.setenv(TRIAL = "azd1222"); COR="D57"; Sys.setenv(VERBOSE = 1) 
#Sys.setenv(TRIAL = "janssen_la_partAsenior"); COR="D29IncludeNotMolecConfirmed"; Sys.setenv(VERBOSE = 1) 
#Sys.setenv(TRIAL = "prevent19"); COR="D35"; Sys.setenv(VERBOSE = 1)
#Sys.setenv(TRIAL = "janssen_la_partA"); COR="D29IncludeNotMolecConfirmed"; Sys.setenv(VERBOSE = 1) 
#Sys.setenv(TRIAL = "hvtn705second"); COR="D210"; Sys.setenv(VERBOSE = 1) 
#Sys.setenv(TRIAL = "janssen_pooled_partA"); COR="D29SevereIncludeNotMolecConfirmed"; Sys.setenv(VERBOSE = 1) 

renv::activate(project = here::here(".."))     
    # There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
    if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))    
source(here::here("..", "_common.R"))
# dat.mock is made

library(kyotil) # p.adj.perm, getFormattedSummary
library(marginalizedRisk)
library(tools) # toTitleCase
library(survey)
library(plotrix) # weighted.hist
library(parallel)
library(forestplot)
library(Hmisc) # wtd.quantile, cut2
library(xtable) # this is a dependency of kyotil
source(here::here("code", "params.R"))
time.start=Sys.time()
myprint(study_name)
myprint(verbose)


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

# some exploratory code
if (config$is_ows_trial) source(here::here("code", "cor_coxph_misc.R"))

#create twophase design object
design.vacc.seroneg<-twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dat.vac.seroneg)
with(dat.vac.seroneg, table(Wstratum, ph2))
    
# create verification object to be populated by the following scripts
rv=list() 
rv$marker.cutpoints=marker.cutpoints

# getting some quantiles
#10**wtd.quantile(dat.vac.seroneg$Day57pseudoneutid50, dat.vac.seroneg$wt, c(0.025, 0.05, seq(.2,.9,by=0.01),seq(.9,.99,by=0.005)))

# table of ph1 and ph2 cases
tab=with(dat.vac.seroneg, table(ph2, EventIndPrimary))
names(dimnames(tab))[2]="Event Indicator"
print(tab)
mytex(tab, file.name="tab1", save2input.only=T, input.foldername=save.results.to)

# for use in competing risk estimation
dat.vac.seroneg.ph2=subset(dat.vac.seroneg, ph2)

begin=Sys.time()
print(date())

# some checks

#with(dat.vac.seroneg.ph2, weighted.mean(Day35bindRBD<log10(100), wt))

#with(dat.vac.seroneg.ph2, table(EventIndPrimary, is.na(seq1.spike.weighted.hamming)))
#with(dat.vac.seroneg.ph2, table(EventIndPrimary, is.na(seq1.log10vl)))

#with(dat.vac.seroneg, table(EventIndPrimary, sieve.status, EventTimePrimary>0))
#with(dat.vac.seroneg, table(EventIndPrimaryD29, sieve.status, EventTimePrimary>0))
#with(dat.vac.seroneg, table(EventIndPrimaryIncludeNotMolecConfirmedD29, sieve.status, EventTimePrimary>0))
#with(dat.vac.seroneg, plot(EventTimePrimary, sieve.time, cex=.2)); abline(0,1)
#subset(dat.vac.seroneg, EventIndPrimary==0 & sieve.status==1)

#with(subset(dat.vac.seroneg, EventIndPrimaryIncludeNotMolecConfirmedD29==1), table(is.na(seq1.variant), EventIndPrimaryHasVLD29))


###################################################################################################
# estimate overall VE in the placebo and vaccine arms
###################################################################################################

source(here::here("code", "cor_coxph_marginalized_risk_no_marker.R"))

if(Sys.getenv("COR_COXPH_NO_MARKER_ONLY")==1) q("no")



###################################################################################################
# run PH models
###################################################################################################
    
source(here::here("code", "cor_coxph_ph.R"))


# unit testing of coxph results
if (Sys.getenv("TRIAL") == "janssen_pooled_EUA" & COR=="D29IncludeNotMolecConfirmedstart1") {
    tmp.1=c(rv$tab.1[,4], rv$tab.2[,"overall.p.0"])
    tmp.2=c("0.162","0.079","0.006",      "0.498","   ","   ","0.162","   ","   ","0.003","   ","   ")
    assertthat::assert_that(all(tmp.1==tmp.2), msg = "failed cor_coxph unit testing")    
    
} else if (attr(config, "config")=="moderna_real" & COR=="D57") {
    assertthat::assert_that(all(abs(p.unadj-c(0.004803168, 0.002172787, 0.000129743, 0.000202068, 0.064569846, 0.005631520, 0.009016447, 0.051800145, 0.011506959, 0.579164657))<1e-6), msg = "failed cor_coxph unit testing")    
    
} else if (attr(config, "config")=="prevent19" & COR=="D35") {
    assertthat::assert_that(all(abs(p.unadj-c(0.000453604, 0.0023274, 0.013258206))<1e-6), msg = "failed cor_coxph unit testing")    
    
}
print("Passed cor_coxph unit testing")    


# forest plots
if(length(config$forestplot_script)==1 & !study_name %in% c("PREVENT19","VAT08m") & !contain(attr(config,"config"),"senior")) {
    tmp=here::here("code", config$forestplot_script)
    if (file.exists(tmp)) source(tmp)
    
    # unit testing 
    if (study_name == "MockCOVE") {
        tmp.1=c(sapply(rv$fr.2[-1], function (x) x[c("HR","p.value"),1])) # concatList(tmp.1, ", ")
        if (tpeak=="29") {
            tmp.2=c(2.19803e-01,3.42813e-06,4.00791e-01,1.55780e-03,2.64497e-01,2.90077e-04,2.52391e-01,3.38292e-04,3.11841e-01,1.09284e-03)
        } else if (tpeak=="57") {
            tmp.2=c(1.17284e-01,4.73761e-11,3.91017e-01,7.49144e-04,2.84943e-01,1.36601e-05,2.44480e-01,9.03454e-06,2.70036e-01,9.12665e-06)
        }
        assertthat::assert_that(
            max(abs(tmp.1-tmp.2)/abs(tmp.2))<1e-5,
            msg = "failed MockCOVE unit testing")    
        print("Passed MockCOVE unit testing")    
    }
}



###################################################################################################
# marginalized risk and controlled VE
###################################################################################################
    
# load ylims.cor[[1]] from D29 analyses, which is a list of two: 1 with placebo lines, 2 without placebo lines.
tmp=paste0(here::here(), paste0("/output/", attr(config,"config"), "/", COR, "/ylims.cor.", study_name, ".Rdata"))
if (file.exists(tmp)) load(tmp) # if it does not exist, the code will find alternative ylim

source(here::here("code", "cor_coxph_marginalized_risk_bootstrap.R"))

source(here::here("code", "cor_coxph_marginalized_risk_plotting.R"))

if (attr(config, "config") %in% c("moderna_real", "janssen_pooled_EUA")) source(here::here("code", "cor_coxph_samplesizeratio.R"))



###################################################################################################
# save verification object rv
save(rv, file=paste0(here::here("verification"), "/", COR, ".rv."%.%study_name%.%".Rdata"))

print(date())
print("cor_coxph run time: "%.%format(Sys.time()-time.start, digits=1))
