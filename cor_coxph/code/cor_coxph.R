#Sys.setenv(TRIAL = "janssen_pooled_partA"); COR="D29IncludeNotMolecConfirmed"; Sys.setenv(VERBOSE = 1) 
#Sys.setenv(TRIAL = "janssen_pooled_partA"); COR="D29SevereIncludeNotMolecConfirmed"; Sys.setenv(VERBOSE = 1) 
#Sys.setenv(TRIAL = "janssen_pooled_partA"); COR="D29ModerateIncludeNotMolecConfirmed"; Sys.setenv(VERBOSE = 1) 
#Sys.setenv(TRIAL = "prevent19_stage2"); COR="D35prevent19_stage2_severe"; Sys.setenv(VERBOSE = 1)
#Sys.setenv(TRIAL = "moderna_mock"); COR="D29"; Sys.setenv(VERBOSE = 1) 
#Sys.setenv(TRIAL = "azd1222"); COR="D29"; Sys.setenv(VERBOSE = 1) 
#Sys.setenv(TRIAL = "moderna_real"); COR="D57over29"; Sys.setenv(VERBOSE = 1) 
#Sys.setenv(TRIAL = "moderna_real"); COR="D57"; Sys.setenv(VERBOSE = 1) 
#Sys.setenv(TRIAL = "azd1222_bAb"); COR="D57"; Sys.setenv(VERBOSE = 1) 
#Sys.setenv(TRIAL = "janssen_pooled_EUA"); COR="D29IncludeNotMolecConfirmedstart1"; Sys.setenv(VERBOSE = 1) 
#Sys.setenv(TRIAL = "profiscov"); COR="D91over43"; Sys.setenv(VERBOSE = 1) 
#Sys.setenv(TRIAL = "profiscov_lvmn"); COR="D43start48"; Sys.setenv(VERBOSE = 1) 
#Sys.setenv(TRIAL = "azd1222"); COR="D57"; Sys.setenv(VERBOSE = 1) 
#Sys.setenv(TRIAL = "janssen_la_partAsenior"); COR="D29IncludeNotMolecConfirmed"; Sys.setenv(VERBOSE = 1) 
#Sys.setenv(TRIAL = "prevent19"); COR="D35"; Sys.setenv(VERBOSE = 1)
#Sys.setenv(TRIAL = "hvtn705secondNonRSA"); COR="D210"; Sys.setenv(VERBOSE = 1) 
#Sys.setenv(TRIAL = "janssen_na_partA"); COR="D29IncludeNotMolecConfirmed"; Sys.setenv(VERBOSE = 1) 
#Sys.setenv(TRIAL = "janssen_la_partA"); COR="D29IncludeNotMolecConfirmed"; Sys.setenv(VERBOSE = 1) 
#Sys.setenv(TRIAL = "janssen_sa_partA"); COR="D29IncludeNotMolecConfirmed"; Sys.setenv(VERBOSE = 1) 
#Sys.setenv(TRIAL = "hvtn705second"); COR="D210"; Sys.setenv(VERBOSE = 1) 
#Sys.setenv(TRIAL = "janssen_partA_VL"); COR="D29"; Sys.setenv(VERBOSE = 1) 
print(date())
renv::activate(project = here::here(".."))     
source(here::here("..", "_common.R")) # dat_proc is made
if (TRIAL %in% c("janssen_partA_VL")) stop("This TRIAL has its own cor_coxph_TRIAL.R script")    


{
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

# append to file names for figures and tables
# defined differently in cor_coxph_xx.R
fname.suffix = ""


# B=1e3 and numPerm=1e4 take 10 min to run with 30 CPUS for one analysis
B <-       config$num_boot_replicates 
numPerm <- config$num_perm_replicates # number permutation replicates 1e4
myprint(B)
myprint(numPerm)

myprint(tfinal.tpeak)
write(tfinal.tpeak, file=paste0(save.results.to, "timepoints_cum_risk_"%.%fname.suffix))
    

# define an alias for EventIndPrimaryDxx
dat_proc$yy=dat_proc[[config.cor$EventIndPrimary]]

dat.vac.seroneg=subset(dat_proc, Trt==1 & ph1)
dat.pla.seroneg=subset(dat_proc, Trt==0 & ph1)



# define trichotomized markers
if (is.null(attr(dat_proc, "marker.cutpoints"))) {
    dat.vac.seroneg = add.trichotomized.markers (dat.vac.seroneg, all.markers, wt.col.name="wt")
    marker.cutpoints=attr(dat.vac.seroneg, "marker.cutpoints")
} else {
    marker.cutpoints=attr(dat_proc, "marker.cutpoints")
}
for (a in all.markers) {       
    # needed for tertile plots
    dat.vac.seroneg[[a%.%"cat"]] = as.factor(dat.vac.seroneg[[a%.%"cat"]])
    
    q.a=marker.cutpoints[[a]]
    if (startsWith(a, "Day")) {
        # not fold change
        write(paste0(labels.axis[1,marker.name.to.assay(a)], " [", concatList(round(q.a, 2), ", "), ")%"), 
              file=paste0(save.results.to, "cutpoints_", a,".txt"))
    } else {
        # fold change
        write(paste0(escape(a), " [", concatList(round(q.a, 2), ", "), ")%"), 
              file=paste0(save.results.to, "cutpoints_", a,".txt"))
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
mytex(tab, file.name="tab1_"%.%fname.suffix, save2input.only=T, input.foldername=save.results.to)

# for use in competing risk estimation
dat.vac.seroneg.ph2=subset(dat.vac.seroneg, ph2)

begin=Sys.time()
}


###################################################################################################
# estimate overall VE in the placebo and vaccine arms
###################################################################################################

cor_coxph_risk_no_marker (
    form.0,
    dat=dat.vac.seroneg,
    fname.suffix, 
    save.results.to,
    config,
    config.cor,
    tfinal.tpeak,
    
    dat.plac = dat.pla.seroneg,
    verbose=FALSE
) 

if(Sys.getenv("COR_COXPH_NO_MARKER_ONLY")==1) q("no")



###################################################################################################
# run PH models
###################################################################################################
    
cor_coxph_coef_1(
    form.0,
    design_or_dat = design.vacc.seroneg,
    fname.suffix="",
    save.results.to,
    config,
    config.cor,
    all.markers,
    all.markers.names.short,
    
    dat.pla.seroneg,
    show.q = T,
    verbose = T
)



# unit testing of coxph results
if (Sys.getenv("TRIAL") == "janssen_pooled_EUA" & COR=="D29IncludeNotMolecConfirmedstart1") {
    tmp.1=c(rv$tab.1[,4], rv$tab.2[,"overall.p.0"])
    tmp.2=c("0.162","0.079","0.006",      "0.498","   ","   ","0.162","   ","   ","0.003","   ","   ")
    assertthat::assert_that(all(tmp.1==tmp.2), msg = "failed cor_coxph unit testing")    
    
} else if (TRIAL=="moderna_real" & COR=="D57") {
    assertthat::assert_that(all(abs(p.unadj-c(0.004803168, 0.002172787, 0.000129743, 0.000202068, 0.064569846, 0.005631520, 0.009016447, 0.051800145, 0.011506959, 0.579164657))<1e-6), msg = "failed cor_coxph unit testing")    
    
} else if (TRIAL=="prevent19" & COR=="D35") {
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
    

cor_coxph_risk_bootstrap(
    form.0,
    dat = dat.vac.seroneg,
    fname.suffix,
    save.results.to,
    config,
    config.cor,
    tfinal.tpeak,
    
    markers = all.markers,
    
    run.Sgts = T # whether to get risk conditional on continuous S>=s
)


cor_coxph_risk_plotting (
    form.0,
    dat = dat.vac.seroneg,
    fname.suffix,
    save.results.to,
    config,
    config.cor,
    tfinal.tpeak,
    
    markers = all.markers,
    markers.names.short = all.markers.names.short,
    markers.names.long = all.markers.names.long,
    marker.cutpoints,
    assay_metadata,
    
    dat.plac = dat.pla.seroneg,
    res.plac.cont,
    prev.plac,
    overall.ve,
    
    show.ve.curves = T,
    plot.geq = T,
    plot.w.plac = T,
    for.title = ""
)



cor_coxph_risk_tertile_incidence_curves (
    form.0,
    dat = dat.vac.seroneg,
    fname.suffix,
    save.results.to,
    config,
    config.cor,
    tfinal.tpeak,
    
    markers = all.markers,
    markers.names.short = all.markers.names.short,
    markers.names.long = all.markers.names.long,
    marker.cutpoints,
    assay_metadata,
    
    dat.plac = dat.pla.seroneg,
    for.title = ""
)


if (TRIAL %in% c("moderna_real", "janssen_pooled_EUA")) source(here::here("code", "cor_coxph_samplesizeratio.R"))



###################################################################################################
# save verification object rv
save(rv, file=paste0(here::here("verification"), "/", COR, ".rv."%.%study_name%.%".Rdata"))

print(date())
print("cor_coxph run time: "%.%format(Sys.time()-time.start, digits=1))
