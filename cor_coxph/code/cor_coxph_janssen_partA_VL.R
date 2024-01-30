{
Sys.setenv(TRIAL = "janssen_partA_VL"); 
COR="D29variant"; 
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

# "Ancestral.Lineage", "Alpha", "Beta", "Delta", "Epsilon", "Gamma", "Lambda", "Mu", "Zeta", "Iota"
variants=lapply(tfinal.tpeak.ls, function(x) names(x))
}


marker.cutpoints=attr(dat.mock, "marker.cutpoints"); marker.cutpoints
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

nevents.df=data.frame(
  region = character(0), 
  variant = character(0), 
  nevents=integer(0), 
  stringsAsFactors=FALSE)

cox.df=data.frame(
  region = character(0), 
  variant = character(0), 
  assay = character(0), 
  est = numeric(0),
  lb = numeric(0),
  ub = numeric(0),
  pval = numeric(0),
  stringsAsFactors=FALSE)



################################################################################
# loop through regions. 1: US, 2: LatAm, 3: RSA

for (iRegion in c(3,1,2)) { # 3 is put first to help debugging b/c there are fewest cases and most likely to throw errors. 
# iRegion=1; variant="Ancestral.Lineage"
# iRegion=3; variant="Beta"
# iRegion=2; variant="Ancestral.Lineage"
  
  region=regions[iRegion]
  
  # subset dataset to region
  dat.vac=subset(dat.vac.seroneg.allregions, Region==iRegion-1)
  dat.pla=subset(dat.pla.seroneg.allregions, Region==iRegion-1)
  
  
  # loop through variants within this region
  for (variant in variants[[iRegion]]) {
    cat("==============================  "); myprint(region, variant, newline=F); cat("  ==============================\n")
    
    # append to file names for figures and tables
    fname.suffix = paste0(region, "_", variant)
    #fname.suffix = study_name
    
    for.title = paste0(variant, " COVID, ", region)
    
    
    ############################
    # count ph1 and ph2 cases
    
    # imputed events of interest
    tabs=sapply(1:10, simplify="array", function (imp) {
      dat.vac$EventIndOfInterest = ifelse(dat.vac$EventIndPrimary==1 & dat.vac[["seq1.variant.hotdeck"%.%imp]]==variant, 1, 0)
      with(dat.vac, table(ph2, EventIndOfInterest))
    })
    tab =apply(tabs, c(1,2), mean)
    names(dimnames(tab))[2]="Event Indicator"
    tab
    mytex(tab, file.name=paste0("tab1_",fname.suffix), save2input.only=T, input.foldername=save.results.to, digits=1)
    
    # save ph1 case count for forest plots
    nevents.df=rbind(nevents.df, list(region = region, variant = variant, nevents=round(sum(tab[,"1"])) ))
                                      
    # imputed competing events
    tabs=sapply(1:10, simplify="array", function (imp) {
      dat.vac$EventIndCompeting = ifelse(dat.vac$EventIndPrimary==1 & dat.vac[["seq1.variant.hotdeck"%.%imp]]!=variant, 1, 0)
      with(dat.vac, table(ph2, EventIndCompeting))
    })
    tab =apply(tabs, c(1,2), mean)
    names(dimnames(tab))[2]="Event Indicator"
    tab
    mytex(tab[,1,drop=F], file.name=paste0("tab1_competing_",fname.suffix), save2input.only=T, input.foldername=save.results.to, digits=1)
    
    # non-imputed events of interest
    dat.vac$EventIndOfInterest = ifelse(dat.vac$EventIndPrimary==1 & dat.vac[["seq1.variant"]]==variant, 1, 0)
    tab = with(dat.vac, table(ph2, EventIndOfInterest)) # NA not counted, which is what we want
    names(dimnames(tab))[2]="Event Indicator"
    mytex(tab, file.name=paste0("tab1_nonimputed_",fname.suffix), save2input.only=T, input.foldername=save.results.to, digits=0)
    
    ############################
    # formula for coxph

    form.0 = update(Surv(EventTimePrimaryD29, EventIndOfInterest) ~ 1, as.formula(config$covariates_riskscore))

    source(here::here("code", "cor_coxph_ph_MI.R"))
    
    
    #####################################
    # formula for competing risk analysis

    # if there are very few competing events, the coxph for competing event may throw warnings

    form.0=list(
      update(Surv(EventTimePrimaryD29, EventIndOfInterest) ~ 1, as.formula(config$covariates_riskscore)),
      update(Surv(EventTimePrimaryD29, EventIndCompeting)  ~ 1, as.formula(config$covariates_riskscore))
    )

    tfinal.tpeak = tfinal.tpeak.ls[[region]][[variant]]
    write(tfinal.tpeak, file=paste0(save.results.to, "timepoints_cum_risk_", fname.suffix))

    # run analyses
    source(here::here("code", "cor_coxph_risk_no_marker.R"))
    source(here::here("code", "cor_coxph_risk_bootstrap.R"))

    # make tables and figures
    source(here::here("code", "cor_coxph_risk_plotting.R"))

  } # for variant
  
} # for iRegion



###################################################################################################
if(verbose) print("forest plots")

# add cox model fits for all COVID cases
trials=c(US="na", LatAm="la", RSA="sa")
for (iRegion in c(1,2,3)) { 
  # iRegion=1
  region=regions[iRegion]
  
  # subset dataset to region
  dat.vac=subset(dat.vac.seroneg.allregions, Region==iRegion-1)
  dat.pla=subset(dat.pla.seroneg.allregions, Region==iRegion-1)
  
  # ph1 case count
  nevents.df=rbind(nevents.df, list(region = region, variant = "All", nevents=sum(dat.vac$EventIndPrimaryIncludeNotMolecConfirmedD29) ))
  
  # load coxph results from regional analyses
  # the coxph results from analysis ready dataset for janssen_partA_VL are different
  # b/c cases without VL are not considered ph2 in the janssen_partA_VL dataset
  load(paste0("/home/yfong/dev/correlates_reporting2-2/cor_coxph/output/janssen_",trials[region],"_partA/D29IncludeNotMolecConfirmed/coxph_fits.Rdata"))
  for (a in all.markers) {
    res=fits.cont.coef.ls[[a]]
    cox.df=rbind(cox.df, list(
      region = region, 
      variant = "All", 
      assay = a, 
      est = exp(res[nrow(res),"HR"]),
      lb =  exp(res[nrow(res),"(lower"]),
      ub =  exp(res[nrow(res),"upper)"])
    ))  
  }
}


row.order=c(
  "LatAm_D614 Ab_All",
  "LatAm_D614 Ab_D614G",
  "LatAm_D614 Ab_Zeta",
  # "LatAm_Zeta Ab_Zeta",
  "LatAm_D614 Ab_Mu",
  # "LatAm_Mu Ab_Mu",
  "LatAm_D614 Ab_Gamma",
  # "LatAm_Gamma Ab_Gamma",
  "LatAm_D614 Ab_Lambda",
  # "LatAm_Lambda Ab_Lambda",
  "US_D614 Ab_All",
  "US_D614 Ab_D614G",
  "RSA_D614 Ab_All",
  "RSA_D614 Ab_Beta"
  # "RSA_Beta Ab_Beta"
)

# add row names
rownames(nevents.df) = paste0(nevents.df$region, "_D614 Ab_", nevents.df$variant)
rownames(nevents.df) = sub("Ancestral.Lineage", "D614G", rownames(nevents.df))
# order by row names
nevents = nevents.df[row.order,"nevents"]


# make forest plots
assay.titles=c("Binding Antibody to Spike", "Binding Antibody to RBD", "PsV Neutralization" )
aa=c("bindSpike",      "bindRBD",        "pseudoneutid50")
names(assay.titles)=aa
for (a in aa) {
  #width and height decide margin
  # to make CI wider, make width bigger and graphwidth larger 
  # onefile has to be F otherwise there will be an empty page inserted
  mypdf(onefile=F, width=10,height=6, file=paste0(save.results.to, "hr_forest_", a)) 
      
    # subset by assay 
    est.ci = subset(cox.df, assay==paste0("Day",tpeak,a))
    # add row names (this needs to be done after subsetting)
    rownames(est.ci) = paste0(est.ci$region, "_D614", if(a=="pseudoneutid50") "G", " Ab_", est.ci$variant)
    rownames(est.ci) = sub("Ancestral.Lineage", "D614G", rownames(est.ci))
    # order by rownames
    est.ci = est.ci[row.order,]
    
    theforestplot(
      nEvents=nevents, 
      point.estimates=est.ci[,"est"], 
      lower.bounds=est.ci[,"lb"], 
      upper.bounds=est.ci[,"ub"], 
      p.values=NA, 
      group=rownames(est.ci), 
      graphwidth=unit(120, "mm"), fontsize=1.2, decimal.places=2, 
      table.labels = c("Group", "HR (95% CI)","No. Events"), 
      title=paste0("Day ", tpeak, " ", assay.titles[a]))
    
  dev.off()
}



###################################################################################################
if(verbose) print("differences in risk curves")

load(paste0("/home/yfong/dev/correlates_reporting2-2/cor_coxph/output/janssen_partA_VL/D29VL/risks.all.1.LatAm.Gamma.Rdata"))
risks.gamma.covid = risks.all.1$Day29bindSpike
risks.gamma.covid = risks.all.1$Day29bindSpike


print(date())
print("cor_coxph run time: "%.%format(Sys.time()-time.start, digits=1))
