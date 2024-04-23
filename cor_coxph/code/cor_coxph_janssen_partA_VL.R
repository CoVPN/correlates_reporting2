{
renv::activate(project = here::here(".."))     
Sys.setenv(TRIAL = "janssen_partA_VL")
Sys.setenv(VERBOSE = 1) 
COR="D29VLvariant" 

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

# COVID lineage: "Ancestral.Lineage", "Alpha", "Beta", "Delta", "Epsilon", "Gamma", "Lambda", "Mu", "Zeta", "Iota"
variants=lapply(tfinal.tpeak.ls, function(x) names(x))
myprint(variants)
}



{
marker.cutpoints=attr(dat_proc, "marker.cutpoints"); marker.cutpoints
for (a in "Day29"%.%assays) {        
  q.a=marker.cutpoints[[a]]
  if (startsWith(a, "Day")) {
    write(paste0(gsub("_", " ", a, fixed = TRUE),     " [", concatList(round(q.a, 2), ", "), ")%"), file=paste0(save.results.to, "cutpoints_", a))
  }
}

# add placebo counterpart
dat.vacc.seroneg.allregions=subset(dat_proc, Trt==1 & ph1)
dat.plac.seroneg.allregions=subset(dat_proc, Trt==0 & ph1)

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
}


################################################################################
# loop through regions. 1: US, 2: LatAm, 3: RSA
# saves nevents.df and cox.df

for (iRegion in c(1,2,3)) {
# iRegion=1; variant="Ancestral.Lineage"
# iRegion=2; variant="Gamma" # Lambda Ancestral.Lineage
# iRegion=3; variant="Beta"
  
  region=regions[iRegion]
  dat.vacc=subset(dat.vacc.seroneg.allregions, Region==iRegion-1)
  dat.plac=subset(dat.plac.seroneg.allregions, Region==iRegion-1)
  
  
  # loop through COVID lineage variants within this region
  for (variant in variants[[iRegion]]) {
    cat("==============================  "); myprint(region, variant, newline=F); cat("  ==============================\n")
    
    # append to file names for figures and tables
    fname.suffix = paste0(region, "_", variant)
    for.title = paste0(variant, " COVID, ", region)
    
    tfinal.tpeak = tfinal.tpeak.ls[[region]][[variant]]
    write(tfinal.tpeak, file=paste0(save.results.to, "timepoints_cum_risk_", fname.suffix))

    ############################
    # count ph1 and ph2 cases
    
    {
    # imputed events of interest and ph1/ph2 for variants Ab ph2
    tabs=sapply(1:10, simplify="array", function (imp) {
      dat.vacc$EventIndOfInterest = ifelse(dat.vacc$EventIndPrimary==1 & dat.vacc[["seq1.variant.hotdeck"%.%imp]]==variant, 1, 0)
      with(dat.vacc, table(ph2, EventIndOfInterest))
    })
    tab =apply(tabs, c(1,2), mean)
    names(dimnames(tab))[2]="Event Indicator"
    tab
    mytex(tab, file.name=paste0("tab1_",fname.suffix), save2input.only=T, input.foldername=save.results.to, digits=1)
    # save for forest plots
    nevents.df=rbind(nevents.df, list(region = region, variant = variant, nevents=round(sum(tab[,"1"])) ))
                                      
    
    # imputed events of interest and ph1/ph2 for ancestral Ab ph2
    tabs=sapply(1:10, simplify="array", function (imp) {
      dat.vacc$EventIndOfInterest = ifelse(dat.vacc$EventIndPrimary==1 & dat.vacc[["seq1.variant.hotdeck"%.%imp]]==variant, 1, 0)
      with(dat.vacc, table(ph2.D29, EventIndOfInterest))
    })
    tab =apply(tabs, c(1,2), mean)
    names(dimnames(tab))[2]="Event Indicator"
    tab
    mytex(tab, file.name=paste0("tab1_",fname.suffix,"_ancestral"), save2input.only=T, input.foldername=save.results.to, digits=1)
    # does not have different nevents.df, so not saved 
    
    
    # imputed competing events
    tabs=sapply(1:10, simplify="array", function (imp) {
      dat.vacc$EventIndCompeting = ifelse(dat.vacc$EventIndPrimary==1 & dat.vacc[["seq1.variant.hotdeck"%.%imp]]!=variant, 1, 0)
      with(dat.vacc, table(ph2, EventIndCompeting))
    })
    tab =apply(tabs, c(1,2), mean)
    names(dimnames(tab))[2]="Event Indicator"
    tab
    # show cases counts only b/c controls numbers are a mix of non-cases and cases of interest
    mytex(tab[,2,drop=F], file.name=paste0("tab1_competing_",fname.suffix), save2input.only=T, input.foldername=save.results.to, digits=1)
    
    
    # non-imputed events of interest
    dat.vacc$EventIndOfInterest = ifelse(dat.vacc$EventIndPrimary==1 & dat.vacc[["seq1.variant"]]==variant, 1, 0)
    tab = with(dat.vacc, table(ph2, EventIndOfInterest)) # NA not counted, which is what we want
    names(dimnames(tab))[2]="Event Indicator"
    mytex(tab, file.name=paste0("tab1_nonimputed_",fname.suffix), save2input.only=T, input.foldername=save.results.to, digits=0)
    
    }
    
    ############################
    # formula for coxph

    form.0 = update(Surv(EventTimePrimaryD29, EventIndOfInterest) ~ 1, as.formula(config$covariates))
  
    # the markers to study depend on region and COVID lineage
    if (iRegion==1) {
      # the only COVID variant is ancestral
      all.markers=c("Day29bindSpike", "Day29pseudoneutid50")
      multivariate_assays = "bindSpike+pseudoneutid50"
      
    } else if (iRegion==2) {
      if (variant=="Ancestral.Lineage") {
        all.markers=c("Day29bindSpike", "Day29pseudoneutid50")
        multivariate_assays = "bindSpike+pseudoneutid50"
        
      } else if (variant=="Gamma") {
        all.markers=c("Day29bindSpike", "Day29bindSpike_P.1", "Day29pseudoneutid50", "Day29pseudoneutid50_Gamma")
        multivariate_assays = c("bindSpike+pseudoneutid50",
                                "bindSpike_P.1+pseudoneutid50_Gamma",
                                "bindSpike+bindSpike_P.1",
                                "pseudoneutid50+pseudoneutid50_Gamma")
        
      } else if (variant=="Lambda") {
        all.markers=c("Day29bindSpike", "Day29bindSpike_C.37", "Day29pseudoneutid50", "Day29pseudoneutid50_Lambda")
        multivariate_assays = c("bindSpike+pseudoneutid50",
                                "bindSpike_C.37+pseudoneutid50_Lambda",
                                "bindSpike+bindSpike_C.37",
                                "pseudoneutid50+pseudoneutid50_Lambda")
        
      } else if (variant=="Mu") {
        all.markers=c("Day29bindSpike", "Day29bindSpike_B.1.621", "Day29pseudoneutid50", "Day29pseudoneutid50_Mu")
        multivariate_assays = c("bindSpike+pseudoneutid50",
                                "bindSpike_B.1.621+pseudoneutid50_Mu",
                                "bindSpike+bindSpike_B.1.621",
                                "pseudoneutid50+pseudoneutid50_Mu")
        
      } else if (variant=="Zeta") {
        all.markers=c("Day29bindSpike", "Day29pseudoneutid50", "Day29pseudoneutid50_Zeta")
        multivariate_assays = c("bindSpike+pseudoneutid50",
                                "pseudoneutid50+pseudoneutid50_Zeta")
        
      }
        
    } else if (iRegion==3) {
      all.markers=c("Day29bindSpike", "Day29bindSpike_B.1.351", "Day29pseudoneutid50", "Day29pseudoneutid50_Beta")
      multivariate_assays = c("bindSpike+pseudoneutid50",
                              "bindSpike+bindSpike_B.1.351",
                              "bindSpike_B.1.351+pseudoneutid50_Beta",
                              "pseudoneutid50+pseudoneutid50_Beta")
    }
    
    all.markers.names.short = assay_metadata$assay_label_short[match(sub("Day29","",all.markers),assays)]
    all.markers.names.long  = assay_metadata$assay_label[match(sub("Day29","",all.markers),assays)]
    names(all.markers.names.short) = all.markers
    names(all.markers.names.long) = all.markers
    
    source(here::here("code", "cor_coxph_ph_MI.R"))
    
    
    #####################################
    # formula for competing risk analysis in risk curves

    # if there are very few competing events, the coxph for competing event may throw warnings

    form.0=list(
      update(Surv(EventTimePrimaryD29, EventIndOfInterest) ~ 1, as.formula(config$covariates)),
      update(Surv(EventTimePrimaryD29, EventIndCompeting)  ~ 1, as.formula(config$covariates))
    )
    
    cor_coxph_risk_no_marker (
      form.0,
      dat=dat.vacc,
      fname.suffix, 
      save.results.to,
      config,
      config.cor,
      tfinal.tpeak,
      
      dat.plac,
      verbose=FALSE
    ) 
    

    cor_coxph_risk_bootstrap (
      form.0,
      dat = dat.vacc,
      fname.suffix,
      save.results.to,
      config,
      config.cor,
      tfinal.tpeak,

      markers=all.markers,

      run.Sgts = F # whether to get risk conditional on continuous S>=s
    )

    cor_coxph_risk_plotting (
      form.0,
      dat = dat.vacc,
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

      dat.plac,
      res.plac.cont,
      prev.plac,
      overall.ve,

      show.ve.curves = T,
      plot.geq = F,
      plot.w.plac = T,
      for.title = for.title
    )


  } # variant for loop
  
} # iRegion for loop

save(nevents.df, cox.df, file=paste0(save.results.to, "cox_results.Rdata"))



###################################################################################################
# Forest plots

# {
#   if(verbose) print("forest plots")
#   
#   load(file=paste0(save.results.to, "cox_results.Rdata")) # nevents.df, cox.df
#   
#   trials=c(US="na", LatAm="la", RSA="sa") # used to load Cox model fits for all COVID cases
#   
#   for (iRegion in c(1,2,3)) {
#     # iRegion=1
#     region=regions[iRegion]
#   
#     # subset dataset to region
#     dat.vacc=subset(dat.vacc.seroneg.allregions, Region==iRegion-1)
#     dat.plac=subset(dat.plac.seroneg.allregions, Region==iRegion-1)
#   
#     # ph1 case count
#     nevents.df=rbind(nevents.df, list(region = region, variant = "All", nevents=sum(dat.vacc$EventIndPrimaryIncludeNotMolecConfirmedD29) ))
#   
#     # load coxph results from regional analyses
#     # the coxph results from analysis ready dataset for janssen_partA_VL are different
#     # b/c cases without VL are not considered ph2 in the janssen_partA_VL dataset
#     load(paste0("output/janssen_",trials[region],"_partA/D29IncludeNotMolecConfirmed/coxph_fits.Rdata"))
#     for (a in all.markers) {
#       res=fits.cont.coef.ls[[a]]
#       cox.df=rbind(cox.df, list(
#         region = region,
#         variant = "All",
#         assay = a,
#         est = exp(res[nrow(res),"HR"]),
#         lb =  exp(res[nrow(res),"(lower"]),
#         ub =  exp(res[nrow(res),"upper)"])
#       ))
#     }
#   }
#   
#   
#   row.order=c(
#     "LatAm_D614 Ab_All",
#     "LatAm_D614 Ab_D614G",
#     "LatAm_D614 Ab_Zeta",
#     # "LatAm_Zeta Ab_Zeta",
#     "LatAm_D614 Ab_Mu",
#     # "LatAm_Mu Ab_Mu",
#     "LatAm_D614 Ab_Gamma",
#     # "LatAm_Gamma Ab_Gamma",
#     "LatAm_D614 Ab_Lambda",
#     # "LatAm_Lambda Ab_Lambda",
#     "US_D614 Ab_All",
#     "US_D614 Ab_D614G",
#     "RSA_D614 Ab_All",
#     "RSA_D614 Ab_Beta"
#     # "RSA_Beta Ab_Beta"
#   )
#   
#   # add row names
#   rownames(nevents.df) = paste0(nevents.df$region, "_D614 Ab_", nevents.df$variant)
#   rownames(nevents.df) = sub("Ancestral.Lineage", "D614G", rownames(nevents.df))
#   # order by row names
#   nevents = nevents.df[row.order,"nevents"]
#   
#   
#   # make forest plots
#   assay.titles=c("Binding Antibody to Spike", "Binding Antibody to RBD", "PsV Neutralization" )
#   aa=c("bindSpike",      "bindRBD",        "pseudoneutid50")
#   names(assay.titles)=aa
#   for (a in aa) {
#     #width and height decide margin
#     # to make CI wider, make width bigger and graphwidth larger
#     # onefile has to be F otherwise there will be an empty page inserted
#     mypdf(onefile=F, width=10,height=6, file=paste0(save.results.to, "hr_forest_", a))
#   
#       # subset by assay
#       est.ci = subset(cox.df, assay==paste0("Day",tpeak,a))
#       # add row names (this needs to be done after subsetting)
#       rownames(est.ci) = paste0(est.ci$region, "_D614", if(a=="pseudoneutid50") "G", " Ab_", est.ci$variant)
#       rownames(est.ci) = sub("Ancestral.Lineage", "D614G", rownames(est.ci))
#       # order by rownames
#       est.ci = est.ci[row.order,]
#   
#       theforestplot(
#         nEvents=nevents,
#         point.estimates=est.ci[,"est"],
#         lower.bounds=est.ci[,"lb"],
#         upper.bounds=est.ci[,"ub"],
#         p.values=NA,
#         group=rownames(est.ci),
#         graphwidth=unit(120, "mm"), fontsize=1.2, decimal.places=2,
#         table.labels = c("Group", "HR (95% CI)","No. Events"),
#         title=paste0("Day ", tpeak, " ", assay.titles[a]))
#   
#     dev.off()
#   }
# 
# }


###################################################################################################
if(verbose) print("differences in risk curves")

# load risk bootstrap results

{
# LatAm Ancestral
# vaccine risk
load(paste0(save.results.to, "/risks.all.1_LatAm_Ancestral.Lineage.Rdata"))
LA_AncestralCOVID.AncestralbAb = risks.all.1[["Day29bindSpike"]]
LA_AncestralCOVID.AncestralnAb = risks.all.1[["Day29pseudoneutid50"]]
# placebo risk
load(paste0(save.results.to, "/marginalized.risk.no.marker.LatAm_Ancestral.Lineage.Rdata"))
LA_AncestralCOVID.plac = res.plac.cont
# VE
bAb.est.Ancestral = 1 - LA_AncestralCOVID.AncestralbAb$prob / LA_AncestralCOVID.plac["est"]
bAb.boot.Ancestral = 1 - t( t(LA_AncestralCOVID.AncestralbAb$boot) / LA_AncestralCOVID.plac[2:(1+B)] )
nAb.est.Ancestral = 1 - LA_AncestralCOVID.AncestralnAb$prob / LA_AncestralCOVID.plac["est"]
nAb.boot.Ancestral = 1 - t( t(LA_AncestralCOVID.AncestralnAb$boot) / LA_AncestralCOVID.plac[2:(1+B)] )

# LatAm Gamma
# vaccine risk
load(paste0(save.results.to, "/risks.all.1_LatAm_Gamma.Rdata"))
LA_GammaCOVID.GammabAb = risks.all.1[["Day29bindSpike_P.1"]]
LA_GammaCOVID.GammanAb = risks.all.1[["Day29pseudoneutid50_Gamma"]]
# placebo risk
load(paste0(save.results.to, "/marginalized.risk.no.marker.LatAm_Gamma.Rdata"))
LA_GammaCOVID.plac = res.plac.cont
# VE
bAb.est.Gamma = 1 - LA_GammaCOVID.GammabAb$prob / LA_GammaCOVID.plac["est"]
bAb.boot.Gamma = 1 - t( t(LA_GammaCOVID.GammabAb$boot) / LA_GammaCOVID.plac[2:(1+B)] )
nAb.est.Gamma = 1 - LA_GammaCOVID.GammanAb$prob / LA_GammaCOVID.plac["est"]
nAb.boot.Gamma = 1 - t( t(LA_GammaCOVID.GammanAb$boot) / LA_GammaCOVID.plac[2:(1+B)] )

# LatAm Lambda
# vaccine risk
load(paste0(save.results.to, "/risks.all.1_LatAm_Lambda.Rdata"))
LA_LambdaCOVID.LambdabAb = risks.all.1[["Day29bindSpike_C.37"]]
LA_LambdaCOVID.LambdanAb = risks.all.1[["Day29pseudoneutid50_Lambda"]]
# placebo risk
load(paste0(save.results.to, "/marginalized.risk.no.marker.LatAm_Lambda.Rdata"))
LA_LambdaCOVID.plac = res.plac.cont
# VE
bAb.est.Lambda = 1 - LA_LambdaCOVID.LambdabAb$prob / LA_LambdaCOVID.plac["est"]
bAb.boot.Lambda = 1 - t( t(LA_LambdaCOVID.LambdabAb$boot) / LA_LambdaCOVID.plac[2:(1+B)] )
nAb.est.Lambda = 1 - LA_LambdaCOVID.LambdanAb$prob / LA_LambdaCOVID.plac["est"]
nAb.boot.Lambda = 1 - t( t(LA_LambdaCOVID.LambdanAb$boot) / LA_LambdaCOVID.plac[2:(1+B)] )

# LatAm Mu
# vaccine risk
load(paste0(save.results.to, "/risks.all.1_LatAm_Mu.Rdata"))
LA_MuCOVID.MubAb = risks.all.1[["Day29bindSpike_B.1.621"]]
LA_MuCOVID.MunAb = risks.all.1[["Day29pseudoneutid50_Mu"]]
# placebo risk
load(paste0(save.results.to, "/marginalized.risk.no.marker.LatAm_Mu.Rdata"))
LA_MuCOVID.plac = res.plac.cont
# VE
bAb.est.Mu = 1 - LA_MuCOVID.MubAb$prob / LA_MuCOVID.plac["est"]
bAb.boot.Mu = 1 - t( t(LA_MuCOVID.MubAb$boot) / LA_MuCOVID.plac[2:(1+B)] )
nAb.est.Mu = 1 - LA_MuCOVID.MunAb$prob / LA_MuCOVID.plac["est"]
nAb.boot.Mu = 1 - t( t(LA_MuCOVID.MunAb$boot) / LA_MuCOVID.plac[2:(1+B)] )

# LatAm Zeta, no Zeta-specific bAb
# vaccine risk
load(paste0(save.results.to, "/risks.all.1_LatAm_Zeta.Rdata"))
LA_ZetaCOVID.ZetanAb = risks.all.1[["Day29pseudoneutid50_Zeta"]]
# placebo risk
load(paste0(save.results.to, "/marginalized.risk.no.marker.LatAm_Zeta.Rdata"))
LA_ZetaCOVID.plac = res.plac.cont
# VE
nAb.est.Zeta = 1 - LA_ZetaCOVID.ZetanAb$prob / LA_ZetaCOVID.plac["est"]
nAb.boot.Zeta = 1 - t( t(LA_ZetaCOVID.ZetanAb$boot) / LA_ZetaCOVID.plac[2:(1+B)] )

}


plot.diff = function(marker.a, est.a, boot.a, label.a, 
                     marker.b, est.b, boot.b, label.b,
                     xlab
                     ) {
  
  # find common support, overlap between two [5% percentile, 95% percentile]
  n=100 # a grid of size n
  min = max(marker.a[" 5.0%"], marker.b[" 5.0%"])
  max = min(marker.a["95.0%"], marker.b["95.0%"])
  seq = seq(min, max, len=n)
  
  # interpolate values on the same grid
  # There may be some NAs in y.imputed (out of range of x), it is okay
  
  # curve 1
  select = marker.a>=min & marker.a<=max
  x= c(seq, marker.a[select])
  # impute point est
  est.a.seq = zoo::na.approx(c(rep(NA,n), est.a[select]), x = x, na.rm = FALSE)
  # impute bootstrap values
  boot.a.seq = sapply(1:500, function(i) zoo::na.approx(c(rep(NA,n), boot.a[select,i]), x = x, na.rm = FALSE))
  
  # curve 2
  select = marker.b>=min & marker.b<=max
  x = c(seq, marker.b[select])
  # impute point est
  est.b.seq = zoo::na.approx(c(rep(NA,n), est.b[select]), x = x, na.rm = FALSE)
  # impute bootstrap values
  boot.b.seq = sapply(1:500, function(i) zoo::na.approx(c(rep(NA,n), boot.b[select,i]), x = x, na.rm = FALSE))
  
  # diff
  ve.diff = est.a.seq[1:n] - est.b.seq[1:n]
  ve.diff.boot = boot.a.seq[1:n,] - boot.b.seq[1:n,]
  # pointwise CI
  ci.band=apply(ve.diff.boot, 1, function (x) quantile(x, c(.025,.975), na.rm=T))
  
  # plot
  plot(seq, ve.diff, type="l", xlab=xlab, ylab=paste0("VE Difference"), main=paste0(label.a, " - ", label.b), ylim=c(-1,1), col=4, xaxt="n")
  axis(1, at = seq[seq(1,100,length=6)], labels = round(10**seq[seq(1,100,length=6)]))
  
  lines(seq, ci.band[1,], lty=2, col=4)
  lines(seq, ci.band[2,], lty=2, col=4)
  abline(h=0)
}

# nAb
mypdf(file=paste0(save.results.to, "VE_diff_nAb_anc_lambda"))
plot.diff (
  LA_AncestralCOVID.AncestralnAb$marker, nAb.est.Ancestral, nAb.boot.Ancestral, "Ancestral",
  LA_LambdaCOVID.LambdanAb$marker, nAb.est.Lambda, nAb.boot.Lambda, "Lambda",
  xlab="Lineage-Specific PsV Neutralization"
) 
dev.off()

# bAb
mypdf(file=paste0(save.results.to, "VE_diff_bAb_anc_gamma"))
  plot.diff (
  LA_AncestralCOVID.AncestralbAb$marker, bAb.est.Ancestral, bAb.boot.Ancestral, "Ancestral",
  LA_GammaCOVID.GammabAb$marker, bAb.est.Gamma, bAb.boot.Gamma, "Gamma",
  xlab="Binding Antibody to Lineage-Specific Spike"
) 
dev.off()
  
mypdf(file=paste0(save.results.to, "VE_diff_bAb_anc_lambda"))
plot.diff (
  LA_AncestralCOVID.AncestralbAb$marker, bAb.est.Ancestral, bAb.boot.Ancestral, "Ancestral",
  LA_LambdaCOVID.LambdabAb$marker, bAb.est.Lambda, bAb.boot.Lambda, "Lambda",
  xlab="Binding Antibody to Lineage-Specific Spike"
) 
dev.off()

mypdf(file=paste0(save.results.to, "VE_diff_bAb_anc_mu"))
  plot.diff (
  LA_AncestralCOVID.AncestralbAb$marker, bAb.est.Ancestral, bAb.boot.Ancestral, "Ancestral",
  LA_MuCOVID.MubAb$marker, bAb.est.Mu, bAb.boot.Mu, "Mu",
  xlab="Binding Antibody to Lineage-Specific Spike"
) 
dev.off()
  
mypdf(file=paste0(save.results.to, "VE_diff_bAb_gamma_lambda"))
  plot.diff (
  LA_GammaCOVID.GammabAb$marker, bAb.est.Gamma, bAb.boot.Gamma, "Gamma",
  LA_LambdaCOVID.LambdabAb$marker, bAb.est.Lambda, bAb.boot.Lambda, "Lambda",
  xlab="Binding Antibody to Lineage-Specific Spike"
) 
dev.off()
  
mypdf(file=paste0(save.results.to, "VE_diff_bAb_gamma_mu"))
  plot.diff (
  LA_GammaCOVID.GammabAb$marker, bAb.est.Gamma, bAb.boot.Gamma, "Gamma",
  LA_MuCOVID.MubAb$marker, bAb.est.Mu, bAb.boot.Mu, "Mu",
  xlab="Binding Antibody to Lineage-Specific Spike"
) 
dev.off()
  
mypdf(file=paste0(save.results.to, "VE_diff_bAb_lambda_mu"))
  plot.diff (
  LA_LambdaCOVID.LambdabAb$marker, bAb.est.Lambda, bAb.boot.Lambda, "Lambda",
  LA_MuCOVID.MubAb$marker, bAb.est.Mu, bAb.boot.Mu, "Mu",
  xlab="Binding Antibody to Strain-Specific Spike"
) 
dev.off()
  

################################################################################


print(date())
print("cor_coxph run time: "%.%format(Sys.time()-time.start, digits=1))
