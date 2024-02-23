# COR="D15to181"
# COR="D15to181BA45"
# COR="D15to91"
# COR="D92to181"

renv::activate(project = here::here(".."))
Sys.setenv(TRIAL = "covail"); 
source(here::here("..", "_common.R")) 


{
Sys.setenv(VERBOSE = 1) 
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


marker.cutpoints = attr(dat.mock, "marker.cutpoints")
# save cut points to files
for (a in names(marker.cutpoints)) {        
  write(paste0(gsub("_", "\\_", a, fixed = TRUE),     " [", concatList(round(marker.cutpoints[[a]], 2), ", "), ")%"), 
        file=paste0(save.results.to, "cutpoints_", a))
}

# create scaled version of markers for later use, which is necessary because we do not want to do scaling within naive and non-naive separately
for (a in c("Day15"%.%assays, "B"%.%assays, "Delta15overB"%.%assays)) {
  dat.mock[[a%.%"scaled"]] = scale(dat.mock[[a]], scale=F)
}

assays = c("pseudoneutid50_D614G", "pseudoneutid50_Delta", "pseudoneutid50_Beta", "pseudoneutid50_BA.1", "pseudoneutid50_BA.4.BA.5", "pseudoneutid50_MDW")
dat.sanofi = subset(dat.mock, ph1.D15 & TrtSanofi==1)
dat.sanofi$ph2=1
dat.onedosemRNA = subset(dat.mock, ph1.D15 & TrtonedosemRNA==1) 
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
  form.0 = update(Surv(COVIDtimeD22toD181, COVIDIndD22toD181) ~ 1, as.formula(config$covariates_riskscore))
  dat.sanofi$yy     =dat.sanofi$COVIDIndD22toD181
  dat.onedosemRNA$yy=dat.onedosemRNA$COVIDIndD22toD181
  
} else if (COR=="D15to91") {
  form.0 = update(Surv(COVIDtimeD22toD91, COVIDIndD22toD91) ~ 1, as.formula(config$covariates_riskscore))
  dat.onedosemRNA$yy=dat.onedosemRNA$COVIDIndD22toD91
  dat.sanofi$yy     =dat.sanofi$COVIDIndD22toD91
  
} else if (COR=="D92to181") {
  form.0 = update(Surv(COVIDtimeD92toD181, COVIDIndD92toD181) ~ 1, as.formula(config$covariates_riskscore))
  dat.onedosemRNA$yy=dat.onedosemRNA$COVIDIndD92toD181
  dat.sanofi$yy     =dat.sanofi$COVIDIndD92toD181

} else if (COR=="D15to181BA45") {
  dat.onedosemRNA$EventIndOfInterest = ifelse(dat.onedosemRNA$COVIDIndD22toD181==1 & dat.onedosemRNA$COVIDlineage %in% c("BA.4","BA.5"), 1, 0)
  form.0 = update(Surv(COVIDtimeD22toD181, EventIndOfInterest) ~ 1, as.formula(config$covariates_riskscore))
  dat.sanofi$yy     =dat.sanofi$EventIndOfInterest
  dat.onedosemRNA$yy=dat.onedosemRNA$EventIndOfInterest

} else if (COR=="D15to91BA45") {
  dat.onedosemRNA$EventIndOfInterest = ifelse(dat.onedosemRNA$COVIDIndD22toD91==1 & dat.onedosemRNA$COVIDlineage %in% c("BA.4","BA.5"), 1, 0)
  form.0 = update(Surv(COVIDtimeD22toD91, EventIndOfInterest) ~ 1, as.formula(config$covariates_riskscore))
  dat.onedosemRNA$yy=dat.onedosemRNA$EventIndOfInterest
  dat.sanofi$yy     =dat.sanofi$EventIndOfInterest

} else if (COR=="D92to181BA45") {
  dat.onedosemRNA$EventIndOfInterest = ifelse(dat.onedosemRNA$COVIDIndD92toD181==1 & dat.onedosemRNA$COVIDlineage %in% c("BA.4","BA.5"), 1, 0)
  form.0 = update(Surv(COVIDtimeD92toD181, EventIndOfInterest) ~ 1, as.formula(config$covariates_riskscore))
  dat.onedosemRNA$yy=dat.onedosemRNA$EventIndOfInterest
  dat.sanofi$yy     =dat.sanofi$EventIndOfInterest

} 


form.1 = update(form.0, ~.-naive)

dat.n = subset(dat.onedosemRNA, naive==1)
dat.nn = subset(dat.onedosemRNA, naive==0)

assays = c("pseudoneutid50_D614G", "pseudoneutid50_Delta", "pseudoneutid50_Beta", "pseudoneutid50_BA.1", "pseudoneutid50_BA.4.BA.5", "pseudoneutid50_MDW")

# parameters for ph R scripts
show.q=F # no multiplicity adjustment
use.svy = F 
has.plac = F
show.q = F
}




################################################################################
# Peak Obj 1-3 for mRNA vaccines

# 1 and 11 are both Obj 1, no itxn
# 2 and 21 are both Obj 2, itxn by naive
# 3 and 31 are both Obj 3, itxn by B
for (iObj in c(1,11,12,2,21,3,31)) {
# iObj=1
  
  # define all.markers
  if(iObj==1) {
    all.markers = c("B"%.%assays, "Day15"%.%assays, "Delta15overB"%.%assays)
    all.markers.names.short = assay_metadata$assay_label_short[match(assays,assay_metadata$assay)]
    all.markers.names.short = c("B "%.%all.markers.names.short, "D15 "%.%all.markers.names.short, "D15/B "%.%all.markers.names.short)
    
  } else if(iObj==11){
    # B marker + D15/B
    all.markers = sapply(assays, function (a) paste0("B",a, "scaled + Delta15overB",a, "scaled"))
    all.markers.names.short = sub("Pseudovirus-", "", assay_metadata$assay_label_short[match(assays,assay_metadata$assay)])
    all.markers.names.short = all.markers.names.short
    # parameters for R script
    nCoef=2
    col.headers=c("center(B)", "center(D15/B)")
    
  } else if(iObj==12){
    # B marker + D15 + D15^2
    all.markers = sapply(assays, function (a) paste0("B",a, "scaled + Day15",a, "scaled + I(Day15",a, "scaled^2)"))
    all.markers.names.short = sub("Pseudovirus-", "", assay_metadata$assay_label_short[match(assays,assay_metadata$assay)])
    all.markers.names.short = all.markers.names.short
    # parameters for R script
    nCoef=3
    col.headers=c("center(B)", "center(D15)", "center(D15)**2")
    
  } else if(iObj==2){
    # interaction naive x D15 marker
    all.markers = c(sapply(assays, function (a) paste0("naive * Day15",a, "scaled")),
                    sapply(assays, function (a) paste0("naive * Delta15overB",a, "scaled")))
    all.markers.names.short = sub("Pseudovirus-", "", assay_metadata$assay_label_short[match(assays,assay_metadata$assay)])
    all.markers.names.short = c("D15 "%.%all.markers.names.short, "D15/B "%.%all.markers.names.short)
    # parameters for R script
    nCoef=3
    col.headers=c("naive", "center(D15 or fold)", "naive:center(D15 or fold)")
    
  } else if(iObj==21){
    # interaction (1-naive) x D15 marker
    all.markers = c(sapply(assays, function (a) paste0("I(1-naive) * Day15",a, "scaled")),
                    sapply(assays, function (a) paste0("I(1-naive) * Delta15overB",a, "scaled")))
    all.markers.names.short = sub("Pseudovirus-", "", assay_metadata$assay_label_short[match(assays,assay_metadata$assay)])
    all.markers.names.short = c("D15 "%.%all.markers.names.short, "D15/B "%.%all.markers.names.short)
    nCoef=3
    col.headers=c("I(1-naive)", "center(D15 or fold)", "I(1-naive):center(D15 or fold)")

  } else if(iObj==3){
    # interaction B marker x D15 marker or D15/B
    all.markers = c(sapply(assays, function (a) paste0("B",a, "scaled * Day15",a, "scaled")),
                    sapply(assays, function (a) paste0("B",a, "scaled * Delta15overB",a, "scaled"))
    )
    all.markers.names.short = sub("Pseudovirus-", "", assay_metadata$assay_label_short[match(assays,assay_metadata$assay)])
    all.markers.names.short = c("D15 "%.%all.markers.names.short, "D15/B "%.%all.markers.names.short)
    # parameters for R script
    nCoef=3
    col.headers=c("center(B)", "center(D15 or fold)", "center(B):center(D15 or fold)")
    
  } else if(iObj==31){
    # B_cat marker x D15 marker or D15/B
    all.markers = c(sapply(assays, function (a) paste0("B",a, "cat * Day15",a, "scaled")),
                    sapply(assays, function (a) paste0("B",a, "cat * Delta15overB",a, "scaled")))
    all.markers.names.short = sub("Pseudovirus-", "", assay_metadata$assay_label_short[match(assays,assay_metadata$assay)])
    all.markers.names.short = c("D15 "%.%all.markers.names.short, "D15/B "%.%all.markers.names.short)
    # table col names
    col1="Med B:center(D15 or fold)"
    col2="High B:center(D15 or fold)"
    col3="center(D15 or fold) at Low B"
    col4="center(D15 or fold) at Med B"
    col5="center(D15 or fold) at High B"
    
  }
  
  # same formula, repeat 7 subpopulations
  for (iPop in 1:7) {
    if (iPop==1) {
      dat=dat.onedosemRNA
      fname.suffix = 'mRNA_onedose'
      
    } else if (iPop==2) {
      dat=subset(dat.onedosemRNA, TrtA==1)
      fname.suffix = 'mRNA_Moderna'
    } else if (iPop==3) {
      dat=subset(dat.onedosemRNA, TrtA==0)
      fname.suffix = 'mRNA_Pfizer'
      
    } else if (iPop==4) {
      dat=subset(dat.onedosemRNA, TrtB==1)
      fname.suffix = 'mRNA_Prototype'
    } else if (iPop==5) {
      dat=subset(dat.onedosemRNA, TrtB==0)
      fname.suffix = 'mRNA_Omicron-Containing'
      
    } else if (iPop==6) {
      dat=subset(dat.onedosemRNA, TrtC==1)
      fname.suffix = 'mRNA_Bivalent'
    } else if (iPop==7) {
      dat=subset(dat.onedosemRNA, TrtC==0)
      fname.suffix = 'mRNA_Monovalent'
    } 
    
    if(iObj==1) {
      source(here::here("code", "cor_coxph_ph.R"))
      
    } else if(iObj==11) {
      fname.suffix = paste0(fname.suffix, "_B+D15")
      source(here::here("code", "cor_coxph_ph_coef.R"))
      
      # unit testing 
      if (iPop==1 & COR=="D15to181") {
        # print(pvals.cont)
        assertthat::assert_that(all(
          abs(pvals.cont-c(0.304557, 0.414779, 0.650192, 0.746262, 0.771249, 0.892435))/pvals.cont < 1e-6
        ), msg = "failed cor_coxph unit testing")    
        print("Passed cor_coxph unit testing")    
      }
      
    } else if(iObj==12) {
      fname.suffix = paste0(fname.suffix, "_B+D15^2")
      source(here::here("code", "cor_coxph_ph_coef.R"))
      
    } else if(iObj==2) {
      fname.suffix = paste0(fname.suffix, "_NxD15")
      source(here::here("code", "cor_coxph_ph_coef.R"))
      
    } else if(iObj==21) {
      fname.suffix = paste0(fname.suffix, "_NxD15_2")
      source(here::here("code", "cor_coxph_ph_coef.R"))
      
    } else if(iObj==3) {
      fname.suffix = paste0(fname.suffix, "_BxD15")
      source(here::here("code", "cor_coxph_ph_coef.R"))
      
    } else if(iObj==31) {
      fname.suffix = paste0(fname.suffix, "_BxD15_cat")
      source(here::here("code", "cor_coxph_ph_coef_itxn3.R"))
    }
  }
  
  
  # repeat Obj 1 and 11 in the naive or nonnaive subpopulation
  for (iPop in 1:2) {
    if (iPop==1) {
      dat=subset(dat.onedosemRNA, naive==1)
      fname.suffix = 'mRNA_naive'
    } else if (iPop==2) {
      dat=subset(dat.onedosemRNA, naive==0)
      fname.suffix = 'mRNA_nnaive'
    } 
    
    if(iObj==1) {
      source(here::here("code", "cor_coxph_ph.R"))
    } else if(iObj==11) {
      fname.suffix = paste0(fname.suffix, "_B+D15")
      source(here::here("code", "cor_coxph_ph_coef.R"))
    } else if(iObj==12) {
      fname.suffix = paste0(fname.suffix, "_B+D15^2")
      source(here::here("code", "cor_coxph_ph_coef.R"))
    }
  }
  
  
}


################################################################################
# Peak Obj 4

# same formula, repeated 3 times
for (iPop in 1:3) {
  if (iPop==1) {
    dat=subset(dat.onedosemRNA, TrtA %in% c(1,0))
    fname.suffix = 'mRNA_Mod_Pfi'
    
    all.markers = sapply(assays, function (a) paste0("TrtA * Day15",a, "scaled"))
    all.markers.names.short = sub("Pseudovirus-", "", assay_metadata$assay_label_short[match(assays,assay_metadata$assay)])
    # parameters for R script
    nCoef=3
    col.headers=c("TrtA Moderna~Pfizer", "center(D15)", "TrtA:center(D15)")
    
  } else if (iPop==2) {
    dat=subset(dat.onedosemRNA, TrtB %in% c(1,0))
    fname.suffix = 'mRNA_Pro_Omi'
    
    all.markers = sapply(assays, function (a) paste0("TrtB * Day15",a, "scaled"))
    all.markers.names.short = sub("Pseudovirus-", "", assay_metadata$assay_label_short[match(assays,assay_metadata$assay)])
    # parameters for R script
    nCoef=3
    col.headers=c("TrtB Prot~Omicron", "center(D15)", "TrtB:center(D15)")
    
  } else if (iPop==3) {
    dat=subset(dat.onedosemRNA, TrtC %in% c(1,0))
    fname.suffix = 'mRNA_Bi_Mono'
    
    all.markers = sapply(assays, function (a) paste0("TrtC * Day15",a, "scaled"))
    all.markers.names.short = sub("Pseudovirus-", "", assay_metadata$assay_label_short[match(assays,assay_metadata$assay)])
    # parameters for R script
    nCoef=3
    col.headers=c("TrtC Biv~Monovalent", "center(D15)", "TrtC:center(D15)")
  } 
  
  source(here::here("code", "cor_coxph_ph_coef.R"))

}
  

################################################################################
# Peak Obj 1 for Sanofi vaccines

if (!endsWith(COR,"BA45")) {
  
for (iObj in c(1,11)) {
  
  # define all.markers
  if(iObj==1) {
    all.markers = c("B"%.%assays, "Day15"%.%assays, "Delta15overB"%.%assays)
    all.markers.names.short = assay_metadata$assay_label_short[match(assays,assay_metadata$assay)]
    all.markers.names.short = c("B "%.%all.markers.names.short, "D15 "%.%all.markers.names.short, "D15/B "%.%all.markers.names.short)
    
  } else if(iObj==11){
    # B marker + D15/B
    all.markers = sapply(assays, function (a) paste0("B",a, "scaled + Delta15overB",a, "scaled")
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
  } else if(iObj==11) {
    fname.suffix = paste0(fname.suffix, "_B+D15")
    source(here::here("code", "cor_coxph_ph_coef.R"))
  }
  
}

}


################################################################################
# a quick study of D15 markers in subpopulations defined by baseline titers

res = sapply(assays, function (a) {
      f= update(form.0, as.formula(paste0("~.+", "Day15", a)))
      getFormattedSummary(
        list(coxph(f, dat.onedosemRNA[as.integer(dat.onedosemRNA[[paste0("B",a,"cat")]])==1,]),
             coxph(f, dat.onedosemRNA[as.integer(dat.onedosemRNA[[paste0("B",a,"cat")]])==2,]),
             coxph(f, dat.onedosemRNA[as.integer(dat.onedosemRNA[[paste0("B",a,"cat")]])==3,])),
        type=12, robust=F, exp=T)[4,]
})
tab=t(res)
colnames(tab)=c("L","M","H")
tab
mytex(tab, file.name="CoR_univariable_svycoxph_pretty_Bmarkercat", input.foldername=save.results.to, align="c")



################################################################################
# modeling

llik=NULL; zphglobal=NULL; zphmarker=NULL; coef.ls=list()

for (i in 1:3) { # three different populations: N+NN, N, NN
  if (i==1) {
    f=form.0 
    dat=dat.onedosemRNA
    suffix="N+NN"
  } else if (i==2) {
    f=form.1
    dat=dat.n
    suffix="N"
  } else {
    f=form.1
    dat=dat.nn
    suffix="NN"
  }
  
  model.names=c("B", "D15", "B+D15","B*D15","B+D15^2","B^2+D15","Bcat*D15cat") 
  fits=list(
     coxph(update(f, ~.+Bpseudoneutid50_MDWscaled), dat)
    ,coxph(update(f, ~.+Day15pseudoneutid50_MDWscaled), dat)
    ,coxph(update(f, ~.+Bpseudoneutid50_MDWscaled + Day15pseudoneutid50_MDWscaled), dat)
    ,coxph(update(f, ~.+Bpseudoneutid50_MDWscaled * Day15pseudoneutid50_MDWscaled), dat)
    ,coxph(update(f, ~.+Bpseudoneutid50_MDWscaled + Day15pseudoneutid50_MDWscaled + I(Day15pseudoneutid50_MDWscaled^2)), dat)
    ,coxph(update(f, ~.+Bpseudoneutid50_MDWscaled + I(Bpseudoneutid50_MDWscaled^2) + Day15pseudoneutid50_MDWscaled), dat)
    ,coxph(update(f, ~.+Bpseudoneutid50_MDWcat * Day15pseudoneutid50_MDWcat), dat)
  )
  names(fits)=model.names
  
  llik = rbind(llik, sapply(fits, function(fit) fit$loglik[2]))
  
  zphglobal = rbind(zphglobal, sapply(fits, function(fit) cox.zph(fit)$table["GLOBAL","p"]))
  
  zphmarker = rbind(zphmarker, sapply(fits[1:2], function(fit) cox.zph(fit)$table[length(coef(fit)),"p"]))
  
  tab=getFormattedSummary(fits[1:6], type=6, robust=F, exp=T)[-(1:ifelse(i==1,3,2)),,drop=F]
  mytex(tab, file.name="mdw_models_"%.%suffix, input.foldername=save.results.to, align="c")
  
}

rownames(llik)<-rownames(zphglobal)<-rownames(zphmarker)<-c("N+NN","N","NN")
print(llik)
print(zphglobal)
print(zphmarker)

mytex(llik, file.name="mdw_llik", input.foldername=save.results.to, align="c")
mytex(zphglobal, file.name="mdw_zphglobal", input.foldername=save.results.to, align="c")
mytex(zphmarker, file.name="mdw_zphmarker", input.foldername=save.results.to, align="c")


# myboxplot(Day15pseudoneutid50_MDW~naive, dat.onedosemRNA)
# tmp = aggregate(Day15pseudoneutid50_MDW~naive, dat.onedosemRNA, summary)
# (tmp[1,]-tmp[2,])[4]


###################################
# discrete at both baseline and D15

for (i in 1:3) { # three different populations: N+NN, N, NN
  if (i==1) {
    dat=dat.onedosemRNA
    suffix="N+NN"
    f=form.0
  } else if (i==2) {
    dat=dat.n
    suffix="N"
    f=form.1
  } else {
    dat=dat.nn
    suffix="NN"
    f=form.1
  }
  
  f=update(f, ~ . 
           + I(B_MDW_L * Day15_MDW_M_H) + I(B_MDW_M_H * Day15_MDW_L)
           + I(B_MDW_M * Day15_MDW_M) + I(B_MDW_M * Day15_MDW_H)
           + I(B_MDW_H * Day15_MDW_M) + I(B_MDW_H * Day15_MDW_H)
  )
  
  fit=coxph(f, dat)
  
  tab=getFormattedSummary(list(fit), robust=F, exp=T); tab
  mytex(tab, file.name="mdw_discrete_discrete_itxn_model_"%.%suffix, input.foldername=save.results.to, align="c")

}





print(date())
print("cor_coxph run time: "%.%format(Sys.time()-time.start, digits=1))
