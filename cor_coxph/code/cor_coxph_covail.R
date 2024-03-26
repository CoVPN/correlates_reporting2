# COR="D15to181BA45"
# COR="D15to181"
# COR="D15to91"
# COR="D92to181"

renv::activate(project = here::here(".."))
Sys.setenv(TRIAL = "covail")
source(here::here("..", "_common.R")) 
source(here::here("code", "params.R"))
source(here::here("code", "cor_coxph_risk_bootstrap.R"))
source(here::here("code", "cor_coxph_risk_plotting.R"))


{
Sys.setenv(VERBOSE = 1) 
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

multi.imp=FALSE

# path for figures and tables etc
save.results.to = here::here("output"); if (!dir.exists(save.results.to))  dir.create(save.results.to)
save.results.to = paste0(save.results.to, "/", attr(config,"config")); if (!dir.exists(save.results.to))  dir.create(save.results.to)
save.results.to = paste0(save.results.to, "/", COR, "/"); if (!dir.exists(save.results.to))  dir.create(save.results.to)


# B=1e3 and numPerm=1e4 take 10 min to run with 30 CPUS for one analysis
B <-       config$num_boot_replicates 
numPerm <- config$num_perm_replicates # number permutation replicates 1e4
myprint(B, numPerm)


marker.cutpoints = attr(dat_proc, "marker.cutpoints")
# save cut points to files
for (a in names(marker.cutpoints)) {        
  write(paste0(escape(a),     " [", concatList(round(marker.cutpoints[[a]], 2), ", "), ")%"), 
        file=paste0(save.results.to, "cutpoints_", a))
}

# create centered version of markers for later use, which is necessary because we do not want to do scaling within naive and non-naive separately
for (a in c("Day15"%.%assays, "B"%.%assays, "Delta15overB"%.%assays)) {
  dat_proc[[a%.%"centered"]] = scale(dat_proc[[a]], scale=F)
}

assays = c("pseudoneutid50_D614G", "pseudoneutid50_Delta", "pseudoneutid50_Beta", "pseudoneutid50_BA.1", "pseudoneutid50_BA.4.BA.5", "pseudoneutid50_MDW")
dat.sanofi = subset(dat_proc, ph1.D15 & TrtSanofi==1)
dat.sanofi$ph2=1
dat.onedosemRNA = subset(dat_proc, ph1.D15 & TrtonedosemRNA==1) 
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
  form.0 = update(Surv(COVIDtimeD22toD181, COVIDIndD22toD181) ~ 1, as.formula(config$covariates))
  dat.onedosemRNA$yy=dat.onedosemRNA$COVIDIndD22toD181
  dat.sanofi$yy     =dat.sanofi$COVIDIndD22toD181
  
} else if (COR=="D15to91") {
  form.0 = update(Surv(COVIDtimeD22toD91, COVIDIndD22toD91) ~ 1, as.formula(config$covariates))
  dat.onedosemRNA$yy=dat.onedosemRNA$COVIDIndD22toD91
  dat.sanofi$yy     =dat.sanofi$COVIDIndD22toD91
  
} else if (COR=="D92to181") {
  form.0 = update(Surv(COVIDtimeD92toD181, COVIDIndD92toD181) ~ 1, as.formula(config$covariates))
  dat.onedosemRNA$yy=dat.onedosemRNA$COVIDIndD92toD181
  dat.sanofi$yy     =dat.sanofi$COVIDIndD92toD181

} else if (COR=="D15to181BA45") {
  form.0 = update(Surv(COVIDtimeD22toD181, EventIndOfInterest) ~ 1, as.formula(config$covariates))
  dat.onedosemRNA$EventIndOfInterest = ifelse(dat.onedosemRNA$COVIDIndD22toD181==1 &  dat.onedosemRNA$COVIDlineage %in% c("BA.4","BA.5"), 1, 0)
  dat.onedosemRNA$EventIndCompeting  = ifelse(dat.onedosemRNA$COVIDIndD22toD181==1 & !dat.onedosemRNA$COVIDlineage %in% c("BA.4","BA.5"), 1, 0)
  dat.onedosemRNA$yy=dat.onedosemRNA$EventIndOfInterest
  
} else if (COR=="D15to91BA45") {
  form.0 = update(Surv(COVIDtimeD22toD91, EventIndOfInterest) ~ 1, as.formula(config$covariates))
  dat.onedosemRNA$EventIndOfInterest = ifelse(dat.onedosemRNA$COVIDIndD22toD91==1 &  dat.onedosemRNA$COVIDlineage %in% c("BA.4","BA.5"), 1, 0)
  dat.onedosemRNA$EventIndCompeting  = ifelse(dat.onedosemRNA$COVIDIndD22toD91==1 & !dat.onedosemRNA$COVIDlineage %in% c("BA.4","BA.5"), 1, 0)
  dat.onedosemRNA$yy=dat.onedosemRNA$EventIndOfInterest
  
} else if (COR=="D92to181BA45") {
  form.0 = update(Surv(COVIDtimeD92toD181, EventIndOfInterest) ~ 1, as.formula(config$covariates))
  dat.onedosemRNA$EventIndOfInterest = ifelse(dat.onedosemRNA$COVIDIndD92toD181==1 &  dat.onedosemRNA$COVIDlineage %in% c("BA.4","BA.5"), 1, 0)
  dat.onedosemRNA$EventIndCompeting  = ifelse(dat.onedosemRNA$COVIDIndD92toD181==1 & !dat.onedosemRNA$COVIDlineage %in% c("BA.4","BA.5"), 1, 0)
  dat.onedosemRNA$yy=dat.onedosemRNA$EventIndOfInterest
  
} else stop("Wrong COR")


form.1 = update(form.0, ~.-naive)

dat.n = subset(dat.onedosemRNA, naive==1)
dat.nn = subset(dat.onedosemRNA, naive==0)

prev.vacc = get.marginalized.risk.no.marker(form.0, dat.onedosemRNA, tfinal.tpeak)
myprint(prev.vacc)

assays = c("pseudoneutid50_D614G", "pseudoneutid50_Delta", "pseudoneutid50_Beta", "pseudoneutid50_BA.1", "pseudoneutid50_BA.4.BA.5", "pseudoneutid50_MDW")

# parameters for ph R scripts
show.q=F # no multiplicity adjustment
use.svy = F 
has.plac = F
}




################################################################################
# Peak Obj 1-3 for mRNA vaccines

# 1, 11 and 12 are main effects models
# 2 and 21 are itxn by naive
# 3 and 31 are itxn by baseline markers
for (iObj in c(1,11,12,2,21,3,31)) {
# iObj=1; iPop=1  
  
  # define the list of all.markers to work on
  # an item in the list need not be a single marker but is more like a formula

  if(iObj==1) {
    all.markers = c("B"%.%assays, "Day15"%.%assays, "Delta15overB"%.%assays)
    all.markers.names.short = assay_metadata$assay_label_short[match(assays,assay_metadata$assay)]
    all.markers.names.short = c("B "%.%all.markers.names.short, "D15 "%.%all.markers.names.short, "D15/B "%.%all.markers.names.short)
    
  } else if(iObj==11){
    # B marker + D15/B
    all.markers = sapply(assays, function (a) paste0("B",a, "centered + Delta15overB",a, "centered"))
    all.markers.names.short = sub("Pseudovirus-", "", assay_metadata$assay_label_short[match(assays,assay_metadata$assay)])
    # parameters for R script
    nCoef=2
    col.headers=c("center(B)", "center(D15/B)")
    
  } else if(iObj==12){
    # B marker + D15 + D15^2
    all.markers = sapply(assays, function (a) paste0("B",a, "centered + Day15",a, "centered + I(Day15",a, "centered^2)"))
    all.markers.names.short = sub("Pseudovirus-", "", assay_metadata$assay_label_short[match(assays,assay_metadata$assay)])
    # parameters for R script
    nCoef=3
    col.headers=c("center(B)", "center(D15)", "center(D15)**2")
    
  } else if(iObj==2){
    # interaction naive x D15 marker
    all.markers = c(sapply(assays, function (a) paste0("naive * Day15",a, "centered")),
                    sapply(assays, function (a) paste0("naive * Delta15overB",a, "centered")))
    all.markers.names.short = sub("Pseudovirus-", "", assay_metadata$assay_label_short[match(assays,assay_metadata$assay)])
    all.markers.names.short = c("D15 "%.%all.markers.names.short, "D15/B "%.%all.markers.names.short)
    # parameters for R script
    nCoef=3
    col.headers=c("naive", "center(D15 or fold)", "naive:center(D15 or fold)")
    
  } else if(iObj==21){
    # interaction (1-naive) x D15 marker
    all.markers = c(sapply(assays, function (a) paste0("I(1-naive) * Day15",a, "centered")),
                    sapply(assays, function (a) paste0("I(1-naive) * Delta15overB",a, "centered")))
    all.markers.names.short = sub("Pseudovirus-", "", assay_metadata$assay_label_short[match(assays,assay_metadata$assay)])
    all.markers.names.short = c("D15 "%.%all.markers.names.short, "D15/B "%.%all.markers.names.short)
    nCoef=3
    col.headers=c("I(1-naive)", "center(D15 or fold)", "I(1-naive):center(D15 or fold)")

  } else if(iObj==3){
    # interaction B marker x D15 marker or D15/B
    all.markers = c(sapply(assays, function (a) paste0("B",a, "centered * Day15",a, "centered")),
                    sapply(assays, function (a) paste0("B",a, "centered * Delta15overB",a, "centered"))
    )
    all.markers.names.short = sub("Pseudovirus-", "", assay_metadata$assay_label_short[match(assays,assay_metadata$assay)])
    all.markers.names.short = c("D15 "%.%all.markers.names.short, "D15/B "%.%all.markers.names.short)
    # parameters for R script
    nCoef=3
    col.headers=c("center(B)", "center(D15 or fold)", "center(B):center(D15 or fold)")
    
  } else if(iObj==31){
    # B_cat marker x D15 marker or D15/B
    all.markers = c(sapply(assays, function (a) paste0("B",a, "cat * Day15",a, "centered")),
                    sapply(assays, function (a) paste0("B",a, "cat * Delta15overB",a, "centered")))
    all.markers.names.short = sub("Pseudovirus-", "", assay_metadata$assay_label_short[match(assays,assay_metadata$assay)])
    all.markers.names.short = c("D15 "%.%all.markers.names.short, "D15/B "%.%all.markers.names.short)
    # table col names
    col1="Med B:center(D15 or fold)"
    col2="High B:center(D15 or fold)"
    col3="center(D15 or fold) at Low B"
    col4="center(D15 or fold) at Med B"
    col5="center(D15 or fold) at High B"
    
  }
  names(all.markers.names.short) = all.markers
  
  # repeat all objectives over several subpopulations, save results with different fname.suffix
  
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
      
      # forest plot
      fits = lapply ("Day15"%.%assays, function (a) coxph(update(form.0, as.formula(paste0("~.+", a))), dat) )
      forest.covail (fits, names=assays, fname.suffix, save.results.to)
      
      # risk curves using competing risk for BA45 COVID for all mRNA
      if (COR=="D15to181BA45" & iPop==1) {
        
        cor_coxph_risk_bootstrap(  
          form.0 = list(form.0, as.formula(sub("EventIndOfInterest", "EventIndCompeting", paste0(deparse(form.0,width.cutoff=500))))),
          dat, 
          fname.suffix, 
          save.results.to,
          
          config,
          config.cor,
          
          tfinal.tpeak,
          all.markers = "Day15"%.%assays,
          
          comp.risk=T, 
          run.Sgts=F # whether to get risk conditional on continuous S>=s
        )
        
        cor_coxph_risk_plotting(
          form.0 = list(form.0, as.formula(sub("EventIndOfInterest", "EventIndCompeting", paste0(deparse(form.0,width.cutoff=500))))),
          dat,
          fname.suffix,
          save.results.to,
          config,
          config.cor,
          
          assay_metadata,
          
          tfinal.tpeak,
          all.markers = "Day15"%.%assays,
          all.markers.names.short,
          all.markers.names.long,
          marker.cutpoints,

          multi.imp=F,
          comp.risk=T, 
          
          dat.pla.seroneg = NULL,
          res.plac.cont = NULL,
          prev.plac=NULL,
          
          variant=NULL,
          
          show.ve.curves=F,
          plot.geq = F,
          plot.no.plac = F,
          for.title=""
        )
        
      }
      
    } else if(iObj==11) {
      
      fname.suffix = paste0(fname.suffix, "_B+D15overB")
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
  
  
  # repeat Obj 1, 11 and 12 in the naive or nonnaive subpopulation
  # it may be better to use form.1 b/c the caption uses formula, but it is harder to implement and 
  # it is okay to use form.0 because coxph handles it gracefully by setting everything related to it to NA
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
      
      # forest plot
      fits = lapply ("Day15"%.%assays, function (a) coxph(update(form.0, as.formula(paste0("~.+", a))), dat) )
      forest.covail (fits, names=assays, fname.suffix, save.results.to)
      
    } else if(iObj==11) {
      fname.suffix = paste0(fname.suffix, "_B+D15overB")
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
    
    all.markers = sapply(assays, function (a) paste0("TrtA * Day15",a, "centered"))
    all.markers.names.short = sub("Pseudovirus-", "", assay_metadata$assay_label_short[match(assays,assay_metadata$assay)])
    # parameters for R script
    nCoef=3
    col.headers=c("TrtA Moderna~Pfizer", "center(D15)", "TrtA:center(D15)")
    
  } else if (iPop==2) {
    dat=subset(dat.onedosemRNA, TrtB %in% c(1,0))
    fname.suffix = 'mRNA_Pro_Omi'
    
    all.markers = sapply(assays, function (a) paste0("TrtB * Day15",a, "centered"))
    all.markers.names.short = sub("Pseudovirus-", "", assay_metadata$assay_label_short[match(assays,assay_metadata$assay)])
    # parameters for R script
    nCoef=3
    col.headers=c("TrtB Prot~Omicron", "center(D15)", "TrtB:center(D15)")
    
  } else if (iPop==3) {
    dat=subset(dat.onedosemRNA, TrtC %in% c(1,0))
    fname.suffix = 'mRNA_Bi_Mono'
    
    all.markers = sapply(assays, function (a) paste0("TrtC * Day15",a, "centered"))
    all.markers.names.short = sub("Pseudovirus-", "", assay_metadata$assay_label_short[match(assays,assay_metadata$assay)])
    # parameters for R script
    nCoef=3
    col.headers=c("TrtC Biv~Monovalent", "center(D15)", "TrtC:center(D15)")
  } 
  
  source(here::here("code", "cor_coxph_ph_coef.R"))

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
tab=t(res); colnames(tab)=c("L","M","H") #tab
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
     coxph(update(f, ~.+Bpseudoneutid50_MDWcentered), dat)
    ,coxph(update(f, ~.+Day15pseudoneutid50_MDWcentered), dat)
    ,coxph(update(f, ~.+Bpseudoneutid50_MDWcentered + Day15pseudoneutid50_MDWcentered), dat)
    ,coxph(update(f, ~.+Bpseudoneutid50_MDWcentered * Day15pseudoneutid50_MDWcentered), dat)
    ,coxph(update(f, ~.+Bpseudoneutid50_MDWcentered + Day15pseudoneutid50_MDWcentered + I(Day15pseudoneutid50_MDWcentered^2)), dat)
    ,coxph(update(f, ~.+Bpseudoneutid50_MDWcentered + I(Bpseudoneutid50_MDWcentered^2) + Day15pseudoneutid50_MDWcentered), dat)
    ,coxph(update(f, ~.+Bpseudoneutid50_MDWcat * Day15pseudoneutid50_MDWcat), dat)
  )
  names(fits)=model.names
  
  llik = rbind(llik, sapply(fits, function(fit) fit$loglik[2]))
  
  zphglobal = rbind(zphglobal, sapply(fits, function(fit) cox.zph(fit)$table["GLOBAL","p"]))
  
  zphmarker = rbind(zphmarker, sapply(fits[1:2], function(fit) cox.zph(fit)$table[length(coef(fit)),"p"]))
  
  tab=getFormattedSummary(fits[1:6], type=6, robust=F, exp=T)[-(1:ifelse(i==1,3,2)),,drop=F]
  # # make row names shorter to fit the page
  # rownames(tab)=sub("centered","",rownames(tab))
  mytex(tab, file.name="mdw_models_"%.%suffix, input.foldername=save.results.to, align="c")
  
}

rownames(llik)<-rownames(zphglobal)<-rownames(zphmarker)<-c("N+NN","N","NN")

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
  
  tab=getFormattedSummary(list(fit), robust=F, exp=T); #tab
  mytex(tab, file.name="mdw_discrete_discrete_itxn_model_"%.%suffix, input.foldername=save.results.to, align="c")

}

# 3-way interaction 

dat=dat.onedosemRNA
suffix="N*LMH"
f=update(form.0, ~ . 
         + I(B_MDW_L * Day15_MDW_M_H) * naive + I(B_MDW_M_H * Day15_MDW_L) * naive
         + I(B_MDW_M * Day15_MDW_M) * naive + I(B_MDW_M * Day15_MDW_H) * naive
         + I(B_MDW_H * Day15_MDW_M) * naive + I(B_MDW_H * Day15_MDW_H) * naive
)

fit=coxph(f, dat)

tab=getFormattedSummary(list(fit), robust=F, exp=T); #tab
mytex(tab, file.name="mdw_discrete_discrete_itxn_model_"%.%suffix, input.foldername=save.results.to, align="c")


# naive * D15^2 interaction

dat=dat.onedosemRNA
suffix="N*marker"

fs=list(
  update(form.0, ~ . + Bpseudoneutid50_MDWcentered*naive + Day15pseudoneutid50_MDWcentered*naive+ I(Day15pseudoneutid50_MDWcentered^2)*naive ),
  update(form.0, ~ . + Bpseudoneutid50_MDWcentered + Day15pseudoneutid50_MDWcentered * naive + I(Day15pseudoneutid50_MDWcentered^2) * naive),
  update(form.0, ~ . + Bpseudoneutid50_MDWcentered*naive + Day15pseudoneutid50_MDWcentered+ I(Day15pseudoneutid50_MDWcentered^2) ))
fits=lapply(fs, function (f) coxph(f, dat))
tab=getFormattedSummary(fits, robust=F, exp=T, type=6); rownames(tab)=names(coef(fits[[1]])); 
colnames(tab)=c("","N*D15","N*B"); #tab

mytex(tab, file.name="mdw_discrete_discrete_itxn_model_"%.%suffix, input.foldername=save.results.to, align="l")





################################################################################
print(date())
print("cor_coxph run time: "%.%format(Sys.time()-time.start, digits=1))
