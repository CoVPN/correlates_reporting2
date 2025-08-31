# COR="C0iliad_ib202p"
# renv::activate(project = here::here("..")) #  in .Rprofile
Sys.setenv(TRIAL = "iliad_ib202p")
Sys.setenv(VERBOSE = 1)
source(here::here("..", "_common.R")) 


{
library(kyotil) # p.adj.perm, getFormattedSummary
library(Hmisc) # wtd.quantile, cut2
library(xtable) # this is a dependency of kyotil
library(sandwich)
library(lmtest)
  
  
source(here::here("code", "params.R"))
time.start=Sys.time()
myprint(study_name)
myprint(verbose)

# redefine form.0
form.0 = update (EventIndPrimary~1, as.formula(config$covariates))
print(form.0)


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

# define an alias for EventIndPrimaryDxx
dat_proc$yy=dat_proc[[config.cor$EventIndPrimary]]

# save cut points to files
marker.cutpoints = attr(dat_proc, "marker.cutpoints")
for (a in names(marker.cutpoints)) {        
  write(paste0(escape(a),     " [", concatList(round(marker.cutpoints[[a]], 2), ", "), ")%"), 
        file=paste0(save.results.to, "cutpoints_", a,".txt"))
}

all.markers = c("Day"%.%timepoints%.%assays, "B"%.%assays, "Delta28overB"%.%assays)
all.markers.names.short=c(
  "D28 "%.%labels.assays.short,
  "B "%.%labels.assays.short,
  "Fold-rise "%.% sub("\\(.+\\)", "", labels.assays.short) 
)

show.q=F

begin=Sys.time()
}


trts=1:0
for (trt in trts) {
  
  if (trt==1) {
    dat.ph1=subset(dat_proc, ph1==1 & Trt=="BPZE1")
    fname.suffix <- trt.label <- "BPZE1"
  } else if (trt==0) {
    dat.ph1=subset(dat_proc, ph1==1 & Trt=="PBO")
    fname.suffix <- trt.label <- "PBO"
  }
  
  
  # table of ph1 and ph2 cases
  tab=with(dat.ph1, table(ph2, EventIndPrimary))
  names(dimnames(tab))[2]="Event Indicator"
  print(tab)
  mytex(tab, file.name="tab1_" %.% fname.suffix, save2input.only=T, input.foldername=save.results.to)
  
    
    
  
  ################################################################################
  # get RR on continuous markers
  ################################################################################
  
  
  fits=list()
  fits.scaled=list()
  for (i in 1:2) { # 1: not scaled, 2: scaled
    for (a in all.markers) {
      f=update (form.0,  as.formula('~. + '%.%ifelse(i==2,"scale("%.%a%.%")", a)))
      
      fit <- glm(f, family = poisson(link = "log"), data = dat.ph1)
      
      # coeftest(fit, vcov. = sandwich)    # Robust (Huber-White) standard errors
        
      if (i==1) fits[[a]]=fit else fits.scaled[[a]]=fit
    }
  }
  
  natrisk=nrow(dat.ph1)
  nevents=sum(dat.ph1$yy==1)
  
  # make pretty table
  {
  rows=length(coef(fits[[1]]))
  est=getFormattedSummary(fits, exp=T, robust=F, rows=rows, type=1)
  ci= getFormattedSummary(fits, exp=T, robust=F, rows=rows, type=7)
  p=  getFormattedSummary(fits, exp=T, robust=F, rows=rows, type=10)
  est.scaled=getFormattedSummary(fits.scaled, exp=T, robust=F, rows=rows, type=1)
  ci.scaled= getFormattedSummary(fits.scaled, exp=T, robust=F, rows=rows, type=7)
  }
  
  pvals.cont = sapply(fits, function(x) {
      tmp=getFixedEf(x)
      p.val.col=which(startsWith(tolower(colnames(tmp)),"p"))
      tmp[nrow(tmp),p.val.col]
  })
  
  
  
  
  
  ###################################################################################################
  if(verbose) print("regression for trichotomized markers")
  
  marker.levels = sapply(all.markers, function(a) length(table(dat.ph1[[a%.%"cat"]])))
  
  fits.tri=list()
  for (a in all.markers) {
    if(verbose) myprint(a)
    f= update(form.0, as.formula(paste0("~.+", a, "cat")))
    
    fit <- glm(f, family = poisson(link = "log"), data = dat.ph1)
    
    fits.tri[[a]]=fit
  }
  
  fits.tri.coef.ls= lapply(fits.tri, function (fit) getFixedEf(fit, robust=F))
  
  
  # get generalized Wald p values
  overall.p.tri=sapply(all.markers, function(a) {
      fit=fits.tri[[a]]
      rows=length(fit$coef) - (marker.levels[a]-2):0
      
      if (length(fit)==1) NA else {
          stat=coef(fit)[rows] %*% solve(vcov(fit,robust=F)[rows,rows]) %*% coef(fit)[rows]
          pchisq(stat, length(rows), lower.tail = FALSE)
      }
  })
  #
  overall.p.0=formatDouble(c(rbind(overall.p.tri, NA,NA)), digits=3, remove.leading0 = F);   overall.p.0=sub("0.000","<0.001",overall.p.0)
  
  
  
  ###################################################################################################
  if(verbose) print("# multitesting adjustment for continuous and trichotomized markers together")
  
  # If primary_assays is not defined in config, multitesting adjustment is over all assays. 
  # If primary_assays defined, multitesting adjustment is only over this subset. If this set is empty, then no multitesting adjustment is done
  
  p.unadj=c(cont=pvals.cont, tri=overall.p.tri)
  # save a copy for later use
  p.unadj.1 = p.unadj 
  # pick out a subset based on config
  if (!is.null(config$primary_assays)) {
      if (length(config$primary_assays)>0) {
          p.unadj = p.unadj[c(paste0("cont.",DayPrefix,tpeak,primary_assays), paste0("tri.",DayPrefix,tpeak,primary_assays))]
      } else {
          p.unadj=c()
      }
  }
  
  #### Holm and FDR adjustment
  pvals.adj.fdr=p.adjust(p.unadj, method="fdr")
  pvals.adj.hol=p.adjust(p.unadj, method="holm")
  
  print("not doing Westfall and Young")
  pvals.adj=cbind(p.unadj, p.FWER=pvals.adj.hol, p.FDR=pvals.adj.fdr)
  write(NA, file=paste0(save.results.to, "permutation_replicates_"%.%study_name))     # so the rmd file can compile
  
  
  # since we take ID80 out earlier, we may need to add it back for the table and we do it with the help of p.unadj.1
  pvals.adj = cbind(p.unadj=p.unadj.1, pvals.adj[match(names(p.unadj.1), rownames(pvals.adj)),2:3, drop=F])
  
  
  
  ###################################################################################################
  # make continuous markers table
  
  p.1=formatDouble(pvals.adj["cont."%.%names(pvals.cont),"p.FWER"], 3, remove.leading0=F); p.1=sub("0.000","<0.001",p.1)
  p.2=formatDouble(pvals.adj["cont."%.%names(pvals.cont),"p.FDR" ], 3, remove.leading0=F); p.2=sub("0.000","<0.001",p.2)
  
  
  tab.1=cbind(paste0(nevents, "/", format(natrisk, big.mark=",")), t(est), t(ci), t(p), if(show.q) p.2, if(show.q) p.1)
  rownames(tab.1)=all.markers.names.short
  tab.1
  
  if (show.q) {
    header=paste0("\\hline\n 
           \\multicolumn{1}{l}{", escape(fname.suffix), "}    & \\multicolumn{1}{c}{No. cases /}   & \\multicolumn{2}{c}{RR per 10-fold incr.}                    
              & \\multicolumn{1}{c}{P-value}   & \\multicolumn{1}{c}{q-value}   & \\multicolumn{1}{c}{FWER} \\\\ 
           \\multicolumn{1}{l}{Immunologic Marker}  & \\multicolumn{1}{c}{No. at-risk**} & \\multicolumn{1}{c}{Pt. Est.} & \\multicolumn{1}{c}{95\\% CI} 
              & \\multicolumn{1}{c}{(2-sided)} & \\multicolumn{1}{c}{***}       & \\multicolumn{1}{c}{} \\\\ 
           \\hline\n 
        ")
  } else {
    header=paste0("\\hline\n 
           \\multicolumn{1}{l}{", escape(fname.suffix), "}    & \\multicolumn{1}{c}{No. cases /}   & \\multicolumn{2}{c}{RR per 10-fold incr.}                    
              & \\multicolumn{1}{c}{P-value}   \\\\ 
           \\multicolumn{1}{l}{Immunologic Marker}  & \\multicolumn{1}{c}{No. at-risk**} & \\multicolumn{1}{c}{Pt. Est.} & \\multicolumn{1}{c}{95\\% CI} 
              & \\multicolumn{1}{c}{(2-sided)} \\\\ 
           \\hline\n 
        ")
  }
  mytex(tab.1, file.name="CoR_univariable_poisson_pretty_"%.%fname.suffix, align="c", include.colnames = F, save2input.only=T, input.foldername=save.results.to,
        col.headers=header,
        longtable=T, 
        label=paste0("tab:CoR_univariable_poisson_pretty"), 
        caption.placement = "top", 
        caption=paste0("Inference for ", DayPrefix, tpeak, " antibody marker covariate-adjusted correlates of risk of ", config.cor$txt.endpoint, " pooled over treatment arms: Odds ratios per 10-fold increment in the marker*"))
  tab.cont=tab.1
  
  tab.1.nop12=cbind(paste0(nevents, "/", format(natrisk, big.mark=",")), t(est), t(ci), t(p))
  rownames(tab.1.nop12)=all.markers.names.short
  
  # scaled markers
  tab.1.scaled=cbind(paste0(nevents, "/", format(natrisk, big.mark=",")), t(est.scaled), t(ci.scaled), t(p), if(show.q) p.2, if(show.q) p.1)
  rownames(tab.1.scaled)=all.markers.names.short
  tab.1.scaled
  
  if (show.q) {
    header=paste0("\\hline\n 
           \\multicolumn{1}{l}{", escape(fname.suffix), "} & \\multicolumn{1}{c}{No. cases /}   & \\multicolumn{2}{c}{RR per SD incr.}                     & \\multicolumn{1}{c}{P-value}   & \\multicolumn{1}{c}{q-value}   & \\multicolumn{1}{c}{FWER} \\\\ 
           \\multicolumn{1}{l}{Immunologic Marker}            & \\multicolumn{1}{c}{No. at-risk**} & \\multicolumn{1}{c}{Pt. Est.} & \\multicolumn{1}{c}{95\\% CI} & \\multicolumn{1}{c}{(2-sided)} & \\multicolumn{1}{c}{***} & \\multicolumn{1}{c}{} \\\\ 
           \\hline\n 
      ")
  } else {
    header=paste0("\\hline\n 
           \\multicolumn{1}{l}{", escape(fname.suffix), "} & \\multicolumn{1}{c}{No. cases /}   & \\multicolumn{2}{c}{RR per SD incr.}                     & \\multicolumn{1}{c}{P-value}   \\\\ 
           \\multicolumn{1}{l}{Immunologic Marker}            & \\multicolumn{1}{c}{No. at-risk**} & \\multicolumn{1}{c}{Pt. Est.} & \\multicolumn{1}{c}{95\\% CI} & \\multicolumn{1}{c}{(2-sided)}  \\\\ 
           \\hline\n 
      ")
  }
  
  mytex(tab.1.scaled, file.name="CoR_univariable_poisson_pretty_scaled_"%.%fname.suffix, align="c", include.colnames = F, save2input.only=T, input.foldername=save.results.to,
        col.headers=header,
        longtable=T, 
        label=paste0("tab:CoR_univariable_poisson_pretty_scaled"), 
        caption.placement = "top", 
        caption=paste0("Inference for ", DayPrefix, tpeak, " antibody marker covariate-adjusted correlates of risk of ", config.cor$txt.endpoint, " pooled over treatment groups: Odds ratios per SD increment in the marker*"))
  tab.cont.scaled=tab.1.scaled
  
  
  
  ###################################################################################################
  # make trichotomized markers table
  
  
  get.est=function(a) {
    fit=fits.tri[[a]]
    rows=length(fit$coef) - (marker.levels[a]-2):0
    out = getFormattedSummary(list(fit), exp=T, robust=F, rows=rows, type=1)
    if (length(out)==1) c(NA,out) else out
  }
  get.ci =function(a) {
    fit=fits.tri[[a]]
    rows=length(fit$coef) - (marker.levels[a]-2):0
    out = getFormattedSummary(list(fit), exp=T, robust=F, rows=rows, type=7)
    if (length(out)==1) c(NA,out) else out
  }
  get.p  =function(a) {
    fit=fits.tri[[a]]
    rows=length(fit$coef) - (marker.levels[a]-2):0
    out = getFormattedSummary(list(fit), exp=T, robust=F, rows=rows, type=10)
    if (length(out)==1) c(NA,out) else out
  }
  get.est(all.markers[1]); get.ci (all.markers[1]); get.p  (all.markers[1])
  get.est(all.markers[6]); get.ci (all.markers[6]); get.p  (all.markers[6])
  
  # regression parameters
  est=c(rbind(1.00,  sapply(all.markers, function (a) get.est(a))))
  ci= c(rbind(NA, sapply(all.markers, function (a) get.ci (a))))
  p=  c(rbind(NA, sapply(all.markers, function (a) get.p  (a))))
  
  overall.p.1=formatDouble(pvals.adj["tri."%.%names(pvals.cont),"p.FWER"], 3, remove.leading0=F);   overall.p.1=sub("0.000","<0.001",overall.p.1)
  overall.p.2=formatDouble(pvals.adj["tri."%.%names(pvals.cont),"p.FDR" ], 3, remove.leading0=F);   overall.p.2=sub("0.000","<0.001",overall.p.2)
  
  # add space
  overall.p.1=c(rbind(overall.p.1, NA,NA))
  overall.p.2=c(rbind(overall.p.2, NA,NA))
  
  
  # n cases and n at risk
  natrisk = round(c(sapply (all.markers%.%"cat", function(a) {
    out = aggregate(subset(dat.ph1,ph2==1)        [["wt"]], subset(dat.ph1,ph2==1        )[a], sum, na.rm=T, drop=F)[,2]
    if (length(out)==3) out else c(out[1],NA,out[2])
  } )))
  nevents = round(c(sapply (all.markers%.%"cat", function(a) {
    out = aggregate(subset(dat.ph1,yy==1 & ph2==1)[["wt"]], subset(dat.ph1,yy==1 & ph2==1)[a], sum, na.rm=T, drop=F)[,2]
    out[is.na(out)]=0
    if (length(out)==3) out else c(out[1],NA,out[2])
  } )))
  
  tab=cbind(
      rep(c("Lower","Middle","Upper"), length(p)/3),
      paste0(nevents, "/", format(natrisk, big.mark=",",digit=0, scientific=F)), 
      formatDouble(nevents/natrisk, digit=4, remove.leading0=F),
      est, ci, p, overall.p.0, if(show.q) overall.p.2, if(show.q) overall.p.1
  )
  
  tmp=rbind(all.markers.names.short, "", "")
  rownames(tab)=c(tmp)
  tab
  
  if(show.q) {
    header=paste0("\\hline\n 
           \\multicolumn{1}{l}{", escape(fname.suffix), "} & \\multicolumn{1}{c}{Tertile}   & \\multicolumn{1}{c}{No. cases /}   & \\multicolumn{1}{c}{Attack}   & \\multicolumn{2}{c}{RR}                     & \\multicolumn{1}{c}{P-value}   & \\multicolumn{1}{c}{Overall P-}      & \\multicolumn{1}{c}{Overall q-}   & \\multicolumn{1}{c}{Overall} \\\\ 
           \\multicolumn{1}{l}{Immunologic Marker}            & \\multicolumn{1}{c}{}          & \\multicolumn{1}{c}{No. at-risk**} & \\multicolumn{1}{c}{rate}   & \\multicolumn{1}{c}{Pt. Est.} & \\multicolumn{1}{c}{95\\% CI} & \\multicolumn{1}{c}{(2-sided)} & \\multicolumn{1}{c}{value***} & \\multicolumn{1}{c}{value $\\dagger$} & \\multicolumn{1}{c}{FWER} \\\\ 
           \\hline\n 
        ")
  } else {
    header=paste0("\\hline\n 
           \\multicolumn{1}{l}{", escape(fname.suffix), "} & \\multicolumn{1}{c}{Tertile}   & \\multicolumn{1}{c}{No. cases /}   & \\multicolumn{1}{c}{Attack}   & \\multicolumn{2}{c}{RR}                     & \\multicolumn{1}{c}{P-value}   & \\multicolumn{1}{c}{Overall P-}       \\\\ 
           \\multicolumn{1}{l}{Immunologic Marker}            & \\multicolumn{1}{c}{}          & \\multicolumn{1}{c}{No. at-risk**} & \\multicolumn{1}{c}{rate}   & \\multicolumn{1}{c}{Pt. Est.} & \\multicolumn{1}{c}{95\\% CI} & \\multicolumn{1}{c}{(2-sided)} & \\multicolumn{1}{c}{value***} \\\\ 
           \\hline\n 
        ")
  }
  
  # use longtable because this table could be long, e.g. in hvtn705second
  mytex(tab[1:(nrow(tab)),], file.name="CoR_univariable_poisson_cat_pretty_"%.%fname.suffix, align="c", include.colnames = F, save2input.only=T, input.foldername=save.results.to,
        col.headers=header,        
        longtable=T, 
        label=paste0("tab:CoR_univariable_poisson_cat_pretty_", study_name), 
        caption.placement = "top", 
        caption=paste0("Inference for ", DayPrefix, tpeak, " antibody marker covariate-adjusted correlates of risk of ", config.cor$txt.endpoint, " in the vaccine group: Hazard ratios for Middle vs. Upper tertile vs. Lower tertile*"))
  
  
  tab.nop12=cbind(
      rep(c("Lower","Middle","Upper"), length(p)/3), 
      paste0(nevents, "/", format(natrisk, big.mark=",",digit=0, scientific=F)), 
      formatDouble(nevents/natrisk, digit=4, remove.leading0=F),
      est, ci, p, overall.p.0
  )
  rownames(tab.nop12)=c(rbind(all.markers.names.short, "", ""))
  
  
}

print(date())
print("cor_poisson run time: "%.%format(Sys.time()-time.start, digits=1))
