# Output:
# two tables for continuous markers 
# one table for discrete markers

# possible errors:
# Error in .local(x, i, j, ..., value) : not-yet-implemented 'Matrix[<-' method
#   likely caused by missing values


cor_coxph_coef_1_mi = function(
  form.0,
  dat, # either design when tps is T, or data frame otherwise
  fname.suffix, #used in the file names to save results
  save.results.to,
  config,
  config.cor,
  markers,
  markers.names.short,
  
  dat.pla.seroneg = NULL,
  show.q=TRUE, # whether to show fwer and q values in tables
  
  forestplot.markers=list(1:length(markers)), # make forestplot only for a list of subsets of markers
  for.title="",
  
  verbose=FALSE
) {
  
  
  if(verbose) print(paste0("Running cor_coxph_coef_1_mi: ", fname.suffix))
  
  ###################################################################################################
  if(verbose) print("Regression for continuous markers")
  # Report point and 95\% confidence interval estimates for the hazard ratio per 10-fold change in the antibody marker, 
  # for the entire baseline negative vaccine cohort
  
  fits=list()
  fits.scaled=list()
  for (i in 1:2) { # 1: not scaled, 2: scaled
    for (a in markers) {
      myprint(a)
      
      models=lapply(1:10, function(imp) {
        # when this script becomes a function, mclapply runs into error:   unable to fork, possible reason: Cannot allocate memory
        # models=mclapply(1:10, mc.cores = 10, FUN=function(imp) {
          # imp=1
        
        if (TRIAL=="janssen_partA_VL" | 
            TRIAL=="vat08_combined" & endsWith(COR,"st1.nAb.batch0and1")) {
          # for janssen_partA_VL, form.0 is not a list because this is for Cox regression and not risk
          f = update(form.0, 
                     as.formula(paste0("~.+", if(i==2) "scale", "(", a, "_"%.%imp, ")"))
                     )
        } else if (TRIAL=="vat08_combined") {
          f = update(form.0, 
                     as.formula(paste0("~.+", if(i==2) "scale", "(", a, ")"))
                     )
        } else {
          stop("wrong TRIAL")
        }
        
        
        # set event indicator and time
        if (TRIAL=="janssen_partA_VL") {
          dat$EventIndOfInterest = ifelse(dat$EventIndPrimary==1 & dat[["seq1.variant.hotdeck"%.%imp]]==variant, 1, 0)
        } else if (TRIAL=="vat08_combined") {
          # not competing risk, imputation only
          dat$EventIndOfInterest  = dat[[config.cor$EventIndPrimary  %.% imp]]
          dat$EventTimeOfInterest = dat[[config.cor$EventTimePrimary %.% imp]]
        } else stop('wrong TRIAL: '%.%TRIAL)
        
        
        if (TRIAL=="janssen_partA_VL" & a %in% c("Day29bindSpike","Day29pseudoneutid50")) {
          dat$ph2a = dat$ph2.D29
        } else {
          dat$ph2a = dat$ph2
        }
        
        design.vac<-twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2a, data=dat)
        svycoxph(f, design=design.vac) 
      
      })
      
      betas<-MIextract(models, fun=coef)
      vars<- MIextract(models, fun=vcov)
      
      capture.output({
        res<-summary(MIcombine(betas,vars)) # MIcombine prints the results, there is no way to silent it
      }, type="output") # type = message captures stderr, type=output is for stdout
      
      if (i==1) {
        fits[[a]]=res
      } else {
        fits.scaled[[a]]=res
      }
    }
  }
  

  # make continuous markers tables
  # one for per 10 fold inc and one for per SD increase
  for (i in 1:2) { # 1: not scaled, 2: scaled
    
    # remove missInfo and cast as matrix to get a numeric matrix
    res=as.matrix(mysapply(if(i==1) fits else fits.scaled, function (fit) as.matrix(fit[,names(fit)!="missInfo"])[nrow(fit),]))
    
    
    # make table
    tab.1=cbind(
      # est
      formatDouble(exp(res[,1]), 2, remove.leading0=F),
      # CI
      paste0("(", formatDouble(exp(res[,'(lower']), 2, remove.leading0=F), "-", formatDouble(exp(res[,'upper)']), 2, remove.leading0=F), ")"),
      # compute p value based on Gaussian
      formatDouble(2*pnorm(abs(res[,1])/res[,"se"], lower.tail=F), 3, remove.leading0=F)
    )
    rownames(tab.1)=markers.names.short
    tab.1
    
    mytex(tab.1, file.name=paste0("CoR_univariable_svycoxph_pretty_",ifelse(i==1,"","scaled_"),fname.suffix), align="c", include.colnames = F, save2input.only=T, input.foldername=save.results.to,
          col.headers=paste0("\\hline\n 
           \\multicolumn{1}{l}{} & \\multicolumn{2}{c}{HR per ",ifelse(i==1,"10-fold","SD")," incr.}                     & \\multicolumn{1}{c}{P-value}   \\\\ 
           \\multicolumn{1}{l}{Immunologic Marker}            & \\multicolumn{1}{c}{Pt. Est.} & \\multicolumn{1}{c}{95\\% CI} & \\multicolumn{1}{c}{(2-sided)}  \\\\ 
           \\hline\n 
      "),
          longtable=T, 
          label=paste0("tab:CoR_univariable_svycoxph_pretty",ifelse(i==1,"","_scaled")), 
          caption.placement = "top", 
          caption=paste0("Inference for Day ", config.cor$tpeak, " antibody marker covariate-adjusted correlates of risk of ", 
                         config.cor$txt.endpoint, 
                         " in the vaccine group: Hazard ratios per ",ifelse(i==1,"10-fold","SD")," increment in the marker*")
    )
    
    # make forest plots
    
    if (!is.list(forestplot.markers)) forestplot.markers=list(forestplot.markers)
    
    for (iM in 1:length(forestplot.markers)) {
      
      est.ci = rbind(exp(res[forestplot.markers[[iM]], 1]), 
                     exp(res[forestplot.markers[[iM]], '(lower']), 
                     exp(res[forestplot.markers[[iM]], 'upper)']), 
                     2*pnorm(abs(res[,1])/res[,"se"], lower.tail=F)[forestplot.markers[[iM]]]
                     )
      colnames(est.ci)=markers.names.short[forestplot.markers[[iM]]]
      est.ci
      
      nevents=rep(NA, ncol(est.ci))
      
      # make two versions, one log and one antilog
      
      mypdf(onefile=F, width=10,height=4, file=paste0(save.results.to, "hr_forest_", ifelse(i==1,"","scaled_"), fname.suffix, if (iM>1) iM)) 
        theforestplot(point.estimates=est.ci[1,], lower.bounds=est.ci[2,], upper.bounds=est.ci[3,], group=colnames(est.ci), 
                    nEvents=nevents, 
                    p.values=formatDouble(est.ci[4,], 3, remove.leading0=F), 
                    decimal.places=2, graphwidth=unit(120, "mm"), fontsize=1.2, 
                    table.labels = c("", "  HR (95% CI)",""), 
                    xlog=F,
                    title=for.title, 
                    x.ticks = get.forestplot.ticks(est.ci, forestplot.xlog=F)  # controls the limit
      )
      dev.off()
      
      mypdf(onefile=F, width=10,height=4, file=paste0(save.results.to, "hr_forest_log_", ifelse(i==1,"","scaled_"), fname.suffix, if (iM>1) iM)) 
        theforestplot(point.estimates=est.ci[1,], lower.bounds=est.ci[2,], upper.bounds=est.ci[3,], group=colnames(est.ci), 
                    nEvents=nevents, p.values=formatDouble(est.ci[4,], 3, remove.leading0=F), 
                    decimal.places=2, graphwidth=unit(120, "mm"), fontsize=1.2, 
                    table.labels = c("", "  HR (95% CI)",""), 
                    xlog=T,
                    title=for.title, 
                    x.ticks = get.forestplot.ticks(est.ci, forestplot.xlog=T)  # controls the limit
      )
      dev.off()
    
    }
    
    
  }
  
  
  
  ###################################################################################################
  if(verbose) print("regression for trichotomized markers")
  
  marker.levels = sapply(markers, function(a) length(table(dat[[a%.%"cat"]]))); marker.levels
  
  fits.tri=list()
  overall.p.tri=c()
  for (a in markers) {
    if(verbose) myprint(a)
    
    models=lapply(1:10, function(imp) {
      # models=mclapply(1:10, mc.cores = 10, FUN=function(imp) {
        # imp=1
      f = update(form.0, as.formula(paste0("~.+", a, if(TRIAL=="janssen_partA_VL" |
                                                        TRIAL=="vat08_combined" & endsWith(COR,"st1.nAb.batch0and1")) "_"%.%imp, "cat"))); f
  
      if (TRIAL=="janssen_partA_VL") {
        dat$EventIndOfInterest = ifelse(dat$EventIndPrimary==1 & dat[["seq1.variant.hotdeck"%.%imp]]==variant, 1, 0)
        
      } else if (TRIAL=="vat08_combined") {
        dat$EventIndOfInterest  = dat[[config.cor$EventIndPrimary %.% imp]]
        dat$EventTimeOfInterest = dat[[config.cor$EventTimePrimary %.% imp]]
        
      } else stop('wrong TRIAL: '%.%TRIAL)
      
      
      if (TRIAL=="janssen_partA_VL" & a %in% c("Day29bindSpike","Day29pseudoneutid50")) {
        dat$ph2a = dat$ph2.D29
      
      } else {
        dat$ph2a = dat$ph2
      }
      
      design.vac<-twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2a, data=dat)
      svycoxph(f, design=design.vac) 
    })
    
    betas<-MIextract(models, fun=coef)
    vars<-MIextract(models, fun=vcov)
    
    
    capture.output({
      combined.model<-MIcombine(betas,vars) # MIcombine prints the results, there is no way to silent it
      fits.tri[[a]]=summary(combined.model)
    }, type="output") # type = message captures stderr, type=output is for stdout

    
    # get generalized Wald p values
    
    rows=nrow(fits.tri[[a]])- (marker.levels[a]-2):0
    if (length(rows)>1) {
      stat=coef(combined.model)[rows] %*% solve(vcov(combined.model)[rows,rows]) %*% coef(combined.model)[rows]
    } else {
      stat=coef(combined.model)[rows] * 1/(vcov(combined.model)[rows,rows]) * coef(combined.model)[rows]
    }
    overall.p.tri=c(overall.p.tri, pchisq(stat, length(rows), lower.tail = FALSE))
  }
  names(overall.p.tri) = markers
  
  overall.p.0=formatDouble(c(rbind(overall.p.tri, NA,NA)), digits=3, remove.leading0 = F);   
  overall.p.0=sub("0.000","<0.001",overall.p.0)
  
  save(overall.p.tri, file=paste0(save.results.to, "svycoxph_cat_overall_pvalues_",fname.suffix,".Rdata"))
        
  
  ###################################################################################################
  # make trichotomized markers table
  
  get.est=function(a) {
    fit=fits.tri[[a]]
    if (length(fit)==1) return (rep(NA,2))
    res=fit[,names(fit)!="missInfo"][nrow(fit)-(marker.levels[a]-2):0,1]
    out=formatDouble(exp(res), 2, remove.leading0=F)
    if (length(out)==1) c(NA,out) else out
  }
  get.ci =function(a) {
    fit=fits.tri[[a]]
    if (length(fit)==1) return (rep(NA,2))
    res=fit[,names(fit)!="missInfo"][nrow(fit)-(marker.levels[a]-2):0,,drop=F] 
    out=paste0("(", formatDouble(exp(res[,'(lower']), 2, remove.leading0=F), "-", formatDouble(exp(res[,'upper)']), 2, remove.leading0=F), ")")
    if (length(out)==1) c(NA,out) else out
  }
  get.p  =function(a) {
    fit=fits.tri[[a]]
    if (length(fit)==1) return (rep(NA,2))
    res=fit[,names(fit)!="missInfo"][nrow(fit)-(marker.levels[a]-2):0,,drop=F] 
    out=formatDouble(2*pnorm(abs(res[,1])/res[,"se"], lower.tail=F), 3, remove.leading0=F)
    if (length(out)==1) c(NA,out) else out
  }
  get.est(markers[1])
  
  # regression parameters
  est=c(rbind(1.00,  sapply(markers, function (a) get.est(a))))
  ci= c(rbind("N/A", sapply(markers, function (a) get.ci (a))))
  p=  c(rbind("N/A", sapply(markers, function (a) get.p  (a))))
  
  tab=cbind(
    rep(c("Lower","Middle","Upper"), length(p)/3), 
    est, ci, p, overall.p.0
  )
  tmp=rbind(markers.names.short, "", "")
  rownames(tab)=c(tmp)
  tab
  
  # use longtable because this table could be long, e.g. in hvtn705second
  mytex(tab[1:(nrow(tab)),], file.name=paste0("CoR_univariable_svycoxph_cat_pretty_",fname.suffix), align="c", include.colnames = F, save2input.only=T, input.foldername=save.results.to,
        col.headers=paste0("\\hline\n 
           \\multicolumn{1}{l}{} & \\multicolumn{1}{c}{Tertile}   & \\multicolumn{2}{c}{Haz. Ratio}                     & \\multicolumn{1}{c}{P-value}   & \\multicolumn{1}{c}{Overall P-}       \\\\ 
           \\multicolumn{1}{l}{Immunologic Marker}            & \\multicolumn{1}{c}{}          & \\multicolumn{1}{c}{Pt. Est.} & \\multicolumn{1}{c}{95\\% CI} & \\multicolumn{1}{c}{(2-sided)} & \\multicolumn{1}{c}{value***}\\\\ 
           \\hline\n 
      "),       
        longtable=T, 
        label=paste0("tab:CoR_univariable_svycoxph_cat_pretty_", config$study_name), 
        caption.placement = "top", 
        caption=paste0("Inference for Day ", config.cor$tpeak, " antibody marker covariate-adjusted correlates of risk of ", 
                       config.cor$txt.endpoint, 
                       " in the vaccine group: Hazard ratios for Middle vs. Upper tertile vs. Lower tertile*")
  )
  
  
  
  
  ###################################################################################################
  # multivariate_assays models
  
  if (!is.null(config$multivariate_assays)) {
  
  for (a in config$multivariate_assays) {
    aa=trim(strsplit(a, "\\+")[[1]])
    
    for (i in 1:2) { # 1: per SD; 2: per 10-fold
      
      # the function also depends on DayPrefix, tpeak, TRIAL, aa, a
      get.f=function(i, imp) {
        a.tmp=a
        for (x in aa[!contain(aa, "\\*")]) {
          # replace x with, e.g., Day210x, with scale if i==1 and without scale if i==2
          # marker may also be multiple imputed
          
          # changed from gsub to sub at some point, which may cause a bug
          # the order of items in aa is important, e.g. 
          # a="bindSpike+bindSpike_B.1.351" works, but a="bindSpike_B.1.351+bindSpike" does not
          # since bindSpike_B.1.351 contains bindSpike
          
          a.tmp=sub(x, paste0(if(i==1) "scale",  "(", DayPrefix, config.cor$tpeak, x, if(TRIAL=="janssen_partA_VL") "_"%.%imp, ")"), a.tmp) 
        }
        
        f = update(form.0, as.formula(paste0("~.+", a.tmp)))
        f
      }
      
      models = lapply(1:10, function (imp) {
        dat$EventIndOfInterest = ifelse(dat$EventIndPrimary==1 & dat[["seq1.variant.hotdeck"%.%imp]]==variant, 1, 0)
        
        if (TRIAL=="janssen_partA_VL" & a=="bindSpike+pseudoneutid50") {
          dat$ph2a = dat$ph2.D29
        } else {
          dat$ph2a = dat$ph2
        }
        design.vac<-twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2a, data=dat)
        svycoxph(get.f(i, imp), design=design.vac) 
      })
      betas<-MIextract(models, fun=coef)
      vars<-MIextract(models, fun=vcov)

      capture.output({
        combined.model<-MIcombine(betas,vars) # MIcombine prints the results, there is no way to silent it
        res<-summary(combined.model) # MIcombine prints the results, there is no way to silent it
      }, type="output") # type = message captures stderr, type=output is for stdout

      
      est=formatDouble(exp(res[,1]), 2, remove.leading0=F)
      ci= paste0("(", formatDouble(exp(res[,'(lower']), 2, remove.leading0=F), "-", formatDouble(exp(res[,'upper)']), 2, remove.leading0=F), ")")
      est = paste0(est, " ", ci)
      p=  formatDouble(2*pnorm(abs(res[,1])/res[,"se"], lower.tail=F), 3, remove.leading0=F)
      
      # get generalized Wald p values
      rows=length(betas[[1]]) - length(aa):1 + 1
      stat=coef(combined.model)[rows] %*% solve(vcov(combined.model)[rows,rows]) %*% coef(combined.model)[rows]
      p.gwald=pchisq(stat, length(rows), lower.tail = FALSE)
      
      tab=cbind(est, p)[rows,,drop=F]
      ## set row names
      # tmp=match(aa, colnames(labels.axis))
      # tmp[is.na(tmp)]=1 # otherwise, labels.axis["Day"%.%tpeak, tmp] would throw an error when tmp is NA
      # rownames(tab)=ifelse(aa %in% colnames(labels.axis), labels.axis["Day"%.%config.cor$tpeak, tmp], aa)
      colnames(tab)=c(paste0("HR per ",ifelse(i==1,"sd","10 fold")," incr."), "P value")
      tab
      tab=rbind(tab, "Generalized Wald Test"=c("", formatDouble(p.gwald,3, remove.leading0 = F)))
  
      mytex(tab, file.name=paste0("CoR_multivariable_svycoxph_pretty", 
                                  match(a, config$multivariate_assays), 
                                  if(i==2) "_per10fold",
                                  fname.suffix), 
            align="c", include.colnames = T, save2input.only=T, input.foldername=save.results.to)
    }
  }
  
  }
  
  
  write(NA, file=paste0(save.results.to, "permutation_replicates_"%.%config$study_name))     # so the rmd file can compile

  
}