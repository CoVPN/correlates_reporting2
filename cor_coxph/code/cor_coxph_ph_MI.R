# multiple imputation version of cor_coxph

###################################################################################################
if(verbose) print("Regression for continuous markers")
# Report point and 95\% confidence interval estimates for the hazard ratio per 10-fold change in the antibody marker, 
# for the entire baseline negative vaccine cohort

fits=list()
fits.scaled=list()
for (i in 1:2) { # 1: not scaled, 2: scaled
  for (a in all.markers) {
    if(i==1) {
      f= update(form.0, as.formula(paste0("~.+", a)))
    } else {
      f= update(form.0, as.formula(paste0("~.+scale(", a, ")")))
    }
    
    models = lapply(1:10, function (imp) {
      # imp=1
      if (TRIAL=="janssen_partA_VL") {
        dat.vac$EventIndOfInterest = ifelse(dat.vac$EventIndPrimary==1 & dat.vac[["seq1.variant.hotdeck"%.%imp]]==variant, 1, 0)
      } else if (TRIAL=="vat08_combined") {
        dat.vac$EventIndOfInterest  = dat.vac[[paste0("EventIndOmicronD",tpeak,"hotdeck",imp)]]
        dat.vac$EventTimeOfInterest = dat.vac[[paste0("EventTimeOmicronD",tpeak,"hotdeck",imp)]]
      } else stop('wrong TRIAL: '%.%TRIAL)
      
      design.vac<-twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dat.vac)
      svycoxph(f, design=design.vac) 
    })
    betas<-MIextract(models, fun=coef)
    vars<-MIextract(models, fun=vcov)
    res<-summary(MIcombine(betas,vars)) # MIcombine prints the results, there is no way to silent it
    
    if (i==1) {
      fits[[a]]=res
      
      # save for forest plots
      if (TRIAL=='janssen_partA_VL') {
        cox.df=rbind(cox.df, list(
          region = region, 
          variant = variant, 
          assay = a, 
          est = exp(res[nrow(res),"results"]),
          lb =exp(res[nrow(res),"(lower"]),
          ub = exp(res[nrow(res),"upper)"])
        ))
      }

    } else {
      fits.scaled[[a]]=res
    }
  }
}

# # sanity check
# i=1
# a="Day29pseudoneutid50"
# # remove cases with missing variants info and cases with competing types
# dat.tmp=dat.vac
# dat.tmp=subset(dat.tmp, !(EventIndPrimary==1 & (is.na(seq1.variant) | seq1.variant!=variant)))
# design.vac<-twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dat.tmp)
# svycoxph(Surv(EventTimePrimaryD29, EventIndPrimary) ~ risk_score + Day29pseudoneutid50, design=design.vac) 


# make continuous markers tables
# one for per 10 fold inc and one for per SD increase

for (i in 1:2) { # 1: not scaled, 2: scaled
  # remove missInfo and cast as matrix to get a numeric matrix
  res=as.matrix(mysapply(if(i==1) fits else fits.scaled, function (fit) as.matrix(subset(fit, select=-missInfo))[nrow(fit),]))
  # exp HR
  tab.1=cbind(
    # est
    formatDouble(exp(res[,1]), 2, remove.leading0=F),
    # CI
    glue("({formatDouble(exp(res[,'(lower']), 2, remove.leading0=F)}-{formatDouble(exp(res[,'upper)']), 2, remove.leading0=F)})"),
    # compute p value based on Gaussian
    formatDouble(2*pnorm(abs(res[,1])/res[,"se"], lower.tail=F), 3, remove.leading0=F)
  )
  rownames(tab.1)=all.markers.names.short
  tab.1
  
  mytex(tab.1, file.name=paste0("CoR_univariable_svycoxph_pretty_",ifelse(i==1,"","scaled_"),fname.suffix), align="c", include.colnames = F, save2input.only=T, input.foldername=save.results.to,
        col.headers=paste0("\\hline\n 
         \\multicolumn{1}{l}{", toTitleCase(study_name), "} & \\multicolumn{2}{c}{HR per ",ifelse(i==1,"10-fold","SD")," incr.}                     & \\multicolumn{1}{c}{P-value}   \\\\ 
         \\multicolumn{1}{l}{Immunologic Marker}            & \\multicolumn{1}{c}{Pt. Est.} & \\multicolumn{1}{c}{95\\% CI} & \\multicolumn{1}{c}{(2-sided)}  \\\\ 
         \\hline\n 
    "),
        longtable=T, 
        label=paste0("tab:CoR_univariable_svycoxph_pretty",ifelse(i==1,"","_scaled")), 
        caption.placement = "top", 
        caption=paste0("Inference for Day ", tpeak, "antibody marker covariate-adjusted correlates of risk of ", 
                       config.cor$txt.endpoint, 
                       " in the vaccine group: Hazard ratios per ",ifelse(i==1,"10-fold","SD")," increment in the marker*")
  )
}



###################################################################################################
if(verbose) print("regression for trichotomized markers")

marker.levels = sapply(all.markers, function(a) length(table(dat.vac[[a%.%"cat"]])))

fits.tri=list()
overall.p.tri=c()
for (a in all.markers) {
  if(verbose) myprint(a)
  f= update(form.0, as.formula(paste0("~.+", a, "cat")))
  
  models = lapply(1:10, function (imp) {
    # imp=1
    if (TRIAL=="janssen_partA_VL") {
      dat.vac$EventIndOfInterest = ifelse(dat.vac$EventIndPrimary==1 & dat.vac[["seq1.variant.hotdeck"%.%imp]]==variant, 1, 0)
    } else if (TRIAL=="vat08_combined") {
      dat.vac$EventIndOfInterest  = dat.vac[[paste0("EventIndOmicronD",tpeak,"hotdeck",imp)]]
      dat.vac$EventTimeOfInterest = dat.vac[[paste0("EventTimeOmicronD",tpeak,"hotdeck",imp)]]
    } else stop('wrong TRIAL: '%.%TRIAL)
    
    design.vac<-twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dat.vac)
    svycoxph(f, design=design.vac) 
  })
  
  betas<-MIextract(models, fun=coef)
  vars<-MIextract(models, fun=vcov)
  combined.model = MIcombine(betas,vars)
  fits.tri[[a]]=summary(combined.model)
  
  # get generalized Wald p values
  
  rows=nrow(fits.tri[[a]])- (marker.levels[a]-2):0
  if (length(rows)>1) {
    stat=coef(combined.model)[rows] %*% solve(vcov(combined.model)[rows,rows]) %*% coef(combined.model)[rows]
  } else {
    stat=coef(combined.model)[rows] * 1/(vcov(combined.model)[rows,rows]) * coef(combined.model)[rows]
  }
  overall.p.tri=c(overall.p.tri, pchisq(stat, length(rows), lower.tail = FALSE))
}
names(overall.p.tri) = all.markers

overall.p.0=formatDouble(c(rbind(overall.p.tri, NA,NA)), digits=3, remove.leading0 = F);   
overall.p.0=sub("0.000","<0.001",overall.p.0)



###################################################################################################
# make trichotomized markers table

get.est=function(a) {
  fit=fits.tri[[a]]
  if (length(fit)==1) return (rep(NA,2))
  res=subset(fit, select=-missInfo)[nrow(fit)-(marker.levels[a]-2):0,1]
  out=formatDouble(exp(res), 2, remove.leading0=F)
  if (length(out)==1) c(NA,out) else out
}
get.ci =function(a) {
  fit=fits.tri[[a]]
  if (length(fit)==1) return (rep(NA,2))
  res=subset(fit, select=-missInfo)[nrow(fit)-(marker.levels[a]-2):0,,drop=F] 
  out=glue("({formatDouble(exp(res[,'(lower']), 2, remove.leading0=F)}-{formatDouble(exp(res[,'upper)']), 2, remove.leading0=F)})")
  if (length(out)==1) c(NA,out) else out
}
get.p  =function(a) {
  fit=fits.tri[[a]]
  if (length(fit)==1) return (rep(NA,2))
  res=subset(fit, select=-missInfo)[nrow(fit)-(marker.levels[a]-2):0,,drop=F] 
  out=formatDouble(2*pnorm(abs(res[,1])/res[,"se"], lower.tail=F), 3, remove.leading0=F)
  if (length(out)==1) c(NA,out) else out
}
get.est(all.markers[1])

# regression parameters
est=c(rbind(1.00,  sapply(all.markers, function (a) get.est(a))))
ci= c(rbind("N/A", sapply(all.markers, function (a) get.ci (a))))
p=  c(rbind("N/A", sapply(all.markers, function (a) get.p  (a))))

tab=cbind(
  rep(c("Lower","Middle","Upper"), length(p)/3), 
  est, ci, p, overall.p.0
)
tmp=rbind(all.markers.names.short, "", "")
rownames(tab)=c(tmp)
tab

# use longtable because this table could be long, e.g. in hvtn705second
mytex(tab[1:(nrow(tab)),], file.name=paste0("CoR_univariable_svycoxph_cat_pretty_",fname.suffix), align="c", include.colnames = F, save2input.only=T, input.foldername=save.results.to,
      col.headers=paste0("\\hline\n 
         \\multicolumn{1}{l}{", toTitleCase(study_name), "} & \\multicolumn{1}{c}{Tertile}   & \\multicolumn{2}{c}{Haz. Ratio}                     & \\multicolumn{1}{c}{P-value}   & \\multicolumn{1}{c}{Overall P-}       \\\\ 
         \\multicolumn{1}{l}{Immunologic Marker}            & \\multicolumn{1}{c}{}          & \\multicolumn{1}{c}{Pt. Est.} & \\multicolumn{1}{c}{95\\% CI} & \\multicolumn{1}{c}{(2-sided)} & \\multicolumn{1}{c}{value***}\\\\ 
         \\hline\n 
    "),       
      longtable=T, 
      label=paste0("tab:CoR_univariable_svycoxph_cat_pretty_", study_name), 
      caption.placement = "top", 
      caption=paste0("Inference for Day ", tpeak, "antibody marker covariate-adjusted correlates of risk of ", 
                     config.cor$txt.endpoint, 
                     " in the vaccine group: Hazard ratios for Middle vs. Upper tertile vs. Lower tertile*")
)




###################################################################################################
# multivariate_assays models

if (!is.null(config$multivariate_assays)) {
  if(verbose) print("Multiple regression")
  
  for (a in config$multivariate_assays) {
    for (i in 1:2) {
      # 1: per SD; 2: per 10-fold
      a.tmp=a
      aa=trim(strsplit(a, "\\+")[[1]])
      for (x in aa[!contain(aa, "\\*")]) {
        # replace every variable with Day210x, with or without scale
        a.tmp=gsub(x, paste0(if(i==1) "scale","(Day",tpeak,x,")"), a.tmp) 
      }
      f= update(form.0, as.formula(paste0("~.+", a.tmp)))
      
      ## todo
      models = lapply(1:10, function (imp) {
        dat.vac$EventIndOfInterest = ifelse(dat.vac$EventIndPrimary==1 & dat.vac[["seq1.variant.hotdeck"%.%imp]]==variant, 1, 0)
        design.vac<-twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dat.vac)
        svycoxph(f, design=design.vac) 
      })
      betas<-MIextract(models, fun=coef)
      vars<-MIextract(models, fun=vcov)
      combined.model = MIcombine(betas,vars)
      res<-summary(combined.model) # MIcombine prints the results, there is no way to silent it
      
      est=formatDouble(exp(res[,1]), 2, remove.leading0=F)
      ci= glue("({formatDouble(exp(res[,'(lower']), 2, remove.leading0=F)}-{formatDouble(exp(res[,'upper)']), 2, remove.leading0=F)})")
      est = paste0(est, " ", ci)
      p=  formatDouble(2*pnorm(abs(res[,1])/res[,"se"], lower.tail=F), 3, remove.leading0=F)
      
      # get generalized Wald p values
      rows=length(betas[[1]]) - length(aa):1 + 1
      stat=coef(combined.model)[rows] %*% solve(vcov(combined.model)[rows,rows]) %*% coef(combined.model)[rows]
      p.gwald=pchisq(stat, length(rows), lower.tail = FALSE)
      
      tab=cbind(est, p)[rows,,drop=F]
      tmp=match(aa, colnames(labels.axis))
      tmp[is.na(tmp)]=1 # otherwise, labels.axis["Day"%.%tpeak, tmp] would throw an error when tmp is NA
      rownames(tab)=ifelse(aa %in% colnames(labels.axis), labels.axis["Day"%.%tpeak, tmp], aa)
      colnames(tab)=c(paste0("HR per ",ifelse(i==1,"sd","10 fold")," incr."), "P value")
      tab
      tab=rbind(tab, "Generalized Wald Test"=c("", formatDouble(p.gwald,3, remove.leading0 = F)))

      mytex(tab, file.name=paste0("CoR_multivariable_svycoxph_pretty", match(a, config$multivariate_assays), if(i==2) "_per10fold",fname.suffix), align="c", include.colnames = T, save2input.only=T, 
            input.foldername=save.results.to)
    }
  }
  
}



write(NA, file=paste0(save.results.to, "permutation_replicates_"%.%study_name))     # so the rmd file can compile
