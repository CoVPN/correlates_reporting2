# multiple imputation version of cor_coxph

###################################################################################################
if(verbose) print("Regression for continuous markers")
# Report point and 95\% confidence interval estimates for the hazard ratio per 10-fold change in the antibody marker, 
# for the entire baseline negative vaccine cohort

fits=list()
for (a in all.markers) {
  f= update(form.0, as.formula(paste0("~.+", a)))
  models=lapply(datasets, function (dataset) {
    design.vacc.seroneg<-twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dataset)
    svycoxph(f, design=design.vacc.seroneg) 
  })
  betas<-MIextract(models, fun=coef)
  vars<-MIextract(models, fun=vcov)
  fits[[a]]<-summary(MIcombine(betas,vars))
}

# scaled marker
fits.scaled=list()
for (a in all.markers) {
  f= update(form.0, as.formula(paste0("~.+scale(", a, ")")))
  models=lapply(datasets, function (dataset) {
    design.vacc.seroneg<-twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dataset)
    svycoxph(f, design=design.vacc.seroneg) 
  })
  betas<-MIextract(models, fun=coef)
  vars<-MIextract(models, fun=vcov)
  fits.scaled[[a]]<-summary(MIcombine(betas,vars))
}


###################################################################################################
# make two continuous markers tables
# one for per 10 fold inc and one for per SD increase

for (i in 1:2) {
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
  
  mytex(tab.1, file.name=paste0("CoR_univariable_svycoxph_pretty_",ifelse(i==1,"","scaled_"),region,"_",variant), align="c", include.colnames = F, save2input.only=T, input.foldername=save.results.to,
        col.headers=paste0("\\hline\n 
         \\multicolumn{1}{l}{", toTitleCase(study_name), "} & \\multicolumn{2}{c}{HR per ",ifelse(i==1,"10-fold","SD")," incr.}                     & \\multicolumn{1}{c}{P-value}   \\\\ 
         \\multicolumn{1}{l}{Immunologic Marker}            & \\multicolumn{1}{c}{Pt. Est.} & \\multicolumn{1}{c}{95\\% CI} & \\multicolumn{1}{c}{(2-sided)}  \\\\ 
         \\hline\n 
    "),
        longtable=T, 
        label=paste0("tab:CoR_univariable_svycoxph_pretty",ifelse(i==1,"","_scaled")), 
        caption.placement = "top", 
        caption=paste0("Inference for Day ", tpeak, "antibody marker covariate-adjusted correlates of risk of ", 
                       variant, "-", config.cor$txt.endpoint, " in ", region,
                       " in the vaccine group: Hazard ratios per ",ifelse(i==1,"10-fold","SD")," increment in the marker*")
  )
}



###################################################################################################
if(verbose) print("regression for trichotomized markers")

fits.tri=list()
overall.p.tri=c()
for (a in all.markers) {
  if(verbose) myprint(a)
  f= update(form.0, as.formula(paste0("~.+", a, "cat")))
  models=lapply(datasets, function (dataset) {
    design.vacc.seroneg<-twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dataset)
    svycoxph(f, design=design.vacc.seroneg) 
  })
  betas<-MIextract(models, fun=coef)
  vars<-MIextract(models, fun=vcov)
  combined.model = MIcombine(betas,vars)
  fits.tri[[a]]=summary(combined.model)
  
  # get generalized Wald p values
  rows=nrow(fits.tri[[1]])-1:0
  stat=coef(combined.model)[rows] %*% solve(vcov(combined.model)[rows,rows]) %*% coef(combined.model)[rows]
  overall.p.tri=c(overall.p.tri, pchisq(stat, length(rows), lower.tail = FALSE))
}
names(overall.p.tri) = all.markers

overall.p.0=formatDouble(c(rbind(overall.p.tri, NA,NA)), digits=3, remove.leading0 = F);   
overall.p.0=sub("0.000","<0.001",overall.p.0)



###################################################################################################
# make trichotomized markers table

get.est=function(fit) {
  res=subset(fit, select=-missInfo)[nrow(fit)-1:0,1]
  formatDouble(exp(res), 2, remove.leading0=F)
}
get.ci =function(fit) {
  res=subset(fit, select=-missInfo)[nrow(fit)-1:0,] 
  glue("({formatDouble(exp(res[,'(lower']), 2, remove.leading0=F)}-{formatDouble(exp(res[,'upper)']), 2, remove.leading0=F)})")
}
get.p  =function(fit) {
  res=subset(fit, select=-missInfo)[nrow(fit)-1:0,] 
  formatDouble(2*pnorm(abs(res[,1])/res[,"se"], lower.tail=F), 3, remove.leading0=F)
}
get.p(fits.tri[[1]])

# regression parameters
est=c(rbind(1.00,  sapply(fits.tri, function (fit) if(length(fit)==1) rep(NA,2) else get.est(fit))))
ci= c(rbind("N/A", sapply(fits.tri, function (fit) if(length(fit)==1) rep(NA,2) else get.ci (fit))))
p=  c(rbind("N/A", sapply(fits.tri, function (fit) if(length(fit)==1) rep(NA,2) else get.p  (fit))))

tab=cbind(
  rep(c("Lower","Middle","Upper"), length(p)/3), 
  est, ci, p, overall.p.0
)
tmp=rbind(all.markers.names.short, "", "")
rownames(tab)=c(tmp)
tab

# use longtable because this table could be long, e.g. in hvtn705second
mytex(tab[1:(nrow(tab)),], file.name=paste0("CoR_univariable_svycoxph_cat_pretty_",region,"_",variant), align="c", include.colnames = F, save2input.only=T, input.foldername=save.results.to,
      col.headers=paste0("\\hline\n 
         \\multicolumn{1}{l}{", toTitleCase(study_name), "} & \\multicolumn{1}{c}{Tertile}   & \\multicolumn{2}{c}{Haz. Ratio}                     & \\multicolumn{1}{c}{P-value}   & \\multicolumn{1}{c}{Overall P-}       \\\\ 
         \\multicolumn{1}{l}{Immunologic Marker}            & \\multicolumn{1}{c}{}          & \\multicolumn{1}{c}{Pt. Est.} & \\multicolumn{1}{c}{95\\% CI} & \\multicolumn{1}{c}{(2-sided)} & \\multicolumn{1}{c}{value***}\\\\ 
         \\hline\n 
    "),       
      longtable=T, 
      label=paste0("tab:CoR_univariable_svycoxph_cat_pretty_", study_name), 
      caption.placement = "top", 
      caption=paste0("Inference for Day ", tpeak, "antibody marker covariate-adjusted correlates of risk of ", 
                     variant, "-", config.cor$txt.endpoint, " in ", region,
                     " in the vaccine group: Hazard ratios for Middle vs. Upper tertile vs. Lower tertile*")
)
