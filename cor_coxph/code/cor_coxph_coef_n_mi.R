# Make a table to report HR, CI and P value of the last nCoef coefficients. One model per row
# nCoef=1, used for model with one continuous effect of interest
# nCoef=2, used for model with two continuous effects of interest
# nCoef=3, used for interaction between two continuous models

# no multitesting b/c we usually do multitesting across continuous and discrete variables

# all.markers, which defines the formula and can be flexible, e.g., Bmarkercat * scale(Day15marker,scale=F)

cor_coxph_coef_n_mi = function(
  form.0,
  dat, 
  fname.suffix, #used in the file names to save results
  save.results.to,
  config,
  config.cor,
  
  all.markers,
  all.markers.names.short, 
  
  nCoef,
  col.headers, # length n coef, e.g., c("TrtA Moderna~Pfizer", "center(D15)", "TrtA:center(D15)")
  
  verbose=FALSE
) {
  
  
  if(verbose) {
    cat("Running cor_coxph_coef_n_mi \n")
    myprint(col.headers)
  }
  
  # multiple imputation
  fits=list()
  for (a in all.markers) {
  
    models=lapply(1:10, function(imp) {
    # for some unknown reason, mclapply leads to out of memory fault here. change f to as.formula() does not help 
    # models=mclapply(1:10, mc.cores = 10, FUN=function(imp) {
      # imp=1
      
      f = update(form.0, as.formula(paste0("~.+", a, if(TRIAL=="janssen_partA_VL" |
                                                        TRIAL=="vat08_combined" & endsWith(COR,"st1.nAb.batch0and1")) "_"%.%imp))); f

      if (TRIAL %in% c("vat08_combined")) {
        dat$EventIndOfInterest  = dat[[config.cor$EventIndPrimary  %.% imp]]
        dat$EventTimeOfInterest = dat[[config.cor$EventTimePrimary %.% imp]]
      }
      
      design.vac<-twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dat)
      svycoxph(f, design=design.vac) 
      
    })
    
    betas<-MIextract(models, fun=coef)
    vars<- MIextract(models, fun=vcov)
    
    # calling summary on the output of MIcombine prints the results
    # getFormattedSummary calls getFixedEf.MIresult, which uses capture.output, which apparently silences it
    
    fits[[a]]=MIcombine(betas,vars)  
  
  }
  
  rows=length(coef(fits[[1]]))
  rows=(rows-(nCoef-1)):rows
  # robust=tps b/c not an option to use robust=T for coxph, but it is a required argument for getFormattedSummary
  tps=T
  est=getFormattedSummary(fits, exp=T, robust=tps, rows=rows, type=1)
  ci= getFormattedSummary(fits, exp=T, robust=tps, rows=rows, type=7)
  p=  getFormattedSummary(fits, exp=T, robust=tps, rows=rows, type=10)

  
  # # include number of cases in table
  # natrisk=nrow(dat)
  # nevents=sum(dat$yy==1)
  # tab.1=rep(paste0(nevents, "/", format(natrisk, big.mark=",")), ncol(est))
  # header=paste0("\\hline\n 
  #      \\multicolumn{1}{l}{", '', "} & \\multicolumn{1}{c}{}", concatList(paste0("& \\multicolumn{3}{c}{",col.headers,"} ")), "   \\\\ 
  #      \\multicolumn{1}{l}{", '', "} & \\multicolumn{1}{c}{No. cases /}   ", concatList(rep("& \\multicolumn{2}{c}{HR per 10-fold incr.} & \\multicolumn{1}{c}{P-value}",nCoef)), "  \\\\ 
  #      \\multicolumn{1}{l}{Immunologic Marker}  & \\multicolumn{1}{c}{No. at-risk**} ", concatList(rep("& \\multicolumn{1}{c}{Pt. Est.} & \\multicolumn{1}{c}{95\\% CI} & \\multicolumn{1}{c}{}", nCoef)), "  \\\\ 
  #      \\hline\n 
  # ")
  
  # skip number of cases in table
  tab.1=NULL
  header=paste0("\\hline\n 
       \\multicolumn{1}{l}{", '', "}", concatList(paste0("& \\multicolumn{3}{c}{",col.headers,"} ")), "   \\\\ 
       \\multicolumn{1}{l}{", '', "}", concatList(rep("& \\multicolumn{2}{c}{HR per 10-fold incr.} & \\multicolumn{1}{c}{P-value}",nCoef)), "  \\\\ 
       \\multicolumn{1}{l}{Immunologic Marker}", concatList(rep("& \\multicolumn{1}{c}{Pt. Est.} & \\multicolumn{1}{c}{95\\% CI} & \\multicolumn{1}{c}{}", nCoef)), "  \\\\ 
       \\hline\n 
  ")
  
  for (i in 1:nCoef) tab.1=cbind(tab.1, est[i,], ci[i,], p[i,])
  rownames(tab.1)=all.markers.names.short
  tab.1
  mytex(tab.1, file.name="CoR_univariable_svycoxph_pretty_"%.%fname.suffix, align="c", include.colnames = F, save2input.only=T, input.foldername=save.results.to,
        col.headers=header,
        longtable=T, 
        label=paste0("tab:CoR_univariable_svycoxph_pretty"), 
        caption.placement = "top", 
        caption=paste0("Inference for Day ", config.cor$tpeak, " antibody marker covariate-adjusted correlates of risk of ", config.cor$txt.endpoint, " in the ", escape(fname.suffix), " group: Hazard ratios per 10-fold increment in the marker*")
  )
  
  
  # optional
  # code return p values
}

