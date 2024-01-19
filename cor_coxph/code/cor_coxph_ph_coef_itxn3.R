# can be updated to deal with two phase

# for datasets where markers are measured for everyone in the cohort

# for interaction between a continuous marker and a trichotomized marker
# make two tables: 
# one showing the interaction terms and a generalized wald test p value
# two showing the effects of the continuous marker in each level of the trichotomized marker

# control whether fwer and q values are shown in tables
if (is.null(show.q)) show.q=T

###################################################################################################
if(verbose) print("Regression for continuous markers")

# Report point and 95\% confidence interval estimates for the hazard ratio per 10-fold change in the antibody marker, 
# for the entire baseline negative vaccine cohort


fits=list()
p.gwalds=c()
# for making a table of effects of continuous variables in the 3 levels
est.2=NULL
ci.2=NULL
p.2=NULL
for (a in all.markers) {
    f = update(form.0, as.formula(paste0("~.+", a)))
    
    if (use.svy) {
      fits[[a]]=svycoxph(f, design=design.dat) 
    } else {
      fits[[a]]=coxph(f, dat) 
    }
    
    fit=fits[[a]]
    
    if (use.svy) {
      # generalized Wald test for interaction terms
      var.ind=length(coef(fit))-1:0
      stat=coef(fit)[var.ind] %*% solve(vcov(fit)[var.ind,var.ind]) %*% coef(fit)[var.ind] 
      p.gwald=pchisq(stat, length(var.ind), lower.tail = FALSE); p.gwald
      p.gwalds=c(p.gwalds, p.gwald)
      
    } else {
      # likelihood ratio test if we are not doing two phase
      # still use p.gwalds to hold results
      f.1 = update(form.0, as.formula(paste0("~.+", sub("\\*","\\+",a))))
      p.likra= anova(fit, coxph(f.1, dat))$P[2]
      p.gwalds=c(p.gwalds, p.likra)
    
    }

    # get effects in l/m/h
    res=sapply(1:3, function (i) {
      if (i==1) var.ind=length(coef(fit))-c(2)
      if (i==2) var.ind=length(coef(fit))-c(2,1)
      if (i==3) var.ind=length(coef(fit))-c(2,0)
      ef=rep(1,length(var.ind)) %*% coef(fit)[var.ind];
      sd=sqrt(rep(1,length(var.ind)) %*% vcov(fit)[var.ind,var.ind] %*% rep(1,length(var.ind))); 
      p = 2*(1-pnorm(abs(ef)/sd))
      lb=formatDouble(exp(ef-sd*1.96), 2, remove.leading0=F) 
      ub=formatDouble(exp(ef+sd*1.96), 2, remove.leading0=F) 
      ci = "(" %.% lb %.% ", " %.% ub %.% ")" 
      list(exp(ef), ci, p)
    })
    
    est.2=cbind(est.2, unlist(res[1,]))
    ci.2=cbind(ci.2, unlist(res[2,]))
    p.2=cbind(p.2, unlist(res[3,]))
}


natrisk=nrow(dat)
nevents=sum(dat$yy==1)

# make a table for interaction terms, plus p value for generalized wald

rows=length(coef(fits[[1]]))
rows=(rows-1):rows
# robust=F b/c not an option to use robust=T for coxph, but it is a required argument for getFormattedSummary
est=getFormattedSummary(fits, exp=T, robust=F, rows=rows, type=1)
ci= getFormattedSummary(fits, exp=T, robust=F, rows=rows, type=13)
p=  getFormattedSummary(fits, exp=T, robust=F, rows=rows, type=10)

tab.1=cbind(paste0(nevents, "/", format(natrisk, big.mark=",")), 
            est[1,], ci[1,], p[1,],
            est[2,], ci[2,], p[2,],
          formatDouble(p.gwalds, 3, remove.leading0 = F))
rownames(tab.1)=all.markers.names.short
tab.1

header=paste0("\\hline\n 
       \\multicolumn{1}{l}{", '', "} & \\multicolumn{1}{c}{}                         & \\multicolumn{3}{c}{",col1,"}   & \\multicolumn{3}{c}{",col2,"}   & \\multicolumn{1}{c}{",ifelse(use.svy,"Generalized Wald","Lik Ratio"),"}    \\\\ 
       \\multicolumn{1}{l}{", '', "} & \\multicolumn{1}{c}{No. cases /}              & \\multicolumn{2}{c}{Ratio of HRs}                      & \\multicolumn{1}{c}{P-value}   & \\multicolumn{2}{c}{Ratio of HRs}                     & \\multicolumn{1}{c}{P-value}   & \\multicolumn{1}{c}{P-value}    \\\\ 
       \\multicolumn{1}{l}{Immunologic Marker}  & \\multicolumn{1}{c}{No. at-risk**} & \\multicolumn{1}{c}{Pt. Est.} & \\multicolumn{1}{c}{95\\% CI}  & \\multicolumn{1}{c}{(2-sided)} & \\multicolumn{1}{c}{Pt. Est.} & \\multicolumn{1}{c}{95\\% CI} & \\multicolumn{1}{c}{(2-sided)} & \\multicolumn{1}{c}{}   \\\\ 
       \\hline\n 
  ")

# keep svycoxph in the name so that we can reuse rmd file
mytex(tab.1, file.name="CoR_univariable_svycoxph_pretty_"%.%fname.suffix, align="c", include.colnames = F, save2input.only=T, input.foldername=save.results.to,
    col.headers=header,
    longtable=T, 
    label=paste0("tab:CoR_univariable_svycoxph_pretty"), 
    caption.placement = "top", 
    caption=paste0("Inference for Day ", tpeak, " antibody marker covariate-adjusted correlates of risk of ", config.cor$txt.endpoint, ": Interaction terms and generalized Wald test P values over interaction terms*")
)


# table of effects sizes in different strata
est.2=formatDouble(est.2, digits=2, remove.leading0 = F)
p.2=formatDouble(p.2, digits=3, remove.leading0 = F)
tab.2=cbind(paste0(nevents, "/", format(natrisk, big.mark=",")), 
            est.2[1,], ci.2[1,], p.2[1,],
            est.2[2,], ci.2[2,], p.2[2,],
            est.2[3,], ci.2[3,], p.2[3,])
rownames(tab.2)=all.markers.names.short
tab.2

header=paste0("\\hline\n 
       \\multicolumn{1}{l}{", '', "} & \\multicolumn{1}{c}{}              & \\multicolumn{3}{c}{",col3,"}   & \\multicolumn{3}{c}{",col4,"}   & \\multicolumn{3}{c}{",col5,"}    \\\\ 
       \\multicolumn{1}{l}{", '', "} & \\multicolumn{1}{c}{No. cases /}              & \\multicolumn{2}{c}{HR per 10-fold incr.} & \\multicolumn{1}{c}{P-value}   & \\multicolumn{2}{c}{HR per 10-fold incr.} & \\multicolumn{1}{c}{P-value}   & \\multicolumn{2}{c}{HR per 10-fold incr.} & \\multicolumn{1}{c}{P-value}    \\\\ 
       \\multicolumn{1}{l}{Immunologic Marker}  & \\multicolumn{1}{c}{No. at-risk**} & \\multicolumn{1}{c}{Pt. Est.} & \\multicolumn{1}{c}{95\\% CI} & \\multicolumn{1}{c}{(2-sided)} & \\multicolumn{1}{c}{Pt. Est.} & \\multicolumn{1}{c}{95\\% CI} & \\multicolumn{1}{c}{(2-sided)} & \\multicolumn{1}{c}{Pt. Est.} & \\multicolumn{1}{c}{95\\% CI} & \\multicolumn{1}{c}{(2-sided)}  \\\\ 
       \\hline\n 
  ")


# keep svycoxph in the name so that we can reuse rmd file
mytex(tab.2, file.name="CoR_univariable_svycoxph_pretty_"%.%fname.suffix%.%"_2", align="c", include.colnames = F, save2input.only=T, input.foldername=save.results.to,
      col.headers=header,
      longtable=T, 
      label=paste0("tab:CoR_univariable_svycoxph_pretty"), 
      caption.placement = "top", 
      caption=paste0("Inference for Day ", tpeak, " antibody marker covariate-adjusted correlates of risk of ", config.cor$txt.endpoint, ": Hazard ratios per 10-fold increment in the marker in separate subgroups*")
)
