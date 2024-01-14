# for datasets where markers are measured for everyone in the cohort

# control whether fwer and q values are shown in tables
if (is.null(show.q)) show.q=T

###################################################################################################
if(verbose) print("Regression for continuous markers")

# Report point and 95\% confidence interval estimates for the hazard ratio per 10-fold change in the antibody marker, 
# for the entire baseline negative vaccine cohort


fits=list()
for (a in all.markers) {
    f= update(form.0, as.formula(paste0("~.+", a)))
    fits[[a]]=coxph(f, dat) 
}


natrisk=nrow(dat)
nevents=sum(dat$yy==1)

# make pretty table
rows=length(coef(fits[[1]]))
rows=(rows-2):rows
# robust=F b/c not an option to use robust=T for coxph, but it is a required argument for getFormattedSummary
est=getFormattedSummary(fits, exp=T, robust=F, rows=rows, type=1)
ci= getFormattedSummary(fits, exp=T, robust=F, rows=rows, type=13)
p=  getFormattedSummary(fits, exp=T, robust=F, rows=rows, type=10)
est.scaled=getFormattedSummary(fits.scaled, exp=T, robust=F, rows=rows, type=1)
ci.scaled= getFormattedSummary(fits.scaled, exp=T, robust=F, rows=rows, type=13)

pvals.cont = sapply(fits, function(x) {
    tmp=getFixedEf(x)
    p.val.col=which(startsWith(tolower(colnames(tmp)),"p"))
    tmp[nrow(tmp),p.val.col]
})



tab.1=cbind(paste0(nevents, "/", format(natrisk, big.mark=",")), 
            est[1,], ci[1,], p[1,],
            est[2,], ci[2,], p[2,],
            est[3,], ci[3,], p[3,])
rownames(tab.1)=all.markers.names.short
tab.1

header=paste0("\\hline\n 
       \\multicolumn{1}{l}{", '', "} & \\multicolumn{1}{c}{}              & \\multicolumn{3}{c}{",col1,"}   & \\multicolumn{3}{c}{",col2,"}   & \\multicolumn{3}{c}{",col3,"}    \\\\ 
       \\multicolumn{1}{l}{", '', "} & \\multicolumn{1}{c}{No. cases /}              & \\multicolumn{2}{c}{HR per 10-fold incr.} & \\multicolumn{1}{c}{P-value}   & \\multicolumn{2}{c}{HR per 10-fold incr.} & \\multicolumn{1}{c}{P-value}   & \\multicolumn{2}{c}{HR per 10-fold incr.} & \\multicolumn{1}{c}{P-value}    \\\\ 
       \\multicolumn{1}{l}{Immunologic Marker}  & \\multicolumn{1}{c}{No. at-risk**} & \\multicolumn{1}{c}{Pt. Est.} & \\multicolumn{1}{c}{95\\% CI} & \\multicolumn{1}{c}{(2-sided)} & \\multicolumn{1}{c}{Pt. Est.} & \\multicolumn{1}{c}{95\\% CI} & \\multicolumn{1}{c}{(2-sided)} & \\multicolumn{1}{c}{Pt. Est.} & \\multicolumn{1}{c}{95\\% CI} & \\multicolumn{1}{c}{(2-sided)}  \\\\ 
       \\hline\n 
  ")


# keep svycoxph in the name so that we can reuse rmd file
mytex(tab.1, file.name="CoR_univariable_svycoxph_pretty_"%.%fname.suffix, align="c", include.colnames = F, save2input.only=T, input.foldername=save.results.to,
    col.headers=header,
    longtable=T, 
    label=paste0("tab:CoR_univariable_svycoxph_pretty"), 
    caption.placement = "top", 
    caption=paste0("Inference for Day ", tpeak, "antibody marker covariate-adjusted correlates of risk of ", config.cor$txt.endpoint, " in the vaccine group: Hazard ratios per 10-fold increment in the marker*")
)
