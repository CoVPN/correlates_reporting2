
 
load.data=function(TRIAL, COR, trt=1) {
    Sys.setenv("TRIAL"=TRIAL)
    source(here::here("..", "_common.R"), local=T)
    # uloq censoring    
    for (a in assays) {
        tmp="Day"%.%config.cor$tpeak %.% a
        dat_proc[[tmp]] <- ifelse(dat_proc[[tmp]] > log10(uloqs[a]), log10(uloqs[a]), dat_proc[[tmp]])
    }
    if (trt==1) {
        dat=subset(dat_proc, Trt==1 & ph1)
        dat = add.trichotomized.markers (dat, tpeak, wt.col.name="wt")
    } else {
        dat=subset(dat_proc, Trt==0 & ph1)        
    }
    dat
}


get.trial=function(x, assay) {
    if (startsWith(x,"janssen")) {
        if(startsWith(assay,"bind")) x=paste0(x, "")
        if(startsWith(assay,"pseudo")) x=paste0(x, "")
        if(startsWith(assay,"ADCP")) x=paste0(x, "")
    } else if (startsWith(x,"azd1222")) {
        if(startsWith(assay,"bind")) x=paste0(x, "_bAb")
        if(startsWith(assay,"pseudo")) x=paste0(x, "")
    } 
    x
}






# need this function b/c svycoxh may error due to singularity if, e.g. all cases have the same marker value
run.svycoxph=function(f, design) {
    fit=try(svycoxph(f, design=design), silent=T)
    if (class(fit)[1]=="try-error") NA else fit
}





# forest.covail=function (fits, names, fname.suffix, save.results.to) {
#   
#   names(fits)=names
#   est.ci = sapply(fits, function (fit) {
#     if (length(fit)==1) return (rep(NA,4))
#     tmp=getFixedEf(fit, exp=T, robust=F)
#     tmp[nrow(tmp),c("HR", "(lower", "upper)", "p.value")]
#   })
#   
#   #nevents=rep(sum(dat$yy==1), length(fits))
#   # not showing nevents 
#   nevents=rep(NA, length(fits))
#   
#   mypdf(onefile=F, width=10,height=4, file=paste0(save.results.to, "hr_forest_", fname.suffix)) 
#   theforestplot(point.estimates=est.ci[1,], lower.bounds=est.ci[2,], upper.bounds=est.ci[3,], group=colnames(est.ci), 
#                 nEvents=nevents, title=paste0(""), p.values=NA, 
#                 decimal.places=2, graphwidth=unit(120, "mm"), fontsize=1.2, 
#                 table.labels = c("", "  HR (95% CI)",""), 
#                 x.ticks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4) # controls the limit
#   )
#   dev.off()
#   
# 
#   mypdf(onefile=F, width=10,height=4, file=paste0(save.results.to, "hr_forest_log_", fname.suffix)) 
#   theforestplot(point.estimates=est.ci[1,], lower.bounds=est.ci[2,], upper.bounds=est.ci[3,], group=colnames(est.ci), 
#                 nEvents=nevents, title=paste0(""), p.values=NA, 
#                 decimal.places=2, graphwidth=unit(120, "mm"), fontsize=1.2, 
#                 table.labels = c("", "  HR (95% CI)",""), 
#                 xlog=T,
#                 x.ticks = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4) # controls the limit
#   )
#   dev.off()
#   
# }




# for some subpopulations, some strata may have empty ph2 
# this function removes those strata
get.dat.with.no.empty=function(dat.tmp) {
    tab=with(dat.tmp, table(Wstratum, ph2))
    subset(dat.tmp, !Wstratum %in% as.integer(rownames(tab)[which(tab[,"TRUE"]==0)]))
}



# a more comprehensive, slower version of copcor::escape
escape_latex <- function(text) {
  # Define the special characters and their LaTeX escaped equivalents
  special_chars <- c("\\", "%", "$", "#", "_", "{", "}", "&", "^", "~", "<", ">", "|", "\"")
  latex_escapes <- c("\\textbackslash{}", "\\%", "\\$", "\\#", "\\_", "\\{", "\\}", "\\&", "\\textasciicircum{}", "\\textasciitilde{}", "\\textless{}", "\\textgreater{}", "\\textbar{}", "\\textquotedbl{}")
  
  # Replace each special character in the text with its escaped version
  for (i in seq_along(special_chars)) {
    text <- gsub(special_chars[i], latex_escapes[i], text, fixed = TRUE)
  }
  
  return(text)
}
