renv::activate(project = here::here(".."))    
    # There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
    if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))
    
library(kyotil)
library(tools) # toTitleCase
library(xtable) # this is a dependency of kyotil


trials=c("janssen_pooled_real", "janssen_pooled_realPsV", "janssen_pooled_realADCP")

for (COR in c("D29IncludeNotMolecConfirmed", "D29IncludeNotMolecConfirmedstart1")) {
    out=lapply(trials, function (trial) {
        load (file=paste0(here::here("output", trial, COR), "/coxph_slopes.Rdata")) 
        list(tab.cont, tab.cat, save.s.1, save.s.2)
    })
        
    tab.cont=do.call(rbind, sapply(out, function (res) res[[1]]))
    tab.cat=do.call(rbind, sapply(out, function (res) res[[2]]))
    
    mytex(tab.cont, file.name="CoR_univariable_svycoxph_pretty_ENSEMBLE_"%.%COR, align="c", include.colnames = F, save2input.only=T, input.foldername=here::here("output"),
        col.headers=paste0("\\hline\n 
             \\multicolumn{1}{l}{} & \\multicolumn{1}{c}{No. cases /}   & \\multicolumn{2}{c}{HR per 10-fold incr.}                     & \\multicolumn{1}{c}{P-value}   & \\multicolumn{1}{c}{q-value}   & \\multicolumn{1}{c}{FWER} \\\\ 
             \\multicolumn{1}{l}{Immunologic Marker}            & \\multicolumn{1}{c}{No. at-risk**} & \\multicolumn{1}{c}{Pt. Est.} & \\multicolumn{1}{c}{95\\% CI} & \\multicolumn{1}{c}{(2-sided)} & \\multicolumn{1}{c}{***} & \\multicolumn{1}{c}{} \\\\ 
             \\hline\n 
        ")
    )    
    
    mytex(tab.cat, file.name="CoR_univariable_svycoxph_cat_pretty_ENSEMBLE_"%.%COR, align="c", include.colnames = F, save2input.only=T, input.foldername=here::here("output"),
        col.headers=paste0("\\hline\n 
             \\multicolumn{1}{l}{} & \\multicolumn{1}{c}{Tertile}   & \\multicolumn{1}{c}{No. cases /}   & \\multicolumn{1}{c}{Attack}   & \\multicolumn{2}{c}{Haz. Ratio}                     & \\multicolumn{1}{c}{P-value}   & \\multicolumn{1}{c}{Overall P-}      & \\multicolumn{1}{c}{Overall q-}   & \\multicolumn{1}{c}{Overall} \\\\ 
             \\multicolumn{1}{l}{Immunologic Marker}            & \\multicolumn{1}{c}{}          & \\multicolumn{1}{c}{No. at-risk**} & \\multicolumn{1}{c}{rate}   & \\multicolumn{1}{c}{Pt. Est.} & \\multicolumn{1}{c}{95\\% CI} & \\multicolumn{1}{c}{(2-sided)} & \\multicolumn{1}{c}{value***} & \\multicolumn{1}{c}{value $\\dagger$} & \\multicolumn{1}{c}{FWER} \\\\ 
             \\hline\n 
        "),        
        add.to.row=list(list(nrow(tab.cat)), # insert at the beginning of table, and at the end of, say, the first table
            c(paste0(" \n \\multicolumn{8}{l}{} \\\\ \n", 
                      "\n \\multicolumn{2}{l}{Placebo} & ", 
                     out[[1]][[3]], "&",  
                     out[[1]][[4]], "&",  
                     "\\multicolumn{4}{l}{}  \\\\ \n")
             )
        )
    )    
    
}
