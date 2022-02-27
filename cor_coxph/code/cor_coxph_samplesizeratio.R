renv::activate(project = here::here(".."))     
library(kyotil)
library(Hmisc)
library(xtable) # this is a dependency of kyotil

Args <- commandArgs(trailingOnly=TRUE) 
if (length(Args)!=2) {
    stop("Two command line arguments are required: study_name, COR. For example, COVE D57 or ENSEMBLE D29IncludeNotMolecConfirmedstart1")
}
myprint(Args)
study.name <- Args[1]
COR=Args[2]

if (study.name=="COVE") {
    load(paste0("output/moderna_real/", COR, "/marginalized.risk.Rdata"))
} else if (study.name=="ENSEMBLE") {
    
} 
    
