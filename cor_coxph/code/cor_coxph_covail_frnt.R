# COR="D15to91covail_frnt"
# COR="D15to181covail_frnt"

renv::activate(project = here::here(".."))
Sys.setenv(TRIAL = "covail_frnt")
source(here::here("..", "_common.R")) 
source(here::here("code", "params.R"))


{
Sys.setenv(VERBOSE = 1) 
library(kyotil) # p.adj.perm, getFormattedSummary
library(xtable) # this is a dependency of kyotil
library(marginalizedRisk)
library(tools) # toTitleCase
library(survey)
library(plotrix) # weighted.hist
library(parallel)
library(Hmisc) # wtd.quantile, cut2
library(mitools)
library(glue)

time.start=Sys.time()
print(date())

multi.imp=FALSE

# path for figures and tables etc
save.results.to = here::here("output"); if (!dir.exists(save.results.to))  dir.create(save.results.to)
save.results.to = paste0(save.results.to, "/", attr(config,"config")); if (!dir.exists(save.results.to))  dir.create(save.results.to)
save.results.to = paste0(save.results.to, "/", COR, "/"); if (!dir.exists(save.results.to))  dir.create(save.results.to)


# B=1e3 and numPerm=1e4 take 10 min to run with 30 CPUS for one analysis
B <-       config$num_boot_replicates 
numPerm <- config$num_perm_replicates # number permutation replicates 1e4
myprint(B, numPerm)


# save cut points to files
marker.cutpoints = attr(dat_proc, "marker.cutpoints")
for (a in names(marker.cutpoints)) {        
  write(paste0(escape(a),     " [", concatList(round(marker.cutpoints[[a]], 2), ", "), ")%"), 
        file=paste0(save.results.to, "cutpoints_", a,".txt"))
}

# create centered version of markers for later use, which is necessary because we do not want to do scaling within naive and non-naive separately
for (a in c("Day15"%.%assays, "B"%.%assays, "Delta15overB"%.%assays)) {
  dat_proc[[a%.%"centered"]] = scale(dat_proc[[a]], scale=F)
}

dat.sanofi =      subset(dat_proc, ph1.D15.frnt & TrtSanofi==1)
dat.onedosemRNA = subset(dat_proc, ph1.D15.frnt & TrtonedosemRNA==1) 


# all cases have covid lineage observed

if (COR=="D15to181covail_frnt") {
  form.0 = update(Surv(COVIDtimeD22toD181, COVIDIndD22toD181) ~ 1, as.formula(config$covariates))
  dat.onedosemRNA$yy=dat.onedosemRNA$COVIDIndD22toD181
  dat.sanofi$yy     =dat.sanofi$COVIDIndD22toD181
  
} else if (COR=="D15to91covail_frnt") {
  form.0 = update(Surv(COVIDtimeD22toD91, COVIDIndD22toD91) ~ 1, as.formula(config$covariates))
  dat.onedosemRNA$yy=dat.onedosemRNA$COVIDIndD22toD91
  dat.sanofi$yy     =dat.sanofi$COVIDIndD22toD91
  
} else if (COR=="D92to181covail_frnt") {
  form.0 = update(Surv(COVIDtimeD92toD181, COVIDIndD92toD181) ~ 1, as.formula(config$covariates))
  dat.onedosemRNA$yy=dat.onedosemRNA$COVIDIndD92toD181
  dat.sanofi$yy     =dat.sanofi$COVIDIndD92toD181

} else stop("Wrong COR: "%.% COR)


form.1 = update(form.0, ~.-naive)

dat.n = subset(dat.onedosemRNA, naive==1)
dat.nn = subset(dat.onedosemRNA, naive==0)

prev.vacc = get.marginalized.risk.no.marker(form.0, dat.onedosemRNA, tfinal.tpeak)
myprint(prev.vacc)

# parameters for ph R scripts
show.q=F # no multiplicity adjustment
use.svy = F 
has.plac = F

frnt=assays[startsWith(assays, "frnt")]; 
assays = frnt

}




################################################################################
# Peak Obj 1-3 for mRNA vaccines

# 1, 11 and 12 are main effects models
# 2 and 21 are itxn models: marker * naive
# 3 and 31 are itxn models: D15 marker * baseline marker
# Note that 4 and 5 here are not peak object 4 and 5
# 4: like 1, but subset to naive
# 5: like 1, but subset to nnaive

for (iObj in c(4,5)) {
#for (iObj in c(1,11,12,2,21,3,31,4,5)) {
  # iObj=5; iPop=7
  
  # define the list of all.markers to work on
  # an item in the list need not be a single marker but is more like a formula

  if(iObj %in% c(1,4,5)) {
    all.markers = c("B"%.%assays, "Day15"%.%assays, "Delta15overB"%.%assays)
    all.markers.names.short = sub("Pseudovirus-", "", assay_metadata$assay_label_short[match(assays,assay_metadata$assay)])
    all.markers.names.short = c("B "%.%all.markers.names.short, "D15 "%.%all.markers.names.short, "D15/B "%.%all.markers.names.short)
    
  } else if(iObj==11){
    # B marker + D15/B
    all.markers = sapply(assays, function (a) paste0("B",a, "centered + Delta15overB",a, "centered"))
    all.markers.names.short = sub("Pseudovirus-", "", assay_metadata$assay_label_short[match(assays,assay_metadata$assay)])
    # parameters for R script
    nCoef=2
    col.headers=c("center(B)", "center(D15/B)")
    
  } else if(iObj==12){
    # B marker + D15 + D15^2
    all.markers = sapply(assays, function (a) paste0("B",a, "centered + Day15",a, "centered + I(Day15",a, "centered^2)"))
    all.markers.names.short = sub("Pseudovirus-", "", assay_metadata$assay_label_short[match(assays,assay_metadata$assay)])
    # parameters for R script
    nCoef=3
    col.headers=c("center(B)", "center(D15)", "center(D15)**2")
    
  } else if(iObj==2){
    # interaction naive x D15 marker
    all.markers = c(sapply(assays, function (a) paste0("naive * Day15",a, "centered")),
                    sapply(assays, function (a) paste0("naive * Delta15overB",a, "centered")))
    all.markers.names.short = sub("Pseudovirus-", "", assay_metadata$assay_label_short[match(assays,assay_metadata$assay)])
    all.markers.names.short = c("D15 "%.%all.markers.names.short, "D15/B "%.%all.markers.names.short)
    # parameters for R script
    nCoef=3
    col.headers=c("naive", "center(D15 or fold)", "naive:center(D15 or fold)")
    
  } else if(iObj==21){
    # interaction (1-naive) x D15 marker
    all.markers = c(sapply(assays, function (a) paste0("I(1-naive) * Day15",a, "centered")),
                    sapply(assays, function (a) paste0("I(1-naive) * Delta15overB",a, "centered")))
    all.markers.names.short = sub("Pseudovirus-", "", assay_metadata$assay_label_short[match(assays,assay_metadata$assay)])
    all.markers.names.short = c("D15 "%.%all.markers.names.short, "D15/B "%.%all.markers.names.short)
    nCoef=3
    col.headers=c("I(1-naive)", "center(D15 or fold)", "I(1-naive):center(D15 or fold)")

  } else if(iObj==3){
    # interaction B marker x D15 marker or D15/B
    all.markers = c(sapply(assays, function (a) paste0("B",a, "centered * Day15",a, "centered")),
                    sapply(assays, function (a) paste0("B",a, "centered * Delta15overB",a, "centered"))
    )
    all.markers.names.short = sub("Pseudovirus-", "", assay_metadata$assay_label_short[match(assays,assay_metadata$assay)])
    all.markers.names.short = c("D15 "%.%all.markers.names.short, "D15/B "%.%all.markers.names.short)
    # parameters for R script
    nCoef=3
    col.headers=c("center(B)", "center(D15 or fold)", "center(B):center(D15 or fold)")
    
  } else if(iObj==31){
    # B_cat marker x D15 marker or D15/B
    all.markers = c(sapply(assays, function (a) paste0("B",a, "cat * Day15",a, "centered")),
                    sapply(assays, function (a) paste0("B",a, "cat * Delta15overB",a, "centered")))
    all.markers.names.short = sub("Pseudovirus-", "", assay_metadata$assay_label_short[match(assays,assay_metadata$assay)])
    all.markers.names.short = c("D15 "%.%all.markers.names.short, "D15/B "%.%all.markers.names.short)
    # table col names
    col1="Med B:center(D15 or fold)"
    col2="High B:center(D15 or fold)"
    col3="center(D15 or fold) at Low B"
    col4="center(D15 or fold) at Med B"
    col5="center(D15 or fold) at High B"
    
  }
  names(all.markers.names.short) = all.markers
  
  # repeat all objectives over several subpopulations, save results with different fname.suffix
  for (iPop in 1:4) {
    cat("\n")
    myprint(iObj, iPop)
    
    if (iPop==1) {
      dat=dat.onedosemRNA
      fname.suffix = 'mRNA_onedose'
    } else if (iPop==2) {
      dat=subset(dat.onedosemRNA, TrtA==1)
      fname.suffix = 'mRNA_Moderna'
    } else if (iPop==3) {
      dat=subset(dat.onedosemRNA, TrtA==0)
      fname.suffix = 'mRNA_Pfizer'
    } else if (iPop==4) {
      dat=dat.sanofi
      fname.suffix = 'mRNA_Pfizer'
      
    # } else if (iPop==4) {
    #   dat=subset(dat.onedosemRNA, TrtB==1)
    #   fname.suffix = 'mRNA_Prototype'
    # } else if (iPop==5) {
    #   dat=subset(dat.onedosemRNA, TrtB==0)
    #   fname.suffix = 'mRNA_Omicron-Containing'
    #   
    # } else if (iPop==6) {
    #   dat=subset(dat.onedosemRNA, TrtC==1)
    #   fname.suffix = 'mRNA_Bivalent'
    # } else if (iPop==7) {
    #   dat=subset(dat.onedosemRNA, TrtC==0)
    #   fname.suffix = 'mRNA_Monovalent'
    } 
    
    if(iObj==4) dat=subset(dat, naive==1)
    if(iObj==5) dat=subset(dat, naive==0)
    
    # table of ph1 and ph2 cases
    tab1 = with(dat, table(ph2, EventIndPrimary))
    names(dimnames(tab1))[2] = "Event Indicator"; print(tab1)
    mytex(tab1, file.name = "tab1_" %.% fname.suffix%.%ifelse(iObj==4, "_N", "_NN"), save2input.only = T, input.foldername = save.results.to)
    
    # if(iObj==4 | iObj==5) {
    #   cor_coxph_coef_1(
    #     form.0,
    #     design_or_dat = dat,
    #     fname.suffix=fname.suffix%.%ifelse(iObj==4, "_N", "_NN"),
    #     save.results.to,
    #     config,
    #     config.cor,
    #     
    #     all.markers,
    #     all.markers.names.short,
    #     
    #     dat.plac = NULL,
    #     show.q=F, # whether to show fwer and q values in tables
    #     run.trichtom=T,
    #     verbose = T)
    #   
    #   # trichotomitized curves
    #   cor_coxph_risk_tertile_incidence_curves(
    #     # need to remove naive from formula. otherwise risk will be NA
    #     form.0 = as.formula(
    #       sub("naive", 1, # if naive is the first covariate, there will be no +, so the next line won't work
    #         sub("\\+ naive", "", paste0(deparse(form.0,width.cutoff=500)))
    #       )
    #     ),
    #     dat,
    #     fname.suffix=fname.suffix%.%ifelse(iObj==4, "_N", "_NN"),
    #     save.results.to,
    #     config,
    #     config.cor,
    #     tfinal.tpeak,
    #     
    #     markers = all.markers,
    #     markers.names.short = all.markers.names.short,
    #     markers.names.long = all.markers.names.long,
    #     marker.cutpoints,
    #     assay_metadata,
    #     
    #     dat.plac = NULL,
    #     for.title="", 
    #     verbose=T
    #   )
    # } 
  }
  
  
  
  
}


################################################################################
# Peak Obj 4

# same formula, repeated 3 times
# for (iPop in 1:3) {
#   if (iPop==1) {
#     dat=subset(dat.onedosemRNA, TrtA %in% c(1,0))
#     fname.suffix = 'mRNA_Mod_Pfi'
#     
#     all.markers = sapply(assays, function (a) paste0("TrtA * Day15",a, "centered"))
#     all.markers.names.short = sub("Pseudovirus-", "", assay_metadata$assay_label_short[match(assays,assay_metadata$assay)])
#     # parameters for R script
#     nCoef=3
#     col.headers=c("TrtA Moderna~Pfizer", "center(D15)", "TrtA:center(D15)")
#     
#   } else if (iPop==2) {
#     dat=subset(dat.onedosemRNA, TrtB %in% c(1,0))
#     fname.suffix = 'mRNA_Pro_Omi'
#     
#     all.markers = sapply(assays, function (a) paste0("TrtB * Day15",a, "centered"))
#     all.markers.names.short = sub("Pseudovirus-", "", assay_metadata$assay_label_short[match(assays,assay_metadata$assay)])
#     # parameters for R script
#     nCoef=3
#     col.headers=c("TrtB Prot~Omicron", "center(D15)", "TrtB:center(D15)")
#     
#   } else if (iPop==3) {
#     dat=subset(dat.onedosemRNA, TrtC %in% c(1,0))
#     fname.suffix = 'mRNA_Bi_Mono'
#     
#     all.markers = sapply(assays, function (a) paste0("TrtC * Day15",a, "centered"))
#     all.markers.names.short = sub("Pseudovirus-", "", assay_metadata$assay_label_short[match(assays,assay_metadata$assay)])
#     # parameters for R script
#     nCoef=3
#     col.headers=c("TrtC Biv~Monovalent", "center(D15)", "TrtC:center(D15)")
#   } 
#   
#   cor_coxph_coef_n (
#     form.0,
#     design_or_dat = dat,
#     fname.suffix,
#     save.results.to,
#     config,
#     config.cor,
#     all.markers,
#     all.markers.names.short,
#     
#     nCoef,
#     col.headers,
#     verbose = T)
# }
  




################################################################################
print(date())
print("cor_coxph run time: "%.%format(Sys.time()-time.start, digits=1))
