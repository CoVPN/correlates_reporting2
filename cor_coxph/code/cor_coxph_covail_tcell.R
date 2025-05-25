# COR="D15to91covail_tcell"
# COR="D15to181covail_tcell"

# renv::activate(project = here::here(".."))
Sys.setenv(TRIAL = "covail_tcell"); source(here::here("..", "_common.R")); source(here::here("code", "params.R"))

# hack
# source("~/copcor/R/cor_coxph_coef_1.R")

marker_sets = c("primary", "secondary", "exploratory")
trts=1:8 # onedosemRNA, etc
# trt=8; marker_set = 2

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
myprint(save.results.to)

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

# # experiment with truncating t cell markers at 0.02 or -1.699
# for (a in c(primary,secondary)) {
#   if (endsWith(a, ".N")) {
#     dat_proc[[a]] = ifelse(dat_proc[[a]] < -2, -2, dat_proc[[a]])
#   } else {
#     dat_proc[[a]] = ifelse(dat_proc[[a]] < -1.699, -1.699, dat_proc[[a]])
#   }
# }
# summary(dat_proc[primary])
# # filter out ones with no variablity
# secondary=setdiff(secondary, c("Bcd4_IL4.IL5.IL13.154_Wuhan.N"))

# note that arm 16 and 17 are excluded, because no T cells are done for them. This is different from the antibody correlates where 16 and 17 are included.
dat.onedosemRNA  =    subset(dat_proc, ph1.D15.tcell & TrtonedosemRNA==1 & !arm %in% c(16,17)) 
dat.onedoseModerna = subset(dat_proc, ph1.D15.tcell & arm %in% c(1,2,5,6))
dat.onedosePfizer  = subset(dat_proc, ph1.D15.tcell & arm %in% c(7,9,12))
dat.sanofi = subset(dat_proc, ph1.D15.tcell & TrtSanofi==1)

# all cases have covid lineage observed

if (startsWith(COR,"D15to181")) {
  form.0 = update(Surv(COVIDtimeD22toD181, COVIDIndD22toD181) ~ 1, as.formula(config$covariates))

} else if (startsWith(COR,"D15to91")) {
  form.0 = update(Surv(COVIDtimeD22toD91, COVIDIndD22toD91) ~ 1, as.formula(config$covariates))

} else if (startsWith(COR,"D92to181")) {
  form.0 = update(Surv(COVIDtimeD92toD181, COVIDIndD92toD181) ~ 1, as.formula(config$covariates))

} else stop("Wrong COR: "%.% COR)

# updated config to remove naive
# form.0 = update(form.0, ~.-naive)

prev.vacc = get.marginalized.risk.no.marker(form.0, dat.onedosemRNA, tfinal.tpeak)
myprint(prev.vacc)

nAbs=subset(assay_metadata, panel=="id50", assay, drop=T)

# parameters for ph R scripts
show.q=F # no multiplicity adjustment
use.svy = F 
has.plac = F
}


################################################################################
# time periods (vaccine proximal, distal and all) are implemented through COR


for (trt in trts) {
  # naive
  if (trt==1) {
    dat.1=subset(dat.onedosemRNA, naive==1);    fname.suffix.0 <- trt.label <- "OneDosemRNA_Naive"
  } else if (trt==2) {
    dat.1=subset(dat.onedoseModerna, naive==1); fname.suffix.0 <- trt.label <- "OneDoseModerna_Naive"
  } else if (trt==3) {
    dat.1=subset(dat.onedosePfizer, naive==1);  fname.suffix.0 <- trt.label <- "OneDosePfizer_Naive"
  } else if (trt==4) {
    dat.1=subset(dat.sanofi, naive==1);         fname.suffix.0 <- trt.label <- "Sanofi_Naive"
    
  # nnaive  
  } else if (trt==5) {
    dat.1=subset(dat.onedosemRNA, naive==0);    fname.suffix.0 <- trt.label <- "OneDosemRNA_NNaive"
  } else if (trt==6) {
    dat.1=subset(dat.onedoseModerna, naive==0); fname.suffix.0 <- trt.label <- "OneDoseModerna_NNaive"
  } else if (trt==7) {
    dat.1=subset(dat.onedosePfizer, naive==0);  fname.suffix.0 <- trt.label <- "OneDosePfizer_NNaive"
  } else if (trt==8) {
    dat.1=subset(dat.sanofi, naive==0);         fname.suffix.0 <- trt.label <- "Sanofi_NNaive"
  }
  
  design.1 <- twophase(id = list( ~ 1,  ~ 1), strata = list(NULL,  ~ Wstratum), subset =  ~ ph2, data = dat.1)
  dat.0=NULL
  
  # table of ph1 and ph2 cases
  tab1 = with(dat.1, table(ph2, EventIndPrimary))
  names(dimnames(tab1))[2] = "Event Indicator"; print(tab1)
  
  for (marker_set in marker_sets) {
    
    fname.suffix = fname.suffix.0%.%"_"%.%marker_set
    all.markers=get(marker_set)
    all.markers.names.short <- all.markers.names.long <- all.markers
    names(all.markers.names.short) = all.markers
    names(all.markers.names.long) = all.markers

    # need to save tab1 for each distinct fname.suffix
    mytex(tab1, file.name = "tab1_" %.% fname.suffix, save2input.only = T, input.foldername = save.results.to)
    
    # # create centered version of markers for later use, which is necessary because we do not want to center within naive and non-naive separately
    # for (a in c("Day15"%.%all.markers, "B"%.%all.markers, "Delta15overB"%.%all.markers)) {
    #   dat_proc[[a%.%"centered"]] = scale(dat_proc[[a]], scale=F)
    # }
    
    cat("\n\n")
    cor_coxph_coef_1 (
      form.0 = if(trt==1 | trt==5) update(form.0, ~.+ strata(stage)) else form.0, 
      design_or_dat = design.1,
      fname.suffix,
      save.results.to,
      config,
      config.cor,
      markers = all.markers,
      markers.names.short = all.markers.names.short,
      
      dat.plac = dat.0,
      show.q = marker_set=="primary",
      
      forestplot.markers=NULL, 
      for.title="",
      run.trichtom=TRUE,
      verbose = T
    )
    
    cor_coxph_risk_tertile_incidence_curves (
      form.0 = if(trt==1 | trt==5) update(form.0, ~.+ strata(stage)) else form.0, 
      dat = dat.1,
      fname.suffix,
      save.results.to,
      config,
      config.cor,
      tfinal.tpeak,
      
      markers = all.markers,
      markers.names.short = all.markers.names.short,
      markers.names.long = all.markers.names.long,
      marker.cutpoints,
      assay_metadata,
      
      dat.plac = dat.0,
      for.title = "",
      trt.label = trt.label, 
      verbose=T
    )
    
    
  }
  
}



################################################################################
print(date())
print("cor_coxph run time: "%.%format(Sys.time()-time.start, digits=1))
