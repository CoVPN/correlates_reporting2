# COR="D15to91covail_frnt"
# COR="D15to181covail_frnt"
# COR="D92to181covail_frnt"

renv::activate(project = here::here(".."))
Sys.setenv(TRIAL = "covail_frnt")
source(here::here("..", "_common.R")) 
source(here::here("code", "params.R"))

# hack
# source("~/copcor/R/cor_coxph_coef_1.R")
# source("~/copcor/R/cor_coxph_risk_tertile_incidence_curves.R")

{
Sys.setenv(VERBOSE = 1) 
quiet_library("kyotil") # p.adj.perm, getFormattedSummary
quiet_library("xtable") # this is a dependency of kyotil
quiet_library("marginalizedRisk")
quiet_library("tools") # toTitleCase
quiet_library("survey")
quiet_library("plotrix") # weighted.hist
quiet_library("parallel")
quiet_library("Hmisc") # wtd.quantile, cut2
quiet_library("mitools")
quiet_library("glue")

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

# all cases have covid lineage observed

if (COR=="D15to181covail_frnt") {
  form.0 = update(Surv(COVIDtimeD22toD181, COVIDIndD22toD181) ~ 1, as.formula(config$covariates))
  dat.sanofi =      subset(dat_proc, ph1.D15.frnt & TrtSanofi==1)
  dat.onedosemRNA = subset(dat_proc, ph1.D15.frnt & TrtonedosemRNA==1) 
  dat.onedosemRNA$yy=dat.onedosemRNA$COVIDIndD22toD181
  dat.sanofi$yy     =dat.sanofi$COVIDIndD22toD181
  
  trts=1:8 
  
} else if (COR=="D15to91covail_frnt") {
  form.0 = update(Surv(COVIDtimeD22toD91, COVIDIndD22toD91) ~ 1, as.formula(config$covariates))
  dat.sanofi =      subset(dat_proc, ph1.D15.frnt & TrtSanofi==1)
  dat.onedosemRNA = subset(dat_proc, ph1.D15.frnt & TrtonedosemRNA==1) 
  dat.onedosemRNA$yy=dat.onedosemRNA$COVIDIndD22toD91
  dat.sanofi$yy     =dat.sanofi$COVIDIndD22toD91
  
  trts=1:8 
  
} else if (COR=="D92to181covail_frnt") {
  form.0 = update(Surv(COVIDtimeD92toD181, COVIDIndD92toD181) ~ 1, as.formula(config$covariates))
  dat.sanofi =      subset(dat_proc, ph1.D92.frnt & TrtSanofi==1)
  dat.onedosemRNA = subset(dat_proc, ph1.D92.frnt & TrtonedosemRNA==1) 
  dat.onedosemRNA$yy=dat.onedosemRNA$COVIDIndD92toD181
  dat.sanofi$yy     =dat.sanofi$COVIDIndD92toD181

  trts=c(1,5:8) 
  
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
id50=c("pseudoneutid50_D614G", "pseudoneutid50_BA.1", "pseudoneutid50_BA.4.BA.5")

marker_sets = c("frnt", "id50")

}


for (trt in trts) {

  if (trt==1) {
    dat=subset(dat.onedosemRNA, naive==0)
    fname.suffix.0 = 'mRNA_onedose_NN'
  } else if (trt==2) {
    dat=subset(dat.onedosemRNA, TrtA==1 & naive==0)
    fname.suffix.0 = 'mRNA_Moderna_NN'
  } else if (trt==3) {
    dat=subset(dat.onedosemRNA, TrtA==0 & naive==0)
    fname.suffix.0 = 'mRNA_Pfizer_NN'
  } else if (trt==4) {
    dat=subset(dat.sanofi, naive==0)
    fname.suffix.0 = 'mRNA_Sanofi_NN'
    # no baseline id50 cat marker defined for this cohort
    
  } else if (trt==5) {
    dat=subset(dat.onedosemRNA, naive==1)
    fname.suffix.0 = 'mRNA_onedose_N'
  } else if (trt==6) {
    dat=subset(dat.onedosemRNA, TrtA==1 & naive==1)
    fname.suffix.0 = 'mRNA_Moderna_N'
  } else if (trt==7) {
    dat=subset(dat.onedosemRNA, TrtA==0 & naive==1)
    fname.suffix.0 = 'mRNA_Pfizer_N'
  } else if (trt==8) {
    dat=subset(dat.sanofi, naive==1)
    fname.suffix.0 = 'mRNA_Sanofi_N'
    # no baseline id50 cat marker defined for this cohort
  } 
    
  cat("\n",fname.suffix.0,"\n")
  
  
  # table of ph1 and ph2 cases
  tab1 = with(dat, table(ph2, EventIndPrimary))
  names(dimnames(tab1))[2] = "Event Indicator"; print(tab1)
  
  design <- twophase(id = list( ~ 1,  ~ 1), strata = list(NULL,  ~ Wstratum), subset =  ~ ph2, data = dat)
  
  for (marker_set in marker_sets) {
    cat("\n")
    myprint(marker_set)
    
    fname.suffix = fname.suffix.0%.%"_"%.%marker_set
    assays=get(marker_set)
    all.markers = c("B"%.%assays, "Day15"%.%assays, "Delta15overB"%.%assays)
    # all.markers.names.short = sub("Pseudovirus-", "", assay_metadata$assay_label_short[match(assays,assay_metadata$assay)])
    # all.markers.names.short = c("B "%.%all.markers.names.short, "D15 "%.%all.markers.names.short, "D15/B "%.%all.markers.names.short)
    # names(all.markers.names.short) = all.markers
    
    all.markers.names.short <- all.markers.names.long <- all.markers
    names(all.markers.names.short) = all.markers
    names(all.markers.names.long) = all.markers
    
    mytex(tab1, file.name = "tab1_" %.% fname.suffix, save2input.only = T, input.foldername = save.results.to)
    
    cat("\n")
    
    cor_coxph_coef_1(
      form.0,
      design_or_dat = design,
      fname.suffix,
      save.results.to,
      config,
      config.cor,

      markers = all.markers,
      markers.names.short = all.markers.names.short,
      
      run.trichtom=F,

      dat.plac = NULL,
      show.q=F, 
      verbose = T)


    # cor_coxph_risk_tertile_incidence_curves(
    #   form.0,
    #   dat,
    #   fname.suffix,
    #   save.results.to,
    #   config,
    #   config.cor,
    #   tfinal.tpeak,
    # 
    #   markers = all.markers,
    #   markers.names.short = all.markers.names.short,
    #   markers.names.long = all.markers.names.long,
    #   marker.cutpoints,
    #   assay_metadata,
    # 
    #   dat.plac = NULL,
    #   for.title="",
    #   verbose=T
    # )
  
  } # end loop marker_set

} # end loop trt




################################################################################
print(date())
print("cor_coxph run time: "%.%format(Sys.time()-time.start, digits=1))
