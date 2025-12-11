#COR="D15toEOS_cov2008_tcell_symp"
#COR="D15toEOS_cov2008_tcell_anyinf"

renv::activate(project = here::here(".."))
Sys.setenv(TRIAL = "cov2008_tcell")
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

dat_proc$yy=dat_proc$EventIndPrimary

prev.vacc = get.marginalized.risk.no.marker(form.0, dat_proc[dat_proc$ph1==1,], tfinal.tpeak)
myprint(prev.vacc)

}




################################################################################

for (trt in 1) {
  
  dat.1=subset(dat_proc, ph1.D15.tcell==1);    fname.suffix.0 <- trt.label <- "D15"

  design.1 <- twophase(id = list( ~ 1,  ~ 1), strata = list(NULL,  ~ Wstratum.tcell), subset =  ~ ph2.D15.tcell, data = dat.1)
  
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
