library(methods)
library(dplyr)
library(kyotil)
library(copcor)
?library(marginalizedRisk)
library(survival)

# disable lower level parallelization in favor of higher level of parallelization
library(RhpcBLASctl)
blas_get_num_procs()
blas_set_num_threads(1L)
#stopifnot(blas_get_num_procs() == 1L) # Commented this out as it does not work as expected any more!
omp_set_num_threads(1L)
stopifnot(omp_get_max_threads() == 1L)
    
set.seed(98109)
    
if(!exists("verbose")) verbose=0
if (Sys.getenv("VERBOSE") %in% c("T","TRUE")) verbose=1
if (Sys.getenv("VERBOSE") %in% c("1", "2", "3")) verbose=as.integer(Sys.getenv("VERBOSE"))
    
# COR defines the analysis to be done, e.g. D14
if(!exists("COR")) {
    if(!exists("Args")) Args <- commandArgs(trailingOnly=TRUE)
    if (length(Args)>0) {
        COR=Args[1]
        myprint(COR)
    } else {
        warning("No COR. This is okay if _common.R is sourced just to load common functions. If needed, COR can be defined through command line argument or in R script before _common.R is sourced.")
    }
}

# if DESCRIPTIVE env variable is set, then we are doing descriptive analyses, e.g. immuno_tabular, or cor_graphical
DESCRIPTIVE = Sys.getenv("DESCRIPTIVE") %in% c("1", "T", "TRUE")
myprint(DESCRIPTIVE)

# if EXPOSUREPROXIMAL env variable is set, then we are doing exposure-proximal analyses
EXPOSUREPROXIMAL = Sys.getenv("EXPOSUREPROXIMAL") %in% c("1", "T", "TRUE")
myprint(EXPOSUREPROXIMAL)


if(Sys.getenv("TRIAL")=="") {
  stop(" *************************************  environmental variable TRIAL not defined  *************************************")
}

config <- config::get(config = Sys.getenv("TRIAL"))
for(opt in names(config)) eval(parse(text = paste0(names(config[opt])," <- config[[opt]]")))
TRIAL=attr(config, "config")


DayPrefix = switch (TRIAL, 
                    'moderna_boost' = "BD", 
                    'id27hpv' = "M",
                    'id27hpvnAb' = "M",
                    "Day")

if (is.null(config$threshold_grid_size)) {
  # Should be 15 at least for the plots of the threshold-response and its inverse to be representative of the true functions.
  threshold_grid_size <- 30 
} else {
  threshold_grid_size = config$threshold_grid_size
}


STOP("need to filter out bindSpike_D614 form assays")
###################################################################################################
# assay metadata

if (!is.null(config$assay_metadata)) {
  
  # created named lists for assay metadata to easier access, e.g. assay_labels_short["bindSpike"]
  assay_metadata = read.csv(paste0(dirname(attr(config,"file")),"/",config$assay_metadata))
  
  if(any(is.na(assay_metadata$uloq))) stop('uloq cannot be NA, set it to Inf if not needed')
  
  # remove bindN
  assay_metadata=subset(assay_metadata, assay!="bindN")
  
  if (TRIAL=='vat08_combined') {
    if (exists('COR')) {
      if (contain(COR, "nAb") | endsWith(COR,'original2')) {
        # only keeps ID50 markers
        assay_metadata = subset(assay_metadata, panel=='id50')
      } else {
        # assay_metadata = subset(assay_metadata, panel=='bindSpike')
        if (!DESCRIPTIVE) {
          # change lloq for bAb to min(...) 
          lloq_min = min (subset(assay_metadata, panel=='bindSpike' & assay!="bindSpike_mdw", lloq))
          assay_metadata[assay_metadata$panel=='bindSpike' & assay_metadata$assay!="bindSpike_mdw",'lloq'] = lloq_min
        }
      }
    }
    
  } else if (TRIAL=='id27hpv') {
    if (exists('COR')) {
      assay_metadata = subset(assay_metadata, panel=='bind')
    }
    
  } else if (TRIAL=='id27hpvnAb') {
    if (exists('COR')) {
      assay_metadata = subset(assay_metadata, panel=='id50')
    }
    
  }
  

  assays=assay_metadata$assay
  
  labels.assays=assay_metadata$assay_label; names(labels.assays)=assays
  labels.assays.short=assay_metadata$assay_label_short; names(labels.assays.short)=assays

  llox_labels=assay_metadata$llox_label; names(llox_labels)=assays
  lloqs=assay_metadata$lloq; names(lloqs)=assays
  uloqs=assay_metadata$uloq; names(uloqs)=assays
  lods=assay_metadata$lod; names(lods)=assays
  lloxs=ifelse(llox_labels=="lloq", lloqs, lods)
  lloxs=ifelse(llox_labels=="pos", assay_metadata$pos.cutoff, lloxs)
  
  
} else {
  
  lloxs=NULL
  
  if(length(config$llox_label)==1) {
    llox_labels = rep(config$llox_label, length(config$assays))
  } else {
    stopifnot(length(config$llox_label)==length(config$assays))
    llox_labels = config$llox_label
  }
  names(llox_labels)=config$assays

  # assays labels. This needs to come before all.markers
  labels.assays = config$assay_labels
  names(labels.assays) = config$assays
  
  if (is.null(config$assay_labels_short)) {
    labels.assays.short=labels.assays
  } else {
    labels.assays.short = config$assay_labels_short
    names(labels.assays.short) = config$assays
  }
  
  names(assays)=assays # add names so that lapply results will have names
  
  # if the following part changes, make sure to copy to _common.R in the processing repo
  
  # uloqs etc are hardcoded for ows trials but driven by config for other trials
  # For bAb, IU and BAU are the same thing
  # all values on BAU or IU
  # LOQ can not be NA, it is needed for computing delta
  pos.cutoffs<-llods<-lloqs<-uloqs<-c()
  if (study_name %in% c("COVE", "MockCOVE", "MockENSEMBLE")) {
    tmp=list(
      bindSpike=c(
        pos.cutoff=10.8424,
        LLOD = 0.3076,
        ULOD = 172226.2,
        LLOQ = 1.7968,
        ULOQ = 10155.95)
      ,
      bindRBD=c(
        pos.cutoff=14.0858,
        LLOD = 1.593648,
        ULOD = 223074,
        LLOQ = 3.4263,
        ULOQ = 16269.23)
      ,
      bindN=c( 
        pos.cutoff=23.4711,
        LLOD = 0.093744,
        ULOD = 52488,
        LLOQ = 4.4897,
        ULOQ = 574.6783)
      ,
      pseudoneutid50=c( 
        pos.cutoff=2.42,# as same lod
        LLOD = 2.42,
        ULOD = NA,
        LLOQ = 4.477,
        ULOQ = 10919)
      ,
      pseudoneutid80=c( 
        pos.cutoff=15.02,# as same lod
        LLOD = 15.02,
        ULOD = NA,
        LLOQ = 21.4786,
        ULOQ = 15368)
      ,
      liveneutmn50=c( 
        pos.cutoff=82.1*0.276,# as same lod
        LLOD = 82.11*0.276,
        ULOD = NA,
        LLOQ =  159.79*0.276,
        ULOQ = 11173.21*0.276)
    )
    
    pos.cutoffs=sapply(tmp, function(x) unname(x["pos.cutoff"]))
    llods=sapply(tmp, function(x) unname(x["LLOD"]))
    lloqs=sapply(tmp, function(x) unname(x["LLOQ"]))
    uloqs=sapply(tmp, function(x) unname(x["ULOQ"]))        
    
    
  } else if(study_name=="ENSEMBLE") {
    
    # data less than pos cutoff is set to pos.cutoff/2
    llods["bindSpike"]=NA 
    lloqs["bindSpike"]=1.7968 
    uloqs["bindSpike"]=238.1165 
    pos.cutoffs["bindSpike"]=10.8424
    
    # data less than pos cutoff is set to pos.cutoff/2
    llods["bindRBD"]=NA                 
    lloqs["bindRBD"]=3.4263                 
    uloqs["bindRBD"]=172.5755    
    pos.cutoffs["bindRBD"]=14.0858
    
    # data less than lod is set to lod/2
    llods["ADCP"]=11.57
    lloqs["ADCP"]=8.87
    uloqs["ADCP"]=211.56
    pos.cutoffs["ADCP"]=11.57# as same lod
    
    llods["bindN"]=0.093744
    lloqs["bindN"]=4.4897
    uloqs["bindN"]=574.6783
    pos.cutoffs["bindN"]=23.4711
    
    # the limits below are different for EUA and Part A datasets
    if (contain(TRIAL, "EUA")) {
      # EUA data
      
      # data less than lloq is set to lloq/2
      llods["pseudoneutid50"]=NA  
      lloqs["pseudoneutid50"]=42*0.0653  #2.7426
      uloqs["pseudoneutid50"]=9484*0.0653 # 619.3052
      pos.cutoffs["pseudoneutid50"]=lloqs["pseudoneutid50"]
      
      # repeat for two synthetic markers that are adapted to SA and LA
      llods["pseudoneutid50sa"]=NA  
      lloqs["pseudoneutid50sa"]=42*0.0653  #2.7426
      uloqs["pseudoneutid50sa"]=9484*0.0653 # 619.3052
      pos.cutoffs["pseudoneutid50sa"]=lloqs["pseudoneutid50sa"]
      
      llods["pseudoneutid50la"]=NA  
      lloqs["pseudoneutid50la"]=42*0.0653  #2.7426
      uloqs["pseudoneutid50la"]=9484*0.0653 # 619.3052
      pos.cutoffs["pseudoneutid50la"]=lloqs["pseudoneutid50la"]
      
      
    } else if (contain(TRIAL, "partA")) {
      # complete part A data
      
      # data less than lloq is set to lloq/2
      llods["pseudoneutid50"]=NA  
      lloqs["pseudoneutid50"]=75*0.0653  #4.8975
      uloqs["pseudoneutid50"]=12936*0.0653 # 844.7208
      pos.cutoffs["pseudoneutid50"]=lloqs["pseudoneutid50"]
      
      # repeat for two synthetic markers that are adapted to SA and LA
      llods["pseudoneutid50sa"]=NA  
      lloqs["pseudoneutid50sa"]=75*0.0653  #4.8975
      uloqs["pseudoneutid50sa"]=12936*0.0653 # 844.7208
      pos.cutoffs["pseudoneutid50sa"]=lloqs["pseudoneutid50sa"]
      
      llods["pseudoneutid50la"]=NA  
      lloqs["pseudoneutid50la"]=75*0.0653  #4.8975
      uloqs["pseudoneutid50la"]=12936*0.0653 # 844.7208
      pos.cutoffs["pseudoneutid50la"]=lloqs["pseudoneutid50la"]
    }
    
    # data less than lod is set to lod/2
    llods["pseudoneutid50uncensored"]=40*0.0653 #2.612
    lloqs["pseudoneutid50uncensored"]=40*0.0653  
    uloqs["pseudoneutid50uncensored"]=12936*0.0653 # 844.7208
    pos.cutoffs["pseudoneutid50uncensored"]=lloqs["pseudoneutid50uncensored"]
    
  } else if(study_name=="PREVENT19") {
    # Novavax
    
    # data less than lloq is set to lloq/2 in the raw data
    llods["bindSpike"]=NA 
    lloqs["bindSpike"]=150.4*0.0090 # 1.3536
    uloqs["bindSpike"]=770464.6*0.0090 # 6934.181
    pos.cutoffs["bindSpike"]=10.8424 # use same as COVE
    
    # data less than lloq is set to lloq/2
    llods["bindRBD"]=NA  
    lloqs["bindRBD"]=1126.7*0.0272  #30.6
    uloqs["bindRBD"]=360348.7*0.0272 # 9801
    pos.cutoffs["bindRBD"]=lloqs["bindRBD"]
    
    # data less than lod is set to lod/2 in the raw data
    llods["pseudoneutid50"]=2.612 # 40 * 0.0653
    lloqs["pseudoneutid50"]=51*0.0653 # 3.3303
    uloqs["pseudoneutid50"]=127411*0.0653 # 8319.938
    pos.cutoffs["pseudoneutid50"]=llods["pseudoneutid50"]
    
    llods["bindN"]=0.093744
    lloqs["bindN"]=4.4897
    uloqs["bindN"]=574.6783
    pos.cutoffs["bindN"]=23.4711
    
  } else if(TRIAL=="azd1222") {
    
    # data less than lod is set to lod/2
    llods["pseudoneutid50"]=2.612  
    lloqs["pseudoneutid50"]=56*0.0653 # 3.6568
    uloqs["pseudoneutid50"]=47806*0.0653 # 3121.732
    pos.cutoffs["pseudoneutid50"]=llods["pseudoneutid50"]
    
    # bindN info missing in SAP
    
  } else if(TRIAL=="azd1222_bAb") {
    
    # data less than lloq is set to lloq/2 in the raw data, Nexelis
    llods["bindSpike"]=NA 
    lloqs["bindSpike"]=62.8*0.0090 # 0.5652
    uloqs["bindSpike"]=238528.4*0.0090 # 2146.756
    pos.cutoffs["bindSpike"]=10.8424 # use same as COVE
    
  } else if(study_name=="HVTN705") {
    
    # get uloqs and lloqs from config
    # config$uloqs is a list before this processing
    if (!is.null(config$uloqs)) uloqs=sapply(config$uloqs, function(x) ifelse(is.numeric(x), x, Inf))  else uloqs=sapply(assays, function(a) Inf)
    if (!is.null(config$lloxs)) lloxs=sapply(config$lloxs, function(x) ifelse(is.numeric(x), x, NA))   else lloxs=sapply(assays, function(a) NA)
    lloqs=lloxs
    llods=lloxs
    names(uloqs)=assays # this is necessary because config$uloqs does not have names
    names(lloxs)=assays
    names(lloqs)=assays
    
  } else if(study_name=="PROFISCOV") { # Butantan
    
    # lod and lloq are the same
    # data less than lod is set to lloq/2
    
    #SARS-CoV-2 Spike           49 70,000 696 49
    #SARS-CoV-2 Spike (P.1)     32 36,000 463 32
    #SARS-CoV-2 Spike (B.1.351) 72 21,000 333 72
    #SARS-CoV-2 Spike (B.1.1.7) 70 47,000 712 70
    
    lloqs["bindSpike"] <- llods["bindSpike"] <- 49*0.0090 # 0.441
    uloqs["bindSpike"]=70000*0.0090 # 630
    pos.cutoffs["bindSpike"]=696*0.0090 # 15.0
    
    lloqs["bindSpike_P.1"] <- llods["bindSpike_P.1"] <- 32*0.0090 
    uloqs["bindSpike_P.1"]=36000*0.0090 
    pos.cutoffs["bindSpike_P.1"]=463*0.0090 
    
    lloqs["bindSpike_B.1.351"] <- llods["bindSpike_B.1.351"] <- 72*0.0090 
    uloqs["bindSpike_B.1.351"]=21000*0.0090 
    pos.cutoffs["bindSpike_B.1.351"]=333*0.0090 
    
    lloqs["bindSpike_B.1.1.7"] <- llods["bindSpike_B.1.1.7"] <- 70*0.0090 
    uloqs["bindSpike_B.1.1.7"]=47000*0.0090 
    pos.cutoffs["bindSpike_B.1.1.7"]=712*0.0090 
    
    #SARS-CoV-2 S1 RBD           35  30,000 1264 35
    #SARS-CoV-2 S1 RBD (P.1)     91  10,000 572  91
    #SARS-CoV-2 S1 RBD (B.1.351) 53  6,300  368  53
    #SARS-CoV-2 S1 RBD (B.1.1.7) 224 20,000 1111 224
    
    lloqs["bindRBD"] <- llods["bindRBD"] <- 35*0.0272 
    uloqs["bindRBD"]=30000*0.0272 # 630
    pos.cutoffs["bindRBD"]=1264*0.0272 # 15.0
    
    lloqs["bindRBD_P.1"] <- llods["bindRBD_P.1"] <- 91*0.0272 
    uloqs["bindRBD_P.1"]=10000*0.0272 
    pos.cutoffs["bindRBD_P.1"]=572*0.0272 
    
    lloqs["bindRBD_B.1.351"] <- llods["bindRBD_B.1.351"] <- 53*0.0272 
    uloqs["bindRBD_B.1.351"]=6300*0.0272 
    pos.cutoffs["bindRBD_B.1.351"]=368*0.0272 
    
    lloqs["bindRBD_B.1.1.7"] <- llods["bindRBD_B.1.1.7"] <- 224*0.0272 
    uloqs["bindRBD_B.1.1.7"]=20000*0.0272 
    pos.cutoffs["bindRBD_B.1.1.7"]=1111*0.0272 
    
    #SARS-CoV-2 Nucleocapsid 46 80,000 7015 46
    
    lloqs["bindN"] <- llods["bindN"] <- 46*0.00236 
    uloqs["bindN"]=80000*0.00236 
    pos.cutoffs["bindN"]=7015*0.00236 
    
    #LVMN
    llods["liveneutmn50"]=27.56 
    lloqs["liveneutmn50"]=27.84
    uloqs["liveneutmn50"]=20157.44 
    pos.cutoffs["liveneutmn50"]=llods["liveneutmn50"] 
    
  } else stop("unknown study_name 1")
  
  
  # llox is for plotting and can be either llod or lloq depending on trials
  if (is.null(lloxs)) {
    lloxs=ifelse(llox_labels=="LOD", llods[names(llox_labels)], lloqs[names(llox_labels)])
    lloxs=ifelse(llox_labels=="POS", pos.cutoffs[names(llox_labels)], lloxs)
  }
  
  # create assay_metadata from llods etc 
  assay_metadata = data.frame(assay=names(lloqs), lod=llods, lloq=lloqs, uloq=uloqs)
  # llox_label is treated differently b/c it may not contain bindN
  assay_metadata = cbinduneven(list(assay_metadata, llox_label=data.frame(llox_labels)))
  

}




# hacky fix for tabular, since unclear who else is using
# the truncated labels.assays.short later
labels.assays.short.tabular <- labels.assays.short

labels.time = config$time_labels
names(labels.time) = config$times

# axis labeling
labels.axis <- outer(rep("", length(times)), labels.assays.short[assays], "%.%")
labels.axis <- as.data.frame(labels.axis)
rownames(labels.axis) <- times

# title labeling
labels.title <- outer(labels.assays[assays], ": " %.% labels.time, paste0)
labels.title <- as.data.frame(labels.title)
colnames(labels.title) <- times
# NOTE: hacky solution to deal with changes in the number of markers
rownames(labels.title)[seq_along(assays)] <- assays
labels.title <- as.data.frame(t(labels.title))


do.fold.change=F
do.fold.change.overB=F #TRIAL %in% c("vat08_combined")


# if this flag is true, then the N IgG binding antibody is reported 
# in the immuno report (but is not analyzed in the cor or cop reports).
include_bindN <- !study_name %in% c("PREVENT19","AZD1222","VAT08m")



# COR-related config
if (exists("COR")) {
    myprint(COR)
    # making sure we are inadvertently using the wrong COR
    if(study_name=="ENSEMBLE") {
        if (contain(TRIAL, "EUA")) {
            # EUA datasets
            if (COR %in% c("D29","D29start1")) stop("For ENSEMBLE, we should not use D29 or D29start1")
        } 
    } 
    
    config.cor <- config::get(config = COR)
    stopifnot(!is.null(config.cor))
    
    if (startsWith(config.cor$tpeak%.%"","Delta")) { 
        tpeak = as.integer(strsplit(sub("Delta","",config.cor$tpeak), "over")[[1]][1])    
        do.fold.change=T    
    } else {
        tpeak=as.integer(paste0(config.cor$tpeak))
    }
    
    tpeaklag=as.integer(paste0(config.cor$tpeaklag))
    tinterm=as.integer(paste0(config.cor$tinterm))
    myprint(tpeak, tpeaklag, tinterm)
    # some config may not have all fields
    if (length(tpeak)==0 | length(tpeaklag)==0) stop("config "%.%COR%.%" misses some fields")

    if (do.fold.change) {
        all.markers=paste0(config.cor$tpeak, assays)
        names(all.markers)=all.markers
        
        all.markers.names.short=sub("\\(.+\\)", config.cor$tpeak, labels.assays.short)    # e.g. "Pseudovirus-nAb ID50 (IU50/ml)" => "Pseudovirus-nAb ID50 D57over29"
        names(all.markers.names.short)=all.markers
        
        all.markers.names.long=as.matrix(labels.title)[config.cor$tpeak, assays]
        names(all.markers.names.long)=all.markers
        
    } else {
        all.markers=paste0(DayPrefix, tpeak, assays)
        if (do.fold.change.overB) all.markers=c(all.markers, paste0("Delta", tpeak, "overB", assays))
        names(all.markers)=all.markers
        
        all.markers.names.short=c(
            labels.assays.short,
            if (do.fold.change.overB) sub("\\(.+\\)", "fold change", labels.assays.short) # e.g. "Pseudovirus-nAb ID50 (IU50/ml)" => "Pseudovirus-nAb ID50 fold change"
        )
        names(all.markers.names.short)=all.markers
        
        all.markers.names.long=c(
          as.matrix(labels.title)[DayPrefix%.%tpeak, assays],
          if (do.fold.change.overB) as.matrix(labels.title)["Delta"%.%tpeak%.%"overB", assays]
        )
        names(all.markers.names.long)=all.markers
    }
    
}
    
# to be deprecated
has57 = study_name %in% c("COVE","MockCOVE")
has29 = study_name %in% c("COVE","ENSEMBLE", "MockCOVE","MockENSEMBLE")



###################################################################################################
# read data

if (startsWith(tolower(study_name), "mock")) {
    data_name_updated <- paste0(TRIAL, "_data_processed_with_riskscore.csv")
    # the path depends on whether _common.R is sourced from Rmd or from R scripts in modules
    path_to_data = ifelse (endsWith(here::here(), "correlates_reporting2"), here::here("data_clean", data_name_updated), here::here("..", "data_clean", data_name_updated))
    data_name = data_name_updated    
} else {
    # the path depends on whether _common.R is sourced from Rmd or from R scripts in modules
    path_to_data = data_cleaned
    data_name = path_to_data
    # if path is relative, needs to do some processing
    if(endsWith(here::here(), "correlates_reporting2") & startsWith(path_to_data,"..")) path_to_data=substr(path_to_data, 4, nchar(path_to_data))
}
cat("Analysis-ready data: ", path_to_data, "\n")
# if this is run under _reporting level, it will not load. Thus we only warn and not stop
if (!file.exists(path_to_data)) stop ("_common.R: dataset with risk score not available ===========================================")

dat.mock <- read.csv(path_to_data)


###################################################################################################
# some more data processing

if (!DESCRIPTIVE & !EXPOSUREPROXIMAL) {
  
  # uloq censoring when it is for peak correlates analyses, not for descriptive analyses, or for exposure proximal correlates where decay model uses uncensored values
  for (a in assays) {
    uloq=uloqs[a]
    for (t in c(DayPrefix%.%timepoints)  ) {
      if ('t'%.%a %in% names(dat.mock)) {
        dat.mock[[t %.% a]] <- ifelse(dat.mock[[t %.% a]] > log10(uloq), log10(uloq), dat.mock[[t %.% a]])
      }
    }
    # process baseline marker if exists
    if ('B'%.%a %in% names(dat.mock)) {
      dat.mock[['B' %.% a]] <- ifelse(dat.mock[['B' %.% a]] > log10(uloq), log10(uloq), dat.mock[['B' %.% a]])
    }
  }    
  
}


if(TRIAL %in% c("janssen_pooled_partA", "janssen_na_partA", "janssen_la_partA", "janssen_sa_partA", "janssen_partA_VL")) {
  # make endpointDate.Bin a factor variable
  dat.mock$endpointDate.Bin = as.factor(dat.mock$endpointDate.Bin)

  
} else if (TRIAL %in% c("hvtn705secondRSA", "hvtn705secondNonRSA")) {
  # subset to RSA or non-RSA
  dat.mock = subset(dat.mock, RSA==ifelse(TRIAL=="hvtn705secondRSA", 1, 0))
  
  
} else if (study_name=='VAT08') {
  if (DESCRIPTIVE) {
    # censor bAb markers by lloq when the DESCRIPTIVE flag is set (data in analysis ready dataset is censored by lloq_min)
    bAb_markers = subset(assay_metadata, panel=='bindSpike' & assay!='bindSpike_mdw', assay, drop=T)
    for (a in bAb_markers) {
      for (t in c("B", paste0(DayPrefix, timepoints)) ) {
        dat.mock[[t%.%a]] = ifelse(dat.mock[[t%.%a]] < log10(lloqs[a]), log10(lloqs[a]/2), dat.mock[[t%.%a]])
      }
    }
    
    # recompute delta using censored markers
    # need to censor by uloq first as in data processing
    tmp=list()
    for (a in bAb_markers) {
      for (t in c("B", paste0(DayPrefix, timepoints)) ) {
        tmp[[t %.% a]] <- ifelse(dat.mock[[t %.% a]] > log10(uloqs[a]), log10(uloqs[a]), dat.mock[[t %.% a]])
      }
    }
    tmp=as.data.frame(tmp) # cannot subtract list from list, but can subtract data frame from data frame
    for (tp in rev(timepoints)) {
      dat.mock["Delta"%.%tp%.%"overB" %.% bAb_markers] <- tmp[DayPrefix%.%tp %.% bAb_markers] - tmp["B" %.% bAb_markers]
    }   
    dat.mock["Delta"%.%timepoints[2]%.%"over"%.%timepoints[1] %.% bAb_markers] <- tmp[DayPrefix%.% timepoints[2]%.% bAb_markers] - tmp[DayPrefix%.%timepoints[1] %.% bAb_markers]
  }
  
}


# converts discrete markers to factors from strings
if (TRIAL=="covail") {
  assays1 = c("pseudoneutid50_D614G", "pseudoneutid50_Delta", "pseudoneutid50_Beta", "pseudoneutid50_BA.1", "pseudoneutid50_BA.4.BA.5", "pseudoneutid50_MDW")
  assays1 = assays1%.%"cat"
  all.markers1 = c("B"%.%assays1, "Day15"%.%assays1, "Delta15overB"%.%assays1)
  for (a in all.markers1) {
    dat.mock[[a]] = as.factor(dat.mock[[a]])
  }
}



###################################################################################################
# For correlates reports

if (exists("COR")) {       
    
    # subset according to baseline seronegative status
    if (!is.null(config$Bserostatus)) {
        dat.mock=subset(dat.mock, Bserostatus==config$Bserostatus)
    }
    
    # for Novavax trial only, subset to US for the correlates modules
    # this is redundant in a way because only US participants have non-NA risk scores, but good to add
    if (study_name=="PREVENT19") dat.mock=subset(dat.mock, Country==0)
    
    # formula
    if (TRIAL %in% c("janssen_partA_VL")) {
      # will be defined in cor_coxph_ensemble_variant.R
      # form.0 is different for cox model and risk estimate
      # for risk estimate, it uses competing risk 

    } else {
      form.s = Surv(EventTimePrimary, EventIndPrimary) ~ 1
      form.0 = update (form.s, as.formula(config$covariates_riskscore))
      print(form.0)
    }
    
    ###########################################################
    # single time point COR config such as D29
    if (is.null(config.cor$tinterm)) {    
    
        dat.mock$ph1=dat.mock[[config.cor$ph1]]
        dat.mock$ph2=dat.mock[[config.cor$ph2]]
        dat.mock$Wstratum=dat.mock[[config.cor$WtStratum]]
        dat.mock$wt=dat.mock[[config.cor$wt]]
        dat.mock$EventIndPrimary =dat.mock[[config.cor$EventIndPrimary]]
        if (!is.null(config.cor$EventTimePrimary)) dat.mock$EventTimePrimary=dat.mock[[config.cor$EventTimePrimary]]
        if (!is.null(config.cor$tpsStratum)) dat.mock$tps.stratum=dat.mock[[config.cor$tpsStratum]]
        if (!is.null(config.cor$Earlyendpoint)) dat.mock$Earlyendpoint=dat.mock[[config.cor$Earlyendpoint]]
        
        # this day may be different from tpeak. it is the origin of followup days
        tpeak1 = as.integer(sub(".*[^0-9]+", "", config.cor$EventTimePrimary))
        
        # subset to require risk_score
        # check to make sure that risk score is not missing in ph1
        if(!is.null(dat.mock$risk_score)) {
            if (!TRIAL %in% c("janssen_pooled_EUA","janssen_na_EUA","janssen_na_partA")) { 
                # check this for backward compatibility
                stopifnot(nrow(subset(dat.mock, ph1 & is.na(risk_score)))==0)
            }
            dat.mock=subset(dat.mock, !is.na(risk_score))
        }        
    
        # data integrity checks
        if (!is.null(dat.mock$ph1)) {
            # missing values in variables that should have no missing values
            variables_with_no_missing <- paste0(c("ph2", "EventIndPrimary", "EventTimePrimary"))
            ans=sapply(variables_with_no_missing, function(a) all(!is.na(dat.mock[dat.mock$ph1==1, a])))
            if(!all(ans)) stop(paste0("Unexpected missingness in: ", paste(variables_with_no_missing[!ans], collapse = ", ")))   
            
            # ph1 should not have NA in Wstratum
            ans=with(subset(dat.mock,ph1==1), all(!is.na(Wstratum)))
            if(!ans) stop("Some Wstratum in ph1 are NA")
        } else {
            # may not be defined if COR is not provided in command line and used the default value
        }
        
        
        # define tfinal.tpeak
        if (TRIAL == "moderna_boost") {
            tfinal.tpeak = 92 # as computed in reporting3 repo
        } else if (TRIAL == "moderna_real" & COR == "D57a") {
          tfinal.tpeak = 92 # for comparing with stage 2 
        } else if (TRIAL == "janssen_na_EUA") {
            tfinal.tpeak=53
        } else if (TRIAL == "janssen_la_EUA") { # from day 48 to 58, risk jumps from .008 to .027
            tfinal.tpeak=48 
        } else if (TRIAL == "janssen_sa_EUA") {
            tfinal.tpeak=40            
        } else if (TRIAL == "janssen_pooled_EUA") {
            tfinal.tpeak=54
            
        } else if (startsWith(TRIAL, "janssen_") & endsWith(TRIAL, "partA")) {
          # smaller of the two: 1) last case in ph2 in vaccine, 2) last time to have 15 at risk in subcohort vaccine arm
          tfinal.tpeak=min(
            with(subset(dat.mock, Trt==1 & ph2 & EventIndPrimary==1), max(EventTimePrimary)),
            with(subset(dat.mock, Trt==1 & ph2 & SubcohortInd==1),    sort(EventTimePrimary, decreasing=T)[15]-1)
          )
          # for moderate, we choose to use the same tfinal.tpeak as the overall COVID
          if (COR=="D29ModerateIncludeNotMolecConfirmed") {
            tfinal.tpeak=min(
              with(subset(dat.mock, Trt==1 & ph2 & ModerateEventIndPrimaryIncludeNotMolecConfirmedD29==1), max(EventTimePrimary)),
              with(subset(dat.mock, Trt==1 & ph2 & SubcohortInd==1),    sort(EventTimePrimary, decreasing=T)[15]-1)
            )
          }
          
        } else if (TRIAL %in% c("janssen_partA_VL")) {
          # variant-specific tfinal.tpeak. set it to NULL so that it is not inadverdently used
          tfinal.tpeak = NULL 
          # smaller of the two: 1) last case in ph2 in vaccine, 2) last time to have 15 at risk in subcohort vaccine arm
          tfinal.tpeak.ls=list(
            US=list(
              # All=min(with(subset(dat.mock, Trt==1 & ph2 & EventIndPrimaryIncludeNotMolecConfirmedD29==1 & Region==0), max(EventTimePrimary)),
              #         with(subset(dat.mock, Trt==1 & ph2 & SubcohortInd==1 & Region==0),    sort(EventTimePrimary, decreasing=T)[15]-1)),
              
              Ancestral.Lineage=min(with(subset(dat.mock, Trt==1 & ph2 & EventIndPrimary==1 & Region==0 & seq1.variant=="Ancestral.Lineage"), max(EventTimePrimary)),
                                    with(subset(dat.mock, Trt==1 & ph2 & SubcohortInd==1 & Region==0),    sort(EventTimePrimary, decreasing=T)[15]-1))
            ),
            
            LatAm=list(
              # All=min(with(subset(dat.mock, Trt==1 & ph2 & EventIndPrimaryIncludeNotMolecConfirmedD29==1 & Region==1), max(EventTimePrimary)),
              #         with(subset(dat.mock, Trt==1 & ph2 & SubcohortInd==1 & Region==1),    sort(EventTimePrimary, decreasing=T)[15]-1)),
              
              Ancestral.Lineage = min(with(subset(dat.mock, Trt==1 & ph2 & EventIndPrimary==1 & Region==1 & seq1.variant=="Ancestral.Lineage"), max(EventTimePrimary)),
                        with(subset(dat.mock, Trt==1 & ph2 & SubcohortInd==1 & Region==1),    sort(EventTimePrimary, decreasing=T)[15]-1)),
              
              Gamma = min(with(subset(dat.mock, Trt==1 & ph2 & EventIndPrimary==1 & Region==1 & seq1.variant=="Gamma"), max(EventTimePrimary)),
                          with(subset(dat.mock, Trt==1 & ph2 & SubcohortInd==1 & Region==1),    sort(EventTimePrimary, decreasing=T)[15]-1)), 
              
              Lambda =min(with(subset(dat.mock, Trt==1 & ph2 & EventIndPrimary==1 & Region==1 & seq1.variant=="Lambda"), max(EventTimePrimary)),
                          with(subset(dat.mock, Trt==1 & ph2 & SubcohortInd==1 & Region==1),    sort(EventTimePrimary, decreasing=T)[15]-1)), 
              
              Mu    = min(with(subset(dat.mock, Trt==1 & ph2 & EventIndPrimary==1 & Region==1 & seq1.variant=="Mu"), max(EventTimePrimary)),
                        with(subset(dat.mock, Trt==1 & ph2 & SubcohortInd==1 & Region==1),    sort(EventTimePrimary, decreasing=T)[15]-1)),
              
              Zeta  = min(with(subset(dat.mock, Trt==1 & ph2 & EventIndPrimary==1 & Region==1 & seq1.variant=="Zeta"), max(EventTimePrimary)),
                          with(subset(dat.mock, Trt==1 & ph2 & SubcohortInd==1 & Region==1),    sort(EventTimePrimary, decreasing=T)[15]-1))
            ),
            
            RSA=list(
              # All=min(with(subset(dat.mock, Trt==1 & ph2 & EventIndPrimaryIncludeNotMolecConfirmedD29==1 & Region==2), max(EventTimePrimary)),
              #         with(subset(dat.mock, Trt==1 & ph2 & SubcohortInd==1 & Region==2),    sort(EventTimePrimary, decreasing=T)[15]-1)),
              
              Beta  = min(with(subset(dat.mock, Trt==1 & ph2 & EventIndPrimary==1 & Region==2 & seq1.variant=="Beta"), max(EventTimePrimary)),
                          with(subset(dat.mock, Trt==1 & ph2 & SubcohortInd==1 & Region==2),    sort(EventTimePrimary, decreasing=T)[15]-1))
            )
          )
          
        } else if (TRIAL %in% c("profiscov", "profiscov_lvmn")) {
            if (COR=="D91") tfinal.tpeak=66 else if(COR=="D43") tfinal.tpeak= 91+66-43 else stop("no tfinal.tpeak")
            
        } else if (study_name=="HVTN705") {
          tfinal.tpeak=550
          
        } else if (study_name=="VAT08") {
          # hardcode 180 days post dose 2
          if (COR=="D22M6omi" | COR=="D22M6ominAb") {
            tfinal.tpeak=180 # tpeak is dose 2
          } else if (COR=="D43M6omi" | COR=="D43M6ominAb") {
            tfinal.tpeak=180-21 # tpeak is 21 days post dose 2
          } # else is M12
          
        } else if (study_name=="IARCHPV") {
          tfinal.tpeak=NULL
          
        } else if (study_name=="COVAIL") {
          if (COR %in% c("D15to181","D92to181")) {
            tfinal.tpeak=181
          } else if (COR %in% c("D15to91")) {
            tfinal.tpeak=91
          } else {
            stop("COVAIL, wrong COR")
          }
          
        } else {
          # default rule for followup time is the last case in ph2 in vaccine arm
          tfinal.tpeak=with(subset(dat.mock, Trt==1 & ph2), max(EventTimePrimary[EventIndPrimary==1]))
        }

                
        if (!TRIAL %in% c("janssen_partA_VL", "vat08_combined", "id27hpv", "id27hpvnAb", "covail")) {
          # this block depends on tfinal.tpeak. For variants analysis, there is not just one tfinal.tpeak
          prev.vacc = get.marginalized.risk.no.marker(form.0, subset(dat.mock, Trt==1 & ph1), tfinal.tpeak)
          prev.plac = get.marginalized.risk.no.marker(form.0, subset(dat.mock, Trt==0 & ph1), tfinal.tpeak)   
          overall.ve = c(1 - prev.vacc/prev.plac) 
          myprint(prev.plac, prev.vacc, overall.ve)
        }

#        # get VE in the first month or two of followup
#        dat.tmp=dat.mock
#        t.tmp=30
#        # censor at t.tmp 
#        dat.tmp$EventIndPrimary =ifelse(dat.tmp$EventTimePrimary<=t.tmp, dat.tmp$EventIndPrimary, 0)
#        dat.tmp$EventTimePrimary=ifelse(dat.tmp$EventTimePrimary<=t.tmp, dat.tmp$EventTimePrimary, t.tmp)
#        prev.vacc = get.marginalized.risk.no.marker(form.0, subset(dat.tmp, Trt==1 & ph1), t.tmp)
#        prev.plac = get.marginalized.risk.no.marker(form.0, subset(dat.tmp, Trt==0 & ph1), t.tmp)
#        overall.ve = c(1 - prev.vacc/prev.plac)    
#        myprint(prev.plac, prev.vacc, overall.ve)        
        
    } else {
        # subset to require risk_score
        # note that it is assumed there no risk_score is missing for anyone in the analysis population. 
        # in the case of single time point COR, we do a check for that after definining ph1, 
        # which is why we have to do subset separately here again
        if(!is.null(dat.mock$risk_score)) dat.mock=subset(dat.mock, !is.na(risk_score))
        
    }
    
    
}

## wt can be computed from ph1, ph2 and Wstratum. See config for redundancy note
#wts_table <- dat.mock %>% dplyr::filter(ph1==1) %>% with(table(Wstratum, ph2))
#wts_norm <- rowSums(wts_table) / wts_table[, 2]
#dat.mock$wt <- wts_norm[dat.mock$Wstratum %.% ""]
#dat.mock$wt = ifelse(with(dat.mock, ph1), dat.mock$wt, NA) # the step above assigns weights for some subjects outside ph1. the next step makes them NA





###################################################################################################

# some common graphing parameters
if(config$is_ows_trial) {
    # maxed over Spike, RBD, N, restricting to Day 29 or 57
    if("bindSpike" %in% assays & "bindRBD" %in% assays) {
        if(has29) MaxbAbDay29 = max(dat.mock[,paste0("Day29", c("bindSpike", "bindRBD", "bindN"))], na.rm=T)
        if(has29) MaxbAbDelta29overB = max(dat.mock[,paste0("Delta29overB", c("bindSpike", "bindRBD", "bindN"))], na.rm=T)
        if(has57) MaxbAbDay57 = max(dat.mock[,paste0("Day57", c("bindSpike", "bindRBD", "bindN"))], na.rm=T)
        if(has57) MaxbAbDelta57overB = max(dat.mock[,paste0("Delta57overB", c("bindSpike", "bindRBD", "bindN"))], na.rm=T)
    }
        
    # maxed over ID50 and ID80, restricting to Day 29 or 57
    if("pseudoneutid50" %in% assays & "pseudoneutid80" %in% assays) {
        if(has29) MaxID50ID80Day29 = max(dat.mock[,paste0("Day29", c("pseudoneutid50", "pseudoneutid80"))], na.rm=T)
        if(has29) MaxID50ID80Delta29overB = max(dat.mock[,paste0("Delta29overB", c("pseudoneutid50", "pseudoneutid80"))], na.rm=TRUE)
        if(has57) MaxID50ID80Day57 = max(dat.mock[,paste0("Day57", c("pseudoneutid50", "pseudoneutid80"))], na.rm=T)        
        if(has57) MaxID50ID80Delta57overB = max(dat.mock[,paste0("Delta57overB", c("pseudoneutid50", "pseudoneutid80"))], na.rm=TRUE)
    }
    
            
}     


## map tps.stratum to stratification variables
#tps.stratums=sort(unique(dat.mock$tps.stratum)); names(tps.stratums)=tps.stratums
#decode.tps.stratum=t(sapply(tps.stratums, function(i) unlist(subset(dat.mock, tps.stratum==i)[1,
#    if (study_name=="COVE" | study_name=="MockCOVE" ) {
#        c("Senior", "HighRiskInd", "URMforsubcohortsampling")
#    } else if (study_name=="ENSEMBLE" | study_name=="MockENSEMBLE" ) {
#        c("Senior", "HighRiskInd", "Region", "URMforsubcohortsampling")
#    } else {
#        NA
#    }
#])))

    


###############################################################################
# figure labels and titles for markers
###############################################################################

# markers <- c(outer(times[which(times %in% c("B", "Day29", "Day57"))], assays, "%.%"))

# race labeling
labels.race <- c(
  "White", 
  "Black or African American",
  "Asian", 
  if ((study_name=="ENSEMBLE" | study_name=="MockENSEMBLE") & startsWith(TRIAL,"janssen_la")) "Indigenous South American" else "American Indian or Alaska Native",
  "Native Hawaiian or Other Pacific Islander", 
  "Multiracial",
  if ((study_name=="COVE" | study_name=="MockCOVE")) "Other", 
  "Not reported and unknown"
)

# ethnicity labeling
labels.ethnicity <- c(
  "Hispanic or Latino", "Not Hispanic or Latino",
  "Not reported and unknown"
)


#labels.assays <- c("Binding Antibody to Spike", 
#                   "Binding Antibody to RBD",
#                   "PsV Neutralization 50% Titer",
#                   "PsV Neutralization 80% Titer",
#                   "WT LV Neutralization 50% Titer")
#
#names(labels.assays) <- c("bindSpike", 
#                          "bindRBD", 
#                          "pseudoneutid50",
#                          "pseudoneutid80",
#                          "liveneutmn50")

#labels.assays.short <- c("Anti N IgG (BAU/ml)", 
#                         "Anti Spike IgG (BAU/ml)", 
#                         "Anti RBD IgG (BAU/ml)", 
#                         "Pseudovirus-nAb cID50", 
#                         "Pseudovirus-nAb cID80", 
#                         "Live virus-nAb cMN50")
#names(labels.assays.short) <- c("bindN",
#  "bindSpike",
#  "bindRBD",
#  "pseudoneutid50",
#  "pseudoneutid80",
#  "liveneutmn50")

#labels.time=c()
#for (t in times) {
#    labels.time=c(labels.time, "Day "%.%ifelse(t=="B", 1, t))
#}
#if (length(timepoints)==2) {
#    labels.time <- c("Day 1", "Day "%.%timepoints[1], "Day "%.%timepoints[2], 
#                     "D"%.%timepoints[1]%.%" fold-rise over D1", 
#                     "D"%.%timepoints[2]%.%" fold-rise over D1", 
#                     "D"%.%timepoints[2]%.%" fold-rise over D"%.%timepoints[1])
#    names(labels.time) <- c("B", DayPrefix%.%timepoints[1], DayPrefix%.%timepoints[2], 
#                        "Delta"%.%timepoints[1]%.%"overB", "Delta"%.%timepoints[2]%.%"overB", "Delta"%.%timepoints[2]%.%"over"%.%timepoints[1])
#} else {
#    labels.time <- c("Day 1", "Day "%.%timepoints[1], "D"%.%timepoints[1]%.%" fold-rise over D1")
#    names(labels.time) <- c("B", DayPrefix%.%timepoints[1], "Delta"%.%timepoints[1]%.%"overB")
#}




# baseline stratum labeling
if (study_name %in% c("COVE", "MockCOVE", "COVEBoost")) {
    Bstratum.labels <- c(
      "Age >= 65",
      "Age < 65, At risk",
      "Age < 65, Not at risk"
    )
    
} else if (study_name=="ENSEMBLE" | study_name=="MockENSEMBLE") {
    Bstratum.labels <- c(
      "Age < 60, Not at risk",
      "Age < 60, At risk",
      "Age >= 60, Not at risk",
      "Age >= 60, At risk"
    )
    
} else if (study_name %in% c("PREVENT19","AZD1222")) {
    Bstratum.labels <- c(
      "Age >= 65",
      "Age < 65"
    )

} else if (study_name %in% c("VAT08")) {
    Bstratum.labels <- c(
      "Age >= 60",
      "Age < 60"
    )

} else if (study_name=="HVTN705") {
    # do nothing

} else if (study_name %in% c("PROFISCOV")) {
    Bstratum.labels <- c("All")

} else if (study_name == 'IARCHPV') {
  Bstratum.labels <- c(
    "Age > 14",
    "Age <= 14"
  )
  
} else if (study_name=="COVAIL") {
  # do nothing
  
} else stop("unknown study_name 2")



# baseline stratum labeling
if (study_name %in% c("COVE", "MockCOVE", "COVEBoost")) {
  demo.stratum.labels <- c(
      "Age >= 65, URM",
      "Age < 65, At risk, URM",
      "Age < 65, Not at risk, URM",
      "Age >= 65, White non-Hisp",
      "Age < 65, At risk, White non-Hisp",
      "Age < 65, Not at risk, White non-Hisp"
    )
    
} else if (study_name=="ENSEMBLE" | study_name=="MockENSEMBLE") {
    demo.stratum.labels <- c(
      "US URM, Age 18-59, Not at risk",
      "US URM, Age 18-59, At risk",
      "US URM, Age >= 60, Not at risk",
      "US URM, Age >= 60, At risk",
      "US White non-Hisp, Age 18-59, Not at risk",
      "US White non-Hisp, Age 18-59, At risk",
      "US White non-Hisp, Age >= 60, Not at risk",
      "US White non-Hisp, Age >= 60, At risk",
      "Latin America, Age 18-59, Not at risk",
      "Latin America, Age 18-59, At risk",
      "Latin America, Age >= 60, Not at risk",
      "Latin America, Age >= 60, At risk",
      "South Africa, Age 18-59, Not at risk",
      "South Africa, Age 18-59, At risk",
      "South Africa, Age >= 60, Not at risk",
      "South Africa, Age >= 60, At risk"
    )
    
} else if (study_name=="PREVENT19") {
    demo.stratum.labels <- c(
      "US White non-Hisp, Age 18-64, Not at risk",
      "US White non-Hisp, Age 18-64, At risk",
      "US White non-Hisp, Age >= 65, Not at risk",
      "US White non-Hisp, Age >= 65, At risk",
      "US URM, Age 18-64, Not at risk",
      "US URM, Age 18-64, At risk",
      "US URM, Age >= 65, Not at risk",
      "US URM, Age >= 65, At risk",
      "Mexico, Age 18-64",
      "Mexico, Age >= 65"
    )

} else if (study_name=="AZD1222") {
    demo.stratum.labels <- c(
      "US White non-Hisp, Age 18-64",
      "US White non-Hisp, Age >= 65",
      "US URM, Age 18-64",
      "US URM, Age >= 65",
      "Non-US, Age 18-64",
      "Non-US, Age >= 65"
    )

} else if (study_name=="VAT08") {
  #    Stage 1, Not HND, Not senior
  #    Stage 1, Not HND, senior
  #    Stage 1, HND, Not senior
  #    Stage 1, HND, senior
  #    Stage 2, Not senior
  #    Stage 2, senior
  
    demo.stratum.labels <- c(
      "Stage 1 Not HND, Age 18-59",
      "Stage 1 Not HND, Age >= 60",
      "Stage 1 HND, Age 18-59",
      "Stage 1 HND, Age >= 60",
      "Stage 2, Age 18-59",
      "Stage 2, Age >= 60"
    )

} else if (study_name=="HVTN705") {
  # do nothing
  
} else if (study_name=="COVAIL") {
  # do nothing
  
} else if (study_name=="PROFISCOV") {
    demo.stratum.labels <- c("All")

} else if (study_name == 'IARCHPV') {
  demo.stratum.labels <- c(
    "Age > 14",
    "Age <= 14"
  )
  
} else stop("unknown study_name 3")

labels.regions.ENSEMBLE =c("0"="Northern America", "1"="Latin America", "2"="Southern Africa")
regions.ENSEMBLE=0:2
names(regions.ENSEMBLE)=labels.regions.ENSEMBLE

labels.countries.ENSEMBLE=c("0"="United States", "1"="Argentina", "2"="Brazil", "3"="Chile", "4"="Columbia", "5"="Mexico", "6"="Peru", "7"="South Africa")
countries.ENSEMBLE=0:7
names(countries.ENSEMBLE)=labels.countries.ENSEMBLE


###############################################################################
# theme options
###############################################################################

# fixed knitr chunk options
knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  out.width = "80%",
  out.extra = "",
  fig.pos = "H",
  fig.show = "hold",
  fig.align = "center",
  fig.width = 6,
  fig.asp = 0.618,
  fig.retina = 0.8,
  dpi = 600,
  echo = FALSE,
  message = FALSE,
  warning = FALSE
)

# global options
options(
  digits = 6,
  #scipen = 999,
  dplyr.print_min = 6,
  dplyr.print_max = 6,
  crayon.enabled = FALSE,
  bookdown.clean_book = TRUE,
  knitr.kable.NA = "NA",
  repos = structure(c(CRAN = "https://cran.rstudio.com/"))
)

# no complaints from installation warnings
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")

# overwrite options by output type
if (knitr:::is_html_output()) {
  #options(width = 80)

  # automatically create a bib database for R packages
  knitr::write_bib(c(
    .packages(), "bookdown", "knitr", "rmarkdown"
  ), "packages.bib")
}
if (knitr:::is_latex_output()) {
  #knitr::opts_chunk$set(width = 67)
  #options(width = 67)
  options(cli.unicode = TRUE)

  # automatically create a bib database for R packages
  knitr::write_bib(c(
    .packages(), "bookdown", "knitr", "rmarkdown"
  ), "packages.bib")
}

# create and set global ggplot theme
# borrowed from https://github.com/tidymodels/TMwR/blob/master/_common.R
theme_transparent <- function(...) {
  # use black-white theme as base
  ret <- ggplot2::theme_bw(...)

  # modify with transparencies
  trans_rect <- ggplot2::element_rect(fill = "transparent", colour = NA)
  ret$panel.background  <- trans_rect
  ret$plot.background   <- trans_rect
  ret$legend.background <- trans_rect
  ret$legend.key        <- trans_rect

  # always have legend below
  ret$legend.position <- "bottom"
  return(ret)
}

library(ggplot2)
theme_set(theme_transparent())
theme_update(
  text = element_text(size = 25),
  axis.text.x = element_text(colour = "black", size = 30),
  axis.text.y = element_text(colour = "black", size = 30)
)

# custom ggsave function with updated defaults
ggsave_custom <- function(filename = default_name(plot),
                          height= 15, width = 21, ...) {
  ggsave(filename = filename, height = height, width = width, ...)
}




