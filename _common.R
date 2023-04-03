#if (exists(".DEF.COMMON")) stop ("_common.R has already been loaded") else .DEF.COMMON=TRUE
library(methods)
library(dplyr)
library(kyotil)
library(marginalizedRisk)
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
    } else {
        warning("No COR. This is okay if _common.R is sourced just to load common functions. If needed, COR can be defined through command line argument or in R script before _common.R is sourced.")
    }
}


###################################################################################################
# read config

if(Sys.getenv("TRIAL")=="") stop(" *************************************  environmental variable TRIAL not defined  *************************************")

# TRIAL-related config
config <- config::get(config = Sys.getenv("TRIAL"))
if(length(config$llox_label)==1) {
    config$llox_label=rep(config$llox_label, length(config$assays))
} else {
    stopifnot(length(config$assays)==length(config$llox_label))
}
names(config$llox_label)=config$assays
for(opt in names(config)){
  eval(parse(text = paste0(names(config[opt])," <- config[[opt]]")))
}



# assays labels. This needs to come before all.markers
labels.assays = config$assay_labels
names(labels.assays) = config$assays

if (is.null(config$assay_labels_short)) {
    labels.assays.short=labels.assays
} else {
    labels.assays.short = config$assay_labels_short
    names(labels.assays.short) = config$assays
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

# creating short and long labels
#labels.assays.short <- labels.axis[1, ] # should not create this again
labels.assays.long <- labels.title

do.fold.change.overB=attr(config, "config") %in% c("vat08m_nonnaive")
do.fold.change=F

# if this flag is true, then the N IgG binding antibody is reported 
# in the immuno report (but is not analyzed in the cor or cop reports).
include_bindN <- !study_name %in% c("PREVENT19","AZD1222","VAT08m")



# COR-related config
if (exists("COR")) {
    myprint(COR)
    # making sure we are inadvertently using the wrong COR
    if(study_name=="ENSEMBLE") {
        if (contain(attr(config, "config"), "EUA")) {
            # EUA datasets
            if (COR %in% c("D29","D29start1")) stop("For ENSEMBLE, we should not use D29 or D29start1")
        } 
    } 
    
    config.cor <- config::get(config = COR)
    
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
        
        all.markers.names.long=as.matrix(labels.assays.long)[config.cor$tpeak, assays]
        names(all.markers.names.long)=all.markers
        
    } else {
        all.markers=paste0("Day", tpeak, assays)
        if (do.fold.change.overB) all.markers=c(all.markers, paste0("Delta", tpeak, "overB", assays))
        names(all.markers)=all.markers
        
        all.markers.names.short=c(
            labels.assays.short,
            if (do.fold.change.overB) sub("\\(.+\\)", "fold change", labels.assays.short) # e.g. "Pseudovirus-nAb ID50 (IU50/ml)" => "Pseudovirus-nAb ID50 fold change"
        )
        names(all.markers.names.short)=all.markers
        
        all.markers.names.long=c(
            as.matrix(labels.assays.long)["Day"%.%tpeak, assays],
            if (do.fold.change.overB) as.matrix(labels.assays.long)["Delta"%.%tpeak%.%"overB", assays]
        )
        names(all.markers.names.long)=all.markers
    }
    
}
    
# to be deprecated
has57 = study_name %in% c("COVE","MockCOVE")
has29 = study_name %in% c("COVE","ENSEMBLE", "MockCOVE","MockENSEMBLE")



###################################################################################################
# read data

data_name = paste0(attr(config, "config"), "_data_processed.csv")
if (startsWith(tolower(study_name), "mock")) {
    data_name_updated <- sub(".csv", "_with_riskscore.csv", data_name)
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

if(attr(config, "config") %in% c("janssen_pooled_partA", "janssen_na_partA", "janssen_la_partA", "janssen_sa_partA")) {
    # make endpointDate.Bin a factor variable
    dat.mock$endpointDate.Bin = as.factor(dat.mock$endpointDate.Bin)
}


###################################################################################################
# get marginalized risk without marker

get.marginalized.risk.no.marker=function(formula, dat, day){
    if (!is.list(formula)) {
        # model=T is required because the type of prediction requires it, see Note on ?predict.coxph
        fit.risk = coxph(formula, dat, model=T) 
        dat$EventTimePrimary=day
        risks = 1 - exp(-predict(fit.risk, newdata=dat, type="expected"))
        mean(risks)
    } else {
        # competing risk estimation
        out=pcr2(formula, dat, day)
        mean(out)
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
    
    # formulae
    form.s = Surv(EventTimePrimary, EventIndPrimary) ~ 1
    form.0 = update (form.s, as.formula(config$covariates_riskscore))
    print(form.0)
    
    comp.risk=FALSE
    
####Deprecated since we now take a simpler approach for the severe disease paper to treat severe endpoints as marginal survival endpoint
#    if (COR=="D29SevereIncludeNotMolecConfirmed") {
#        # formulae for competing risk
#        comp.risk=TRUE
#        form.0=list(
#            update(Surv(EventTimePrimaryIncludeNotMolecConfirmedD29, SevereEventIndPrimaryIncludeNotMolecConfirmedD29) ~ 1, 
#                as.formula(config$covariates_riskscore )),
#            update(Surv(EventTimePrimaryIncludeNotMolecConfirmedD29, ModerateEventIndPrimaryIncludeNotMolecConfirmedD29) ~ 1, 
#                as.formula(config$covariates_riskscore ))
#        )
#    }
    
#### Deprecated since we now use the hotdeck imputation approach
#    if (COR=="D29VL") {
#        # formulae for competing risk
#        comp.risk=TRUE
#        form.0=list(
#            update(Surv(EventTimePrimaryIncludeNotMolecConfirmedD29, EventIndPrimaryHasVLD29) ~ 1, 
#                as.formula(config$covariates_riskscore )),
#            update(Surv(EventTimePrimaryIncludeNotMolecConfirmedD29, EventIndPrimaryHasnoVLD29) ~ 1, 
#                as.formula(config$covariates_riskscore ))
#        )
#   }
    
    ###########################################################
    # single time point COR config such as D29
    if (is.null(config.cor$tinterm)) {    
    
        dat.mock$ph1=dat.mock[[config.cor$ph1]]
        dat.mock$ph2=dat.mock[[config.cor$ph2]]
        dat.mock$EventIndPrimary =dat.mock[[config.cor$EventIndPrimary]]
        dat.mock$EventTimePrimary=dat.mock[[config.cor$EventTimePrimary]]
        dat.mock$Wstratum=dat.mock[[config.cor$WtStratum]]
        dat.mock$wt=dat.mock[[config.cor$wt]]
        if (!is.null(config.cor$tpsStratum)) dat.mock$tps.stratum=dat.mock[[config.cor$tpsStratum]]
        if (!is.null(config.cor$Earlyendpoint)) dat.mock$Earlyendpoint=dat.mock[[config.cor$Earlyendpoint]]
        
        # this day may be different from tpeak. it is the origin of followup days
        tpeak1 = as.integer(sub(".*[^0-9]+", "", config.cor$EventTimePrimary))
        
        # subset to require risk_score
        # check to make sure that risk score is not missing in ph1
        if(!is.null(dat.mock$risk_score)) {
            if (!attr(config, "config") %in% c("janssen_na_EUA","janssen_na_partA")) { 
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
        
        # default rule for followup time is the last case in ph2 in vaccine arm
        tfinal.tpeak=with(subset(dat.mock, Trt==1 & ph2), max(EventTimePrimary[EventIndPrimary==1]))
        
        # exceptions
        if (attr(config, "config") == "janssen_na_EUA") {
            tfinal.tpeak=53
        } else if (attr(config, "config") == "janssen_la_EUA") { # from day 48 to 58, risk jumps from .008 to .027
            tfinal.tpeak=48 
        } else if (attr(config, "config") == "janssen_sa_EUA") {
            tfinal.tpeak=40            
        } else if (attr(config, "config") == "janssen_pooled_EUA") {
            tfinal.tpeak=54
            
        } else if (startsWith(attr(config, "config"), "janssen_") & contain(attr(config, "config"), "partA")) {
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

        } else if (attr(config, "config") %in% c("profiscov", "profiscov_lvmn")) {
            if (COR=="D91") tfinal.tpeak=66 else if(COR=="D43") tfinal.tpeak= 91+66-43
            
        } else if (study_name=="HVTN705") {
            tfinal.tpeak=550

        }
        
        prev.vacc = get.marginalized.risk.no.marker(form.0, subset(dat.mock, Trt==1 & ph1), tfinal.tpeak)
        prev.plac = get.marginalized.risk.no.marker(form.0, subset(dat.mock, Trt==0 & ph1), tfinal.tpeak)   
        overall.ve = c(1 - prev.vacc/prev.plac) 
        myprint(prev.plac, prev.vacc, overall.ve)
        

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

    



###################################################################################################

names(assays)=assays # add names so that lapply results will have names

# if the following part changes, make sure to copy to _common.R in the processing repo

# uloqs etc are hardcoded for ows trials but driven by config for other trials
# For bAb, IU and BAU are the same thing
# all values on BAU or IU
# LOQ can not be NA, it is needed for computing delta
pos.cutoffs<-llods<-lloqs<-uloqs<-c()
lloxs=NULL
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
    if (contain(attr(config, "config"), "EUA")) {
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
    
        
    } else if (contain(attr(config, "config"), "partA")) {
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
    
} else if(study_name=="AZD1222") {
       
    # data less than lloq is set to lloq/2 in the raw data, Nexelis
    llods["bindSpike"]=NA 
    lloqs["bindSpike"]=62.8*0.0090 # 0.5652
    uloqs["bindSpike"]=238528.4*0.0090 # 2146.756
    pos.cutoffs["bindSpike"]=10.8424 # use same as COVE
    
    # data less than lod is set to lod/2
    llods["pseudoneutid50"]=2.612  
    lloqs["pseudoneutid50"]=56*0.0653 # 3.6568
    uloqs["pseudoneutid50"]=47806*0.0653 # 3121.732
    pos.cutoffs["pseudoneutid50"]=llods["pseudoneutid50"]
    
    # bindN info missing in SAP
    
} else if(study_name=="VAT08m") { # Sanofi
       
    # data less than lod is set to lod/2
    llods["pseudoneutid50"]=2.612  
    lloqs["pseudoneutid50"]=95*0.0653 # 3.6568
    uloqs["pseudoneutid50"]=191429*0.0653 # 3121.732
    pos.cutoffs["pseudoneutid50"]=llods["pseudoneutid50"]
    
    llods["bindN"]=0.093744
    lloqs["bindN"]=4.4897
    uloqs["bindN"]=574.6783
    pos.cutoffs["bindN"]=23.4711
    
} else if(study_name=="HVTN705") {
    
    # get uloqs and lloqs from config
    # config$uloqs is a list before this processing
    if (!is.null(config$uloqs)) uloqs=sapply(config$uloqs, function(x) ifelse(is.numeric(x), x, Inf))  else uloqs=sapply(assays, function(a) Inf)
    if (!is.null(config$lloxs)) lloxs=sapply(config$lloxs, function(x) ifelse(is.numeric(x), x, NA))   else lloxs=sapply(assays, function(a) NA)
    names(uloqs)=assays # this is necessary because config$uloqs does not have names
    names(lloxs)=assays
    
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
if (is.null(lloxs)) lloxs=ifelse(config$llox_label=="LOD", llods[names(config$llox_label)], lloqs[names(config$llox_label)])




###############################################################################
# figure labels and titles for markers
###############################################################################

markers <- c(outer(times[which(times %in% c("B", "Day29", "Day57"))], assays, "%.%"))

# race labeling
labels.race <- c(
  "White", 
  "Black or African American",
  "Asian", 
  if ((study_name=="ENSEMBLE" | study_name=="MockENSEMBLE") & startsWith(attr(config, "config"),"janssen_la")) "Indigenous South American" else "American Indian or Alaska Native",
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
#    names(labels.time) <- c("B", "Day"%.%timepoints[1], "Day"%.%timepoints[2], 
#                        "Delta"%.%timepoints[1]%.%"overB", "Delta"%.%timepoints[2]%.%"overB", "Delta"%.%timepoints[2]%.%"over"%.%timepoints[1])
#} else {
#    labels.time <- c("Day 1", "Day "%.%timepoints[1], "D"%.%timepoints[1]%.%" fold-rise over D1")
#    names(labels.time) <- c("B", "Day"%.%timepoints[1], "Delta"%.%timepoints[1]%.%"overB")
#}




# baseline stratum labeling
if (study_name=="COVE" | study_name=="MockCOVE") {
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

} else if (study_name %in% c("VAT08m")) {
    Bstratum.labels <- c(
      "Age >= 60",
      "Age < 60"
    )

} else if (study_name=="HVTN705") {
    # do nothing

} else if (study_name %in% c("PROFISCOV")) {
    Bstratum.labels <- c("All")

} else stop("unknown study_name 2")



# baseline stratum labeling
if (study_name=="COVE" | study_name=="MockCOVE") {
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

} else if (study_name=="VAT08m") {
#    demo.stratum.labels <- c(
#      "Not HND, Age 18-59",
#      "Not HND, Age >= 60",
#      "HND, Age 18-59",
#      "HND, Age >= 60",
#      "USA, Age 18-59",
#      "USA, Age >= 60",
#      "JPN, Age 18-59",
#      "JPN, Age >= 60"
#    )

    # in this partial dataset, we need to collapse "Not HND, US or JPN, senior" and "HND, senior" due to sparsity
    demo.stratum.labels <- c(
      "Not HND, Age 18-59",
      "Not USA or JPN, Age >= 60",
      "HND, Age 18-59",
      "USA, Age 18-59",
      "USA, Age >= 60",
      "JPN, Age 18-59",
      "JPN, Age >= 60"
    )

} else if (study_name=="HVTN705") {
    # do nothing

} else if (study_name=="PROFISCOV") {
    demo.stratum.labels <- c("All")

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




############## Utility func

# e.g. Day22pseudoneutid50 => pseudoneutid50, Delta22overBpseudoneutid50 => pseudoneutid50
get.assay.from.name=function(a) {
    if (startsWith(a,"Day")) {
        sub("Day[[0123456789]+", "", a)
    } else if (contain(a,"overB")) {
        sub("Delta[[0123456789]+overB", "", a)
    } else if (contain(a,"over")) {
        sub("Delta[[0123456789]+over[[0123456789]+", "", a)
    } else stop("get.assay.from.name: not sure what to do")
}


get.range.cor=function(dat, assay, time) {
    if(assay %in% c("bindSpike", "bindRBD") & all(c("pseudoneutid50", "pseudoneutid80") %in% assays)) {
        ret=range(dat[["Day"%.%time%.%"bindSpike"]], 
                  dat[["Day"%.%time%.%"bindRBD"]], 
                  log10(lloxs[c("bindSpike","bindRBD")]/2), na.rm=T)
        
    } else if(assay %in% c("pseudoneutid50", "pseudoneutid80") & all(c("pseudoneutid50", "pseudoneutid80") %in% assays)) {
        ret=range(dat[["Day"%.%time%.%"pseudoneutid50"]], 
                  dat[["Day"%.%time%.%"pseudoneutid80"]], 
                  #log10(uloqs[c("pseudoneutid50","pseudoneutid80")]),
                  log10(lloxs[c("pseudoneutid50","pseudoneutid80")]/2), na.rm=T) 
    } else {
        ret=range(dat[["Day"%.%time%.%assay]], 
        log10(lloxs[assay]/2), na.rm=T)        
    }
    delta=(ret[2]-ret[1])/20     
    c(ret[1]-delta, ret[2]+delta)
}

draw.x.axis.cor=function(xlim, llox, llox.label){
        
    xx=seq(ceiling(xlim[1]), floor(xlim[2]))        
    if (is.na(llox)) {
        for (x in xx) {
            axis(1, at=x, labels=if (x>=3) bquote(10^.(x)) else 10^x )    
        }
    } else if (llox.label=="delta") {
        for (x in xx) {
            axis(1, at=x, labels=if (x>=3 | x<=-3) bquote(10^.(x)) else 10^x )    
        }    
    } else {
        axis(1, at=log10(llox), labels=llox.label)
        for (x in xx[xx>log10(llox*1.8)]) {
            axis(1, at=x, labels= if(x>=3) bquote(10^.(x)) else 10^x)
        }
    }
    
    # add e.g. 30 between 10 and 100
    if (length(xx)<=3 & length(xx)>1) { 
        # a hack for prevent19 ID50 to not draw 3 b/c it is too close to LOD
        tmp=2:length(xx)
        if (study_name=="PREVENT19") tmp=3:length(xx)
        for (i in tmp) {
            x=xx[i-1]
            axis(1, at=x+log10(3), labels=if (x>=3) bquote(3%*%10^.(x)) else 3*10^x )
        }
    }
    
}

##### Copy of draw.x.axis.cor but returns the x-axis ticks and labels
# This is necessary if one works with ggplot as the "axis" function does not work.
get.labels.x.axis.cor=function(xlim, llox){
  xx=seq(floor(xlim[1]), ceiling(xlim[2]))
  if (!is.na(llox)) xx=xx[xx>log10(llox*2)]
  x_ticks <- xx
  if (is.na(llox)) {
      labels <- sapply(xx, function(x) {
        if (x>=3) bquote(10^.(x)) else 10^x
      })
  } else {
      labels <- sapply(xx, function(x) {
        if (log10(llox)==x) config$llox_label else if (x>=3) bquote(10^.(x)) else 10^x
      })
      #if(!any(log10(llox)==x_ticks)){
        x_ticks <- c(log10(llox), x_ticks)
        labels <- c(config$llox_label, labels)
      #}
  }
  return(list(ticks = x_ticks, labels = labels))
}


# bootstrap from case control studies is done by resampling cases, ph2 controls, and non-ph2 controls separately. 
# Across bootstrap replicates, the number of cases does not stay constant, neither do the numbers of ph2 controls by demographics strata. 
# Specifically,
# 1) sample with replacement to get dat.b. From this dataset, take the cases and count ph2 and non-ph2 controls by strata
# 2) sample with replacement ph2 and non-ph2 controls by strata
bootstrap.case.control.samples=function(dat.ph1, seed, delta.name="EventIndPrimary", strata.name="tps.stratum", ph2.name="ph2", min.cell.size=1) {
#dat.ph1=dat.tmp; delta.name="EventIndPrimary"; strata.name="tps.stratum"; ph2.name="ph2"; min.cell.size=0
    
    set.seed(seed)
    
    dat.tmp=data.frame(ptid=1:nrow(dat.ph1), delta=dat.ph1[,delta.name], strata=dat.ph1[,strata.name], ph2=dat.ph1[,ph2.name])
    
    nn.ph1=with(dat.tmp, table(strata, delta))
    strat=rownames(nn.ph1); names(strat)=strat
    # ctrl.ptids is a list of lists
    ctrl.ptids = with(subset(dat.tmp, delta==0), lapply(strat, function (i) list(ph2=ptid[strata==i & ph2], nonph2=ptid[strata==i & !ph2])))
    
    # 1. resample dat.ph1 to get dat.b, but only take the cases 
    dat.b=dat.tmp[sample.int(nrow(dat.tmp), r=TRUE),]
    
    # re-do resampling if the bootstrap dataset has too few samples in a cell in nn.ctrl.b
    while(TRUE) {   
        nn.ctrl.b=with(subset(dat.b, !delta), table(strata, ph2))
        if (min(nn.ctrl.b)<min.cell.size | ncol(nn.ctrl.b)<2) dat.b=dat.tmp[sample.int(nrow(dat.tmp), r=TRUE),] else break
    }
    
    # take the case ptids
    case.ptids.b = dat.b$ptid[dat.b$delta==1]
    
    # 2. resample controls in dat.ph1 (numbers determined by dat.b) stratified by strata and ph2/nonph2
    # ph2 and non-ph2 controls by strata
    nn.ctrl.b=with(subset(dat.b, !delta), table(strata, ph2))
    # sample the control ptids
    ctrl.ptids.by.stratum.b=lapply(strat, function (i) {
        c(sample(ctrl.ptids[[i]]$ph2, nn.ctrl.b[i,2], r=T),
          sample(ctrl.ptids[[i]]$nonph2, nn.ctrl.b[i,1], r=T))
    })
    ctrl.ptids.b=do.call(c, ctrl.ptids.by.stratum.b)    
    
    # return data frame
    dat.ph1[c(case.ptids.b, ctrl.ptids.b), ]
}

## testing
#dat.b=bootstrap.case.control.samples(dat.vac.seroneg)
#with(dat.vac.seroneg, table(ph2, tps.stratum, EventIndPrimary))
#with(dat.b, table(ph2, tps.stratum, EventIndPrimary))
#> with(dat.vac.seroneg, table(ph2, tps.stratum, EventIndPrimary))
#, , EventIndPrimary = 0
#
#       tps.stratum
#ph2       33   34   35   36   37   38   39   40   41   42   43   44   45   46   47   48
#  FALSE 1483  915  759  439 1677 1138  894  591 3018 1973 1559 1051 1111  693  511  329
#  TRUE    57   53   55   57   56   57   57   56   58   55   55   57   57   56   56   56
#
#, , EventIndPrimary = 1
#
#       tps.stratum
#ph2       33   34   35   36   37   38   39   40   41   42   43   44   45   46   47   48
#  FALSE    1    0    0    1    0    1    0    0    2    1    2    1    0    0    0    1
#  TRUE     3    7    7   10    8   11    2   13   17   23   15   23    5    6    4    6
#
#> with(dat.b, table(ph2, tps.stratum, EventIndPrimary))
#, , EventIndPrimary = 0
#
#       tps.stratum
#ph2       33   34   35   36   37   38   39   40   41   42   43   44   45   46   47   48
#  FALSE 1487  911  750  462 1675 1181  884  570 3058 2023 1499 1034 1094  694  487  329
#  TRUE    47   57   65   62   50   53   50   64   55   61   65   53   64   53   54   60
#
#, , EventIndPrimary = 1
#
#       tps.stratum
#ph2       33   34   35   36   37   38   39   40   41   42   43   44   45   46   47   48
#  FALSE    0    0    0    0    0    2    0    0    1    1    3    3    0    0    0    2
#  TRUE     2    6    8    5    9   13    0   11   20   26   10   20    4    3    4    5


# for bootstrap use
get.ptids.by.stratum.for.bootstrap = function(data) {
    strat=sort(unique(data$tps.stratum))
    ptids.by.stratum=lapply(strat, function (i) 
        list(subcohort=subset(data, tps.stratum==i & SubcohortInd==1, Ptid, drop=TRUE), nonsubcohort=subset(data, tps.stratum==i & SubcohortInd==0, Ptid, drop=TRUE))
    )    
    # add a pseudo-stratum for subjects with NA in tps.stratum (not part of Subcohort). 
    # we need this group because it contains some cases with missing tps.stratum
    # if data is ph2 only, then this group is only cases because ph2 = subcohort + cases
    tmp=list(subcohort=subset(data, is.na(tps.stratum), Ptid, drop=TRUE),               nonsubcohort=NULL)
    ptids.by.stratum=append(ptids.by.stratum, list(tmp))    
    ptids.by.stratum
}


# bootstrap case cohort samples
# data is assumed to contain only ph1 ptids
get.bootstrap.data.cor = function(data, ptids.by.stratum, seed) {
    set.seed(seed)    
    
    # For each sampling stratum, bootstrap samples in subcohort and not in subchort separately
    tmp=lapply(ptids.by.stratum, function(x) c(sample(x$subcohort, r=TRUE), sample(x$nonsubcohort, r=TRUE)))
    
    dat.b=data[match(unlist(tmp), data$Ptid),]
    
    # compute weights
    tmp=with(dat.b, table(Wstratum, ph2))
    weights=rowSums(tmp)/tmp[,2]
    dat.b$wt=weights[""%.%dat.b$Wstratum]
    # we assume data only contains ph1 ptids, thus weights is defined for every bootstrapped ptids
    
    dat.b
}

# extract assay from marker name such as Day57pseudoneutid80, Bpseudoneutid80
marker.name.to.assay=function(marker.name) {
    if(endsWith(marker.name, "bindSpike")) {
        "bindSpike"
    } else if(endsWith(marker.name, "bindRBD")) {
        "bindRBD"
    } else if(endsWith(marker.name, "bindN")) {
        "bindN"
    } else if(endsWith(marker.name, "pseudoneutid50")) {
        "pseudoneutid50"
    } else if(endsWith(marker.name, "pseudoneutid80")) {
        "pseudoneutid80"
    } else if(endsWith(marker.name, "liveneutmn50")) {
        "liveneutmn50"
    } else stop("marker.name.to.assay: wrong marker.name")
}


# x is the marker values
# assay is one of assays, e.g. pseudoneutid80
report.assay.values=function(x, assay){
    lars.quantiles=seq(0,1,length.out=30) [round(seq.int(1, 30, length.out = 10))]
    sens.quantiles=c(0.15, 0.85)
    # cannot have different lengths for different assays, otherwise downstream code may break
    fixed.values = log10(c("500"=500, "1000"=1000))
    # if we want to add "llox/2"=unname(lloxs[assay]/2))) to fixed.values, we have to get assay right, which will take some thought because marker.name.to.assay is hardcoded
    out=sort(c(quantile(x, c(lars.quantiles,sens.quantiles), na.rm=TRUE), fixed.values[fixed.values<max(x, na.rm=T) & fixed.values>min(x, na.rm=T)]))    
    out
    #out[!duplicated(out)] # unique strips away the names. But don't take out duplicates because 15% may be needed and because we may want the same number of values for each assay
}
#report.assay.values (dat.vac.seroneg[["Day57pseudoneutid80"]], "pseudoneutid80")


add.trichotomized.markers=function(dat, markers, wt.col.name) {
    
    if(verbose) print("add.trichotomized.markers ...")
    
    marker.cutpoints <- list()    
    for (a in markers) {
        if (verbose) myprint(a, newline=F)
        tmp.a=dat[[a]]
        
        # if we estimate cutpoints using all non-NA markers, it may have an issue when a lot of subjects outside ph2 have non-NA markers
        # since that leads to uneven distribution of markers between low/med/high among ph2
        # this issue did not affect earlier trials much, but it is a problem with vat08m. We are changing the code for trials after vat08m
        if (attr(config, "config") %in% c("hvtn705","hvtn705V1V2","hvtn705second","hvtn705secondprimary","moderna_real","moderna_mock","prevent19",
                "janssen_pooled_EUA","janssen_na_EUA","janssen_la_EUA","janssen_sa_EUA")) {
            flag=rep(TRUE, length(tmp.a))
        } else {
            flag=dat$ph2
        }
    
        if(startsWith(a, "Day")) {
            # not fold change
            uppercut=log10(uloqs[get.assay.from.name(a)]); uppercut=uppercut*ifelse(uppercut>0,.9999,1.0001)
            lowercut=min(tmp.a, na.rm=T)*1.0001; lowercut=lowercut*ifelse(lowercut>0,1.0001,.9999)
            if (mean(tmp.a>uppercut, na.rm=T)>1/3) {
                # if more than 1/3 of vaccine recipients have value > ULOQ, let q.a be (median among those < ULOQ, ULOQ)
                if (verbose) cat("more than 1/3 of vaccine recipients have value > ULOQ\n")
                q.a=c(wtd.quantile(tmp.a[dat[[a]]<=uppercut & flag], weights = dat[[wt.col.name]][tmp.a<=uppercut & flag], probs = c(1/2)),  uppercut)
            } else if (mean(tmp.a<lowercut, na.rm=T)>1/3) {
                # if more than 1/3 of vaccine recipients have value at min, let q.a be (min, median among those > LLOQ)
                if (verbose) cat("more than 1/3 of vaccine recipients have at min\n")
                q.a=c(lowercut, wtd.quantile(tmp.a[dat[[a]]>=lowercut & flag], weights = dat[[wt.col.name]][tmp.a>=lowercut & flag], probs = c(1/2))  )
            } else {
                # this implementation uses all non-NA markers, which include a lot of subjects outside ph2, and that leads to uneven distribution of markers between low/med/high among ph2
                #q.a <- wtd.quantile(tmp.a, weights = dat[[wt.col.name]], probs = c(1/3, 2/3))
                q.a <- wtd.quantile(tmp.a[flag], weights = dat[[wt.col.name]][flag], probs = c(1/3, 2/3))
            }
        } else {
            # fold change
            q.a <- wtd.quantile(tmp.a[flag], weights = dat[[wt.col.name]][flag], probs = c(1/3, 2/3))
        }
        tmp=try(factor(cut(tmp.a, breaks = c(-Inf, q.a, Inf))), silent=T)
 
        do.cut=FALSE # if TRUE, use cut function which does not use weights
        # if there is a huge point mass, an error would occur, or it may not break into 3 groups
        if (inherits(tmp, "try-error")) do.cut=TRUE else if(length(table(tmp)) != 3) do.cut=TRUE
        
        if(!do.cut) {
            dat[[a %.% "cat"]] <- tmp
            marker.cutpoints[[a]] <- q.a
        } else {
            cat("\nfirst cut fails, call cut again with breaks=3 \n")
            # cut is more robust but it does not incorporate weights
            tmp=cut(tmp.a, breaks=3)
            stopifnot(length(table(tmp))==3)
            dat[[a %.% "cat"]] = tmp
            # extract cut points from factor level labels
            tmpname = names(table(tmp))[2]
            tmpname = substr(tmpname, 2, nchar(tmpname)-1)
            marker.cutpoints[[a]] <- as.numeric(strsplit(tmpname, ",")[[1]])
        }
        stopifnot(length(table(dat[[a %.% "cat"]])) == 3)
        if(verbose) {
            print(table(dat[[a %.% "cat"]]))
            cat("\n")
        }
    }
    
    attr(dat, "marker.cutpoints")=marker.cutpoints
    dat
    
}



# a function to print tables of cases counts with different marker availability
# note that D57 cases and intercurrent cases may add up to more than D29 cases because ph1.D57 requires EarlyendpointD57==0 while ph1.D29 requires EarlyendpointD29==0
make.case.count.marker.availability.table=function(dat) {
    if (study_name=="COVE" | study_name=="MockCOVE" ) {
        idx.trt=1:0
        names(idx.trt)=c("vacc","plac")
        cnts = sapply (idx.trt, simplify="array", function(trt) {
             idx=1:3
             names(idx)=c("Day 29 Cases", "Day 57 Cases", "Intercurrent Cases")
             tab=t(sapply (idx, function(i) {           
                tmp.1 = with(subset(dat, Trt==trt & Bserostatus==0 & (if(i==2) EventIndPrimaryD57 else EventIndPrimaryD29) &   (if(i==2) ph1.D57 else if(i==1) ph1.D29 else ph1.intercurrent.cases)), is.na(BbindSpike)     | is.na(BbindRBD) )
                tmp.2 = with(subset(dat, Trt==trt & Bserostatus==0 & (if(i==2) EventIndPrimaryD57 else EventIndPrimaryD29) &   (if(i==2) ph1.D57 else if(i==1) ph1.D29 else ph1.intercurrent.cases)), is.na(Day29bindSpike) | is.na(Day29bindRBD))
                tmp.3 = with(subset(dat, Trt==trt & Bserostatus==0 & (if(i==2) EventIndPrimaryD57 else EventIndPrimaryD29) &   (if(i==2) ph1.D57 else if(i==1) ph1.D29 else ph1.intercurrent.cases)), is.na(Day57bindSpike) | is.na(Day57bindRBD))    
                
                c(sum(tmp.1 & tmp.2 & tmp.3), sum(tmp.1 & tmp.2 & !tmp.3), sum(tmp.1 & !tmp.2 & tmp.3), sum(tmp.1 & !tmp.2 & !tmp.3), 
                  sum(!tmp.1 & tmp.2 & tmp.3), sum(!tmp.1 & tmp.2 & !tmp.3), sum(!tmp.1 & !tmp.2 & tmp.3), sum(!tmp.1 & !tmp.2 & !tmp.3))
            }))
            colnames(tab)=c("---", "--+", "-+-", "-++", "+--", "+-+", "++-", "+++")
            tab
        })
        cnts
    } else if (study_name=="ENSEMBLE" | study_name=="MockENSEMBLE" ) {
        idx.trt=1:0
        names(idx.trt)=c("vacc","plac")
        cnts = sapply (idx.trt, simplify="array", function(trt) {
             idx=1:1
             tab=t(sapply (idx, function(i) {           
                tmp.1 = with(subset(dat, Trt==trt & Bserostatus==0 & if(i==2) EventIndPrimaryD57 else EventIndPrimaryD29 &   if(i==2) ph1.D57 else if(i==1) ph1.D29 else ph1.intercurrent.cases), is.na(BbindSpike)     | is.na(BbindRBD) )
                tmp.2 = with(subset(dat, Trt==trt & Bserostatus==0 & if(i==2) EventIndPrimaryD57 else EventIndPrimaryD29 &   if(i==2) ph1.D57 else if(i==1) ph1.D29 else ph1.intercurrent.cases), is.na(Day29bindSpike) | is.na(Day29bindRBD))
                
                c(sum(tmp.1 & tmp.2), sum(!tmp.1 & tmp.2), sum(tmp.1 & !tmp.2), sum(!tmp.1 & !tmp.2))
             }))
             colnames(tab)=c("--", "+-", "-+", "++")
             tab
        })
        t(drop(cnts))
    } else {
        NA
    }
}
#make.case.count.marker.availability.table(dat.mock)


# get histogram object to add to VE plots etc
get.marker.histogram=function(marker, wt, trial, marker.break=marker) {
    # first call hist to get breaks, then call weighted.hist
    tmp.1=hist(marker.break,breaks=ifelse(trial=="moderna_real",25,15),plot=F)  # 15 is treated as a suggestion and the actual number of breaks is determined by pretty()
    tmp=weighted.hist(marker,wt, breaks=tmp.1$breaks, plot=F)
    attr(tmp,"class")="histogram" 
    tmp
}
