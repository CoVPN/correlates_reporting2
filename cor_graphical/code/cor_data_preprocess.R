#Sys.setenv(TRIAL = "vat08_combined") D22D43omi
#D22vat08_combined_M12_bAb, D22vat08_combined_M12_nAb, D43vat08_combined_M12_bA, D43vat08_combined_M12_nAb
#Sys.setenv(TRIAL = "id27hpv") id27hpvnAb
#Sys.setenv(TRIAL = "janssen_partA_VL")
#Sys.setenv(TRIAL = "janssen_pooled_partA")
#Sys.setenv(TRIAL = "prevent19_stage2") # D35prevent19_stage2_delta, D35prevent19_stage2_severe
#Sys.setenv(TRIAL = "azd1222_stage2") # D57azd1222_stage2_delta_nAb, D57azd1222_stage2_delta_bAb, D57azd1222_stage2_severe_nAb, D57azd1222_stage2_severe_bAb
#Sys.setenv(TRIAL = "nvx_uk302") # D35nvx_uk302
#Sys.setenv(TRIAL = "prevent19nvx") # D35prevent19nvx
#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
Sys.setenv(DESCRIPTIVE = 1)
source(here::here("..", "_common.R"))
#-----------------------------------------------

source(here::here("code", "cor_process_function.R"))
library(here)
library(dplyr)
library(tidyverse)
library(stringr)
if (!is.null(config$assay_metadata)) {pos.cutoffs = assay_metadata$pos.cutoff; names(pos.cutoffs) <- assays}

## COR has a set of analysis-specific parameters defined in the config file
config.cor <- config::get(config = COR)
if (study_name=="VAT08") {dat_proc = dat_proc %>% filter(Trialstage == 1)} # need manually update "Trialstage" and line 69 in report.Rmd

# forcing this is not a good idea. ~ Youyi
# set wt.DXX missingness to 0
wt.vars <- colnames(dat_proc)[grepl("wt.D", colnames(dat_proc))]
for (a in wt.vars) dat_proc[a][is.na(dat_proc[a])]<-0

if(T){ # for ENSEMBLE SA and LA reports only
  # copy Bpseudoneutid50 to Bpseudoneutid50la & calculate delta value if Day29pseudoneutid50la exists and is required for reporting
  # copy Bpseudoneutid50 to Bpseudoneutid50sa & calculate delta value if Day29pseudoneutid50sa exists and is required for reporting
  if ("Day29pseudoneutid50la" %in% colnames(dat_proc) & "pseudoneutid50la" %in% assays) {
    dat_proc$Bpseudoneutid50la = dat_proc$Bpseudoneutid50
    dat_proc$Delta29overBpseudoneutid50la = pmin(log10(uloqs["pseudoneutid50la"]), dat_proc$Day29pseudoneutid50la) - pmin(log10(uloqs["pseudoneutid50la"]), dat_proc$Bpseudoneutid50la)
    }
  if ("Day29pseudoneutid50sa" %in% colnames(dat_proc) & "pseudoneutid50sa" %in% assays) {
    dat_proc$Bpseudoneutid50sa = dat_proc$Bpseudoneutid50
    dat_proc$Delta29overBpseudoneutid50sa = pmin(log10(uloqs["pseudoneutid50sa"]), dat_proc$Day29pseudoneutid50sa) - pmin(log10(uloqs["pseudoneutid50sa"]), dat_proc$Bpseudoneutid50sa)
    }
}

if(F){# for simulated figure Peter requested on 8/2/2022, hardly lower the readouts for AZ PsV in dat_proc
  set.seed(12345)
  dat_proc <- dat_proc %>%
    mutate(randnum=runif(min=1, max=2, n()),
           Day57pseudoneutid50 = ifelse(ph2.D57==1 & EventIndPrimaryD57==1 & Day57pseudoneutid50>2,  # cases
                                        Day57pseudoneutid50-randnum, 
                                        ifelse(ph2.D57==1 & EventIndPrimaryD57==1 & Day57pseudoneutid50>0.8,  # cases
                                               Day57pseudoneutid50-randnum/2,
                                          ifelse(ph2.D57==1 & AnyinfectionD1==0 & EventIndPrimaryD1==0 & Day57pseudoneutid50<1.2, # non-cases
                                                 Day57pseudoneutid50+randnum, 
                                                 ifelse(ph2.D57==1 & AnyinfectionD1==0 & EventIndPrimaryD1==0 & Day57pseudoneutid50<1.8,
                                                        Day57pseudoneutid50+randnum/2,
                                                        Day57pseudoneutid50)))),
           Day57pseudoneutid50 = ifelse(Day57pseudoneutid50<log10(2.612/2), log10(2.612/2), Day57pseudoneutid50),
           # hard code here to remove one data point
           Day57pseudoneutid50 = ifelse(ph2.D57==1 & EventIndPrimaryD57==1 & Day57pseudoneutid50>1.2, 1.03, Day57pseudoneutid50)
    )
}

# load parameters
source(here("code", "params.R"))


################################################
dat <- as.data.frame(dat_proc); #dat$ph2.D43 = dat$ph2.D43.original; dat$wt.D22 = dat$wt.D22.original; dat$wt.D43 = dat$wt.D43.original

Args <- commandArgs(trailingOnly=TRUE)
COR=Args[1]
if (grepl("IncludeNotMolecConfirmed", COR)) {incNotMol <- "IncludeNotMolecConfirmed"
} else {incNotMol <- ""}

# set EventIndTimePrimary to EventIndTimeOmicron if study_name=="VAT08_combined" & COR=="D22D43omi"
if (study_name=="VAT08" & grepl("omi", COR)){
  # All COVID endpoint cases of observed non-Omicron lineages, or with unknown lineage before January 17, 2022, are excluded
  dat$EventIndPrimaryD1 = as.numeric(dat$EventIndKnownLineageOmicronOrMissingLineageD1 & dat$Omi_or_NA_after_cutoff==1) # used by cohort_event def
  dat$EventIndPrimaryD22 = as.numeric(dat$EventIndKnownLineageOmicronOrMissingLineageD22 & dat$Omi_or_NA_after_cutoff==1) 
  dat$EventIndPrimaryD43 = as.numeric(dat$EventIndKnownLineageOmicronOrMissingLineageD43 & dat$Omi_or_NA_after_cutoff==1) # used by cohort_event def
  dat$EventTimePrimaryD1 = dat$EventTimeKnownLineageOmicronOrMissingLineageD1
  dat$EventTimePrimaryD22 = dat$EventTimeKnownLineageOmicronOrMissingLineageD22 # used by scatter plot
  dat$EventTimePrimaryD43 = dat$EventTimeKnownLineageOmicronOrMissingLineageD43
}

# tpeak for this study is not set upstream for some reason, so set it here
if (study_name=="IARCHPV") {tpeak=18}

## label the subjects according to their case-control status
if (study_name=="IARCHPV"){
  dat = dat %>%
    filter(enrolltype!="Cohort" & !!as.name(config.cor$ph2) == 1) %>%
    mutate(cohort_event = factor(ifelse(enrolltype=="Case", "Any HPV Cases", ifelse(enrolltype=="Control", "Controls", NA)),
      levels = c("Any HPV Cases", "Controls")))
} else if(study_name=="ENSEMBLE" | (study_name=="PREVENT19" & COR %in% c("D35","D35prevent19nvx")) | study_name=="UK302")   { # one timepoint stage 1 studies: ENSEMBLE, PREVENT19 301 stage 1, UK302, PREVENT 301 NVX
  
  #intcur2 <- paste0("Day 15-", 28+tpeaklag, " Cases")
  dat = dat %>%
    mutate(cohort_event = factor(
      #ifelse(Perprotocol==1 & Bserostatus==0 & TwophasesampIndD29==1 & (!!as.name(paste0("EventIndPrimary", incNotMol, "D1")))==1  & (!!as.name(paste0("EventTimePrimary", incNotMol, "D1"))) <= 13, "Day 2-14 Cases",
      #       ifelse(Perprotocol==1 & Bserostatus==0 & TwophasesampIndD29==1 & (!!as.name(paste0("EventIndPrimary", incNotMol, "D1")))==1  & (!!as.name(paste0("EventTimePrimary", incNotMol, "D1"))) > 13 & (!!as.name(paste0("EventTimePrimary", incNotMol, "D1"))) <= tpeaklag-1 + NumberdaysD1toD29, intcur2,
      case_when(!!as.name(config.cor$ph2)==1 &
                  !!as.name(config.cor$EventIndPrimary)==1 ~ "Post-Peak Cases",
                !!as.name(config.cor$ph2)==1 & 
                  !!as.name(gsub(tpeak, "1", config.cor$EventIndPrimary))==0 & # remove HasVL because the there is no D1 endpoint indicator with "HasVL" in the variable name
                  AnyinfectionD1==0 ~ "Non-Cases"),
      levels = c(#"Day 2-14 Cases", intcur2, 
        "Post-Peak Cases", "Non-Cases"))
    )
} else if (study_name=="PREVENT19" & grepl("stage2", COR)) {# one timepoint stage 2 studies: PREVENT19 301 stage 2
  
  dat = dat %>%
    mutate(cohort_event = factor(
      case_when(!!as.name(config.cor$ph2)==1 &
                  !!as.name(config.cor$EventIndPrimary)==1 ~ paste0(config.cor$txt.endpoint, " Cases"),
                !!as.name(config.cor$ph2)==1 & 
                  AnyInfectionD1to22Mar26==0 ~ "Non-Cases"),
      levels = c(paste0(config.cor$txt.endpoint, " Cases"), "Non-Cases"))
    )
} else if (study_name=="AZD1222" & grepl("stage2", COR)) {# one timepoint stage 2 studies: AZD1222 stage2
  
  dat = dat %>%
    mutate(cohort_event = factor(
      case_when(!!as.name(config.cor$ph2)==1 &
                  !!as.name(config.cor$EventIndPrimary)==1 ~ paste0(config.cor$txt.endpoint, " Cases"),
                !!as.name(config.cor$ph2)==1 & 
                  AnyInfectionD1toD360==0 ~ "Non-Cases"),
      levels = c(paste0(config.cor$txt.endpoint, " Cases"), "Non-Cases"))
    )
}else if (study_name=="COVE" | study_name=="MockCOVE") {
  # for COVE, can't use ph2.tinterm=1 for now because non-case requires EarlyendpointD57==0 instead of EarlyendpointD29, may replace it with AnyinfectionD1 later
  
  dat <- dat %>%
    mutate(cohort_event = factor(
      case_when(ph2.intercurrent.cases==1 ~ "Intercurrent Cases",
                Perprotocol==1 & (!!as.name(paste0("EarlyendpointD", tpeak)))==0 & 
                  (!!as.name(paste0("TwophasesampIndD", tinterm)))==1 & 
                  (!!as.name(paste0("EventIndPrimaryD", tpeak)))==1 ~ "Post-Peak Cases", 
                    # definition for post-peak cases include people with and without D57 marker data for downstream plotting
                    # will filter out those without D57 marker data in the D57 panels
                Perprotocol==1 & 
                  (!!as.name(paste0("EarlyendpointD", tpeak)))==0 & 
                  (!!as.name(paste0("TwophasesampIndD", tpeak)))==1 & 
                  EventIndPrimaryD1==0 ~ "Non-Cases"),
      levels = c("Intercurrent Cases", "Post-Peak Cases", "Non-Cases"))
      )
} else if (study_name=="AZD1222" & !grepl("stage2", COR)){ # for two timepoints studies requiring D29 marker for D29 set, and D57 for D57 set, such as AZ
  # for AZ stage 1, can't use ph2.tinterm=1 for now because non-case requires EarlyendpointD57==0 instead of EarlyendpointD29
  
  dat <- dat %>%
    mutate(cohort_event = factor(
      case_when(ph2.intercurrent.cases==1 ~ "Intercurrent Cases",
                Perprotocol==1 & (!!as.name(paste0("EarlyendpointD", tpeak)))==0 & 
                  (!!as.name(paste0("TwophasesampIndD", tinterm)))==1 & 
                  (!!as.name(paste0("EventIndPrimaryD", tpeak)))==1 ~ "Post-Peak Cases", 
                # definition for post-peak cases include people with and without D57 marker data for downstream plotting
                # will filter out those without D57 marker data in the D57 panels
                Perprotocol==1 & 
                  #AnyinfectionD1==0 & # use EarlyendpointD29/57==0 per discussion on 8/5/2022
                  (!!as.name(paste0("EarlyendpointD", tinterm)))==0 &
                  # will filter out those with EarlyendpointD57==1 in the D57 panels
                  (!!as.name(paste0("TwophasesampIndD", tinterm)))==1 & 
                  # definition for non-cases include people with and without D57 marker data for downstream plotting
                  # will filter out those without D57 marker data in the D57 panels
                  EventIndPrimaryD1==0 ~ "Non-Cases"),
      levels = c("Intercurrent Cases", "Post-Peak Cases", "Non-Cases"))
    )
} else if (study_name=="VAT08"){
  
  dat <- dat %>%
    mutate(cohort_event = factor(
      case_when(Perprotocol==1 & (!!as.name(paste0("EarlyinfectionD", tinterm)))==0 & 
                  (ph2.D22.bAb==1 | ph2.D22.nAb==1 | ph2.D43.nAb | ph2.D43.bAb) & 
                  (!!as.name(paste0("EventIndPrimaryD", tinterm)))==1 & 
                  (!!as.name(paste0("EventTimePrimaryD", tinterm))) >= 7 &
                  (!!as.name(paste0("EventTimePrimaryD", tpeak))) < 7 ~ "7-27 days PD2 cases",
                Perprotocol==1 & (!!as.name(paste0("EarlyinfectionD", tpeak)))==0 & 
                  (ph2.D22.bAb==1 | ph2.D22.nAb==1 | ph2.D43.nAb | ph2.D43.bAb) & 
                  (!!as.name(paste0("EventIndPrimaryD", tpeak)))==1 &
                  (!!as.name(paste0("EventTimePrimaryD", tpeak))) >= 7 &
                  (!!as.name(paste0("EventTimePrimaryD", tpeak))) <= 180 ~ "28-180 days PD2 cases", 
                #Perprotocol==1 & (!!as.name(paste0("EarlyinfectionD", tpeak)))==0 & 
                #  (ph2.D22.bAb==1 | ph2.D22.nAb==1 | ph2.D43.nAb | ph2.D43.bAb) & 
                #  (!!as.name(paste0("EventIndPrimaryD", tpeak)))==1 &
                #  (!!as.name(paste0("EventTimePrimaryD", tpeak))) >= 181 &
                #  (!!as.name(paste0("EventTimePrimaryD", tpeak))) <= 365 ~ "181-365 days PD2 cases", 
                # definition for post-peak cases include people with and without D57 marker data for downstream plotting
                # will filter out those without D57 marker data in the D57 panels
                Perprotocol==1 & (!!as.name(paste0("EarlyinfectionD", tpeak)))==0 & 
                  (ph2.D22.bAb==1 | ph2.D22.nAb==1 | ph2.D43.nAb | ph2.D43.bAb) & 
                  EventIndPrimaryD1==0 ~ "Non-Cases"),
      levels = c("7-27 days PD2 cases", "28-180 days PD2 cases", #"181-365 days PD2 cases", 
                 "Non-Cases"))
    )
  
  } else {# keep other two timepoint studies except for Moderna, AZ and Sanofi here
  
  dat <- dat %>%
    mutate(cohort_event = factor(
      case_when(ph2.intercurrent.cases==1 ~ "Intercurrent Cases",
                Perprotocol==1 & (!!as.name(paste0("EarlyendpointD", tpeak)))==0 & 
                  (!!as.name(paste0("TwophasesampIndD", tinterm)))==1 & 
                  (!!as.name(paste0("EventIndPrimaryD", tpeak)))==1 ~ "Post-Peak Cases", 
                # definition for post-peak cases include people with and without D57 marker data for downstream plotting
                # will filter out those without D57 marker data in the D57 panels
                Perprotocol==1 & 
                  AnyinfectionD1==0 & 
                  (!!as.name(paste0("TwophasesampIndD", tpeak)))==1 & 
                  EventIndPrimaryD1==0 ~ "Non-Cases"),
      levels = c("Intercurrent Cases", "Post-Peak Cases", "Non-Cases"))
    )
  
}

if(study_name=="ENSEMBLE" & COR=="D29VLvariant"){
  
  # filter to baseline negative, vaccine ppt and
  # split "peak-cases" by variant type
  dat <- dat %>%
    filter(Trt==1 & Bserostatus==0) %>%
    mutate(cohort_event2 = case_when(as.character(cohort_event)=="Post-Peak Cases" & EventIndPrimaryIncludeNotMolecConfirmedD1_Anc==1 ~ "Post-Peak Cases-Reference",
                                     as.character(cohort_event)=="Post-Peak Cases" & EventIndPrimaryIncludeNotMolecConfirmedD1_Delta==1 ~ "Post-Peak Cases-Delta",
                                     as.character(cohort_event)=="Post-Peak Cases" & EventIndPrimaryIncludeNotMolecConfirmedD1_Beta==1 ~ "Post-Peak Cases-Beta",
                                     as.character(cohort_event)=="Post-Peak Cases" & EventIndPrimaryIncludeNotMolecConfirmedD1_Zeta==1 ~ "Post-Peak Cases-Zeta",
                                     as.character(cohort_event)=="Post-Peak Cases" & EventIndPrimaryIncludeNotMolecConfirmedD1_Mu==1 ~ "Post-Peak Cases-Mu",
                                     as.character(cohort_event)=="Post-Peak Cases" & EventIndPrimaryIncludeNotMolecConfirmedD1_Gamma==1 ~ "Post-Peak Cases-Gamma",
                                     as.character(cohort_event)=="Post-Peak Cases" & EventIndPrimaryIncludeNotMolecConfirmedD1_Lambda==1 ~ "Post-Peak Cases-Lambda",
                                     TRUE ~ as.character(cohort_event)))
  
}

dat <- dat[!is.na(dat$cohort_event),]

if (length(timepoints)==1) {
  ph2.indicator = config.cor$ph2
} else if (study_name != "VAT08"){
  ph2.indicator = paste0("ph2.D", tpeak) # for example: no config.cor$ph2 when COR=D29D57
}

if (study_name == "VAT08"){
  dat.cor.subset <- dat
} else {
  dat.cor.subset <- dat %>%
    dplyr::filter(!!as.name(ph2.indicator)==1)
}

write.csv(dat.cor.subset, file = here::here("data_clean", "cor_data.csv"), row.names=F)
saveRDS(dat.cor.subset, file = here::here("data_clean", "cor_data.rds"))

## arrange the dataset in the long form, expand by assay types
## dat.long.subject_level is the subject level covariates;
## dat.long.assay_value is the long-form time variables that differ by the assay type
dat.long.subject_level <- dat %>%
  replicate(length(assays),., simplify = FALSE) %>%
  bind_rows()

dat.long.assay_value.names <- times_
dat.long.assay_value <- as.data.frame(matrix(
  nrow = nrow(dat) * length(assays),
  ncol = length(dat.long.assay_value.names)
))
colnames(dat.long.assay_value) <- dat.long.assay_value.names

for (tt in seq_along(times_)) {
  dat_mock_col_names <- paste(times_[tt], assays, sep = "")
  dat.long.assay_value[, dat.long.assay_value.names[tt]] <- unlist(lapply(
    # B, Day29, Delta29overB
    dat_mock_col_names,
    # BbindSpike, BbindRBD
    function(nn) {
      if (nn %in% colnames(dat)) {
        dat[, nn]
      } else {
        rep(NA, nrow(dat))
      }
    }
  ))
}

dat.long.assay_value$assay <- rep(assays, each = nrow(dat))

dat.long <- cbind(dat.long.subject_level, dat.long.assay_value)


## change the labels of the factors for plot labels
if (study_name=="IARCHPV") {dat.long$Trt <- factor(dat.long$Trt, levels = c(1, 2, 3, 4), labels = trt.labels)
} else {dat.long$Trt <- factor(dat.long$Trt, levels = c(0, 1), labels = trt.labels)}

# no baseline serostatus for the IARCHPV study, set all Bserostatus to 0
if (study_name=="IARCHPV") {
  dat.long$Bserostatus=0
  dat.long$Bserostatus <- factor(dat.long$Bserostatus,
                                 levels = c(0),
                                 labels = bstatus.labels
  )
} else{
  dat.long$Bserostatus <- factor(dat.long$Bserostatus,
  levels = c(0, 1),
  labels = bstatus.labels
)
}


# add Hispanic or Latino vs. Not Hispanic or Latino variable
if ("EthnicityHispanic" %in% colnames(dat.long)){ # IARCHPV doesn't have this variable
  dat.long$Dich_RaceEthnic = with(dat.long,
                                ifelse(EthnicityHispanic==1, "Hispanic or Latino",
                                       ifelse(EthnicityHispanic==0 & EthnicityNotreported==0 & EthnicityUnknown==0, "Not Hispanic or Latino", NA)))
}

# add label = LLoD / poscutoff, uloq values to show in the plot
dat.long$LLoD = with(dat.long, log10(lods[assay]))
dat.long$LLoQ = with(dat.long, log10(lloqs[assay]))
dat.long$pos.cutoffs = with(dat.long, log10(pos.cutoffs[assay]))

dat.long$assay <- factor(dat.long$assay, levels = assays, labels = assays)

if ((study_name=="ENSEMBLE" | study_name=="MockENSEMBLE") & COR!="D29VLvariant"){ # for ENSEMBLE, ID50 uses LLOQ, ADCP uses LLOD
  dat.long$lb = with(dat.long, ifelse(grepl("bind", assay), "Pos.Cut", ifelse(assay=="ADCP", "LoD", "LoQ"))) 
  dat.long$lbval =  with(dat.long, ifelse(grepl("bind", assay), pos.cutoffs, ifelse(assay=="ADCP", LLoD, LLoQ))) 
} else if ((study_name=="ENSEMBLE" | study_name=="MockENSEMBLE") & COR=="D29VLvariant"){
  dat.long$lb = with(dat.long, ifelse(assay=="bindSpike", "Pos.Cut", ifelse(grepl("bind", assay), "LoQ", "LoD"))) 
  dat.long$lbval =  with(dat.long, ifelse(assay=="bindSpike", pos.cutoffs, ifelse(grepl("bind", assay), LLoQ, LLoD)))
} else if (study_name=="IARCHPV" | study_name=="UK302"){
  dat.long$lb = with(dat.long, "Pos.Cut")
  dat.long$lbval =  with(dat.long, pos.cutoffs)
} else if ((study_name %in% c("PREVENT19","AZD1222") & grepl("stage2", COR)) | study_name == "VAT08"){
  dat.long$lb = with(dat.long, ifelse(grepl("bind", assay), "LoQ", "LoD"))
  dat.long$lbval =  with(dat.long, ifelse(grepl("bind", assay), LLoQ, LLoD))
} else { # e.g. prevent19nvx
  dat.long$lb = with(dat.long, ifelse(grepl("bind", assay), "Pos.Cut", "LoD"))
  dat.long$lbval =  with(dat.long, ifelse(grepl("bind", assay), pos.cutoffs, LLoD))
}

if (study_name=="IARCHPV") {
  dat.long$ULoQ = with(dat.long, log10(uloqs[as.character(assay)]))
  
  # plot lloq instead of uloq for IARCHPV
  dat.long$lb2 = "LoQ"
  dat.long$lbval2 =  dat.long$LLoQ
} else {
  dat.long$ULoQ = with(dat.long, log10(uloqs[as.character(assay)]))
  
  #dat.long$lb2 = with(dat.long, ifelse(grepl("bind", assay) | !study_name %in% c("COVE","MockCOVE","ENSEMBLE","MockENSEMBLE"), "ULoQ", ""))
  #dat.long$lbval2 =  with(dat.long, ifelse(grepl("bind", assay) | !study_name %in% c("COVE","MockCOVE","ENSEMBLE","MockENSEMBLE"), ULoQ, -99))
  dat.long$lb2 = "ULoQ"
  dat.long$lbval2 =  with(dat.long, ifelse(!is.na(ULoQ) & ULoQ!="Inf", ULoQ, -99))
}

# assign values above the uloq to the uloq
for (t in times_[!grepl("Delta", times_)]) {
  dat.long[[t]] <- ifelse(dat.long[[t]] > dat.long$ULoQ, dat.long$ULoQ, dat.long[[t]])
}

# reset Delta29overB & Delta57overB for response call later using LLoD & ULoQ truncated data at Day 1, Day 29, Day 57
if (study_name!="IARCHPV") { # IARCHPV doesn't have delta assay variables
  for (t in timepoints) { # unique(gsub("Day", "", times_[!grepl("Delta|B", times_)]))
    tp = paste0("Day", t) # ifelse(grepl("[A-Za-z]", t), t, paste0("Day", t))
    if (!grepl("stage2", COR) & !study_name %in% c("UK302")) {dat.long[, "Delta"%.%t%.%"overB"] = dat.long[, tp] - dat.long[, "B"]
    }
  }
}

# age threshold
if (study_name=="ENSEMBLE" | study_name=="MockENSEMBLE") {age_thres=60; younger_age="Age 18 - 59"; older_age="Age >= 60"
} else {age_thres=65; younger_age="Age < 65"; older_age="Age >= 65"}
dat.long$age.geq.65 = as.integer(dat.long$Age >= age_thres)

# labels of the demographic strata for the subgroup plotting
dat.long$age_geq_65_label <-
  with(
    dat.long,
    factor(age.geq.65,
           levels = c(0, 1),
           labels = c(younger_age, older_age)
    )
  )

if (!study_name %in% c("IARCHPV","UK302")) { # IARCHPV doesn't have the high risk, sex, ethnicity variable, UK302 doesn't need these
  dat.long$highrisk_label <-
    with(
      dat.long,
      factor(HighRiskInd,
             levels = c(0, 1),
             labels = c("Not at risk", "At risk")
      )
    )


dat.long$age_risk_label <-
  with(
    dat.long,
    factor(paste0(age.geq.65, HighRiskInd),
           levels = c("00", "01", "10", "11"),
           labels = c(
             paste(younger_age, "not at risk"),
             paste(younger_age, "at risk"),
             paste(older_age, "not at risk"),
             paste(older_age, "at risk")
           )
    )
  )


dat.long$sex_label <-
  with(
    dat.long,
    factor(Sex,
           levels = c(1, 0),
           labels = c("Female", "Male")
    )
  )

  dat.long$ethnicity_label <-
    with(
      dat.long,
      ifelse(EthnicityHispanic == 1,
             "Hispanic or Latino",
             ifelse(
               EthnicityNotreported == 0 & EthnicityUnknown == 0,
               "Not Hispanic or Latino",
               "Not reported and unknown"
             ))
    ) %>% factor(
      levels = c("Hispanic or Latino", "Not Hispanic or Latino", "Not reported and unknown")
    )
}

if (study_name=="ENSEMBLE" | study_name=="MockENSEMBLE") {minor_var = "URMforsubcohortsampling"} else {minor_var = "MinorityInd"}

if (study_name!="IARCHPV") { # IARCHPV doesn't have the minority variable
  dat.long$minority_label <-
      factor(dat.long[, minor_var],
             levels = c(0, 1),
             labels = c("White Non-Hispanic", "Comm. of Color")
      )
}

# For immunogenicity characterization, complete ignore any information on cases
# vs. non-cases.  The goal is to characterize immunogenicity in the random
# subcohort, which is a stratified sample of enrolled participants. So,
# immunogenicity analysis is always done in ppts that meet all of the criteria.

# Here, only filter based on ph2.D29==1. Filtering by ph2.D57 will occur downstream,
# since it should only happen for D57-related figures.
if(length(timepoints)==1 | study_name == "VAT08"){ # one timepoint study: ph2.tpeak
  dat.long.cor.subset <- dat.long #%>%
    #dplyr::filter(!!as.name(paste0("ph2.D", tpeak, ifelse(grepl("start1", COR), "start1","")))==1)
} else {
  dat.long.cor.subset <- dat.long %>%
    dplyr::filter(!!as.name(paste0("ph2.D", tpeak, ifelse(grepl("start1", COR), "start1", ifelse(grepl("variant", COR), "variant",""))))==1)
}

if (study_name !="VAT08") {
  write.csv(dat.long.cor.subset, file = here::here("data_clean", "long_cor_data.csv"), row.names=F)
  saveRDS(dat.long.cor.subset, file = here::here("data_clean", "long_cor_data.rds"))
}

# long to longer format by time
dat.longer.cor.subset <- dat.long.cor.subset %>%
  pivot_longer(cols = all_of(times_), names_to = "time", values_to = "value")

# phase 2 filters: 
#    include both +++ and ++- at D29 for intercurrent cases and Post-Peak Cases
#    include only +++ at D57 for Post-Peak Cases
#    non-cases is defined as +++ only for Moderna, but ++-/+++ at D29/57 for AZ and Sanofi
#    for intercurrent cases at D57, Day 2-14 Cases & Day 15-35 Cases at D29, can't use ph2.D57/ph2.D29 because they are before D57/D29
if(length(timepoints) > 1 & study_name !="VAT08") {
  dat.longer.cor.subset <- dat.longer.cor.subset %>% 
    filter(!(time == paste0("Day", tpeak) & (!!as.name(paste0("ph2.D", tpeak)))==0))  # set "Day 57" in the ph2.D57 cohort  
}

if(study_name =="VAT08") {
  
  if (nrow(subset(dat.longer.cor.subset, cohort_event=="7-27 days PD2 cases"))!=0) {
    # append with a pooled case group: union of 7-27 days PD2 cases, 28-180 days PD2 cases at both day 1 and day 22 as cases for D22 marker correlates analyses
    dat.longer.cor.subset <- dat.longer.cor.subset %>%
      bind_rows(dat.longer.cor.subset %>% 
                  filter(cohort_event %in% c("7-27 days PD2 cases","28-180 days PD2 cases") & time %in% c("B", "Day22", "Delta22overB")) %>%
                  mutate(cohort_event = "7-180 days PD2 cases")) %>%
      mutate(cohort_event = factor(cohort_event, 
                                   levels = c("7-27 days PD2 cases","28-180 days PD2 cases","7-180 days PD2 cases", "Non-Cases"),
                                   labels = c("C1","C2","C3","Non-Cases")))
  } else {
    
    dat.longer.cor.subset <- dat.longer.cor.subset %>%
      mutate(cohort_event = factor(cohort_event, 
                                   levels = c("28-180 days PD2 cases","Non-Cases"),
                                   labels = c("C2","Non-Cases")))
   
  }
  
  dat.longer.cor.subset <- dat.longer.cor.subset %>% 
    filter(!(time %in% c("B", "Day22", "Delta22overB") & grepl("bind", assay) & ph2.D22.bAb == 0)) %>%
    filter(!(time %in% c("B", "Day22", "Delta22overB") & grepl("pseudoneutid50", assay) & ph2.D22.nAb == 0)) %>%
    filter(!(time %in% c("Day43", "Delta43overB", "Delta43over22") & grepl("bind", assay) & ph2.D43.bAb == 0)) %>%
    filter(!(time %in% c("Day43", "Delta43overB", "Delta43over22") & grepl("pseudoneutid50", assay) & ph2.D43.nAb == 0))
  
  # exclude EarlyinfectionD43==1 at Day 43 and D43 fold-rise over D1 for all case groups
  stopifnot(nrow(subset(dat.longer.cor.subset, cohort_event %in% c("7-27 days PD2 cases","28-180 days PD2 cases","7-180 days PD2 cases") & time %in% c("Day43", "Delta43overB", "Delta43over22") & EarlyinfectionD43==1))==0)
  
}

if (study_name=="AZD1222" & !grepl("stage2", COR)) {
  # for studies like AZ, exclude non-cases with EarlyendpointD57==1 for Day 57 panel
  # non_cases_d57 <- subset(dat.longer.cor.subset, cohort_event=="Non-Cases" & time %in% c("Day57","Delta57over29","Delta57overB"))
  # table(non_cases_d57$time, non_cases_d57$EarlyendpointD57)
  dat.longer.cor.subset <-  dat.longer.cor.subset %>%
    filter(! (cohort_event=="Non-Cases" & EarlyendpointD57==1 & time=="Day57"))
}

# define response rates
if (attr(config,"config")=="janssen_pooled_partA"){
  dat_proc$Day71pseudoneutid50uncensored=NA; dat_proc$Mon6pseudoneutid50uncensored=NA;
}

resp <- getResponder(dat_proc, cutoff.name="llox", times=grep(ifelse(study_name!="IARCHPV", "Day|Mon|^DD|^C", "M"), times_, value=T), # IARCHPV uses "M" instead of "Day" in the assay variables
               assays=assays, pos.cutoffs = pos.cutoffs)
resp_by_time_assay <- resp[, c("Ptid", colnames(resp)[grepl("Resp", colnames(resp))])] %>%
  pivot_longer(!Ptid, names_to = "category", values_to = "response")
  
dat.longer.cor.subset <- dat.longer.cor.subset %>%
  mutate(category=paste0(time, assay, "Resp")) %>%
  left_join(resp_by_time_assay, by=c("Ptid", "category")) %>%
  mutate(
    time = labels.time[time])

# define severe: severe case or non-case
# comment out on 4/7/2023 because only ENSEMBLE partA primary manuscript needs to be looped through "sev" 
#if (study_name=="ENSEMBLE" | study_name=="MockENSEMBLE") {
#  dat.longer.cor.subset <- dat.longer.cor.subset %>%
#    mutate(severe = case_when((time=="Day 1" & cohort_event != "Non-Cases" & (!!as.name(paste0("SevereEventIndPrimary", incNotMol, "D1")))==1) ~ 1,
#                              (time=="Day 29" & cohort_event != "Non-Cases" & (!!as.name(paste0("SevereEventIndPrimary", incNotMol, "D29")))==1) ~ 1,
#                              cohort_event == "Non-Cases" ~ 1,
#                              TRUE ~ 0)
#           )
#} else {dat.longer.cor.subset$severe = NA}

# only keep fold change over B for do.fold.change.overB=1: e.g. vat08
if (do.fold.change.overB==1 | study_name %in% c("VAT08")){
  dat.longer.cor.subset <- dat.longer.cor.subset %>% filter(!grepl(paste0("over D", tinterm), time))
} else if (grepl("stage2", COR)){ # keep all timepoints in times_ for stage 2 studies, such as prevent19 stage2
  dat.longer.cor.subset <- dat.longer.cor.subset
} else ( # exclude delta timepoints
  dat.longer.cor.subset <- dat.longer.cor.subset %>% filter(grepl(ifelse(study_name!="IARCHPV", "Day|Mon", "M"), time)) # IARCHPV uses "M" instead of "Day" in the assay variables
)


##### change cohort_event value for PROFISCOV because both groups in PROFISCOV are post-peak case groups
if (study_name=="PROFISCOV"){
  dat.longer.cor.subset$cohort_event = factor(with(dat.longer.cor.subset,
                                                 case_when(as.character(cohort_event)=="Intercurrent Cases" ~ "Early Post-Peak Cases",
                                                           as.character(cohort_event)=="Post-Peak Cases" ~ "Late Post-Peak Cases",
                                                           TRUE ~ as.character(cohort_event))),
                                            levels = c(if(length(timepoints)!=1)"Early Post-Peak Cases", "Late Post-Peak Cases", "Non-Cases"))
}



if ((study_name=="ENSEMBLE" | study_name=="MockENSEMBLE") & COR=="D29VLvariant"){
  # remove post-peak case rows for variant strain assay but don't belong to a variant strain (e.g. EventIndPrimaryIncludeNotMolecConfirmedD1_Beta) 
  # The only data points included in 'Zeta Cases' are from ptids with Zeta strain COVID-19
  dat.longer.cor.subset$keep_day29 = 1
  for (a in c("Anc","Delta|DeltaMDW","Beta|B.1.351","Zeta","Mu|B.1.621","Gamma|P.1","Lambda|C.37")){
    
    dat.longer.cor.subset <- dat.longer.cor.subset %>%
      mutate(assay = as.character(assay)) %>%
      mutate(assay = case_when(assay=="pseudoneutid50" ~ "pseudoneutid50_Anc", 
                               assay=="bindSpike" ~ "bindSpike_Anc",
                               TRUE ~ assay)) %>%
      mutate(keep_day29 = ifelse(grepl(a, assay) & cohort_event=="Post-Peak Cases" & 
                                   time=="Day 29" & get(paste0("EventIndPrimaryIncludeNotMolecConfirmedD1_", unlist(strsplit(a, split = "|", fixed = TRUE))[1])) == 0, 0, 
                                 # flag the peak case rows for variant assay, but is not a corresponding variant case to 0
                                 ifelse(grepl(a, assay) & cohort_event=="Non-Cases" & time=="Day 29" & is.na(value), 0, 
                                        # flag all non-case rows without the data for the corresponding variant to 0
                                        ifelse(time!="Day 29", 0,
                                               # flag all non day 29 rows to 0
                                               keep_day29)))) %>%
      mutate(assay = factor(gsub("_Anc", "", assay), levels = assays, labels = assays)) %>%
      filter(keep_day29==1)
    
    table(dat.longer.cor.subset$keep_day29, dat.longer.cor.subset$assay, dat.longer.cor.subset$cohort_event)
    dat.longer.cor.subset %>% 
      group_by(Trt, Bserostatus, cohort_event, time, assay) %>%
      summarise(n=n(), n_keep_day29 = sum(keep_day29, na.rm=T))
  }
  
  # check if rows for cases not of the variant are not kept for that assay
  for (a in c("Delta","Beta","Zeta","Mu","Gamma","Lambda")){
    # check nAb
    stopifnot(nrow(subset(dat.longer.cor.subset, cohort_event=="Post-Peak Cases" & time=="Day29" & assay==paste0("pseudoneutid50_", a) & get(paste0("EventIndPrimaryIncludeNotMolecConfirmedD1_", a)) == 0 & keep_day29 == 1)) == 0)
    
    # check bAb
    bind_var = case_when(a=="Delta" ~ "DeltaMDW",
                         a=="Beta" ~"B.1.351",
                         a=="Mu" ~ "B.1.621",
                         a=="Gamma" ~ "P.1",
                         a=="Lambda" ~ "C.37",
                         TRUE ~ "")
    if (a!="Zeta") {stopifnot(nrow(subset(dat.longer.cor.subset, cohort_event=="Post-Peak Cases" & time=="Day29" & assay==paste0("bindSpike_", bind_var) & get(paste0("EventIndPrimaryIncludeNotMolecConfirmedD1_", a)) == 0 & keep_day29 == 1)) == 0)}
  }

  stopifnot(nrow(subset(dat.longer.cor.subset, cohort_event=="Post-Peak Cases" & time=="Day29" & assay %in% c("pseudoneutid50","bindSpike") & get(paste0("EventIndPrimaryIncludeNotMolecConfirmedD1_Anc")) == 0 & keep_day29 == 1)) == 0)
}

# subsets for violin/line plots
#### figure specific data prep
# 1. define response rate:
# 2. make subsample datasets such that the violin plot only shows <= 100 non-case data points

#### for Figure 1. intercurrent vs pp, case vs non-case, (Day 1), Day 29 Day 57
groupby_vars1=c("Trt", "Bserostatus", "cohort_event", "time", "assay")

# define response rate
# for studies like IARCHPV, pooled violin plots are requested, so stack the dataset by pooling all arms thus the statistics are calculated based on the pooled arm as well
if (study_name=="IARCHPV") {
  dat.longer.cor.subset_ = dat.longer.cor.subset %>% 
    mutate(Trt="pooled") %>%
    bind_rows(dat.longer.cor.subset)
} else if (attr(config,"config")=="janssen_pooled_partA"){ # remove day 71 and mon 6 for this set of figures 1 of janssen_pooled_partA
  dat.longer.cor.subset_ = dat.longer.cor.subset %>% filter(time %in% c("Day 1","Day 29") & !is.na(value))
} else {
  dat.longer.cor.subset_ = dat.longer.cor.subset %>% filter(!is.na(value))
}

dat.longer.cor.subset.plot1 <- get_resp_by_group(dat.longer.cor.subset_, groupby_vars1)
dat.longer.cor.subset.plot1 <- dat.longer.cor.subset.plot1 %>%
  mutate(N_RespRate = ifelse(grepl("Day|M", time) && !is.na(pos.cutoffs), N_RespRate, ""),
         lb = ifelse(grepl("Day|M", time), lb, ""),
         lbval = ifelse(grepl("Day|M", time), lbval, NA),
         lb2 = ifelse(grepl("Day|M", time), lb2, ""),
         lbval2 = ifelse(grepl("Day|M", time), lbval2, NA)) # set fold-rise resp to ""
write.csv(dat.longer.cor.subset.plot1, file = here("data_clean", "longer_cor_data_plot1.csv"), row.names=F)
saveRDS(dat.longer.cor.subset.plot1, file = here("data_clean", "longer_cor_data_plot1.rds"))

# make subsample
plot.25sample1 <- get_sample_by_group(dat.longer.cor.subset.plot1, groupby_vars1)
write.csv(plot.25sample1, file = here("data_clean", "plot.25sample1.csv"), row.names=F)
saveRDS(plot.25sample1, file = here("data_clean", "plot.25sample1.rds"))

if (attr(config,"config")=="janssen_pooled_partA") { 
  # for adhoc plots showing moderate, severe and non-cases at 3 timepoints (day 29, day 71, month 6) 
  # 1 ptid with both severe and moderate cases post day 29 (moderate first), VAC31518COV3001-3010000, so moderate post-peak cases can be derived from config.cor$EventIndPrimary, i.e., EventIndPrimaryIncludeNotMolecConfirmedD29
  groupby_vars1_adhoc=c("Trt", "Bserostatus", "cohort_event_adhoc", "time", "assay")
  
  dat.longer.cor.subset.adhoc <- subset(dat.longer.cor.subset, ModerateEventIndPrimaryIncludeNotMolecConfirmedD29==1) %>%
    mutate(cohort_event_adhoc = "Moderate Post-Peak Cases") %>%
    
    bind_rows(subset(dat.longer.cor.subset, SevereEventIndPrimaryIncludeNotMolecConfirmedD29==1) %>%
                mutate(Ptid = ifelse(Ptid=="VAC31518COV3001-3010000", "VAC31518COV3001-3010000_severe", Ptid),
                       cohort_event_adhoc = "Severe Post-Peak Cases")) %>%
    
    bind_rows(subset(dat.longer.cor.subset, as.character(cohort_event) == "Non-Cases") %>%
                mutate(cohort_event_adhoc = "Non-Cases")) %>%
    mutate(cohort_event_adhoc = factor(cohort_event_adhoc, levels = c("Moderate Post-Peak Cases","Severe Post-Peak Cases","Non-Cases")))
  
  #table(unique(dat.longer.cor.subset.adhoc[,c("Ptid","cohort_event_adhoc")])$cohort_event_adhoc)
  #table(unique(dat.longer.cor.subset[,c("Ptid","cohort_event")])$cohort_event)
  dat.longer.cor.subset.adhoc = dat.longer.cor.subset.adhoc %>% filter(!is.na(value))
  
  dat.longer.cor.subset.plot1_adhoc <- get_resp_by_group(dat.longer.cor.subset.adhoc, groupby_vars1_adhoc)
  dat.longer.cor.subset.plot1_adhoc <- dat.longer.cor.subset.plot1_adhoc %>%
    mutate(N_RespRate = ifelse(grepl("Day|M", time) && !is.na(pos.cutoffs), N_RespRate, ""),
           lb = ifelse(grepl("Day|M", time), lb, ""),
           lbval = ifelse(grepl("Day|M", time), lbval, NA),
           lb2 = ifelse(grepl("Day|M", time), lb2, ""),
           lbval2 = ifelse(grepl("Day|M", time), lbval2, NA)) # set fold-rise resp to ""
  saveRDS(dat.longer.cor.subset.plot1_adhoc, file = here("data_clean", "longer_cor_data_plot1_adhoc.rds"))
}

#### for Figure 3. intercurrent vs pp, case vs non-case, (Day 1) Day 29 Day 57, by if Age >=65 and if at risk
if (!study_name %in% c("IARCHPV","UK302")) { # IARCHPV, UK302 doesn't have high risk variable
  
  if (attr(config,"config")=="janssen_pooled_partA"){ # remove day 71 and mon 6 for this set of figures 1 of janssen_pooled_partA
    dat.longer.cor.subset_ = dat.longer.cor.subset %>% filter(time %in% c("Day 1","Day 29") & !is.na(value))
  } else {
    dat.longer.cor.subset_ = dat.longer.cor.subset %>% filter(!is.na(value))
  }
  
  groupby_vars3 <- c("Trt", "Bserostatus", "cohort_event", "time", "assay", "age_geq_65_label", "highrisk_label")
  
  # define response rate
  dat.longer.cor.subset.plot3 <- get_resp_by_group(dat.longer.cor.subset_, groupby_vars3)
  dat.longer.cor.subset.plot3 <- dat.longer.cor.subset.plot3 %>%
    mutate(N_RespRate = ifelse(grepl("Day", time), N_RespRate, ""),
           lb = ifelse(grepl("Day", time), lb, ""),
           lbval = ifelse(grepl("Day", time), lbval, NA),
           lb2 = ifelse(grepl("Day", time), lb2, ""),
           lbval2 = ifelse(grepl("Day", time), lbval2, NA)) # set fold-rise resp to ""
  saveRDS(dat.longer.cor.subset.plot3, file = here("data_clean", "longer_cor_data_plot3.rds"))
  
  # make subsample
  plot.25sample3 <- get_sample_by_group(dat.longer.cor.subset.plot3, groupby_vars3)
  saveRDS(plot.25sample3, file = here("data_clean", "plot.25sample3.rds"))
}

saveRDS(as.data.frame(dat.longer.cor.subset),
        file = here("data_clean", "longer_cor_data.rds"))

if ((study_name=="ENSEMBLE" | study_name=="MockENSEMBLE") & COR=="D29VLvariant") {
  #### for Figure (variant). post peak case of specific strain vs non-case, (Day 1), Day 29 Day 57
  groupby_vars_variant=c("Trt", "Bserostatus", "cohort_event2", "time", "assay", "Region") # diff from figure 1, uses cohort_event2 instead of cohort_event, add region
  
  # define response rate
  dat.longer.cor.subset = dat.longer.cor.subset %>% filter(!is.na(value))
  
  dat.longer.cor.subset.plot.variant <- get_resp_by_group(dat.longer.cor.subset, groupby_vars_variant)
  dat.longer.cor.subset.plot.variant <- dat.longer.cor.subset.plot.variant %>%
    mutate(N_RespRate = ifelse(grepl("Day", time), N_RespRate, ""),
           lb = ifelse(grepl("Day", time), lb, ""),
           lbval = ifelse(grepl("Day", time), lbval, NA),
           lb2 = ifelse(grepl("Day", time), lb2, ""),
           lbval2 = ifelse(grepl("Day", time), lbval2, NA)) # set fold-rise resp to ""
  # unique(dat.longer.cor.subset.plot4[, c("N_RespRate","cohort_event2","assay","counts","Region")])
  write.csv(dat.longer.cor.subset.plot.variant, file = here("data_clean", "longer_cor_data_plot_variant.csv"), row.names=F)
  saveRDS(dat.longer.cor.subset.plot.variant, file = here("data_clean", "longer_cor_data_plot_variant.rds"))
}
