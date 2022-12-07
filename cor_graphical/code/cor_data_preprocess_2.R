#Sys.setenv(TRIAL = "janssen_pooled_real")
#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
source(here::here("..", "_common.R"))
#-----------------------------------------------

library(here)
library(dplyr)
library(tidyverse)
library(stringr)
#dat.mock <- read.csv(here("..", "data_clean", data_name))# superceded by _common.R data read

# forcing this is not a good idea. ~ Youyi
# set wt.DXX missingness to 0
wt.vars <- colnames(dat.mock)[grepl("wt.D", colnames(dat.mock))]
for (a in wt.vars) dat.mock[a][is.na(dat.mock[a])]<-0

if(T){ # for ENSEMBLE SA and LA reports only
  # copy Bpseudoneutid50 to Bpseudoneutid50la & calculate delta value if Day29pseudoneutid50la exists and is required for reporting
  # copy Bpseudoneutid50 to Bpseudoneutid50sa & calculate delta value if Day29pseudoneutid50sa exists and is required for reporting
  if ("Day29pseudoneutid50la" %in% colnames(dat.mock) & "pseudoneutid50la" %in% assays) {
    dat.mock$Bpseudoneutid50la = dat.mock$Bpseudoneutid50
    dat.mock$Delta29overBpseudoneutid50la = pmin(log10(uloqs["pseudoneutid50la"]), dat.mock$Day29pseudoneutid50la) - pmin(log10(uloqs["pseudoneutid50la"]), dat.mock$Bpseudoneutid50la)
  }
  if ("Day29pseudoneutid50sa" %in% colnames(dat.mock) & "pseudoneutid50sa" %in% assays) {
    dat.mock$Bpseudoneutid50sa = dat.mock$Bpseudoneutid50
    dat.mock$Delta29overBpseudoneutid50sa = pmin(log10(uloqs["pseudoneutid50sa"]), dat.mock$Day29pseudoneutid50sa) - pmin(log10(uloqs["pseudoneutid50sa"]), dat.mock$Bpseudoneutid50sa)
  }
}

# load parameters
source(here("code", "params.R"))
################################################
dat <- as.data.frame(dat.mock)

Args <- commandArgs(trailingOnly=TRUE)
COR=Args[1]
if (grepl("IncludeNotMolecConfirmed", COR)) {incNotMol <- "IncludeNotMolecConfirmed"
} else {incNotMol <- ""}

# set EventIndTimePrimary to EventIndTimeOmicron if study_name=="VAT08m" & COR %in% c("D22omi","D43omi")
if (study_name=="VAT08m" & grepl("omi", COR)){
  dat$EventIndPrimaryD1 = dat$EventIndOmicronD1 # used by cohort_event def
  dat$EventIndPrimaryD22 = dat$EventIndOmicronD22
  dat$EventIndPrimaryD43 = dat$EventIndOmicronD43
  dat$EventTimePrimaryD1 = dat$EventTimeOmicronD1
  dat$EventTimePrimaryD22 = dat$EventTimeOmicronD22
  dat$EventTimePrimaryD43 = dat$EventTimeOmicronD43
}


## label the subjects according to their case-control status
## add case vs non-case indicators
if(study_name=="COVE" | study_name=="MockCOVE")  {
  # only Moderna uses TwophasesampIndD57==1 for non-cases at both D29 and D57
  
  dat = dat %>%
    mutate(cohort_event = factor(case_when(
      #Perprotocol==1 &
        #!!as.name(config.cor$Earlyendpoint)==0 & 
        #!!as.name(paste0("TwophasesampIndD", config.cor$tpeak))==1 &
        #ph2.Dxx = EarlyendpointDxx==0 & Perprotocol==1 & EventTimePrimaryDxx>=7 & TwophasesampIndDxx
      !!as.name(config.cor$ph2)==1 &
        !!as.name(config.cor$EventIndPrimary)==1 ~ "Cases",
      Perprotocol==1 & # do not use config.cor$ph2 here for now to avoid confusion (config is D29 but need D57 args for non-cases)
        !!as.name(paste0("EarlyendpointD", timepoints[length(timepoints)]))==0 &
        !!as.name(paste0("TwophasesampIndD", timepoints[length(timepoints)]))==1 & 
        EventIndPrimaryD1==0 ~ "Non-Cases"),
      levels = c(#"Day 2-14 Cases", intcur2, 
        "Cases", "Non-Cases")
    ))
  
} else if (study_name=="AZD1222") {
  # for whatever reason very different sets of ptids had D29 or D57 measured for AZ, so for D29, include a ptid if have D29 data ignoring D57 availability; for D57, include a ptid if have D57 data ignoing D29 availability
  # keep Sanofi here as well
  
  dat = dat %>%
    mutate(cohort_event = factor(case_when(
      #Perprotocol==1 &
        #!!as.name(config.cor$Earlyendpoint)==0 & 
        #!!as.name(paste0("TwophasesampIndD", config.cor$tpeak))==1 & 
        #ph2.Dxx = EarlyendpointDxx==0 & Perprotocol==1 & EventTimePrimaryDxx>=7 & TwophasesampIndDxx
      !!as.name(config.cor$ph2)==1 &
        !!as.name(config.cor$EventIndPrimary)==1 ~ "Cases",
      Perprotocol==1 & # do not use config.cor$ph2 here for now because it doesn't include AnyinfectionD1
        AnyinfectionD1==0 &
        !!as.name(paste0("TwophasesampIndD", config.cor$tpeak))==1 & 
        !!as.name(paste0("EventIndPrimary", incNotMol, "D1"))==0 ~ "Non-Cases"),
      levels = c(#"Day 2-14 Cases", intcur2, 
        "Cases", "Non-Cases")
    ))
} else {# keep Sanofi and other two timepoint studies except for AZ and Moderna here
  
  dat = dat %>%
    mutate(cohort_event = factor(case_when(
      #Perprotocol==1 &
      #!!as.name(config.cor$Earlyendpoint)==0 & 
      #!!as.name(paste0("TwophasesampIndD", config.cor$tpeak))==1 & 
      #ph2.Dxx = EarlyendpointDxx==0 & Perprotocol==1 & EventTimePrimaryDxx>=7 & TwophasesampIndDxx
      !!as.name(config.cor$ph2)==1 &
        !!as.name(config.cor$EventIndPrimary)==1 ~ "Cases",
      Perprotocol==1 & # do not use config.cor$ph2 here for now because it doesn't include AnyinfectionD1
        AnyinfectionD1==0 &
        !!as.name(paste0("TwophasesampIndD", timepoints[length(timepoints)]))==1 & 
        !!as.name(paste0("EventIndPrimary", incNotMol, "D1"))==0 ~ "Non-Cases"),
      levels = c(#"Day 2-14 Cases", intcur2, 
        "Cases", "Non-Cases")
    ))
}

dat <- dat[!is.na(dat$cohort_event),]

## arrange the dataset in the long form, expand by assay types
## dat.long.subject_level is the subject level covariates;
## dat.long.assay_value is the long-form time variables that differ by the assay type
dat.long.subject_level <- dat %>%
  replicate(length(assays),., simplify = FALSE) %>%
  bind_rows()

dat.long.assay_value.names <- times
dat.long.assay_value <- as.data.frame(matrix(
  nrow = nrow(dat) * length(assays),
  ncol = length(dat.long.assay_value.names)
))
colnames(dat.long.assay_value) <- dat.long.assay_value.names

for (tt in seq_along(times)) {
  dat_mock_col_names <- paste(times[tt], assays, sep = "")
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
dat.long$Trt <- factor(dat.long$Trt, levels = c(0, 1), labels = trt.labels)
dat.long$Bserostatus <- factor(dat.long$Bserostatus,
  levels = c(0, 1),
  labels = bstatus.labels
)
dat.long$assay <- factor(dat.long$assay, levels = assays, labels = assays)

# add LLoQ pos.cutoffs, and ULoQ value for response call and censoring - log10 scales
#dat.long$LLoQ = with(dat.long, log10(lloqs[as.character(assay)]))
#dat.long$pos.cutoffs = with(dat.long, log10(pos.cutoffs[as.character(assay)]))
dat.long$ULoQ = with(dat.long, log10(uloqs[as.character(assay)]))

# add label = LLoD / poscutoff, uloq values to show in the plot
#dat.long$LLoD = with(dat.long, log10(llods[as.character(assay)]))
#dat.long$lb = with(dat.long, ifelse(grepl("bind", assay), "Pos.Cut", "LoD")) 
#dat.long$lbval =  with(dat.long, ifelse(grepl("bind", assay), pos.cutoffs, LLoD))
#dat.long$lb2 = with(dat.long, ifelse(grepl("bind", assay), "ULoQ", "")) 
#dat.long$lbval2 =  with(dat.long, ifelse(grepl("bind", assay), ULoQ, -99))

# assign values above the uloq to the uloq
for (t in times[!grepl("Delta", times)]) {
  dat.long[[t]] <- ifelse(dat.long[[t]] > dat.long$ULoQ, dat.long$ULoQ, dat.long[[t]])
}

# reset Delta29overB & Delta57overB for response call later using LLoD & ULoQ truncated data at Day 1, Day 29, Day 57
for (t in unique(gsub("Day", "", times[!grepl("Delta|B", times)]))) {
  dat.long[, "Delta"%.%t%.%"overB"] = dat.long[, "Day"%.%t] - dat.long[, "B"]
}

# age threshold
if (study_name=="ENSEMBLE" | study_name=="MockENSEMBLE") {age_thres=60; younger_age="Age 18 - 59"; older_age="Age >= 60"
} else {age_thres=65; younger_age="Age < 65"; older_age="Age >= 65"}
dat.long$age.geq.65 = as.integer(dat.long$Age >= age_thres)

# # matrix to decide the sampling strata
dat.long$demo_lab <-
  with(dat.long, factor(paste0(age.geq.65, HighRiskInd),
    levels = c("00", "01", "10", "11"),
    labels = c(
      paste(younger_age, "not at risk"),
      paste(younger_age, "at risk"),
      paste(older_age, "not at risk"),
      paste(older_age, "at risk")
    )
  ))

# labels of the demographic strata for the subgroup plotting
#dat.long$trt_bstatus_label <-
#  with(
#    dat.long,
#    factor(paste0(as.numeric(Trt), as.numeric(Bserostatus)),
#      levels = c("11", "12", "21", "22"),
#      labels = c(
#        "Placebo, Baseline Neg",
#        "Placebo, Baseline Pos",
#        "Vaccine, Baseline Neg",
#        "Vaccine, Baseline Pos"
#      )
#    )
#  )

# For immunogenicity characterization, complete ignore any information on cases
# vs. non-cases.  The goal is to characterize immunogenicity in the random
# subcohort, which is a stratified sample of enrolled participants. So,
# immunogenicity analysis is always done in ppts that meet all of the criteria.

dat.cor.subset <- dat %>%
  dplyr::filter(!!as.name(paste0("ph2.D", tpeak, ifelse(grepl("start1", COR), "start1","")))==1)
dat.long.cor.subset <- dat.long %>%
  dplyr::filter(!!as.name(paste0("ph2.D", tpeak, ifelse(grepl("start1", COR), "start1","")))==1)


saveRDS(as.data.frame(dat.long.cor.subset),
        file = here("data_clean", "long_cor_data.rds")
)
saveRDS(as.data.frame(dat.cor.subset),
        file = here("data_clean", "cor_data.rds")
)

