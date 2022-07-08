#Sys.setenv(TRIAL = "janssen_pooled_real")
#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
source(here::here("..", "_common.R"))
#-----------------------------------------------

source(here::here("code", "cor_process_function.R"))
library(here)
library(dplyr)
library(tidyverse)
library(stringr)
#dat.mock <- read.csv(here("..", "data_clean", data_name))# superceded by _common.R data read

## moved to _common.R
## COR defines the analysis to be done, e.g. D29, D57, D29start1
#Args <- commandArgs(trailingOnly=TRUE)
#if (length(Args)==0) Args=c(COR="D29D57") 
#COR=Args[1]; myprint(COR)
#
## COR has a set of analysis-specific parameters defined in the config file
config.cor <- config::get(config = COR)
#tpeak=as.integer(paste0(config.cor$tpeak))
#tpeaklag=as.integer(paste0(config.cor$tpeaklag))
#myprint(tpeak, tpeaklag)
#if (length(tpeak)==0 | length(tpeaklag)==0) stop("config "%.%COR%.%" misses some fields")

# forcing this is not a good idea. ~ Youyi
# set wt.DXX missingness to 0
wt.vars <- colnames(dat.mock)[grepl("wt.D", colnames(dat.mock))]
for (a in wt.vars) dat.mock[a][is.na(dat.mock[a])]<-0

# load parameters
source(here("code", "params.R"))


################################################
dat <- as.data.frame(dat.mock)

Args <- commandArgs(trailingOnly=TRUE)
COR=Args[1]
if (grepl("IncludeNotMolecConfirmed", COR)) {incNotMol <- "IncludeNotMolecConfirmed"
} else {incNotMol <- ""}

# set EventIndTimePrimary to EventIndTimeOmicron if study_name=="VAT08m" & COR=="D22D43omi"
if (study_name=="VAT08m" & grepl("omi", COR)){
  dat$EventIndPrimaryD1 = dat$EventIndOmicronD1 # used by cohort_event def
  dat$EventIndPrimaryD22 = dat$EventIndOmicronD22
  dat$EventIndPrimaryD43 = dat$EventIndOmicronD43 # used by cohort_event def
  dat$EventTimePrimaryD1 = dat$EventTimeOmicronD1
  dat$EventTimePrimaryD22 = dat$EventTimeOmicronD22 # used by scatter plot
  dat$EventTimePrimaryD43 = dat$EventTimeOmicronD43
}

## label the subjects according to their case-control status
## add case vs non-case indicators
if(study_name=="ENSEMBLE" | study_name=="MockENSEMBLE" | study_name=="HVTN705" | study_name=="PREVENT19")  {
  
  #intcur2 <- paste0("Day 15-", 28+tpeaklag, " Cases")
dat = dat %>%
    mutate(cohort_event = factor(
      #ifelse(Perprotocol==1 & Bserostatus==0 & TwophasesampIndD29==1 & (!!as.name(paste0("EventIndPrimary", incNotMol, "D1")))==1  & (!!as.name(paste0("EventTimePrimary", incNotMol, "D1"))) <= 13, "Day 2-14 Cases",
      #       ifelse(Perprotocol==1 & Bserostatus==0 & TwophasesampIndD29==1 & (!!as.name(paste0("EventIndPrimary", incNotMol, "D1")))==1  & (!!as.name(paste0("EventTimePrimary", incNotMol, "D1"))) > 13 & (!!as.name(paste0("EventTimePrimary", incNotMol, "D1"))) <= tpeaklag-1 + NumberdaysD1toD29, intcur2,
      case_when(!!as.name(config.cor$ph2)==1 & !!as.name(config.cor$EventIndPrimary)==1 ~ "Post-Peak Cases",
                !!as.name(config.cor$ph2)==1 & !!as.name(paste0("EventIndPrimary", incNotMol, "D1"))==0  & AnyinfectionD1==0 ~ "Non-Cases"),
      levels = c(#"Day 2-14 Cases", intcur2, 
        "Post-Peak Cases", "Non-Cases"))
    )
} else if (study_name=="COVE" | study_name=="MockCOVE") {
  # for COVE, can't use ph2.tinterm=1 for now because non-case requires EarlyendpointD57==0 instead of EarlyendpointD29, may replace it with AnyinfectionD1 later
  
  dat <- dat %>%
    mutate(cohort_event = factor(
      case_when(ph2.intercurrent.cases==1 ~ "Intercurrent Cases",
                Perprotocol==1 & (!!as.name(paste0("EarlyendpointD", tpeak)))==0 & (!!as.name(paste0("TwophasesampIndD", tinterm)))==1 & (!!as.name(paste0("EventIndPrimaryD", tpeak)))==1 ~ "Post-Peak Cases", 
                    # definition for post-peak cases include people with and without D57 marker data for downstream plotting
                    # will filter out those without D57 marker data in the D57 panels
                Perprotocol==1 & (!!as.name(paste0("EarlyendpointD", tpeak)))==0 & (!!as.name(paste0("TwophasesampIndD", tpeak)))==1 & EventIndPrimaryD1==0 ~ "Non-Cases"),
      levels = c("Intercurrent Cases", "Post-Peak Cases", "Non-Cases"))
      )
} else {
  # for AZ, can't use ph2.tinterm=1 for now because non-case requires EarlyendpointD57==0 instead of EarlyendpointD29
  
  dat <- dat %>%
    mutate(cohort_event = factor(
      case_when(ph2.intercurrent.cases==1 ~ "Intercurrent Cases",
                Perprotocol==1 & (!!as.name(paste0("EarlyendpointD", tpeak)))==0 & (!!as.name(paste0("TwophasesampIndD", tinterm)))==1 & (!!as.name(paste0("EventIndPrimaryD", tpeak)))==1 ~ "Post-Peak Cases", 
                    # definition for post-peak cases include people with and without D57 marker data for downstream plotting
                    # will filter out those without D57 marker data in the D57 panels
                Perprotocol==1 & AnyinfectionD1==0 & (!!as.name(paste0("TwophasesampIndD", tpeak)))==1 & EventIndPrimaryD1==0 ~ "Non-Cases"),
      levels = c("Intercurrent Cases", "Post-Peak Cases", "Non-Cases"))
    )
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


# add Hispanic or Latino vs. Not Hispanic or Latino variable
dat.long$Dich_RaceEthnic = with(dat.long,
                                ifelse(EthnicityHispanic==1, "Hispanic or Latino",
                                       ifelse(EthnicityHispanic==0 & EthnicityNotreported==0 & EthnicityUnknown==0, "Not Hispanic or Latino", NA)))

# add label = LLoD / poscutoff, uloq values to show in the plot
dat.long$LLoD = with(dat.long, log10(llods[as.character(assay)]))
dat.long$pos.cutoffs = with(dat.long, log10(pos.cutoffs[as.character(assay)]))
dat.long$lb = with(dat.long, ifelse(grepl("bind", assay), "Pos.Cut", "LoD")) 
dat.long$lbval =  with(dat.long, ifelse(grepl("bind", assay), pos.cutoffs, LLoD))
dat.long$ULoQ = with(dat.long, log10(uloqs[as.character(assay)]))
dat.long$lb2 = with(dat.long, ifelse(grepl("bind", assay) | !study_name %in% c("COVE","MockCOVE","ENSEMBLE","MockENSEMBLE"), "ULoQ", ""))
dat.long$lbval2 =  with(dat.long, ifelse(grepl("bind", assay) | !study_name %in% c("COVE","MockCOVE","ENSEMBLE","MockENSEMBLE"), ULoQ, -99))

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

# labels of the demographic strata for the subgroup plotting
dat.long$age_geq_65_label <-
  with(
    dat.long,
    factor(age.geq.65,
           levels = c(0, 1),
           labels = c(younger_age, older_age)
    )
  )

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

if (study_name=="ENSEMBLE" | study_name=="MockENSEMBLE") {minor_var = "URMforsubcohortsampling"} else {minor_var = "MinorityInd"}
dat.long$minority_label <-
    factor(dat.long[, minor_var],
           levels = c(0, 1),
           labels = c("White Non-Hispanic", "Comm. of Color")
    )

# For immunogenicity characterization, complete ignore any information on cases
# vs. non-cases.  The goal is to characterize immunogenicity in the random
# subcohort, which is a stratified sample of enrolled participants. So,
# immunogenicity analysis is always done in ppts that meet all of the criteria.

# Here, only filter based on ph2.D29==1. Filtering by ph2.D57 will occur downstream,
# since it should only happen for D57-related figures.
if(study_name=="ENSEMBLE" | study_name=="MockENSEMBLE" | study_name=="HVTN705" | study_name=="PREVENT19"){ # one timepoint study: ph2.tpeak
  dat.long.cor.subset <- dat.long #%>%
    #dplyr::filter(!!as.name(paste0("ph2.D", tpeak, ifelse(grepl("start1", COR), "start1","")))==1)
} else {# two timepoints study: ph2.tinterm
  dat.long.cor.subset <- dat.long %>%
    dplyr::filter(!!as.name(paste0("ph2.D", tinterm, ifelse(grepl("start1", COR), "start1","")))==1)
}


# long to longer format by time
dat.longer.cor.subset <- dat.long.cor.subset %>%
  pivot_longer(cols = all_of(times), names_to = "time", values_to = "value")

# phase 2 filters: 
#    include both +++ and ++- at D29 for intercurrent cases and Post-Peak Cases
#    include only +++ at D57 for intercurrent cases and Post-Peak Cases
#    non-cases is defined as +++
#    for intercurrent cases at D57, Day 2-14 Cases & Day 15-35 Cases at D29, can't use ph2.D57/ph2.D29 because they are before D57/D29
if(!(study_name=="ENSEMBLE" | study_name=="MockENSEMBLE" | study_name=="PREVENT19")) {
  dat.longer.cor.subset <- dat.longer.cor.subset %>% 
    filter(!(time == paste0("Day", tpeak) & (!!as.name(paste0("ph2.D", tpeak)))==0))  # set "Day 57" in the ph2.D57 cohort  
}

# define response rates
resp <- getResponder(dat.mock, cutoff.name="llox", times=grep("Day", times, value=T), 
               assays=assays, pos.cutoffs = pos.cutoffs)
resp_by_time_assay <- resp[, c("Ptid", colnames(resp)[grepl("Resp", colnames(resp))])] %>%
  pivot_longer(!Ptid, names_to = "category", values_to = "response")
  
dat.longer.cor.subset <- dat.longer.cor.subset %>%
  mutate(category=paste0(time, assay, "Resp")) %>%
  left_join(resp_by_time_assay, by=c("Ptid", "category")) %>%
  mutate(
    time = labels.time[time])

# define severe: severe case or non-case
if (study_name=="ENSEMBLE" | study_name=="MockENSEMBLE") {
  dat.longer.cor.subset <- dat.longer.cor.subset %>%
    mutate(severe = case_when((time=="Day 1" & cohort_event != "Non-Cases" & (!!as.name(paste0("SevereEventIndPrimary", incNotMol, "D1")))==1) ~ 1,
                              (time=="Day 29" & cohort_event != "Non-Cases" & (!!as.name(paste0("SevereEventIndPrimary", incNotMol, "D29")))==1) ~ 1,
                              cohort_event == "Non-Cases" ~ 1,
                              TRUE ~ 0)
           )
} else {dat.longer.cor.subset$severe = NA}

# only keep fold change for do.fold.change=1: e.g. vat08m_nonnaive
if (do.fold.change==1){
  dat.longer.cor.subset <- dat.longer.cor.subset %>% filter(!grepl(paste0("over D", tinterm), time))
} else (
  dat.longer.cor.subset <- dat.longer.cor.subset %>% filter(grepl("Day", time))
)


# subsets for violin/line plots
#### figure specific data prep
# 1. define response rate:
# 2. make subsample datasets such that the violin plot only shows <= 100 non-case data points

#### for Figure 1. intercurrent vs pp, case vs non-case, (Day 1), Day 29 Day 57
groupby_vars1=c("Trt", "Bserostatus", "cohort_event", "time", "assay")

# define response rate
dat.longer.cor.subset.plot1 <- get_resp_by_group(dat.longer.cor.subset, groupby_vars1)
dat.longer.cor.subset.plot1 <- dat.longer.cor.subset.plot1 %>%
  mutate(N_RespRate = ifelse(grepl("Day", time), N_RespRate, ""),
         lb = ifelse(grepl("Day", time), lb, ""),
         lbval = ifelse(grepl("Day", time), lbval, NA),
         lb2 = ifelse(grepl("Day", time), lb2, ""),
         lbval2 = ifelse(grepl("Day", time), lbval2, NA)) # set fold-rise resp to ""
write.csv(dat.longer.cor.subset.plot1, file = here("data_clean", "longer_cor_data_plot1.csv"), row.names=F)
saveRDS(dat.longer.cor.subset.plot1, file = here("data_clean", "longer_cor_data_plot1.rds"))

# make subsample
plot.25sample1 <- get_sample_by_group(dat.longer.cor.subset.plot1, groupby_vars1)
write.csv(plot.25sample1, file = here("data_clean", "plot.25sample1.csv"), row.names=F)
saveRDS(plot.25sample1, file = here("data_clean", "plot.25sample1.rds"))

#### for Figure 3. intercurrent vs pp, case vs non-case, (Day 1) Day 29 Day 57, by if Age >=65 and if at risk
groupby_vars3 <- c("Trt", "Bserostatus", "cohort_event", "time", "assay", "age_geq_65_label", "highrisk_label")

# define response rate
dat.longer.cor.subset.plot3 <- get_resp_by_group(dat.longer.cor.subset, groupby_vars3)
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


saveRDS(as.data.frame(dat.longer.cor.subset),
        file = here("data_clean", "longer_cor_data.rds"))
