#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
source(here::here("..", "_common.R"))
source(here::here("code", "params.R")) # load parameters
#-----------------------------------------------

library(here)
library(dplyr)
library(stringr)
if (F){
  # adhoc for AZ, pair plot with spike and pseudovirus side by side
  # 1. add azd1222_all with both assays in config.yml, Sys.setenv(TRIAL="azd1222_all")
  # azd1222_all: &azd1222_all
  # <<: *azd1222_base
  # assays: [bindSpike, pseudoneutid50]
  # llox_label: [LLOQ, LOD]
  # assay_labels: [Binding Antibody to Spike, PsV Neutralization 50% Titer]
  # assay_labels_short: [Anti Spike IgG (BAU/ml), Pseudovirus-nAb ID50 (IU50/ml)]
  # 2. create azd1222_all_data_processed_with_riskscore 
  # by combining azd1222_data_processed_with_riskscore.csv and azd1222_bAb_data_processed_with_riskscore.csv and assign to dat.mock
  azd1222_bAb <- read.csv(here("..", "data_clean", "azd1222_bAb_data_processed_with_riskscore.csv"), header = TRUE)
  azd1222 <- read.csv(here("..", "data_clean", "azd1222_data_processed_with_riskscore.csv"), header = TRUE)
  azd1222_bAb$Bpseudoneutid50=NULL
  azd1222_bAb$Day29pseudoneutid50=NULL
  azd1222_bAb$Day57pseudoneutid50=NULL
  azd1222_bAb$wt.subcohort=NULL
  dat.mock <- azd1222_bAb %>%
    left_join(azd1222[,c("Ptid","Bpseudoneutid50","Day29pseudoneutid50","Day57pseudoneutid50",
                         "Delta29overBpseudoneutid50","Delta57overBpseudoneutid50","Delta57over29pseudoneutid50",
                         "ph2.immuno","wt.subcohort","TwophasesampIndD29","TwophasesampIndD57")], by="Ptid")
  # wt.subcohort is from nAb dataset based on email discussion with Youyi on 7/22/2022:
  # Youyi: one based on ID50 weights because we have less ID50 samples than bAb samples
  table(dat.mock$ph2.immuno.x, dat.mock$ph2.immuno.y)
  dat.mock$ph2.immuno = with(dat.mock, ph2.immuno.x==1 & ph2.immuno.y==1, 1, 0) # 628
  dat.mock$TwophasesampIndD29 = with(dat.mock, TwophasesampIndD29.x==1 & TwophasesampIndD29.y==1, 1, 0) # 828
  dat.mock$TwophasesampIndD57 = with(dat.mock, TwophasesampIndD57.x==1 & TwophasesampIndD57.y==1, 1, 0) # 659
  
  dim(subset(dat.mock, EarlyendpointD57==0 & Perprotocol==1 & SubcohortInd==1 & 
               !is.na(Bpseudoneutid50) & !is.na(Day29pseudoneutid50) & 
               !is.na(BbindSpike) & !is.na(Day29bindSpike))) # 773, 
  # if change to EarlyendpointD29==0, save three participants by using EarlyendpointD29 instead of EarlyendpointD57 for Day 29 plots
  dim(subset(dat.mock, EarlyendpointD57==0 & Perprotocol==1 & SubcohortInd==1 & 
               !is.na(Bpseudoneutid50) & !is.na(Day29pseudoneutid50) & !is.na(Day57pseudoneutid50) & 
               !is.na(BbindSpike) & !is.na(Day29bindSpike) & !is.na(Day57bindSpike))) # 628

  # subsetting on vaccine recipients with ID50 value > LOD and with IgG spike > positivity cut-off at Day 57
  dat.mock <- subset(dat.mock, Day57bindSpike > log10(pos.cutoffs["bindSpike"]) & Day57pseudoneutid50 > log10(llods["pseudoneutid50"]))
}
  
dat.mock <- read.csv(data_name, header = TRUE)

# for unknown reason, the Senior variable has no value in the PROFISCOV dataset
if (study_name=="PROFISCOV") {
  dat.mock$Senior <- as.numeric(with(dat.mock, Age >= age.cutoff, 1, 0))
  
  dat.mock <- dat.mock %>%
    mutate(
      race = labels.race[1],
      race = case_when(
        Black == 1 ~ labels.race[2],
        Asian == 1 ~ labels.race[3],
        NatAmer == 1 ~ labels.race[4],
        PacIsl == 1 ~ labels.race[5],
        Multiracial == 1 ~ labels.race[6],
        Notreported == 1 | Unknown == 1 ~ labels.race[7],
        TRUE ~ labels.race[1]
      ),
      race = factor(race, levels = labels.race)
    )
}

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

dat <- dat.mock

print("Data preprocess")

# For immunogenicity characterization, complete ignore any information on cases
# vs. non-cases.  The goal is to characterize immunogenicity in the random
# subcohort, which is a stratified sample of enrolled participants. So,
# immunogenicity analysis is always done in ppts that meet all of the criteria.
dat.twophase.sample <- dat %>%
  filter(ph2.immuno == 1)
twophase_sample_id <- dat.twophase.sample$Ptid

important.columns <- c("Ptid", "Trt", "MinorityInd", "HighRiskInd", "Age", "Sex",
  "Bserostatus", "Senior", "Bstratum", "wt.subcohort", 
  "race","EthnicityHispanic","EthnicityNotreported", 
  "EthnicityUnknown", "WhiteNonHispanic", if (study_name !="COVE" & study_name!="MockCove") "HIVinfection", 
  if (study_name !="COVE" & study_name !="MockCove" & study_name !="PROFISCOV") "Country")

## arrange the dataset in the long form, expand by assay types
## dat.long.subject_level is the subject level covariates;
## dat.long.assay_value is the long-form time variables that differ by the assay type
dat.long.subject_level <- dat[, important.columns] %>%
  replicate(length(assay_immuno), ., simplify = FALSE) %>%
  bind_rows()


dat.long.assay_value.names <- times
dat.long.assay_value <- as.data.frame(matrix(
  nrow = nrow(dat) * length(assay_immuno),
  ncol = length(dat.long.assay_value.names)
))
colnames(dat.long.assay_value) <- dat.long.assay_value.names

for (tt in seq_along(times)) {
  dat_mock_col_names <- paste(times[tt], assay_immuno, sep = "")
  dat.long.assay_value[, dat.long.assay_value.names[tt]] <- unlist(lapply(
    dat_mock_col_names,
    function(nn) {
      if (nn %in% colnames(dat)) {
        dat[, nn]
      } else {
        rep(NA, nrow(dat))
      }
    }
  ))
}

dat.long.assay_value$assay <- rep(assay_immuno, each = nrow(dat))

dat.long <- cbind(dat.long.subject_level, dat.long.assay_value)


## change the labels of the factors for plot labels
dat.long$Trt <- factor(dat.long$Trt, levels = c(0, 1), labels = trt.labels)
dat.long$Bserostatus <- factor(dat.long$Bserostatus,
  levels = c(0, 1),
  labels = bstatus.labels
)
dat.long$assay <- factor(dat.long$assay, levels = assay_immuno, labels = assay_immuno)

dat.long.twophase.sample <- dat.long[dat.long$Ptid %in% twophase_sample_id, ]
dat.twophase.sample <- subset(dat, Ptid %in% twophase_sample_id)


# labels of the demographic strata for the subgroup plotting
dat.long.twophase.sample$trt_bstatus_label <-
  with(
    dat.long.twophase.sample,
    factor(paste0(as.numeric(Trt), as.numeric(Bserostatus)),
      levels = c("11", "12", "21", "22"),
      labels = c(
        "Placebo, Baseline Neg",
        "Placebo, Baseline Pos",
        "Vaccine, Baseline Neg",
        "Vaccine, Baseline Pos"
      )
    )
  )

dat.long.twophase.sample$age_geq_65_label <-
  with(
    dat.long.twophase.sample,
    factor(Senior,
      levels = c(0, 1),
      labels = paste0(c("Age < ", "Age >= "), age.cutoff)
    )
  )

dat.long.twophase.sample$highrisk_label <-
  with(
    dat.long.twophase.sample,
    factor(HighRiskInd,
      levels = c(0, 1),
      labels = c("Not at risk", "At risk")
    )
  )

dat.long.twophase.sample$age_risk_label <-
  with(
    dat.long.twophase.sample,
    factor(paste0(Senior, HighRiskInd),
      levels = c("00", "01", "10", "11"),
      labels = c(
        paste0("Age < ", age.cutoff, " not at risk"),
        paste0("Age < ", age.cutoff, " at risk"),
        paste0("Age >= ", age.cutoff, " not at risk"),
        paste0("Age >= ", age.cutoff, " at risk")
      )
    )
  )

if (study_name!="ENSEMBLE" & study_name!="MockENSEMBLE") {

  dat.long.twophase.sample$sex_label <-
    with(
      dat.long.twophase.sample,
      factor(Sex,
        levels = c(1, 0),
        labels = c("Female", "Male")
      )
    )
  
  dat.long.twophase.sample$age_sex_label <-
    with(
      dat.long.twophase.sample,
      factor(paste0(Senior, Sex),
        levels = c("00", "01", "10", "11"),
        labels = c(
          paste0("Age < ", age.cutoff, " male"),
          paste0("Age < ", age.cutoff, " female"),
          paste0("Age >= ", age.cutoff, " male"),
          paste0("Age >= ", age.cutoff, " female")
        )
      )
    )

} else if (study_name=="ENSEMBLE" | study_name=="MockENSEMBLE") {
  
  dat.long.twophase.sample$sex_label <-
    with(
      dat.long.twophase.sample,
      factor(Sex,
             levels = c(0, 1, 2, 3),
             labels = c("Male", "Female", "Undifferentiated", "Unknown")
      )
    )
  
  dat.long.twophase.sample$age_sex_label <-
    with(
      dat.long.twophase.sample,
      factor(paste0(Senior, Sex),
             levels = c("00", "01", "02", "03", "10", "11", "12", "13"),
             labels = c(
               paste0("Age < ", age.cutoff, " male"),
               paste0("Age < ", age.cutoff, " female"),
               paste0("Age < ", age.cutoff, " undifferentiated"),
               paste0("Age < ", age.cutoff, " unknown"),
               paste0("Age >= ", age.cutoff, " male"),
               paste0("Age >= ", age.cutoff, " female"),
               paste0("Age >= ", age.cutoff, " undifferentiated"),
               paste0("Age >= ", age.cutoff, " unknown")
             )
      )
    )
  
  # Ignore undifferentiated participants
  dat.long.twophase.sample$sex_label[dat.long.twophase.sample$sex_label == "Undifferentiated"] <- NA
  
  dat.long.twophase.sample$age_sex_label[endsWith(as.character(dat.long.twophase.sample$age_sex_label), "undifferentiated")] <- NA
  
  # Remove factor levels that aren't present in the data
  dat.long.twophase.sample$sex_label <- droplevels.factor(dat.long.twophase.sample$sex_label)
  
  dat.long.twophase.sample$age_sex_label <- droplevels.factor(dat.long.twophase.sample$age_sex_label)
  
}

dat.long.twophase.sample$ethnicity_label <-
  with(
    dat.long.twophase.sample,
    ifelse(
      EthnicityHispanic == 1,
      "Hispanic or Latino",
      ifelse(
        EthnicityNotreported == 0 & EthnicityUnknown == 0,
        "Not Hispanic or Latino",
        "Not reported and unknown"
      )
    )
  ) %>% factor(
    levels = c("Hispanic or Latino", "Not Hispanic or Latino", "Not reported and unknown", "Others")
  )



dat.long.twophase.sample$minority_label <-
  with(
    dat.long.twophase.sample,
    factor(MinorityInd,
      levels = c(0, 1),
      labels = c("White Non-Hispanic", "Comm. of Color")
    )
  )

dat.long.twophase.sample$age_minority_label <-
  with(
    dat.long.twophase.sample,
    factor(paste0(Senior, MinorityInd),
      levels = c("01", "00", "11", "10"),
      labels = c(
        paste0("Age < ", age.cutoff, " Comm. of Color"),
        paste0("Age < ", age.cutoff, " White Non-Hispanic"),
        paste0("Age >= ", age.cutoff, " Comm. of Color"),
        paste0("Age >= ", age.cutoff, " White Non-Hispanic")
      )
    )
  )

if(study_name=="ENSEMBLE" | study_name=="MockENSEMBLE") {
  dat.long.twophase.sample$country_label <- factor(sapply(dat.long.twophase.sample$Country, function(x) {
    names(countries.ENSEMBLE)[countries.ENSEMBLE==x]
  }), levels = names(countries.ENSEMBLE))
}

if(study_name!="COVE" & study_name!="MockCOVE") {dat.long.twophase.sample$hiv_label <- factor(sapply(dat.long.twophase.sample$HIVinfection, function(x) {
  ifelse(x,
         "HIV Positive",
         "HIV Negative")
}), levels=c("HIV Negative", "HIV Positive"))
}

dat.long.twophase.sample$race <- as.factor(dat.long.twophase.sample$race)
dat.twophase.sample$race <- as.factor(dat.twophase.sample$race)

dat.long.twophase.sample$Ptid <- as.character(dat.long.twophase.sample$Ptid) 
dat.twophase.sample$Ptid <- as.character(dat.twophase.sample$Ptid) 


dat.long.twophase.sample <- filter(dat.long.twophase.sample, assay %in% assay_immuno)


saveRDS(as.data.frame(dat.long.twophase.sample),
  file = here("data_clean", "long_twophase_data.rds")
)
saveRDS(as.data.frame(dat.twophase.sample),
  file = here("data_clean", "twophase_data.rds")
)
