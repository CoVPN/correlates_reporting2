##################################################
# obligatory to append to the top of each script #
renv::activate(project = here::here(".."))
source(here::here("..", "_common.R")) #
##################################################

library(survey)
library(tidyverse)
library(PropCIs)
library(dplyr, warn.conflicts = FALSE)
# Suppress summarise info
options(dplyr.summarise.inform = FALSE)
# For stratum with 1 ppt
options(survey.lonely.psu="adjust")
source(here::here("code", "make_functions.R"))


###################################################
#                  Parameters                     #
###################################################
# To select which tables are included in the report.
# Also to modify the headers and footers for each table.

randomsubcohort <- case_when(study_name %in% c("COVE", "MockCOVE") ~ "This table summarizes the 
      random subcohort, which was randomly sampled from the per-protocol cohort. The 
      sampling was stratified by 24 strata defined by enrollment characteristics: Assigned 
      treatment arm $\\\\times$ Baseline SARS-CoV-2 naive vs. non-naive status 
      (defined by serostatus and NAAT testing) $\\\\times$ Randomization strata 
      (Age < 65 and at-risk, Age < 65 and not at-risk, Age $\\\\geq 65)\\\\times$ 
      Communities of color (Yes/No) defined by White Non-Hispanic vs. all 
      others (following the primary COVE trial paper).",
                             
      study_name %in% c("ENSEMBLE", "MockENSEMBLE") ~ "This table summarizes characteristics of 
      per-protocol participants in the immunogenicity subcohort, which was randomly 
      sampled from the study cohort. The sampling was stratified by 
      strata defined by enrollment characteristics: Assigned randomization arm $\\\\times$ 
      Baseline SARS-CoV-2 seronegative vs. seropositive $\\\\times$ Randomization strata. 
      The U.S. subcohort includes 8 baseline demographic strata; the Latin America 
      and South Africa subcohorts each include 4 baseline demographic strata.",
      
      TRUE~"This table summarizes characteristics of 
      per-protocol participants in the immunogenicity subcohort, which was randomly 
      sampled from the study cohort.")


tlf <-
  list(
    tab_dm_neg_ph1 = list(
      table_header = "Demographic and Clinical Characteristics at Baseline in 
      the Per-Protocol Cohort (Correlates Cohort)",
      table_footer = randomsubcohort,
      deselect = "subgroup",
      pack_row = "subgroup",
      col1="7cm"
    ),
    
    
    tab_dm_neg = list(
      table_header = "Demographic and Clinical Characteristics at Baseline in 
      the Per-Protocol Cohort (Phase 2 Correlates Cohort)",
      table_footer = randomsubcohort,
      deselect = "subgroup",
      pack_row = "subgroup",
      col1="7cm"
    ),
    
    tab_strtm1 = NULL,
    
    case_vacc_neg = list(
      table_header = "Antibody levels in the per-protocol cohort (BPZE1 recipients)",
      table_footer =c(
        paste(paste(sprintf("Cases for Day %s markers are per-protocol BPZE1 recipients 
      with the symptomatic primary endpoint diagnosed starting %s day(s) 
      after the Day %s study visit.", config.cor$tpeak, config.cor$tpeaklag, config.cor$tpeak), collapse=" "),
          "Non-cases/Controls are per-protocol BPZE1 recipients with no endpoint diagnosis by the time of data-cut.")
        ),

      col_name = c("Visit", "Marker", "N", "Resp rate", "GMT/GMC", "N",
                   "Resp rate", "GMT/GMC", "Resp Rate\nDifference", "GMTR/GMCR"),
      header_above1 = c(" "=2, "Cases*" = 3, "Non-Cases/Control" = 3,
                        "Comparison" = 2),
      header_above2 = c(" "=2,
                        "BPZE1 Recipients" = 8),
      col1="1cm"),
    
    case_plcb_neg = list(
      table_header = "Antibody levels in the per-protocol cohort (PBO recipients)",
      table_footer =c(
        paste(paste(sprintf("Cases for Day %s markers are per-protocol PBO recipients 
      with the symptomatic infection primary endpoint diagnosed starting %s day(s) 
      after the Day %s study visit.", config.cor$tpeak, config.cor$tpeaklag, config.cor$tpeak), collapse=" "),
              "Non-cases/Controls are per-protocol PBO recipients with no endpoint diagnosis by the time of data-cut.")
        ),
      
      col_name = c("Visit", "Marker", "N", "Resp rate", "GMT/GMC", "N",
                   "Resp rate", "GMT/GMC", "Resp Rate\nDifference", "GMTR/GMCR"),
      header_above1 = c(" "=2, "Cases*" = 3, "Non-Cases/Control" = 3,
                        "Comparison" = 2),
      header_above2 = c(" "=2,
                        "PBO Recipients" = 8),
      col1="1cm")

    )

tlf <- lapply(tlf, function(x){
    x$table_header <- paste0(COR,": ", x$table_header)
    x
  }
)

timepoints <- config$timepoints

labels.age <- case_when(study_name %in% c("ENSEMBLE", "MockENSEMBLE") ~ c("Age 18 - 59", "Age $\\geq$ 60"), 
                        TRUE~ c("Age $<$ 65", "Age $\\geq$ 65"))

labels.minor <- c("Communities of Color", "White Non-Hispanic")

labels.BMI <- c("Underweight BMI < 18.5", "Normal 18.5 $\\leq$ BMI < 25", 
                "Overweight 25 $\\leq$ BMI < 30", "Obese BMI $\\geq$ 30")

labels.time <- labels.time[times]

if ("BbindN" %in% names(dat_proc) & any(grepl("bind", assays))) assays <- union(assays, "bindN")


bindN <- "Anti N IgG (BAU/ml)"
names(bindN) <- "bindN"
labels.assays.short.tabular <- c(labels.assays.short.tabular, bindN)

labels.assays.short <- labels.assays.short.tabular[assays]

# redefines what is in _common.R to use shorter names
labels.assays.long <- data.frame (purrr::imap_dfc(labels.assays.short, ~ paste0(labels.assays.short[.y], ": ", labels.time)))
rownames(labels.assays.long) <- names(labels.time)

visits <- names(labels.time)[!grepl("Delta", names(labels.time))]
assays_col <- as.vector(outer(visits, assays, paste0))

labels.assays <- expand.grid(
  time = rownames(labels.assays.long),
  marker = colnames(labels.assays.long),
  stringsAsFactors = FALSE
) %>%
  rowwise() %>%
  mutate(
    label.long = labels.assays.long[time, marker],
    label.short = sapply(labels.assays.short, as.character)[marker],
    Marker = strsplit(as.character(label.long), ": ", fixed = T)[[1]][1],
    Visit = strsplit(as.character(label.long), ": ", fixed = T)[[1]][2],
    colname = paste0(time, marker)
  )

resp.lb <- expand.grid(
  time = visits, marker = assays,
  ind = c("Resp", "FR2", "FR4", "2lloq", "4lloq", "2llod", "4llod"), stringsAsFactors = F
) %>%
  mutate(Ind = case_when(
    ind == "FR2" ~ "% 2-Fold Rise",
    ind == "FR4" ~ "% 4-Fold Rise",
    ind == "Resp" ~ "Responder",
    ind == "2lloq" ~ "% Greater than 2xLLOQ",
    ind == "4lloq" ~ "% Greater than 4xLLOQ",
    ind == "2llod" ~ "% Greater than 2xLLOD",
    ind == "4llod" ~ "% Greater than 4xLLOD"
  )) 

labels_all <- full_join(labels.assays, resp.lb, by = c("time", "marker")) %>% 
  mutate(mag_cat = colname, resp_cat = paste0(colname, ind))



###################################################
#                Clean the Data                   #
###################################################

### Table 1. Demographics 
# Output: tab_dm
# Select the covariates to be summarised.
# num_v are columns from ds_long;
# cat_v are rows of `subgroup`


# dat_proc was made in _common.R
dat <- dat_proc


# The stratified random cohort for immunogenicity
if (study_name=="ILIAD_IB202P") {
    ds_s <- dat %>%
    mutate(
      raceC = case_when(Black==1 ~ "Black or African American",
                        Asian==1 ~ "Asian",
                        Multiracial==1 ~ "Multiracail",
                        Other==1 ~ "Other"),
      # HighRiskC = ifelse(HighRiskInd == 1, "At-risk", "Not at-risk"),
      AgeC = ifelse(Age >= 65, labels.age[2], labels.age[1]),
      SexC = ifelse(Sex == 1, "Female", "Male"),
      # AgeRiskC = paste(AgeC, HighRiskC),
      AgeSexC = paste(AgeC, SexC),
      # AgeMinorC = ifelse(is.na(MinorityC), NA, paste(AgeC, MinorityC)),
      `Baseline SARS-CoV-2` = "Negative",
      Arm = factor(Trt, levels = c("BPZE1", "PBO")),
      All = "All participants")
  
} else {

  ds_s <- dat %>%
  mutate(
    raceC = as.character(race),
    ethnicityC = case_when(EthnicityHispanic==1 ~ "Hispanic or Latino",
                           EthnicityHispanic==0 & EthnicityNotreported==0 & 
                           EthnicityUnknown==0 ~ "Not Hispanic or Latino",
                           EthnicityNotreported==1 | 
                           EthnicityUnknown==1 ~ "Not reported and unknown "),
    RaceEthC = case_when(
      WhiteNonHispanic==1 ~ "White Non-Hispanic ",
      TRUE ~ raceC
    ),
    MinorityC = case_when(
      MinorityInd == 1 ~ "Communities of Color",
      MinorityInd == 0 ~ "White Non-Hispanic"
    ),
    HighRiskC = ifelse(HighRiskInd == 1, "At-risk", "Not at-risk"),
    AgeC = ifelse(Senior == 1, labels.age[2], labels.age[1]),
    SexC = ifelse(Sex == 1, "Female", "Male"),
    AgeRiskC = paste(AgeC, HighRiskC),
    AgeSexC = paste(AgeC, SexC),
    AgeMinorC = ifelse(is.na(MinorityC), NA, paste(AgeC, MinorityC)),
    `Baseline SARS-CoV-2` = factor(ifelse(Bserostatus == 1, "Positive", "Negative"),
                                   levels = c("Negative", "Positive")
    ),
    Arm = factor(ifelse(Trt == 1, "Vaccine", "Placebo"), 
                 levels = c("Vaccine", "Placebo")),
    
    All = "All participants"
    )
}



# Step2: Responders
# Post baseline visits
pos.cutoffs <- assay_metadata$pos.cutoff
names(pos.cutoffs) <- assay_metadata$assay

ds <- getResponder(ds_s, times=grep("Day", times, value=T), 
                   assays=assays, pos.cutoffs = pos.cutoffs)

subgrp <- c(
  All = "All participants", 
  AgeC = "Age",
  BMI="BMI",
  HighRiskC = "Risk for Severe Covid-19",
  AgeRiskC = "Age, Risk for Severe Covid-19",
  AgeRisk1 = paste0(labels.age[1], ", Risk for Severe Covid-19"),
  AgeRisk2 = paste0(labels.age[2], ", Risk for Severe Covid-19"),
  SexC = "Sex", 
  AgeSexC = "Age, sex",
  ethnicityC = "Hispanic or Latino ethnicity", 
  RaceEthC = "Race",
  MinorityC = "Underrepresented minority status",
  AgeMinorC = "Age, Communities of color",
  URMC = "Underrepresented Minority Status in the U.S.",
  AgeURM = "Age, Underrepresented Minority Status in the U.S.",
  CountryC = "Country",
  HIVC = "HIV Infection"
)


###################################################
#             Generating the Tables               #
###################################################

# Setup empty tables 
for (i in names(tlf)){
  assign(i, NULL)
}


if (study_name %in% c("COVE", "MockCOVE")) {
  num_v1 <- c("Age") # Summaries - Mean & Range
  num_v2 <- c("BMI") # Summaries - Mean & St.d
  cat_v <- c("AgeC", "SexC", "raceC", "ethnicityC", "HighRiskC", "AgeRiskC", "MinorityC")
} else if (study_name %in% c("ENSEMBLE", "MockENSEMBLE")) {
  num_v1 <- c("Age") # Summaries - Mean & Range
  num_v2 <- NULL # Summaries - Mean & St.d
  cat_v <- c("AgeC", "SexC", "raceC", "ethnicityC", 
             "HighRiskC", "AgeRiskC", "URMC",  "CountryC", "HIVC", "BMI")
} else if (study_name %in% c("AZD1222")) {
  num_v1 <- c("Age") # Summaries - Mean & Range
  num_v2 <- NULL # Summaries - Mean & St.d
  cat_v <- c("AgeC", "SexC", "CountryC", "raceC", "ethnicityC", "HighRiskC", "AgeRiskC")
} else if (study_name %in% c("UK302")) {
  num_v1 <- c("Age") # Summaries - Mean & Range
  num_v2 <- c("BMI") # Summaries - Mean & St.d
  cat_v <- c("AgeC", "SexC", "raceC", "ethnicityC")
} else if (study_name %in% c("ILIAD_IB202P")) {
    num_v1 <- c("Age") # Summaries - Mean & Range
    num_v2 <- NULL # Summaries - Mean & St.d
    cat_v <- c("AgeC", "SexC", "raceC")
} else{ # Keeping the minimal
  num_v1 <- c("Age") # Summaries - Mean & Range
  num_v2 <- NULL # Summaries - Mean & St.d
  cat_v <- c("AgeC", "SexC", "raceC", "ethnicityC", "HighRiskC", "AgeRiskC")
} 

ds_long_ttl <- ds %>%
  dplyr::filter(!!as.name(config.cor$ph2)==1) %>%
  # dplyr::filter(!!as.name(paste0("ph2.D", tpeak))) %>% 
  bind_rows(mutate(., Arm="Total")) %>% 
  # mutate(AgeRiskC = ifelse(grepl("$\\geq$ 65", AgeRiskC, fixed=T), "Age $\\geq$ 65 ", AgeRiskC)) %>% 
  mutate_all(as.character) %>% 
  pivot_longer(all_of(c(num_v1, num_v2, cat_v)), names_to="subgroup", values_to="subgroup_cat")

ds_long_ttl_ph1 <- ds %>%
  dplyr::filter(Perprotocol==1) %>%
  # dplyr::filter(!!as.name(paste0("ph2.D", tpeak))) %>% 
  bind_rows(mutate(., Arm="Total")) %>% 
  # mutate(AgeRiskC = ifelse(grepl("$\\geq$ 65", AgeRiskC, fixed=T), "Age $\\geq$ 65 ", AgeRiskC)) %>% 
  mutate_all(as.character) %>% 
  pivot_longer(all_of(c(num_v1, num_v2, cat_v)), names_to="subgroup", values_to="subgroup_cat")


# Calculate % for categorical covariates
dm_cat <- inner_join(
  ds_long_ttl %>%
    group_by(`Baseline SARS-CoV-2`, Arm, subgroup, subgroup_cat) %>%
    summarise(n = n(), .groups = 'drop'),
  ds_long_ttl %>%
    group_by(`Baseline SARS-CoV-2`, Arm, subgroup) %>%
    summarise(N = n(), .groups = 'drop'),
  by = c("Baseline SARS-CoV-2", "Arm", "subgroup")
) %>%
  mutate(pct = n / N,
         rslt1 = sprintf("%s (%.1f%%)", n, n / N * 100), 
         rslt2 = sprintf("%s/%s = %.1f%%", n, N, n / N * 100)) %>% 
  dplyr::filter(subgroup %in% cat_v) 

dm_cat_ph1 <- inner_join(
  ds_long_ttl_ph1 %>%
    group_by(`Baseline SARS-CoV-2`, Arm, subgroup, subgroup_cat) %>%
    summarise(n = n(), .groups = 'drop'),
  ds_long_ttl_ph1 %>%
    group_by(`Baseline SARS-CoV-2`, Arm, subgroup) %>%
    summarise(N = n(), .groups = 'drop'),
  by = c("Baseline SARS-CoV-2", "Arm", "subgroup")
) %>%
  mutate(pct = n / N,
         rslt1 = sprintf("%s (%.1f%%)", n, n / N * 100), 
         rslt2 = sprintf("%s/%s = %.1f%%", n, N, n / N * 100)) %>% 
  dplyr::filter(subgroup %in% cat_v) 

# Calculate mean and range for numeric covariates
dm_num <- ds_long_ttl %>%
  dplyr::filter(subgroup %in% c(num_v1, num_v2)) %>% 
  mutate(subgroup_cat=as.numeric(subgroup_cat)) %>%
  group_by(`Baseline SARS-CoV-2`, Arm, subgroup) %>%
  summarise(
    min = min(subgroup_cat, na.rm = T), 
    max = max(subgroup_cat, na.rm = T),
    mean = mean(subgroup_cat, na.rm = T),
    sd = sd(subgroup_cat, na.rm = T), 
    rslt1 = sprintf("%.1f (%.1f, %.1f)", mean, min, max),
    rslt2 = sprintf("%.1f $\\pm$ %.1f", mean, sd),
    N = n(),
    .groups = 'drop'
  ) %>% 
  mutate(subgroup_cat = case_when(subgroup %in% num_v1 ~ "Mean (Range)",
                                  subgroup %in% num_v2 ~ "Mean $\\pm$ SD"),
         subgroup=ifelse(subgroup=="Age", "AgeC", subgroup))

dm_num_ph1 <- ds_long_ttl_ph1 %>%
  dplyr::filter(subgroup %in% c(num_v1, num_v2)) %>% 
  mutate(subgroup_cat=as.numeric(subgroup_cat)) %>%
  group_by(`Baseline SARS-CoV-2`, Arm, subgroup) %>%
  summarise(
    min = min(subgroup_cat, na.rm = T), 
    max = max(subgroup_cat, na.rm = T),
    mean = mean(subgroup_cat, na.rm = T),
    sd = sd(subgroup_cat, na.rm = T), 
    rslt1 = sprintf("%.1f (%.1f, %.1f)", mean, min, max),
    rslt2 = sprintf("%.1f $\\pm$ %.1f", mean, sd),
    N = n(),
    .groups = 'drop'
  ) %>% 
  mutate(subgroup_cat = case_when(subgroup %in% num_v1 ~ "Mean (Range)",
                                  subgroup %in% num_v2 ~ "Mean $\\pm$ SD"),
         subgroup=ifelse(subgroup=="Age", "AgeC", subgroup))

char_lev <- c(labels.age, "Mean (Range)","Mean $\\pm$ SD",
              "Female","Male", "White", "Black or African American",
              "Asian", "American Indian or Alaska Native",
              "Native Hawaiian or Other Pacific Islander", "Multiracial",
              "Other", "Not reported and unknown", 
              "White Non-Hispanic", "Communities of Color",
              "Hispanic or Latino","Not Hispanic or Latino",
              "Not reported and unknown ","At-risk","Not at-risk",
              paste(labels.age[1],"At-risk"), paste(labels.age[1], "Not at-risk"), 
              paste(labels.age[2],"At-risk"), paste(labels.age[2], "Not at-risk"),
              paste(labels.age[2], ""), 
              # "Communities of Color", "White Non-Hispanic",
              labels.countries.ENSEMBLE,
              "Negative", "Positive", labels.BMI)

tab_dm <- bind_rows(dm_cat, dm_num) %>%
  mutate(rslt = case_when(subgroup %in% cat_v ~ rslt1,
                          subgroup %in% num_v1 ~ rslt1,
                          subgroup %in% num_v2 ~ rslt2)) %>%
  mutate(subgroup=ifelse(subgroup %in% c("MinorityC", "raceC"), "RaceEthC", subgroup)) %>% 
  dplyr::filter(subgroup_cat %in% char_lev) %>% 
  inner_join(ds_long_ttl %>% 
               distinct(`Baseline SARS-CoV-2`, Arm, Ptid) %>% 
               group_by(`Baseline SARS-CoV-2`, Arm) %>%
               summarise(tot = n()),
             by = c("Baseline SARS-CoV-2", "Arm")) %>% 
  mutate(Arm = paste0(Arm, "\n(N = ", tot, ")"), subgroup=subgrp[subgroup]) %>%
  pivot_wider(id_cols=c(`Baseline SARS-CoV-2`, subgroup, subgroup_cat),
              names_from = Arm, 
              names_sort = T,
              values_from = c(rslt)) %>%
  mutate(Characteristics = factor(subgroup_cat, levels=char_lev),
         subgroup=factor(subgroup, levels=subgrp)) %>%
  arrange(`Baseline SARS-CoV-2`, subgroup, Characteristics)

tab_dm_ph1 <- bind_rows(dm_cat_ph1, dm_num_ph1) %>%
  mutate(rslt = case_when(subgroup %in% cat_v ~ rslt1,
                          subgroup %in% num_v1 ~ rslt1,
                          subgroup %in% num_v2 ~ rslt2)) %>%
  mutate(subgroup=ifelse(subgroup %in% c("MinorityC", "raceC"), "RaceEthC", subgroup)) %>% 
  dplyr::filter(subgroup_cat %in% char_lev) %>% 
  inner_join(ds_long_ttl_ph1 %>% 
               distinct(`Baseline SARS-CoV-2`, Arm, Ptid) %>% 
               group_by(`Baseline SARS-CoV-2`, Arm) %>%
               summarise(tot = n()),
             by = c("Baseline SARS-CoV-2", "Arm")) %>% 
  mutate(Arm = paste0(Arm, "\n(N = ", tot, ")"), subgroup=subgrp[subgroup]) %>%
  pivot_wider(id_cols=c(`Baseline SARS-CoV-2`, subgroup, subgroup_cat),
              names_from = Arm, 
              names_sort = T,
              values_from = c(rslt)) %>%
  mutate(Characteristics = factor(subgroup_cat, levels=char_lev),
         subgroup=factor(subgroup, levels=subgrp)) %>%
  arrange(`Baseline SARS-CoV-2`, subgroup, Characteristics)

if ("Negative" %in% tab_dm$`Baseline SARS-CoV-2`){
  tab_dm_neg <- tab_dm %>% 
    dplyr::filter(`Baseline SARS-CoV-2` == "Negative") %>% 
    select_if(~ !all(is.na(.))) %>% 
    select_at(c("subgroup", "Characteristics", 
                grep("Vaccine" ,names(.), value = T),
                grep("Placebo" ,names(.), value = T),
                grep("Total" ,names(.), value = T)))
  
  tab_dm_neg_ph1 <- tab_dm_ph1 %>% 
    dplyr::filter(`Baseline SARS-CoV-2` == "Negative") %>% 
    select_if(~ !all(is.na(.))) %>% 
    select_at(c("subgroup", "Characteristics", 
                grep("Vaccine" ,names(.), value = T),
                grep("Placebo" ,names(.), value = T),
                grep("Total" ,names(.), value = T)))
}

if ("Positive" %in% tab_dm$`Baseline SARS-CoV-2`){
  tab_dm_pos <- tab_dm %>% 
    dplyr::filter(`Baseline SARS-CoV-2` == "Positive") %>% 
    select_if(~ !all(is.na(.))) %>% 
    select_at(c("subgroup", "Characteristics", 
                grep("Vaccine" ,names(.), value = T),
                grep("Placebo" ,names(.), value = T),
                grep("Total" ,names(.), value = T)))
  
  tab_dm_pos_ph1 <- tab_dm_ph1 %>% 
    dplyr::filter(`Baseline SARS-CoV-2` == "Positive") %>% 
    select_if(~ !all(is.na(.))) %>% 
    select_at(c("subgroup", "Characteristics", 
                grep("Vaccine" ,names(.), value = T),
                grep("Placebo" ,names(.), value = T),
                grep("Total" ,names(.), value = T)))
}
print("Done with table 1") 


# Cases & Non-cases
# 1. Moderna 
# !!as.name(paste0("EarlyendpointD",nonCaseD))==0
# nonCaseD always 57
# 2. AZ
# !!as.name(paste0("EarlyendpointD",nonCaseD))==0
# nonCaseD depends on timepoint, 29 or 57 
# 3. All other
# AnyinfectionD1==0
# if (study_name %in% c("COVE", "MockCOVE", "MockENSEMBLE", "ENSEMBLE", "PREVENT19", "VAT08m")){

if (study_name=="AZD1222"){
  nonCaseD <- tpeak
} else{ #study_name %in% c("COVE", "MockCOVE", "MockENSEMBLE", "ENSEMBLE", "PREVENT19", "VAT08m")
  nonCaseD <- timepoints[length(timepoints)]
}


ds <- ds %>% 
  mutate(Case = case_when(!!as.name(config.cor$EventIndPrimary)==1 ~ "Cases",
                          TRUE ~ "Non-Cases"))

# Generate a full table with all the estimates: response rates, GMT, GMTR, etc.
# (Per Peter's email Feb 5, 2021)
# Cases vs Non-cases

  sub.by <- c("Arm")
  ds.i <- filter(ds, !!as.name(config.cor$ph1)==1)
  resp.v <- intersect(grep("Resp", names(ds), value = T), 
                      grep(config.cor$tpeak, names(ds), value = T))
  
  gm.v <- intersect(assays_col, names(ds))
  resp.v <- resp.v[sapply(resp.v, function(x)any(!is.na(ds.i[,x])))]
  gm.v <- gm.v[sapply(gm.v, function(x)any(!is.na(ds.i[,x])))]
  
  subs <- "Case"
  comp.i <- c("Cases", "Non-Cases")


  ds.l.resp <- filter(ds, !!as.name(config.cor$ph1)==1) %>% 
    pivot_longer(cols=all_of(resp.v), names_to = "resp_cat", values_to = "resp_value") %>% 
    mutate(assay=gsub("Resp", "", resp_cat)) 
  
  ds.l.mag <- filter(ds, !!as.name(config.cor$ph1)==1) %>% 
    # select_at(c("Ptid", gm.v)) %>% 
    pivot_longer(cols=all_of(gm.v), names_to = "mag_cat", values_to = "mag_value") %>% 
    mutate(assay=mag_cat) %>% 
    left_join(labels_all %>% distinct(mag_cat, Marker)) %>% 
    group_by(Marker, Ptid) %>% 
    mutate(mag_filt = sum(!is.na(mag_value))) %>% 
    filter(mag_filt!=1) %>% 
    ungroup()
  
    
  n_GM <- 2
  rpcnt_case <- ds.l.resp %>% 
    group_by(Arm, resp_cat, Case) %>% 
    summarise(N=sum(!is.na(resp_value)), Pos=sum(resp_value, na.rm = T), .groups="drop") %>% 
    rowwise() %>% 
    mutate(response=Pos/N, ci_l=exactci(Pos, N, .95)$conf.int[1], ci_u=exactci(Pos, N, .95)$conf.int[2], 
           rslt=sprintf("%s/%s = %.1f%%\n(%.1f%%, %.1f%%)", Pos, N, response*100, ci_l*100, ci_u*100)) %>% 
    inner_join(distinct(labels_all, resp_cat, Visit, Marker, Ind), by = "resp_cat") %>% 
    filter(N!=0)
  
  rgm_case <- ds.l.mag %>% 
    group_by(Arm, mag_cat, Case) %>% 
    summarise(N=sum(!is.na(mag_value)), lgmt=mean(mag_value, na.rm=T), lsd=sd(mag_value, na.rm=T), .groups="drop") %>% 
    rowwise() %>% 
    mutate(gmt=ifelse(grepl("FS", mag_cat), lgmt, 10^lgmt), 
           ci_l=ifelse(grepl("FS", mag_cat), lgmt+qt(.025, N-1)*lsd/sqrt(N), 10^(lgmt+qt(.025, N-1)*lsd/sqrt(N))), 
           ci_u=ifelse(grepl("FS", mag_cat), lgmt+qt(.975, N-1)*lsd/sqrt(N), 10^(lgmt+qt(.975, N-1)*lsd/sqrt(N))), 
           `GMT/GMC`=sprintf(sprintf("%%.%sf\n(%%.%sf, %%.%sf)", n_GM, n_GM, n_GM), gmt, ci_l, ci_u)) %>% 
    inner_join(distinct(labels_all, mag_cat, Visit, Marker), by = "mag_cat") 
  # filter(N!=0)
  
  rgmt_excl <- rgm_case %>% 
    filter(N==0) %>% 
    select(mag_cat)
  
  rgmt_case <- ds.l.mag %>%
    anti_join(rgmt_excl) %>% 
    mutate(Case=factor(Case)) %>% 
    mutate(Case=relevel(Case, ref="Cases")) %>%
    filter(!is.na(mag_value)) %>% 
    group_by(Arm, mag_cat) %>% 
    summarise(N=sum(!is.na(mag_value)), 
              est=glm(mag_value~Case)$coefficients[2], 
              ci_l=unique(ifelse(class(try(confint(glm(mag_value~Case))))=="try-error", NA, confint(glm(mag_value~Case))[2,1])),
              ci_u=unique(ifelse(class(try(confint(glm(mag_value~Case))))=="try-error", NA, confint(glm(mag_value~Case))[2,2])),.groups="drop") %>% 
    mutate(comp="Non-Cases vs Cases", 
           `Ratios of GMT/GMC`=ifelse(grepl("FS", mag_cat),
                                      sprintf("%.2f\n(%.2f, %.2f)", est, ci_l, ci_u),
                                      sprintf("%.2f\n(%.2f, %.2f)", 10^est, 10^ci_l, 10^ci_u))) %>% 
    inner_join(distinct(labels_all, mag_cat, Visit, Marker))
  
  
  rrdiff_case <- rpcnt_case %>%
    # dplyr::filter(subgroup %in% subs & grepl("Resp",resp_cat)) %>%
    mutate(groupn = 2-match(Case, comp.i)%%2) %>%
    pivot_wider(id_cols = c(Arm, Visit, Marker, Ind),
                names_from = groupn, values_from = c(response, ci_l, ci_u), names_sep = "")
  
  responseNA <- setdiff(as.vector(outer(c("response", "ci_l", "ci_u"), 1:2, paste0)), names(rrdiff_case))
  rrdiff_case[, responseNA] <- NA
  
  rrdiff_case <- rrdiff_case %>%
    mutate(Estimate = response1-response2,
           ci_l = Estimate-sqrt((response1-ci_l1)^2+(response2-ci_u2)^2),
           ci_u = Estimate+sqrt((response1-ci_u1)^2+(response2-ci_l2)^2),
           rrdiff = ifelse(!is.na(Estimate),
                           sprintf("%s\n(%s, %s)", round(Estimate, 3), round(ci_l, 3), round(ci_u, 3)),
                           "-"))
  
  tab_case <- full_join(rpcnt_case, rgm_case,
                        by = c("Case", "N", "Marker", "Visit", "Arm")) %>%
    filter(N!=0) %>% 
    pivot_wider(id_cols = c(Arm, Marker, Visit),
                names_from = Case,
                values_from = c(N, rslt, `GMT/GMC`)) %>%
    full_join(rrdiff_case, by = c("Arm", "Marker", "Visit")) %>%
    full_join(rgmt_case, by = c("Arm", "Marker", "Visit")) %>% 
    mutate(Marker=factor(Marker, levels=assay_metadata$assay_label_short))
  
  
  if(length(comp_NA <- setdiff(comp.i, rpcnt_case$Case))!=0){
    tab_case <- tab_case %>%
      mutate(!!paste0("N_", comp_NA) := 0,
             !!paste0("rslt_", comp_NA) := "-",
             !!paste0("GMT/GMC_", comp_NA) :="-",
             `Ratios of GMT/GMC`=replace_na(`Ratios of GMT/GMC`, "-"))
  }else{
    tab_case <- tab_case %>%
      mutate_at(vars(starts_with("N_")), replace_na, replace=0) %>%
      mutate_at(vars(starts_with("rslt_")), replace_na, replace="-") %>%
      mutate_at(vars(starts_with("GMT/GMC_")), replace_na, replace="-") %>%
      mutate(`Ratios of GMT/GMC`=replace_na(`Ratios of GMT/GMC`, "-")) %>% 
      mutate(rrdiff=replace_na(rrdiff, "-")) 
  }
  
  tab_case <- tab_case %>%
    select(Arm, Visit, Marker, `N_Cases`, `rslt_Cases`,
           `GMT/GMC_Cases`, `N_Non-Cases`, `rslt_Non-Cases`, `GMT/GMC_Non-Cases`,
           rrdiff, `Ratios of GMT/GMC`) 
  
  tab_gmtr <- tab_case %>% 
    filter(grepl("over", Visit)) %>% 
    select(Arm, Visit, Marker, `N_Cases`, `GMT/GMC_Cases`, `N_Non-Cases`, `GMT/GMC_Non-Cases`,`Ratios of GMT/GMC`) %>% 
    arrange(Arm, Visit, Marker)
  
  tab_case <- tab_case %>% 
    filter(!grepl("over", Visit)) %>% 
    arrange(Arm, Visit, Marker)
  


print("Done with table6")

case_vacc_neg <- tab_case %>%
  dplyr::filter(Arm == "BPZE1") %>% 
  select(-Arm) %>% 
  arrange(Visit, Marker)

case_plcb_neg <- tab_case %>%
  dplyr::filter(Arm == "PBO") %>% 
  select(-Arm) %>% 
  arrange(Visit, Marker)


print("Done with all tables") 


if(tpeak!=timepoints[1]){
  tlf <- tlf[!names(tlf) %in% c("tab_dm_neg_ph1", "tab_dm_pos_ph1", "tab_dm_neg", "tab_dm_pos")]
}


# path for tables
save.results.to <- here::here("output")
if (!dir.exists(save.results.to))  dir.create(save.results.to)
save.results.to <- paste0(here::here("output"), "/", attr(config,"config"))

if (!dir.exists(save.results.to))  dir.create(save.results.to)
print(paste0("save.results.to equals ", save.results.to))

save(list = c("tlf", names(tlf)), file = file.path(save.results.to, sprintf("Tables%s.Rdata", ifelse(exists("COR"), gsub("D15", "D22", COR), ""))))
