##################################################
# obligatory to append to the top of each script #
renv::activate(project = here::here(".."))
source(here::here("..", "_common.R")) #
##################################################

library(survey)
library(tidyverse)
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
      the Baseline SARS-CoV-2 Negative Per-Protocol Cohort (Immunogenicity Cohort)",
      table_footer = randomsubcohort,
      deselect = "subgroup",
      pack_row = "subgroup",
      col1="7cm"
    ),
    
    tab_dm_pos_ph1 = list(
      table_header = "Demographic and Clinical Characteristics at Baseline in 
      the Baseline SARS-CoV-2 Positive Per-Protocol Cohort (Immunogenicity Cohort)",
      table_footer = randomsubcohort,
      deselect = "subgroup",
      pack_row = "subgroup",
      col1="7cm"
    ),
    
    tab_dm_neg = list(
      table_header = "Demographic and Clinical Characteristics at Baseline in 
      the Baseline SARS-CoV-2 Negative Per-Protocol Cohort (Phase 2 Immunogenicity Cohort)",
      table_footer = randomsubcohort,
      deselect = "subgroup",
      pack_row = "subgroup",
      col1="7cm"
    ),
    
    tab_dm_pos = list(
      table_header = "Demographic and Clinical Characteristics at Baseline in 
      the Baseline SARS-CoV-2 Positive Per-Protocol Cohort (Phase 2 Immunogenicity Cohort)",
      table_footer = randomsubcohort,
      deselect = "subgroup",
      pack_row = "subgroup",
      col1="7cm"
    ),
    
    tab_strtm1 = list(
      table_header = "",
      deselect = "Arm",
      pack_row = "Arm"
    ),
    
    tab_strtm2 = list(
      table_header = "",
      deselect = "Arm",
      pack_row = "Arm"
    ),
    
    tab_case_cnt = list(
      table_header = "Availability of immunogenicity data by case status",
      deselect = "Arm",
      pack_row = "Arm",
      table_footer = c("The $+$ (available) and $-$ (unavailable) in the column 
                       labels refer to the availability of the baseline, D29 and D57 markers, respectively."),
      col1="7cm"),  
    
    tab_days = list(
      table_header = sprintf("Duration from vaccination to D%s visit in the 
                             baseline SARS-CoV-2 negative per-protocol cohort", config.cor$tpeak),
      deselect = "Arm",
      pack_row = "Arm"
    ),
    
    case_vacc_neg = list(
      table_header = "Antibody levels in the baseline SARS-CoV-2 negative
      per-protocol cohort (vaccine recipients)",
      table_footer =c(
        paste(paste(sprintf("Cases for Day %s markers are baseline negative per-protocol vaccine recipients 
      with the symptomatic infection COVID-19 primary endpoint diagnosed starting %s day(s) 
      after the Day %s study visit.", config.cor$tpeak, config.cor$tpeaklag, config.cor$tpeak), collapse=" "),
          "Non-cases/Controls are baseline negative per-protocol vaccine recipients sampled into the random subcohort 
      with no COVID-19 endpoint diagnosis by the time of data-cut."),
          "N is the number of cases sampled into the subcohort within baseline covariate strata.",
          "The denominator in Resp Rate is the number of participants in the whole per-protocol cohort within baseline 
      covariate strata, calculated using inverse probability weighting."),

      col_name = c("Visit", "Marker", "N", "Resp rate", "GMT/GMC", "N",
                   "Resp rate", "GMT/GMC", "Resp Rate\nDifference", "GMTR/GMCR"),
      header_above1 = c(" "=2, "Cases*" = 3, "Non-Cases/Control" = 3,
                        "Comparison" = 2),
      header_above2 = c(" "=2,
                        "Baseline SARS-CoV-2 Negative Vaccine Recipients" = 8),
      col1="1cm"),
    
    case_plcb_neg = list(
      table_header = "Antibody levels in the baseline SARS-CoV-2 negative
      per-protocol cohort (placebo recipients)",
      table_footer =c(
        paste(paste(sprintf("Cases for Day %s markers are baseline negative per-protocol placebo recipients 
      with the symptomatic infection COVID-19 primary endpoint diagnosed starting %s day(s) 
      after the Day %s study visit.", config.cor$tpeak, config.cor$tpeaklag, config.cor$tpeak), collapse=" "),
              "Non-cases/Controls are baseline negative per-protocol placebo recipients sampled into the random subcohort 
      with no COVID-19 endpoint diagnosis by the time of data-cut."),
        "N is the number of cases sampled into the subcohort within baseline covariate strata.",
        "The denominator in Resp Rate is the number of participants in the whole per-protocol cohort within baseline 
      covariate strata, calculated using inverse probability weighting."),
      
      col_name = c("Visit", "Marker", "N", "Resp rate", "GMT/GMC", "N",
                   "Resp rate", "GMT/GMC", "Resp Rate\nDifference", "GMTR/GMCR"),
      header_above1 = c(" "=2, "Cases*" = 3, "Non-Cases/Control" = 3,
                        "Comparison" = 2),
      header_above2 = c(" "=2,
                        "Baseline SARS-CoV-2 Negative Placebo Recipients" = 8),
      col1="1cm"),
    
    case_vacc_pos = list(
      table_header = "Antibody levels in the baseline SARS-CoV-2 positive
      per-protocol cohort (vaccine recipients)",
      table_footer =c(
        paste(paste(sprintf("Cases for Day %s markers are baseline positive per-protocol vaccine recipients 
      with the symptomatic infection COVID-19 primary endpoint diagnosed starting %s day(s) 
      after the Day %s study visit.", config.cor$tpeak, config.cor$tpeaklag, config.cor$tpeak), collapse=" "),
              "Non-cases/Controls are baseline positive per-protocol vaccine recipients sampled into the random subcohort 
      with no COVID-19 endpoint diagnosis by the time of data-cut."),
        "N is the number of cases sampled into the subcohort within baseline covariate strata.",
        "The denominator in Resp Rate is the number of participants in the whole per-protocol cohort within baseline 
      covariate strata, calculated using inverse probability weighting."),
      
      col_name = c("Visit", "Marker", "N", "Resp rate", "GMT/GMC", "N",
                   "Resp rate", "GMT/GMC", "Resp Rate\nDifference", "GMTR/GMCR"),
      header_above1 = c(" "=2, "Cases*" = 3, "Non-Cases/Control" = 3,
                        "Comparison" = 2),
      header_above2 = c(" "=2,
                        "Baseline SARS-CoV-2 Positive Vaccine Recipients" = 8),
      col1="1cm"),
    
    case_plcb_pos = list(
      table_header = "Antibody levels in the baseline SARS-CoV-2 positive
      per-protocol cohort (placebo recipients)",
      table_footer =c(
        paste(paste(sprintf("Cases for Day %s markers are baseline positive per-protocol placebo recipients 
      with the symptomatic infection COVID-19 primary endpoint diagnosed starting %s day(s) 
      after the Day %s study visit.", config.cor$tpeak, config.cor$tpeaklag, config.cor$tpeak), collapse=" "),
              "Non-cases/Controls are baseline positive per-protocol placebo recipients sampled into the random subcohort 
      with no COVID-19 endpoint diagnosis by the time of data-cut."),
        "N is the number of cases sampled into the subcohort within baseline covariate strata.",
        "The denominator in Resp Rate is the number of participants in the whole per-protocol cohort within baseline 
      covariate strata, calculated using inverse probability weighting."),
      
      col_name = c("Visit", "Marker", "N", "Resp rate", "GMT/GMC", "N",
                   "Resp rate", "GMT/GMC", "Resp Rate\nDifference", "GMTR/GMCR"),
      header_above1 = c(" "=2, "Cases*" = 3, "Non-Cases/Control" = 3,
                        "Comparison" = 2),
      header_above2 = c(" "=2,
                        "Baseline SARS-CoV-2 Positive Placebo Recipients" = 8),
      col1="1cm"))


    
# cutoff.name <- config$llox_label

timepoints <- config$timepoints

labels.age <- case_when(study_name %in% c("ENSEMBLE", "MockENSEMBLE") ~ c("Age 18 - 59", "Age $\\geq$ 60"), 
                        TRUE~ c("Age $<$ 65", "Age $\\geq$ 65"))

labels.minor <- c("Communities of Color", "White Non-Hispanic")

labels.BMI <- c("Underweight BMI < 18.5", "Normal 18.5 $\\leq$ BMI < 25", 
                "Overweight 25 $\\leq$ BMI < 30", "Obese BMI $\\geq$ 30")

labels.time <- labels.time[times]

if ("BbindN" %in% names(dat.mock) & any(grepl("bind", assays))) assays <- union(assays, "bindN")


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


# dat.mock was made in _common.R
dat <- dat.mock


# The stratified random cohort for immunogenicity
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
    
    demo.stratum.ordered=case_when(study_name=="VAT08m" & (max(demo.stratum) > length(demo.stratum.labels)) ~ demo.stratum-as.numeric(demo.stratum>4),
                                   !is.na(demo.stratum) ~ as.numeric(demo.stratum), 
                                   age.geq.65 == 1 ~ 7, 
                                   age.geq.65 == 0 & HighRiskInd==1 ~ 8,
                                   age.geq.65 == 0 & HighRiskInd==0 ~ 9), 
    
    AgeRisk1 = ifelse(AgeC==labels.age[1], AgeRiskC, NA),
    AgeRisk2 = ifelse(AgeC==labels.age[2], AgeRiskC, NA),
    All = "All participants"
    )

if(study_name %in% c("ENSEMBLE", "MockENSEMBLE")){
  ds_s <- ds_s %>% 
    mutate(CountryC=labels.countries.ENSEMBLE[Country+1],
           RegionC=labels.regions.ENSEMBLE[Region+1],
           URMC = case_when(URMforsubcohortsampling == 1 & Country ==0 ~ "Communities of Color",
                            URMforsubcohortsampling == 0 & Country ==0 ~ "White Non-Hispanic", 
                            TRUE ~ as.character(NA)),
           AgeURM = case_when(is.na(URMC) ~ as.character(NA), 
                              TRUE ~ paste(AgeC, URMC)),
           demo.stratum.ordered=demo.stratum,
           HIVC = c("Positive", "Negative")[2-HIVinfection],
           BMI = case_when(max(BMI, na.rm=T) < 5 ~ labels.BMI[BMI],
                           BMI>=30 ~ "Obese BMI $\\geq$ 30", 
                           BMI>=25 ~ "Overweight 25 $\\leq$ BMI < 30",
                           BMI>=18.5 ~ "Normal 18.5 $\\leq$ BMI < 25",
                           BMI<18.5 ~ "Underweight BMI < 18.5")
           )
}

if(study_name %in% c("AZD1222")){
  ds_s <- ds_s %>% 
    mutate(CountryC=c("Chile", "Peru", "United States")[Country+1])
}

# Step2: Responders
# Post baseline visits
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
tab_dm_neg <- tab_dm_pos <- tab_dm_neg_ph1 <- tab_dm_pos_ph1 <- NULL
tab_strtm1 <- tab_strtm2 <- tab_strtm2_1 <- tab_strtm2_2 <- tab_case_cnt <- NULL
rpcnt_case <- rgm_case <- rgmt_case <- NULL
case_vacc_neg <-case_plcb_neg <-case_vacc_pos <-case_plcb_pos <- NULL

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
  num_v1 <- c("Age") # Summaries - Mean &ßß Range
  num_v2 <- NULL # Summaries - Mean & St.d
  cat_v <- c("AgeC", "SexC", "CountryC", "raceC", "ethnicityC", "HighRiskC", "AgeRiskC")
} else{ # Keeping the minimal
  num_v1 <- c("Age") # Summaries - Mean & Range
  num_v2 <- NULL # Summaries - Mean & St.d
  cat_v <- c("AgeC", "SexC", "raceC", "ethnicityC", "HighRiskC", "AgeRiskC")
} 

ds_long_ttl <- ds %>%
  dplyr::filter(ph2.immuno) %>%
  # dplyr::filter(!!as.name(paste0("ph2.D", tpeak))) %>% 
  bind_rows(mutate(., Arm="Total")) %>% 
  mutate(AgeRiskC = ifelse(grepl("$\\geq$ 65", AgeRiskC, fixed=T), "Age $\\geq$ 65 ", AgeRiskC)) %>% 
  mutate_all(as.character) %>% 
  pivot_longer(all_of(c(num_v1, num_v2, cat_v)), names_to="subgroup", values_to="subgroup_cat")

ds_long_ttl_ph1 <- ds %>%
  dplyr::filter(Perprotocol & SubcohortInd==1) %>%
  # dplyr::filter(!!as.name(paste0("ph2.D", tpeak))) %>% 
  bind_rows(mutate(., Arm="Total")) %>% 
  mutate(AgeRiskC = ifelse(grepl("$\\geq$ 65", AgeRiskC, fixed=T), "Age $\\geq$ 65 ", AgeRiskC)) %>% 
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
  pivot_wider(c(`Baseline SARS-CoV-2`, Arm, subgroup, subgroup_cat, rslt),
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
  pivot_wider(c(`Baseline SARS-CoV-2`, Arm, subgroup, subgroup_cat, rslt),
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
if (study_name %in% c("COVE", "MockCOVE", "MockENSEMBLE", "PREVENT19", "VAT08m")){
  nonCaseD <- timepoints[length(timepoints)]
} else {
  nonCaseD <- tpeak
}

ds <- ds %>% 
  mutate(
    EventIndPrimaryD1 = ifelse(study_name=="VAT08m" & grepl("omi", COR), EventIndOmicronD1, EventIndPrimaryD1),
    Case = case_when(Perprotocol==1 & 
                            !!as.name(config.cor$Earlyendpoint)==0 & 
                            !!as.name(paste0("TwophasesampIndD", config.cor$tpeak))==1 & 
                            !!as.name(config.cor$EventIndPrimary)==1 ~ "Cases",
                          Perprotocol==1 & 
                            # !!as.name(ifelse(length(timepoints)>1, paste0("EarlyendpointD",timepoints[length(timepoints)]), config.cor$Earlyendpoint))==0 &
                            # AnyinfectionD1==0 & 
                            !!as.name(paste0("EarlyendpointD",nonCaseD))==0 &
                            !!as.name(paste0("TwophasesampIndD", nonCaseD))==1 & 
                            EventIndPrimaryD1==0 ~ "Non-Cases"))


# Added table: 
demo.stratum.ordered <- gsub(">=", "$\\\\geq$", demo.stratum.labels, fixed=T)

if (study_name %in% c("COVE", "MockCOVE")){
  demo.stratum.ordered <- gsub("URM", "Minority", demo.stratum.ordered)
  demo.stratum.ordered <- gsub("White non-Hisp", "Non-Minority", demo.stratum.ordered)
  demo.stratum.ordered[7:9] <- c("Age $\\\\geq$ 65, Unknown", "Age < 65, At risk, Unknown", "Age < 65, Not at risk, Unknown")
} else if (study_name %in% c("ENSEMBLE", "MockENSEMBLE")) {
  demo.stratum.ordered <- gsub("URM", "Communities of Color", demo.stratum.ordered)
  demo.stratum.ordered <- gsub("At risk", "Presence of comorbidities", demo.stratum.ordered)
  demo.stratum.ordered <- gsub("Not at risk", "Absence of comorbidities", demo.stratum.ordered)
}

strtm_cutoff <- ifelse(study_name %in% c("ENSEMBLE", "MockENSEMBLE"), length(demo.stratum.ordered)/2, length(demo.stratum.ordered))

tab_strtm <- ds %>% 
  filter(!!as.name(config.cor$ph2)) %>% 
  group_by(demo.stratum.ordered, Arm, `Baseline SARS-CoV-2`) %>%
  summarise("Day {tpeak} Cases":=sum(Case=="Cases", na.rm=T), 
            `Non-Cases`=sum(Case=="Non-Cases", na.rm=T)) %>% 
  pivot_longer(cols=c(!!as.name(paste("Day", tpeak, "Cases")), `Non-Cases`)) %>% 
  arrange(`Baseline SARS-CoV-2`, demo.stratum.ordered) %>% 
  pivot_wider(id_cols=c(Arm, name), 
              names_from = c(`Baseline SARS-CoV-2`, demo.stratum.ordered), 
              values_from=value) 


tab_strtm1 <- tab_strtm %>% select(Arm, name, any_of(paste0("Negative_", 1:strtm_cutoff)), 
                                   any_of(paste0("Positive_", 1:strtm_cutoff)))
tab_strtm2 <- tab_strtm %>% select(Arm, name, any_of(paste0("Negative_", (strtm_cutoff+1):(strtm_cutoff*2))), 
                                   any_of(paste0("Positive_", (strtm_cutoff+1):(strtm_cutoff*2))))

ls_strtm <- list(tab_strtm1, tab_strtm2)

for (i in 1:2){
  if ((n_strtm.i <- ceiling(ncol(ls_strtm[[i]])/2-1))!=0) {
  tlf[[paste0("tab_strtm", i)]]$col_name <- colnames(ls_strtm[[i]])[-1] %>%
    gsub("name", " ", .) %>% 
    gsub("Negative_", "", .) %>% 
    gsub("Positive_", "", .) 
  
  ds.i <- filter(ds, demo.stratum.ordered %in% ((i-1)*strtm_cutoff+1):(i*strtm_cutoff))
  
  tlf[[paste0("tab_strtm", i)]]$table_header <- 
    sprintf("Sample Sizes of Random Subcohort Strata (with antibody markers data at D%s) Plus All Other Cases Outside the Random Subcohort %s",
            config.cor$tpeak, ifelse(is.null(ds.i$RegionC), "", paste("in", paste(sort(unique(ds.i$RegionC))))))
  
  tlf[[paste0("tab_strtm", i)]]$header_above1 <- c(" "=1, "Baseline SARS-CoV-2 Negative" = sum(grepl("Negative", colnames(ls_strtm[[i]]))), 
                                    "Baseline SARS-CoV-2 Positive" = sum(grepl("Positive", colnames(ls_strtm[[i]]))))
  
  tlf[[paste0("tab_strtm", i)]]$header_above1 <- tlf[[paste0("tab_strtm", i)]]$header_above1[tlf[[paste0("tab_strtm", i)]]$header_above1!=0]
  
  tab_strtm_header2 <- ncol(ls_strtm[[i]])-1
  names(tab_strtm_header2) <- sprintf("%s\nSample Sizes (N=%s Participants) (%s Trial)", 
                                      tlf[[paste0("tab_strtm", i)]]$table_header,
                                      sum(ds[ds$demo.stratum.ordered%in%1:strtm_cutoff, paste0("ph2.D", tpeak)]), 
                                      stringr::str_to_title(study_name))
  tlf[[paste0("tab_strtm", i)]]$header_above2 <- tab_strtm_header2
  tlf[[paste0("tab_strtm", i)]]$table_footer <- c("Demographic covariate strata:",
                                   paste(sort(unique(ds.i$demo.stratum.ordered)), 
                                         demo.stratum.ordered[sort(unique(ds.i$demo.stratum.ordered))], 
                                         sep=". "),
                                   " ",
                                   "Minority includes Blacks or African Americans, Hispanics or Latinos, American Indians or
                   Alaska Natives, Native Hawaiians, and other Pacific Islanders."[study_name %in% c("COVE", "MockCOVE")],
                                   "Non-Minority includes all other races with observed race (Asian, Multiracial, White, Other) and observed ethnicity Not Hispanic or Latino.
                   Participants not classifiable as Minority or Non-Minority because of unknown, unreported or missing were not included."[study_name %in% c("COVE", "MockCOVE")],
                                   " "[study_name %in% c("COVE", "MockCOVE")],
                                   "Observed = Numbers of participants sampled into the subcohort within baseline covariate strata.",
                                   "Estimated = Estimated numbers of participants in the whole per-protocol cohort within baseline 
  covariate strata, calculated using inverse probability weighting.")
  } 
}

if (ncol(tab_strtm1)==2) tab_strtm1 <- NULL
if (ncol(tab_strtm2)==2) tab_strtm2 <- NULL


# median (interquartile range) days from vaccination to the tpeak visit

if ((Numberdays <- paste0("NumberdaysD1toD", config.cor$tpeak)) %in% names(ds)) {
  tab_days <- ds %>% 
    filter(!!as.name(config.cor$ph2), !is.na(Case)) %>% 
    mutate(Visit = paste("Day", config.cor$tpeak)) %>% 
    bind_rows(., mutate(., Case = "Total")) %>% 
    group_by(Visit, Arm, Case) %>% 
    summarise(N=n(), dmed=median(!!as.name(Numberdays), na.rm=T), iqr=IQR(!!as.name(Numberdays))) %>% 
    select(Visit, Arm, ` `=Case, 
           N, `Median\n(days)`=dmed, `Interquatile Range\n(days)`=iqr)
} else {
  tab_days <- NULL
}

# Case counts by availability of markers at baseline, d29, d57

if (study_name %in% c("COVE", "MockCOVE")){
  tab_case_cnt <- make.case.count.marker.availability.table(dat) %>% 
    data.frame(check.names = F) %>% 
    rename_all(gsub, pattern=".", replacement="_", fixed=T) %>% 
    rownames_to_column("Case") %>% 
    pivot_longer(cols = !Case,
                 names_to = c(".value", "Arm"),
                 names_pattern = "(.*)_(.*)") %>% 
    mutate(Arm = factor(ifelse(Arm=="vacc", "Vaccine", "Placebo"), levels=c("Vaccine", "Placebo"))) %>%
    arrange(Arm, Case) %>% 
    rename_at(-c(1:2), function(x)paste0("$",x,"$"))
} 

# Generate a full table with all the estimates: response rates, GMT, GMTR, etc.
# (Per Peter's email Feb 5, 2021)
# Cases vs Non-cases

  sub.by <- c("Arm", "`Baseline SARS-CoV-2`")
  ds.i <- filter(ds, !!as.name(config.cor$ph1))
  resp.v <- intersect(grep("Resp", names(ds), value = T), 
                      grep(config.cor$tpeak, names(ds), value = T))
  gm.v <- intersect(assays_col, grep(config.cor$tpeak, names(ds), value = T))
  
  subs <- "Case"
  comp.i <- c("Cases", "Non-Cases")
  
  rpcnt_case <- get_rr(ds.i, resp.v, subs, sub.by, strata=config.cor$WtStratum, weights=config.cor$wt, subset=config.cor$ph2) 
  rgm_case <- get_gm(ds.i, gm.v, subs, sub.by, strata=config.cor$WtStratum, weights=config.cor$wt, subset=config.cor$ph2) 
  rgmt_case <- get_rgmt(ds.i, gm.v, subs, comp_lev=comp.i, sub.by, strata=config.cor$WtStratum, weights=config.cor$wt, subset=config.cor$ph2) 
  
  print("Done with table 2b & 3b") 


rrdiff_case <- rpcnt_case %>% 
  # dplyr::filter(subgroup %in% subs & grepl("Resp",resp_cat)) %>% 
  mutate(groupn = 2-match(Group, comp.i)%%2) %>%
  pivot_wider(id_cols = c(subgroup, `Baseline SARS-CoV-2`, Arm, Visit, Marker, Ind),
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
  
print("Done with table6")

tab_case <- full_join(rpcnt_case, rgm_case,
                      by = c("Group", "Arm", "Baseline SARS-CoV-2", 
                             "N", "Marker", "Visit")) %>% 
  pivot_wider(id_cols = c(Arm, `Baseline SARS-CoV-2`, Marker, Visit),
              names_from = Group, 
              values_from = c(N, rslt, `GMT/GMC`)) %>% 
  full_join(rrdiff_case, by = c("Arm", "Baseline SARS-CoV-2", "Marker", "Visit")) %>% 
  full_join(rgmt_case, by = c("Arm", "Baseline SARS-CoV-2", "Marker", "Visit"))

if(length(comp_NA <- setdiff(comp.i, rpcnt_case$Group))!=0){
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
      mutate(`Ratios of GMT/GMC`=replace_na(`Ratios of GMT/GMC`, "-"))
  }

tab_case <- tab_case %>% 
  select(Arm, `Baseline SARS-CoV-2`, Visit, Marker, `N_Cases`, `rslt_Cases`, 
         `GMT/GMC_Cases`, `N_Non-Cases`, `rslt_Non-Cases`, `GMT/GMC_Non-Cases`,  
         rrdiff, `Ratios of GMT/GMC`) %>% 
  arrange(Arm, `Baseline SARS-CoV-2`, Visit) 
  
case_vacc_neg <- tab_case %>% 
  dplyr::filter(Arm == "Vaccine" & `Baseline SARS-CoV-2` == "Negative") %>% 
  select(-c(Arm, `Baseline SARS-CoV-2`))

case_plcb_neg <- tab_case %>% 
  dplyr::filter(Arm == "Placebo" & `Baseline SARS-CoV-2` == "Negative") %>% 
  select(-c(Arm, `Baseline SARS-CoV-2`))

case_vacc_pos <- tab_case %>% 
  dplyr::filter(Arm == "Vaccine" & `Baseline SARS-CoV-2` == "Positive") %>% 
  select(-c(Arm, `Baseline SARS-CoV-2`))

case_plcb_pos <- tab_case %>% 
  dplyr::filter(Arm == "Placebo" & `Baseline SARS-CoV-2` == "Positive") %>% 
  select(-c(Arm, `Baseline SARS-CoV-2`))

print("Done with all tables") 

if(study_name %in% c("PREVENT19") & all(ds$Country==0)){
  for (i in 1:length(tlf)){
    if(!is.null(tlf[[i]]$table_header)){
      tlf[[i]]$table_header <- paste0(tlf[[i]]$table_header, " in U.S. only")
    }
  }
}

if(tpeak!=timepoints[1]){
  tlf <- tlf[!names(tlf) %in% c("tab_dm_neg_ph1", "tab_dm_pos_ph1", "tab_dm_neg", "tab_dm_pos")]
}


# path for tables
save.results.to <- here::here("output")
if (!dir.exists(save.results.to))  dir.create(save.results.to)
save.results.to <- paste0(here::here("output"), "/", attr(config,"config"))

if (!dir.exists(save.results.to))  dir.create(save.results.to)
print(paste0("save.results.to equals ", save.results.to))

save(tlf, tab_dm_neg, tab_dm_pos, tab_dm_neg_ph1, tab_dm_pos_ph1, 
     tab_strtm1, tab_strtm2, tab_strtm2_1, tab_strtm2_2, 
     tab_case_cnt, tab_days, 
     case_vacc_neg, case_plcb_neg,
     case_vacc_pos, case_plcb_pos,
     file = file.path(save.results.to, sprintf("Tables%s.Rdata", ifelse(exists("COR"), COR, ""))))
