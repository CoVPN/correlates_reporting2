##################################################
# obligatory to append to the top of each script #
renv::activate(project = here::here("..")) #
source(here::here("..", "_common.R")) #
##################################################

source(here::here("code", "make_functions.R"))
library(survey)
library(tidyverse)
library(dplyr, warn.conflicts = FALSE)
# Suppress summarise info
options(dplyr.summarise.inform = FALSE)
# For stratum with 1 ppt
options(survey.lonely.psu="adjust")
# To select which tables are included in the report.
# Also to modify the headers, footers, etc. for each table

# The stratified random cohort for immunogenicity
if (grepl("prevent19", COR) & grepl("stage2", COR)){
  immuno_timepoint <- c("D35", "C1")
} 

randomsubcohort <- case_when(study_name=="COVE" ~ "This table summarizes the 
      random subcohort, which was randomly sampled from the per-protocol cohort. The 
      sampling was stratified by 24 strata defined by enrollment characteristics: Assigned 
      treatment arm $\\\\times$ Baseline SARS-CoV-2 naïve vs. non-naïve status 
      (defined by serostatus and NAAT testing) $\\\\times$ Randomization strata 
      (Age < 65 and at-risk, Age < 65 and not at-risk, Age $\\\\geq 65)\\\\times$ 
      Communities of color (Yes/No) defined by White Non-Hispanic vs. all 
      others (following the primary COVE trial paper).",
                             
                             study_name=="ENSEMBLE" ~ "This table summarizes characteristics of 
      per-protocol participants in the immunogenicity subcohort, which was randomly 
      sampled from the study cohort. The sampling was The sampling was stratified by 
      strata defined by enrollment characteristics: Assigned randomization arm $\\\\times$ 
      Baseline SARS-CoV-2 seronegative vs. seropositive $\\\\times$ Randomization strata. 
      The U.S. subcohort includes 8 baseline demographic strata; the Latin America 
      and South Africa subcohorts each include 4 baseline demographic strata.",
                             
                             study_name=="PREVENT19" ~ "This table summarizes characteristics of 
      per-protocol participants in the immunogenicity subcohort, which was randomly 
      sampled from the study cohort. The sampling was The sampling was stratified by 
      strata defined by enrollment characteristics: Assigned randomization arm $\\\\times$ 
      Baseline SARS-CoV-2 seronegative vs. seropositive $\\\\times$ Randomization strata. 
      The U.S. subcohort includes 8 baseline demographic strata; the Mexico subcohort includes 2 baseline demographic strata.",
                             
                             TRUE~ "This table summarizes characteristics of 
      per-protocol participants in the immunogenicity subcohort, which was randomly 
      sampled from the study cohort.")

tlf <-
  list(
    tab_dm_neg = list(
      table_header = "Demographic and Clinical Characteristics at Baseline in 
      the Baseline SARS-CoV-2 Negative Per-Protocol Cohort",
      table_footer = randomsubcohort,
      deselect = c("subgroup", "Timepoint"),
      pack_row = "subgroup",
      loop = "Timepoint",
      col1="7cm"
    ),
    
    tab_dm_pos = list(
      table_header = "Demographic and Clinical Characteristics at Baseline in 
      the Baseline SARS-CoV-2 Positive Per-Protocol Cohort",
      table_footer = randomsubcohort,
      deselect = c("subgroup", "Timepoint"),
      pack_row = "subgroup",
      loop = "Timepoint",
      col1="7cm"
    ),
    
    tab_strtm1 = list(
      table_header = "Sample Sizes of Random Subcohort Strata for Measuring Antibody Markers",
      loop = "Timepoint",
      deselect = c("Arm", "Timepoint"),
      pack_row = "Arm"
    ),
    
    tab_strtm2 = list(
      table_header = "Sample Sizes of Random Subcohort Strata for Measuring Antibody Markers",
      loop = "Timepoint",
      deselect =  c("Arm", "Timepoint"),
      pack_row = "Arm"
    ),
    
    tab_bind1 = list(
      table_header = "Percentage of responders, and participants
      with 2-fold rise, and participants with 4-fold rise for binding antibody
      markers",
      table_footer = c(
        sprintf("Binding Antibody Responders are defined as participants with concentration 
        above the specified positivity cut-off."),
      "Percentages are calculated for the whole per-protocol group/subgroup, 
        using inverse probability weighting."),
      loop = "subgroup",
      group_table_col = c("Rx", "Group", "Baseline", "Visit", "N", "Marker"),
      deselect = "subgroup",
      pack_row = "subgroup", 
      font_size=8
    ),
    
    tab_bind2 = list(
      table_header = "Percentage of responders, and participants
      with 2-fold rise, and participants with 4-fold rise for binding antibody
      markers",
      table_footer = c(
        sprintf("Binding Antibody Responders are defined as participants with concentration 
        above the specified positivity cut-off"),
        "Percentages are calculated for the whole per-protocol group/subgroup, 
        using inverse probability weighting."),
      loop = "subgroup",
      group_table_col = c("Rx", "Group", "Baseline", "Visit", "N", "Marker"),
      deselect = "subgroup",
      pack_row = "subgroup", 
      font_size=8),
    
    tab_pseudo = list(
      table_header = "Percentage of responders, and participants
      participants with 2-fold rise, and participants with 4-fold rise for 
      ID50 pseudo-virus neutralization antibody markers",
      table_footer = c(
        "Neutralization Responders are defined as participants who had baseline
        values below the lower limit of detection (LLOQ) with detectable
        ID50 neutralization titer above the assay LLOQ, or as participants with
        baseline values above the LLOQ with a 4-fold increase in ID50.",
        "Percentages are calculated for the whole per-protocol group/subgroup, 
        using inverse probability weighting."
      ),
      loop = "subgroup",
      group_table_col = c("Rx", "Group", "Baseline", "Visit", "N", "Marker"),
      deselect = "subgroup",
      pack_row = "subgroup", 
      font_size=8),
    
    tab_wt = list(
      table_header = "Percentage of responders, and participants
      participants with 2-fold rise, and participants with 4-fold rise
      for MN50 WT live virus neutralization antibody markers",
      table_footer = c(
        "Neutralization Responders are defined as participants who had baseline
        values below the lower limit of detection (LLOQ) with detectable
        ID50 neutralization titer above the assay LLOQ, or as participants with
        baseline values above the LLOQ with a 4-fold increase in ID50.",
        "Percentages are calculated for the whole per-protocol group/subgroup, 
        using inverse probability weighting."
      ),
      loop = "subgroup",
      group_table_col = c("Rx", "Group", "Baseline", "Visit", "N", "Marker"),
      deselect = "subgroup",
      pack_row = "subgroup", 
      font_size=8),
    
    tab_gm = list(
      table_header = "Geometric mean titers (GMTs) and geometric mean
      concentrations (GMCs)",
      table_footer = "",
      loop = "subgroup",
      group_table_col = c("Rx", "Group", "Baseline", "Visit", "N", "Marker"),
      deselect = "subgroup",
      pack_row = "subgroup",
      col1="7cm", 
      font_size=8),
    
    tab_rgmt = list(
      table_header = "The ratios of GMTs/GMCs between groups",
      table_footer = " ",
      loop = "subgroup",
      pack_row = "subgroup",
      deselect = "subgroup",
      group_table_col = c("subgroup","Rx", "Baseline", "Visit"),
      col1="4cm", 
      font_size=8),
    
    tab_rrdiff = list(
      table_header = "Differences in the responder rates, 2FRs, 4FRs between 
      the groups",
      table_footer = "Percentages are calculated for the whole per-protocol 
      group/subgroup, using inverse probability weighting.",
      loop = "subgroup",
      pack_row = "subgroup",
      group_table_col = c( "Group", "Baseline","Visit", "Marker"),
      deselect = "subgroup",
      col1="4cm", 
      font_size=8))
   

# Depends on the Incoming data
if(include_bindN & !"bindN" %in% assays & study_name!="PROFISCOV"){
  assays <- sort(c("bindN", assays))
}

labels.age <- case_when(study_name %in% c("ENSEMBLE", "MockENSEMBLE") ~ c("Age 18 - 59", "Age $\\geq$ 60"), 
                        TRUE~ c("Age $<$ 65", "Age $\\geq$ 65"))

labels.minor <- c("Communities of Color", "White Non-Hispanic")

labels.BMI <- c("Underweight BMI < 18.5", "Normal 18.5 $\\leq$ BMI < 25", 
                "Overweight 25 $\\leq$ BMI < 30", "Obese BMI $\\geq$ 30")

labels.time <- labels.time[times]
# hacky fix
labels.assays.short <- labels.assays.short.tabular[assays]

# redefines what is in _common.R to use shorter names
labels.assays.long <- data.frame (purrr::imap_dfc(labels.assays.short, ~ paste0(labels.assays.short[.y], ": ", labels.time)))
rownames(labels.assays.long) <- names(labels.time)

assays_col <- as.vector(outer(gsub("D", "Day", immuno_timepoint), assays, paste0))

labels.assays <- expand.grid(
  time = gsub("D", "Day", immuno_timepoint),
  marker = assays,
  stringsAsFactors = FALSE
) %>%
  left_join(assay_metadata, by=c("marker"="assay")) %>% 
  mutate(label.short = assay_label_short, 
         Marker = assay_label_short,
         Visit=time, 
         colname = paste0(time, marker))
 
resp.lb <- expand.grid(
  time = gsub("D", "Day", immuno_timepoint), marker = assays,
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


# Immunogenicity Tables

# Read in original data
dat <- dat_proc
pos.cutoffs <- assay_metadata$pos.cutoff
names(pos.cutoffs) <- assay_metadata$assay

# The stratified random cohort for immunogenicity

# ph1.immuno.D35 is same as ph1.immuno.C1
ds_s <- dat %>%
    dplyr::filter(!!as.name(paste0("ph1.immuno.", immuno_timepoint))) %>%
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
      AgeC = ifelse(is.na(Senior), ifelse(age.geq.65 == 1, labels.age[2], labels.age[1]), ifelse(Senior == 1, labels.age[2], labels.age[1])),
      SexC = ifelse(Sex == 1, "Female", "Male"),
      AgeRiskC = paste(AgeC, HighRiskC),
      AgeSexC = paste(AgeC, SexC),
      AgeMinorC = ifelse(is.na(MinorityC), NA, paste(AgeC, MinorityC)),
      `Baseline SARS-CoV-2` = factor(ifelse(Bserostatus == 1, "Positive", "Negative"),
                                     levels = c("Negative", "Positive")
      ),
      Arm = factor(ifelse(Trt == 1, "Vaccine", "Placebo"), 
                   levels = c("Vaccine", "Placebo")),
      demo.stratum.ordered=case_when(!is.na(demo.stratum) ~ as.numeric(demo.stratum), 
                                     age.geq.65 == 1 ~ 7, 
                                     age.geq.65 == 0 & HighRiskInd==1 ~ 8,
                                     age.geq.65 == 0 & HighRiskInd==0 ~ 9), 
      AgeRisk1 = ifelse(AgeC==labels.age[1], AgeRiskC, NA),
      AgeRisk2 = ifelse(AgeC==labels.age[2], AgeRiskC, NA),
      All = "All participants"
    ) 
  
  if(study_name %in% c("ENSEMBLE", "MockENSEMBLE", "PREVENT19")){
    ds_s <- ds_s %>% 
      mutate(CountryC = labels.countries.ENSEMBLE[Country+1],
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
  
  if(study_name %in% c("ENSEMBLE", "MockENSEMBLE")){
    ds_s <- ds_s %>% 
      mutate(RegionC = labels.regions.ENSEMBLE[Region+1])
  } 
  
  if(study_name %in% c("PROFISCOV")){
    ds_s <- ds_s %>% 
      mutate(URMC = case_when(URMforsubcohortsampling == 1 ~ "Communities of Color",
                              URMforsubcohortsampling == 0 ~ "White Non-Hispanic", 
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
  # Step2: Responders, % >=2FR, % >=4FR, % >=2lloq, % >=4lloq
  # Post baseline visits
  
  ds <- getResponder(ds_s, times=gsub("D", "Day", immuno_timepoint), 
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
  MinorityC = "Underrepresented Minority Status",
  AgeMinorC = "Age, Communities of color",
  URMC = "Underrepresented Minority Status in the U.S.",
  AgeURM = "Age, Underrepresented Minority Status in the U.S.",
  CountryC = "Country",
  HIVC = "HIV Infection"
)

grplev <- c("", labels.age, "At-risk", "Not at-risk", 
            paste(labels.age[1], c("At-risk", "Not at-risk")),
            paste(labels.age[2], c("At-risk", "Not at-risk")),
            "Male", "Female", 
            paste(labels.age[1], c("Female", "Male")),
            paste(labels.age[2], c("Female", "Male")),
            "Hispanic or Latino", "Not Hispanic or Latino", "Not reported and unknown ", 
            "White Non-Hispanic ", "Black or African American", "Asian", 
            "American Indian or Alaska Native", 
            "Native Hawaiian or Other Pacific Islander", 
            "Multiracial", "Other", "Not reported and unknown",  
            labels.minor, 
            paste(labels.age[1], labels.minor),  
            paste(labels.age[2], labels.minor),
            labels.countries.ENSEMBLE,
            "Negative", "Positive")

names(grplev) <- c("All participants", grplev[-1])



# Immunogenicity Tables

###################################################
#             Generating the Tables               #
###################################################

### Table 1. Demographics 
# Output: tab_dm
# Select the covariates to be summarised.
# num_v are columns from ds_long;
# cat_v are rows of `subgroup`

for (i in names(tlf)){
  assign(i, NULL)
}

if (study_name=="COVE") {
  num_v1 <- c("Age") # Summaries - Mean & Range
  num_v2 <- c("BMI") # Summaries - Mean & St.d
  cat_v <- c("AgeC", "SexC", "raceC", "ethnicityC", "HighRiskC", "AgeRiskC", "MinorityC")
} else { #if (study_name %in% c("ENSEMBLE", "PREVENT19")) {
  num_v1 <- c("Age") # Summaries - Mean & Range
  num_v2 <- NULL # Summaries - Mean & St.d
  cat_v <- c("AgeC", "SexC", "raceC", "ethnicityC", 
             "HighRiskC", "AgeRiskC", "URMC",  "CountryC", "HIVC", "BMI")
  
  if (study_name %in% c("PROFISCOV")) {
    cat_v <- c("AgeC", "SexC", "raceC", "ethnicityC", 
               "HighRiskC", "AgeRiskC", "URMC", "HIVC", "BMI")
  }
} 

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
              paste(labels.age[2], ""), "URM", "Non-URM", labels.countries.ENSEMBLE,
              "Negative", "Positive", labels.BMI)



tab_dm <- lapply(immuno_timepoint, function(x){
  ds_long_ttl <- ds %>%
  dplyr::filter(!!as.name(paste0("ph2.immuno.", x))) %>% 
  # bind_rows(mutate(., Arm ="Total")) %>%
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
  dplyr::filter(as.character(subgroup) %in% cat_v)


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
  pivot_wider(id_cols = c(`Baseline SARS-CoV-2`,subgroup, subgroup_cat),
              names_from = Arm, 
              names_sort = T,
              values_from = c(rslt)) %>%
  mutate(Characteristics = factor(subgroup_cat, levels=char_lev),
         subgroup=factor(subgroup, levels=subgrp), Timepoint=x) %>%
  arrange(`Baseline SARS-CoV-2`, subgroup, Characteristics)


return(tab_dm)
})

tab_dm <- bind_rows(tab_dm)

tab_dm_pos <- tab_dm %>% 
  dplyr::filter(`Baseline SARS-CoV-2` == "Positive") %>% 
  select_at(c("subgroup", "Characteristics", "Timepoint",
                         grep("Vaccine" ,names(.), value = T),
                         grep("Placebo" ,names(.), value = T),
                         grep("Total" ,names(.), value = T))) %>% 
  select_if(~ !all(is.na(.)))

tab_dm_neg <- tab_dm %>% 
  dplyr::filter(`Baseline SARS-CoV-2` == "Negative") %>% 
  select_at(c("subgroup", "Characteristics", "Timepoint",
              grep("Vaccine" ,names(.), value = T),
              grep("Placebo" ,names(.), value = T),
              grep("Total" ,names(.), value = T))) %>% 
  select_if(~ !all(is.na(.))) 


print("Done with table 1") 

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

tab_strtm <- lapply(immuno_timepoint, function(x){
  
  tab_strtm <- ds %>% 
    dplyr::filter(!!as.name(paste0("ph2.immuno.", x))) %>% 
    group_by(demo.stratum.ordered, !!as.name(paste0("wt.immuno.", x)), Arm, `Baseline SARS-CoV-2`) %>% 
    summarise(Observed=n(), Estimated=round(sum(!!as.name(paste0("wt.immuno.", x))), 0)) %>% 
    pivot_longer(cols=c(Observed, Estimated)) %>% 
    arrange(`Baseline SARS-CoV-2`, demo.stratum.ordered) %>% 
    pivot_wider(id_cols=c(Arm, name), 
                names_from = c(`Baseline SARS-CoV-2`, demo.stratum.ordered), 
                values_from=value) %>% 
    mutate(Timepoint=x)


  tab_strtm1 <- tab_strtm %>% select(Arm, Timepoint, name, any_of(paste0("Negative_", 1:strtm_cutoff)), 
                                     any_of(paste0("Positive_", 1:strtm_cutoff)))
  
  tab_strtm2 <- tab_strtm %>% select(Arm, Timepoint, name, any_of(paste0("Negative_", (strtm_cutoff+1):(strtm_cutoff*2))), 
                                     any_of(paste0("Positive_", (strtm_cutoff+1):(strtm_cutoff*2))))
  
  ls_strtm <- list(tab_strtm1, tab_strtm2)

  for (i in 1:2){
    if ((n_strtm.i <- ceiling(ncol(ls_strtm[[i]])/2-1))!=0) {
      tlf[[paste0("tab_strtm", i)]]$col_name <- colnames(ls_strtm[[i]])[-(1:length(tlf[[paste0("tab_strtm", i)]]$deselect))] %>%
        gsub("name", " ", .) %>% 
        gsub("Negative_", "", .) %>% 
        gsub("Positive_", "", .) 
      
      ds.i <- filter(ds, demo.stratum.ordered %in% ((i-1)*strtm_cutoff+1):(i*strtm_cutoff))
      
      tlf[[paste0("tab_strtm", i)]]$header_above1 <- c(" "=1, "Baseline SARS-CoV-2 Negative" = sum(grepl("Negative", colnames(ls_strtm[[i]]))), 
                                                       "Baseline SARS-CoV-2 Positive" = sum(grepl("Positive", colnames(ls_strtm[[i]]))))
    
      tlf[[paste0("tab_strtm", i)]]$header_above1 <- tlf[[paste0("tab_strtm", i)]]$header_above1[tlf[[paste0("tab_strtm", i)]]$header_above1!=0]
    
    
    
      tab_strtm_header2 <- ncol(ls_strtm[[i]])-length(tlf[[paste0("tab_strtm", i)]]$deselect)
      names(tab_strtm_header2) <- sprintf("%sRandom Subcohort Sample Sizes (N=%s Participants) (%s Trial)",
                                          case_when(study_name=="COVE" ~ "", 
                                                    study_name=="ENSEMBLE" ~ 
                                                      paste0(paste(c("U.S.", "Latin America", "South Africa")[sort(unique(ds.i$Region))+1], collapse=" and "), " "),
                                                    TRUE ~ ""),
                                          sum(ds.i[,paste0("ph2.immuno.", x)]), # 
                                          study_name)
      
      tlf[[paste0("tab_strtm", i)]]$header_above2 <- tab_strtm_header2
      tlf[[paste0("tab_strtm", i)]]$table_footer <- c("Demographic covariate strata:",
                                     paste(sort(unique(ds.i$demo.stratum.ordered)), 
                                           demo.stratum.ordered[sort(unique(ds.i$demo.stratum.ordered))], 
                                           sep=". "),
                                     " ",
                                     "Minority includes Blacks or African Americans, Hispanics or Latinos, American Indians or
                   Alaska Natives, Native Hawaiians, and other Pacific Islanders."[study_name=="COVE"],
                                     "Non-Minority includes all other races with observed race (Asian, Multiracial, White, Other) and observed ethnicity Not Hispanic or Latino.
                   Participants not classifiable as Minority or Non-Minority because of unknown, unreported or missing were not included."[study_name=="COVE"],
                                     " "[study_name=="COVE"],
                                     "Observed = Numbers of participants sampled into the subcohort within baseline covariate strata.",
                                     "Estimated = Estimated numbers of participants in the whole per-protocol cohort within baseline 
                                     covariate strata, calculated using inverse probability weighting.")
    }
  }
  tlf <<- tlf

  if (ncol(tab_strtm1)==3) tab_strtm1 <- NULL
  if (ncol(tab_strtm2)==3) tab_strtm2 <- NULL
  
  return(list(tab_strtm1 = tab_strtm1, tab_strtm2 = tab_strtm2, 
              tab_strtm1_header2 = tlf$tab_strtm1$header_above2, tab_strtm2_header2 = tlf$tab_strtm2$header_above2))
})


tab_strtm1 <- tab_strtm %>% 
  map("tab_strtm1") %>%
  bind_rows()


tab_strtm2 <- tab_strtm %>% 
  map("tab_strtm2") %>%
  bind_rows()

tlf$tab_strtm1$header_above2  <- tab_strtm %>% map("tab_strtm1_header2")

tlf$tab_strtm2$header_above2 <- tab_strtm %>% map("tab_strtm2_header2")


### Table 2. Responder Rates & Proportions of Magnitudes >= 2FR, 4FR
# For each binding antibody marker, the estimated percentage of participants
# defined as responders, and with concentrations >= 2x LLOQ or >=
# 4 x LLOQ, will be provided with the corresponding 95% CIs
# 
# Output: tab_bind


# Variables used for stratification in the tables
# subgroup: SAP Table 6: Baseline Subgroups
# Arm and Baseline: Assigned treatment Arms * Baseline SARS-CoV-2-19 Status
# Group: Category in each subgroup

sub.by <- c("Arm", "`Baseline SARS-CoV-2`")


if (study_name=="COVE") {
  subs <- c("All", "AgeC", "HighRiskC", "AgeRiskC", "AgeRisk1", "AgeRisk2", "SexC",
            "AgeSexC", "ethnicityC", "RaceEthC", "MinorityC", "AgeMinorC")
} else if (study_name %in% c("ENSEMBLE", "PREVENT19")) {
  subs <- c("All", "AgeC", "HIVC", "CountryC", "HighRiskC", "AgeRiskC", "AgeRisk1", "AgeRisk2", "SexC",
            "AgeSexC", "ethnicityC", "RaceEthC", "URMC", "AgeURM")
} else if (study_name %in% c("PROFISCOV")) {
  subs <- c("All", "AgeC", "HIVC", "HighRiskC", "AgeRiskC", "AgeRisk1", "AgeRisk2", "SexC",
            "AgeSexC", "ethnicityC", "RaceEthC", "URMC", "AgeURM")
}

rpcnt <- lapply(immuno_timepoint, function(x){

  resp.v <- grep(gsub("D", "Day", x), grep("Resp|2lloq|4lloq|FR2|FR4", names(ds), value = T), value = T)
  resp.v <- resp.v[sapply(resp.v, function(x)!all(is.na(ds[, x])))]
  
  get_rr(dat=ds, v=resp.v, subs=subs, sub.by=sub.by, strata="tps.stratum",
         weights=paste0("wt.immuno.", x), subset=paste0("ph2.immuno.", x))
  
  })

rpcnt <- bind_rows(rpcnt)

tab_rr <- rpcnt %>% 
  dplyr::filter(!subgroup %in% c("AgeRisk1", "AgeRisk2") & Visit != "Day 1" & Group %in% names(grplev)) %>% 
  mutate(subgroup=factor(subgrp[subgroup], levels=subgrp), Group=factor(grplev[Group], levels=grplev)) %>% 
  pivot_wider(
    id_cols = c(subgroup, Group, Arm, `Baseline SARS-CoV-2`, Marker, Visit),
    names_from = Ind, values_from = rslt) %>% 
  arrange(subgroup, Group, Visit, Arm, `Baseline SARS-CoV-2`, Marker) 

if(any(grepl("bind", assays))){
  
  if (all(c("% Greater than 2xLLOQ", "% Greater than 4xLLOQ") %in% names(tab_rr))){
    tab_bind1 <- tab_rr %>% 
      dplyr::filter(Marker %in% labels_all$Marker[grep("bind", labels_all$marker)]) %>% 
      select(subgroup, Group, Visit, Arm, `Baseline SARS-CoV-2`, Marker, Responder, 
             c("% Greater than 2xLLOQ", "% Greater than 4xLLOQ"))
  }
  
  tab_bind2 <- tab_rr %>% 
    dplyr::filter(Marker %in% labels_all$Marker[grep("bind", labels_all$marker)]) %>% 
    select(subgroup, Group, Visit, Arm, `Baseline SARS-CoV-2`, Marker, Responder, 
          `% 2-Fold Rise`, `% 4-Fold Rise`)
}

print("Done with table2") 

# Table 3 & 4. For the ID50 pseudo-virus & MN50 WT live virus neutralization 
# antibody marker, the estimated percentage of participants defined as 
# responders, participants with % 2-Fold Rise (2FR), and participants with 4-fold 
# rise (4FR) will be provided with the corresponding 95% CIs 
# 
# Output: tab_pseudo & tab_wt

if(any(grepl("pseudoneutid50", assays))){

tab_pseudo <- tab_rr %>% 
  dplyr::filter(Marker %in% labels_all$Marker[grep("pseudoneutid50", labels_all$marker)]) %>% 
  select(subgroup, Group, Visit, Arm, `Baseline SARS-CoV-2`, Marker, Responder,
         `% 2-Fold Rise`, `% 4-Fold Rise`)
}

if(any(grepl("liveneutmn50", assays))){
tab_wt <- tab_rr %>% 
  dplyr::filter(Marker %in% labels_all$Marker[grep("liveneutmn50", labels_all$marker)]) %>% 
  select(subgroup, Group, Visit, Arm, `Baseline SARS-CoV-2`, Marker, Responder,
         `% 2-Fold Rise`, `% 4-Fold Rise`)
}
print("Done with table3 & 4") 

# Table 5. Geometric mean titers (GMTs) and geometric mean concentrations (GMCs)
# will be summarized along with their 95% CIs using the t-distribution
# approximation of log-transformed concentrations/titers (for each of the 5
# Spike-targeted marker types including pseudovirus-nAb ID50 and ID80
# and WT live virus-nAb MN50, as well as for binding Ab to N).
# 


rgm <- lapply(immuno_timepoint, function(x){
  
  gm.v <- grep(gsub("D", "Day", x), assays_col, value = T)
  gm.v <- gm.v[sapply(gm.v, function(x)!all(is.na(ds[, x])))]
  
  rgm <- get_gm(dat=ds, v=gm.v, subs=subs, sub.by=sub.by, strata="tps.stratum",
              weights=paste0("wt.immuno.", x), subset=paste0("ph2.immuno.", x))

  
  return(rgm)
})

rgm <- bind_rows(rgm)

tab_gm <- rgm %>% 
  dplyr::filter(!subgroup %in% c("AgeRisk1", "AgeRisk2") & !grepl("Delta", mag_cat) & Group %in% names(grplev)) %>% 
  mutate(subgroup=factor(subgrp[subgroup], levels=subgrp), Group=factor(grplev[Group], levels=grplev)) %>% 
  arrange(subgroup, Group, Visit, Arm, `Baseline SARS-CoV-2`, Marker) %>% 
  select(subgroup, Group, Visit, Arm, `Baseline SARS-CoV-2`, Marker, N, `GMT/GMC`)

print("Done with table5") 

### Table 6. GMTRs/GMCRs will be summarized with 95% CI (t-distribution 
# approximation) for any post-baseline values compared to baseline, and
# post-Day 57 values compared to Day 57
# 
# Output: tab_gmr



### Table 7. The ratios of GMTs/GMCs will be estimated between groups with the
# two-sided 95% CIs calculated using t-distribution approximation of 
# log-transformed titers/concentrations
# Output: tab_rgmt
# 
# Ratios of GMT/GMC between subgroups among vacinees

comp_lev <- c(labels.age[2:1],
              "At-risk", "Not at-risk",
              paste(labels.age[1], c("At-risk", "Not at-risk")),
              paste(labels.age[2], c("At-risk", "Not at-risk")),
              "Male", "Female",
              "Hispanic or Latino", "Not Hispanic or Latino",
              labels.minor,
              "Positive", "Negative")

groups <- c("AgeC", "HighRiskC", "AgeRisk1", "AgeRisk2",
            "SexC", "ethnicityC",
            case_when(study_name=="COVE"~"MinorityC",
                      study_name%in%c("ENSEMBLE", "PREVENT19", "PROFISCOV")~"URMC"),
            "HIVC"[study_name %in% c("ENSEMBLE", "PROFISCOV")])


rgmt <- lapply(immuno_timepoint, function(x){
  
  mag_groups <- grep(gsub("D", "Day", x), assays_col, value = T)
  mag_groups <- mag_groups[sapply(mag_groups, function(x)!all(is.na(ds[,x])))]

  rgmt <- get_rgmt(ds, mag_groups, groups, comp_lev=comp_lev, 
                   sub.by, "tps.stratum", paste0("wt.immuno.", x), paste0("ph2.immuno.", x))

  return(rgmt)
  
  })

rgmt <- bind_rows(rgmt)

  rgmt_gm <- rgm %>% 
    dplyr::filter(!grepl("Delta", mag_cat) & Group %in% names(grplev)) %>% 
    mutate(subgroup=factor(subgrp[subgroup], levels=subgrp), Group=factor(grplev[Group], levels=grplev)) %>% 
    dplyr::filter(subgroup %in% subgrp[groups]) %>% 
    mutate(groupn = 2-match(Group, comp_lev)%%2) %>% 
    pivot_wider(id_cols = c(subgroup, Arm, `Baseline SARS-CoV-2`, Visit, Marker),
                names_from = groupn, values_from = `GMT/GMC`, 
                names_prefix = "Group")
  
  tab_rgmt <- inner_join(rgmt_gm, 
                         rgmt %>% mutate(subgroup=factor(subgrp[subgroup], levels=subgrp)),  
                         c("Baseline SARS-CoV-2", "Arm", "subgroup", "Marker", "Visit")) %>% 
    rename(`Group 1 vs 2` = comp, 
           `Group 1 GMT/GMC` = `Group1`, 
           `Group 2 GMT/GMC` = `Group2`) %>% 
    select(`Group 1 vs 2`, subgroup, Visit, Arm, `Baseline SARS-CoV-2`, Marker, 
           `Group 1 GMT/GMC`, `Group 2 GMT/GMC`, `Ratios of GMT/GMC`) %>% 
    arrange(subgroup, Visit, Arm, `Baseline SARS-CoV-2`, Marker) 

print("Done with table7") 

### Table 8. The differences in the responder rates, 2FRs, 4FRs between groups 
# will be computed along with the two-sided 95% CIs by the Wilson-Score
# method without continuity correction (Newcombe, 1998).
# Output: tab_rrdiff

# tab_rrdiff <- bind_rows(rpcnt %>% 
#                           dplyr::filter(subgroup=="All") %>% 
#                           mutate(Group=Arm, Arm="-"),
#                         rpcnt %>% 
#                           dplyr::filter(subgroup=="All") %>% 
#                           mutate(Group=`Baseline SARS-CoV-2`, `Baseline SARS-CoV-2`="-"),
#                         rpcnt)%>% 
tab_rrdiff <- rpcnt %>%
  dplyr::filter(subgroup %in% c(groups, "All") & grepl("Resp|FR2|FR4",resp_cat)) %>% 
  mutate(groupn = 2-match(Group, c(comp_lev, "Vaccine", "Placebo", "Positive", "Negative"))%%2) %>% 
  pivot_wider(id_cols = c(subgroup, Arm, `Baseline SARS-CoV-2`, Marker, Visit, Ind),
              names_from = groupn, values_from = c(response, ci_l, ci_u), names_sep = "") %>% 
  full_join(distinct(rgmt, subgroup, comp)) %>% 
  mutate(Comparison = case_when(Arm=="-"~"Vaccine vs Placebo",
                          `Baseline SARS-CoV-2`=="-"~"Positive vs Negative",
                          TRUE~comp),
         subgroup = factor(case_when(Arm=="-" ~ "Arm",
                           `Baseline SARS-CoV-2`=="-" ~ "Baseline SARS-CoV-2",
                            TRUE~subgrp[subgroup]), levels=c("Arm", "Baseline SARS-CoV-2", subgrp)),
         Estimate = response1-response2,
         ci_l = Estimate-sqrt((response1-ci_l1)^2+(response2-ci_u2)^2),
         ci_u = Estimate+sqrt((response1-ci_u1)^2+(response2-ci_l2)^2),
         rslt = ifelse(is.na(Estimate), "-", 
                       sprintf("%s\n(%s, %s)", round(Estimate,2), round(ci_l,2), round(ci_u,2)))) %>%
  dplyr::filter(!is.na(Comparison)) %>%
  pivot_wider(id_cols=c(Comparison, subgroup, `Baseline SARS-CoV-2`, Arm, Visit, Marker), 
              names_from = Ind, values_from = rslt) %>%
  select(Comparison, subgroup, `Baseline SARS-CoV-2`, Arm, Visit, Marker, Responder, `% 2-Fold Rise`, `% 4-Fold Rise`) %>% 
  arrange(subgroup, Visit, `Baseline SARS-CoV-2`, Marker, Comparison) 

print("Done with table8") 

# Generate a full table with all the estimates: response rates, GMT, GMTR, etc.
# (Per Peter's email Feb 5, 2021) 
# Vaccine vs Placebo within Baseline status
# groups <- c("Arm")
# comp_i <- c("Vaccine", "Placebo")
# mag_groups <- assays_col
# mag_groups <- mag_groups[sapply(mag_groups, function(x)!all(is.na(ds[,x])))]
# 
# rgmt_Rx <- get_rgmt(ds, mag_groups, groups, comp_lev=comp_i, sub.by="`Baseline SARS-CoV-2`", 
#                     "tps.stratum", "wt.subcohort", "ph2.immuno") %>% 
#   mutate(subgroup=factor(subgrp[subgroup], levels=subgrp))
# 
# rrdiff_Rx <- rpcnt %>% 
#   dplyr::filter(subgroup=="All" & grepl("Resp|FR2|FR4",resp_cat)) %>% 
#   mutate(groupn = 2-match(Arm, comp_i)%%2) %>% 
#   pivot_wider(id_cols = c(subgroup, `Baseline SARS-CoV-2`, Visit, Marker, Ind),
#               names_from = groupn, values_from = c(response, ci_l, ci_u), names_sep = "") %>% 
#   mutate(Estimate = response1-response2,
#          ci_l = Estimate-sqrt((response1-ci_l1)^2+(response2-ci_u2)^2),
#          ci_u = Estimate+sqrt((response1-ci_u1)^2+(response2-ci_l2)^2),
#          rrdiff = sprintf("%s\n(%s, %s)", round(Estimate, 2), 
#                           round(ci_l, 2), round(ci_u, 2)),
#          subgroup=factor(subgrp[subgroup], levels=subgrp)) 
# 
# tab_Rx <- full_join(tab_rr, tab_gm,
#                     c("subgroup", "Arm", "Baseline SARS-CoV-2", 
#                       "Group", "Visit", "N", "Marker")) %>%
#   dplyr::filter(as.character(subgroup) == "All participants") %>%
#   pivot_wider(id_cols = c(subgroup, Group, `Baseline SARS-CoV-2`, Marker, Visit),  
#               names_from = Arm, 
#               values_from = c(N, Responder, `GMT/GMC`)) %>% 
#   inner_join(rrdiff_Rx %>% dplyr::filter(Ind=="Responder"), 
#              by = c("subgroup", "Baseline SARS-CoV-2", "Visit", "Marker")) %>% 
#   inner_join(rgmt_Rx, 
#              by = c("Baseline SARS-CoV-2", "Visit", "Marker")) %>% 
#   select(`Baseline SARS-CoV-2`, Visit, Marker, `N_Vaccine`, `Responder_Vaccine`, 
#          `GMT/GMC_Vaccine`, `N_Placebo`, `Responder_Placebo`, `GMT/GMC_Placebo`, 
#          rrdiff, `Ratios of GMT/GMC`) %>% 
#   arrange(`Baseline SARS-CoV-2`, Visit, Marker)
# 
# tab_neg <- dplyr::filter(tab_Rx, `Baseline SARS-CoV-2` == "Negative") %>% 
#   select(-c(`Baseline SARS-CoV-2`))
# tab_pos <- dplyr::filter(tab_Rx, `Baseline SARS-CoV-2` == "Positive") %>% 
#   select(-c(`Baseline SARS-CoV-2`))
# 
# print("Done with table13-14") 
# 
# ###################################################
# 
# groups <- "`Baseline SARS-CoV-2`"
# comp_i <- c("Positive", "Negative")
# resp_v <- grep("Resp", names(ds), value=T)
# 
# rrdiff_bl <- rpcnt %>% 
#   dplyr::filter(subgroup == "All" & grepl("Resp", resp_cat)) %>% 
#   mutate(groupn = match(as.character(`Baseline SARS-CoV-2`), comp_i), 
#          comp = paste(comp_i, collapse = " vs ")) %>% 
#   pivot_wider(id_cols = c(subgroup, Arm, Visit, Marker),
#               names_from = groupn, 
#               values_from = c(response, ci_l, ci_u), names_sep = "") %>% 
#   mutate(Estimate = response1-response2,
#          ci_l = Estimate-sqrt((response1-ci_l1)^2+(response2-ci_u2)^2),
#          ci_u = Estimate+sqrt((response1-ci_u1)^2+(response2-ci_l2)^2),
#          rrdiff = sprintf("%s\n(%s, %s)", round(Estimate, 2), 
#                           round(ci_l, 2), round(ci_u, 2)),
#          subgroup=factor(subgrp[subgroup], levels=subgrp))  
# 
# rgmt_bl <- get_rgmt(ds, mag_groups, groups, comp_lev = comp_i, sub.by="Arm",
#                     "tps.stratum", "wt.subcohort", "ph2.immuno") %>% 
#   mutate(subgroup=factor(subgrp[subgroup], levels=subgrp))
# 
# 
# tab_bl <- full_join(tab_rr, tab_gm, 
#                     by = c("Arm", "subgroup", "Group",
#                            "Baseline SARS-CoV-2", "Visit", "Marker", "N")) %>%
#   dplyr::filter(as.character(subgroup) == "All participants" & Visit!="Day 1") %>%
#   pivot_wider(id_cols = c(Arm, subgroup, Group, Marker, Visit),
#               names_from=`Baseline SARS-CoV-2`, 
#               values_from=c(N, Responder, `GMT/GMC`)) %>% 
#   inner_join(rrdiff_bl, c("Arm", "subgroup", "Marker", "Visit")) %>% 
#   inner_join(rgmt_bl, c("Arm", "Marker", "Visit")) %>% 
#   select(Arm, Visit, Marker, `N_Positive`, `Responder_Positive`, 
#          `GMT/GMC_Positive`, `N_Negative`, `Responder_Negative`, 
#          `GMT/GMC_Negative`, rrdiff, `Ratios of GMT/GMC`) %>% 
#   arrange(Arm, Visit, Marker)
# 
# tab_vacc <- tab_bl %>% dplyr::filter(Arm == "Vaccine") %>% select(-Arm)
# tab_plcb <- tab_bl %>% dplyr::filter(Arm == "Placebo") %>% select(-Arm)

print("Done with all tables") 

save.results.to <- here::here("output")
if (!dir.exists(save.results.to))  dir.create(save.results.to)
save.results.to <- paste0(here::here("output"), "/", attr(config,"config"))

if (!dir.exists(save.results.to))  dir.create(save.results.to)
print(paste0("save.results.to equals ", save.results.to))


save(list = c("tlf", names(tlf)), file = file.path(save.results.to, sprintf("Tables%s.Rdata", ifelse(exists("COR"), gsub("D15", "D22", COR), ""))))
