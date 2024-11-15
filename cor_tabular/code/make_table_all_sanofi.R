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


tlf <-
  list(
    tab_demo = list(
      table_header = sprintf("Demographic and Clinical Characteristics at Baseline in the Per-Protocol Cohort for D%s %s (Phase 1)", tpeak, substr(COR, nchar(COR)-2, nchar(COR))),
      loop="Stage",
      deselect = c("Stage", "subgroup"),
      pack_row = "subgroup",
      col1="7cm"
    ),
    
    tab_demo_naive = list(
      table_header = sprintf("Demographic and Clinical Characteristics at Baseline in the Naive Per-Protocol Cohort for D%s %s (Phase 1)", tpeak, substr(COR, nchar(COR)-2, nchar(COR))),
      loop="Stage",
      deselect = c("Stage", "subgroup"),
      pack_row = "subgroup",
      col1="7cm"
    ),
    
    tab_demo_nnaive = list(
      table_header = sprintf("Demographic and Clinical Characteristics at Baseline in the Non-Naive Per-Protocol Cohort for D%s %s (Phase 1)", tpeak, substr(COR, nchar(COR)-2, nchar(COR))),
      loop="Stage",
      deselect = c("Stage", "subgroup"),
      pack_row = "subgroup",
      col1="7cm"
    ),
  
    tab_demo_ph2 = list(
      table_header = sprintf("Demographic and Clinical Characteristics at Baseline in the Per-Protocol Cohort for D%s %s (Phase 2)", tpeak, substr(COR, nchar(COR)-2, nchar(COR))),
      # table_footer = randomsubcohort,
      loop="Stage",
      deselect = c("Stage", "subgroup"),
      pack_row = "subgroup",
      col1="7cm"
    ),
    
    tab_demo_naive_ph2 = list(
      table_header = sprintf("Demographic and Clinical Characteristics at Baseline in the Naive Per-Protocol Cohort for D%s %s (Phase 2)", tpeak, substr(COR, nchar(COR)-2, nchar(COR))),
      loop="Stage",
      deselect = c("Stage", "subgroup"),
      pack_row = "subgroup",
      col1="7cm"
    ),
    
    tab_demo_nnaive_ph2 = list(
      table_header = sprintf("Demographic and Clinical Characteristics at Baseline in the Non-Naive Per-Protocol Cohort for D%s %s (Phase 2)", tpeak, substr(COR, nchar(COR)-2, nchar(COR))),
      # table_footer = randomsubcohort,
      loop="Stage",
      deselect = c("Stage", "subgroup"),
      pack_row = "subgroup",
      col1="7cm"
    ),
    
    tab_strtm1 = NULL,
    
    tab_case = list(
      table_header = "Antibody levels in the per-protocol cohort",
      table_footer =c("Stage 1 Day 22: 7-180 days PD2 cases. Stage 2 Day 22: 7-150 days PD2 cases.",
                      "Stage 1 Day 43: 28-180 days PD2 cases. Stage 2 Day 43: 28-150 days PD2 cases."),
      loop="Stage",
      deselect = c("Stage"),
      col_name = c("Arm", "Visit", "Marker", "N", "Resp rate", "GMT/GMC", "N",
                   "Resp rate", "GMT/GMC", "Resp Rate\nDifference", "GMTR/GMCR"),
      header_above1 = c(" "=3, "Cases" = 3, "Non-Cases/Control" = 3,
                        "Comparison" = 2),
      font_size=7.5,
      col1="1cm"),
    
    tab_case_naive = list(
      table_header = "Antibody levels in the naive per-protocol cohort",
      table_footer =c("Stage 1 Day 22: 7-180 days PD2 cases. Stage 2 Day 22: 7-150 days PD2 cases.",
                      "Stage 1 Day 43: 28-180 days PD2 cases. Stage 2 Day 43: 28-150 days PD2 cases."),
      loop="Stage",
      deselect = c("Stage", "Naive"),
      col_name = c("Arm", "Visit", "Marker", "N", "Resp rate", "GMT/GMC", "N",
                   "Resp rate", "GMT/GMC", "Resp Rate\nDifference", "GMTR/GMCR"),
      header_above1 = c(" "=3, "Naive Cases" = 3, "Naive Non-Cases/Control" = 3,
                        "Comparison" = 2),
      font_size=7.5,
      col1="1cm"),

    tab_case_nnaive = list(
      table_header = "Antibody levels in the non-naive per-protocol cohort",
      table_footer =c("Stage 1 Day 22: 7-180 days PD2 cases. Stage 2 Day 22: 7-150 days PD2 cases.",
                      "Stage 1 Day 43: 28-180 days PD2 cases. Stage 2 Day 43: 28-150 days PD2 cases."),
      loop="Stage",
      deselect = c("Stage", "Naive"),
      col_name = c("Arm", "Visit", "Marker", "N", "Resp rate", "GMT/GMC", "N",
                   "Resp rate", "GMT/GMC", "Resp Rate\nDifference", "GMTR/GMCR"),
      header_above1 = c(" "=3, "Non-Naive Cases" = 3, "Non-Naive Non-Cases/Control" = 3,
                        "Comparison" = 2),
      font_size=7.5,
      col1="1cm"),
    
    tab_gmtr = list(
      table_header = sprintf("Geometric mean/concentration titer ratios (GMTRs/GMCRs) for the D%s value compared to the D1 value in the per-protocol cohort", tpeak),
      
      table_footer =c("Stage 1 Day 22: 7-180 days PD2 cases. Stage 2 Day 22: 7-150 days PD2 cases.",
                      "Stage 1 Day 43: 28-180 days PD2 cases. Stage 2 Day 43: 28-150 days PD2 cases."),
      loop="Stage",
      deselect = c("Stage"),
      col_name = c("Arm", "Visit", "Marker", "N", "GMT/GMC", "N", "GMT/GMC", "GMTR/GMCR"),
      header_above1 = c(" "=4, "Cases" = 2, "Non-Cases/Control" = 2),
      col1="5cm"),

    tab_gmtr_naive = list(
      table_header = sprintf("Geometric mean/concentration titer ratios (GMTRs/GMCRs) for D%s value compared to the D1 value in the naive per-protocol cohort", tpeak),
      table_footer =c("Stage 1 Day 22: 7-180 days PD2 cases. Stage 2 Day 22: 7-150 days PD2 cases.",
                      "Stage 1 Day 43: 28-180 days PD2 cases. Stage 2 Day 43: 28-150 days PD2 cases."),
      loop="Stage",
      deselect = c("Stage", "Naive"),
      col_name = c("Arm", "Visit", "Marker", "N", "GMT/GMC", "N", "GMT/GMC", "GMTR/GMCR"),
      header_above1 = c(" "=3, "Naive Cases" = 2, "Naive Non-Cases/Control" = 2),
      col1="5cm"),

    tab_gmtr_nnaive = list(
      table_header = sprintf("Geometric mean/concentration titer ratios (GMTRs/GMCRs) for the D%s value compared to the D1 value in the non-naive per-protocol cohort",tpeak),
      table_footer =c("Stage 1 Day 22: 7-180 days PD2 cases. Stage 2 Day 22: 7-150 days PD2 cases.",
                      "Stage 1 Day 43: 28-180 days PD2 cases. Stage 2 Day 43: 28-150 days PD2 cases."),
      loop="Stage",
      deselect = c("Stage", "Naive"),
      col_name = c("Arm", "Visit", "Marker", "N", "GMT/GMC", "N", "GMT/GMC", "GMTR/GMCR"),
      header_above1 = c(" "=3, "Non-Naive Cases" = 2, "Non-Naive Non-Cases/Control" = 2),
      col1="5cm")

  )


    
timepoints <- config$timepoints

labels.age <-  c("Age $<$ 60", "Age $\\geq$ 60")

labels.minor <- c("Communities of Color", "White Non-Hispanic")

labels.BMI <- c("Underweight BMI < 18.5", "Normal 18.5 $\\leq$ BMI < 25", 
                "Overweight 25 $\\leq$ BMI < 30", "Obese BMI $\\geq$ 30")


labels.assays.short <- labels.assays.short.tabular[assays]

# redefines what is in _common.R to use shorter names
labels.assays.long <- data.frame (purrr::imap_dfc(labels.assays.short, ~ paste0(labels.assays.short[.y], ": ", labels.time)))
rownames(labels.assays.long) <- names(labels.time)

visits <- names(labels.time)#[!grepl("Delta", names(labels.time))]
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
dat <- dat_proc

if(COR %in% c("D22vat08_combined_M6_st1.nAb.batch0and1", "D43vat08_combined_M6_st1.nAb.batch0and1")){
  dat <- dat %>% filter(Trialstage==1)
  Stgn <- 1
  assays_10 <- paste0(assays, "_10")
  rename_col <-  grep(paste0("(", paste(assays_10, collapse="|"), ")"), names(dat), value=T)

  for (icol_10 in rename_col){
    icol <- gsub("_10", "", icol_10)
    if (icol_10 %in% names(dat) & icol %in% names(dat)){
      dat[ ,icol] <- dat[ ,icol_10]
    }
  }

} else if(COR %in% c("D22vat08_combined_M5_nAb", "D43vat08_combined_M5_nAb")) {
  dat <- dat %>% filter(Trialstage==2)
  Stgn <- 2
} else {
  Stgn <- 1:2
}
 
# The stratified random cohort for immunogenicity
if (is.null(dat$EthnicityUnknown)) dat$EthnicityUnknown <- 0
ds_s <- dat %>%
  mutate(
    Trialstage = paste("Stage", Trialstage),
    Arm = factor(ifelse(Trt==1, "Vaccine", "Placebo"), levels=c("Vaccine", "Placebo")),
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
    Naive = ifelse(Bserostatus==0, "Naive", "Non-naive"),
    All = "All participants"
    )

# Step2: Responders
# Post baseline visits
pos.cutoffs <- assay_metadata$pos.cutoff
names(pos.cutoffs) <- assay_metadata$assay
ds <- getResponder(ds_s, times=times[times!="B"], responderFR = 4, #c("B","Day15", "Day29", "Day91"), 
                   assays=assays[!grepl("_mdw", assays)], pos.cutoffs = pos.cutoffs, na.rm = T)



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

###################################################
#             Generating the Tables               #
###################################################
# Setup empty tables 

num_v1 <- c("Age") # Summaries - Mean & Range
num_v2 <- c("BMI") # Summaries - Mean & St.d
cat_v <- c("AgeC", "SexC", "raceC", "ethnicityC", "HighRiskC", "AgeRiskC", "MinorityC")

for (i in names(tlf)){
  assign(i, NULL)
}

tab_dm_ph1 <- lapply(paste("Stage", Stgn), function(x){

  ds_long_ttl_ph1 <- ds %>%
    dplyr::filter(!!as.name(config.cor$ph1)==1) %>%
    dplyr::filter(Trialstage==x) %>%
    bind_rows(mutate(., Arm:="Total")) %>% 
    mutate_all(as.character) %>% 
    pivot_longer(all_of(c(num_v1, num_v2, cat_v)), names_to="subgroup", values_to="subgroup_cat")
  
  # Calculate % for categorical covariates
  dm_cat_ph1 <- inner_join(
    ds_long_ttl_ph1 %>%
      group_by(Trialstage, Arm, subgroup, subgroup_cat) %>%
      summarise(n = n(), .groups = 'drop'),
    ds_long_ttl_ph1 %>%
      group_by(Trialstage, Arm, subgroup) %>%
      summarise(N = n(), .groups = 'drop'),
    by = c("Trialstage", "Arm","subgroup")
  ) %>%
    mutate(pct = n / N,
           rslt1 = sprintf("%s (%.1f%%)", n, n / N * 100), 
           rslt2 = sprintf("%s/%s = %.1f%%", n, N, n / N * 100)) %>% 
    dplyr::filter(subgroup %in% cat_v) 
  
  
  # Calculate mean and range for numeric covariates
  dm_num_ph1 <- ds_long_ttl_ph1 %>%
    dplyr::filter(subgroup %in% c(num_v1, num_v2)) %>% 
    mutate(subgroup_cat=as.numeric(subgroup_cat)) %>%
    group_by(Trialstage, Arm, subgroup) %>%
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
  

  tab_dm_ph1 <- bind_rows(dm_cat_ph1, dm_num_ph1) %>%
    mutate(rslt = case_when(subgroup %in% cat_v ~ rslt1,
                            subgroup %in% num_v1 ~ rslt1,
                            subgroup %in% num_v2 ~ rslt2)) %>%
    mutate(subgroup=ifelse(subgroup %in% c("MinorityC", "raceC"), "RaceEthC", subgroup),
           subgroup=subgrp[subgroup]) %>% 
    dplyr::filter(subgroup_cat %in% char_lev) %>% 
    inner_join(ds_long_ttl_ph1 %>%
                 distinct(Trialstage, Arm, Ptid) %>%
                 group_by(Trialstage, Arm) %>%
                 summarise(tot = n()),
               by = c("Trialstage", "Arm")) %>%
    mutate(Arm = paste0(Arm, "\n(N = ", tot, ")")) %>%
    pivot_wider(c(Trialstage, subgroup, subgroup_cat),
                names_from = Arm,
                names_sort = T,
                values_from = c(rslt)) %>%
    mutate(Characteristics = factor(subgroup_cat, levels=char_lev),
           subgroup=factor(subgroup, levels=subgrp)) %>%
    arrange(Trialstage, subgroup, Characteristics) %>% 
    select(Stage=Trialstage, subgroup, Characteristics, contains("Vaccine"), contains("Placebo"),contains("Total"))
  
  cols <- names(tab_dm_ph1)
  names(tab_dm_ph1) <- str_split(cols, "\n", simplify = T)[,1]
  
  return(list(tab=tab_dm_ph1, col=c(tab_dm_ph1$Stage[1], cols)))
})

tab_demo_col <- tab_dm_ph1 %>% 
  map("col") %>% 
  bind_cols()

tab_demo_col <- structure(tab_demo_col[2:nrow(tab_demo_col), ], .Names=tab_demo_col[1,]) %>% 
  data.frame(check.names = F)

tab_demo <- tab_dm_ph1 %>% 
  map("tab") %>% 
  bind_rows()

tlf$tab_demo$col_name <- tab_demo_col



tab_dm_status <- lapply(paste("Stage", Stgn), function(x){
  
  ds_long_ttl_ph1 <- ds %>%
    dplyr::filter(!!as.name(config.cor$ph1)==1) %>%
    dplyr::filter(Trialstage==x) %>%
    bind_rows(mutate(., Arm:="Total")) %>% 
    mutate_all(as.character) %>% 
    pivot_longer(all_of(c(num_v1, num_v2, cat_v)), names_to="subgroup", values_to="subgroup_cat")
  
  # Calculate % for categorical covariates
  dm_cat_ph1 <- inner_join(
    ds_long_ttl_ph1 %>%
      group_by(Trialstage, Naive, Arm, subgroup, subgroup_cat) %>%
      summarise(n = n(), .groups = 'drop'),
    ds_long_ttl_ph1 %>%
      group_by(Trialstage, Naive, Arm, subgroup) %>%
      summarise(N = n(), .groups = 'drop'),
    by = c("Trialstage", "Naive", "Arm","subgroup")
  ) %>%
    mutate(pct = n / N,
           rslt1 = sprintf("%s (%.1f%%)", n, n / N * 100), 
           rslt2 = sprintf("%s/%s = %.1f%%", n, N, n / N * 100)) %>% 
    dplyr::filter(subgroup %in% cat_v) 
  
  
  # Calculate mean and range for numeric covariates
  dm_num_ph1 <- ds_long_ttl_ph1 %>%
    dplyr::filter(subgroup %in% c(num_v1, num_v2)) %>% 
    mutate(subgroup_cat=as.numeric(subgroup_cat)) %>%
    group_by(Trialstage, Naive, Arm, subgroup) %>%
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
  
 
  tab_dm_ph1 <- bind_rows(dm_cat_ph1, dm_num_ph1) %>%
    mutate(rslt = case_when(subgroup %in% cat_v ~ rslt1,
                            subgroup %in% num_v1 ~ rslt1,
                            subgroup %in% num_v2 ~ rslt2)) %>%
    mutate(subgroup=ifelse(subgroup %in% c("MinorityC", "raceC"), "RaceEthC", subgroup),
           subgroup=subgrp[subgroup]) %>% 
    dplyr::filter(subgroup_cat %in% char_lev) %>% 
    inner_join(ds_long_ttl_ph1 %>%
                 distinct(Trialstage, Naive, Arm, Ptid) %>%
                 group_by(Trialstage, Naive, Arm) %>%
                 summarise(tot = n()),
               by = c("Trialstage", "Naive", "Arm")) %>%
    mutate(Arm = paste0(Arm, "\n(N = ", tot, ")")) %>%
    pivot_wider(c(Trialstage, subgroup, subgroup_cat),
                names_from = c(Naive, Arm),
                names_sort = T,
                values_from = c(rslt)) %>%
    mutate(Characteristics = factor(subgroup_cat, levels=char_lev),
           subgroup=factor(subgroup, levels=subgrp)) %>%
    arrange(Trialstage, subgroup, Characteristics) # %>% 
  
  tab_dm_ph1_naive <- tab_dm_ph1 %>% 
    select(Stage=Trialstage, subgroup, Characteristics, contains("Naive_Vaccine", ignore.case = F), contains("Naive_Placebo", ignore.case = F),contains("Naive_Total", ignore.case = F))
 
  tab_dm_ph1_nnaive <- tab_dm_ph1 %>% 
    select(Stage=Trialstage, subgroup, Characteristics, contains("Non-naive_Vaccine", ignore.case = F), contains("Non-naive_Placebo", ignore.case = F),contains("Non-naive_Total", ignore.case = F))
  
   
  cols_naive <- names(tab_dm_ph1_naive)
  names(tab_dm_ph1_naive) <- str_split(cols_naive, "\n", simplify = T)[,1]
  
  cols_nnaive <- names(tab_dm_ph1_nnaive)
  names(tab_dm_ph1_nnaive) <- str_split(cols_nnaive, "\n", simplify = T)[,1]
  
  return(list(tab_naive=tab_dm_ph1_naive, col_naive=c(tab_dm_ph1_naive$Stage[1], cols_naive),
              tab_nnaive=tab_dm_ph1_nnaive, col_nnaive=c(tab_dm_ph1_nnaive$Stage[1], cols_nnaive)))
})

# Naive Tables
tab_demo_col_naive <- tab_dm_status %>% 
  map("col_naive") %>% 
  bind_cols() 

tab_demo_col_naive <- structure(tab_demo_col_naive[2:nrow(tab_demo_col_naive), ], .Names=tab_demo_col_naive[1,]) %>% 
  data.frame(check.names = F) %>% 
  mutate_all(function(x)gsub("Naive_", "", x))
  

tab_demo_naive <- tab_dm_status %>% 
  map("tab_naive") %>% 
  bind_rows()

tlf$tab_demo_naive$col_name <- tab_demo_col_naive

# Non-Naive Tables
tab_demo_col_nnaive <- tab_dm_status %>% 
  map("col_nnaive") %>% 
  bind_cols()

tab_demo_col_nnaive <- structure(tab_demo_col_nnaive[2:nrow(tab_demo_col_nnaive), ], .Names=tab_demo_col_nnaive[1,]) %>% 
  data.frame(check.names = F) %>% 
  mutate_all(function(x)gsub("Non-naive_", "", x))


tab_demo_nnaive <- tab_dm_status %>% 
  map("tab_nnaive") %>% 
  bind_rows()

tlf$tab_demo_nnaive$col_name <- tab_demo_col_nnaive



tab_dm_ph2 <- lapply(paste("Stage", Stgn), function(x){
  
  ds_long_ttl_ph1 <- ds %>%
    dplyr::filter(!!as.name(config.cor$ph2)==1) %>%
    dplyr::filter(Trialstage==x) %>%
    bind_rows(mutate(., Arm:="Total")) %>% 
    mutate_all(as.character) %>% 
    pivot_longer(all_of(c(num_v1, num_v2, cat_v)), names_to="subgroup", values_to="subgroup_cat")
  
  # Calculate % for categorical covariates
  dm_cat_ph1 <- inner_join(
    ds_long_ttl_ph1 %>%
      group_by(Trialstage, Arm, subgroup, subgroup_cat) %>%
      summarise(n = n(), .groups = 'drop'),
    ds_long_ttl_ph1 %>%
      group_by(Trialstage, Arm, subgroup) %>%
      summarise(N = n(), .groups = 'drop'),
    by = c("Trialstage", "Arm","subgroup")
  ) %>%
    mutate(pct = n / N,
           rslt1 = sprintf("%s (%.1f%%)", n, n / N * 100), 
           rslt2 = sprintf("%s/%s = %.1f%%", n, N, n / N * 100)) %>% 
    dplyr::filter(subgroup %in% cat_v) 
  
  
  # Calculate mean and range for numeric covariates
  dm_num_ph1 <- ds_long_ttl_ph1 %>%
    dplyr::filter(subgroup %in% c(num_v1, num_v2)) %>% 
    mutate(subgroup_cat=as.numeric(subgroup_cat)) %>%
    group_by(Trialstage, Arm, subgroup) %>%
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
  
  
  tab_dm_ph1 <- bind_rows(dm_cat_ph1, dm_num_ph1) %>%
    mutate(rslt = case_when(subgroup %in% cat_v ~ rslt1,
                            subgroup %in% num_v1 ~ rslt1,
                            subgroup %in% num_v2 ~ rslt2)) %>%
    mutate(subgroup=ifelse(subgroup %in% c("MinorityC", "raceC"), "RaceEthC", subgroup),
           subgroup=subgrp[subgroup]) %>% 
    dplyr::filter(subgroup_cat %in% char_lev) %>% 
    inner_join(ds_long_ttl_ph1 %>%
                 distinct(Trialstage, Arm, Ptid) %>%
                 group_by(Trialstage, Arm) %>%
                 summarise(tot = n()),
               by = c("Trialstage", "Arm")) %>%
    mutate(Arm = paste0(Arm, "\n(N = ", tot, ")")) %>%
    pivot_wider(c(Trialstage, subgroup, subgroup_cat),
                names_from = Arm,
                names_sort = T,
                values_from = c(rslt)) %>%
    mutate(Characteristics = factor(subgroup_cat, levels=char_lev),
           subgroup=factor(subgroup, levels=subgrp)) %>%
    arrange(Trialstage, subgroup, Characteristics) %>% 
    select(Stage=Trialstage, subgroup, Characteristics, contains("Vaccine"), contains("Placebo"),contains("Total"))
  
  cols <- names(tab_dm_ph1)
  names(tab_dm_ph1) <- str_split(cols, "\n", simplify = T)[,1]
  
  return(list(tab=tab_dm_ph1, col=c(tab_dm_ph1$Stage[1], cols)))
})

tab_demo_col_ph2 <- tab_dm_ph2 %>% 
  map("col") %>% 
  bind_cols()

tab_demo_col_ph2 <- structure(tab_demo_col_ph2[2:nrow(tab_demo_col_ph2), ], .Names=tab_demo_col_ph2[1,]) %>% 
  data.frame(check.names = F)

tab_demo_ph2 <- tab_dm_ph2 %>% 
  map("tab") %>% 
  bind_rows()

tlf$tab_demo_ph2$col_name <- tab_demo_col_ph2



tab_dm_status_ph2 <- lapply(paste("Stage", Stgn), function(x){
  ds_long_ttl_ph1 <- ds %>%
    dplyr::filter(!!as.name(config.cor$ph2)==1) %>%
    dplyr::filter(Trialstage==x) %>%
    bind_rows(mutate(., Arm:="Total")) %>% 
    mutate_all(as.character) %>% 
    pivot_longer(all_of(c(num_v1, num_v2, cat_v)), names_to="subgroup", values_to="subgroup_cat")
  
  # Calculate % for categorical covariates
  dm_cat_ph1 <- inner_join(
    ds_long_ttl_ph1 %>%
      group_by(Trialstage, Naive, Arm, subgroup, subgroup_cat) %>%
      summarise(n = n(), .groups = 'drop'),
    ds_long_ttl_ph1 %>%
      group_by(Trialstage, Naive, Arm, subgroup) %>%
      summarise(N = n(), .groups = 'drop'),
    by = c("Trialstage", "Naive", "Arm","subgroup")
  ) %>%
    mutate(pct = n / N,
           rslt1 = sprintf("%s (%.1f%%)", n, n / N * 100), 
           rslt2 = sprintf("%s/%s = %.1f%%", n, N, n / N * 100)) %>% 
    dplyr::filter(subgroup %in% cat_v) 
  
  
  # Calculate mean and range for numeric covariates
  dm_num_ph1 <- ds_long_ttl_ph1 %>%
    dplyr::filter(subgroup %in% c(num_v1, num_v2)) %>% 
    mutate(subgroup_cat=as.numeric(subgroup_cat)) %>%
    group_by(Trialstage, Naive, Arm, subgroup) %>%
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
  
  
  tab_dm_ph1 <- bind_rows(dm_cat_ph1, dm_num_ph1) %>%
    mutate(rslt = case_when(subgroup %in% cat_v ~ rslt1,
                            subgroup %in% num_v1 ~ rslt1,
                            subgroup %in% num_v2 ~ rslt2)) %>%
    mutate(subgroup=ifelse(subgroup %in% c("MinorityC", "raceC"), "RaceEthC", subgroup),
           subgroup=subgrp[subgroup]) %>% 
    dplyr::filter(subgroup_cat %in% char_lev) %>% 
    inner_join(ds_long_ttl_ph1 %>%
                 distinct(Trialstage, Naive, Arm, Ptid) %>%
                 group_by(Trialstage, Naive, Arm) %>%
                 summarise(tot = n()),
               by = c("Trialstage", "Naive", "Arm")) %>%
    mutate(Arm = paste0(Arm, "\n(N = ", tot, ")")) %>%
    pivot_wider(c(Trialstage, subgroup, subgroup_cat),
                names_from = c(Naive, Arm),
                names_sort = T,
                values_from = c(rslt)) %>%
    mutate(Characteristics = factor(subgroup_cat, levels=char_lev),
           subgroup=factor(subgroup, levels=subgrp)) %>%
    arrange(Trialstage, subgroup, Characteristics) # %>% 
  
  tab_dm_ph1_naive <- tab_dm_ph1 %>% 
    select(Stage=Trialstage, subgroup, Characteristics, contains("Naive_Vaccine", ignore.case = F), contains("Naive_Placebo", ignore.case = F),contains("Naive_Total", ignore.case = F))
  
  tab_dm_ph1_nnaive <- tab_dm_ph1 %>% 
    select(Stage=Trialstage, subgroup, Characteristics, contains("Non-naive_Vaccine", ignore.case = F), contains("Non-naive_Placebo", ignore.case = F),contains("Non-naive_Total", ignore.case = F))
  
  
  cols_naive <- names(tab_dm_ph1_naive)
  names(tab_dm_ph1_naive) <- str_split(cols_naive, "\n", simplify = T)[,1]
  
  cols_nnaive <- names(tab_dm_ph1_nnaive)
  names(tab_dm_ph1_nnaive) <- str_split(cols_nnaive, "\n", simplify = T)[,1]
  
  return(list(tab_naive=tab_dm_ph1_naive, col_naive=c(tab_dm_ph1_naive$Stage[1], cols_naive),
              tab_nnaive=tab_dm_ph1_nnaive, col_nnaive=c(tab_dm_ph1_nnaive$Stage[1], cols_nnaive)))
})

# Naive
tab_demo_col_naive_ph2 <- tab_dm_status_ph2 %>% 
  map("col_naive") %>% 
  bind_cols()

tab_demo_col_naive_ph2 <- structure(tab_demo_col_naive_ph2[2:nrow(tab_demo_col_naive_ph2), ], .Names=tab_demo_col_naive_ph2[1,]) %>% 
  data.frame(check.names = F) %>% 
  mutate_all(function(x)gsub("Naive_", "", x))

tab_demo_naive_ph2 <- tab_dm_status_ph2 %>% 
  map("tab_naive") %>% 
  bind_rows()

tlf$tab_demo_naive_ph2$col_name <- tab_demo_col_naive_ph2

# Non-naive

tab_demo_col_nnaive_ph2 <- tab_dm_status_ph2 %>% 
  map("col_nnaive") %>% 
  bind_cols()

tab_demo_col_nnaive_ph2 <- structure(tab_demo_col_nnaive_ph2[2:nrow(tab_demo_col_nnaive_ph2), ], .Names=tab_demo_col_nnaive_ph2[1,]) %>% 
  data.frame(check.names = F) %>% 
  mutate_all(function(x)gsub("Non-naive_", "", x))

tab_demo_nnaive_ph2 <- tab_dm_status_ph2 %>% 
  map("tab_nnaive") %>% 
  bind_rows()

tlf$tab_demo_nnaive_ph2$col_name <- tab_demo_col_nnaive_ph2

print("Done with table 1") 



# Cases & Non-cases
# To do these analyses, eligible cases will have EventIndOmicronD43M5hotdeck10==1.  
# (Note the ‘M5’ in the file name).  Would you please re-make your cor_graphical 
# and cor_tabular Stage 2 reports based on this new way to define a case

# From Peter's email on 8/8/2024:
# The Stage 2 trial cor_graphical and cor_tabular reports need to define primary COVID-19 endpoints as 
# Omicron COVID-19 through 150 days post D22, NOT through 180 days post D22.  
# This change is needed because overall vaccine efficacy in non-naives is about ~44% through 5 months and only ~27% through 6 months.
# CoP analyses work better with higher VE.  
# To do these analyses, eligible cases will have EventIndOmicronD43M5hotdeck10==1.  
# (Note the ‘M5’ in the file name).  

# Day 22: Stage 1, 7-180 days PD2 cases; Stage 2,  7-150 days PD2 cases
# Day 43: Stage 1, 28-180 days PD2 cases; Stage 2,  28-150 days PD2 cases


ds <- ds %>% 
  mutate(Case = 
    case_when(Perprotocol ==1 & 
              Trialstage == "Stage 1" &
                EarlyinfectionD43 == 0 & 
                EventIndOmicronD43M6hotdeck10 ==1 &
                EventTimeOmicronD43M6hotdeck10 >= 7 ~  "Cases",
                # EventTimeOmicronD22M6hotdeck10 <= 180 ~  "28-180 days PD2 cases", 
              
              Perprotocol ==1 & 
              Trialstage == "Stage 2" &
                EarlyinfectionD43 == 0 & 
                EventIndOmicronD43M5hotdeck10 ==1 &
                EventTimeOmicronD43M5hotdeck10 >= 7 ~  "Cases",
                # EventTimeOmicronD22M5hotdeck10 <= 150 ~  "28-150 days PD2 cases", 
              
              Perprotocol ==1 & 
              Trialstage == "Stage 1" &
                EarlyinfectionD43 == 0 & 
                EventIndOmicronD22M6hotdeck10==0 ~ "Non-Cases",
              
              Perprotocol ==1 & 
              Trialstage == "Stage 2" &
                EarlyinfectionD43 == 0 & 
                EventIndOmicronD22M5hotdeck10==0 ~ "Non-Cases"))

if (tpeak==22){
  ds <- ds %>% 
    mutate(Case=case_when( 
             # Case == "28-180 days PD2 cases" ~ "7-180 days PD2 cases",
             # Case == "28-150 days PD2 cases" ~ "7-150 days PD2 cases",
             Case == "Cases" ~ "Cases",
             
             Perprotocol ==1 & 
             Trialstage == "Stage 1" &
              EarlyinfectionD22 == 0 & 
              EventIndOmicronD22M6hotdeck10 ==1 &
              EventTimeOmicronD22M6hotdeck10 >= 7 &
              # EventTimeOmicronD43M6hotdeck10 < 7 ~  "7-180 days PD2 cases", 
              EventTimeOmicronD43M6hotdeck10 < 7 ~  "Cases", 
             
             Perprotocol ==1 & 
             Trialstage == "Stage 2" &
               EarlyinfectionD22 == 0 & 
               EventIndOmicronD22M5hotdeck10 ==1 &
               EventTimeOmicronD22M5hotdeck10 >= 7 &
               # EventTimeOmicronD43M5hotdeck10 < 7 ~  "7-150 days PD2 cases", 
               EventTimeOmicronD43M5hotdeck10 < 7 ~  "Cases", 
             
             TRUE ~ Case))
  
}

caseName <- setdiff(unique(ds$Case), c("Non-Cases", NA))

  
  

# Generate a full table with all the estimates: response rates, GMT, GMTR, etc.
# (Per Peter's email Feb 5, 2021)
# Cases vs Non-cases

  
  resp.v <- setdiff(grep("Resp", names(ds), value = T), grep("mdw", names(ds), value = T))
  mag.v <- intersect(assays_col, names(ds))
  gm.v <- intersect(assays_col, grep(config.cor$tpeak, names(ds), value = T))
  gm.v <- gm.v[sapply(gm.v, function(x)any(!is.na(ds[,x])))]
 

  subs=c("Case")
  sub.by=c("Trialstage", "Arm", "Naive")
  comp.i=c("Non-Cases", caseName)
  
  ds <- ds %>% filter(!!as.name(config.cor$ph1))
  rpcnt_case <- get_rr(ds, resp.v, subs, sub.by, strata=config.cor$WtStratum, weights=config.cor$wt, subset=config.cor$ph2) 
  rgm_case <- get_gm(ds, gm.v, subs, sub.by, strata=config.cor$WtStratum, weights=config.cor$wt, subset=config.cor$ph2) 
  rgmt_case <- get_rgmt(ds, gm.v, subs, comp_lev=comp.i, sub.by, strata=config.cor$WtStratum, weights=config.cor$wt, subset=config.cor$ph2) 

  
  rrdiff_case <- rpcnt_case %>%
    mutate(groupn = 2-match(Group, comp.i)%%2) %>%
    pivot_wider(id_cols = c(Trialstage, Naive, Arm, Visit, Marker, Ind),
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
                      by = c("Trialstage", "Naive", "Group","Arm", "N", "Marker", "Visit")) %>%
  filter(N!=0) %>% 
  pivot_wider(id_cols = c(Trialstage, Naive, Arm, Marker, Visit),
              names_from =  Group,
              values_from = c(N, rslt, `GMT/GMC`)) %>%
  full_join(rrdiff_case, by = c("Trialstage", "Naive",  "Arm", "Marker", "Visit")) %>%
  full_join(rgmt_case, by = c("Trialstage", "Naive",  "Arm", "Marker", "Visit"))


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
      mutate(`Ratios of GMT/GMC`=replace_na(`Ratios of GMT/GMC`, "-")) %>% 
      mutate(rrdiff=replace_na(rrdiff, "-")) 
  }


tab_case <- tab_case %>%
  select(Stage=Trialstage, Naive, Arm, Visit, Marker, paste0("N_", caseName), paste0("rslt_", caseName),
         paste0("GMT/GMC_", caseName), `N_Non-Cases`, `rslt_Non-Cases`, `GMT/GMC_Non-Cases`,
         rrdiff, `Ratios of GMT/GMC`) 


tab_gmtr <- tab_case %>% 
  # mutate(Stage=paste("Stage", Stage)) %>% 
  filter(grepl("over D1", Visit), grepl(paste0("D", tpeak), Visit)) %>% 
  select(Stage, Naive, Arm, Visit, Marker, paste0("N_", caseName), paste0("GMT/GMC_", caseName), `N_Non-Cases`, `GMT/GMC_Non-Cases`,`Ratios of GMT/GMC`) %>% 
  arrange(Stage, Naive, Arm, Visit, Marker)


tab_gmtr_naive <- tab_gmtr %>% filter(Naive=="Naive")
tab_gmtr_nnaive <- tab_gmtr %>% filter(Naive=="Non-naive")



tab_case <- tab_case %>% 
  # mutate(Stage=paste("Stage", Stage)) %>% 
  filter(!grepl("over", Visit), Visit==paste("Day", tpeak)) %>% 
  arrange(Stage, Naive, Arm, Visit, Marker)

tab_case_naive <- tab_case %>% filter(Naive=="Naive")
tab_case_nnaive <- tab_case %>% filter(Naive=="Non-naive")



# Total

sub.by=c("Trialstage", "Arm")

rpcnt_case <- get_rr(ds, resp.v, subs, sub.by, strata=config.cor$WtStratum, weights=config.cor$wt, subset=config.cor$ph2)
rgm_case <- get_gm(ds, gm.v, subs, sub.by, strata=config.cor$WtStratum, weights=config.cor$wt, subset=config.cor$ph2)
rgmt_case <- get_rgmt(ds, gm.v, subs, comp_lev=comp.i, sub.by, strata=config.cor$WtStratum, weights=config.cor$wt, subset=config.cor$ph2)


rrdiff_case <- rpcnt_case %>%
  mutate(groupn = 2-match(Group, comp.i)%%2) %>%
  pivot_wider(id_cols = c(Trialstage, Arm, Visit, Marker, Ind),
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
                      by = c("Trialstage", "Group","Arm", "N", "Marker", "Visit")) %>%
  filter(N!=0) %>%
  pivot_wider(id_cols = c(Trialstage, Arm, Marker, Visit),
              names_from =  Group,
              values_from = c(N, rslt, `GMT/GMC`)) %>%
  full_join(rrdiff_case, by = c("Trialstage", "Arm", "Marker", "Visit")) %>%
  full_join(rgmt_case, by = c("Trialstage", "Arm", "Marker", "Visit"))


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
    mutate(`Ratios of GMT/GMC`=replace_na(`Ratios of GMT/GMC`, "-")) %>%
    mutate(rrdiff=replace_na(rrdiff, "-"))
}

tab_case <- tab_case %>%
  select(Stage=Trialstage, Arm, Visit, Marker, paste0("N_", caseName), paste0("rslt_", caseName),
         paste0("GMT/GMC_", caseName), `N_Non-Cases`, `rslt_Non-Cases`, `GMT/GMC_Non-Cases`,
         rrdiff, `Ratios of GMT/GMC`)

tab_gmtr <- tab_case %>%
  # mutate(Stage=paste("Stage", Stage)) %>%
  filter(grepl("over D1", Visit), grepl(paste0("D", tpeak), Visit)) %>%
  select(Stage, Arm, Visit, Marker, paste0("N_", caseName), paste0("GMT/GMC_", caseName), `N_Non-Cases`, `GMT/GMC_Non-Cases`,`Ratios of GMT/GMC`) %>%
  arrange(Stage, Arm, Visit, Marker)

tab_case <- tab_case %>%
  # mutate(Stage=paste("Stage", Stage)) %>%
  filter(!grepl("over", Visit), Visit==paste("Day", tpeak)) %>%
  arrange(Stage, Arm, Visit, Marker)


print("Done with all tables")

# path for tables
save.results.to <- here::here("output")
if (!dir.exists(save.results.to))  dir.create(save.results.to)
save.results.to <- paste0(here::here("output"), "/", attr(config,"config"))

if (!dir.exists(save.results.to))  dir.create(save.results.to)
print(paste0("save.results.to equals ", save.results.to))

save(list = c("tlf", names(tlf)), file = file.path(save.results.to, sprintf("Tables%s.Rdata", ifelse(exists("COR"), gsub("D15", "D22", COR), ""))))
