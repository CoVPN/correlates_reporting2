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
source(here::here("code", "make_functions_covail.R"))


###################################################
#                  Parameters                     #
###################################################
# To select which tables are included in the report.
# Also to modify the headers and footers for each table.


tlf <-
  list(
    tab_demo = list(
      table_header = "Demographic and Clinical Characteristics at Baseline in the Per-Protocol Immunogenicity Cohort",
      # table_footer = randomsubcohort,
      loop="Arms",
      deselect = c("Arms", "subgroup"),
      pack_row = "subgroup",
      col1="7cm"
    ),

    tab_case = list(
      table_header = "Antibody levels in the per-protocol cohort",
      table_footer =c(),
      loop="Arms",
      deselect = c("Arms"),
      col_name = c("Visit", "Marker", "N", "Resp rate", "GMT/GMC", "N",
                   "Resp rate", "GMT/GMC", "Resp Rate\nDifference", "GMTR/GMCR"),
      header_above1 = c(" "=2, "Cases" = 3, "Non-Cases/Control" = 3,
                        "Comparison" = 2),
      # header_above2 = c(" "=2,
      #                   "Baseline SARS-CoV-2 Negative Vaccine Recipients" = 8),
      font_size=9,
      col1="1cm"),
    
    tab_case_naive = list(
      table_header = "Antibody levels in the naive per-protocol cohort",
      table_footer =c(),
      loop="Arms",
      deselect = c("Arms"),
      col_name = c("Visit", "Marker", "N", "Resp rate", "GMT/GMC", "N",
                   "Resp rate", "GMT/GMC", "Resp Rate\nDifference", "GMTR/GMCR"),
      header_above1 = c(" "=2, "Naive Cases" = 3, "Naive Non-Cases/Control" = 3,
                        "Comparison" = 2),
      # header_above2 = c(" "=2,
      #                   "Baseline SARS-CoV-2 Negative Vaccine Recipients" = 8),
      font_size=9,
      col1="1cm"),
    
    tab_case_nnaive = list(
      table_header = "Antibody levels in the non-naive per-protocol cohort",
      table_footer =c(),
      loop="Arms",
      deselect = c("Arms"),
      col_name = c("Visit", "Marker", "N", "Resp rate", "GMT/GMC", "N",
                   "Resp rate", "GMT/GMC", "Resp Rate\nDifference", "GMTR/GMCR"),
      header_above1 = c(" "=2, "Non-Naive Cases" = 3, "Non-Naive Non-Cases/Control" = 3,
                        "Comparison" = 2),
      # header_above2 = c(" "=2,
      #                   "Baseline SARS-CoV-2 Negative Vaccine Recipients" = 8),
      font_size=9,
      col1="1cm"),
    
    tab_gmtr = list(
      table_header = "Geometric mean titer ratios (GMTRs) of ID50 for the D15 or D29 or D91 value compared to the D1 value 
      in the per-protocol cohort",
      
      table_footer =c(),
      loop="Arms",
      deselect = c("Arms"),
      col_name = c("Visit", "Marker", "N", "GMT/GMC", "N", "GMT/GMC", "GMTR/GMCR"),
      header_above1 = c(" "=2, "Cases" = 2, "Non-Cases/Control" = 2),
      col1="5cm"),
    
    
    tab_gmtr_naive = list(
      table_header = "Geometric mean titer ratios (GMTRs) of ID50 for the D15 or 
      D29 or D91 value compared to the D1 value in the naive per-protocol cohort",
      table_footer =c(),
      loop="Arms",
      deselect = c("Arms"),
      col_name = c("Visit", "Marker", "N", "GMT/GMC", "N", "GMT/GMC", "GMTR/GMCR"),
      header_above1 = c(" "=2, "Naive Cases" = 2, "Naive Non-Cases/Control" = 2),
      col1="5cm"),
    
    tab_gmtr_nnaive = list(
      table_header = "Geometric mean titer ratios (GMTRs) of ID50 for the D15 or 
      D29 or D91 value compared to the D1 value in the non-naive per-protocol cohort",
      table_footer =c(),
      loop="Arms",
      deselect = c("Arms"),
      col_name = c("Visit", "Marker", "N", "GMT/GMC", "N", "GMT/GMC", "GMTR/GMCR"),
      header_above1 = c(" "=2, "Non-Naive Cases" = 2, "Non-Naive Non-Cases/Control" = 2),
      col1="5cm")

  )


    
timepoints <- config$timepoints

labels.age <-  c("Age $<$ 65", "Age $\\geq$ 65")

labels.minor <- c("Communities of Color", "White Non-Hispanic")

labels.BMI <- c("Underweight BMI < 18.5", "Normal 18.5 $\\leq$ BMI < 25", 
                "Overweight 25 $\\leq$ BMI < 30", "Obese BMI $\\geq$ 30")

# labels.time <- labels.time[times]

labels.time <- c(B="Day 1", 
                 Day15="Day 15",
                 Delta15overB="D15 fold-rise over D1", 
                 Day29="Day 29", 
                 Delta29overB="D29 fold-rise over D1",
                 Day91="Day 91", 
                 Delta91overB="D91 fold-rise over D1")

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
    # ind == "FR2" ~ "% 2-Fold Rise",
    # ind == "FR4" ~ "% 4-Fold Rise",
    ind == "Resp" ~ "Responder"
    # ind == "2lloq" ~ "% Greater than 2xLLOQ",
    # ind == "4lloq" ~ "% Greater than 4xLLOQ",
    # ind == "2llod" ~ "% Greater than 2xLLOD",
    # ind == "4llod" ~ "% Greater than 4xLLOD"
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
if (is.null(dat$EthnicityUnknown)) dat$EthnicityUnknown <- 0
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
    AgeC = ifelse(Senior == 1, labels.age[2], labels.age[1]),
    SexC = ifelse(Sex == 1, "Female", "Male"),
    AgeSexC = paste(AgeC, SexC),
    AgeMinorC = ifelse(is.na(MinorityC), NA, paste(AgeC, MinorityC)),
    Naive = ifelse(naive==1, "Naive", "Non-naive"),
    All = "All participants"
    )

  ncol1 <- ncol(ds_s)

  ds_s <- ds_s %>% 
    mutate(Trt1=ifelse(TrtonedosemRNA==1, "One-dose mRNA booster", ""),
           Trt2=ifelse(arm %in% c(13, 14, 15), "pooled Sanofi", ""),
           Trt3=ifelse(TrtA==1, "mRNA Moderna", ""),
           Trt4=ifelse(TrtA==0, "mRNA Pfizer", ""),
           Trt5=ifelse(TrtB==1, "mRNA Prototype", ""),
           Trt6=ifelse(TrtB==1 & grepl("moderna", tolower(treatment_actual)), "mRNA Prototype Moderna", ""),
           Trt7=ifelse(TrtB==1 & grepl("pfizer", tolower(treatment_actual)), "mRNA Prototype Pfizer", ""),
           Trt8=ifelse(TrtB==0, "mRNA Omicron-Containing", ""),
           Trt9=ifelse(TrtB==0 & grepl("moderna", tolower(treatment_actual)), "mRNA Omicron-Containing Moderna", ""),
           Trt10=ifelse(TrtB==0 & grepl("pfizer", tolower(treatment_actual)), "mRNA Omicron-Containing Pfizer", ""),
           Trt11=ifelse(TrtC==1, "mRNA Bivalent", ""),
           Trt12=ifelse(TrtC==0, "mRNA Monovalent", "")
           )
  ncol2 <- ncol(ds_s)
  Trtn <- ncol2 - ncol1
# Step2: Responders
# Post baseline visits
pos.cutoffs <- assay_metadata$pos.cutoff
names(pos.cutoffs) <- assay_metadata$assay
ds <- getResponder(ds_s, times=c("Day15", "Day29", "Day91"), 
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
num_v2 <- NULL # Summaries - Mean & St.d
cat_v <- c("AgeC", "SexC", "raceC", "ethnicityC")

for (i in names(tlf)){
  assign(i, NULL)
}

tab_dm_ph1 <- lapply(1:Trtn, function(x){
  Trti <- paste0("Trt",x)
  
  ds_long_ttl_ph1 <- ds %>%
    dplyr::filter(Perprotocol==1, !!as.name(config.cor$ph1)==1) %>%
    dplyr::filter(!!as.name(Trti)!="") %>%
    bind_rows(mutate(., Naive:="Total")) %>% 
    mutate_all(as.character) %>% 
    pivot_longer(all_of(c(num_v1, num_v2, cat_v)), names_to="subgroup", values_to="subgroup_cat")
  
  # Calculate % for categorical covariates
  dm_cat_ph1 <- inner_join(
    ds_long_ttl_ph1 %>%
      group_by(!!as.name(Trti), Naive, subgroup, subgroup_cat) %>%
      summarise(n = n(), .groups = 'drop'),
    ds_long_ttl_ph1 %>%
      group_by(!!as.name(Trti), Naive, subgroup) %>%
      summarise(N = n(), .groups = 'drop'),
    by = c(Trti, "Naive","subgroup")
  ) %>%
    mutate(pct = n / N,
           rslt1 = sprintf("%s (%.1f%%)", n, n / N * 100), 
           rslt2 = sprintf("%s/%s = %.1f%%", n, N, n / N * 100)) %>% 
    dplyr::filter(subgroup %in% cat_v) 
  
  
  # Calculate mean and range for numeric covariates
  dm_num_ph1 <- ds_long_ttl_ph1 %>%
    dplyr::filter(subgroup %in% c(num_v1, num_v2)) %>% 
    mutate(subgroup_cat=as.numeric(subgroup_cat)) %>%
    group_by(!!as.name(Trti), Naive, subgroup) %>%
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
    mutate(subgroup=ifelse(subgroup %in% c("MinorityC", "raceC"), "RaceEthC", subgroup)) %>% 
    dplyr::filter(subgroup_cat %in% char_lev) %>% 
    inner_join(ds_long_ttl_ph1 %>% 
                 distinct(!!as.name(Trti), Naive, Ptid) %>% 
                 group_by(!!as.name(Trti), Naive) %>%
                 summarise(tot = n()),
               by = c(Trti, "Naive")) %>% 
    mutate(Naive = paste0(Naive, "\n(N = ", tot, ")"), subgroup=subgrp[subgroup]) %>%
    pivot_wider(c(!!as.name(Trti), subgroup, subgroup_cat, rslt),
                names_from = Naive,
                names_sort = T,
                values_from = c(rslt)) %>%
    mutate(Characteristics = factor(subgroup_cat, levels=char_lev),
           subgroup=factor(subgroup, levels=subgrp)) %>%
    arrange(!!as.name(Trti), subgroup, Characteristics) %>% 
    select(Arms=Trti, subgroup, Characteristics, contains("Naive"), contains("Non-naive"),contains("Total"))
  
  cols <- names(tab_dm_ph1)
  names(tab_dm_ph1) <- str_split(cols, "\n", simplify = T)[,1]
  
  return(list(tab=tab_dm_ph1, col=cols))
})

tab_demo_col <- tab_dm_ph1 %>% 
  map("col")

tab_demo <- tab_dm_ph1 %>% 
  map("tab") %>% 
  bind_rows()

tlf$tab_demo$col_name <- tab_demo_col

print("Done with table 1") 

# Cases & Non-cases



ds <- ds %>%
  mutate(
    # EventIndPrimaryD1 = ifelse(study_name=="VAT08m" & grepl("omi", COR), EventIndOmicronD1, EventIndPrimaryD1),
    Case = case_when(Perprotocol==1 &
                            # Immunemarkerset==1 &
                            # !!as.name(config.cor$Earlyendpoint)==0 &
                            !!as.name(config.cor$EventIndPrimary)==1 ~ "Cases",
                    Perprotocol==1 &
                            # NAntibody=="NEGATIVE" &
                            !!as.name(config.cor$EventIndPrimary)==0 ~ "Non-Cases")
    )

if (any(grepl("SevereEventIndPrimary", names(ds)))) {
  ds <- ds %>%
    mutate(Case = case_when(Case=="Cases" & !!as.name(paste0("SevereEventIndPrimaryIncludeNotMolecConfirmedD", config.cor$tpeak))==1 ~ "Severe Cases",
                            Case=="Cases" & !!as.name(paste0("SevereEventIndPrimaryIncludeNotMolecConfirmedD", config.cor$tpeak))==0 ~ "Non-Cases",
                            Case=="Non-Cases" ~ "Non-Cases"))
}


# Generate a full table with all the estimates: response rates, GMT, GMTR, etc.
# (Per Peter's email Feb 5, 2021)
# Cases vs Non-cases

  sub.by <- c("All")
  comp.i <- c("Non-Cases", "Cases")
  # resp.v <- intersect(grep("Resp", names(ds), value = T),
  #                     grep(config.cor$tpeak, names(ds), value = T))
  # mag.v <- intersect(assays_col, grep(config.cor$tpeak, names(ds), value = T))
  
  resp.v <- grep("Resp", names(ds), value = T)
  mag.v <- intersect(assays_col, names(ds))
  
  ds.l.resp <- filter(ds, !!as.name(config.cor$ph1)==1) %>% 
    pivot_longer(cols=all_of(resp.v), names_to = "resp_cat", values_to = "resp_value") %>% 
    mutate(assay=gsub("Resp", "", resp_cat))
  
  ds.l.mag <- filter(ds, !!as.name(config.cor$ph1)==1) %>% 
    pivot_longer(cols=all_of(mag.v), names_to = "mag_cat", values_to = "mag_value") %>% 
    mutate(assay=mag_cat)
  
  ds.l <- full_join(ds.l.mag, ds.l.resp)
    

tab_assay <- lapply(1:Trtn, function(x, dat=ds.l){
  Trti <- paste0("Trt",x)
  ds.i <- dat %>% dplyr::filter(!!as.name(Trti)!="")
  rpcnt_case <- ds.i %>% 
    group_by(!!as.name(Trti), resp_cat, Case) %>% 
    summarise(N=sum(!is.na(resp_value)), Pos=sum(resp_value, na.rm = T), .groups="drop") %>% 
    rowwise() %>% 
    mutate(response=Pos/N, ci_l=exactci(Pos, N, .95)$conf.int[1], ci_u=exactci(Pos, N, .95)$conf.int[2], 
           rslt=sprintf("%s/%s = %.1f%%\n(%.1f%%, %.1f%%)", Pos, N, response*100, ci_l*100, ci_u*100)) %>% 
    inner_join(distinct(labels_all, resp_cat, Visit, Marker, Ind), by = "resp_cat") %>% 
    filter(N!=0)
  
  rgm_case <- ds.i %>% 
    group_by(!!as.name(Trti), mag_cat, Case) %>% 
    summarise(N=sum(!is.na(mag_value)), lgmt=mean(mag_value, na.rm=T), lsd=sd(mag_value, na.rm=T), .groups="drop") %>% 
    rowwise() %>% 
    mutate(gmt=10^lgmt, ci_l=10^(lgmt+qt(.025, N-1)*lsd/sqrt(N)), ci_u=10^(lgmt+qt(.975, N-1)*lsd/sqrt(N)), 
           `GMT/GMC`=sprintf("%.1f\n(%.1f, %.1f)", gmt, ci_l, ci_u)) %>% 
    inner_join(distinct(labels_all, mag_cat, Visit, Marker), by = "mag_cat") 
    # filter(N!=0)
  
  rgmt_excl <- rgm_case %>% 
    filter(N==0) %>% 
    select(!!as.name(Trti), mag_cat)
  
  rgmt_case <- ds.i %>% 
    anti_join(rgmt_excl) %>% 
    mutate(Case=relevel(factor(Case), ref="Cases")) %>%
    filter(!is.na(mag_value)) %>% 
    group_by(!!as.name(Trti), mag_cat) %>% 
    summarise(N=sum(!is.na(mag_value)), 
              est=glm(mag_value~Case)$coefficients[2], 
              ci_l=confint(glm(mag_value~Case))[2,1], 
              ci_u=confint(glm(mag_value~Case))[2,2], .groups="drop") %>% 
    mutate(`Ratios of GMT/GMC`=sprintf("%.2f\n(%.2f, %.2f)", 10^est, 10^ci_l, 10^ci_u), comp="Non-Cases vs Cases") %>% 
    inner_join(distinct(labels_all, mag_cat, Visit, Marker))

  
rrdiff_case <- rpcnt_case %>%
  # dplyr::filter(subgroup %in% subs & grepl("Resp",resp_cat)) %>%
  mutate(groupn = 2-match(Case, comp.i)%%2) %>%
  pivot_wider(id_cols = c(Trti, Visit, Marker, Ind),
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
                      by = c(Trti, "Case", "N", "Marker", "Visit")) %>%
  filter(N!=0) %>% 
  pivot_wider(id_cols = c(!!as.name(Trti), Marker, Visit),
              names_from = Case,
              values_from = c(N, rslt, `GMT/GMC`)) %>%
  full_join(rrdiff_case, by = c(Trti, "Marker", "Visit")) %>%
  full_join(rgmt_case, by = c(Trti, "Marker", "Visit"))


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
  select(Arms=Trti, Visit, Marker, `N_Cases`, `rslt_Cases`,
         `GMT/GMC_Cases`, `N_Non-Cases`, `rslt_Non-Cases`, `GMT/GMC_Non-Cases`,
         rrdiff, `Ratios of GMT/GMC`) 

tab_gmtr <- tab_case %>% 
  filter(grepl("over", Visit)) %>% 
  select(Arms, Visit, Marker, `N_Cases`, `GMT/GMC_Cases`, `N_Non-Cases`, `GMT/GMC_Non-Cases`,`Ratios of GMT/GMC`) %>% 
  arrange(Arms, Visit, Marker)

tab_case <- tab_case %>% 
  filter(!grepl("over", Visit)) %>% 
  arrange(Arms, Visit, Marker)


return(list(tab_case=tab_case, tab_gmtr=tab_gmtr))
})

tab_case <- tab_assay %>% 
  map("tab_case") %>% 
  bind_rows()

tab_gmtr <- tab_assay %>% 
  map("tab_gmtr") %>% 
  bind_rows()



tab_assay_status <- lapply(1:Trtn, function(x, dat=ds.l){
  Trti <- paste0("Trt",x)
  ds.i <- dat %>% dplyr::filter(!!as.name(Trti)!="")
  rpcnt_case <- ds.i %>% 
    group_by(!!as.name(Trti), resp_cat, Naive, Case) %>% 
    summarise(N=sum(!is.na(resp_value)), Pos=sum(resp_value, na.rm = T), .groups="drop") %>% 
    rowwise() %>% 
    mutate(response=Pos/N, ci_l=exactci(Pos, N, .95)$conf.int[1], ci_u=exactci(Pos, N, .95)$conf.int[2], 
           rslt=sprintf("%s/%s = %.1f%%\n(%.1f%%, %.1f%%)", Pos, N, response*100, ci_l*100, ci_u*100)) %>% 
    inner_join(distinct(labels_all, resp_cat, Visit, Marker, Ind), by = "resp_cat") %>% 
    filter(N!=0)
  
  rgm_case <- ds.i %>% 
    mutate(Case=factor(Case), Naive=factor(Naive)) %>% 
    group_by(!!as.name(Trti), mag_cat, Naive, Case, .drop = F) %>% 
    summarise(N=sum(!is.na(mag_value)), lgmt=mean(mag_value, na.rm=T), lsd=sd(mag_value, na.rm=T), .groups="drop") %>% 
    rowwise() %>% 
    mutate(gmt=10^lgmt, ci_l=10^(lgmt+qt(.025, N-1)*lsd/sqrt(N)), ci_u=10^(lgmt+qt(.975, N-1)*lsd/sqrt(N)), 
           `GMT/GMC`=sprintf("%.1f\n(%.1f, %.1f)", gmt, ci_l, ci_u)) %>% 
    inner_join(distinct(labels_all, mag_cat, Visit, Marker), by = "mag_cat") #%>% 
    # filter(N!=0)
  
  rgmt_excl <- rgm_case %>% 
    filter(N==0) %>% 
    select(!!as.name(Trti), mag_cat, Naive)
  
  rgmt_case <- ds.i %>% 
    anti_join(rgmt_excl) %>% 
    mutate(Case=relevel(factor(Case), ref="Cases")) %>%
    filter(!is.na(mag_value)) %>% 
    # filter(mag_cat%in%c('Bpseudoneutid50_BA.1')) %>% 
    group_by(!!as.name(Trti), Naive, mag_cat) %>% 
    summarise(N=sum(!is.na(mag_value)), 
              est=glm(mag_value~Case)$coefficients[2], 
              ci_l=confint(glm(mag_value~Case))[2,1], 
              ci_u=confint(glm(mag_value~Case))[2,2], .groups="drop") %>% 
    mutate(`Ratios of GMT/GMC`=sprintf("%.2f\n(%.2f, %.2f)", 10^est, 10^ci_l, 10^ci_u), comp="Non-Cases vs Cases") %>% 
    inner_join(distinct(labels_all, mag_cat, Visit, Marker))
  
  
  rrdiff_case <- rpcnt_case %>%
    # dplyr::filter(subgroup %in% subs & grepl("Resp",resp_cat)) %>%
    mutate(groupn = 2-match(Case, comp.i)%%2) %>%
    pivot_wider(id_cols = c(Trti, Visit, Marker, Ind, Naive),
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
                        by = c(Trti, "Case", "N", "Marker", "Visit", "Naive")) %>%
    filter(N!=0) %>% 
    pivot_wider(id_cols = c(!!as.name(Trti), Marker, Visit, Naive),
                names_from = Case,
                values_from = c(N, rslt, `GMT/GMC`)) %>%
    full_join(rrdiff_case, by = c(Trti, "Marker", "Visit", "Naive")) %>%
    full_join(rgmt_case, by = c(Trti, "Marker", "Visit", "Naive"))
  
  
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

  tab_case_naive <- tab_case %>%
    filter(Naive=="Naive", !grepl("over", Visit)) %>% 
    select(Arms=Trti, Visit, Marker, `N_Cases`, `rslt_Cases`,
           `GMT/GMC_Cases`, `N_Non-Cases`, `rslt_Non-Cases`, `GMT/GMC_Non-Cases`,
           rrdiff, `Ratios of GMT/GMC`) %>% 
    arrange(Arms, Visit, Marker)
  
  tab_case_nnaive <- tab_case %>%
    filter(Naive=="Non-naive", !grepl("over", Visit)) %>% 
    select(Arms=Trti, Visit, Marker, `N_Cases`, `rslt_Cases`,
           `GMT/GMC_Cases`, `N_Non-Cases`, `rslt_Non-Cases`, `GMT/GMC_Non-Cases`,
           rrdiff, `Ratios of GMT/GMC`)  %>% 
    arrange(Arms, Visit, Marker)
  
  tab_gmtr_naive <- tab_case %>% 
    filter(Naive=="Naive", grepl("over", Visit)) %>% 
    select(Arms=Trti, Visit, Marker, `N_Cases`, `GMT/GMC_Cases`, `N_Non-Cases`, `GMT/GMC_Non-Cases`,`Ratios of GMT/GMC`) %>% 
    arrange(Arms, Visit, Marker)
  

  tab_gmtr_nnaive <- tab_case %>% 
    filter(Naive=="Non-naive", grepl("over", Visit)) %>% 
    select(Arms=Trti, Visit, Marker, `N_Cases`, `GMT/GMC_Cases`, `N_Non-Cases`, `GMT/GMC_Non-Cases`,`Ratios of GMT/GMC`) %>% 
    arrange(Arms, Visit, Marker)
  
  return(list(tab_case_naive=tab_case_naive, tab_case_nnaive=tab_case_nnaive, 
              tab_gmtr_naive=tab_gmtr_naive, tab_gmtr_nnaive=tab_gmtr_nnaive))
})

tab_case_naive <- tab_assay_status %>% 
  map("tab_case_naive") %>% 
  bind_rows()

tab_case_nnaive <- tab_assay_status %>% 
  map("tab_case_nnaive") %>% 
  bind_rows()

tab_gmtr_naive <- tab_assay_status %>% 
  map("tab_gmtr_naive") %>% 
  bind_rows()

tab_gmtr_nnaive <- tab_assay_status %>% 
  map("tab_gmtr_nnaive") %>% 
  bind_rows()

print("Done with all tables")

# path for tables
save.results.to <- here::here("output")
if (!dir.exists(save.results.to))  dir.create(save.results.to)
save.results.to <- paste0(here::here("output"), "/", attr(config,"config"))

if (!dir.exists(save.results.to))  dir.create(save.results.to)
print(paste0("save.results.to equals ", save.results.to))

save(list = c("tlf", names(tlf)), file = file.path(save.results.to, sprintf("Tables%s.Rdata", ifelse(exists("COR"), gsub("D15", "D22", COR), ""))))
