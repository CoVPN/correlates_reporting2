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
      table_header = "Demographic and Clinical Characteristics at Baseline in the Per-Protocol Immunogenicity Cohort",
      # table_footer = randomsubcohort,
      loop="Arms",
      deselect = c("Arms", "subgroup"),
      pack_row = "subgroup",
      col1="7cm"
    ),
  
    tab_strtm1 = NULL,
	
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
      font_size=7,
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
      font_size=7,
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
      font_size=7,
      col1="1cm"),
    
    tab_case_nevercase = list(
      table_header = "Antibody levels by Case and Never-Case in the per-protocol cohort",
      table_footer =c("Never-Case are participants who were never identified as cases at any timepoints."),
      loop="Arms",
      deselect = c("Arms"),
      col_name = c("Visit", "Marker", "N", "Resp rate", "GMT/GMC", "N",
                   "Resp rate", "GMT/GMC", "Resp Rate\nDifference", "GMTR/GMCR"),
      header_above1 = c(" "=2, "Cases" = 3, "Never-Cases/Control" = 3,
                        "Comparison" = 2),
      # header_above2 = c(" "=2,
      #                   "Baseline SARS-CoV-2 Negative Vaccine Recipients" = 8),
      font_size=7,
      col1="1cm"),
    
    tab_case_naive_nevercase = list(
      table_header = "Antibody levels by Case and Never-Case in the naive per-protocol cohort",
      table_footer =c("Never-Case are participants who were never identified as cases at any timepoints."),
      loop="Arms",
      deselect = c("Arms"),
      col_name = c("Visit", "Marker", "N", "Resp rate", "GMT/GMC", "N",
                   "Resp rate", "GMT/GMC", "Resp Rate\nDifference", "GMTR/GMCR"),
      header_above1 = c(" "=2, "Naive Cases" = 3, "Naive Never-Cases/Control" = 3,
                        "Comparison" = 2),
      # header_above2 = c(" "=2,
      #                   "Baseline SARS-CoV-2 Negative Vaccine Recipients" = 8),
      font_size=7,
      col1="1cm"),
    
    tab_case_nnaive_nevercase = list(
      table_header = "Antibody levels by Case and Never-Case in the non-naive per-protocol cohort",
      table_footer =c("Never-Case are participants who were never identified as cases at any timepoints."),
      loop="Arms",
      deselect = c("Arms"),
      col_name = c("Visit", "Marker", "N", "Resp rate", "GMT/GMC", "N",
                   "Resp rate", "GMT/GMC", "Resp Rate\nDifference", "GMTR/GMCR"),
      header_above1 = c(" "=2, "Non-Naive Cases" = 3, "Never-Naive Non-Cases/Control" = 3,
                        "Comparison" = 2),
      # header_above2 = c(" "=2,
      #                   "Baseline SARS-CoV-2 Negative Vaccine Recipients" = 8),
      font_size=7,
      col1="1cm"),
    
    tab_gmtr = list(
      table_header = "Geometric mean titer ratios (GMTRs) of Responses for the D15 or D91 or D181 value compared to the D1 value 
      in the per-protocol cohort",
      
      table_footer =c(),
      loop="Arms",
      deselect = c("Arms"),
      col_name = c("Visit", "Marker", "N", "GMT/GMC", "N", "GMT/GMC", "GMTR/GMCR"),
      header_above1 = c(" "=2, "Cases" = 2, "Non-Cases/Control" = 2),
	  font_size=7,
      col1="5cm"),
    
    
    tab_gmtr_naive = list(
      table_header = "Geometric mean titer ratios (GMTRs) of Responses for the D15 or 
      D91 or D181 value compared to the D1 value in the naive per-protocol cohort",
      table_footer =c(),
      loop="Arms",
      deselect = c("Arms"),
      col_name = c("Visit", "Marker", "N", "GMT/GMC", "N", "GMT/GMC", "GMTR/GMCR"),
      header_above1 = c(" "=2, "Naive Cases" = 2, "Naive Non-Cases/Control" = 2),
	  font_size=7,
      col1="5cm"),
    
    tab_gmtr_nnaive = list(
      table_header = "Geometric mean titer ratios (GMTRs) of Responses for the D15 or 
      D91 or D181 value compared to the D1 value in the non-naive per-protocol cohort",
      table_footer =c(),
      loop="Arms",
      deselect = c("Arms"),
      col_name = c("Visit", "Marker", "N", "GMT/GMC", "N", "GMT/GMC", "GMTR/GMCR"),
      header_above1 = c(" "=2, "Non-Naive Cases" = 2, "Non-Naive Non-Cases/Control" = 2),
	  font_size=7,
      col1="5cm")

  )

tlf <- lapply(tlf, function(x){
  x$table_header <- paste0(COR,": ", x$table_header)
  x
})


    
# timepoints <- config$timepoints

labels.age <-  c("Age $<$ 65", "Age $\\geq$ 65")

labels.minor <- c("Communities of Color", "White Non-Hispanic")

labels.BMI <- c("Underweight BMI < 18.5", "Normal 18.5 $\\leq$ BMI < 25", 
                "Overweight 25 $\\leq$ BMI < 30", "Obese BMI $\\geq$ 30")

# labels.time <- labels.time[times]

labels.time <- c(B="Day 1",
                 Day15="Day 15",
                 Delta15overB="D15 fold-rise over D1",
                 Day91="Day 91",
                 Delta91overB="D91 fold-rise over D1",
                 Day181="Day 181",
                 Delta181overB="D181 fold-rise over D1")


if(grepl("frnt", COR)) {
  assays <- assays[substr(assays, 1, 4)=="frnt"]
  } else {
  assays <- assays[substr(assays, 1, 2)!="cd"]
  }

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
  # ind = c("Resp", "FR2", "FR4", "2lloq", "4lloq", "2llod", "4llod"), 
  ind = "Resp", 
  stringsAsFactors = F
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

# labels_all <- labels_all %>% filter(!substr(marker, nchar(marker)-1, nchar(marker)) %in% c("S1", "S2"))

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
    mutate(Trt1=ifelse(TrtonedosemRNA==1 & stage %in% c(1:2), "One-dose mRNA arms in Stages 1 and 2", ""),
           Trt2=ifelse(TrtA==1, "One-dose mRNA arms in Stage 1", ""),
           Trt3=ifelse(TrtA==0, "One-dose mRNA arms in Stage 2", ""),
           Trt4=ifelse(arm %in% c(13, 14, 15), "One-dose recombinant protein arms in Stage 3", ""),
           Trt5=ifelse(arm %in% c(1:2, 5:9, 12:15), "One-dose mRNA and recombinant protein arms", "")
           )
  
  ncol2 <- ncol(ds_s)
  Trtn <- ncol2 - ncol1

# Step2: Responders
# Post baseline visits

# assay_metadata <- assay_metadata %>% filter(panel=="tcell")
  
pos.cutoffs <- assay_metadata$pos.cutoff
names(pos.cutoffs) <- assay_metadata$assay

ds <- getResponder(ds_s, times=c("B", paste0("Day", c(15, 91, 181))),
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
    dplyr::filter(!!as.name(paste0("ph1.D15", ".tcell"[grepl("tcell", COR)]))==1) %>%
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
    mutate(subgroup=ifelse(subgroup %in% c("MinorityC", "raceC"), "RaceEthC", subgroup),
           subgroup=subgrp[subgroup]) %>% 
    dplyr::filter(subgroup_cat %in% char_lev) %>% 
    inner_join(ds_long_ttl_ph1 %>%
                 distinct(!!as.name(Trti), Naive, Ptid) %>%
                 group_by(!!as.name(Trti), Naive) %>%
                 summarise(tot = n()),
               by = c(Trti, "Naive")) %>%
    mutate(Naive = paste0(Naive, "\n(N = ", tot, ")")) %>%
    pivot_wider(c(!!as.name(Trti), subgroup, subgroup_cat),
                names_from = Naive,
                names_sort = T,
                values_from = c(rslt)) %>%
    mutate(Characteristics = factor(subgroup_cat, levels=char_lev),
           subgroup=factor(subgroup, levels=subgrp)) %>%
    arrange(!!as.name(Trti), subgroup, Characteristics) %>% 
    select(Arms=Trti, subgroup, Characteristics, contains("Naive"), contains("Non-naive"),contains("Total"))
  
  cols <- names(tab_dm_ph1)
  names(tab_dm_ph1) <- str_split(cols, "\n", simplify = T)[,1]
  
  return(list(tab=tab_dm_ph1, col=c(tab_dm_ph1$Arms[1], cols)))
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

print("Done with table 1") 

# Cases & Non-cases


# Per-protocol is defined by ph1.D15
ds <- ds %>%
  mutate(
    Case = case_when(ph1.D15.frnt==1 &
                            !!as.name(config.cor$EventIndPrimary)==1 ~ "Cases",
                     ph1.D15.frnt==1 &
                            !!as.name(config.cor$EventIndPrimary)==0 ~ "Non-Cases"), 
    CaseNevercase = case_when(ph1.D15.frnt==1 &
                                !!as.name(config.cor$EventIndPrimary)==1 ~ "Cases",
                              ph1.D15.frnt==1 & COVIDIndD22toD181==0 ~ "Never-Cases")
    )

# FS_col <- intersect(grep("FS", names(ds), value=T), assays_col)

# ds_lgFS <- ds %>% mutate_at(FS_col, log10)
# Generate a full table with all the estimates: response rates, GMT, GMTR, etc.
# (Per Peter's email Feb 5, 2021)
# Cases vs Non-cases


  # subs <- c("Case")
  resp.v <- intersect(paste0(assays_col, "Resp"), grep("Resp", names(ds), value = T))
  mag.v <- intersect(assays_col, names(ds))
  
tab_assay_list <- function(x, subs, sub.by.add=NULL, comp.i=c("Non-Cases", "Cases")){
    print(x)
    Trti <- paste0("Trt",x)
    sub.by <- c(Trti, sub.by.add)
    n_GM <- ifelse(grepl("tcell", COR), 3, 1) 
    
    ds.i <- ds %>% dplyr::filter(ph1.D15.frnt==1, !!as.name(Trti)!="")
    
    resp.v <- resp.v[sapply(resp.v, function(x)any(!is.na(ds.i[,x])))]
    mag.v <- mag.v[sapply(mag.v, function(x)any(!is.na(ds.i[,x])))]

    rpcnt_case <- get_rr(ds.i, resp.v, subs, sub.by, strata=config.cor$WtStratum, weights=config.cor$wt, subset="ph2.D15.frnt") %>%
      rename(!!subs:=Group)

    
    rgm_case <- get_gm(ds.i, mag.v, subs, sub.by, strata=config.cor$WtStratum, weights=config.cor$wt, subset="ph2.D15.frnt", digits=n_GM) %>%
      rename(!!subs:=Group) %>% 
      mutate(`GMT/GMC` =ifelse(grepl("FS", mag_cat), sprintf("%.3f\n(%.3f, %.3f)", mag, ci_l, ci_u), `GMT/GMC`))
    
    
    rgmt_excl <- rgm_case %>% 
      filter(N==0) %>% 
      select(all_of(sub.by), mag_cat)

    # gmt.v <- setdiff(mag.v, grep("Delta", mag.v, value=T))
    gmt.v <- mag.v[sapply(mag.v, function(x)!any(is.infinite(ds.i[, x]), na.rm = T))]
    rgmt_case <- get_rgmt(ds.i, gmt.v, subs, comp_lev=comp.i, sub.by, strata=config.cor$WtStratum, weights=config.cor$wt, subset="ph2.D15.frnt") %>% 
    mutate(`Ratios of GMT/GMC`=ifelse(grepl("FS", mag_cat), sprintf("%.3f\n(%.3f, %.3f)", Estimate, ci_l.Estimate, ci_u.Estimate), `Ratios of GMT/GMC`))
    

    rrdiff_case <- rpcnt_case %>% 
      mutate(groupn = 2-match(!!as.name(subs), comp.i)%%2) %>%
      pivot_wider(id_cols = c("Visit", "Marker", "Ind", all_of(sub.by)),
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
                          by = c(sub.by, subs, "N", "Marker", "Visit")) %>%
      filter(N!=0) %>% 
      mutate(`GMT/GMC`=gsub("(NaN, NaN)", "", `GMT/GMC`)) %>% 
      # filter(!is.na(rslt) | grepl("functionality", Marker)) %>% 
      pivot_wider(id_cols = c(all_of(sub.by), "Marker", "Visit"),
                  names_from = all_of(subs),
                  values_from = c(N, rslt, `GMT/GMC`)) %>%
      full_join(rrdiff_case, by = c(sub.by, "Marker", "Visit")) %>%
      full_join(rgmt_case, by = c(sub.by, "Marker", "Visit")) %>% 
      filter(!is.na(N_Cases))
      
    
    
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
      mutate(Visit=factor(Visit, levels=labels.time)) %>% 
      arrange(Trti, Visit, Marker) %>% 
      select(Arms=Trti, Visit, Marker, 
             paste(c("N", "rslt", "GMT/GMC"), comp.i[2], sep="_"),
             paste(c("N", "rslt", "GMT/GMC"), comp.i[1], sep="_"), 
             rrdiff, `Ratios of GMT/GMC`, all_of(sub.by.add)) 
    
    tab_gmtr <- tab_case %>% 
      filter(grepl("over", Visit)) %>% 
      select(Arms, Visit, Marker, 
             paste(c("N", "GMT/GMC"), comp.i[2], sep="_"),
             paste(c("N", "GMT/GMC"), comp.i[1], sep="_"),
             `Ratios of GMT/GMC`, 
             all_of(sub.by.add)) %>% 
      arrange(Arms, Visit, Marker)
    
    tab_case <- tab_case %>% 
      filter(!grepl("over", Visit)) %>% 
      arrange(Arms, Visit, Marker)
    
    
    return(list(tab_case=tab_case, tab_gmtr=tab_gmtr))
  }
  
  
tab_assay <- lapply(1:Trtn, tab_assay_list, subs="Case", comp.i=c("Non-Cases", "Cases"))


tab_case <- tab_assay %>% 
  map("tab_case") %>% 
  bind_rows()

tab_gmtr <- tab_assay %>% 
  map("tab_gmtr") %>% 
  bind_rows()

tab_assay_status <- lapply(1:Trtn, tab_assay_list, subs="Case", sub.by.add="Naive", comp.i=c("Non-Cases", "Cases"))


tab_case_naive <- tab_assay_status %>% 
  map("tab_case") %>% 
  bind_rows() %>% 
  filter(Naive=="Naive") %>% 
  select(-Naive)

tab_case_nnaive <- tab_assay_status %>% 
  map("tab_case") %>% 
  bind_rows() %>% 
  filter(Naive=="Non-naive") %>% 
  select(-Naive)

tab_gmtr_naive <- tab_assay_status %>% 
  map("tab_gmtr") %>% 
  bind_rows() %>% 
  filter(Naive=="Naive") %>% 
  select(-Naive)

tab_gmtr_nnaive <- tab_assay_status %>% 
  map("tab_gmtr") %>% 
  bind_rows() %>% 
  filter(Naive=="Non-naive") %>% 
  select(-Naive)


print("Done with all tables")

# path for tables
save.results.to <- here::here("output")
if (!dir.exists(save.results.to))  dir.create(save.results.to)
save.results.to <- paste0(here::here("output"), "/", attr(config,"config"))

if (!dir.exists(save.results.to))  dir.create(save.results.to)
print(paste0("save.results.to equals ", save.results.to))

save(list = c("tlf", names(tlf)), file = file.path(save.results.to, sprintf("Tables%s.Rdata", ifelse(exists("COR"), gsub("D15", "D22", COR), ""))))
