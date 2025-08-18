# Sys.setenv(TRIAL = "hvtn705second")
# Sys.setenv(TRIAL = "moderna_real")
# Sys.setenv(TRIAL = "janssen_pooled_partA")
# Sys.setenv(TRIAL = "covail_tcell")
# # COR = "D15to91covail_tcell"
# COR = "D15to181covail_tcell"
#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))

# There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))

source(here::here("..", "_common.R"))
#-----------------------------------------------
## ----load-all-SLobjects, message=FALSE, error=FALSE, warning=FALSE----------------------------------------------------------------------------------------------------

library("conflicted")
library("tidyr")
library("purrr")
library("dplyr")
library("cowplot")
library("ggplot2")
library("vimp")
library("kyotil")
library("cvAUC")
library(gridExtra)
library(forcats)
library(cowplot)
library(here)
# filter() is used in both dplyr and stats, so need to set the preference to dplyr
conflict_prefer("filter", "dplyr")
conflict_prefer("summarise", "dplyr")
conflict_prefer("load", "base")
source(here("code", "utils.R"))
method <- "method.CC_nloglik" # since SuperLearner relies on this to be in GlobalEnv
ggplot2::theme_set(theme_cowplot())
load(here(file_path, "objects_for_running_SL.rda"))

# # Get vimps in one dataframe: ----------------------------------------------------------
# pooled_ests_lst <- list.files(here(paste0("output/", Sys.getenv("TRIAL"))), pattern = "pooled_ests_*") %>%
#     tibble(file = .) %>%
#     mutate(listdat = lapply(paste0("output/", Sys.getenv("TRIAL"), "/", file), readRDS))
# 
# all_estimates = as.data.frame(do.call("rbind", pooled_ests_lst$listdat))
# 
# # add on the variable set name
# vim_estimates <- all_estimates %>%
#     mutate(variable_set = rep(varset_names, each = 2), .before = "s")
# 
# # save the output
# saveRDS(vim_estimates, file = paste0("output/", Sys.getenv("TRIAL"), "/vim_estimates.rds"))
# # ----------------------------------------------------------------------------------------
# read in the results; note that createRDAfiles_fromSLobjects has to be run prior to this

if (study_name %in% c("COVE", "MockCOVE")) {
  cvaucs_vacc <- readRDS(file = here::here("output", Sys.getenv("TRIAL"), "cvaucs_vacc_EventIndPrimaryD57.rds"))
  vim_estimates <- readRDS(file = here::here("output", Sys.getenv("TRIAL"), "vim_estimates.rds"))  %>%
    mutate(group = ifelse(variable_set %in% c("2_bAbSpike_D57", "3_bAbRBD_D57", "13_bAbSpike_D29"), TRUE, group))
}
if (study_name == "HVTN705") {
  cvaucs_vacc <- readRDS(file = here::here("output", Sys.getenv("TRIAL"), "cvaucs_vacc_EventIndPrimaryD210.rds"))
  vim_estimates <- readRDS(file = here::here("output", Sys.getenv("TRIAL"), "vim_estimates.rds")) %>%
    mutate(group = ifelse(variable_set %in% c("3_M7Tcells", "36_M7_2_3_4_5_12_13"), TRUE, group))
}
if (study_name == "ENSEMBLE") {
  cvaucs_vacc <- readRDS(file = here::here("output", Sys.getenv("TRIAL"), "cvaucs_vacc_EventIndPrimaryIncludeNotMolecConfirmedD29.rds"))
  vim_estimates <- readRDS(file = here::here("output", Sys.getenv("TRIAL"), "vim_estimates.rds")) 
}
if (study_name == "COVAIL") {
  if(COR %in% c("D15to91covail_tcell", "D15to91covail_xassays")){
    cvaucs_vacc <- readRDS(file = here(file_path, "cvaucs_vacc_COVIDIndD22toD91.rds"))
    #vim_estimates <- readRDS(file = here::here("output", Sys.getenv("TRIAL"), "vim_estimates.rds")) 
  } else if(COR %in% c("D15to181covail_tcell", "D15to181covail_xassays")){
    cvaucs_vacc <- readRDS(file = here(file_path, "cvaucs_vacc_COVIDIndD22toD181.rds"))
    #vim_estimates <- readRDS(file = here::here("output", Sys.getenv("TRIAL"), "vim_estimates.rds")) 
  }
}

ph2_vacc_ptids <- readRDS(file = here(file_path, "ph2_vacc_ptids.rds"))

# Select the random seed from which to display results
if(study_name %in% c("COVE", "MockCOVE", "HVTN705")){
  rseed = 1
} else if(study_name %in% c("ENSEMBLE", "COVAIL")){
  rseed = 2
}

# Create tables ------------------------------------------------------------------------
# Table of learner/screen combinations
# caption <- "All learner-screen combinations (16 in total) used as input to the Superlearner."

tab <- cvaucs_vacc %>%
  filter(!Learner %in% c("SL", "Discrete SL")) %>%
  select(Learner, Screen) %>%
  mutate(Screen = fct_relevel(Screen, c("all", #"glmnet", 
                                        "univar_logistic_pval",
                                        "highcor_random")),
         Learner = as.factor(Learner)) %>%
  arrange(Learner, Screen) %>%
  distinct(Learner, Screen) %>%
  rename("Screen*" = Screen)

if (!grepl("Mock", study_name) & study_name == "COVE") {
  tab <- tab %>%
    mutate(Learner = fct_relevel(Learner, c("SL.mean", "SL.glmnet.0", "SL.glmnet.1", "SL.xgboost.2.no", "SL.xgboost.4.no",
                                            "SL.xgboost.2.yes", "SL.xgboost.4.yes", "SL.ranger.yes", "SL.ranger.no", "SL.glm"))) %>%
    arrange(Learner, `Screen*`)
} else if (!grepl("Mock", study_name) & study_name == "HVTN705") {
  tab <- tab %>%
    mutate(Learner = fct_relevel(Learner, c("SL.mean", "SL.gam", #"SL.bayesglm", 
					    "SL.glm", #"SL.glm.interaction",
					    "SL.glmnet.0",
					    "SL.glmnet.1",
                                            #"SL.ksvm.polydot", "SL.ksvm.rbfdot",
                                            #"SL.polymars",
                                            "SL.xgboost.2.no",
                                            "SL.xgboost.4.no",
                                            "SL.xgboost.2.yes",
					    "SL.xgboost.4.yes",
					    "SL.ranger.no",
                                            "SL.ranger.yes"))) %>%
    arrange(Learner, `Screen*`)
} else {
  tab <- tab %>%
    mutate(Learner = fct_relevel(Learner, c("SL.mean", "SL.glm"))) %>%
    arrange(Learner, `Screen*`)
}

tab %>% write.csv(here(file_path, "learner-screens.csv"))

# Table of variable set definitions
if (study_name %in% c("COVE", "MockCOVE")) {
  caption <- "The 34 variable sets on which an estimated optimal surrogate was built."
  
  only_varsets <- varset_names[1:34]

  tab <- data.frame(`Variable Set Name` = varset_names[1:34],
                    `Variables included in the set` = c("Baseline risk factors only (Reference model)",
                                                        "Baseline risk factors + Day 57 bAb anti-Spike markers",
                                                        "Baseline risk factors + Day 57 bAb anti-RBD markers",
                                                        "Baseline risk factors + Day 57 p-nAb ID50 markers",
                                                        "Baseline risk factors + Day 57 p-nAb ID80 markers",
                                                        "Baseline risk factors + Day 57 l-nAb MN50 markers",
                                                        "Baseline risk factors + Day 57 bAb markers and p-nAb ID50 markers",
                                                        "Baseline risk factors + Day 57 bAb markers and p-nAb ID80 markers",
                                                        "Baseline risk factors + Day 57 bAb markers and l-nAb MN50 markers",
                                                        "Baseline risk factors + Day 57 bAb markers and combination scores across the five markers [PCA1, PCA2, FSDAM1/FSDAM2 (the first two
components of nonlinear PCA), and the maximum signal diversity score]",
                                                        "Baseline risk factors + all individual Day 57 marker variables",
                                                        "Baseline risk factors + all individual Day 57 marker variables and their combination scores (Full model of Day 57 markers)",

                                                        "Baseline risk factors + Day 29 bAb anti-Spike markers",
                                                        "Baseline risk factors + Day 29 bAb anti-RBD markers",
                                                        "Baseline risk factors + Day 29 p-nAb ID50 markers",
                                                        "Baseline risk factors + Day 29 p-nAb ID80 markers",
                                                        "Baseline risk factors + Day 29 l-nAb MN50 markers",
                                                        "Baseline risk factors + Day 29 bAb markers and p-nAb ID50 markers",
                                                        "Baseline risk factors + Day 29 bAb markers and p-nAb ID80 markers",
                                                        "Baseline risk factors + Day 29 bAb markers and l-nAb MN50 markers",
                                                        "Baseline risk factors + Day 29 bAb markers and combination scores across the five markers [PCA1, PCA2, FSDAM1/FSDAM2 (the first two
components of nonlinear PCA), and the maximum signal diversity score]",
                                                        "Baseline risk factors + all individual Day 29 marker variables",
                                                        "Baseline risk factors + all individual Day 29 marker variables and their combination scores (Full model of Day 29 markers)",

                                                        "Baseline risk factors + Day 29 and Day 57 bAb anti-Spike markers",
                                                        "Baseline risk factors + Day 29 and Day 57 bAb anti-RBD markers",
                                                        "Baseline risk factors + Day 29 and Day 57 p-nAb ID50 markers",
                                                        "Baseline risk factors + Day 29 and Day 57 p-nAb ID80 markers",
                                                        "Baseline risk factors + Day 29 and Day 57 l-nAb MN50 markers",
                                                        "Baseline risk factors + Day 29 and Day 57 bAb markers and p-nAb ID50 markers",
                                                        "Baseline risk factors + Day 29 and Day 57 bAb markers and p-nAb ID80 markers",
                                                        "Baseline risk factors + Day 29 and Day 57 bAb markers and l-nAb MN50 markers",
                                                        "Baseline risk factors + Day 29 and Day 57 bAb markers and combination scores across the ten markers [PCA1, PCA2, FSDAM1/FSDAM2 (the first two
components of nonlinear PCA), and the maximum signal diversity score]",
                                                        "Baseline risk factors + all individual Day 29 and Day 57 marker variables",
                                                        "Baseline risk factors + all individual Day 29 and Day 57 marker variables and their combination scores (Full model of Day 29 and Day 57 markers)"))
}
if (study_name == "HVTN705") {
  total_varsets = 36
  caption <- "The 36 variable sets on which an estimated optimal surrogate was built."

  tab <- data.frame(`Variable Set Name` = c("1_baselineFactors", "2_M7ELISA", "3_M7Tcells","4_M7ADCP","5_M7IgG3",
                                            "6_M7IgG3_gp140","7_M7IgG3_gp120","8_M7IgG3_V1V2","9_M7IgG3_gp41",
                                            "10_M7IgG3_breadthScore","11_M7IgG3_multi","12_M7IgGt_V1V2","13_M7ADCC",
                                            "14_M7_combPrimaryMarkers", "15_M7_combAllMarkers","16_M7_overallScore", 
                                            "17_M7_2_3", "18_M7_2_4", "19_M7_2_5_12", "20_M7_2_13", "21_M7_3_4", 
                                            "22_M7_3_5_12", "23_M7_3_13", "24_M7_4_5_12", "25_M7_4_13", 
                                            "26_M7_5_12_13","27_M7_2_3_4", "28_M7_2_3_5_12","29_M7_2_3_13","30_M7_3_4_5_12",
                                            "31_M7_3_4_13","32_M7_4_5_12_13","33_M7_2_3_4_5_12","34_M7_2_3_4_13",
                                            "35_M7_3_4_5_12_13","36_M7_2_3_4_5_12_13"),
                    `Variables included in the set` = c(
                      "Baseline Factors South Africa, Age, BMI, Behavioral risk score",
                      "1.+M7 ELISA ELISA markers ELISA VT-C and ELISA VT-M controlling for baseline factors (Cont-base)",
                      "1.+M7 T cells ELISpot PTE Env, CD4/CD8 T cells Cont-base",
                      "1.+M7 ADCP ADCP markers ADCP C97ZA, ADCP Mosaic Cont-base",
                      "1.+M7 IgG3 All BAMA IgG3 markers Cont-base",
                      "1.+M7 IgG3 gp140 All BAMA IgG3 gp140 markers Cont-base",
                      "1.+M7 IgG3 gp120 All BAMA IgG3 gp120 markers Cont-base",
                      "1.+M7 IgG3 V1V2 All BAMA IgG3 V1V2 markers Cont-base",
                      "1.+M7 IgG3 gp41 All BAMA IgG3 gp41 markers Cont-base",
                      "1.+M7 IgG3 breadth scores All BAMA antigen-specific breadth scores Cont-base",
                      "1.+M7 IgG3 Multi-epitope breadth BAMA IgG3 multi-epitope cross-reactivity score Cont-base",
                      "1.+M7 IgGt V1V2 All BAMA IgG V1V2 markers Cont-base",
                      "1.+M7 ADCC All ADCC markers Cont-base",
                      "1.+M7 comb scores over primary markers All combination scores derived over the 6 primary markers [PCA1, PCA2, FSDAM1/FSDAM2 (the first two components of nonlinear PCA), and the maximum signal diversity
                                                          score He and Fong (2019)]",
                      "1.+M7 comb scores over all markers All combination scores derived over all markers [PCA1, PCA2, FSDAM1/FSDAM2, and max.sig.div.score].",
                      "1.+M7 Overall scores across assays Two overall combined scores Cont-base",
                      "1.+2.+3. Do ELISA + T cells together improve classification?",
                      "1.+2.+4. Do ELISA + ADCP together improve classification?",
                      "1.+2.+5.+12. Do ELISA + IgG3/IgGt together improve classification?",
                      "1.+2.+13. Do ELISA + ADCC together improve classification?",
                      "1.+3.+4. Do T cells + ADCP together improve classification?",
                      "1.+3.+5.+12. Do T cells + IgG3/IgGt together improve classification?",
                      "1.+3.+13. Do T cells + ADCC together improve classification?",
                      "1.+4.+5.+12. Do ADCP + IgG3/IgGt together improve classification?",
                      "1.+4.+13. Do ADCP + ADCC together improve classification?",
                      "1.+5.+12.+13. Do IgG3/IgGt + ADCC together improve classification?",
                      "1.+2.+3.+4. Do ELISA + T cells + ADCP together improve classification?",
                      "1.+2+3.+5.+12. Do ELISA + T cells + IgG3/IgGt together improve classification?",
                      "1.+2+3.+13. Do ELISA + T cells + ADCC together improve classification?",
                      "1.+3.+4.+5.+12. Do T cells + ADCP + IgG3/IgGt together improve classification?",
                      "1.+3.+4.+13. Do T cells + ADCP + ADCC together improve classification?",
                      "1.+4.+5.+12.+13. Do ADCP + IgG3/IgGt + ADCC together improve classification?",
                      "1.+2.+3.+4.+5.+12. Do ELISA + T cells + ADCP + IgG3/IgGt together improve classification?",
                      "1.+2.+3.+4.+13. Do ELISA + T cells + ADCP + ADCC together improve classification?",
                      "1.+3.+4.+5.+12.+13. Do T cells + ADCP + IgG3/IgGt + ADCC together improve classification?",
                      "1.+2.+3.+4.+5.+12.+13. (All) Do all the immune marker sets together improve classification?"
                    ))

}
if (study_name %in% c("ENSEMBLE")) {
  caption <- "The 16 variable sets on which an estimated optimal surrogate was built."
  
  only_varsets <- varset_names[1:16]
  
  tab <- data.frame(`Variable Set Name` = varset_names[1:16],
                    `Variables included in the set` = c("Baseline risk factors only (Reference model)",
                                                        "Baseline risk factors + Day 29 bAb anti-Spike markers",
                                                        "Baseline risk factors + Day 29 bAb anti-RBD markers",
                                                        "Baseline risk factors + Day 29 p-nAb ID50 markers",
                                                        "Baseline risk factors + Day 29 ADCP markers",
                                                        "Baseline risk factors + Day 29 Functional markers minus bAb markers",
                                                        "Baseline risk factors + Day 29 bAb and p-nAb ID50 markers",
                                                        "Baseline risk factors + Day 29 bAb and ADCP markers",
                                                        "Baseline risk factors + Day 29 bAb and p-nAb ID50 markers and difference between them",
                                                        "Baseline risk factors + Day 29 bAb and ADCP markers and difference between them",
                                                        "Baseline risk factors + Day 29 bAb and p-nAb ID50 and ADCP markers",
                                                        "Baseline risk factors + Day 29 bAb and p-nAb ID50 and ADCP markers and functional markers minus bAb markers",
                                                        "Baseline risk factors + Day 29 bAb markers and combination scores across the four markers [PCA1, PCA2, FSDAM1/FSDAM2 (the first two
components of nonlinear PCA), and the maximum signal diversity score]",
                                                        "Baseline risk factors + Day 29 p-nAb ID50 markers and combination scores across the four markers",
                                                        "Baseline risk factors + Day 29 ADCP markers and combination scores across the four markers",
                                                        "Baseline risk factors + all individual Day 29 marker variables and their combination scores (Full model of Day 29 markers)"))
}
if (TRIAL == "covail_tcell" & non_naive == FALSE) {
  caption <- "The 13 variable sets on which an estimated optimal surrogate was built."
  
  only_varsets <- varset_names[1:13]
  
  tab <- data.frame(`Variable Set Name` = varset_names[1:13],
                    `Variables included in the set` = c(
                      "Baseline risk factors only (Reference model): includes baseline risk score and the standardized FOI score only",
                      #"Baseline risk factors only (Reference model): includes baseline risk score, standardized FOI score, insert and stage (vaccine manufacturer) info",
                                                        "Baseline risk factors + D1 antibody markers against BA.1",
                                                        "Baseline risk factors + D15 antibody markers against BA.1",
                                                        "Baseline risk factors + D1 and D15 antibody markers against BA.1",
                                                        "Baseline risk factors + D1 primary and secondary T cell markers",
                                                        "Baseline risk factors + D15 primary and secondary T cell markers",
                                                        "Baseline risk factors + D1 and D15 primary and secondary T cell markers",
                                                        "Baseline risk factors + D1 antibody markers against BA.1 + D1 primary and secondary T cell markers ",
                                                        "Baseline risk factors + D1 antibody markers against BA.1 + D15 primary and secondary T cell markers ",
                                                        "Baseline risk factors + D15 antibody markers against BA.1 + D1 primary and secondary T cell markers ",
                                                        "Baseline risk factors + D15 antibody markers against BA.1 + D15 primary and secondary T cell markers ",
                                                        "Baseline risk factors + D1 and D15 antibody markers against BA.1 + D1 and D15 primary and secondary T cell markers",
                                                        "Baseline risk factors + D1 and D15 antibody markers against BA.1 + D1 and D15 [primary, secondary, and exploratory T cell markers"))
}
if (TRIAL == "covail_tcell" & non_naive == TRUE) {
  caption <- "The 28 variable sets on which an estimated optimal surrogate was built."
  
  only_varsets <- varset_names[1:28]
  
  tab <- data.frame(`Variable Set Name` = varset_names[1:28],
                    `Variables included in the set` = c(
                      "Baseline risk factors only (Reference model): includes baseline risk score only",
                      #"Baseline risk factors only (Reference model): includes baseline risk score, standardized FOI score, insert and stage (vaccine manufacturer) info",
                      "Baseline risk factors + D1 antibody markers against BA.1",
                      "Baseline risk factors + D15 antibody markers against BA.1",
                      "Baseline risk factors + D1 and D15 antibody markers against BA.1",
                      "Baseline risk score + D1 primary CD4+ T cell marker IFN-g/IL-2 Spike BA.4/5",
                      "Baseline risk score + D1 primary CD4+ T cell marker Functionality Score (FS) Spike BA.4/5", 
                      "Baseline risk score + D1 primary CD4+ T cell marker IFN-g/IL-2 N Index", 
                      "Baseline risk score + D1 primary CD4+ T cell marker FS N Index", 
                      "Baseline risk score + D15 primary CD4+ T cell marker IFN-g/IL-2 Spike BA.4/5", 
                      "Baseline risk score + D15 primary CD4+ T cell marker FS Spike BA.4/5", 
                      "Baseline risk score + D1 primary CD4+ T cell marker IFN-g/IL-2 Spike BA.4/5 + D15 primary CD4+ T cell marker IFN-g/IL-2 Spike BA.4/5", 
                      "Baseline risk score + D1 primary CD4+ T cell marker FS Spike BA.4/5 + D15 primary CD4+ T cell marker FS Spike BA.4/5",
                      "Baseline risk score + D1 primary CD4+ T cell marker IFN-g/IL-2 N Index + D15 primary CD4+ T cell marker IFN-g/IL-2 Spike BA.4/5",
                      "Baseline risk score + D1 primary CD4+ T cell marker FS N Index + D15 primary CD4+ T cell marker FS Spike BA.4/5",
                      "Baseline risk score + D1 antibody markers against BA.1 + D1 primary CD4+ T cell marker IFN-g/IL-2 Spike BA.4/5",
                      "Baseline risk score + D1 antibody markers against BA.1 + D1 primary CD4+ T cell marker FS Spike BA.4/5",
                      "Baseline risk score + D1 antibody markers against BA.1 + D1 primary CD4+ T cell marker IFN-g/IL-2 N Index",
                      "Baseline risk score + D1 antibody markers against BA.1 + D1 primary CD4+ T cell marker FS N Index",
                      "Baseline risk score + D1 antibody markers against BA.1 + D15 primary CD4+ T cell marker IFN-g/IL-2 Spike BA.4/5",
                      "Baseline risk score + D1 antibody markers against BA.1 + D15 primary CD4+ T cell marker FS Spike BA.4/5",
                      "Baseline risk score + D15 antibody markers against BA.1 + D1 primary CD4+ T cell marker IFN-g/IL-2 Spike BA.4/5",
                      "Baseline risk score + D15 antibody markers against BA.1 + D1 primary CD4+ T cell marker FS Spike BA.4/5",
                      "Baseline risk score + D15 antibody markers against BA.1 + D1 primary CD4+ T cell marker IFN-g/IL-2 N Index",
                      "Baseline risk score + D15 antibody markers against BA.1 + D1 primary CD4+ T cell marker FS N Index",
                      "Baseline risk score + D15 antibody markers against BA.1 + D15 primary CD4+ T cell marker IFN-g/IL-2 Spike BA.4/5",
                      "Baseline risk score + D15 antibody markers against BA.1 + D15 primary CD4+ T cell marker FS Spike BA.4/5",
                      "Baseline risk score + D15 antibody markers against BA.1 + D15 primary CD4+ T cell marker IFN-g/IL-2 N Index",
                      "Baseline risk score + D15 antibody markers against BA.1 + D15 primary CD4+ T cell marker FS N Index"))
}
if (TRIAL == "covail_xassays" & non_naive == FALSE) {
  #caption <- "The 82 variable sets on which an estimated optimal surrogate was built."
  
  only_varsets <- varset_names[1:82]
  
  tab <- data.frame(`Variable Set Name` = varset_names[1:82],
                    `Variables included in the set` = c(
                      # Dynamic first entry
                      if (setequal(briskfactors, "risk_score")) {
                        "Baseline risk factors only (Reference model): includes baseline risk score only"
                      } else if (all(c("risk_score", "FOIstandardized") %in% briskfactors) &&
                                 !any(grepl("stage|Trtgrp", briskfactors))) {
                        "Baseline risk factors only (Reference model): includes baseline risk score and standardized FOI score only"
                      } else if (all(c("risk_score", "FOIstandardized") %in% briskfactors) &&
                                 any(grepl("stage|Trtgrp", briskfactors))) {
                        "Baseline risk factors only (Reference model): includes baseline risk score, standardized FOI score, insert and stage (vaccine manufacturer) info"
                      },
                      
                      #"Baseline risk factors only (Reference model): includes baseline risk score and the standardized FOI score only",
                      # "Baseline risk factors only (Reference model): includes baseline risk score, standardized FOI score, insert and stage (vaccine manufacturer) info",
                      "Baseline factors + D1 nAb titers set",
                      "Baseline factors + D15 nAb titers set",
                      "Baseline factors + D1 and D15 nAb titers set",
                      "Baseline factors + D1 CD4 T cell set",
                      "Baseline factors + D15 CD4 T cell set",
                      "Baseline factors + D1 and D15 CD4 T cell set",
                      "Baseline factors + D1 CD8 T cell set",
                      "Baseline factors + D15 CD8 T cell set",
                      "Baseline factors + D1 and D15 CD8 T cell set",
                      "Baseline factors + D1 Fc-spike set",
                      "Baseline factors + D15 Fc-spike set",
                      "Baseline factors + D1 and D15 Fc-Spike set",
                      "Baseline factors + D1 bAb-spike set",
                      "Baseline factors + D15 bAb-spike set",
                      "Baseline factors + D1 and D15 bAb-Spike set",
                      "Baseline factors + D1 combination markers set",
                      "Baseline factors + D15 combination markers set",
                      "Baseline factors + D1 and D15 combination markers set",
                      "Baseline factors + D1 nAb titers set + D1 CD4 T cell set",
                      "Baseline factors + D1 nAb titers set + D1 CD8 T cell set",
                      "Baseline factors + D1 nAb titers set + D1 Fc-spike set",
                      "Baseline factors + D1 nAb titers set + D1 bAb-spike set",
                      "Baseline factors + D1 nAb titers set + D1 combination markers set",
                      "Baseline factors + D1 nAb titers set + D15 CD4 T cell set",
                      "Baseline factors + D1 nAb titers set + D15 CD8 T cell set",
                      "Baseline factors + D1 nAb titers set + D15 Fc-spike set",
                      "Baseline factors + D1 nAb titers set + D15 bAb-spike set",
                      "Baseline factors + D1 nAb titers set + D15 combination markers set",
                      "Baseline factors + D15 nAb titers set + D1 CD4 T cell set",
                      "Baseline factors + D15 nAb titers set + D1 CD8 T cell set",
                      "Baseline factors + D15 nAb titers set + D1 Fc-spike set",
                      "Baseline factors + D15 nAb titers set + D1 bAb-spike set",
                      "Baseline factors + D15 nAb titers set + D1 combination markers set",
                      "Baseline factors + D15 nAb titers set + D15 CD4 T cell set",
                      "Baseline factors + D15 nAb titers set + D15 CD8 T cell set",
                      "Baseline factors + D15 nAb titers set + D15 Fc-spike set",
                      "Baseline factors + D15 nAb titers set + D15 bAb-spike set",
                      "Baseline factors + D15 nAb titers set + D15 combination markers set",
                      "Baseline factors + D1 CD4 T cell set + D1 CD8 T cell set",
                      "Baseline factors + D1 CD4 T cell set + D1 Fc-spike set",
                      "Baseline factors + D1 CD4 T cell set + D1 bAb-spike set",
                      "Baseline factors + D1 CD4 T cell set + D1 combination markers set",
                      "Baseline factors + D1 CD4 T cell set + D15 CD8 T cell set",
                      "Baseline factors + D1 CD4 T cell set + D15 Fc-spike set",
                      "Baseline factors + D1 CD4 T cell set + D15 bAb-spike set",
                      "Baseline factors + D1 CD4 T cell set + D15 combination markers set",
                      "Baseline factors + D15 CD4 T cell set + D1 CD8 T cell set",
                      "Baseline factors + D15 CD4 T cell set + D1 Fc-spike set",
                      "Baseline factors + D15 CD4 T cell set + D1 bAb-spike set",
                      "Baseline factors + D15 CD4 T cell set + D1 combination markers set",
                      "Baseline factors + D15 CD4 T cell set + D15 CD8 T cell set",
                      "Baseline factors + D15 CD4 T cell set + D15 Fc-spike set",
                      "Baseline factors + D15 CD4 T cell set + D15 bAb-spike set",
                      "Baseline factors + D15 CD4 T cell set + D15 combination markers set",
                      "Baseline factors + D1 CD8 T cell set + D1 Fc-spike set",
                      "Baseline factors + D1 CD8 T cell set + D1 bAb-spike set",
                      "Baseline factors + D1 CD8 T cell set + D1 combination markers set",
                      "Baseline factors + D1 CD8 T cell set + D15 Fc-spike set",
                      "Baseline factors + D1 CD8 T cell set + D15 bAb-spike set",
                      "Baseline factors + D1 CD8 T cell set + D15 combination markers set",
                      "Baseline factors + D15 CD8 T cell set + D1 Fc-spike set",
                      "Baseline factors + D15 CD8 T cell set + D1 bAb-spike set",
                      "Baseline factors + D15 CD8 T cell set + D1 combination markers set",
                      "Baseline factors + D15 CD8 T cell set + D15 Fc-spike set",
                      "Baseline factors + D15 CD8 T cell set + D15 bAb-spike set",
                      "Baseline factors + D15 CD8 T cell set + D15 combination markers set",
                      "Baseline factors + D1 Fc-spike set + D1 bAb-spike set",
                      "Baseline factors + D1 Fc-spike set + D1 combination markers set",
                      "Baseline factors + D1 Fc-spike set + D15 bAb-spike set",
                      "Baseline factors + D1 Fc-spike set + D15 combination markers set",
                      "Baseline factors + D15 Fc-spike set + D1 bAb-spike set",
                      "Baseline factors + D15 Fc-spike set + D1 combination markers set",
                      "Baseline factors + D15 Fc-spike set + D15 bAb-spike set",
                      "Baseline factors + D15 Fc-spike set + D15 combination markers set",
                      "Baseline factors + D1 bAb-spike set + D1 combination markers set",
                      "Baseline factors + D1 bAb-spike set + D15 combination markers set",
                      "Baseline factors + D15 bAb-spike set + D1 combination markers set",
                      "Baseline factors + D15 bAb-spike set + D15 combination markers set",
                      "Baseline factors + All eight immunoassay sets at D1",
                      "Baseline factors + All six immunoassay sets at D15 (always excluding D15 anti-N markers)",
                      "Baseline factors + All six immunoassay sets at both D1 and at D15 (always excluding D15 anti-N markers)"))
} else if (TRIAL == "covail_xassays" & non_naive == TRUE) {
  #caption <- "The 109 variable sets on which an estimated optimal surrogate was built."
  
  only_varsets <- varset_names[1:109]
  
  tab <- data.frame(`Variable Set Name` = varset_names[1:109],
                    `Variables included in the set` = c(
                      # Dynamic first entry
                      if (setequal(briskfactors, "risk_score")) {
                        "Baseline risk factors only (Reference model): includes baseline risk score only"
                      } else if (all(c("risk_score", "FOIstandardized") %in% briskfactors) &&
                                 !any(grepl("stage|Trtgrp", briskfactors))) {
                        "Baseline risk factors only (Reference model): includes baseline risk score and standardized FOI score only"
                      } else if (all(c("risk_score", "FOIstandardized") %in% briskfactors) &&
                                 any(grepl("stage|Trtgrp", briskfactors))) {
                        "Baseline risk factors only (Reference model): includes baseline risk score, standardized FOI score, insert and stage (vaccine manufacturer) info"
                      },
                      
                      #"Baseline risk factors only (Reference model): includes baseline risk score and the standardized FOI score only",
                      #"Baseline risk factors only (Reference model): includes baseline risk score, standardized FOI score, insert and stage (vaccine manufacturer) info",
                      "Baseline factors + D1 nAb titers set",
                      "Baseline factors + D15 nAb titers set",
                      "Baseline factors + D1 and D15 nAb titers set",
                      "Baseline factors + D1 CD4 T cell set",
                      "Baseline factors + D15 CD4 T cell set",
                      "Baseline factors + D1 and D15 CD4 T cell set",
                      "Baseline factors + D1 CD8 T cell set",
                      "Baseline factors + D15 CD8 T cell set",
                      "Baseline factors + D1 and D15 CD8 T cell set",
                      "Baseline factors + D1 Fc-spike set",
                      "Baseline factors + D15 Fc-spike set",
                      "Baseline factors + D1 and D15 Fc-Spike set",
                      "Baseline factors + D1 bAb-spike set",
                      "Baseline factors + D15 bAb-spike set",
                      "Baseline factors + D1 and D15 bAb-Spike set",
                        "Baseline factors + D1 Fc-N set (only included in non-naive cohort)",
                        "Baseline factors + D1 bAb-N set (only included in non-naive cohort)",
                        "Baseline factors + D1 combination markers set",
                        "Baseline factors + D15 combination markers set",
                        "Baseline factors + D1 and D15 combination markers set",
                        "Baseline factors + D1 nAb titers set + D1 CD4 T cell set",
                        "Baseline factors + D1 nAb titers set + D1 CD8 T cell set",
                        "Baseline factors + D1 nAb titers set + D1 Fc-spike set",
                        "Baseline factors + D1 nAb titers set + D1 bAb-spike set",
                        "Baseline factors + D1 nAb titers set + D1 Fc-N set (this variable set only included for the non-naive cohort)",
                        "Baseline factors + D1 nAb titers set + D1 bAb-N set (this variable set only included for the non-naive cohort)",
                        "Baseline factors + D1 nAb titers set + D1 combination markers set (combination markers include anti-N markers only in non-naïve cohort)",
                        "Baseline factors + D1 nAb titers set + D15 CD4 T cell set",
                        "Baseline factors + D1 nAb titers set + D15 CD8 T cell set",
                        "Baseline factors + D1 nAb titers set + D15 Fc-spike set",
                        "Baseline factors + D1 nAb titers set + D15 bAb-spike set",
                        "Baseline factors + D1 nAb titers set + D15 combination markers set",
                        "Baseline factors + D15 nAb titers set + D1 CD4 T cell set",
                        "Baseline factors + D15 nAb titers set + D1 CD8 T cell set",
                        "Baseline factors + D15 nAb titers set + D1 Fc-spike set",
                        "Baseline factors + D15 nAb titers set + D1 bAb-spike set",
                        "Baseline factors + D15 nAb titers set + D1 Fc-N set (this variable set only included for the non-naive cohort)",
                        "Baseline factors + D15 nAb titers set + D1 bAb-N set (this variable set only included for the non-naive cohort)",
                        "Baseline factors + D15 nAb titers set + D1 combination markers set (combination markers include anti-N markers only in non-naïve cohort)",
                        "Baseline factors + D15 nAb titers set + D15 CD4 T cell set",
                        "Baseline factors + D15 nAb titers set + D15 CD8 T cell set",
                        "Baseline factors + D15 nAb titers set + D15 Fc-spike set",
                        "Baseline factors + D15 nAb titers set + D15 bAb-spike set",
                        "Baseline factors + D15 nAb titers set + D15 combination markers set",
                        "Baseline factors + D1 CD4 T cell set + D1 CD8 T cell set",
                        "Baseline factors + D1 CD4 T cell set + D1 Fc-spike set",
                        "Baseline factors + D1 CD4 T cell set + D1 bAb-spike set",
                        "Baseline factors + D1 CD4 T cell set + D1 Fc-N set (this variable set only included for the non-naive cohort)",
                        "Baseline factors + D1 CD4 T cell set + D1 bAb-N set (this variable set only included for the non-naive cohort)",
                        "Baseline factors + D1 CD4 T cell set + D1 combination markers set (combination markers include anti-N markers only in non-naïve cohort)",
                        "Baseline factors + D1 CD4 T cell set + D15 CD8 T cell set",
                        "Baseline factors + D1 CD4 T cell set + D15 Fc-spike set",
                        "Baseline factors + D1 CD4 T cell set + D15 bAb-spike set",
                        "Baseline factors + D1 CD4 T cell set + D15 combination markers set",
                        "Baseline factors + D15 CD4 T cell set + D1 CD8 T cell set",
                        "Baseline factors + D15 CD4 T cell set + D1 Fc-spike set",
                        "Baseline factors + D15 CD4 T cell set + D1 bAb-spike set",
                        "Baseline factors + D15 CD4 T cell set + D1 Fc-N set (this variable set only included for the non-naive cohort)",
                        "Baseline factors + D15 CD4 T cell set + D1 bAb-N set (this variable set only included for the non-naive cohort)",
                        "Baseline factors + D15 CD4 T cell set + D1 combination markers set (combination markers include anti-N markers only in non-naïve cohort)",
                        "Baseline factors + D15 CD4 T cell set + D15 CD8 T cell set",
                        "Baseline factors + D15 CD4 T cell set + D15 Fc-spike set",
                        "Baseline factors + D15 CD4 T cell set + D15 bAb-spike set",
                        "Baseline factors + D15 CD4 T cell set + D15 combination markers set",
                        "Baseline factors + D1 CD8 T cell set + D1 Fc-spike set",
                        "Baseline factors + D1 CD8 T cell set + D1 bAb-spike set",
                        "Baseline factors + D1 CD8 T cell set + D1 Fc-N set (this variable set only included for the non-naive cohort)",
                        "Baseline factors + D1 CD8 T cell set + D1 bAb-N set (this variable set only included for the non-naive cohort)",
                        "Baseline factors + D1 CD8 T cell set + D1 combination markers set (combination markers include anti-N markers only in non-naïve cohort)",
                        "Baseline factors + D1 CD8 T cell set + D15 Fc-spike set",
                        "Baseline factors + D1 CD8 T cell set + D15 bAb-spike set",
                        "Baseline factors + D1 CD8 T cell set + D15 combination markers set",
                        "Baseline factors + D15 CD8 T cell set + D15 Fc-spike set",
                        "Baseline factors + D15 CD8 T cell set + D15 bAb-spike set",
                        "Baseline factors + D15 CD8 T cell set + D15 combination markers set",
                        "Baseline factors + D15 CD8 T cell set + D1 Fc-spike set",
                        "Baseline factors + D15 CD8 T cell set + D1 bAb-spike set",
                        "Baseline factors + D15 CD8 T cell set + D1 Fc-N set (this variable set only included for the non-naive cohort)",
                        "Baseline factors + D15 CD8 T cell set + D1 bAb-N set (this variable set only included for the non-naive cohort)",
                        "Baseline factors + D15 CD8 T cell set + D1 combination markers set combination markers include anti-N markers only in non-naïve cohort)",
                        "Baseline factors + D1 Fc-spike set + D1 bAb-spike set",
                        "Baseline factors + D1 Fc-spike set + D1 Fc-N set (this variable set only included for the non-naive cohort)",
                        "Baseline factors + D1 Fc-spike set + D1 bAb-N set (this variable set only included for the non-naive cohort)",
                        "Baseline factors + D1 Fc-spike set + D1 combination markers set (combination markers include anti-N markers only in non-naïve cohort)",
                        "Baseline factors + D1 Fc-spike set + D15 bAb-spike set",
                        "Baseline factors + D1 Fc-spike set + D15 combination markers set",
                        "Baseline factors + D15 Fc-spike set + D1 bAb-spike set",
                        "Baseline factors + D15 Fc-spike set + D1 Fc-N set (this variable set only included for the non-naive cohort)",
                        "Baseline factors + D15 Fc-spike set + D1 bAb-N set (this variable set only included for the non-naive cohort)",
                        "Baseline factors + D15 Fc-spike set + D1 combination markers set (combination markers include anti-N markers only in non-naïve cohort)",
                        "Baseline factors + D15 Fc-spike set + D15 bAb-spike set",
                        "Baseline factors + D15 Fc-spike set + D15 combination markers set",
                        "Baseline factors + D1 bAb-spike set + D1 Fc-N set (this variable set only included for the non-naive cohort)",
                        "Baseline factors + D1 bAb-spike set + D1 bAb-N set (this variable set only included for the non-naive cohort)",
                        "Baseline factors + D1 bAb-spike set + D1 combination markers set",
                        "Baseline factors + D1 bAb-spike set + D15 combination markers set",
                        "Baseline factors + D15 bAb-spike set + D1 Fc-N set (this variable set only included for the non-naive cohort)",
                        "Baseline factors + D15 bAb-spike set + D1 bAb-N set (this variable set only included for the non-naive cohort)",
                        "Baseline factors + D15 bAb-spike set + D1 combination markers set (combination markers include anti-N markers only in non-naïve cohort)",
                        "Baseline factors + D15 bAb-spike set + D15 combination markers set",
                        "Baseline factors + D1 Fc-N set + D1 bAb-N set (this variable set only included for the non-naive cohort)",
                        "Baseline factors + D1 Fc-N set + D1 combination markers set (this variable set only included for the non-naive cohort)",
                        "Baseline factors + D1 Fc-N set + D15 combination markers set (this variable set only included for the non-naive cohort)",
                        "Baseline factors + D1 bAb-N set + D1 combination markers set (this variable set only included for the non-naive cohort)",
                        "Baseline factors + D1 bAb-N set + D15 combination markers set (this variable set only included for the non-naive cohort)",
                        "Baseline factors + All eight immunoassay sets at D1 (markers including anti-N markers are only included in non-naïve cohort)",
                        "Baseline factors + All eight immunoassay sets at D15 (always excluding D15 anti-N markers)",
                        "Baseline factors + All eight immunoassay sets at both D1 and at D15 (markers including anti-N markers are only included in non-naïve cohort , and always excluding D15 anti-N markers for both the naïve and non-naïve cohorts)"))
}


tab %>% write.csv(here(file_path, "varsets.csv"))

# Create figures ---------------------------------------------------------------
# Forest Plots
options(bitmapType = "cairo")
# All Superlearners
learner.choice = "SL"
png(file = here(file_path, "figs", paste0("forest_vacc_cvaucs_all", learner.choice, "s.png")), width=1100, height=1100)
top_learner <- make_forest_plot_SL_allVarSets(cvaucs_vacc %>% filter(!is.na(varsetNo)), 
                                              varsets_to_display = 50,
                                              learner.choice,
                                              learner.plot.margin = unit(c(1.0,0.2,0.8,-0.15),"cm"), # Adjusts trbl for forest plot
                                              names.plot.margin = unit(c(-1.5,-0.15,0.45,-0.05),"cm"), # Adjusts trbl for learner-screen names plot
                                              y_at = 51.7) # Add y-coordinate at which the header for names needs to be
grid.arrange(top_learner$top_learner_nms_plot, top_learner$top_learner_plot, ncol=2)
dev.off()

# All discrete.SLs
learner.choice = "Discrete SL"
png(file = here(file_path, "figs", paste0("forest_vacc_cvaucs_all", learner.choice, "s.png")), width=1100, height=1100)
top_learner <- make_forest_plot_SL_allVarSets(cvaucs_vacc %>% filter(!is.na(varsetNo)), 
                                              varsets_to_display = 50,
                                              learner.choice,
                                              learner.plot.margin = unit(c(1.0,0.2,0.8,-0.15),"cm"), # Adjusts trbl for forest plot
                                              names.plot.margin = unit(c(-1.5,-0.15,0.45,-0.05),"cm"), # Adjusts trbl for learner-screen names plot
                                              y_at = 51.7) # Add y-coordinate at which the header for names needs to be
grid.arrange(top_learner$top_learner_nms_plot, top_learner$top_learner_plot, ncol=2)
dev.off()

# Get SL & Discrete SL performance for all variable sets
cvaucs_vacc %>% filter(!is.na(varsetNo)) %>%
  arrange(-AUC) %>% filter(Learner == "SL") %>%
  select(varset, AUCstr) %>%
  write.csv(here(file_path, "SLperformance_allvarsets.csv"))

cvaucs_vacc %>% filter(!is.na(varsetNo)) %>%
  arrange(-AUC) %>% filter(Learner == "Discrete SL") %>%
  select(varset, AUCstr) %>%
  write.csv(here(file_path, "DiscreteSLperformance_allvarsets.csv"))

# Based of Discrete SL performance, choose top 5 varsets and the reference model varset
# Plot forest plots, ROC and pred prob plots only for these 6 varsets. 
cvaucs_vacc %>% filter(!is.na(varsetNo)) %>%
  arrange(-AUC) %>% filter(Learner == "Discrete SL")


top5_plus_baseline <- cvaucs_vacc %>%
  filter(!is.na(varsetNo), Learner == "Discrete SL") %>%
  arrange(desc(AUC)) %>%
  slice_head(n = 5) %>%
  bind_rows(
    cvaucs_vacc %>%
      filter(varset == "1_baselineRiskFactors", Learner == "Discrete SL") 
  ) 



# # Forest plots for individual vaccine models
# for(i in 1:(cvaucs_vacc %>% filter(!is.na(varsetNo)) %>% distinct(varset) %>% nrow())) {
#   variableSet = unique(cvaucs_vacc$varset)[i]
#   png(file = here("figs", Sys.getenv("TRIAL"), paste0("forest_vacc_cvaucs_", variableSet, ".png")), width=1000, height=1100)
#   top_learner <- make_forest_plot(cvaucs_vacc %>% filter(varset==variableSet),
#                                   PLOT.MARGIN = unit(c(3.9,-0.15,1,-0.15),"cm"), # Adjusts trbl for forest plot
#                                   NAMES.PLOT.MARGIN = unit(c(2.01,-0.15,1.0,-0.15),"cm"), # Adjusts trbl for learner-screen names plot
#                                   y_at = 31) # Add y-coordinate at which the header for names needs to be
#   grid.arrange(top_learner$top_learner_nms_plot, top_learner$top_learner_plot, ncol=2)
#   dev.off()
# }

# Forest plots for individual vaccine models
if(basename(file_path) != "nonnaive_1dosemRNA_allcases_briskscore"){
  for(i in 1:6) {
    variableSet = top5_plus_baseline$varset[i]
    png(file = here(file_path, "figs", paste0("forest_vacc_cvaucs_", variableSet, ".png")), width=1000, height=1100)
    top_learner <- make_forest_plot(cvaucs_vacc %>% filter(varset==variableSet),
                                    PLOT.MARGIN = unit(c(3.9,-0.15,1,-0.15),"cm"), # Adjusts trbl for forest plot
                                    NAMES.PLOT.MARGIN = unit(c(2.01,-0.15,1.0,-0.15),"cm"), # Adjusts trbl for learner-screen names plot
                                    y_at = 31) # Add y-coordinate at which the header for names needs to be
    grid.arrange(top_learner$top_learner_nms_plot, top_learner$top_learner_plot, ncol=2)
    dev.off()
  }
} else {
  for(i in 1:6) {
    variableSet = top5_plus_baseline$varset[i]
    png(file = here(file_path, "figs", paste0("forest_vacc_cvaucs_", variableSet, ".png")), width=1000, height=1100)
    top_learner <- make_forest_plot(cvaucs_vacc %>% filter(varset==variableSet),
                                    PLOT.MARGIN = unit(c(2.2,-0.15,1,-0.15),"cm"), # Adjusts trbl for forest plot
                                    NAMES.PLOT.MARGIN = unit(c(0.01,-0.15,1.2,-0.15),"cm"), # Adjusts trbl for learner-screen names plot
                                    y_at = 21) # Add y-coordinate at which the header for names needs to be
    grid.arrange(top_learner$top_learner_nms_plot, top_learner$top_learner_plot, ncol=2)
    dev.off()
  }
}




#################################################################################################################################
#################################################################################################################################

# # plot ROC curve and pred.Prob with SL, Discrete SL and top 2 best-performing individual Learners for all 12 variable sets
# for(i in 1:(cvaucs_vacc %>% filter(!is.na(varsetNo)) %>% distinct(varset) %>% nrow())) {
#   variableSet = unique(cvaucs_vacc$varset)[i]
#   print(paste0("Processing ", variableSet))
#   dat <- cvaucs_vacc %>% filter(varset==variableSet)
# 
#   top2 <- bind_rows(
#     dat %>%
#       arrange(-AUC) %>%
#       filter(!Learner %in% c("SL", "Discrete SL")) %>%
#       dplyr::slice(1:2),
#     dat %>%
#       filter(Learner == "SL"),
#     dat %>%
#       filter(Learner == "Discrete SL")
#   ) %>%
#     mutate(LearnerScreen = ifelse(Learner == "SL", "Super Learner",
#                                   ifelse(Learner == "Discrete SL", Learner,
#                                          paste0(Learner, "_", Screen_fromRun))))
# 
#   # Get cvsl fit and extract cv predictions
#   if(study_name %in% c("COVE", "MockCOVE")){
#     cvfit <- readRDS(file = here("output", Sys.getenv("TRIAL"), paste0("CVSLfits_vacc_EventIndPrimaryD57_", variableSet, ".rds")))
#   } else if(study_name == "HVTN705"){
#     cvfit<- readRDS(file = here("output", Sys.getenv("TRIAL"), paste0("CVSLfits_vacc_Delta.D210_", variableSet, ".rds")))
#   } else if(study_name == "ENSEMBLE"){
#     cvfit<- readRDS(file = here("output", Sys.getenv("TRIAL"), paste0("CVSLfits_vacc_EventIndPrimaryIncludeNotMolecConfirmedD29_", variableSet, ".rds")))
#   } else if(study_name == "COVAIL"){
#     
#     if(COR == "D15to91covail_tcell"){
#       cvfit<- readRDS(file = here("output", Sys.getenv("TRIAL"), paste0("CVSLfits_vacc_COVIDIndD22toD91_", variableSet, ".rds")))
#     } else if(COR == "D15to181covail_tcell"){
#       cvfit<- readRDS(file = here("output", Sys.getenv("TRIAL"), paste0("CVSLfits_vacc_COVIDIndD22toD181_", variableSet, ".rds")))
#     } else if(COR == "D15to91covail_xassays"){
#       cvfit<- readRDS(file = here("output", Sys.getenv("TRIAL"), paste0("CVSLfits_vacc_COVIDIndD22toD91_", variableSet, ".rds")))
#     } else if(COR == "D15to181covail_xassays"){
#       cvfit<- readRDS(file = here("output", Sys.getenv("TRIAL"), paste0("CVSLfits_vacc_COVIDIndD22toD181_", variableSet, ".rds")))
#     }
#     
#   }
#   
#   pred <- get_cv_predictions(cv_fit = cvfit[[rseed]], cvaucDAT = top2, markerDAT = NULL)
#   # #Take average of predictions from the 10 random seeds
#   # pred <- get_cv_predictions(cv_fit = cvfits[[1]], cvaucDAT = top2) %>% rename(pred1 = pred) %>%
#   #   bind_cols(get_cv_predictions(cv_fit = cvfits[[2]], cvaucDAT = top2) %>% select(pred) %>% rename(pred2 = pred))  %>%
#   #   bind_cols(get_cv_predictions(cv_fit = cvfits[[3]], cvaucDAT = top2) %>% select(pred) %>% rename(pred3 = pred)) %>%
#   #   bind_cols(get_cv_predictions(cv_fit = cvfits[[4]], cvaucDAT = top2) %>% select(pred) %>% rename(pred4 = pred)) %>%
#   #   bind_cols(get_cv_predictions(cv_fit = cvfits[[5]], cvaucDAT = top2) %>% select(pred) %>% rename(pred5 = pred)) %>%
#   #   bind_cols(get_cv_predictions(cv_fit = cvfits[[6]], cvaucDAT = top2) %>% select(pred) %>% rename(pred6 = pred)) %>%
#   #   bind_cols(get_cv_predictions(cv_fit = cvfits[[7]], cvaucDAT = top2) %>% select(pred) %>% rename(pred7 = pred)) %>%
#   #   bind_cols(get_cv_predictions(cv_fit = cvfits[[8]], cvaucDAT = top2) %>% select(pred) %>% rename(pred8 = pred)) %>%
#   #   bind_cols(get_cv_predictions(cv_fit = cvfits[[9]], cvaucDAT = top2) %>% select(pred) %>% rename(pred9 = pred)) %>%
#   #   bind_cols(get_cv_predictions(cv_fit = cvfits[[10]], cvaucDAT = top2) %>% select(pred) %>% rename(pred10 = pred)) %>%
#   #   mutate(pred = rowMeans(select(., starts_with("pred")), na.rm = TRUE)) %>%
#   #   #left_join(top2 %>% select(Screen_fromRun, Learner, Screen, AUC, LearnerScreen), by = c("algo" = "LearnerScreen")) %>%
#   #   mutate(
#   #     learnerScreen = paste0(Learner, "_", Screen),
#   #     learnerScreen = ifelse(Learner %in% c("SL", "Discrete SL"), algo, learnerScreen),
#   #     AUCchar = format(round(AUC, 3), nsmall = 3),
#   #     learnerScreen = paste0(learnerScreen, " (", AUCchar, ")"),
#   #     learnerScreen = reorder(learnerScreen, -AUC)
#   #   )
# 
# 
#   # plot ROC curve
#   options(bitmapType = "cairo")
#   png(file = here("figs", paste0(Sys.getenv("TRIAL"), "/ROCcurve_", variableSet, ".png")),
#       width = 1000, height = 1000)
#   if(study_name %in% c("COVE", "MockCOVE")){
#     p1 <- plot_roc_curves(predict = pred, cvaucDAT = top2, weights = ph2_vacc_ptids %>% pull(wt.D57))
#     } else if(study_name == "HVTN705"){
#       p1 <- plot_roc_curves(predict = pred, cvaucDAT = top2, weights = ph2_vacc_ptids %>% pull(wt.D210))
#     } else if(study_name == "ENSEMBLE"){
#       p1 <- plot_roc_curves(predict = pred, cvaucDAT = top2, weights = ph2_vacc_ptids %>% pull(wt.D29))
#     } else if(study_name == "COVAIL" & TRIAL == "covail_tcell"){
#       p1 <- plot_roc_curves(predict = pred, cvaucDAT = top2, weights = ph2_vacc_ptids %>% pull(wt.D15.tcell))
#     } else if(study_name == "COVAIL" & TRIAL == "covail_xassays"){
#       p1 <- plot_roc_curves(predict = pred, cvaucDAT = top2, weights = ph2_vacc_ptids %>% pull(wt.D15.xassays))
#     }
#   print(p1)
#   dev.off()
# 
#   # plot pred prob plot
#   options(bitmapType = "cairo")
#   png(file = here("figs", paste0(Sys.getenv("TRIAL"), "/predProb_", variableSet, ".png")),
#       width = 1000, height = 1000)
#   p2 <- plot_predicted_probabilities(pred)
#   print(p2)
#   dev.off()
# }



# plot ROC curve and pred.Prob with SL, Discrete SL and top 2 best-performing individual Learners for all 12 variable sets
for(i in 1:6) {
  variableSet = top5_plus_baseline$varset[i]
  print(paste0("Processing ", variableSet))
  dat <- cvaucs_vacc %>% filter(varset==variableSet)
  
  top2 <- bind_rows(
    dat %>%
      arrange(-AUC) %>%
      filter(!Learner %in% c("SL", "Discrete SL")) %>%
      dplyr::slice(1:2),
    dat %>%
      filter(Learner == "SL"),
    dat %>%
      filter(Learner == "Discrete SL")
  ) %>%
    mutate(LearnerScreen = ifelse(Learner == "SL", "Super Learner",
                                  ifelse(Learner == "Discrete SL", Learner,
                                         paste0(Learner, "_", Screen_fromRun))))
  
  # Get cvsl fit and extract cv predictions
  if(study_name %in% c("COVE", "MockCOVE")){
    cvfit <- readRDS(file = here("output", Sys.getenv("TRIAL"), paste0("CVSLfits_vacc_EventIndPrimaryD57_", variableSet, ".rds")))
  } else if(study_name == "HVTN705"){
    cvfit<- readRDS(file = here("output", Sys.getenv("TRIAL"), paste0("CVSLfits_vacc_Delta.D210_", variableSet, ".rds")))
  } else if(study_name == "ENSEMBLE"){
    cvfit<- readRDS(file = here("output", Sys.getenv("TRIAL"), paste0("CVSLfits_vacc_EventIndPrimaryIncludeNotMolecConfirmedD29_", variableSet, ".rds")))
  } else if(study_name == "COVAIL"){
    
    if(COR == "D15to91covail_tcell"){
      cvfit<- readRDS(file = here(file_path, paste0("CVSLfits_vacc_COVIDIndD22toD91_", variableSet, ".rds")))
    } else if(COR == "D15to181covail_tcell"){
      cvfit<- readRDS(file = here(file_path, paste0("CVSLfits_vacc_COVIDIndD22toD181_", variableSet, ".rds")))
    } else if(COR == "D15to91covail_xassays"){
      cvfit<- readRDS(file = here(file_path, paste0("CVSLfits_vacc_COVIDIndD22toD91_", variableSet, ".rds")))
    } else if(COR == "D15to181covail_xassays"){
      cvfit<- readRDS(file = here(file_path, paste0("CVSLfits_vacc_COVIDIndD22toD181_", variableSet, ".rds")))
    }
    
  }
  
  pred <- get_cv_predictions(cv_fit = cvfit[[rseed]], cvaucDAT = top2, markerDAT = NULL)
  # #Take average of predictions from the 10 random seeds
  # pred <- get_cv_predictions(cv_fit = cvfits[[1]], cvaucDAT = top2) %>% rename(pred1 = pred) %>%
  #   bind_cols(get_cv_predictions(cv_fit = cvfits[[2]], cvaucDAT = top2) %>% select(pred) %>% rename(pred2 = pred))  %>%
  #   bind_cols(get_cv_predictions(cv_fit = cvfits[[3]], cvaucDAT = top2) %>% select(pred) %>% rename(pred3 = pred)) %>%
  #   bind_cols(get_cv_predictions(cv_fit = cvfits[[4]], cvaucDAT = top2) %>% select(pred) %>% rename(pred4 = pred)) %>%
  #   bind_cols(get_cv_predictions(cv_fit = cvfits[[5]], cvaucDAT = top2) %>% select(pred) %>% rename(pred5 = pred)) %>%
  #   bind_cols(get_cv_predictions(cv_fit = cvfits[[6]], cvaucDAT = top2) %>% select(pred) %>% rename(pred6 = pred)) %>%
  #   bind_cols(get_cv_predictions(cv_fit = cvfits[[7]], cvaucDAT = top2) %>% select(pred) %>% rename(pred7 = pred)) %>%
  #   bind_cols(get_cv_predictions(cv_fit = cvfits[[8]], cvaucDAT = top2) %>% select(pred) %>% rename(pred8 = pred)) %>%
  #   bind_cols(get_cv_predictions(cv_fit = cvfits[[9]], cvaucDAT = top2) %>% select(pred) %>% rename(pred9 = pred)) %>%
  #   bind_cols(get_cv_predictions(cv_fit = cvfits[[10]], cvaucDAT = top2) %>% select(pred) %>% rename(pred10 = pred)) %>%
  #   mutate(pred = rowMeans(select(., starts_with("pred")), na.rm = TRUE)) %>%
  #   #left_join(top2 %>% select(Screen_fromRun, Learner, Screen, AUC, LearnerScreen), by = c("algo" = "LearnerScreen")) %>%
  #   mutate(
  #     learnerScreen = paste0(Learner, "_", Screen),
  #     learnerScreen = ifelse(Learner %in% c("SL", "Discrete SL"), algo, learnerScreen),
  #     AUCchar = format(round(AUC, 3), nsmall = 3),
  #     learnerScreen = paste0(learnerScreen, " (", AUCchar, ")"),
  #     learnerScreen = reorder(learnerScreen, -AUC)
  #   )
  
  
  # plot ROC curve
  options(bitmapType = "cairo")
  png(file = here(file_path, "figs", paste0("ROCcurve_", variableSet, ".png")),
      width = 1000, height = 1000)
  if(study_name %in% c("COVE", "MockCOVE")){
    p1 <- plot_roc_curves(predict = pred, cvaucDAT = top2, weights = ph2_vacc_ptids %>% pull(wt.D57))
  } else if(study_name == "HVTN705"){
    p1 <- plot_roc_curves(predict = pred, cvaucDAT = top2, weights = ph2_vacc_ptids %>% pull(wt.D210))
  } else if(study_name == "ENSEMBLE"){
    p1 <- plot_roc_curves(predict = pred, cvaucDAT = top2, weights = ph2_vacc_ptids %>% pull(wt.D29))
  } else if(study_name == "COVAIL" & TRIAL == "covail_tcell"){
    p1 <- plot_roc_curves(predict = pred, cvaucDAT = top2, weights = ph2_vacc_ptids %>% pull(wt.D15.tcell))
  } else if(study_name == "COVAIL" & TRIAL == "covail_xassays"){
    p1 <- plot_roc_curves(predict = pred, cvaucDAT = top2, weights = ph2_vacc_ptids %>% pull(wt.D15.xassays))
  }
  print(p1)
  dev.off()
  
  # plot pred prob plot
  options(bitmapType = "cairo")
  png(file = here(file_path, "figs", paste0("predProb_", variableSet, ".png")),
      width = 1000, height = 1000)
  p2 <- plot_predicted_probabilities(pred)
  print(p2)
  dev.off()
}


# # For all variable sets, get predictors and coefficients for the DiscreteSL selected in each of the 5 outer folds from the 1st random seed!
# if(!study_name %in% c("HVTN705")){
#   
#   for(i in 1:length(only_varsets)) {
#     print(i)
#     variableSet = only_varsets[i]
#     print(variableSet)
# 
#     # Get cvsl fit and extract cv predictions
#     if(Sys.getenv("TRIAL") == "moderna_real"){
#       cvfits <- readRDS(file = here("output", Sys.getenv("TRIAL"), paste0("CVSLfits_vacc_EventIndPrimaryD57_", variableSet, ".rds")))
#     } else if(Sys.getenv("TRIAL") %in% c("janssen_pooled_partA", "janssen_la_partA")){
#       cvfits <- readRDS(file = here("output", Sys.getenv("TRIAL"), paste0("CVSLfits_vacc_EventIndPrimaryIncludeNotMolecConfirmedD29_", variableSet, ".rds")))
#     } else if(Sys.getenv("TRIAL") %in% c("covail_tcell")){
#       if(COR == "D15to91covail_tcell"){
#         cvfits <- readRDS(file = here("output", Sys.getenv("TRIAL"), paste0("CVSLfits_vacc_COVIDIndD22toD91_", variableSet, ".rds")))
#       } else if(COR == "D15to181covail_tcell"){
#         cvfits <- readRDS(file = here("output", Sys.getenv("TRIAL"), paste0("CVSLfits_vacc_COVIDIndD22toD181_", variableSet, ".rds")))
#       }
#     } else if(Sys.getenv("TRIAL") %in% c("covail_xassays")){
#       if(COR == "D15to91covail_xassays"){
#         cvfits <- readRDS(file = here("output", Sys.getenv("TRIAL"), paste0("CVSLfits_vacc_COVIDIndD22toD91_", variableSet, ".rds")))
#       } else if(COR == "D15to181covail_xassays"){
#         cvfits <- readRDS(file = here("output", Sys.getenv("TRIAL"), paste0("CVSLfits_vacc_COVIDIndD22toD181_", variableSet, ".rds")))
#       }
#     }
#     
#     # For selected random seed (rseed variable), get predictors and coefficients for the DiscreteSL selected in each of the 5 outer folds
#     for (j in seq_along(cvfits[[rseed]]$whichDiscreteSL)) {
#       print(paste0("j = ", j))
#       if(cvfits[[rseed]]$whichDiscreteSL[[j]] %in% c("SL.glm_screen_univariate_logistic_pval", 
#                                                  "SL.glm.interaction_screen_highcor_random",
#                                                  "SL.glm_screen_all",
#                                                  "SL.glm_screen_glmnet",
#                                                  "SL.glm_screen_highcor_random",
#                                                  "SL.glm.interaction_screen_glmnet",
#                                                  "SL.glm.interaction_screen_univariate_logistic_pval",
#                                                  "SL.gam_screen_glmnet",
#                                                  "SL.gam_screen_univariate_logistic_pval",
#                                                  "SL.gam_screen_highcor_random",
#                                                  "SL.bayesglm_screen_all",
#                                                  "SL.bayesglm_screen_glmnet",
#                                                  "SL.bayesglm_screen_univariate_logistic_pval",
#                                                  "SL.bayesglm_screen_highcor_random")) {
#         
#         model <- cvfits[[rseed]]$AllSL[[j]][["fitLibrary"]][[cvfits[[rseed]]$whichDiscreteSL[[j]]]]$object$coefficients %>%
#           as.data.frame() %>%
#           tibble::rownames_to_column(var = "Predictors") %>%
#           rename(`Coefficient` = ".") %>%
#           mutate(
#             `Odds Ratio` = exp(`Coefficient`),
#             Learner = cvfits[[rseed]]$whichDiscreteSL[[j]],
#             fold = j)
#       }
#       
#       # For SL.polymars
#       if(cvfits[[rseed]]$whichDiscreteSL[[j]] %in% c("SL.polymars_screen_glmnet", 
#                                                      "SL.polymars_screen_univariate_logistic_pval",
#                                                      "SL.polymars_screen_highcor_random")) {
#         
#         model <- flevr::extract_importance_polymars(polymars.obj$fit, feature_names = cvfits[[rseed]]$AllSL[[j]]$varNames[cvfits[[rseed]]$AllSL[[j]]$whichScreen[grepl(gsub("^[^_]*_", "", cvfits[[rseed]]$whichDiscreteSL[[j]]), rownames(cvfits[[rseed]]$AllSL[[j]]$whichScreen)),]]) %>%
#           filter(!is.na(importance)) %>%
#           as.data.frame() %>%
#           select(feature, importance) %>%
#           rename(`Coefficient` = `importance`,
#                  Predictors = feature) %>%
#           mutate(
#             Learner = cvfits[[rseed]]$whichDiscreteSL[[j]],
#             fold = j)
#       }
#       
#       # For SL.ksvm
#       if(cvfits[[rseed]]$whichDiscreteSL[[j]] %in% c("SL.polymars_screen_glmnet", 
#                                                      "SL.polymars_screen_univariate_logistic_pval",
#                                                      "SL.polymars_screen_highcor_random")) {
#         
#         model <- flevr::extract_importance_polymars(polymars.obj$fit, feature_names = cvfits[[rseed]]$AllSL[[j]]$varNames[cvfits[[rseed]]$AllSL[[j]]$whichScreen[grepl(gsub("^[^_]*_", "", cvfits[[rseed]]$whichDiscreteSL[[j]]), rownames(cvfits[[rseed]]$AllSL[[j]]$whichScreen)),]]) %>%
#           filter(!is.na(importance)) %>%
#           as.data.frame() %>%
#           select(feature, importance) %>%
#           rename(`Coefficient` = `importance`) %>%
#           mutate(
#             Learner = cvfits[[rseed]]$whichDiscreteSL[[j]],
#             fold = j)
#       }
#       
#       # For SL.nnet
#       if(cvfits[[rseed]]$whichDiscreteSL[[j]] %in% c("SL.nnet.2_screen_univariate_logistic_pval",
#                                                      "SL.nnet.2_screen_glmnet", 
#                                                      "SL.nnet.2_screen_highcor_random")) {
#         
#         model <- flevr::extract_importance_polymars(polymars.obj$fit, feature_names = cvfits[[rseed]]$AllSL[[j]]$varNames[cvfits[[rseed]]$AllSL[[j]]$whichScreen[grepl(gsub("^[^_]*_", "", cvfits[[rseed]]$whichDiscreteSL[[j]]), rownames(cvfits[[rseed]]$AllSL[[j]]$whichScreen)),]]) %>%
#           filter(!is.na(importance)) %>%
#           as.data.frame() %>%
#           select(feature, importance) %>%
#           rename(`Coefficient` = `importance`) %>%
#           mutate(
#             Learner = cvfits[[rseed]]$whichDiscreteSL[[j]],
#             fold = j)
#       }
#       
#       if (cvfits[[rseed]]$whichDiscreteSL[[j]] %in% c("SL.glmnet.0_screen_all", "SL.glmnet.0.33_screen_all", "SL.glmnet.0.67_screen_all", "SL.glmnet.1_screen_all")) {
#         
#         # model <- coef(cvfits[[rseed]]$AllSL[[j]][["fitLibrary"]][[cvfits[[rseed]]$whichDiscreteSL[[j]]]]$object, s = "lambda.min") %>%
#         #   as.matrix() %>%
#         #   as.data.frame() %>%
#         #   tibble::rownames_to_column(var = "Predictors") %>%
#         #   rename(`Coefficient` = "s1") %>%
#         #   mutate(`Odds Ratio` = exp(`Coefficient`),
#         #          Learner = cvfits[[rseed]]$whichDiscreteSL[[j]], 
#         #          fold = j)
#         
#         beta_matrix <- cvfits[[rseed]]$AllSL[[j]][["fitLibrary"]][[cvfits[[rseed]]$whichDiscreteSL[[j]]]]$object$glmnet.fit$beta
#         lambda_selected <- cvfits[[rseed]]$AllSL[[j]][["fitLibrary"]][[cvfits[[rseed]]$whichDiscreteSL[[j]]]]$object$lambda.min
#         lambda_index <- which(cvfits[[rseed]]$AllSL[[j]][["fitLibrary"]][[cvfits[[rseed]]$whichDiscreteSL[[j]]]]$object$lambda == lambda_selected)
#         beta_selected <- beta_matrix[, lambda_index]
#         
#         model <- beta_selected %>%
#           as.matrix() %>%
#           as.data.frame() %>%
#           tibble::rownames_to_column(var = "Predictors") %>%
#           rename(`Coefficient` = "V1") %>%
#           mutate(`Odds Ratio` = exp(`Coefficient`),
#                  Learner = cvfits[[rseed]]$whichDiscreteSL[[j]], 
#                  fold = j)
#       }
#       
#       if (cvfits[[rseed]]$whichDiscreteSL[[j]] %in% c("SL.xgboost.2.no_screen_all",
#                                                   "SL.xgboost.4.no_screen_all",
#                                                   "SL.xgboost.2.yes_screen_all",
#                                                   "SL.xgboost.4.yes_screen_all")) {
#         model <- xgboost::xgb.importance(model = cvfits[[rseed]]$AllSL[[j]][["fitLibrary"]][[cvfits[[rseed]]$whichDiscreteSL[[j]]]]$object) %>%
#           as.data.frame() %>%
#           mutate(Learner = cvfits[[rseed]]$whichDiscreteSL[[j]], 
#                  fold = j)
#       }
#       
#       if (cvfits[[rseed]]$whichDiscreteSL[[j]] %in% c("SL.ranger.yes_screen_all", "SL.ranger.no_screen_all")) {
#         model <- cvfits[[rseed]]$AllSL[[j]][["fitLibrary"]][[cvfits[[rseed]]$whichDiscreteSL[[j]]]]$object$variable.importance %>%
#           as.data.frame() %>%
#           rename(Importance = ".") %>%
#           tibble::rownames_to_column(var = "Predictors") %>%
#           mutate(Learner = cvfits[[rseed]]$whichDiscreteSL[[j]], 
#                  fold = j)
#       }
#       
#       if (cvfits[[rseed]]$whichDiscreteSL[[j]] == "SL.mean_screen_all"){
#         model <- data.frame(Learner="SL.mean_screen_all", 
#                             fold = j) 
#       }
#       
#       if (j == 1) {
#         all_models <- model
#       } else {
#         all_models <- bind_rows(all_models, model)
#       }
#     }
#     
#     if (i == 1) {
#       all_varsets_models <- all_models %>% mutate(varset = variableSet)
#     } else {
#       all_varsets_models <- bind_rows(all_varsets_models, 
#                                       all_models %>% mutate(varset = variableSet))  
#     }
#     
#     rm(variableSet, all_models, model)
#   }
# }




# For all variable sets, get predictors and coefficients for the DiscreteSL selected in each of the 5 outer folds from the 1st random seed!
if(!study_name %in% c("HVTN705")){
  
  for(i in 1:6) {
    print(i)
    variableSet = top5_plus_baseline$varset[i]
    print(variableSet)
    
    # Get cvsl fit and extract cv predictions
    if(Sys.getenv("TRIAL") == "moderna_real"){
      cvfits <- readRDS(file = here("output", Sys.getenv("TRIAL"), paste0("CVSLfits_vacc_EventIndPrimaryD57_", variableSet, ".rds")))
    } else if(Sys.getenv("TRIAL") %in% c("janssen_pooled_partA", "janssen_la_partA")){
      cvfits <- readRDS(file = here("output", Sys.getenv("TRIAL"), paste0("CVSLfits_vacc_EventIndPrimaryIncludeNotMolecConfirmedD29_", variableSet, ".rds")))
    } else if(Sys.getenv("TRIAL") %in% c("covail_tcell")){
      if(COR == "D15to91covail_tcell"){
        cvfits <- readRDS(file = here("output", Sys.getenv("TRIAL"), paste0("CVSLfits_vacc_COVIDIndD22toD91_", variableSet, ".rds")))
      } else if(COR == "D15to181covail_tcell"){
        cvfits <- readRDS(file = here("output", Sys.getenv("TRIAL"), paste0("CVSLfits_vacc_COVIDIndD22toD181_", variableSet, ".rds")))
      }
    } else if(Sys.getenv("TRIAL") %in% c("covail_xassays")){
      if(COR == "D15to91covail_xassays"){
        cvfits <- readRDS(file = here(file_path, paste0("CVSLfits_vacc_COVIDIndD22toD91_", variableSet, ".rds")))
      } else if(COR == "D15to181covail_xassays"){
        cvfits <- readRDS(file = here(file_path, paste0("CVSLfits_vacc_COVIDIndD22toD181_", variableSet, ".rds")))
      }
    }
    
    # For selected random seed (rseed variable), get predictors and coefficients for the DiscreteSL selected in each of the 5 outer folds
    for (j in seq_along(cvfits[[rseed]]$whichDiscreteSL)) {
      print(paste0("j = ", j))
      if(cvfits[[rseed]]$whichDiscreteSL[[j]] %in% c("SL.glm_screen_univariate_logistic_pval", 
                                                     "SL.glm.interaction_screen_highcor_random",
                                                     "SL.glm_screen_all",
                                                     "SL.glm_screen_glmnet",
                                                     "SL.glm_screen_highcor_random",
                                                     "SL.glm.interaction_screen_glmnet",
                                                     "SL.glm.interaction_screen_univariate_logistic_pval",
                                                     "SL.gam_screen_glmnet",
                                                     "SL.gam_screen_univariate_logistic_pval",
                                                     "SL.gam_screen_highcor_random",
                                                     "SL.bayesglm_screen_all",
                                                     "SL.bayesglm_screen_glmnet",
                                                     "SL.bayesglm_screen_univariate_logistic_pval",
                                                     "SL.bayesglm_screen_highcor_random")) {
        
        model <- cvfits[[rseed]]$AllSL[[j]][["fitLibrary"]][[cvfits[[rseed]]$whichDiscreteSL[[j]]]]$object$coefficients %>%
          as.data.frame() %>%
          tibble::rownames_to_column(var = "Predictors") %>%
          rename(`Coefficient` = ".") %>%
          mutate(
            `Odds Ratio` = exp(`Coefficient`),
            Learner = cvfits[[rseed]]$whichDiscreteSL[[j]],
            fold = j)
      }
      
      # For SL.polymars
      if(cvfits[[rseed]]$whichDiscreteSL[[j]] %in% c("SL.polymars_screen_glmnet", 
                                                     "SL.polymars_screen_univariate_logistic_pval",
                                                     "SL.polymars_screen_highcor_random")) {
        
        model <- flevr::extract_importance_polymars(polymars.obj$fit, feature_names = cvfits[[rseed]]$AllSL[[j]]$varNames[cvfits[[rseed]]$AllSL[[j]]$whichScreen[grepl(gsub("^[^_]*_", "", cvfits[[rseed]]$whichDiscreteSL[[j]]), rownames(cvfits[[rseed]]$AllSL[[j]]$whichScreen)),]]) %>%
          filter(!is.na(importance)) %>%
          as.data.frame() %>%
          select(feature, importance) %>%
          rename(`Coefficient` = `importance`,
                 Predictors = feature) %>%
          mutate(
            Learner = cvfits[[rseed]]$whichDiscreteSL[[j]],
            fold = j)
      }
      
      # For SL.ksvm
      if(cvfits[[rseed]]$whichDiscreteSL[[j]] %in% c("SL.polymars_screen_glmnet", 
                                                     "SL.polymars_screen_univariate_logistic_pval",
                                                     "SL.polymars_screen_highcor_random")) {
        
        model <- flevr::extract_importance_polymars(polymars.obj$fit, feature_names = cvfits[[rseed]]$AllSL[[j]]$varNames[cvfits[[rseed]]$AllSL[[j]]$whichScreen[grepl(gsub("^[^_]*_", "", cvfits[[rseed]]$whichDiscreteSL[[j]]), rownames(cvfits[[rseed]]$AllSL[[j]]$whichScreen)),]]) %>%
          filter(!is.na(importance)) %>%
          as.data.frame() %>%
          select(feature, importance) %>%
          rename(`Coefficient` = `importance`) %>%
          mutate(
            Learner = cvfits[[rseed]]$whichDiscreteSL[[j]],
            fold = j)
      }
      
      # For SL.nnet
      if(cvfits[[rseed]]$whichDiscreteSL[[j]] %in% c("SL.nnet.2_screen_univariate_logistic_pval",
                                                     "SL.nnet.2_screen_glmnet", 
                                                     "SL.nnet.2_screen_highcor_random")) {
        
        model <- flevr::extract_importance_polymars(polymars.obj$fit, feature_names = cvfits[[rseed]]$AllSL[[j]]$varNames[cvfits[[rseed]]$AllSL[[j]]$whichScreen[grepl(gsub("^[^_]*_", "", cvfits[[rseed]]$whichDiscreteSL[[j]]), rownames(cvfits[[rseed]]$AllSL[[j]]$whichScreen)),]]) %>%
          filter(!is.na(importance)) %>%
          as.data.frame() %>%
          select(feature, importance) %>%
          rename(`Coefficient` = `importance`) %>%
          mutate(
            Learner = cvfits[[rseed]]$whichDiscreteSL[[j]],
            fold = j)
      }
      
      if (cvfits[[rseed]]$whichDiscreteSL[[j]] %in% c("SL.glmnet.0_screen_all", "SL.glmnet.0.33_screen_all", "SL.glmnet.0.67_screen_all", "SL.glmnet.1_screen_all")) {
        
        # model <- coef(cvfits[[rseed]]$AllSL[[j]][["fitLibrary"]][[cvfits[[rseed]]$whichDiscreteSL[[j]]]]$object, s = "lambda.min") %>%
        #   as.matrix() %>%
        #   as.data.frame() %>%
        #   tibble::rownames_to_column(var = "Predictors") %>%
        #   rename(`Coefficient` = "s1") %>%
        #   mutate(`Odds Ratio` = exp(`Coefficient`),
        #          Learner = cvfits[[rseed]]$whichDiscreteSL[[j]], 
        #          fold = j)
        
        beta_matrix <- cvfits[[rseed]]$AllSL[[j]][["fitLibrary"]][[cvfits[[rseed]]$whichDiscreteSL[[j]]]]$object$glmnet.fit$beta
        lambda_selected <- cvfits[[rseed]]$AllSL[[j]][["fitLibrary"]][[cvfits[[rseed]]$whichDiscreteSL[[j]]]]$object$lambda.min
        lambda_index <- which(cvfits[[rseed]]$AllSL[[j]][["fitLibrary"]][[cvfits[[rseed]]$whichDiscreteSL[[j]]]]$object$lambda == lambda_selected)
        beta_selected <- beta_matrix[, lambda_index]
        
        model <- beta_selected %>%
          as.matrix() %>%
          as.data.frame() %>%
          tibble::rownames_to_column(var = "Predictors") %>%
          rename(`Coefficient` = "V1") %>%
          mutate(`Odds Ratio` = exp(`Coefficient`),
                 Learner = cvfits[[rseed]]$whichDiscreteSL[[j]], 
                 fold = j)
      }
      
      if (cvfits[[rseed]]$whichDiscreteSL[[j]] %in% c("SL.xgboost.2.no_screen_all",
                                                      "SL.xgboost.4.no_screen_all",
                                                      "SL.xgboost.2.yes_screen_all",
                                                      "SL.xgboost.4.yes_screen_all")) {
        model <- xgboost::xgb.importance(model = cvfits[[rseed]]$AllSL[[j]][["fitLibrary"]][[cvfits[[rseed]]$whichDiscreteSL[[j]]]]$object) %>%
          as.data.frame() %>%
          mutate(Learner = cvfits[[rseed]]$whichDiscreteSL[[j]], 
                 fold = j)
      }
      
      if (cvfits[[rseed]]$whichDiscreteSL[[j]] %in% c("SL.ranger.yes_screen_all", "SL.ranger.no_screen_all")) {
        model <- cvfits[[rseed]]$AllSL[[j]][["fitLibrary"]][[cvfits[[rseed]]$whichDiscreteSL[[j]]]]$object$variable.importance %>%
          as.data.frame() %>%
          rename(Importance = ".") %>%
          tibble::rownames_to_column(var = "Predictors") %>%
          mutate(Learner = cvfits[[rseed]]$whichDiscreteSL[[j]], 
                 fold = j)
      }
      
      if (cvfits[[rseed]]$whichDiscreteSL[[j]] == "SL.mean_screen_all"){
        model <- data.frame(Learner="SL.mean_screen_all", 
                            fold = j) 
      }
      
      if (j == 1) {
        all_models <- model
      } else {
        all_models <- bind_rows(all_models, model)
      }
    }
    
    if (i == 1) {
      all_varsets_models <- all_models %>% mutate(varset = variableSet)
    } else {
      all_varsets_models <- bind_rows(all_varsets_models, 
                                      all_models %>% mutate(varset = variableSet))  
    }
    
    rm(variableSet, all_models, model)
  }
}





if("Feature" %in% colnames(all_varsets_models)){
  all_varsets_models %>% 
    mutate(`Predictors/Features` = ifelse(is.na(Predictors), Feature, Predictors)) %>%
    select(-c(Predictors, Feature)) %>%
    select(varset, fold, Learner, `Predictors/Features`, everything()) %>%
    write.csv(here(file_path, "all_varsets_all_folds_discreteSLmodels.csv"))
} else {
  all_varsets_models %>% 
    mutate(`Predictors/Features` = Predictors) %>%
    select(varset, fold, Learner, `Predictors/Features`, everything()) %>%
    write.csv(here(file_path, "all_varsets_all_folds_discreteSLmodels.csv"))
}

# # Variable importance forest plots ---------------------------------------------
# # save off all variable importance estimates as a table
# vim_estimates %>%
#   filter(quantity == "VIM") %>%
#   #select(-group) %>%
#   write.csv(here("output", Sys.getenv("TRIAL"), "vim_estimates.csv"))
# vim_estimates %>%
#   filter(quantity == "Predictiveness") %>%
#   #select(-group) %>%
#   write.csv(here("output", Sys.getenv("TRIAL"), "vim_predictiveness_estimates.csv"))
# 
# num_digits <- 3
# plot_vim_init <- vim_estimates %>%
#   mutate(text_ci = paste0(round(est, num_digits), " [",
#                         round(ci_ll, num_digits), ", ",
#                         round(ci_ul, num_digits), "]")) 
# 
# group_ests <- plot_vim_init %>%
#   filter(group)
# individual_ests <- plot_vim_init %>%
#   filter(!group)
# 
# # Groups
# plot_group_vim <- group_ests %>%
#   mutate(plot_ord = as.numeric(gsub("_[^_]*", "", variable_set)),
#          plot_name = factor(plot_ord, levels = plot_ord, labels = variable_set))
# est_group_vims <- plot_group_vim %>% filter(quantity == "VIM")
# est_group_predictiveness <- plot_group_vim %>% filter(quantity == "Predictiveness")
# 
# group_vim_text_pos <- round(max(est_group_vims$ci_ul, na.rm = TRUE), 2) + 0.05
# group_vim_forest_plot <- est_group_vims %>%
#   filter(!grepl("base", variable_set)) %>%
#   mutate(plot_name = fct_reorder(plot_name, est, .desc = F)) %>%
#   ggplot(aes(x = est, y = plot_name)) +
#   geom_point(color = "blue") +
#   geom_errorbarh(aes(xmin = ci_ll, xmax = ci_ul), height = 0.3, color = "blue") +
#   geom_text(aes(x = rep(group_vim_text_pos, nrow(est_group_vims)-1), label = text_ci), hjust = "left") +
#   ggtitle("Estimated Importance Relative to Baseline Risk Factors") +
#   xlab("Estimated Difference in CV-AUC [95% CI]") +
#   ylab("Variable Set Name") +
#   xlim(c(0, group_vim_text_pos + 0.1)) +
#   geom_vline(xintercept = 0.5, lty = "dashed") +
#   theme_bw()
# 
# ggsave(
#   group_vim_forest_plot, file = here::here("figs", Sys.getenv("TRIAL"), "group_vim_forest_plot.png"),
#   width = 11.5, height = 10, units = "in", dpi = 300
# )
# 
# group_pred_text_pos <- round(max(est_group_predictiveness$ci_ul, na.rm = TRUE), 2) + 0.05
# group_pred_forest_plot <- est_group_predictiveness %>%
#   mutate(plot_name = fct_reorder(plot_name, est, .desc = F)) %>%
#   ggplot(aes(x = est, y = plot_name)) +
#   geom_point(color = "blue") +
#   geom_errorbarh(aes(xmin = ci_ll, xmax = ci_ul), height = 0.3, color = "blue") +
#   geom_text(aes(x = rep(group_pred_text_pos, nrow(est_group_predictiveness)), label = text_ci), hjust = "left") +
#   ggtitle("Estimated Predictiveness") +
#   xlab("CV-AUC") +
#   ylab("Variable Set Name") +
#   xlim(c(0, group_pred_text_pos + 0.2)) +
#   geom_vline(xintercept = 0.5, lty = "dashed") +
#   theme_bw()
# 
# ggsave(
#   group_pred_forest_plot, file = here::here("figs", Sys.getenv("TRIAL"), "group_pred_forest_plot.png"),
#   width = 11.5, height = 10, units = "in", dpi = 300
# )
# 
# # Individual variables
# plot_individual_vim <- individual_ests
# est_individual_vims <- plot_individual_vim %>% filter(quantity == "VIM")
# est_individual_predictiveness <- plot_individual_vim %>% filter(quantity == "Predictiveness")
# 
# individual_vim_text_pos <- round(max(est_individual_vims$ci_ul, na.rm = TRUE), 2) + 0.05
# individual_vim_forest_plot <- est_individual_vims %>%
#   mutate(variable_set = fct_reorder(variable_set, est, .desc = F)) %>%
#   ggplot(aes(x = est, y = variable_set)) +
#   geom_point(color = "blue") +
#   geom_errorbarh(aes(xmin = ci_ll, xmax = ci_ul), height = 0.3, color = "blue") +
#   geom_text(aes(x = rep(individual_vim_text_pos, nrow(est_individual_vims)), label = text_ci), hjust = "left") +
#   ggtitle("Estimated Importance Relative to Baseline Risk Factors") +
#   xlab("Estimated Difference in CV-AUC") +
#   ylab("Variable Name") +
#   xlim(c(0, individual_vim_text_pos + 0.1)) +
#   geom_vline(xintercept = 0.5, lty = "dashed") +
#   theme_bw()
# 
# ggsave(
#   individual_vim_forest_plot, file = here::here("figs", Sys.getenv("TRIAL"), "individual_vim_forest_plot.png"),
#   width = 11.5, height = 10, units = "in", dpi = 300
# )
# 
# individual_pred_text_pos <- round(max(est_individual_predictiveness$ci_ul, na.rm = TRUE), 2) + 0.05
# individual_pred_forest_plot <- est_individual_predictiveness %>%
#   mutate(variable_set = fct_reorder(variable_set, est, .desc = F)) %>%
#   ggplot(aes(x = est, y = variable_set)) +
#   geom_point(color = "blue") +
#   geom_errorbarh(aes(xmin = ci_ll, xmax = ci_ul), height = 0.3, color = "blue") +
#   geom_text(aes(x = rep(individual_pred_text_pos, nrow(est_individual_predictiveness)), label = text_ci), hjust = "left") +
#   ggtitle("Estimated Predictiveness") +
#   xlab("CV-AUC") +
#   ylab("Variable Name") +
#   xlim(c(0, individual_pred_text_pos + 0.2)) +
#   geom_vline(xintercept = 0.5, lty = "dashed") +
#   theme_bw()
# 
# ggsave(
#   individual_pred_forest_plot, file = here::here("figs", Sys.getenv("TRIAL"), "individual_pred_forest_plot.png"),
#   width = 11.5, height = 10, units = "in", dpi = 300
# )
# 
