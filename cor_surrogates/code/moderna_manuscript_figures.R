Sys.setenv(TRIAL = "moderna_real")
#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))

# There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))

source(here::here("..", "_common.R"))
#-----------------------------------------------
## ----load-all-SLobjects, message=FALSE, error=FALSE, warning=FALSE----------------------------------------------------------------------------------------------------

library("cvAUC")
library("conflicted")
library("tidyr")
library("purrr")
library("dplyr")
library("cowplot")
library("ggplot2")
library("vimp")
library("kyotil")
library(gridExtra)
library(forcats)
library(cowplot)
library(here)
# filter() is used in both dplyr and stats, so need to set the preference to dplyr
conflict_prefer("filter", "dplyr")
conflict_prefer("summarise", "dplyr")
source(here("code", "utils.R"))
method <- "method.CC_nloglik" # since SuperLearner relies on this to be in GlobalEnv
ggplot2::theme_set(theme_cowplot())
load(paste0("output/", Sys.getenv("TRIAL"), "/objects_for_running_SL.rda"))

# # Get vimps in one dataframe: ----------------------------------------------------------
# pooled_ests_lst <- list.files(here(paste0("output/", Sys.getenv("TRIAL"))), pattern = "pooled_ests_*") %>%
#   tibble(file = .) %>%
#   mutate(listdat = lapply(paste0("output/", Sys.getenv("TRIAL"), "/", file), readRDS))
# 
# all_estimates = as.data.frame(do.call("rbind", pooled_ests_lst$listdat))
# 
# # add on the variable set name
# final_estimates <- all_estimates %>%
#   mutate(variable_set = rep(varset_names[-1], each = 2), .before = "s")
# 
# # save the output
# saveRDS(final_estimates, file = paste0("output/", Sys.getenv("TRIAL"), "/vim_estimates.rds"))
# ----------------------------------------------------------------------------------------
# read in the results; note that createRDAfiles_fromSLobjects has to be run prior to this

#if (study_name %in% c("COVE", "MockCOVE")) {
  cvaucs_vacc <- readRDS(file = here::here("output", Sys.getenv("TRIAL"), "cvaucs_vacc_EventIndPrimaryD57.rds"))
#   vim_estimates <- readRDS(file = here::here("output", Sys.getenv("TRIAL"), "vim_estimates.rds"))  %>%
#     mutate(group = ifelse(variable_set %in% c("2_bAbSpike_D57", "3_bAbRBD_D57", "13_bAbSpike_D29"), TRUE, group))
# }
# if (study_name == "HVTN705") {
#   cvaucs_vacc <- readRDS(file = here::here("output", Sys.getenv("TRIAL"), "cvaucs_vacc_EventIndPrimaryD210.rds"))
#   vim_estimates <- readRDS(file = here::here("output", Sys.getenv("TRIAL"), "vim_estimates.rds")) %>%
#     mutate(group = ifelse(variable_set %in% c("3_M7Tcells", "36_M7_2_3_4_5_12_13"), TRUE, group))
# }
ph2_vacc_ptids <- readRDS(file = here::here("output", Sys.getenv("TRIAL"), "ph2_vacc_ptids.rds"))

# # Create tables ------------------------------------------------------------------------
# # Table of learner/screen combinations
# caption <- "All learner-screen combinations (16 in total) used as input to the Superlearner."
# 
# tab <- cvaucs_vacc %>%
#   filter(!Learner %in% c("SL", "Discrete SL")) %>%
#   select(Learner, Screen) %>%
#   mutate(Screen = fct_relevel(Screen, c("all", "glmnet", "univar_logistic_pval",
#                                         "highcor_random")),
#          Learner = as.factor(Learner)) %>%
#   arrange(Learner, Screen) %>%
#   distinct(Learner, Screen) %>%
#   rename("Screen*" = Screen)
# 
# if (!grepl("Mock", study_name) & study_name == "COVE") {
#   tab <- tab %>%
#     mutate(Learner = fct_relevel(Learner, c("SL.mean", "SL.glmnet.0", "SL.glmnet.1", "SL.xgboost.2.no", "SL.xgboost.4.no",
#                                             "SL.xgboost.2.yes", "SL.xgboost.4.yes", "SL.ranger.yes", "SL.ranger.no", "SL.glm"))) %>%
#     arrange(Learner, `Screen*`)
# } else if (!grepl("Mock", study_name) & study_name == "HVTN705") {
#   tab <- tab %>%
#     mutate(Learner = fct_relevel(Learner, c("SL.mean", "SL.gam", #"SL.bayesglm", 
#                                             "SL.glm", #"SL.glm.interaction",
#                                             "SL.glmnet.0",
#                                             "SL.glmnet.1",
#                                             #"SL.ksvm.polydot", "SL.ksvm.rbfdot",
#                                             #"SL.polymars",
#                                             "SL.xgboost.2.no",
#                                             "SL.xgboost.4.no",
#                                             "SL.xgboost.2.yes",
#                                             "SL.xgboost.4.yes",
#                                             "SL.ranger.no",
#                                             "SL.ranger.yes"))) %>%
#     arrange(Learner, `Screen*`)
# } else {
#   tab <- tab %>%
#     mutate(Learner = fct_relevel(Learner, c("SL.mean", "SL.glm"))) %>%
#     arrange(Learner, `Screen*`)
# }
# 
# tab %>% write.csv(here("output", Sys.getenv("TRIAL"), "learner-screens.csv"))
# 
# # Table of variable set definitions
# if (study_name %in% c("COVE", "MockCOVE")) {
#   caption <- "The 34 variable sets on which an estimated optimal surrogate was built."
#   
#   tab <- data.frame(`Variable Set Name` = varset_names[1:34],
#                     `Variables included in the set` = c("Baseline risk factors only (Reference model)",
#                                                         "Baseline risk factors + Day 57 bAb anti-Spike markers",
#                                                         "Baseline risk factors + Day 57 bAb anti-RBD markers",
#                                                         "Baseline risk factors + Day 57 p-nAb ID50 markers",
#                                                         "Baseline risk factors + Day 57 p-nAb ID80 markers",
#                                                         "Baseline risk factors + Day 57 l-nAb MN50 markers",
#                                                         "Baseline risk factors + Day 57 bAb markers and p-nAb ID50 markers",
#                                                         "Baseline risk factors + Day 57 bAb markers and p-nAb ID80 markers",
#                                                         "Baseline risk factors + Day 57 bAb markers and l-nAb MN50 markers",
#                                                         "Baseline risk factors + Day 57 bAb markers and combination scores across the five markers [PCA1, PCA2, FSDAM1/FSDAM2 (the first two
# components of nonlinear PCA), and the maximum signal diversity score]",
#                                                         "Baseline risk factors + all individual Day 57 marker variables",
#                                                         "Baseline risk factors + all individual Day 57 marker variables and their combination scores (Full model of Day 57 markers)",
#                                                         
#                                                         "Baseline risk factors + Day 29 bAb anti-Spike markers",
#                                                         "Baseline risk factors + Day 29 bAb anti-RBD markers",
#                                                         "Baseline risk factors + Day 29 p-nAb ID50 markers",
#                                                         "Baseline risk factors + Day 29 p-nAb ID80 markers",
#                                                         "Baseline risk factors + Day 29 l-nAb MN50 markers",
#                                                         "Baseline risk factors + Day 29 bAb markers and p-nAb ID50 markers",
#                                                         "Baseline risk factors + Day 29 bAb markers and p-nAb ID80 markers",
#                                                         "Baseline risk factors + Day 29 bAb markers and l-nAb MN50 markers",
#                                                         "Baseline risk factors + Day 29 bAb markers and combination scores across the five markers [PCA1, PCA2, FSDAM1/FSDAM2 (the first two
# components of nonlinear PCA), and the maximum signal diversity score]",
#                                                         "Baseline risk factors + all individual Day 29 marker variables",
#                                                         "Baseline risk factors + all individual Day 29 marker variables and their combination scores (Full model of Day 29 markers)",
#                                                         
#                                                         "Baseline risk factors + Day 29 and Day 57 bAb anti-Spike markers",
#                                                         "Baseline risk factors + Day 29 and Day 57 bAb anti-RBD markers",
#                                                         "Baseline risk factors + Day 29 and Day 57 p-nAb ID50 markers",
#                                                         "Baseline risk factors + Day 29 and Day 57 p-nAb ID80 markers",
#                                                         "Baseline risk factors + Day 29 and Day 57 l-nAb MN50 markers",
#                                                         "Baseline risk factors + Day 29 and Day 57 bAb markers and p-nAb ID50 markers",
#                                                         "Baseline risk factors + Day 29 and Day 57 bAb markers and p-nAb ID80 markers",
#                                                         "Baseline risk factors + Day 29 and Day 57 bAb markers and l-nAb MN50 markers",
#                                                         "Baseline risk factors + Day 29 and Day 57 bAb markers and combination scores across the ten markers [PCA1, PCA2, FSDAM1/FSDAM2 (the first two
# components of nonlinear PCA), and the maximum signal diversity score]",
#                                                         "Baseline risk factors + all individual Day 29 and Day 57 marker variables",
#                                                         "Baseline risk factors + all individual Day 29 and Day 57 marker variables and their combination scores (Full model of Day 29 and Day 57 markers)"))
# }
# if (study_name == "HVTN705") {
#   total_varsets = 36
#   caption <- "The 36 variable sets on which an estimated optimal surrogate was built."
#   
#   tab <- data.frame(`Variable Set Name` = c("1_baselineFactors", "2_M7ELISA", "3_M7Tcells","4_M7ADCP","5_M7IgG3",
#                                             "6_M7IgG3_gp140","7_M7IgG3_gp120","8_M7IgG3_V1V2","9_M7IgG3_gp41",
#                                             "10_M7IgG3_breadthScore","11_M7IgG3_multi","12_M7IgGt_V1V2","13_M7ADCC",
#                                             "14_M7_combPrimaryMarkers", "15_M7_combAllMarkers","16_M7_overallScore", 
#                                             "17_M7_2_3", "18_M7_2_4", "19_M7_2_5_12", "20_M7_2_13", "21_M7_3_4", 
#                                             "22_M7_3_5_12", "23_M7_3_13", "24_M7_4_5_12", "25_M7_4_13", 
#                                             "26_M7_5_12_13","27_M7_2_3_4", "28_M7_2_3_5_12","29_M7_2_3_13","30_M7_3_4_5_12",
#                                             "31_M7_3_4_13","32_M7_4_5_12_13","33_M7_2_3_4_5_12","34_M7_2_3_4_13",
#                                             "35_M7_3_4_5_12_13","36_M7_2_3_4_5_12_13"),
#                     `Variables included in the set` = c(
#                       "Baseline Factors South Africa, Age, BMI, Behavioral risk score",
#                       "1.+M7 ELISA ELISA markers ELISA VT-C and ELISA VT-M controlling for baseline factors (Cont-base)",
#                       "1.+M7 T cells ELISpot PTE Env, CD4/CD8 T cells Cont-base",
#                       "1.+M7 ADCP ADCP markers ADCP C97ZA, ADCP Mosaic Cont-base",
#                       "1.+M7 IgG3 All BAMA IgG3 markers Cont-base",
#                       "1.+M7 IgG3 gp140 All BAMA IgG3 gp140 markers Cont-base",
#                       "1.+M7 IgG3 gp120 All BAMA IgG3 gp120 markers Cont-base",
#                       "1.+M7 IgG3 V1V2 All BAMA IgG3 V1V2 markers Cont-base",
#                       "1.+M7 IgG3 gp41 All BAMA IgG3 gp41 markers Cont-base",
#                       "1.+M7 IgG3 breadth scores All BAMA antigen-specific breadth scores Cont-base",
#                       "1.+M7 IgG3 Multi-epitope breadth BAMA IgG3 multi-epitope cross-reactivity score Cont-base",
#                       "1.+M7 IgGt V1V2 All BAMA IgG V1V2 markers Cont-base",
#                       "1.+M7 ADCC All ADCC markers Cont-base",
#                       "1.+M7 comb scores over primary markers All combination scores derived over the 6 primary markers [PCA1, PCA2, FSDAM1/FSDAM2 (the first two components of nonlinear PCA), and the maximum signal diversity
#                                                           score He and Fong (2019)]",
#                       "1.+M7 comb scores over all markers All combination scores derived over all markers [PCA1, PCA2, FSDAM1/FSDAM2, and max.sig.div.score].",
#                       "1.+M7 Overall scores across assays Two overall combined scores Cont-base",
#                       "1.+2.+3. Do ELISA + T cells together improve classification?",
#                       "1.+2.+4. Do ELISA + ADCP together improve classification?",
#                       "1.+2.+5.+12. Do ELISA + IgG3/IgGt together improve classification?",
#                       "1.+2.+13. Do ELISA + ADCC together improve classification?",
#                       "1.+3.+4. Do T cells + ADCP together improve classification?",
#                       "1.+3.+5.+12. Do T cells + IgG3/IgGt together improve classification?",
#                       "1.+3.+13. Do T cells + ADCC together improve classification?",
#                       "1.+4.+5.+12. Do ADCP + IgG3/IgGt together improve classification?",
#                       "1.+4.+13. Do ADCP + ADCC together improve classification?",
#                       "1.+5.+12.+13. Do IgG3/IgGt + ADCC together improve classification?",
#                       "1.+2.+3.+4. Do ELISA + T cells + ADCP together improve classification?",
#                       "1.+2+3.+5.+12. Do ELISA + T cells + IgG3/IgGt together improve classification?",
#                       "1.+2+3.+13. Do ELISA + T cells + ADCC together improve classification?",
#                       "1.+3.+4.+5.+12. Do T cells + ADCP + IgG3/IgGt together improve classification?",
#                       "1.+3.+4.+13. Do T cells + ADCP + ADCC together improve classification?",
#                       "1.+4.+5.+12.+13. Do ADCP + IgG3/IgGt + ADCC together improve classification?",
#                       "1.+2.+3.+4.+5.+12. Do ELISA + T cells + ADCP + IgG3/IgGt together improve classification?",
#                       "1.+2.+3.+4.+13. Do ELISA + T cells + ADCP + ADCC together improve classification?",
#                       "1.+3.+4.+5.+12.+13. Do T cells + ADCP + IgG3/IgGt + ADCC together improve classification?",
#                       "1.+2.+3.+4.+5.+12.+13. (All) Do all the immune marker sets together improve classification?"
#                     ))
#   
# }
# 
# tab %>% write.csv(here("output", Sys.getenv("TRIAL"), "varsets.csv"))
# 
# # Create figures ---------------------------------------------------------------
# # Forest plots for vaccine model
# # vaccine group
# options(bitmapType = "cairo")
# for(i in 1:(cvaucs_vacc %>% filter(!is.na(varsetNo)) %>% distinct(varset) %>% nrow())) {
#   variableSet = unique(cvaucs_vacc$varset)[i]
#   png(file = here("figs", Sys.getenv("TRIAL"), paste0("forest_vacc_cvaucs_", variableSet, ".png")), width=1000, height=1100)
#   top_learner <- make_forest_plot(cvaucs_vacc %>% filter(varset==variableSet))
#   grid.arrange(top_learner$top_learner_nms_plot, top_learner$top_learner_plot, ncol=2)
#   dev.off()
# }
# 
# # All Superlearners
# learner.choice = "SL"
# png(file = here("figs", Sys.getenv("TRIAL"), paste0("forest_vacc_cvaucs_all", learner.choice, "s.png")), width=1000, height=1100)
# top_learner <- make_forest_plot_SL_allVarSets(cvaucs_vacc %>% filter(!is.na(varsetNo)), learner.choice)
# grid.arrange(top_learner$top_learner_nms_plot, top_learner$top_learner_plot, ncol=2)
# dev.off()

# All discrete.SLs
cvaucs_vacc_discreteSL <- cvaucs_vacc %>% filter(Learner == "Discrete SL") %>%
  arrange(-AUC) %>%
  filter(varset %in% c("1_baselineRiskFactors", "2_bAbSpike_D57", 
                       "27_pnabID80_D29_D57", "6_lnabMN50_D57", "34_allMarkers_combScores_D29_D57"))

write.csv(cvaucs_vacc_discreteSL, here("output", Sys.getenv("TRIAL"), "for_manuscript_discreteSL_forestplot.csv"))

#################################################################################################################################
#################################################################################################################################
# plot ROC curve and pred.Prob with SL, Discrete SL and top 2 best-performing individual Learners for all 12 variable sets
for(i in 1:(cvaucs_vacc_discreteSL %>% distinct(varset) %>% nrow())) {
  variableSet = unique(cvaucs_vacc_discreteSL$varset)[i]
  dat <- cvaucs_vacc_discreteSL %>% filter(varset==variableSet)
  
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
  cvfits <- readRDS(file = here("output", Sys.getenv("TRIAL"), paste0("CVSLfits_vacc_EventIndPrimaryD57_", variableSet, ".rds")))

  pred <- get_cv_predictions(cv_fit = cvfits[[1]], cvaucDAT = top2, markerDAT = NULL) %>%
    mutate(varset = variableSet) 
  
  if(i == 1){
    pred_data <- pred
  } else {
    pred_data <- bind_rows(pred_data, pred)
  }

  # # plot ROC curve
  # options(bitmapType = "cairo")
  # png(file = here("figs", paste0(Sys.getenv("TRIAL"), "/ROCcurve_", variableSet, ".png")),
  #     width = 1000, height = 1000)
  # if(study_name %in% c("COVE", "MockCOVE")){
  #   p1 <- plot_roc_curves(predict = pred, cvaucDAT = top2, weights = ph2_vacc_ptids %>% pull(wt.D57))
  # }
  # if(study_name == "HVTN705"){
  #   p1 <- plot_roc_curves(predict = pred, cvaucDAT = top2, weights = ph2_vacc_ptids %>% pull(wt.D210))
  # }
  # print(p1)
  # dev.off()
}

pred_data %>% write.csv(here("output", Sys.getenv("TRIAL"), "for_manuscript_discreteSL_pred_prob.csv"))

#################################################


#################################################

options(scipen=999)

if(run_prod){
  all_models %>%
    left_join(sl_weights, by = "Learner") %>%
    mutate(
      Weight = format(round(Weights, 3), nsmall = 3),
      Coefficient = format(round(Coefficient, 3), nsmall = 3),
      `Odds Ratio` = format(round(`Odds Ratio`, 3), nsmall = 3),
      Importance = format(round(Importance, 3), nsmall = 3),
      Gain = format(round(Gain, 3), nsmall = 3),
      Cover = format(round(Cover, 3), nsmall = 3),
      Frequency = format(round(Frequency, 3), nsmall = 3),
    ) %>%
    mutate(
      Screen = paste0("screen_", sapply(strsplit(Learner, "_screen_"), `[`, 2)),
      Learner = sapply(strsplit(Learner, "_screen"), `[`, 1)
    ) %>%
    select(Learner, Screen, Weight, Predictors, Coefficient, `Odds Ratio`,
           Importance, Feature, Gain, Cover, Frequency) %>%
    write.csv(here("output", Sys.getenv("TRIAL"), "SL_all_models_with_predictors.csv"))
}else{
  all_models %>%
    left_join(sl_weights, by = "Learner") %>%
    mutate(
      Weight = format(round(Weights, 3), nsmall = 3),
      Coefficient = format(round(Coefficient, 3), nsmall = 3),
      `Odds Ratio` = format(round(`Odds Ratio`, 3), nsmall = 3),
    ) %>%
    mutate(
      Screen = paste0("screen_", sapply(strsplit(Learner, "_screen_"), `[`, 2)),
      Learner = sapply(strsplit(Learner, "_screen"), `[`, 1)
    ) %>%
    select(Learner, Screen, Weight, Predictors, Coefficient, `Odds Ratio`) %>%
    write.csv(here("output", Sys.getenv("TRIAL"), "SL_all_models_with_predictors.csv"))
}
