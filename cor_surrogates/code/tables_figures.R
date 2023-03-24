# Sys.setenv(TRIAL = "hvtn705second")
# Sys.setenv(TRIAL = "moderna_real")
# Sys.setenv(TRIAL = "janssen_pooled_partA")
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
load(paste0("output/", Sys.getenv("TRIAL"), "/objects_for_running_SL.rda"))

# Get vimps in one dataframe: ----------------------------------------------------------
pooled_ests_lst <- list.files(here(paste0("output/", Sys.getenv("TRIAL"))), pattern = "pooled_ests_*") %>%
    tibble(file = .) %>%
    mutate(listdat = lapply(paste0("output/", Sys.getenv("TRIAL"), "/", file), readRDS))

all_estimates = as.data.frame(do.call("rbind", pooled_ests_lst$listdat))

# add on the variable set name
vim_estimates <- all_estimates %>%
    mutate(variable_set = rep(varset_names, each = 2), .before = "s")

# save the output
saveRDS(vim_estimates, file = paste0("output/", Sys.getenv("TRIAL"), "/vim_estimates.rds"))
# ----------------------------------------------------------------------------------------
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
ph2_vacc_ptids <- readRDS(file = here::here("output", Sys.getenv("TRIAL"), "ph2_vacc_ptids.rds"))

# Select the random seed from which to display results
if(study_name %in% c("COVE", "MockCOVE", "HVTN705")){
  rseed = 1
} else if(study_name %in% c("ENSEMBLE")){
  rseed = 2
}

# Create tables ------------------------------------------------------------------------
# Table of learner/screen combinations
caption <- "All learner-screen combinations (16 in total) used as input to the Superlearner."

tab <- cvaucs_vacc %>%
  filter(!Learner %in% c("SL", "Discrete SL")) %>%
  select(Learner, Screen) %>%
  mutate(Screen = fct_relevel(Screen, c("all", "glmnet", "univar_logistic_pval",
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

tab %>% write.csv(here("output", Sys.getenv("TRIAL"), "learner-screens.csv"))

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

tab %>% write.csv(here("output", Sys.getenv("TRIAL"), "varsets.csv"))

# Create figures ---------------------------------------------------------------
# Forest plots for vaccine model
# vaccine group
options(bitmapType = "cairo")
for(i in 1:(cvaucs_vacc %>% filter(!is.na(varsetNo)) %>% distinct(varset) %>% nrow())) {
  variableSet = unique(cvaucs_vacc$varset)[i]
  png(file = here("figs", Sys.getenv("TRIAL"), paste0("forest_vacc_cvaucs_", variableSet, ".png")), width=1000, height=1100)
  top_learner <- make_forest_plot(cvaucs_vacc %>% filter(varset==variableSet))
  grid.arrange(top_learner$top_learner_nms_plot, top_learner$top_learner_plot, ncol=2)
  dev.off()
}

# All Superlearners
learner.choice = "SL"
png(file = here("figs", Sys.getenv("TRIAL"), paste0("forest_vacc_cvaucs_all", learner.choice, "s.png")), width=1000, height=1100)
top_learner <- make_forest_plot_SL_allVarSets(cvaucs_vacc %>% filter(!is.na(varsetNo)), learner.choice)
grid.arrange(top_learner$top_learner_nms_plot, top_learner$top_learner_plot, ncol=2)
dev.off()

# All discrete.SLs
learner.choice = "Discrete SL"
png(file = here("figs", Sys.getenv("TRIAL"), paste0("forest_vacc_cvaucs_all", learner.choice, "s.png")), width=1000, height=1100)
top_learner <- make_forest_plot_SL_allVarSets(cvaucs_vacc %>% filter(!is.na(varsetNo)), learner.choice)
grid.arrange(top_learner$top_learner_nms_plot, top_learner$top_learner_plot, ncol=2)
dev.off()

#################################################################################################################################
#################################################################################################################################
# plot ROC curve and pred.Prob with SL, Discrete SL and top 2 best-performing individual Learners for all 12 variable sets
for(i in 1:(cvaucs_vacc %>% filter(!is.na(varsetNo)) %>% distinct(varset) %>% nrow())) {
  variableSet = unique(cvaucs_vacc$varset)[i]
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
  png(file = here("figs", paste0(Sys.getenv("TRIAL"), "/ROCcurve_", variableSet, ".png")),
      width = 1000, height = 1000)
  if(study_name %in% c("COVE", "MockCOVE")){
    p1 <- plot_roc_curves(predict = pred, cvaucDAT = top2, weights = ph2_vacc_ptids %>% pull(wt.D57))
    } else if(study_name == "HVTN705"){
      p1 <- plot_roc_curves(predict = pred, cvaucDAT = top2, weights = ph2_vacc_ptids %>% pull(wt.D210))
    } else if(study_name == "ENSEMBLE"){
      p1 <- plot_roc_curves(predict = pred, cvaucDAT = top2, weights = ph2_vacc_ptids %>% pull(wt.D29))
    }
  print(p1)
  dev.off()

  # plot pred prob plot
  options(bitmapType = "cairo")
  png(file = here("figs", paste0(Sys.getenv("TRIAL"), "/predProb_", variableSet, ".png")),
      width = 1000, height = 1000)
  p2 <- plot_predicted_probabilities(pred)
  print(p2)
  dev.off()
}



# Get SL & Discrete SL performance for all variable sets
cvaucs_vacc %>% filter(!is.na(varsetNo)) %>%
  arrange(-AUC) %>% filter(Learner == "SL") %>%
  select(varset, AUCstr) %>%
  write.csv(here("output", Sys.getenv("TRIAL"), "SLperformance_allvarsets.csv"))

cvaucs_vacc %>% filter(!is.na(varsetNo)) %>%
  arrange(-AUC) %>% filter(Learner == "Discrete SL") %>%
  select(varset, AUCstr) %>%
  write.csv(here("output", Sys.getenv("TRIAL"), "DiscreteSLperformance_allvarsets.csv"))

# # Predicted probability of COVID-19 vs antibody marker (x-axis)
# marker_cvaucs_vacc <- readin_SLobjects_fromFolder(data_folder, file_pattern = "CVSLaucs*", endpoint = "EventIndPrimaryD57", trt = "vaccine") %>%
#   filter(file %in% c(paste0("CVSLaucs_vacc_EventIndPrimaryD57_", varset_names, ".rds"))) %>%
#   filter(!file %in% cvaucs_vacc$file)
# 
# for(i in 1:length(individualMarkers)) {
#   varMarker = marker_cvaucs_vacc %>% filter(file == paste0("CVSLaucs_vacc_EventIndPrimaryD57_", individualMarkers[i], ".rds"))
#   top2 <- bind_rows(
#     varMarker %>%
#       arrange(-AUC) %>%
#       filter(!Learner %in% c("SL", "Discrete SL")) %>%
#       dplyr::slice(1:2),
#     varMarker %>%
#       filter(Learner == "SL"),
#     varMarker %>%
#       filter(Learner == "Discrete SL")
#   ) %>%
#     mutate(LearnerScreen = ifelse(Learner == "SL", "Super Learner",
#                                   ifelse(Learner == "Discrete SL", Learner,
#                                          paste0(Learner, "_", Screen_fromRun))))
# 
#   # Get cvsl fit and extract cv predictions
#   if(study_name %in% c("COVE", "MockCOVE")){
#     cvfits <- readRDS(file = here("output", paste0("CVSLfits_vacc_EventIndPrimaryD57_", individualMarkers[i], ".rds")))
#   }
# 
#   pred <- get_cv_predictions(cv_fit = cvfits[[1]], cvaucDAT = top2, markerDAT = dat.ph2 %>% select(individualMarkers[i]))
#   # plot
#   options(bitmapType = "cairo")
#   png(file = here("figs", paste0("marker_predProb_", individualMarkers[i], ".png")),
#       width = 1000, height = 1000)
#   xMarker = dat.ph2 %>% pull(individualMarkers[i])
#   vecNum = match(str_split(individualMarkers[i], paste0(0:9, collapse = "|"))[[1]][3], assays)
#   xlab = assay_labels_short[vecNum]
#   titlelab = paste0(assay_labels[vecNum], "\n", "Day ", gsub("...([0-9]+).*$", "\\1", individualMarkers[i]))
#   xdat =
#   print(pred %>%
#           ggplot(aes(x=dat.ph2 %>% select(individualMarkers[i]), y=pred)) +
#           facet_wrap(~algo, ncol = 2, scales = "free") +
#           geom_point(size=2, shape=23) +
#           xlim(floor(min(dat.ph2 %>% pull(individualMarkers[i]))), ceiling(max(dat.ph2 %>% pull(individualMarkers[i])))) +
#           labs(y = "Predicted probability of COVID-19 disease",
#                x = xlab,
#                title = titlelab) +
#           scale_x_log10("x",
#                         breaks = scales::trans_breaks("log10", function(y) 10^y),
#                         labels = scales::trans_format("log10", math_format(10^.y))))
#   dev.off()
# }


# For all variable sets, get predictors and coefficients for the DiscreteSL selected in each of the 5 outer folds from the 1st random seed!
if(!study_name %in% c("HVTN705")){
  
  for(i in 1:length(only_varsets)) {
    print(i)
    variableSet = only_varsets[i]

    # Get cvsl fit and extract cv predictions
    if(Sys.getenv("TRIAL") == "moderna_real"){
      cvfits <- readRDS(file = here("output", Sys.getenv("TRIAL"), paste0("CVSLfits_vacc_EventIndPrimaryD57_", variableSet, ".rds")))
    } else if(Sys.getenv("TRIAL") %in% c("janssen_pooled_partA", "janssen_la_partA")){
      cvfits <- readRDS(file = here("output", Sys.getenv("TRIAL"), paste0("CVSLfits_vacc_EventIndPrimaryIncludeNotMolecConfirmedD29_", variableSet, ".rds")))
    }
    
    # For selected random seed (rseed variable), get predictors and coefficients for the DiscreteSL selected in each of the 5 outer folds
    for (j in seq_along(cvfits[[rseed]]$whichDiscreteSL)) {
      print(j)
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
      
      if (cvfits[[rseed]]$whichDiscreteSL[[j]] %in% c("SL.glmnet.0_screen_all", "SL.glmnet.1_screen_all")) {
        
        model <- coef(cvfits[[rseed]]$AllSL[[j]][["fitLibrary"]][[cvfits[[rseed]]$whichDiscreteSL[[j]]]]$object, s = "lambda.min") %>%
          as.matrix() %>%
          as.data.frame() %>%
          tibble::rownames_to_column(var = "Predictors") %>%
          rename(`Coefficient` = "s1") %>%
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
    write.csv(here("output", Sys.getenv("TRIAL"), "all_varsets_all_folds_discreteSLmodels.csv"))
} else {
  all_varsets_models %>% 
    mutate(`Predictors/Features` = Predictors) %>%
    select(varset, fold, Learner, `Predictors/Features`, everything()) %>%
    write.csv(here("output", Sys.getenv("TRIAL"), "all_varsets_all_folds_discreteSLmodels.csv"))
}

# Variable importance forest plots ---------------------------------------------
# save off all variable importance estimates as a table
vim_estimates %>%
  filter(quantity == "VIM") %>%
  #select(-group) %>%
  write.csv(here("output", Sys.getenv("TRIAL"), "vim_estimates.csv"))
vim_estimates %>%
  filter(quantity == "Predictiveness") %>%
  #select(-group) %>%
  write.csv(here("output", Sys.getenv("TRIAL"), "vim_predictiveness_estimates.csv"))

num_digits <- 3
plot_vim_init <- vim_estimates %>%
  mutate(text_ci = paste0(round(est, num_digits), " [",
                        round(ci_ll, num_digits), ", ",
                        round(ci_ul, num_digits), "]")) 

group_ests <- plot_vim_init %>%
  filter(group)
individual_ests <- plot_vim_init %>%
  filter(!group)

# Groups
plot_group_vim <- group_ests %>%
  mutate(plot_ord = as.numeric(gsub("_[^_]*", "", variable_set)),
         plot_name = factor(plot_ord, levels = plot_ord, labels = variable_set))
est_group_vims <- plot_group_vim %>% filter(quantity == "VIM")
est_group_predictiveness <- plot_group_vim %>% filter(quantity == "Predictiveness")

group_vim_text_pos <- round(max(est_group_vims$ci_ul, na.rm = TRUE), 2) + 0.05
group_vim_forest_plot <- est_group_vims %>%
  filter(!grepl("base", variable_set)) %>%
  mutate(plot_name = fct_reorder(plot_name, est, .desc = F)) %>%
  ggplot(aes(x = est, y = plot_name)) +
  geom_point(color = "blue") +
  geom_errorbarh(aes(xmin = ci_ll, xmax = ci_ul), height = 0.3, color = "blue") +
  geom_text(aes(x = rep(group_vim_text_pos, nrow(est_group_vims)-1), label = text_ci), hjust = "left") +
  ggtitle("Estimated Importance Relative to Baseline Risk Factors") +
  xlab("Estimated Difference in CV-AUC [95% CI]") +
  ylab("Variable Set Name") +
  xlim(c(0, group_vim_text_pos + 0.1)) +
  geom_vline(xintercept = 0.5, lty = "dashed") +
  theme_bw()

ggsave(
  group_vim_forest_plot, file = here::here("figs", Sys.getenv("TRIAL"), "group_vim_forest_plot.png"),
  width = 11.5, height = 10, units = "in", dpi = 300
)

group_pred_text_pos <- round(max(est_group_predictiveness$ci_ul, na.rm = TRUE), 2) + 0.05
group_pred_forest_plot <- est_group_predictiveness %>%
  mutate(plot_name = fct_reorder(plot_name, est, .desc = F)) %>%
  ggplot(aes(x = est, y = plot_name)) +
  geom_point(color = "blue") +
  geom_errorbarh(aes(xmin = ci_ll, xmax = ci_ul), height = 0.3, color = "blue") +
  geom_text(aes(x = rep(group_pred_text_pos, nrow(est_group_predictiveness)), label = text_ci), hjust = "left") +
  ggtitle("Estimated Predictiveness") +
  xlab("CV-AUC") +
  ylab("Variable Set Name") +
  xlim(c(0, group_pred_text_pos + 0.2)) +
  geom_vline(xintercept = 0.5, lty = "dashed") +
  theme_bw()

ggsave(
  group_pred_forest_plot, file = here::here("figs", Sys.getenv("TRIAL"), "group_pred_forest_plot.png"),
  width = 11.5, height = 10, units = "in", dpi = 300
)

# Individual variables
plot_individual_vim <- individual_ests
est_individual_vims <- plot_individual_vim %>% filter(quantity == "VIM")
est_individual_predictiveness <- plot_individual_vim %>% filter(quantity == "Predictiveness")

individual_vim_text_pos <- round(max(est_individual_vims$ci_ul, na.rm = TRUE), 2) + 0.05
individual_vim_forest_plot <- est_individual_vims %>%
  mutate(variable_set = fct_reorder(variable_set, est, .desc = F)) %>%
  ggplot(aes(x = est, y = variable_set)) +
  geom_point(color = "blue") +
  geom_errorbarh(aes(xmin = ci_ll, xmax = ci_ul), height = 0.3, color = "blue") +
  geom_text(aes(x = rep(individual_vim_text_pos, nrow(est_individual_vims)), label = text_ci), hjust = "left") +
  ggtitle("Estimated Importance Relative to Baseline Risk Factors") +
  xlab("Estimated Difference in CV-AUC") +
  ylab("Variable Name") +
  xlim(c(0, individual_vim_text_pos + 0.1)) +
  geom_vline(xintercept = 0.5, lty = "dashed") +
  theme_bw()

ggsave(
  individual_vim_forest_plot, file = here::here("figs", Sys.getenv("TRIAL"), "individual_vim_forest_plot.png"),
  width = 11.5, height = 10, units = "in", dpi = 300
)

individual_pred_text_pos <- round(max(est_individual_predictiveness$ci_ul, na.rm = TRUE), 2) + 0.05
individual_pred_forest_plot <- est_individual_predictiveness %>%
  mutate(variable_set = fct_reorder(variable_set, est, .desc = F)) %>%
  ggplot(aes(x = est, y = variable_set)) +
  geom_point(color = "blue") +
  geom_errorbarh(aes(xmin = ci_ll, xmax = ci_ul), height = 0.3, color = "blue") +
  geom_text(aes(x = rep(individual_pred_text_pos, nrow(est_individual_predictiveness)), label = text_ci), hjust = "left") +
  ggtitle("Estimated Predictiveness") +
  xlab("CV-AUC") +
  ylab("Variable Name") +
  xlim(c(0, individual_pred_text_pos + 0.2)) +
  geom_vline(xintercept = 0.5, lty = "dashed") +
  theme_bw()

ggsave(
  individual_pred_forest_plot, file = here::here("figs", Sys.getenv("TRIAL"), "individual_pred_forest_plot.png"),
  width = 11.5, height = 10, units = "in", dpi = 300
)

