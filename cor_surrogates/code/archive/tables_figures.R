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
load(paste0("output/", "objects_for_running_SL.rda"))

# read in the results; note that createRDAfiles_fromSLobjects has to be run prior to this
if (study_name %in% c("COVE", "MockCOVE")) {
  cvaucs_vacc <- readRDS(file = here::here("output", "cvaucs_vacc_EventIndPrimaryD57.rds"))
  vim_estimates <- readRDS(file = here::here("output", "vim_estimates.rds"))  %>%
    mutate(group = ifelse(variable_set %in% c("2_bAbSpike_D57", "3_bAbRBD_D57", "13_bAbSpike_D29"), TRUE, group))
}
if (study_name == "HVTN705") {
  cvaucs_vacc <- readRDS(file = here::here("output", "cvaucs_vacc_EventIndPrimaryD210.rds"))
  vim_estimates <- readRDS(file = here::here("output", "vim_estimates.rds")) %>%
    mutate(group = ifelse(variable_set %in% c("4_M7_ADCP", "11_M7_IgG3multi", "12_M7_IgG3overall"), TRUE, group))
}
ph2_vacc_ptids <- readRDS(file = here::here("output", "ph2_vacc_ptids.rds"))

# Create tables ----------------------------------------------------------------
# Table of learner/screen combinations
caption <- "All learner-screen combinations (28 in total) used as input to the Superlearner."

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
    mutate(Learner = fct_relevel(Learner, c("SL.mean", #"SL.bayesglm", "SL.gam",
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

tab %>% write.csv(here("output", "learner-screens.csv"))

# Table of variable set definitions
if (study_name %in% c("COVE", "MockCOVE")) {
  caption <- "The 34 variable sets on which an estimated optimal surrogate was built."

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
  caption <- "The 15 variable sets on which an estimated optimal surrogate was built."

  tab <- data.frame(`Variable Set Name` = c("1_baselineRiskFactors",
                                            "2_M7_ELISA", #"3_bAbRBD_D57",
                                            "4_M7_ADCP", "5_M7_IgG3", "6_M7_IgG3gp140", "7_M7_IgG3gp120", "8_M7_IgG3V1V2", "9_M7_IgG3gp41", "10_M7_IgG3bScores",
                                            "11_M7_IgG3multi", "12_M7_IgG3overall", #"13_pnabID50_D29",
                                            "14_2+4", "15_2+5", #"16_bAb_pnabID80_D29", "17_bAb_combScores_D29",
                                            "18_4+5", "22_2+4+5"),
                    `Variables included in the set` = c("Baseline risk factors only (Reference model)",
                                                        "Baseline risk factors + M7 ELISA",
                                                        "Baseline risk factors + M7 ADCP",
                                                        "Baseline risk factors + M7 IgG3",
                                                        "Baseline risk factors + M7 IgG3 gp140",
                                                        "Baseline risk factors + M7 IgG3 gp120",
                                                        "Baseline risk factors + M7 IgG3 V1V2",
                                                        "Baseline risk factors + M7 IgG3 gp41",
                                                        "Baseline risk factors + M7 IgG3 Breadth Scores",
                                                        "Baseline risk factors + M7 IgG3 Multi-Epitope breadth",
                                                        "Baseline risk factors + M7 Overall score across assays",
                                                        "Baseline risk factors + M7 ELISA + M7 ADCP",
                                                        "Baseline risk factors + M7 ELISA + M7 IgG3",
                                                        "Baseline risk factors + M7 ADCP + M7 IgG3",
                                                        "Baseline risk factors + M7 ELISA + M7 ADCP + M7 IgG3"))

}

tab %>% write.csv(here("output", "varsets.csv"))

# Create figures ---------------------------------------------------------------
# Forest plots for vaccine model
# vaccine group
options(bitmapType = "cairo")
for(i in 1:length(unique(cvaucs_vacc$varset))) {
  variableSet = unique(cvaucs_vacc$varset)[i]
  png(file = here("figs", paste0("forest_vacc_cvaucs_", variableSet, ".png")), width=1000, height=1100)
  top_learner <- make_forest_plot(cvaucs_vacc %>% filter(varset==variableSet))
  grid.arrange(top_learner$top_learner_nms_plot, top_learner$top_learner_plot, ncol=2)
  dev.off()
}

# All Superlearners
learner.choice = "SL"
png(file = here("figs", paste0("forest_vacc_cvaucs_all", learner.choice, "s.png")), width=1000, height=1100)
top_learner <- make_forest_plot_SL_allVarSets(cvaucs_vacc, learner.choice)
grid.arrange(top_learner$top_learner_nms_plot, top_learner$top_learner_plot, ncol=2)
dev.off()

# All discrete.SLs
learner.choice = "Discrete SL"
png(file = here("figs", paste0("forest_vacc_cvaucs_all", learner.choice, "s.png")), width=1000, height=1100)
top_learner <- make_forest_plot_SL_allVarSets(cvaucs_vacc, learner.choice)
grid.arrange(top_learner$top_learner_nms_plot, top_learner$top_learner_plot, ncol=2)
dev.off()

#################################################################################################################################
#################################################################################################################################
# plot ROC curve and pred.Prob with SL, Discrete SL and top 2 best-performing individual Learners for all 12 variable sets
for(i in 1:length(unique(cvaucs_vacc$varset))) {
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
    cvfits <- readRDS(file = here("output", paste0("CVSLfits_vacc_EventIndPrimaryD57_", variableSet, ".rds")))
  }
  if(study_name == "HVTN705"){
    cvfits <- readRDS(file = here("output", paste0("CVSLfits_vacc_Delta.D210_", variableSet, ".rds")))
  }

  pred <- get_cv_predictions(cv_fit = cvfits[[1]], cvaucDAT = top2, markerDAT = NULL)
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
  png(file = here("figs", paste0("ROCcurve_", variableSet, ".png")),
      width = 1000, height = 1000)
  if(study_name %in% c("COVE", "MockCOVE")){
    p1 <- plot_roc_curves(predict = pred, cvaucDAT = top2, weights = ph2_vacc_ptids %>% pull(wt.D57))
  }
  if(study_name == "HVTN705"){
    p1 <- plot_roc_curves(predict = pred, cvaucDAT = top2, weights = ph2_vacc_ptids %>% pull(wt.D210))
  }
  print(p1)
  dev.off()

  # plot pred prob plot
  options(bitmapType = "cairo")
  png(file = here("figs", paste0("predProb_", variableSet, ".png")),
      width = 1000, height = 1000)
  p2 <- plot_predicted_probabilities(pred)
  print(p2)
  dev.off()
}


# Get SL & Discrete SL performance for all variable sets
cvaucs_vacc %>% arrange(-AUC) %>%
  filter(Learner == "SL") %>%
  select(varset, AUCstr) %>%
  write.csv(here("output", "SLperformance_allvarsets.csv"))

cvaucs_vacc %>% arrange(-AUC) %>%
  filter(Learner == "Discrete SL") %>%
  select(varset, AUCstr) %>%
  write.csv(here("output", "DiscreteSLperformance_allvarsets.csv"))

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
#           scale_x_log10(limits=c(10^0, 10^7), breaks=10^(0:7)))
#   dev.off()
# }
# 



# Variable importance forest plots ---------------------------------------------
# save off all variable importance estimates as a table
vim_estimates %>%
  filter(quantity == "VIM") %>%
  #select(-group) %>%
  write.csv(here("output", "vim_estimates.csv"))
vim_estimates %>%
  filter(quantity == "Predictiveness") %>%
  #select(-group) %>%
  write.csv(here("output", "vim_predictiveness_estimates.csv"))

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
  geom_text(aes(x = rep(group_vim_text_pos, nrow(est_group_vims) - 1), label = text_ci), hjust = "left") +
  ggtitle("Estimated Importance Relative to Baseline Risk Factors") +
  xlab("Estimated Difference in CV-AUC [95% CI]") +
  ylab("Variable Set Name") +
  xlim(c(0, group_vim_text_pos + 0.1)) +
  geom_vline(xintercept = 0.5, lty = "dashed") +
  theme_bw()

ggsave(
  group_vim_forest_plot, file = here::here("figs", "group_vim_forest_plot.png"),
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
  group_pred_forest_plot, file = here::here("figs", "group_pred_forest_plot.png"),
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
  individual_vim_forest_plot, file = here::here("figs", "individual_vim_forest_plot.png"),
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
  individual_pred_forest_plot, file = here::here("figs", "individual_pred_forest_plot.png"),
  width = 11.5, height = 10, units = "in", dpi = 300
)

