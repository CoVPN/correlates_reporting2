# Sys.setenv(TRIAL = "janssen_pooled_realbAb")
renv::activate(here::here(".."))
# There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))
source(here::here("..", "_common.R"))
#-----------------------------------------------

# load required libraries and functions
library(tidyverse)
library(here)
library(methods)
library(SuperLearner)
library(e1071)
library(glmnet)
library(kyotil)
library(argparse)
library(vimp)
library(nloptr)
library(RhpcBLASctl)
library(aucm)
library(mice)
library(conflicted)
library(gam)
library(xgboost)
library(ranger)
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("omp_set_num_threads", "RhpcBLASctl")
load(paste0("output/", Sys.getenv("TRIAL"), "/objects_for_running_SL.rda"))
load(paste0("output/", Sys.getenv("TRIAL"), "/plac_top2learners_SL_discreteSL.rda"))
source(here("code", "sl_screens.R")) # set up the screen/algorithm combinations
source(here("code", "utils.R")) # get CV-AUC for all algs

# ## solve cores issue
# #blas_get_num_procs()
# blas_set_num_threads(1)
# #print(blas_get_num_procs())
# stopifnot(blas_get_num_procs() == 1)
# 
# ## construct superlearner on placebo arm-----------------------
# set.seed(20210216)
# sl_riskscore_slfits <- SuperLearner(
#   Y = Y, X = X_markers_varset, family = familyVar,
#   SL.library = SL_library, method = methodVar,
#   cvControl = list(V = V_outer, stratifyCV = TRUE), 
#   verbose = FALSE,
#   obsWeights = weights
# )
# 
# save(sl_riskscore_slfits, file = here("output", "sl_riskscore_slfits.rda"))

load(here("output", "sl_riskscore_slfits.rda"))

# Get Superlearner weights
sl_weights <- sl_riskscore_slfits$coef %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "Learner") %>%
  rename(`Weights` = ".") %>%
  arrange(-Weights)

################################################################################
topmodel1 <- sl_riskscore_slfits[["fitLibrary"]][[sl_weights$Learner[1]]]$object$variable.importance %>%
  as.data.frame() %>%
  rename(Importance = ".") %>%
  tibble::rownames_to_column(var = "Predictors") %>%
  mutate(Learner = sl_weights$Learner[1])

topmodel1

rangerfit <- sl_riskscore_slfits[["fitLibrary"]][["SL.ranger.yes_screen_all"]]
dat = rangerfit$object$predictions %>% data.frame() %>%
  mutate(risk = `X1`) %>%
  bind_cols(X_markers_varset)

################################################################################
# GAM figures were shared at HVTN-705 presentation at FGM 2022
topmodel2 <- sl_riskscore_slfits[["fitLibrary"]][[sl_weights$Learner[2]]]$object$coefficients %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "Predictors") %>%
  rename(`Coefficient` = ".") %>%
  mutate(
    `Odds Ratio` = exp(`Coefficient`),
    Learner = sl_weights$Learner[2])

topmodel2

gamfit <- sl_riskscore_slfits[["fitLibrary"]][["SL.gam_screen_univariate_logistic_pval"]]
d = gamfit[["object"]][["fitted.values"]] %>% data.frame
colnames(d) <- "fitted.value.Y"
d <- d %>% bind_cols(X_markers_varset %>% select(Day210ADCCWITO_pk, Day210ADCCWITO_pAUC)) %>%
  mutate(Day210ADCCWITO_pk = exp(Day210ADCCWITO_pk), 
         Day210ADCCWITO_pAUC = exp(Day210ADCCWITO_pAUC))

png(filename="figs/gamplot.png", width=1200, height=600)
a = ggplot(d, aes(Day210ADCCWITO_pk, fitted.value.Y)) +
  geom_point(size=4, shape=1) +
  #stat_smooth(method = "gam", formula = y ~ s(x), size = 1, aes(weight=weights), color="red") +
  stat_smooth(method = "gam", formula = y ~ s(x), size = 1, se=FALSE, lwd=4, size=4) +
  theme_bw() +
  theme(text = element_text(size=20)) 


b = ggplot(d, aes(Day210ADCCWITO_pAUC, fitted.value.Y)) +
  geom_point(size=4, shape=1) +
  #stat_smooth(method = "gam", formula = y ~ s(x), size = 1, aes(weight=weights), color="red") +
  stat_smooth(method = "gam", formula = y ~ s(x), size = 1, se=FALSE, lwd=4) +
  theme_bw() +
  theme(text = element_text(size=20)) 

grid.arrange(a, b, ncol=2)
dev.off()

################################################################################
