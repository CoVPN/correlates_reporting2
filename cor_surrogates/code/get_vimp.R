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

# obtain the job id
#args <- commandArgs(trailingOnly = TRUE)
#job_id <- as.numeric(args[2])
job_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# common setup for CV super learners and variable importance
source(here::here("code", "cor_surrogates_setup.R"))

# drop "SL.xgboost.2.yes" and "SL.xgboost.4.yes" from SL_library as class-balancing learners in the variable
# importance computation doesn’t make sense – the regression we’re doing there (to account for the two-phase sampling)
# is based on a continuous outcome, not a binary outcome, so there shouldn’t be any imbalance.
vec = vector()
for (i in 1:length(SL_library)) {
  if(SL_library[[i]][1] %in% c("SL.xgboost.2.yes", "SL.xgboost.4.yes")){
    if(!exists("vec"))
      vec = c(vec, i)
  }
}

if(length(vec) != 0){
  SL_library = SL_library[-vec]
  sl_lib = sl_lib[-vec]
}

rm(vec)

# get pooled VIMs, predictiveness for each comparison of interest --------------
# get the same seeds as used for the initial CV SL fits
set.seed(20210216)
seeds <- round(runif(10, 1000, 10000))
# set up
all_estimates <- NULL
X <- phase_1_data_treatmentDAT %>%
  select(!!ptidvar, !!briskfactors) %>%
  left_join(dat.ph2 %>%
    select(!!ptidvar, all_of(markers)), by = ptidvar
  ) %>%
  select(-!!ptidvar)

# read in the fits for the baseline risk factors
baseline_fits <- readRDS(here("output", paste0(Sys.getenv("TRIAL"), "/CVSLfits_vacc_", endpoint, "_", varset_names[1], ".rds")))
baseline_aucs <- readRDS(here("output", paste0(Sys.getenv("TRIAL"), "/CVSLaucs_vacc_", endpoint, "_", varset_names[1], ".rds")))
if (!use_ensemble_sl) {
  baseline_fits <- lapply(as.list(1:length(baseline_fits)), function(i) {
    make_discrete_sl_auc(cvsl_fit = baseline_fits[[i]], all_aucs = baseline_aucs[[i]])
  })
}

# get the common CV folds
list_of_indices <- as.list(seq_len(length(baseline_fits)))
cf_folds <- lapply(list_of_indices, function(l) vimp::get_cv_sl_folds(baseline_fits[[l]]$folds))
vim_V <- length(unique(cf_folds[[1]])) / 2

# get common sample-splitting folds
sample_splitting_folds <- lapply(list_of_indices, function(l) {
  set.seed(seeds[l])
  these_ss_folds <- vimp::make_folds(unique(cf_folds[[l]]), V = 2)
})
# flip the sample splitting folds for the baseline comparison so that we're comparing the correct splits
sample_splitting_folds_baseline <- lapply(sample_splitting_folds, function(l) {
    (-1) * l + 3
})

# get the naive predictor (outcome mean within each fold)
naive_fits <- lapply(list_of_indices, function(l) {
  these_folds <- cf_folds[[l]]
  all_means <- lapply(seq_len(length(unique(these_folds))), function(i) {
    rep(mean(Y[these_folds == i]), sum(these_folds == i))
  })
  this_mean <- vector("numeric", length(Y))
  for (i in seq_len(length(unique(these_folds)))) {
    this_mean[these_folds == i] <- all_means[[i]]
  }
  this_mean
})

# get VIMs etc.
# obtain the job id
job_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# grab the current variable set based on the job id
this_var_set <- varset_matrix[job_id, ]
cat("\n Running", varset_names[job_id], "variable set \n")
interval_scale <- "logit" # CIs are computed on logit scale; helps with values outside boundary of [0, 1]
ipw_scale <- "identity" # IPW correction is applied on identity scale
final_point_estimate <- "average" # helps with point estimate stability
# the column indices of interest
  if (job_id == 1) {
    this_s <- which(briskfactors %in% names(X))
  } else {
    this_s <- which(varset_matrix[job_id, ]) + length(briskfactors)
  }
  # get the correct CV.SL lists
  full_fits <- readRDS(here("output", paste0(Sys.getenv("TRIAL"), "/CVSLfits_vacc_", endpoint, "_", varset_names[job_id], ".rds")))
  full_aucs <- readRDS(here("output", paste0(Sys.getenv("TRIAL"), "/CVSLaucs_vacc_", endpoint, "_", varset_names[job_id], ".rds")))
  if (!use_ensemble_sl) {
    full_fits <- lapply(as.list(1:length(baseline_fits)), function(i) {
      make_discrete_sl_auc(cvsl_fit = full_fits[[i]], all_aucs = full_aucs[[i]])
    })
  }
  if (job_id == 1) {
    vim_lst <- lapply(list_of_indices, function(l) {
      get_cv_vim(seed = seeds[l], Y = full_y, X = X, full_fit = full_fits[[l]], reduced_fit = naive_fits[[l]],
                 index = this_s, type = "auc", scale = interval_scale, cross_fitting_folds = cf_folds[[l]],
                 sample_splitting_folds = sample_splitting_folds_baseline[[l]], V = vim_V,
                 C = C, Z = c("Y", paste0("X", which(briskfactors %in% names(X)))), sl_lib = sl_lib,
                 ipc_scale = ipw_scale, ipc_est_type = "ipw", 
                 ipc_weights = all_ipw_weights_treatment, baseline = TRUE,
                 use_ensemble = use_ensemble_sl,
                 final_point_estimate = final_point_estimate) # this last argument helps with point estimates
    })
  } else {
    # get variable importance for each fold
    vim_lst <- lapply(list_of_indices, function(l) {
      get_cv_vim(seed = seeds[l], Y = full_y, X = X, full_fit = full_fits[[l]], reduced_fit = baseline_fits[[l]],
                 index = this_s, type = "auc", scale = interval_scale, cross_fitting_folds = cf_folds[[l]],
                 sample_splitting_folds = sample_splitting_folds[[l]], V = vim_V,
                 C = C, Z = c("Y", paste0("X", which(briskfactors %in% names(X)))), sl_lib = sl_lib,
                 ipc_scale = ipw_scale, ipc_est_type = "ipw", 
                 ipc_weights = all_ipw_weights_treatment, 
                 use_ensemble = use_ensemble_sl,
                 final_point_estimate = final_point_estimate)
    })
  }
  # pool variable importance and predictiveness over the list
  pooled_ests <- pool_cv_vim(vim_lst = vim_lst, scale = interval_scale)
  
  # save off the output
  saveRDS(pooled_ests, file = here("output", paste0(Sys.getenv("TRIAL"), "/pooled_ests_", endpoint, "_", varset_names[job_id], ".rds")))

  