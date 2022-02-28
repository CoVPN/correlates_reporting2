# Sys.setenv(TRIAL = "moderna_real")
#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))

# There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))

source(here::here("..", "_common.R"))
#-----------------------------------------------

# common setup for CV super learners and variable importance
source(here::here("code", "cor_surrogates_setup.R"))

# drop "SL.xgboost.2.yes" and "SL.xgboost.4.yes" from SL_library as class-balancing learners in the variable
# importance computation doesn’t make sense – the regression we’re doing there (to account for the two-phase sampling)
# is based on a continuous outcome, not a binary outcome, so there shouldn’t be any imbalance.
for (i in 1:length(SL_library)) {
  if(SL_library[[i]][1] %in% c("SL.xgboost.2.yes", "SL.xgboost.4.yes")){
    if(!exists("vec"))
      vec = vector()
    vec = c(vec, i)
  }
}

SL_library = SL_library[-vec]
sl_lib = sl_lib[-vec]
rm(vec)

# get pooled VIMs, predictiveness for each comparison of interest --------------
# get the same seeds as used for the initial CV SL fits
set.seed(20210216)
seeds <- round(runif(10, 1000, 10000))
# set up
all_estimates <- NULL
full_y <- dat.ph1 %>%
  pull(!!endpoint)
X <- dat.ph1 %>%
  select(!!ptidvar, !!briskfactors) %>%
  left_join(dat.ph2 %>%
    select(!!ptidvar, all_of(markers)), by = ptidvar
  ) %>%
  select(-!!ptidvar)

# read in the fits for the baseline risk factors
baseline_fits <- readRDS(here("output", paste0("CVSLfits_vacc_", endpoint, "_", varset_names[1], ".rds")))
baseline_aucs <- readRDS(here("output", paste0("CVSLaucs_vacc_", endpoint, "_", varset_names[1], ".rds")))
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
for (i in seq_len(nrow(varset_matrix))) {
  # the column indices of interest
  if (i == 1) {
    this_s <- which(briskfactors %in% names(X))
  } else {
    this_s <- which(varset_matrix[i, ]) + length(briskfactors)
  }
  # get the correct CV.SL lists
  full_fits <- readRDS(here("output", paste0("CVSLfits_vacc_", endpoint, "_", varset_names[i], ".rds")))
  full_aucs <- readRDS(here("output", paste0("CVSLaucs_vacc_", endpoint, "_", varset_names[i], ".rds")))
  if (!use_ensemble_sl) {
    full_fits <- lapply(as.list(1:length(baseline_fits)), function(i) {
      make_discrete_sl_auc(cvsl_fit = full_fits[[i]], all_aucs = full_aucs[[i]])
    })
  }
  if (i == 1) {
    vim_lst <- lapply(list_of_indices, function(l) {
      get_cv_vim(seed = seeds[l], Y = full_y, X = X, full_fit = full_fits[[l]], reduced_fit = naive_fits[[l]],
                 index = this_s, type = "auc", scale = "identity", cross_fitting_folds = cf_folds[[l]],
                 sample_splitting_folds = sample_splitting_folds_baseline[[l]], V = vim_V,
                 C = C, Z = c("Y", paste0("X", which(briskfactors %in% names(X)))), sl_lib = sl_lib,
                 ipc_est_type = "ipw", ipc_weights = all_ipw_weights_treatment, baseline = TRUE,
                 use_ensemble = use_ensemble_sl)
    })
  } else {
    # get variable importance for each fold
    vim_lst <- lapply(list_of_indices, function(l) {
      get_cv_vim(seed = seeds[l], Y = full_y, X = X, full_fit = full_fits[[l]], reduced_fit = baseline_fits[[l]],
                 index = this_s, type = "auc", scale = "identity", cross_fitting_folds = cf_folds[[l]],
                 sample_splitting_folds = sample_splitting_folds[[l]], V = vim_V,
                 C = C, Z = c("Y", paste0("X", which(briskfactors %in% names(X)))), sl_lib = sl_lib,
                 ipc_est_type = "ipw", ipc_weights = all_ipw_weights_treatment, use_ensemble = use_ensemble_sl)
    })
  }
  # pool variable importance and predictiveness over the list
  pooled_ests <- pool_cv_vim(vim_lst = vim_lst, scale = "identity")
  all_estimates <- bind_rows(all_estimates, pooled_ests)
}
# add on the variable set name
final_estimates <- all_estimates %>%
  mutate(variable_set = rep(varset_names, each = 2), .before = "s")

# save the output
saveRDS(final_estimates, file = here::here("output", "vim_estimates.rds"))
