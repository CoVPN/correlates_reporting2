# Sys.setenv(TRIAL = "moderna_mock")
#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))

# There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))

source(here::here("..", "_common.R"))
#-----------------------------------------------

# common setup for CV super learners and variable importance
source(here::here("code", "cor_surrogates_setup.R"))

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

# get the common CV folds
list_of_indices <- as.list(seq_len(length(baseline_fits)))
cf_folds <- lapply(list_of_indices, function(l) vimp::get_cv_sl_folds(baseline_fits[[l]]$folds))
vim_V <- length(unique(cf_folds[[1]])) / 2

# get common sample-splitting folds
sample_splitting_folds <- lapply(list_of_indices, function(l) {
  set.seed(seeds[l])
  these_ss_folds <- vimp::make_folds(unique(cf_folds[[l]]), V = 2)
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
  # get variable importance for each fold
  vim_lst <- lapply(list_of_indices, function(l) {
    get_cv_vim(seed = seeds[l], Y = full_y, X = X, full_fit = full_fits[[l]], reduced_fit = baseline_fits[[l]],
               index = this_s, type = "auc", scale = "identity", cross_fitting_folds = cf_folds[[l]],
               sample_splitting_folds = sample_splitting_folds[[l]], V = vim_V,
               C = C, Z = c("Y", paste0("X", which(briskfactors %in% names(X)))), sl_lib = sl_lib,
               ipc_est_type = "ipw", ipc_weights = all_ipw_weights_treatment)
  })
  # pool variable importance and predictiveness over the list
  pooled_ests <- pool_cv_vim(vim_lst = vim_lst, scale = "identity")
  # if baseline risk variables, then vimp isn't meaningful
  if (i == 1) {
    pooled_ests[1, c("est", "se", "ci_ll", "ci_ul", "pval")] <- NA_real_
  }
  all_estimates <- bind_rows(all_estimates, pooled_ests)
}
# add on the variable set name
final_estimates <- all_estimates %>%
  mutate(variable_set = rep(varset_names, each = 2), .before = "s")

# save the output
saveRDS(final_estimates, file = here::here("output", "vim_estimates.rds"))
