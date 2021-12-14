# Sys.setenv(TRIAL = "moderna_mock")
#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))

# There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))

source(here::here("..", "_common.R"))
#-----------------------------------------------

# common setup for CV super learners and variable importance
source(here::here("cor_surrogates_setup.R"))

# get pooled VIMs, predictiveness for each comparison of interest --------------
# get the same seeds as used for the initial CV SL fits
set.seed(20210216)
seeds <- round(runif(10, 1000, 10000))
# set up
all_estimates <- NULL
X <- bind_cols(dat.ph2 %>%
  select(!!briskfactors), markers)
# get VIMs etc.
for (i in seq_len(nrow(varset_matrix) - 1)) {
  # the column indices of interest
  if (i == 1) {
    this_s <- which(briskfactors %in% names(X))
  } else {
    this_s <- which(varset_matrix[i, ]) + length(briskfactors)
  }
  # get the correct CV.SL lists
  # full_fits <-
  # reduced_fits <-
  list_of_indices <- as.list(seq_len(length(full_fits)))
  # get the common CV folds
  cf_folds <- lapply(list_of_indices, function(l) vimp::get_cv_sl_folds(full_fits[[l]]$folds))
  # get variable importance for each fold
  vim_lst <- lapply(list_of_indices, function(l) {
    get_cv_vim(seed = seeds[l], full_fit = full_fits[[l]], reduced_fit = reduced_fits[[l]], index = this_s, type = "auc", scale = "identity", )
  })
  # pool variable importance and predictiveness over the list
  pooled_ests <- pool_cv_vim(vim_lst = vim_lst, scale = "identity")
  all_estimates <- bind_rows(all_estimates, pooled_ests)
}
