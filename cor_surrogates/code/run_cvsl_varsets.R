# Sys.setenv(TRIAL = "hvtn705")
#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))

# There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))

source(here::here("..", "_common.R"))
#-----------------------------------------------

# common setup for CV super learners and variable importance
source(here::here("code", "cor_surrogates_setup.R"))

# obtain the job id
args <- commandArgs(trailingOnly = TRUE)
job_id <- as.numeric(args[2])

# grab the current variable set based on the job id
this_var_set <- varset_matrix[job_id, ]
cat("\n Running", varset_names[job_id], "variable set \n")

# select the markers corresponding to the job id (and variable set)
markers_start <- length(briskfactors) + 1
if (job_id == 1){
  X_markers_varset <- X_covars2adjust_ph2[1:length(briskfactors)] %>%
    select_if(function(x) any(!is.na(x))) # Drop column if it has 0 variance, and returned all NAN's from scale function.
} else{
  X_markers_varset <- bind_cols(X_covars2adjust_ph2[1:length(briskfactors)],
                                X_covars2adjust_ph2[markers_start:length(X_covars2adjust_ph2)][this_var_set]) %>%
    select_if(function(x) any(!is.na(x))) # Drop column if it has 0 variance, and returned all NAN's from scale function.
}

# ---------------------------------------------------------------------------------
# run super learner, with leave-one-out cross-validation and all screens
# do 10 random starts, average over these
# use assay groups as screens
# ---------------------------------------------------------------------------------
# ensure reproducibility
set.seed(20210216)
seeds <- round(runif(10, 1000, 10000)) # average over 10 random starts

# disable parallelization in openBLAS and openMP
library(RhpcBLASctl)
blas_get_num_procs()
blas_set_num_threads(1L)
stopifnot(blas_get_num_procs() == 1L)
omp_set_num_threads(1L)

# run the Super Learners
fits <- parallel::mclapply(seeds, FUN = run_cv_sl_once,
                           Y = Y,
                           X_mat = X_markers_varset,
                           family = "binomial",
                           obsWeights = weights,
                           all_weights = all_ipw_weights_treatment,
                           ipc_est_type = "ipw",
                           sl_lib = sl_lib,
                           method = "method.CC_nloglik",
                           cvControl = list(V = V_outer, stratifyCV = TRUE),
                           innerCvControl = list(list(V = V_inner)),
                           Z = Z_treatmentDAT,
                           C = C,
                           z_lib = "SL.glm",
                           scale = "identity", # new argument
                           vimp = FALSE,
                           mc.cores = num_cores
)

# compile all CV-AUCs and CV-SL fits
cvaucs <- list()
cvfits <- list()
for (i in 1:length(seeds)) {
  cvaucs[[i]] = fits[[i]]$cvaucs$aucs
  cvfits[[i]] = fits[[i]]$cvfits
}

# save off the output
saveRDS(cvaucs, file = here("output", paste0("CVSLaucs_vacc_", endpoint, "_", varset_names[job_id], ".rds")))
saveRDS(cvfits, file = here("output", paste0("CVSLfits_vacc_", endpoint, "_", varset_names[job_id], ".rds")))
# only save these objects once
if (job_id == 1) {
  saveRDS(ph2_vacc_ptids, file = here("output", "ph2_vacc_ptids.rds"))
  save(run_prod, Y, dat.ph1, dat.ph2, weights, dat.mock, briskfactors, endpoint, maxVar,
       V_outer, file = here("output", "objects_for_running_SL.rda"))
}
