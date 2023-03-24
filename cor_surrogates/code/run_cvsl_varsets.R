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

cvControlVar = list(V = V_outer, stratifyCV = TRUE)
cvControl_quote = quote(list(V = V_outer, stratifyCV = TRUE))

if(V_inner == length(Y) - 1){
  V_inner_quote <- paste0("length(Y) - 1 = ", length(Y) - 1)
}
innerCvControlVar = list(list(V = V_inner))
innerCvControl_quote = quote(list(list(V = V_inner)))

# CV.SL inputs
familyVar = "binomial"
methodVar = "method.CC_nloglik"
interval_scaleVar = "logit"
ipc_scaleVar <- "identity"
ipc_est_typeVar = "ipw"
cvsl_args <- data.frame(matrix(ncol = 2, nrow = 10)) %>%
  rename(Argument = X1,
         Value = X2) %>%
  mutate(Argument = as.character(c("Cases/Total Subjects in vaccine group (%)", "family",
                                   "method", "CI scale", "IPW correction scale", "V_outer", "cvControl (outer CV control)",
                                   "V_inner", "innerCvControl", "Weighting")),
         Value = as.character(c(paste0(nv, "/", length(Y), " (", round(nv*100/length(Y), 2), "%)"), familyVar,
                                methodVar, interval_scaleVar, ipc_scaleVar, V_outer, cvControl_quote,
                                V_inner, innerCvControl_quote, ipc_est_typeVar)))

if(V_inner == length(Y) - 1){
  cvsl_args <- cvsl_args %>% mutate(Value = ifelse(Argument == "V_inner", V_inner_quote, Value))
}

cvsl_args %>% add_row(Argument = "vimp package version",
                      Value = paste0(packageVersion("vimp"))) %>%    # Add vimp package version
  write.csv(paste0("output/", Sys.getenv("TRIAL"), "/cvsl_args.csv"))
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
#stopifnot(blas_get_num_procs() == 1L) # Commented this out as it does not work as expected any more!
omp_set_num_threads(1L)
stopifnot(omp_get_max_threads() == 1L)

# run the Super Learners
fits <- parallel::mclapply(seeds, FUN = run_cv_sl_once,
                           Y = Y,
                           X_mat = X_markers_varset,
                           family = familyVar,
                           obsWeights = weights,
                           all_weights = all_ipw_weights_treatment,
                           ipc_est_type = ipc_est_typeVar,
                           sl_lib = sl_lib,
                           method = methodVar,
                           cvControl = cvControlVar,
                           innerCvControl = innerCvControlVar,
                           Z = Z_treatmentDAT,
                           C = C,
                           z_lib = "SL.glm",
                           scale = interval_scaleVar, # scale on which intervals are computed (helps with intervals lying outside (0, 1))
                           ipc_scale = ipc_scaleVar, # scale on which IPW correction is applied
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
saveRDS(cvaucs, file = paste0("output/", Sys.getenv("TRIAL"), paste0("/CVSLaucs_vacc_", endpoint, "_", varset_names[job_id], ".rds")))
saveRDS(cvfits, file = here("output/", Sys.getenv("TRIAL"), paste0("/CVSLfits_vacc_", endpoint, "_", varset_names[job_id], ".rds")))
# only save these objects once
if (job_id == 1) {
  saveRDS(ph2_vacc_ptids, file = paste0("output/", Sys.getenv("TRIAL"), "/ph2_vacc_ptids.rds"))
  save(run_prod, Y, dat.ph1, dat.ph2, weights, dat.mock, briskfactors, endpoint, maxVar,
       V_outer, varset_names, individualMarkers, SL_library, file = paste0("output/", Sys.getenv("TRIAL"), "/objects_for_running_SL.rda"))
}
cat("\n Finished ", varset_names[job_id], "variable set \n") 
