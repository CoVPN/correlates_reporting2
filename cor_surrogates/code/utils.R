#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))

# There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))

source(here::here("..", "_common.R"))
#-----------------------------------------------

# Cross-validated predictiveness -----------------------------------------------

# get the CV-AUC for a single learner's predicted values
# @param preds the fitted values
# @param Y the outcome
# @param full_y the observed outcome (from the entire dataset, for cross-fitted estimates)
# @param scale what scale should the IPCW correction be applied on?
#              Can help with numbers outside of (0, 1)
#             ("identity" denotes the identity scale;
#              "logit" means that AUC is transformed to the logit scale,
#              correction is applied, and then back-transformed)
# @param weights the inverse probability of censoring weights
# @param C the indicator of being observed in phase 2 (1) or not (0)
# @param Z a matrix of predictors observed on all participants in phase 1
#          (can include the outcome)
# @param ... other arguments to measure_auc, but should include at least:
#            a library of learners (using arg "SL.library")
#            and may include control parameters for the super learner
#            (e.g., cvControl = list(V = 5)
#             for 5-fold cross-validated super learner)
one_auc <- function(preds, Y, full_y = NULL, scale = "identity",
                    weights = rep(1, length(Y)), C = rep(1, length(Y)),
                    Z = NULL, ...) {
  auc_lst <- vimp::measure_auc(
    fitted_values = preds, y = Y, full_y = full_y, C = C, Z = Z,
    ipc_weights = weights, scale = scale,
    ipc_fit_type = "SL", ...
  )
  list(auc = auc_lst$point_est, eif = auc_lst$eif)
}

# get the cross-fitted CV-AUC for a single learner's predicted values
# @param preds the fitted values
# @param Y the outcome
# @param folds the different cv folds that the learner was evaluated on
# @param scale what scale should the CI be applied on? Can help with numbers outside of (0, 1)
# @param ipc_scale what scale should the IPCW correction be applied on?
#              Can help with numbers outside of (0, 1)
#             ("identity" denotes the identity scale;
#              "logit" means that AUC is transformed to the logit scale,
#              correction is applied, and then back-transformed)
# @param weights the inverse probability of censoring weights
# @param C the indicator of being observed in phase 2 (1) or not (0)
# @param Z a matrix of predictors observed on all participants in phase 1
#          (can include the outcome)
# @param ... other arguments to measure_auc, but should include at least:
#            a library of learners (using arg "SL.library")
#            and may include control parameters for the super learner
#            (e.g., cvControl = list(V = 5)
#             for 5-fold cross-validated super learner)
cv_auc <- function(preds, Y, folds, scale = "logit",
                   ipc_scale = "identity",
                   weights = rep(1, length(Y)), C = rep(1, length(Y)),
                   Z = NULL, ...) {
  V <- length(folds)
  folds_numeric <- get_cv_sl_folds(folds)
  if (is.null(Z)) {
    folds_z <- folds_numeric
  } else {
    folds_z <- c(folds_numeric, sample(seq_len(V), nrow(Z) - length(folds_numeric),
      replace = TRUE
    ))
  }
  ests_eifs <- lapply(as.list(seq_len(V)), function(v) {
    one_auc(
      preds = preds[folds_numeric == v], Y[folds_numeric == v],
      full_y = Y, scale = ipc_scale,
      weights = weights[folds_z == v], C = C[folds_z == v],
      Z = Z[folds_z == v, , drop = FALSE], ...
    )
  })
  est <- mean(unlist(lapply(ests_eifs, function(l) l$auc)))
  var <- mean(unlist(lapply(ests_eifs, function(l) mean(l$eif ^ 2))))
  se <- sqrt(var / length(Y))
  ci <- vimp::vimp_ci(est, se, scale = scale, level = 0.95)
  return(list(auc = est, se = se, ci = ci))
}

# get the folds from a CV.SL object, make them a vector
# @param cv_sl_folds the CV.SL folds (a named list of row numbers)
# @return a vector with the correct folds
get_cv_sl_folds <- function(cv_sl_folds) {
  folds_with_row_nums <- sapply(1:length(cv_sl_folds),
    function(x) {
      list(
        row_nums = cv_sl_folds[[x]],
        fold = rep(x, length(cv_sl_folds[[x]]))
      )
    },
    simplify = FALSE
  )
  folds_df <- data.table::rbindlist(folds_with_row_nums)
  folds_df$fold[order(folds_df$row_nums)]
}

# get the CV-AUC for an individual learner in the SL
# @param sl_fit the fitted SL
# @param col the column of interest (corresponds to a fitted algorithm)
# @param scale scale for intervals estimation
# @param ipc_scale scale for EIF estimation
# @param weights IPC weights
# @param C the censoring indicator
# @param Z data observed on all participants
# @param ... other arguments to pass to Super Learner (for EIF estimation)
get_individual_auc <- function(sl_fit, col, scale = "logit",
                               ipc_scale = "identity",
                               weights = rep(1, length(sl_fit$Y)),
                               C = rep(1, length(sl_fit$Y)), Z = NULL, ...) {
  if (any(is.na(sl_fit$library.predict[, col]))) {
    return(NULL)
  }
  alg_auc <- cv_auc(
    preds = sl_fit$library.predict[, col], Y = sl_fit$Y,
    scale = scale, ipc_scale = ipc_scale,
    folds = sl_fit$folds, weights = weights,
    C = C, Z = Z, ...
  )
  # get the regexp object
  alg_screen_string <- strsplit(colnames(sl_fit$library.predict)[col], "_",
    fixed = TRUE
  )[[1]]
  alg <- tail(alg_screen_string[grepl(".", alg_screen_string,
    fixed = TRUE
  )], n = 1)
  screen <- paste0(alg_screen_string[!grepl(alg, alg_screen_string,
    fixed = TRUE
  )],
  collapse = "_"
  )
  data.frame(
    Learner = alg, Screen = screen, AUC = alg_auc$auc,
    se = alg_auc$se,
    ci_ll = alg_auc$ci[1], ci_ul = alg_auc$ci[2]
  )
}
# get the CV-AUC for all learners fit with SL
# @param sl_fit the super learner fit object
# @param scale scale for interval estimation
# @param ipc_scale what scale should the IPCW correction be applied on?
#              Can help with numbers outside of (0, 1)
#             ("identity" denotes the identity scale;
#              "logit" means that AUC is transformed to the logit scale,
#              correction is applied, and then back-transformed)
# @param weights the inverse probability of censoring weights
# @param C the indicator of being observed in phase 2 (1) or not (0)
# @param Z a matrix of predictors observed on all participants in phase 1
#          (can include the outcome)
# @param ... other arguments to measure_auc, but should include at least:
#            a library of learners (using arg "SL.library")
#            and may include control parameters for the super learner
#            (e.g., cvControl = list(V = 5)
#             for 5-fold cross-validated super learner)
get_all_aucs <- function(sl_fit, scale = "logit", ipc_scale = "identity",
                         weights = rep(1, length(sl_fit$Y)),
                         C = rep(1, length(sl_fit$Y)),
                         Z = NULL, ...) {
  # get the CV-AUC of the SuperLearner predictions
  sl_auc <- cv_auc(
    preds = sl_fit$SL.predict, Y = sl_fit$Y,
    folds = sl_fit$folds, scale = scale,
    ipc_scale = ipc_scale, weights = weights, C = C, Z = Z, ...
  )
  out <- data.frame(
    Learner = "SL", Screen = "All", AUC = sl_auc$auc,
    se = sl_auc$se, ci_ll = sl_auc$ci[1], ci_ul = sl_auc$ci[2]
  )

  # Get the CV-auc of the Discrete SuperLearner predictions
  discrete_sl_auc <- cv_auc(
    preds = sl_fit$discreteSL.predict, Y = sl_fit$Y,
    folds = sl_fit$folds, scale = scale, ipc_scale = ipc_scale,
    weights = weights, C = C,
    Z = Z, ...
  )
  out <- rbind(out, data.frame(
    Learner = "Discrete SL", Screen = "All",
    AUC = discrete_sl_auc$auc,
    se = discrete_sl_auc$se,
    ci_ll = discrete_sl_auc$ci[1],
    ci_ul = discrete_sl_auc$ci[2]
  ))

  # Get the cvauc of the individual learners in the library
  other_aucs <- plyr::ldply(
    1:ncol(sl_fit$library.predict),
    function(x) {
      get_individual_auc(
        sl_fit = sl_fit,
        col = x,
        scale = scale,
        ipc_scale = ipc_scale,
        weights = weights,
        C = C, Z = Z, ...
      )
    }
  )
  rbind(out, other_aucs)
}

# determine the discrete Super Learner using CV-AUC
# @param cvsl_fit the CV.SL object
# @param all_aucs the estimated CV-AUCs for that CV.SL object
make_discrete_sl_auc <- function(cvsl_fit, all_aucs) {
    individual_aucs <- all_aucs %>%
      filter(!(Learner %in% c("SL", "Discrete SL")))
    which_discrete_sl <- which.max(individual_aucs$AUC)
    cvsl_fit$discreteSL.predict <- cvsl_fit$library.predict[, which_discrete_sl]
    return(cvsl_fit)
}

# Run the CV Super Learner -----------------------------------------------------
# run CV.SuperLearner for one given random seed
# @param seed the random number seed
# @param Y the outcome
# @param X_mat the covariates
# @param family the family, for super learner (e.g., "binomial")
# @param obsWeights the inverse probability of censoring weights
#                   (for super learner)
# @param sl_lib the super learner library (e.g., "SL.ranger")
# @param method the method for determining the optimal combination
#               of base learners
# @param cvControl a list of control parameters to pass to the
#                  outer super learner
# @param innerCvControl a list of control parameters to pass to the
#                       inner super learners
# @param all_weights the IPC weights for variable importance (full data)
# @param Z the entire (phase 1) dataset (required only if analysis
#                  involves data from phase 2)
# @param C the outcome from the entire (phase 1) dataset
# @param z_lib the learner/s to be used for weighted auc calculations
# @param ipc_scale the scale that the IPC correction should be computed on
# @param scale scale for interval estimation
# @param vimp determines whether or not we save the entire SL fit object
#        (for variable importance, we don't need it so can save some memory
#         by excluding large objects)
run_cv_sl_once <- function(seed = 1, Y = NULL, X_mat = NULL,
                           family = "binomial",
                           obsWeights = rep(1, length(Y)),
                           ipc_est_type = "ipw",
                           all_weights = rep(1, nrow(Z)),
                           sl_lib = "SL.ranger",
                           method = "method.CC_nloglik",
                           cvControl = list(V = 5),
                           innerCvControl = list(V = 5),
                           Z = NULL,
                           C = rep(1, length(Y)),
                           z_lib = "SL.glm",
                           scale = "logit",
                           ipc_scale = "identity",
                           vimp = FALSE) {

  set.seed(seed)
  fit <- SuperLearner::CV.SuperLearner(
    Y = Y, X = X_mat, family = family,
    obsWeights = obsWeights,
    SL.library = sl_lib,
    control = list(saveFitLibrary = TRUE),
    method = method, cvControl = cvControl,
    innerCvControl = innerCvControl
    #, verbose = FALSE
  )

  aucs <- get_all_aucs(sl_fit = fit, scale = scale, ipc_scale = ipc_scale, weights = all_weights,
                       C = C, Z = Z, SL.library = z_lib, ipc_est_type = ipc_est_type, family = gaussian())

  ret_lst <- list(fit = fit, folds = fit$folds, aucs = aucs)

  if (vimp) {
    ret_lst <- list(fit = fit$SL.predict, folds = fit$folds, aucs = aucs)
  }
  return(list(cvaucs = ret_lst, cvfits = fit))
}

# Cross-validated variable importance ------------------------------------------
# get the CV vim using pre-computed CV.SL objects
# @param seed the random number seed (only useful if (A)IPW estimates of variable importance are desired)
# @param Y the outcome (measured on all participants, including phase 1)
# @param X the marker and baseline risk variables
# @param full_fit the CV.SL object from a call to run_cv_sl_once (i.e., one cross-validated super learner)
#         based on all covariates
# @param reduced_fit the CV.SL object from a call to run_cv_sl_once (i.e., one cross-validated super learner)
#         based on the reduced set of covariates
# @param index the index of interest
# @param type the type of importance (e.g., "auc")
# @param scale the scale to measure importance/CIs on
# @param cross_fitting_folds the cross-fitting folds used to fit full_fit, reduced_fit
# @param sample_splitting_folds folds for sample splitting (splits up cross-fitting folds into unique groups)
# @param V the number of folds used for variable importance (equal to the number of folds for outer cross-fitting in the CV.SL / 2)
# @param sl_library the Super Learner library, for IPW estimation of the efficient influence function
# @param ipc_scale the scale to apply IPW or AIPW correction on
# @param ipc_est_type whether to use IPW or AIPW for the EIF
# @param ipc_weights the inverse probability weights for coarsened data
# @param Z the predictors of the missingness mechanism (a character vector)
# @param C the indicator of missing (0) or selected into phase 2 (1)
# @param baseline should we get predictiveness of the baseline risk factors (TRUE) or of the risk factors + variables of interest?
# @param use_ensemble should we use the ensemble SL (TRUE) or the discrete SL (FALSE)?
# @param final_point_estimate should the final point estimate returned be based on the sample split for inference ("split"),
#         the average of the split-specific importance ("average"), or on the full dataset ("full")? All three point estimates
#         are valid; sample-splitting is only required for inference.
# @return an object of class "vimp" with the results
get_cv_vim <- function(seed = NULL, Y = NULL, X = NULL, full_fit = NULL, reduced_fit = NULL, index = 1, type = "auc", scale = "logit",
                       cross_fitting_folds = NULL, sample_splitting_folds = NULL, V = 2, ipc_scale = "identity",
                       sl_library = c("SL.glmnet"), ipc_est_type = "ipw", ipc_weights = rep(1, length(cross_fitting_folds)), Z = NULL,
                       C = NULL, baseline = FALSE, use_ensemble = TRUE, sample_splitting = TRUE, 
                       final_point_estimate = "average") {

  # if not passing sample-splitting folds, create them
  if (all(is.null(sample_splitting_folds))) {
    set.seed(seed)
    sample_splitting_folds <- vimp::make_folds(unique(cross_fitting_folds), V = 2)
  }
  full_cv_fit <- switch(as.numeric(use_ensemble) + 1, full_fit$discreteSL.predict, full_fit$SL.predict)
  if (baseline) {
      reduced_cv_fit <- reduced_fit
  } else {
      reduced_cv_fit <- switch(as.numeric(use_ensemble) + 1, reduced_fit$discreteSL.predict, reduced_fit$SL.predict)
  }
  if (is.null(cross_fitting_folds)) {
    cross_fitting_folds <- full_cv_fit$folds
  }
  # add on CV folds for the people only sampled into phase 1
  if (is.null(Z)) {
    cf_folds <- cross_fitting_folds
  } else {
    cf_folds <- c(cross_fitting_folds, sample(seq_len(length(unique(cross_fitting_folds))), length(C) - length(cross_fitting_folds), replace = TRUE))
  }
  set.seed(seed)
  # estimate variable importance
  est_vim <- vimp::cv_vim(Y = Y, X = X, cross_fitted_f1 = full_cv_fit,
                          cross_fitted_f2 = reduced_cv_fit, indx = index,
                          delta = 0, V = V, run_regression = FALSE,
                          sample_splitting = sample_splitting, cross_fitting_folds = cf_folds,
                          sample_splitting_folds = sample_splitting_folds,
                          type = type, scale = scale, ipc_scale = ipc_scale,
                          SL.library = sl_library, ipc_est_type = ipc_est_type,
                          ipc_weights = ipc_weights, Z = Z, C = C,
                          final_point_estimate = final_point_estimate)
  # return the vimp object
  return(est_vim)
}

# pool CV-VIM estimates over multiple random starts
# @param vim_lst a list of vimp objects (from a call to get_cv_vim for each of a list of random number seeds)
# @param scale the scale that intervals should be on (should match the scale used for each individual random start)
# @return a tibble, with: the quantity (predictiveness, VIM), the point estimate, SE, 95% CI, and a p-value for variable importance
pool_cv_vim <- function (vim_lst, scale = "logit") {
  # for variable importance
  vim_point_est <- mean(unlist(lapply(vim_lst, function(l) l$est)), na.rm = TRUE)
  vim_se <- sqrt(mean(unlist(lapply(vim_lst, function(l) l$se ^ 2)), na.rm = TRUE))
  vim_ci <- vimp::vimp_ci(vim_point_est, vim_se, scale = scale, level = 0.95)
  # for predictiveness
  pred_point_est <- mean(unlist(lapply(vim_lst, function(l) l$predictiveness_full)), na.rm = TRUE)
  pred_se <- sqrt(mean(unlist(lapply(vim_lst, function(l) l$se_full ^ 2)), na.rm = TRUE))
  pred_ci <- vimp::vimp_ci(pred_point_est, pred_se, scale = scale, level = 0.95)
  pred_reduced_point_est <- mean(unlist(lapply(vim_lst, function(l) l$predictiveness_reduced)), na.rm = TRUE)
  # for variable importance
  vim_point_est <- max(pred_point_est - pred_reduced_point_est, 0)
  vim_se <- sqrt(mean(unlist(lapply(vim_lst, function(l) l$se ^ 2)), na.rm = TRUE))
  vim_ci <- vimp::vimp_ci(vim_point_est, vim_se, scale = scale, level = 0.95)
  # hypothesis test for variable importance
  vim_pval <- vimp::vimp_hypothesis_test(predictiveness_full = pred_point_est,
    predictiveness_reduced = pred_reduced_point_est, se = vim_se, delta = 0, alpha = 0.05)$p_value
  output <- tibble::tibble(s = vim_lst[[1]]$s, group = grepl(",", s), quantity = c("VIM", "Predictiveness"),
                           est = c(vim_point_est, pred_point_est), se = c(vim_se, pred_se)) %>%
    bind_cols(as_tibble(rbind(vim_ci, pred_ci), .name_repair = ~ paste0("ci_", c("ll", "ul"))), pval = c(vim_pval, NA))
  return(output)
}


##########################################################################################################
get.pca.scores <- function(dat){
  scale_each_var <- function(data){

    for (a in colnames(data)) {
      data[[a]] <- scale(data[[a]],
                         center = mean(data[[a]], na.rm = T),
                         scale = sd(data[[a]], na.rm = T))
    }
    return(data)
  }

  scaled_dat = scale_each_var(dat %>% select(-Ptid))

  dat.cov <- cov(scaled_dat)
  dat.eigen <- eigen(dat.cov)

  (phi <- dat.eigen$vectors[,1:2])

  phi <- -phi
  row.names(phi) <- names(scaled_dat)
  colnames(phi) <- c("PC1", "PC2")

  # Calculate Principal Components scores
  PC1 <- as.matrix(scaled_dat) %*% phi[,1]
  PC2 <- as.matrix(scaled_dat) %*% phi[,2]

  # Create data frame with Principal Components scores
  PC <- data.frame(Ptid = dat$Ptid, PC1, PC2)

  return(PC)
}


create_varsets <- function(markers, vars){
  varset <- rep(FALSE, length(markers))
  varset[match(vars, markers)] = TRUE
  return(varset)
}

library(mdw)
get.maxSignalDivScore <- function(dat, day){
  ## scale markers to have sd 1 for all vars
  for (a in colnames(dat)) {
    dat[[a]] <- scale(dat[[a]],
                      scale = sd(dat[[a]], na.rm = T))
  }
  # compute pairwise correlations between all marker vars
  marker.wts <- cor(dat, method = "spearman") %>%
    mdw::tree.weight(plot=FALSE)
  
  ## multiply marker values with weights
  for (a in colnames(dat)) {
    dat[[a]] <- dat[[a]] * marker.wts[a]
  }
  
  if(day %in% c("Day29", "Day57")){
    dat <- dat %>%
      mutate(max.signal.div.score = rowSums(.[1:4])) 
  }else if(day == "Day57_29"){
    dat <- dat %>%
      mutate(max.signal.div.score = rowSums(.[1:8])) 
  }else if(day == "705-Day210-primaryMarkers"){
    dat <- dat %>%
      mutate(max.signal.div.score = rowSums(.[1:6])) 
  }else if(day == "705-Day210-allMarkers"){
    dat <- dat %>%
      mutate(max.signal.div.score = rowSums(.[1:41])) 
  }
  
  return(dat$max.signal.div.score)
}


get.nonlinearPCA.scores <- function(dat){
  print("running non-linear PCA")
  ptid.vec <- dat %>% pull(Ptid)
  dat <- dat %>% select(-Ptid)

  ## scale markers to have sd 1 for all vars
  for (a in colnames(dat)) {
    dat[[a]] <- scale(dat[[a]],
                      scale = sd(dat[[a]], na.rm = T))
  }
  dat.scaled = data.frame(as.matrix(dat))

  fit=FSDAM::fsdam(dat.scaled, opt_numCode = 2)
  nlPCA1 <- kyotil::INT(fit$code[,1]) # the first component
  nlPCA2 <- kyotil::INT(fit$code[,2])

  # Create data frame with non-linear Principal Components
  nlPCA <- data.frame(Ptid = ptid.vec, nlPCA1, nlPCA2)
  print("non-linear PCA done")
  return(nlPCA)
}



##########################################################################################################

# # remove any binary risk variables with fewer than 10 ptids that have a 0 or 1 for that variable
# # @param dat the phase 1 dataset
# # @param risk_vars the vector of column names of risk variables
# # @return a data frame upon removal of any binary risk variables with fewer than 10 ptids that have a 0 or 1 for that variable
# drop_predVars_with_fewer_0s_or_1s <- function(dat, risk_vars) {
#   for (i in 1:length(risk_vars)) {
#     if ((dat %>% select(matches(risk_vars[i])) %>% unique() %>% dim())[1] == 2) {
#       if ((dim(dat %>% filter(get(risk_vars[i]) == 1))[1] < 10) | (dim(dat %>% filter(get(risk_vars[i]) == 0))[1] < 10)){
#         dat <- dat %>% select(-matches(risk_vars[i]))
#       }
#     }
#   }
#   return(dat)
# }


# @param dat the phase 1 dataset
# @param pred_vars the vector of column names of predictor variables
# @return a data frame upon removal of any binary predictor variables with fewer than 10 ptids that have a 0 or 1 for that variable (COVE analysis)
# @return a data frame upon removal of any binary predictor variables with number of cases in the variable = 1 or 0 subgroup is <= 3 (ENSEMBLE analysis)
drop_predVars_with_fewer_0s_or_1s <- function(dat, pred_vars) {
  
  if(study_name == "COVE"){
    # delete the file drop_predVars_with_fewer_0s_or_1s.csv
    unlink(here("output", Sys.getenv("TRIAL"), "drop_predVars_with_fewer_0s_or_1s.csv"))
    # Remove binary risk variables with fewer than 10 ptids that have a 0 or 1 for that variable
    for (i in 1:length(pred_vars)) {
      if ((dat %>% select(starts_with(pred_vars[i])) %>% unique() %>% dim())[1] == 2) {
        if ((dim(dat %>% filter(get(pred_vars[i]) == 1))[1] < 10) | (dim(dat %>% filter(get(pred_vars[i]) == 0))[1] < 10)){
          dat <- dat %>% select(-starts_with(pred_vars[i]))
          print(paste0(pred_vars[i], " dropped from risk score analysis as it had fewer than 10 1's or 0's."))
          # Also print to file
          paste0(pred_vars[i], " dropped from risk score analysis as it had fewer than 10 1's or 0's.") %>%
            write.table(file = here("output", Sys.getenv("TRIAL"), "drop_predVars_with_fewer_0s_or_1s.csv"), sep=",", append = TRUE, row.names = F, col.names = F)
        }
      }
    }
  }
  
  if(study_name == "COVE" & !file.exists(here("output", Sys.getenv("TRIAL"), "drop_predVars_with_fewer_0s_or_1s.csv"))){
    paste0("No binary input variable had fewer than 10 ptids with a 0 or 1 for that variable.") %>%
      write.table(file = here("output", Sys.getenv("TRIAL"), "drop_predVars_with_fewer_0s_or_1s.csv"), sep=",", append = TRUE, row.names = F, col.names = F)
  }
  
  if(study_name != "COVE"){
    # delete the file drop_predVars_with_fewer_0s_or_1s.csv
    unlink(here("output", Sys.getenv("TRIAL"), "drop_predVars_with_fewer_0s_or_1s.csv"))
    # Remove a variable if the number of cases in the variable = 1 subgroup is <= 3 or the number of cases in the variable = 0 subgroup is <= 3
    for (i in 1:length(pred_vars)) {
      if ((dat %>% select(starts_with(pred_vars[i])) %>% unique())[[1]] == c(0,1) | 
          (dat %>% select(starts_with(pred_vars[i])) %>% unique())[[1]] == c(1,0)) {
        if (dat %>% filter(get(pred_vars[i]) == 1) %>% pull(endpoint) %>% sum() <= 3 | dat %>% filter(get(pred_vars[i]) == 0) %>% pull(endpoint) %>% sum() <= 3){
          dat <- dat %>% select(-starts_with(pred_vars[i]))
          print(paste0(pred_vars[i], " dropped from risk score analysis as the number of cases in the variable = 1 or 0 subgroup is <= 3."))
          # Also print to file
          paste0(pred_vars[i], " dropped from risk score analysis as the number of cases in the variable = 1 or 0 subgroup is <= 3.") %>%
            write.table(file = here("output", Sys.getenv("TRIAL"), "drop_predVars_with_fewer_0s_or_1s.csv"), sep=",", append = TRUE, row.names = F, col.names = F)
        }
      }
    }
  }
  
  if(study_name != "COVE" & !file.exists(here("output", Sys.getenv("TRIAL"), "drop_predVars_with_fewer_0s_or_1s.csv"))){
    paste0("No binary input variable had number of cases in the variable = 1 or 0 subgroup <= 3") %>%
      write.table(file = here("output", Sys.getenv("TRIAL"), "drop_predVars_with_fewer_0s_or_1s.csv"), sep=",", append = TRUE, row.names = F, col.names = F)
  }
  
  return(dat)
}



# # remove any binary risk variables that have more than 5% values missing
# # @param X dataframe containing all risk variables
# # @param predVars the vector of column names of risk variables
# # @return a data frame upon removal of any binary risk variables that have more than 5% values missing
# drop_predVars_with_high_total_missing_values <- function(X, predVars) {
#   covars_highNAvalues <- vector()
#   for (i in 1:length(predVars)) {
#     total_NAs <- sum(is.na(X %>% pull(predVars[i])))
#     percent_NAs <- total_NAs / length(X %>% pull(predVars[i]))
# 
#     if (percent_NAs > 0.05) {
#       print(paste0("WARNING: ", predVars[i], " variable has more than 5% values missing! This variable will be dropped from SuperLearner analysis."))
#       covars_highNAvalues <- predVars[i]
#     }
#   }
# 
#   if(length(covars_highNAvalues) == 0)
#     return(X)
#   else
#     return(X %>% select(-all_of(covars_highNAvalues)))
# }



# remove any binary risk variables that have more than 5% values missing
# @param X dataframe containing all predictor variables
# @param predVars the vector of column names of predictor variables
# @return a data frame upon removal of any binary predictor variables that have more than 5% values missing
drop_predVars_with_high_total_missing_values <- function(X, predVars) {
  # delete the file drop_predVars_with_high_total_missing_values.csv
  unlink(here("output", Sys.getenv("TRIAL"), "drop_predVars_with_high_total_missing_values.csv"))
  covars_highNAvalues <- vector()
  for (i in 1:length(predVars)) {
    total_NAs <- sum(is.na(X %>% pull(predVars[i])))
    percent_NAs <- total_NAs / length(X %>% pull(predVars[i]))
    
    if (percent_NAs > 0.05) {
      print(paste0("WARNING: ", predVars[i], " variable has more than 5% values missing! This variable will be dropped from SuperLearner analysis."))
      # Also print to file
      paste0(predVars[i], " variable has more than 5% values missing and was dropped from risk score analysis.") %>%
        write.table(file = here("output", Sys.getenv("TRIAL"), "drop_predVars_with_high_total_missing_values.csv"), sep=",", append = TRUE, row.names = F, col.names = F)
      covars_highNAvalues <- predVars[i]
    }
  }
  
  if(!file.exists(here("output", Sys.getenv("TRIAL"), "drop_predVars_with_high_total_missing_values.csv"))){
    paste0("No variables had more than 5% values missing.") %>%
      write.table(file = here("output", Sys.getenv("TRIAL"), "drop_predVars_with_high_total_missing_values.csv"), sep=",", append = TRUE, row.names = F, col.names = F)
  }
  
  if(length(covars_highNAvalues) == 0)
    return(X)
  else
    return(X %>% select(-all_of(covars_highNAvalues)))
}




# # impute missing values in the risk variables using mice package
# # @param X dataframe containing all risk variables
# # @param predVars the vector of column names of risk variables
# # @return a data frame updated with imputations for missing values of risk variables
# impute_missing_values <- function(X, predVars) {
#   covars <- vector()
#   # First identify risk demographic variables having missing values
#   for (i in 1:length(predVars)) {
#     total_NAs <- sum(is.na(X %>% pull(predVars[i])))
#     percent_NAs <- total_NAs / length(X %>% pull(predVars[i]))
# 
#     if (percent_NAs > 0.05) {
#       print(paste0("WARNING: ", predVars[i], " variable has more than 5% values missing; all missing values will be imputed!"))
#     }
# 
#     if (total_NAs > 0) {
#       if (i == 1) {
#         covars <- predVars[i]
#       } else {
#         covars <- c(covars, predVars[i])
#       }
#     }
#   }
# 
#   if (length(covars) == 0) {
#     print("No missing values to impute in any risk variables!")
#   } else {
#     print(paste("Imputing missing values in following variables: ", paste(as.character(covars), collapse = ", ")))
#     n.imp <- 1
# 
#     imp <- X %>% select(all_of(covars))
# 
#     # deal with constant variables
#     for (a in names(imp)) {
#       if (all(imp[[a]]==min(imp[[a]], na.rm=TRUE), na.rm=TRUE)) imp[[a]]=min(imp[[a]], na.rm=TRUE)
#     }
# 
#     # diagnostics = FALSE , remove_collinear=F are needed to avoid errors due to collinearity
#     imp <- imp %>%
#       mice(m = n.imp, printFlag = FALSE, seed=1, diagnostics = FALSE, remove_collinear = FALSE)
# 
#     X[, covars] <- mice::complete(imp, action = 1L)
#   }
#   return(X)
# }



# impute missing values in the risk variables using mice package
# @param X dataframe containing all risk variables
# @param predVars the vector of column names of risk variables
# @return a data frame updated with imputations for missing values of risk variables
impute_missing_values <- function(X, predVars) {
  # delete the file impute_missing_values.csv
  unlink(here("output", Sys.getenv("TRIAL"), "impute_missing_values.csv"))
  covars <- vector()
  # First identify risk demographic variables having missing values
  for (i in 1:length(predVars)) {
    total_NAs <- sum(is.na(X %>% pull(predVars[i])))
    percent_NAs <- total_NAs / length(X %>% pull(predVars[i]))
    
    if (percent_NAs > 0.05) {
      print(paste0("WARNING: ", predVars[i], " variable has more than 5% values missing; all missing values will be imputed!"))
    }
    
    if (total_NAs > 0) {
      if (i == 1) {
        covars <- predVars[i]
      } else {
        covars <- c(covars, predVars[i])
      }
    }
  }
  
  if (length(covars) == 0) {
    print("No missing values to impute in any risk variables!")
    # Also print to file
    paste("Imputing missing values in following variables: None") %>%
      write.table(file = here("output", Sys.getenv("TRIAL"), "impute_missing_values.csv"), sep=",", append = TRUE, row.names = F, col.names = F)
  } else {
    print(paste("Imputing missing values in following variables: ", paste(as.character(covars), collapse = ", ")))
    # Also print to file
    paste("Missing values were imputed for the following variables:", paste(as.character(covars), collapse = ", ")) %>%
      write.table(file = here("output", Sys.getenv("TRIAL"), "impute_missing_values.csv"), sep=",", append = TRUE, row.names = F, col.names = F)
    
    n.imp <- 1
    
    impVars <- X %>% select(all_of(covars))
    
    # deal with constant variables
    for (a in names(impVars)) {
      print(a)
      if (all(impVars[[a]]==min(impVars[[a]], na.rm=TRUE), na.rm=TRUE)) {
        X[[a]] = min(impVars[[a]], na.rm=TRUE)
        covars = covars[!covars %in% a]
      }
    }
    
    noimpVars <- X %>% select(-all_of(covars)) %>% colnames()
    
    init = mice(X, maxit=0) 
    meth = init$method
    predM = init$predictorMatrix 
    meth[c(noimpVars)] = ""
    
    # diagnostics = FALSE , remove_collinear=F are needed to avoid errors due to collinearity
    X <- mice(X, method=meth, predictorMatrix=predM, m=n.imp,
              printFlag = FALSE, seed=1, diagnostics = FALSE, remove_collinear = FALSE)
    X <- mice::complete(X, action = 1L)
    
  }
  return(X)
}


# choose learners for getting predictions
# @param cvaucDAT a dataframe containing Learner, Screen, and AUC
# @return dataframe sorted according to descending AUC for individual learners followed by SL and Discrete SL
choose_learners <- function(cvaucDAT) {
  cvaucDAT %>%
    filter(!Learner %in% c("SL", "Discrete SL")) %>%
    arrange(-AUC) %>%
    .[1:2, ] %>%
    bind_rows(cvaucDAT %>%
      filter(Learner %in% c("SL", "Discrete SL")))
}

# get CV predictions
# @param cv_fit fit from running CV.Superlearner
# @param cvaucDAT a dataframe containing Learner, Screen, and AUC
# @return a dataframe with predictions for each subject for each Learner/Learner-Screen combination present in cvaucDAT
get_cv_predictions <- function(cv_fit, cvaucDAT, markerDAT = NULL) {
  top3 <- choose_learners(cvaucDAT)

  predict <- cv_fit[["library.predict"]] %>%
    as.data.frame() %>%
    bind_cols(cv_fit[["discreteSL.predict"]] %>% as.data.frame() %>% `colnames<-`(c("Discrete SL"))) %>%
    bind_cols(cv_fit[["SL.predict"]] %>% as.data.frame() %>% `colnames<-`(c("Super Learner"))) %>%
    bind_cols(cv_fit[["Y"]] %>% as.data.frame() %>% `colnames<-`(c("Y")))

  if(!is.null(markerDAT)){
    predict <- predict %>% bind_cols(markerDAT)
  }

  predict <- predict %>%
    gather("algo", "pred", -c(Y, names(markerDAT))) %>%
    filter(algo %in% c(top3$LearnerScreen))

  predict %>%
    left_join(top3 %>% select(Screen_fromRun, Learner, Screen, AUC, LearnerScreen), by = c("algo" = "LearnerScreen")) %>%
    mutate(
      learnerScreen = paste0(Learner, "_", Screen),
      learnerScreen = ifelse(Learner %in% c("SL", "Discrete SL"), algo, learnerScreen),
      AUCchar = format(round(AUC, 3), nsmall = 3),
      learnerScreen = paste0(learnerScreen, " (", AUCchar, ")"),
      learnerScreen = reorder(learnerScreen, -AUC)
    )
}


# Compute inverse probability weighted true and false positive rate, for CV-ROC curves
# @param pred_obj an ROCR prediction object
# @param weights the inverse probability weights
# @return the IPW CV-ROC curve
ipw_roc_curves <- function(pred_obj, weights = rep(1, length(pred_obj@predictions[[1]]))) {
  y <- pred_obj@labels[[1]]
  pos_label <- levels(pred_obj@labels[[1]])[2]
  neg_label <- levels(pred_obj@labels[[1]])[1]
  # compute weighted num pos and neg
  n_0_weighted <- sum((y == neg_label) * weights)
  n_1_weighted <- sum((y == pos_label) * weights)
  # compute weighted true and false positive rates
  preds <- as.numeric(pred_obj@predictions[[1]])
  pred_order <- order(preds, decreasing = TRUE)
  ordered_preds <- preds[pred_order]
  weighted_tp <- cumsum(weights[pred_order] * (y[pred_order] == pos_label))
  weighted_fp <- cumsum(weights[pred_order] * (y[pred_order] == neg_label))
  dups <- rev(duplicated(rev(ordered_preds)))
  fp <- c(0, weighted_fp[!dups])
  tp <- c(0, weighted_tp[!dups])
  pred_obj@fp[[1]] <- fp
  pred_obj@tp[[1]] <- tp
  pred_obj@n.pos[[1]] <- n_1_weighted
  pred_obj@n.neg[[1]] <- n_0_weighted
  pred_obj
}



# Plot ROC curves for SL, discrete.SL and topRanking learner-screen combinations
# @param predict dataframe returned by get_cv_predictions function
# @param cvaucDAT a dataframe containing Learner, Screen, and AUC
# @param weights the inverse probability weights for the participants observed in phase 2
# @return ggplot object containing the ROC curves
plot_roc_curves <- function(predict, cvaucDAT, weights) {
  top3 <- choose_learners(cvaucDAT)

  roc.obj <- predict %>%
    group_by(algo) %>%
    nest() %>%
    mutate(
      pred.obj = purrr::map(data, ~ ROCR::prediction(.x$pred, .x$Y)),
      weighted_pred_obj = purrr::map(pred.obj, ~ ipw_roc_curves(.x, weights)),
      perf.obj = purrr::map(weighted_pred_obj, ~ ROCR::performance(.x, "tpr", "fpr")),
      roc.dat = purrr::map(perf.obj, ~ tibble(
        xval = .x@x.values[[1]],
        yval = .x@y.values[[1]]
      ))
    )

  roc.obj %>%
    unnest(roc.dat) %>%
    select(algo, xval, yval) %>%
    ungroup() %>%
    left_join(top3 %>% select(Screen_fromRun, Learner, Screen, AUC, LearnerScreen), by = c("algo" = "LearnerScreen")) %>%
    mutate(
      learnerScreen = paste0(Learner, "_", Screen),
      learnerScreen = ifelse(Learner %in% c("SL", "Discrete SL"), algo, learnerScreen),
      AUCchar = format(round(AUC, 3), nsmall = 3),
      learnerScreen = paste0(learnerScreen, " (", AUCchar, ")"),
      learnerScreen = reorder(learnerScreen, -AUC)
    ) %>%
    ggplot(aes(x = xval, y = yval, col = learnerScreen, linetype = learnerScreen)) +
    geom_step(lwd = 2) +
    theme(
      legend.position = "top",
      legend.direction = "vertical",
      legend.box = "horizontal",
      legend.title=element_text(size=20),
      legend.text=element_text(size=20),
      axis.ticks.length = unit(.35, "cm"),
      axis.text = element_text(size = 23),
      axis.title = element_text(size = 30)
    ) +
    labs(x = "Cross-Validated False Positive Rate", y = "Cross-Validated True Positive Rate",
         col = "Model (CV-AUC)", linetype = "Model (CV-AUC)") +
    geom_abline(intercept = 0, slope = 1) +
    scale_linetype_manual(values=c("dashed", "dotted", "dotdash", "twodash"))
}


# Plot predicted probability plots for SL, discrete.SL and topRanking learner-screen combinations
# @param pred dataframe returned by get_cv_predictions function
# @param weights the inverse probability weights
# @return ggplot object containing the predicted probability plots
plot_predicted_probabilities <- function(pred, weights = rep(1, length(pred$pred))) {
  if(study_name %in% c("COVE", "MockCOVE")){
    cases = "Post Day 57 Cases"
    disease_name = "COVID-19"
  }
  if(study_name == "ENSEMBLE"){
    cases = "Post Day 29 Cases"
    disease_name = "COVID-19"
  }
  if(study_name == "HVTN705"){
    cases = "Post Day 210 Cases"
    disease_name = "HIV"
  }
  pred %>%
    mutate(Ychar = ifelse(Y == 0, "Non-Cases", cases)) %>%
    ggplot(aes(x = Ychar, y = pred, color = Ychar)) +
    geom_jitter(width = 0.06, size = 3, shape = 21, fill = "white") +
    geom_violin(alpha = 0.05, color = "black", lwd=1.0) +
    geom_boxplot(alpha = 0.05, width = 0.15, color = "black", outlier.size = NA, outlier.shape = NA, lwd=1.0) +
    theme_bw() +
    #scale_color_manual(values = c("#56B4E9", "#E69F00")) +
    scale_color_manual(values = c("#00468B", "#8B0000")) +
    facet_wrap(vars(learnerScreen), ncol = 1, scales = "free") +
    labs(y = paste0("CV estimated predicted probability of ", disease_name, " disease"), x = "") +
    theme(
      legend.position = "none",
      strip.text.x = element_text(size = 25),
      axis.text = element_text(size = 23),
      axis.ticks.length = unit(.35, "cm"),
      axis.title.y = element_text(size = 30)
    )
}






# Create forest plot for demo run
# @param avgs dataframe containing Screen, Learner, AUC estimates and CIs as columns
# @return list of 2 ggplot objects: one containing forest plot and the other containing labels (Screen, Learner and CV-AUCs)
make_forest_plot_demo <- function(avgs) {
  lowestXTick <- floor(min(avgs$ci_ll) * 10) / 10
  top_learner_plot <- ggplot() +
    geom_pointrange(avgs %>% mutate(LearnerScreen = fct_reorder(LearnerScreen, AUC, .desc = F)), mapping = aes(x = LearnerScreen, y = AUC, ymin = ci_ll, ymax = ci_ul), size = 1, color = "blue", fill = "blue", shape = 20) +
    coord_flip() +
    scale_y_continuous(breaks = seq(lowestXTick, 1, 0.1), labels = seq(lowestXTick, 1, 0.1), limits = c(lowestXTick, 1)) +
    theme_bw() +
    labs(y = "CV-AUC [95% CI]", x = "") +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(size = 16),
      axis.title.x = element_text(size = 16),
      axis.text.y = element_blank(),
      plot.margin = unit(c(1, -0.15, 1, -0.15), "cm"),
      panel.border = element_blank(),
      axis.line = element_line(colour = "black")
    )

  total_learnerScreen_combos <- length(avgs$LearnerScreen)

  avgs_withCoord <- avgs %>%
    select(Learner, Screen, AUCstr) %>%
    gather("columnVal", "strDisplay") %>%
    mutate(
      xcoord = case_when(
        columnVal == "Learner" ~ 1,
        columnVal == "Screen" ~ 1.5,
        columnVal == "AUCstr" ~ 2
      ),
      ycoord = rep(total_learnerScreen_combos:1, 3)
    )

  top_learner_nms_plot <- ggplot(avgs_withCoord, aes(x = xcoord, y = ycoord, label = strDisplay)) +
    geom_text(hjust = 1, vjust = 0, size = 5) +
    xlim(0.7, 2) +
    theme(
      plot.margin = unit(c(2.2, -0.15, 3, -0.15), "cm"),
      axis.line = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 2, color = "white"),
      axis.ticks = element_blank(),
      axis.title = element_blank()
    )

  return(list(top_learner_plot = top_learner_plot, top_learner_nms_plot = top_learner_nms_plot))
}

# Create forest plot for all variable sets for either SL or discrete SL
# @param cvaucs_vacc dataframe containing Variable Sets, Screen, Learner, AUC estimates and CIs as columns
# @return list of 2 ggplot objects: one containing forest plot and the other containing labels (Variable Sets and CV-AUCs)
make_forest_plot_SL_allVarSets <- function(dat, learner.choice = "SL"){
  allSLs <- dat %>% filter(Learner == learner.choice) %>%
    mutate(varsetNo = sapply(strsplit(varset, "_"), `[`, 1),
           varsetNo = as.numeric(varsetNo)) %>%
    arrange(varsetNo) %>%
    mutate(varset = fct_reorder(varset, AUC, .desc = F)) %>%
    arrange(-AUC)

  lowestXTick <- floor(min(allSLs$ci_ll)*10)/10
  highestXTick <- ceiling(max(allSLs$ci_ul)*10)/10
  if(Sys.getenv("TRIAL") == "hvtn705second"){
    learner.plot.margin = unit(c(2.25,0.2,0.8,-0.15),"cm")
    names.plot.margin = unit(c(0.1,-0.15,0.65,-0.15),"cm")
    y_at = 35.2
  } else if(Sys.getenv("TRIAL") == "moderna_real"){
    learner.plot.margin = unit(c(0.8,0.2,0.8,-0.15),"cm")
    names.plot.margin = unit(c(-1.5,-0.15,0.65,-0.15),"cm")
    y_at = 35.2
  } else if(Sys.getenv("TRIAL") %in% c("janssen_pooled_partA", "janssen_la_partA")){
    learner.plot.margin = unit(c(3.25,0.2,0.8,-0.15),"cm")
    names.plot.margin = unit(c(1.6,-0.15,1.39,-0.15),"cm")
    y_at = 16.7
  }
  
  top_learner_plot <- ggplot() +
    geom_pointrange(allSLs %>% mutate(varset = fct_reorder(varset, AUC, .desc = F)), mapping=aes(x=varset, y=AUC, ymin=ci_ll, ymax=ci_ul), size = 1, color="blue", fill="blue", shape=20) +
    coord_flip() +
    scale_y_continuous(breaks = seq(lowestXTick, max(highestXTick, 1), 0.1), labels = seq(lowestXTick, max(highestXTick, 1), 0.1), limits = c(lowestXTick, max(highestXTick, 1))) +
    theme_bw() +
    labs(y = "CV-AUC [95% CI]", x = "") +
    geom_hline(yintercept = 0.5, lty = "dashed") +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size=16),
          axis.title.x = element_text(size=16),
          axis.text.y = element_blank(),
          plot.margin=learner.plot.margin,
          panel.border = element_blank(),
          axis.line = element_line(colour = "black"))

  total_learnerScreen_combos = length(allSLs$LearnerScreen)

  allSLs_withCoord <- allSLs %>%
    select(varset, AUCstr) %>%
    mutate(varset = as.character(varset)) %>%
    tidyr::gather("columnVal", "strDisplay") %>%
    mutate(xcoord = case_when(columnVal=="varset" ~ 1.5,
                              columnVal=="AUCstr" ~ 2),
           ycoord = rep(total_learnerScreen_combos:1, 2))
  
  top_learner_nms_plot <- ggplot(allSLs_withCoord, aes(x = xcoord, y = ycoord, label = strDisplay)) +
    geom_text(hjust=1, vjust=0, size=5) +
    xlim(0.7, 2) +
    theme_bw() +
    theme(plot.margin=names.plot.margin,
          axis.line = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(size = 2, color = "white"),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          panel.border = element_blank()) +
    annotate("text", x = 1.5, y = y_at, size = 5,
             label = "Variable Set",
             fontface = "bold",
             hjust = 1) +
    annotate("text", x = 2, y = y_at, size = 5,
             label = "CV-AUC [95% CI]",
             fontface = "bold",
             hjust = 1)
  
  return(list(top_learner_plot = top_learner_plot, top_learner_nms_plot = top_learner_nms_plot))
}


# Create forest plot for prod run
# @param avgs dataframe containing Screen, Learner, AUC estimates and CIs as columns
# @return list of 2 ggplot objects: one containing forest plot and the other containing labels (Screen, Learner and CV-AUCs)
make_forest_plot_prod <- function(avgs) {
  
  if(Sys.getenv("TRIAL") == "moderna_real"){
    PLOT.MARGIN = unit(c(0.8, -0.15, 0.2, -0.15), "cm")
    NAMES.PLOT.MARGIN = unit(c(1.0, -0.15, 1.7, -0.15), "cm")
  } else if(Sys.getenv("TRIAL") == "janssen_pooled_partA"){
    PLOT.MARGIN = unit(c(0.8, -0.15, 0.2, -0.15), "cm")
    NAMES.PLOT.MARGIN = unit(c(1.0, -0.15, 1.7, -0.15), "cm")
  }
  
  lowestXTick <- floor(min(avgs$ci_ll) * 10) / 10
  top_learner_plot <- ggplot() +
    geom_pointrange(avgs %>% mutate(LearnerScreen = fct_reorder(LearnerScreen, AUC, .desc = F)), mapping = aes(x = LearnerScreen, y = AUC, ymin = ci_ll, ymax = ci_ul), size = 3, color = "blue", fill = "blue", shape = 20) +
    coord_flip() +
    scale_y_continuous(breaks = seq(lowestXTick, 1, 0.1), labels = seq(lowestXTick, 1, 0.1), limits = c(lowestXTick, 1)) +
    theme_bw() +
    labs(y = "CV-AUC [95% CI]", x = "") +
    theme(
      axis.title.y = element_blank(),
      axis.text = element_text(size = 30),
      axis.title = element_text(size = 30),
      axis.ticks.length = unit(.35, "cm"),
      axis.text.y = element_blank(),
      plot.margin = PLOT.MARGIN,
      panel.border = element_blank(),
      axis.line = element_line(colour = "black")
    )

  total_learnerScreen_combos <- length(avgs$LearnerScreen)

  avgs_withCoord <- avgs %>%
    select(Learner, Screen, AUCstr) %>%
    gather("columnVal", "strDisplay") %>%
    mutate(
      xcoord = case_when(
        columnVal == "Learner" ~ 2,
        columnVal == "Screen" ~ 1.5,
        columnVal == "AUCstr" ~ 2
      ),
      ycoord = rep(total_learnerScreen_combos:1, 3)
    )

  top_learner_nms_plot <- ggplot(avgs_withCoord, aes(x = xcoord, y = ycoord, label = strDisplay)) +
    geom_text(hjust = 1, vjust = 0, size = 10) +
    xlim(0.7, 2) +
    theme(
      plot.margin = NAMES.PLOT.MARGIN,
      axis.line = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 2, color = "white"),
      axis.ticks = element_blank(),
      axis.title = element_blank()
    )

  return(list(top_learner_plot = top_learner_plot, top_learner_nms_plot = top_learner_nms_plot))
}


# Create forest plot in general
# @param avgs dataframe containing Screen, Learner, AUC estimates and CIs as columns
# @return list of 2 ggplot objects: one containing forest plot and the other containing labels (Screen, Learner and CV-AUCs)
make_forest_plot <- function(avgs){

  if(Sys.getenv("TRIAL") == "hvtn705second"){
    PLOT.MARGIN = unit(c(3.4,-0.15,1,-0.15),"cm")
    NAMES.PLOT.MARGIN = unit(c(1.6,-0.15,1.5,-0.15),"cm")
    y_at = 18.75
  } else if(Sys.getenv("TRIAL") == "moderna_real"){
    PLOT.MARGIN = unit(c(2.8,-0.15,1,-0.15),"cm")
    NAMES.PLOT.MARGIN = unit(c(1.0,-0.15,1.7,-0.15),"cm")
    y_at = 15.75
  } else if(Sys.getenv("TRIAL") %in% c("janssen_pooled_partA", "janssen_la_partA")){
    PLOT.MARGIN = unit(c(2.8,0.15,1,-0.15),"cm")
    NAMES.PLOT.MARGIN = unit(c(0.6,-0.15,0.85,-0.15),"cm")
    y_at = 31.25
  }
  
  lowestXTick <- floor(min(avgs$ci_ll)*10)/10
  highestXTick <- ceiling(max(avgs$ci_ul)*10)/10

  top_learner_plot <- ggplot() +
    geom_pointrange(avgs %>% mutate(LearnerScreen = fct_reorder(LearnerScreen, AUC, .desc = F)), mapping=aes(x=LearnerScreen, y=AUC, ymin=ci_ll, ymax=ci_ul), size = 1, color="blue", fill="blue", shape=20) +
    coord_flip() +
    scale_y_continuous(breaks = seq(lowestXTick, max(highestXTick, 1), 0.1), labels = seq(lowestXTick, max(highestXTick, 1), 0.1), limits = c(lowestXTick, max(highestXTick, 1))) +
    theme_bw() +
    labs(y = "CV-AUC [95% CI]", x = "") +
    geom_hline(yintercept = 0.5, lty = "dashed") +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size=16),
          axis.title.x = element_text(size=16),
          axis.text.y = element_blank(),
          plot.margin = PLOT.MARGIN,
          panel.border = element_blank(),
          axis.line = element_line(colour = "black"))

  total_learnerScreen_combos = length(avgs$LearnerScreen)

  avgs_withCoord <- avgs %>%
    select(Learner, Screen, AUCstr) %>%
    gather("columnVal", "strDisplay") %>%
    mutate(xcoord = case_when(columnVal=="Learner" ~ 1,
                              columnVal=="Screen" ~ 1.5,
                              columnVal=="AUCstr" ~ 2),
           ycoord = rep(total_learnerScreen_combos:1, 3))


  top_learner_nms_plot <- ggplot(avgs_withCoord, aes(x = xcoord, y = ycoord, label = strDisplay)) +
    geom_text(hjust=1, vjust=0, size=5) +
    xlim(0.7,2) +
    theme(plot.margin = NAMES.PLOT.MARGIN,
          axis.line=element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(size = 2, color = "white"),
          axis.ticks = element_blank(),
          axis.title = element_blank()) +
    annotate("text", x = 1, y = y_at, size = 5,
             label = "Learner",
             fontface = "bold",
             hjust = 1)  +
    annotate("text", x = 1.5, y = y_at, size = 5,
             label = "Screen",
             fontface = "bold",
             hjust = 1) +
    annotate("text", x = 2, y = y_at, size = 5,
             label = "CV-AUC [95% CI]",
             fontface = "bold",
             hjust = 1)

  return(list(top_learner_plot = top_learner_plot, top_learner_nms_plot = top_learner_nms_plot))
}
