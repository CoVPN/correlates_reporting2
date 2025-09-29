#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
    
# There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))
    
source(here::here("..", "_common.R"))
#-----------------------------------------------

## create SL screens, algorithm + screen combinations

## -------------------------------------------------------------------------------------
## SL screens; all models adjust for baseline maternal enrollment variables
## -------------------------------------------------------------------------------------
## screen based on logistic regression univariate p-value < level
rank_univariate_logistic_pval <- function(Y, X, family, obsWeights, id, ...) {
  ## logistic regression of outcome on each variable
  listp <- apply(X, 2, function(x, Y, family) {
    summ <- coef(summary(glm(as.formula(briskfactors_correction),
                             family = family, weights = obsWeights
    )))
    ifelse(dim(summ)[1] > 1, summ[2, 4], 1)
  }, Y = Y, family = family)
  ## rank the p-values
  ranked_vars <- rank(listp, ties = "average")
  # Give briskfactors the lowest rank 
  # (will always set to TRUE anyways)
  ranked_vars[names(X) %in% briskfactors] <- 999
  return(ranked_vars)
}




# no screen (only the top 6 covariates with lowest univariate logistic pvalue)
screen_all <- function(Y, X, family, obsWeights, id, nVar = maxVar, ...) {
  # X contain baseline variables
  ## logistic regression of outcome on each variable
  vars <- rep(TRUE, ncol(X))
  listp <- apply(X, 2, function(x, Y, family) {
    summ <- coef(summary(glm(as.formula(briskfactors_correction),
                             family = family, weights = obsWeights
    )))
    ifelse(dim(summ)[1] > 1, summ[2, 4], 1)
  }, Y = Y, family = family)

  rankedVars <- rank(listp, ties = "average")
  listp <- rep(TRUE, length(listp))
  listp[listp][rankedVars > nVar] <- FALSE
  vars <- listp
  names(vars) <- names(X)
  # always keep the briskfactors columns of X 
  vars[names(X) %in% briskfactors] <- TRUE
  updatedX = X[vars]
  # Weed out highly-correlated vars at random
  low_cor_vars = screen_highcor_random(Y, X = updatedX, family, obsWeights, id, nVar = maxVar)
  
  # Step 1: Identify missing variable names
  missing_names <- setdiff(names(vars), names(low_cor_vars))
  
  # Step 2: Create a logical vector for the missing variables (set to FALSE or NA)
  # Here we set to FALSE:
  missing_vals <- setNames(rep(FALSE, length(missing_names)), missing_names)
  
  # Step 3: Combine the existing and missing values
  low_cor_vars_expanded <- c(low_cor_vars, missing_vals)
  
  # Step 4: Reorder to match the order of vars
  low_cor_vars_expanded <- low_cor_vars_expanded[names(vars)]
  
  # Create data frame
  df_compare <- data.frame(
    variable = names(vars),
    vars = as.logical(vars),
    low_cor_vars = low_cor_vars[names(vars)],  # will introduce NAs if missing
    low_cor_vars_expanded = as.logical(low_cor_vars_expanded),
    row.names = NULL
  )
  
#  print(paste0("vars is length: ", length(vars)))
#  print(vars)
  return(low_cor_vars_expanded)
}



## screen based on lasso
# vars <- screen_glmnet(Y = Y, X = X_mat, family = family, obsWeights = obsWeights, id = NULL, alpha = 1,
#                       minscreen = 2, nfolds = 10, nlambda = 100, nVar = maxVar)
screen_glmnet <- function(Y, X, family, obsWeights, id,
                          alpha = 1,
                          minscreen = 2,
                          nfolds = 10,
                          nlambda = 100,
                          nVar = maxVar,
                          ...) {
  
  # Check for any NULL columns in X
  #stopifnot(all(sapply(X, function(x) !is.null(x))))
  
  set.seed(123)
  
  # Step 1: Correlation-based screening
  highcor_vars <- screen_highcor_random(Y, X, family, obsWeights, id, nVar = nVar)
  
  # Make sure highcor_vars is a logical vector with names
  if (is.null(names(highcor_vars))) {
    names(highcor_vars) <- colnames(X)
  }
  
  # Step 2: Subset X to those selected by screen_highcor_random
  updatedX <- X[, highcor_vars, drop = FALSE]
  
  # Step 3: Run glmnet screen on reduced X
  glmnet_selected <- screen.glmnet(
    Y = Y,
    X = updatedX,
    family = binomial(), 
    obsWeights = obsWeights,
    id = id,
    alpha = alpha,
    minscreen = minscreen,
    nfolds = nfolds,
    nlambda = nlambda,
    ...
  )
  
  # glmnet_selected is a logical vector indicating selected columns of updatedX
  # Build a logical vector of length ncol(X), FALSE by default
  glmnet_vars <- setNames(rep(FALSE, ncol(X)), colnames(X))
  glmnet_vars[names(updatedX)[glmnet_selected]] <- TRUE
  
  # Step 4: Combine glmnet_vars with briskfactors - ensure they are TRUE
  full_vars <- glmnet_vars
  full_vars[briskfactors] <- TRUE
  
  # # Safety: Ensure names and length
  # stopifnot(!is.null(names(full_vars)))
  # stopifnot(length(full_vars) == ncol(X))
  # stopifnot(is.logical(full_vars))
  
  # Optional debugging / comparison data frame (can be returned or printed if desired)
  df_compare <- data.frame(
    variable = colnames(X),
    highcor_vars = highcor_vars,
    glmnet_vars = glmnet_vars,
    full_vars = full_vars,
    stringsAsFactors = FALSE
  )
  
  # Step 5: Limit number of immune markers by ranking p-values
  # Subset X to variables marked TRUE in full_vars
  selected_marker_names <- names(full_vars)[full_vars]
  
  # Include briskfactors regardless
  selected_marker_names <- unique(c(selected_marker_names, briskfactors))
  
  # Guard against missing columns in X
  selected_marker_names <- intersect(selected_marker_names, colnames(X))
  
  X_screened <- X[, selected_marker_names, drop = FALSE]
  
  ranked_vars <- rank_univariate_logistic_pval(Y, X_screened, family, obsWeights, id)
  
  # Add ranked_vars to df_compare (with NA for variables not in ranked_vars)
  df_compare$ranked_vars <- ranked_vars[match(df_compare$variable, names(ranked_vars))]
  
  # Keep only top nVar variables by p-value
  df_compare$final_vars <- with(df_compare, 
                                ifelse(is.na(ranked_vars), FALSE,
                                       ranked_vars == 999 | (rank(ranked_vars, ties.method = "first") <= nVar & ranked_vars != 999))
  )
  
  final_vars <- setNames(df_compare$final_vars, df_compare$variable)
  
  # Final safety checks again
  # stopifnot(length(final_vars) == ncol(X))
  # stopifnot(!is.null(names(final_vars)))
  # stopifnot(is.logical(final_vars))
  
  return(final_vars)
}





# 
# 
# ## screen based on logistic regression univariate p-value < level
# screen_univariate_logistic_pval <- function(Y, X, family, obsWeights, id,
#                                             minPvalue = 0.1, minscreen = 2,
#                                             nVar = maxVar, ...) {
#   
#   # Step 1: Correlation-based screening
#   highcor_vars <- screen_highcor_random(Y, X, family, obsWeights, id, nVar = maxVar)
#   
#   # Step 2: Subset X to those selected by screen_highcor_random
#   updatedX <- X[, highcor_vars, drop = FALSE]
#   
#   ## logistic regression of outcome on each variable
#   listp <- apply(updatedX, 2, function(x, Y, family) {
#     summ <- coef(summary(glm(as.formula(briskfactors_correction),
#                              family = family, weights = obsWeights
#     )))
#     ifelse(dim(summ)[1] > 1, summ[2, 4], 1)
#   }, Y = Y, family = family)
#   vars <- (listp <= minPvalue)
#   # always keep the briskfactors columns of X 
#   vars[names(X) %in% briskfactors] <- TRUE
#   
#   if (sum(vars) < minscreen) {
#     warning("number of variables with p value less than minPvalue is less than minscreen")
#     vars[rank(listp) <= minscreen] <- TRUE
#   }
#   # keep only a max of nVar immune markers; rank by univariate p-value
#   X_initial_screen <- X %>%
#     select(names(X)[vars], all_of(briskfactors))
#   ranked_vars <- rank_univariate_logistic_pval(Y, X_initial_screen, family,
#                                                obsWeights, id)
# 
#   vars[vars][ranked_vars > nVar] <- FALSE
#   # always keep the briskfactors columns of X 
#   vars[names(X) %in% briskfactors] <- TRUE
#   
# #  print(paste0("vars is length: ", length(vars)))
# #  print(vars)
#   return(vars)
# }


screen_univariate_logistic_pval <- function(Y, X, family, obsWeights, id,
                                            minPvalue = 0.1,
                                            minscreen = 2,
                                            nVar = maxVar,
                                            ...) {
  # Step 1: Correlation-based screening
  highcor_vars <- screen_highcor_random(Y, X, family, obsWeights, id, nVar = maxVar)
  
  if (is.null(names(highcor_vars))) {
    names(highcor_vars) <- colnames(X)
  }
  
  # Step 2: Subset X to those selected by screen_highcor_random
  updatedX <- X[, highcor_vars, drop = FALSE]
  
  # Step 3: Univariate logistic regression of outcome on each variable
  listp <- apply(updatedX, 2, function(x, Y, family) {
    summ <- coef(summary(glm(Y ~ x, family = family, weights = obsWeights)))
    if (nrow(summ) > 1) summ[2, 4] else 1
  }, Y = Y, family = family)
  
  vars <- (listp <= minPvalue)
  
  # Always keep briskfactors
  vars[names(updatedX) %in% briskfactors] <- TRUE
  
  # Expand highcor_vars to match colnames(X)
  highcor_vars_expanded <- setNames(rep(FALSE, length(colnames(X))), colnames(X))
  highcor_vars_expanded[names(highcor_vars)] <- highcor_vars
  
  # Expand vars to match colnames(X)
  full_vars <- setNames(rep(FALSE, length(colnames(X))), colnames(X))
  full_vars[names(vars)] <- vars
  
  # Ensure minimum number of variables retained
  if (sum(full_vars) < minscreen) {
    warning("number of variables with p-value less than minPvalue is less than minscreen")
    ranked_all <- rank(listp, ties.method = "first")
    full_vars[ranked_all <= minscreen] <- TRUE
  }
  
  # Create data frame
  df_compare <- data.frame(
    variable = colnames(X),
    highcor_vars = highcor_vars_expanded,
    full_vars = full_vars,
    stringsAsFactors = FALSE
  )
  
  # Step 5: Limit number of immune markers by ranking p-values
  # Subset X to variables marked TRUE in full_vars
  selected_marker_names <- names(full_vars)[full_vars]
  
  # Include briskfactors regardless
  selected_marker_names <- unique(c(selected_marker_names, briskfactors))
  
  # Guard against missing columns in X
  selected_marker_names <- intersect(selected_marker_names, colnames(X))
  
  X_screened <- X[, selected_marker_names, drop = FALSE]
  
  ranked_vars <- rank_univariate_logistic_pval(Y, X_screened, family, obsWeights, id)
  
  # Add ranked_vars to df_compare (with NA for variables not in ranked_vars)
  df_compare$ranked_vars <- ranked_vars[match(df_compare$variable, names(ranked_vars))]
  
  # Keep only top nVar variables by p-value
  df_compare$final_vars <- with(df_compare, 
                                ifelse(is.na(ranked_vars), FALSE,
                                       ranked_vars == 999 | (rank(ranked_vars, ties.method = "first") <= nVar & ranked_vars != 999))
  )
  
  final_vars <- setNames(df_compare$final_vars, df_compare$variable)
  
  # Final safety checks again
  stopifnot(length(final_vars) == ncol(X))
  stopifnot(!is.null(names(final_vars)))
  stopifnot(is.logical(final_vars))
  
  return(final_vars)
}




# 
# 
# ## screen to avoid high-correlation amongst risk variables
# screen_highcor_random <- function(Y, X, family, obsWeights, id, nVar = maxVar, ...) {
#   
#   # Step 1: The get_priorities function
#   if(study_name != "HVTN705"){
#     # Function to select markers based off priorities
#     get_priorities <- function(x){
#       if((str_detect(x[1], "mn50") & str_detect(x[2], "pseudo")) | (str_detect(x[2], "mn50") & str_detect(x[1], "pseudo"))){
#         selectVar = if_else(str_detect(x[1], "mn50"), x[1], x[2]) 
#       } else if((str_detect(x[1], "pseudo") & str_detect(x[2], "bind")) | (str_detect(x[2], "pseudo") & str_detect(x[1], "bind"))){
#         selectVar = if_else(str_detect(x[1], "pseudo"), x[1], x[2])
#       } else if((str_detect(x[1], "Spike") & str_detect(x[2], "RBD")) | (str_detect(x[2], "Spike") & str_detect(x[1], "RBD"))){
#         selectVar = if_else(str_detect(x[1], "RBD"), x[1], x[2])
#       } else if((str_detect(x[1], "id50") & str_detect(x[2], "id80")) | (str_detect(x[2], "id50") & str_detect(x[1], "id80"))){
#         selectVar = if_else(str_detect(x[1], "id50"), x[1], x[2])
#       } else {
#         selectVar = sample(c(x[1], x[2]), 1, replace = TRUE)
#       }
#       return(selectVar)
#     }
#   } else if(study_name == "HVTN705"){
#     # Function to select markers based off priorities
#     get_priorities <- function(x){
#       selectVar = sample(c(x[1], x[2]), 1, replace = TRUE)
#       return(selectVar)
#     }
#   }
#   
#   all_vars = colnames(X)
#   # set all vars to FALSE
#   vars <- rep(FALSE, ncol(X))
#   # compute pairwise correlations between all marker vars
#   cors <- cor(X, method = "spearman")
#   diag(cors) <- NA
#   
#   # Step 2: Get high-correlation pairs (> 0.9), fix factor-to-character issue
#   high_corr_pairs <- as.data.frame(as.table(cors), stringsAsFactors = FALSE) %>%
#     mutate(
#       var1 = as.character(Var1),
#       var2 = as.character(Var2),
#       corr = Freq
#     ) %>%
#     filter(!is.na(corr)) %>%
#     filter(var1 != var2) %>%
#     filter(abs(corr) > 0.9) %>%
#     rowwise() %>%
#     mutate(pair = paste(sort(c(var1, var2)), collapse = "_")) %>%
#     ungroup() %>%
#     distinct(pair, .keep_all = TRUE) %>%
#     select(var1, var2)
#   
#   # Step 3: Process each pair independently
#   set.seed(123)  # for reproducibility
#   
#   selected <- briskfactors
#   dropped <- character()
#   
#   for (i in seq_len(nrow(high_corr_pairs))) {
#     a <- high_corr_pairs$var1[i]
#     b <- high_corr_pairs$var2[i]
#     
#     if (a %in% dropped | b %in% dropped) {
#       next  # one already dropped
#     }
#     
#     if (a %in% selected & !(b %in% selected)) {
#       dropped <- c(dropped, b)
#     } else if (b %in% selected & !(a %in% selected)) {
#       dropped <- c(dropped, a)
#     } else if (a %in% selected & b %in% selected) {
#       # both selected, drop one (arbitrarily or based on get_priorities)
#       keep <- get_priorities(c(a, b))
#       drop <- setdiff(c(a, b), keep)
#       selected <- setdiff(selected, drop)
#       dropped <- c(dropped, drop)
#     } else {
#       keep <- get_priorities(c(a, b))
#       drop <- setdiff(c(a, b), keep)
#       selected <- c(selected, keep)
#       dropped <- c(dropped, drop)
#     }
#   }
#   
#   
#   # Step 4: Add variables not involved in any high-correlation pair
#   involved_vars <- unique(c(high_corr_pairs$var1, high_corr_pairs$var2))
#   uninvolved_vars <- setdiff(all_vars, involved_vars)
#   uninvolved_vars <- setdiff(uninvolved_vars, dropped)
#   
#   selected <- unique(c(selected, uninvolved_vars))
#   selected <- sort(selected)
#   
#   # Output
#   #print(selected)
#   vars <- setNames(all_vars %in% selected, all_vars)
#   
#   # Final check: ensure no pair in `selected` has correlation > 0.9
#   check_matrix <- abs(cors[selected, selected])
#   check_matrix[lower.tri(check_matrix, diag = TRUE)] <- NA
#   violations <- which(check_matrix > 0.9, arr.ind = TRUE)
#   
#   if (nrow(violations) > 0) {
#     msg <- "Correlation violations found among selected variables:\n"
#     viol_str <- apply(violations, 1, function(idx) {
#       sprintf("  %s ~ %s = %.3f",
#               rownames(check_matrix)[idx[1]],
#               colnames(check_matrix)[idx[2]],
#               cors[rownames(check_matrix)[idx[1]], colnames(check_matrix)[idx[2]]])
#     })
#     full_msg <- paste(c(msg, viol_str), collapse = "\n")
#     stop(full_msg)
#   }
#   
#   #  print(paste0("vars is length: ", length(vars)))
#   #  print(vars)
#   return(vars)
# }
# 



## screen to avoid high-correlation amongst risk variables
screen_highcor_random <- function(Y, X, family, obsWeights, id, nVar = ncol(X), corr_threshold = 0.9, ...) {
  
  # ---- 1. Prioritization Function ----
  get_priorities <- if (study_name != "HVTN705") {
    function(x){
      if((str_detect(x[1], "mn50") & str_detect(x[2], "pseudo")) | 
         (str_detect(x[2], "mn50") & str_detect(x[1], "pseudo"))){
        selectVar <- if_else(str_detect(x[1], "mn50"), x[1], x[2]) 
      } else if((str_detect(x[1], "pseudo") & str_detect(x[2], "bind")) | 
                (str_detect(x[2], "pseudo") & str_detect(x[1], "bind"))){
        selectVar <- if_else(str_detect(x[1], "pseudo"), x[1], x[2])
      } else if((str_detect(x[1], "Spike") & str_detect(x[2], "RBD")) | 
                (str_detect(x[2], "Spike") & str_detect(x[1], "RBD"))){
        selectVar <- if_else(str_detect(x[1], "RBD"), x[1], x[2])
      } else if((str_detect(x[1], "id50") & str_detect(x[2], "id80")) | 
                (str_detect(x[2], "id50") & str_detect(x[1], "id80"))){
        selectVar <- if_else(str_detect(x[1], "id50"), x[1], x[2])
      } else {
        selectVar <- sample(c(x[1], x[2]), 1)
      }
      return(selectVar)
    }
  } else {
    function(x) {
      selectVar <- sample(c(x[1], x[2]), 1)
      return(selectVar)
    }
  }
  
  all_vars <- colnames(X)
  cors <- cor(X, method = "spearman", use = "pairwise.complete.obs")
  diag(cors) <- NA
  
  # ---- 2. Identify High Correlation Pairs ----
  high_corr_pairs <- as.data.frame(as.table(cors), stringsAsFactors = FALSE) %>%
    mutate(var1 = as.character(Var1),
           var2 = as.character(Var2),
           corr = Freq) %>%
    filter(!is.na(corr), var1 != var2, abs(corr) > corr_threshold) %>%
    rowwise() %>%
    mutate(pair = paste(sort(c(var1, var2)), collapse = "_")) %>%
    ungroup() %>%
    distinct(pair, .keep_all = TRUE) %>%
    select(var1, var2)
  
  # ---- 3. Selection Logic ----
  set.seed(123)
  selected <- briskfactors
  dropped <- character()
  
  repeat {
    current_pairs <- high_corr_pairs %>%
      filter(!(var1 %in% dropped | var2 %in% dropped))
    
    if (nrow(current_pairs) == 0) break
    
    for (i in seq_len(nrow(current_pairs))) {
      a <- current_pairs$var1[i]
      b <- current_pairs$var2[i]
      
      if (a %in% dropped | b %in% dropped) next
      
      if (a %in% selected & !(b %in% selected)) {
        dropped <- c(dropped, b)
      } else if (b %in% selected & !(a %in% selected)) {
        dropped <- c(dropped, a)
      } else if (a %in% selected & b %in% selected) {
        keep <- get_priorities(c(a, b))
        drop <- setdiff(c(a, b), keep)
        selected <- setdiff(selected, drop)
        dropped <- c(dropped, drop)
      } else {
        keep <- get_priorities(c(a, b))
        drop <- setdiff(c(a, b), keep)
        selected <- c(selected, keep)
        dropped <- c(dropped, drop)
      }
    }
    
    # Check for any remaining violations
    if (length(selected) > 1) {
      check_matrix <- abs(cors[selected, selected])
      check_matrix[lower.tri(check_matrix, diag = TRUE)] <- NA
      violations <- which(check_matrix > corr_threshold, arr.ind = TRUE)
      if (nrow(violations) == 0) break
    } else {
      break
    }
  }
  
  
  #############################################
  # ---- 4. Add Uninvolved Vars ----
  involved_vars <- unique(c(high_corr_pairs$var1, high_corr_pairs$var2))
  uninvolved_vars <- setdiff(all_vars, union(involved_vars, dropped))
  selected <- unique(c(selected, uninvolved_vars))
  
  # ---- 6. Final Check (safety) ----
  if (length(selected) > 1) {
    check_matrix <- abs(cors[selected, selected])
    check_matrix[lower.tri(check_matrix, diag = TRUE)] <- NA
    violations <- which(check_matrix > corr_threshold, arr.ind = TRUE)
    if (nrow(violations) > 0) {
      msg <- "Correlation violations found among selected variables:\n"
      viol_str <- apply(violations, 1, function(idx) {
        sprintf("  %s ~ %s = %.3f",
                rownames(check_matrix)[idx[1]],
                colnames(check_matrix)[idx[2]],
                cors[rownames(check_matrix)[idx[1]], colnames(check_matrix)[idx[2]]])
      })
      stop(paste(c(msg, viol_str), collapse = "\n"))
    }
  }
  
  # ---- 7. Return Named Logical Vector ----
  vars <- setNames(all_vars %in% selected, all_vars)
  return(vars)
}


## -------------------------------------------------------------------------------------
## SL algorithms
## -------------------------------------------------------------------------------------

## --------------------------------------------------------------------------
## define wrappers that are less memory-intensive than the usual SL functions
## --------------------------------------------------------------------------
# skinny glm
SL.glm.skinny <- function(Y, X, newX, family, obsWeights, ...) {
  SL.glm.fit <- SL.glm(Y = Y, X = X, newX = newX, family = family, obsWeights = obsWeights, ...)
  SL.glm.fit$fit$object$y <- NULL
  SL.glm.fit$fit$object$model <- NULL
  SL.glm.fit$fit$object$residuals <- NULL
  SL.glm.fit$fit$object$fitted.values <- NULL
  SL.glm.fit$fit$object$effects <- NULL
  SL.glm.fit$fit$object$qr$qr <- NULL
  SL.glm.fit$fit$object$linear.predictors <- NULL
  SL.glm.fit$fit$object$weights <- NULL
  SL.glm.fit$fit$object$prior.weights <- NULL
  SL.glm.fit$fit$object$data <- NULL
  SL.glm.fit$fit$object$family$variance <- NULL
  SL.glm.fit$fit$object$family$dev.resids <- NULL
  SL.glm.fit$fit$object$family$aic <- NULL
  SL.glm.fit$fit$object$family$validmu <- NULL
  SL.glm.fit$fit$object$family$simulate <- NULL
  attr(SL.glm.fit$fit$object$terms, ".Environment") <- NULL
  attr(SL.glm.fit$fit$object$formula, ".Environment") <- NULL
  return(SL.glm.fit)
}

## skinny glm with interactions
SL.glm.interaction.skinny <- function(Y, X, newX, family, obsWeights, ...) {
  SL.glm.fit <- SL.glm.interaction(Y = Y, X = X, newX = newX, family = family, obsWeights = obsWeights, ...)
  SL.glm.fit$fit$object$y <- NULL
  SL.glm.fit$fit$object$model <- NULL
  SL.glm.fit$fit$object$residuals <- NULL
  SL.glm.fit$fit$object$fitted.values <- NULL
  SL.glm.fit$fit$object$effects <- NULL
  SL.glm.fit$fit$object$qr$qr <- NULL
  SL.glm.fit$fit$object$linear.predictors <- NULL
  SL.glm.fit$fit$object$weights <- NULL
  SL.glm.fit$fit$object$prior.weights <- NULL
  SL.glm.fit$fit$object$data <- NULL
  SL.glm.fit$fit$object$family$variance <- NULL
  SL.glm.fit$fit$object$family$dev.resids <- NULL
  SL.glm.fit$fit$object$family$aic <- NULL
  SL.glm.fit$fit$object$family$validmu <- NULL
  SL.glm.fit$fit$object$family$simulate <- NULL
  attr(SL.glm.fit$fit$object$terms, ".Environment") <- NULL
  attr(SL.glm.fit$fit$object$formula, ".Environment") <- NULL
  return(SL.glm.fit)
}

# skinny stepwise with interactions
SL.step.interaction.skinny <- function(Y, X, newX, family, obsWeights, ...) {
  SL.step.interaction.fit <- SL.step.interaction(
    Y = Y, X = X, newX = newX, family = family,
    obsWeights = obsWeights, direction = "forward", ...
  )
  SL.step.interaction.fit$fit$object$y <- NULL
  SL.step.interaction.fit$fit$object$model <- NULL
  SL.step.interaction.fit$fit$object$residuals <- NULL
  SL.step.interaction.fit$fit$object$fitted.values <- NULL
  SL.step.interaction.fit$fit$object$effects <- NULL
  SL.step.interaction.fit$fit$object$qr$qr <- NULL
  SL.step.interaction.fit$fit$object$linear.predictors <- NULL
  SL.step.interaction.fit$fit$object$weights <- NULL
  SL.step.interaction.fit$fit$object$prior.weights <- NULL
  SL.step.interaction.fit$fit$object$data <- NULL
  SL.step.interaction.fit$fit$object$family$variance <- NULL
  SL.step.interaction.fit$fit$object$family$dev.resids <- NULL
  SL.step.interaction.fit$fit$object$family$aic <- NULL
  SL.step.interaction.fit$fit$object$family$validmu <- NULL
  SL.step.interaction.fit$fit$object$family$simulate <- NULL
  attr(SL.step.interaction.fit$fit$object$terms, ".Environment") <- NULL
  attr(SL.step.interaction.fit$fit$object$formula, ".Environment") <- NULL
  return(SL.step.interaction.fit)
}

# skinny stepwise (forward)
SL.step.skinny <- function(Y, X, newX, family, obsWeights, ...) {
  SL.step.fit <- SL.step(Y = Y, X = X, newX = newX, family = family, obsWeights = obsWeights, direction = "forward", ...)
  SL.step.fit$fit$object$y <- NULL
  SL.step.fit$fit$object$model <- NULL
  SL.step.fit$fit$object$residuals <- NULL
  SL.step.fit$fit$object$fitted.values <- NULL
  SL.step.fit$fit$object$effects <- NULL
  SL.step.fit$fit$object$qr$qr <- NULL
  SL.step.fit$fit$object$linear.predictors <- NULL
  SL.step.fit$fit$object$weights <- NULL
  SL.step.fit$fit$object$prior.weights <- NULL
  SL.step.fit$fit$object$data <- NULL
  SL.step.fit$fit$object$family$variance <- NULL
  SL.step.fit$fit$object$family$dev.resids <- NULL
  SL.step.fit$fit$object$family$aic <- NULL
  SL.step.fit$fit$object$family$validmu <- NULL
  SL.step.fit$fit$object$family$simulate <- NULL
  attr(SL.step.fit$fit$object$terms, ".Environment") <- NULL
  attr(SL.step.fit$fit$object$formula, ".Environment") <- NULL
  return(SL.step.fit)
}

#############################################################################################
#############################################################################################
# boosted decision trees
SL.xgboost.2.no <- function(Y, X, newX, family, obsWeights, max_depth = 2, shrinkage = 0.1, scale_pos_weight = 1, ...) {
  SL.xgboost2(Y = Y, X = X, newX = newX, family = family, obsWeights = obsWeights, max_depth = max_depth, shrinkage = shrinkage,
              scale_pos_weight = scale_pos_weight, ...)
}
SL.xgboost.4.no <- function(Y, X, newX, family, obsWeights, max_depth = 4, shrinkage = 0.1, scale_pos_weight = 1, ...) {
  SL.xgboost2(Y = Y, X = X, newX = newX, family = family, obsWeights = obsWeights, max_depth = max_depth, shrinkage = shrinkage, 
              scale_pos_weight = scale_pos_weight, ...)
}
SL.xgboost.2.yes <- function(Y, X, newX, family, obsWeights, max_depth = 2, shrinkage = 0.1, scale_pos_weight = sum(1 - Y) / sum(Y), ...) {
  SL.xgboost2(Y = Y, X = X, newX = newX, family = family, obsWeights = obsWeights, max_depth = max_depth, shrinkage = shrinkage,
              scale_pos_weight = scale_pos_weight, ...)
}
SL.xgboost.4.yes <- function(Y, X, newX, family, obsWeights, max_depth = 4, shrinkage = 0.1, scale_pos_weight = sum(1 - Y) / sum(Y), ...) {
  SL.xgboost2(Y = Y, X = X, newX = newX, family = family, obsWeights = obsWeights, max_depth = max_depth, shrinkage = shrinkage, 
              scale_pos_weight = scale_pos_weight, ...)
}

# boosted trees
SL.xgboost2 <- function (Y, X, newX, family, obsWeights, id, ntrees = 1000,
                         max_depth = 4, shrinkage = 0.1, minobspernode = 10, scale_pos_weight = 1, params = list(),
                         nthread = 1, verbose = 0, save_period = NULL, ...){
  if (!is.matrix(X)) {
    X = model.matrix(~. - 1, X)
  }
  xgmat = xgboost::xgb.DMatrix(data = X, label = Y, weight = obsWeights)
  
  if (family$family == "gaussian") {
    model = xgboost::xgboost(data = xgmat, objective = "reg:squarederror", ## if xgboost version >=1.1.1.1, changed from reg:linear to reg:squarederror
                             nrounds = ntrees, max_depth = max_depth, min_child_weight = minobspernode, scale_pos_weight = scale_pos_weight,
                             eta = shrinkage, verbose = verbose, nthread = nthread,
                             params = params, save_period = save_period)
  }
  
  if (family$family == "binomial") {
    model = xgboost::xgboost(data = xgmat, objective = "binary:logistic",
                             nrounds = ntrees, max_depth = max_depth, min_child_weight = minobspernode, scale_pos_weight = scale_pos_weight,
                             eta = shrinkage, verbose = verbose, nthread = nthread,
                             params = params, save_period = save_period)
  }
  
  if (family$family == "multinomial") {
    model = xgboost::xgboost(data = xgmat, objective = "multi:softmax",
                             nrounds = ntrees, max_depth = max_depth, min_child_weight = minobspernode, scale_pos_weight = scale_pos_weight,
                             eta = shrinkage, verbose = verbose, num_class = length(unique(Y)),
                             nthread = nthread, params = params, save_period = save_period)
  }
  
  if (!is.matrix(newX)) {
    newX = model.matrix(~. - 1, newX)
  }
  
  pred = predict(model, newdata = newX)
  fit = list(object = model)
  class(fit) = c("SL.xgboost")
  out = list(pred = pred, fit = fit)
  
  return(out)
}

# random forests w/out class balancing
SL.ranger.no <- function(Y, X, newX, family, obsWeights = rep(1, length(Y)), case.weights = NULL, ...) {
  SL.ranger.imp(Y = Y, X = X, newX = newX, family = family, obsWeights = obsWeights, case.weights = case.weights, ...)
}
SL.ranger.yes <- function(Y, X, newX, family, obsWeights = rep(1, length(Y)), case.weights = Y * sum(Y) ^ (-1) + (1 - Y) * sum(1 - Y) ^ (-1), ...) {
  SL.ranger.imp(Y = Y, X = X, newX = newX, family = family, obsWeights = obsWeights, case.weights = case.weights, ...)
}


# neural nets with size 2 and 5
SL.nnet.2 <- function(Y, X, newX, family, obsWeights = rep(1, length(Y)), size = 2, ...) {
  SL.nnet(Y = Y, X = X, newX = newX, family = family, obsWeights = obsWeights, size = size, ...)
}
SL.nnet.5 <- function(Y, X, newX, family, obsWeights = rep(1, length(Y)), size = 5, ...) {
  SL.nnet(Y = Y, X = X, newX = newX, family = family, obsWeights = obsWeights, size = size, ...)
}



# ksvm with kernel = "rbfdot" or "polydot"
SL.ksvm.rbfdot <- function(Y, X, newX, family, obsWeights = rep(1, length(Y)), kernel = "rbfdot", ...) {
  SL.ksvm(Y = Y, X = X, newX = newX, family = family, obsWeights = obsWeights, kernel = kernel, ...)
}
SL.ksvm.polydot <- function(Y, X, newX, family, obsWeights = rep(1, length(Y)), kernel = "polydot", ...) {
  SL.ksvm(Y = Y, X = X, newX = newX, family = family, obsWeights = obsWeights, kernel = kernel, ...)
}


# Lasso with alpha 0 or 1
SL.glmnet.0 <- function(Y, X, newX, family, obsWeights = rep(1, length(Y)), alpha = 0, ...) {
  SL.glmnet(Y = Y, X = X, newX = newX, family = family, obsWeights = obsWeights, alpha = alpha, ...)
}
SL.glmnet.0.33 <- function(Y, X, newX, family, obsWeights = rep(1, length(Y)), alpha = 0.33, ...) {
  SL.glmnet(Y = Y, X = X, newX = newX, family = family, obsWeights = obsWeights, alpha = alpha, ...)
}
SL.glmnet.0.67 <- function(Y, X, newX, family, obsWeights = rep(1, length(Y)), alpha = 0.67, ...) {
  SL.glmnet(Y = Y, X = X, newX = newX, family = family, obsWeights = obsWeights, alpha = alpha, ...)
}
SL.glmnet.1 <- function(Y, X, newX, family, obsWeights = rep(1, length(Y)), alpha = 1, ...) {
  SL.glmnet(Y = Y, X = X, newX = newX, family = family, obsWeights = obsWeights, alpha = alpha, ...)
}


# naive bayes wrapper
SL.naivebayes <- function(Y, X, newX, family, obsWeights, laplace = 0, ...){
  SuperLearner:::.SL.require("e1071")
  if(family$family == "gaussian"){
    stop("SL.naivebayes only works with binary outcomes")
  }else{
    nb <- naiveBayes(y = Y, x = X, laplace = laplace)
    pred <- predict(nb, newX, type = "raw")[,2]
    out <- list(fit = list(object = nb), pred = pred)
    class(out$fit) <- "SL.naivebayes"
    return(out)
  }
}

# predict method for naive bayes wrapper
predict.SL.naivebayes <- function(object, newdata, ...){
  pred <- predict(object$object, newdata = newdata, type = "raw")[,2]
  return(pred)
}

#############################################################################################
#############################################################################################

# random forests
SL.ranger.imp <- function (Y, X, newX, family, obsWeights, case.weights,
                           num.trees = 500, mtry = floor(sqrt(ncol(X))),
                           write.forest = TRUE, probability = family$family == "binomial",
                           min.node.size = ifelse(family$family == "gaussian", 5, 1),
                           replace = TRUE, sample.fraction = ifelse(replace, 1, 0.632),
                           num.threads = 1, verbose = TRUE, ...) {
  SuperLearner:::.SL.require("ranger")
  if (family$family == "binomial") {
    Y = as.factor(Y)
  }
  if (is.matrix(X)) {
    X = data.frame(X)
  }
  fit <- ranger::ranger(`_Y` ~ ., data = cbind(`_Y` = Y, X),
                        num.trees = num.trees, mtry = mtry, min.node.size = min.node.size,
                        replace = replace, sample.fraction = sample.fraction,
                        case.weights = case.weights, write.forest = write.forest,
                        probability = probability, num.threads = num.threads,
                        verbose = verbose, importance = "impurity")
  pred <- predict(fit, data = newX)$predictions
  if (family$family == "binomial") {
    pred = pred[, "1"]
  }
  fit <- list(object = fit, verbose = verbose)
  class(fit) <- c("SL.ranger")
  out <- list(pred = pred, fit = fit)
  return(out)
}


## ----------------------------------------------------------------------------
## Create SL Library
## ----------------------------------------------------------------------------
if (run_prod & !Sys.getenv("TRIAL") %in% c("covail_xassays")) {
  # NOTE: fancier library for production run
  # learners in the method1 are also combined with no screen
  if(Sys.getenv("TRIAL") %in% c("hvtn705second", "moderna_real")){
    methods1 <- c("SL.mean", "SL.glm", 
                  "SL.glmnet.0", "SL.glmnet.1", 
                  "SL.xgboost.2.no", "SL.xgboost.4.no", "SL.xgboost.2.yes", "SL.xgboost.4.yes",
                  "SL.ranger.no", "SL.ranger.yes") 
  } else if(Sys.getenv("TRIAL") %in% c("janssen_pooled_partA", "janssen_la_partA")){
    methods1 <- c("SL.mean", "SL.glm", "SL.bayesglm", "SL.glm.interaction",
                  "SL.glmnet.1", #"SL.glmnet.0", 
                  #"SL.gam", "SL.ksvm.rbfdot", "SL.ksvm.polydot", "SL.polymars",
                  "SL.xgboost.4.no", #"SL.xgboost.2.no", "SL.xgboost.2.yes", "SL.xgboost.4.yes",
                  "SL.ranger.no") #, "SL.ranger.yes"  
  } else if (Sys.getenv("TRIAL") %in% c("covail_tcell", "covail_xassays")){
    methods1 <- c("SL.mean"#, 
                  #"SL.glm" #,  
                  # "SL.bayesglm", 
                  # "SL.glmnet.0", "SL.glmnet.0.33", "SL.glmnet.0.67", "SL.glmnet.1", 
                  # #"SL.xgboost.2.no", 
                  # "SL.xgboost.4.no", 
                  # #"SL.xgboost.2.yes", "SL.xgboost.4.yes",
                  # "SL.ranger.no", "SL.ranger.yes"
                  ) 
  }

  # learners in the method2 are learners that can have screens
  if(Sys.getenv("TRIAL") %in% c("moderna_real")){
    methods2 <- c("SL.glm")
  } else if(Sys.getenv("TRIAL") == "hvtn705second"){
    methods2 <- c("SL.glm", "SL.gam") #, "SL.bayesglm", "SL.glm.interaction", "SL.gam", "SL.ksvm.rbfdot", "SL.ksvm.polydot", "SL.polymars"
  } else if(Sys.getenv("TRIAL") %in% c("janssen_pooled_partA", "janssen_la_partA")){
    methods2 <- c("SL.glm", "SL.bayesglm", "SL.glm.interaction", "SL.gam",
                  "SL.ksvm.rbfdot", "SL.ksvm.polydot", "SL.polymars") 
  } else if(Sys.getenv("TRIAL") %in% c("covail_tcell", "covail_xassays")){
    methods2 <- c("SL.glm",
                  "SL.bayesglm",
                  "SL.glm.interaction", "SL.gam",
                  #"SL.ksvm.rbfdot", "SL.ksvm.polydot", # never included as model fitting code breaks
                  "SL.nnet.2", "SL.nnet.5"
                  #"SL.polymars" # never included as model fitting code breaks
                  ) 
  }
} else {
  # NOTE: smaller library for ~faster~ demo run
  # learners in the method1 are also combined with no screen
  methods1 <- c("SL.mean", 
                "SL.glm")

  # learners in the method2 are learners that can have screens
  methods2 <- c("SL.glm")
}

screens1 <- "screen_all"
screens2 <- c(
  "screen_glmnet",
  "screen_univariate_logistic_pval",
  "screen_highcor_random"
)

SL_library1 <- sapply(
  1:length(methods1),
  function(i) c(methods1[i], screens1),
  simplify = FALSE
)

SL_library2 <- sapply(
  1:length(methods2),
  function(i) c(methods2[i], screens2),
  simplify = FALSE
)

SL_library <- c(SL_library1, SL_library2)

#############################################################################

if(run_prod & Sys.getenv("TRIAL") %in% c("covail_xassays")){
  methods1 <- c("SL.mean")
  
  methods2 <- c(
    "SL.glm",
    "SL.bayesglm",
    "SL.glm.interaction", 
    "SL.glmnet.0", "SL.glmnet.0.33", "SL.glmnet.0.67", "SL.glmnet.1", 
    "SL.gam",
    "SL.nnet.2", "SL.nnet.5",
    "SL.xgboost.4.no", "SL.xgboost.4.yes",
    "SL.ranger.yes", "SL.ranger.no"
  )
  
  # Learners that only use screen_highcor_random
  only_highcor_screen <- c(
    "SL.glmnet.0", "SL.glmnet.0.33", "SL.glmnet.0.67", "SL.glmnet.1", "SL.xgboost.4.no", "SL.xgboost.4.yes"
  )
  
  # Learners that use only glmnet + univariate p-value screens
  only_glmnet_univariate <- c("SL.glm.interaction", "SL.gam")
  
  # SL.mean with screen_all only
  SL_library1 <- lapply(methods1, function(m) c(m, "screen_all"))
  
  # All others with conditionally assigned screeners
  SL_library2 <- lapply(methods2, function(m) {
    if (m %in% only_highcor_screen) {
      c(m, "screen_highcor_random")
    } else if (m %in% only_glmnet_univariate) {
      c(m, "screen_glmnet", "screen_univariate_logistic_pval")
    } else {
      c(m, "screen_glmnet", "screen_univariate_logistic_pval", "screen_highcor_random")
    }
  })
  
  SL_library <- c(SL_library1, SL_library2)
  # SL_library <- SL_library[c(1, 2, 3, 9, 10, 11)]
  # 
  # # Remove "screen_glmnet" from all elements
  # SL_library <- lapply(SL_library, function(x) setdiff(x, "screen_glmnet"))
}

