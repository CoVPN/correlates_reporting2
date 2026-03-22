require(caret)
require(corpcor)

# Taken from Antonio Olivas's repository here: 
# https://github.com/aolivasm/proxi-mtp/

#######################################################
########## NON-PARAMETRIC ESTIMATION ##################

### Gaussian kernel
rbfkernel <- function(X, sigma = 1, Y = NULL){
  ## sigma is sigma^2
  
  # test if X is a matrix
  if(!is.matrix(X))
  {
    print("X must be a matrix containing samples in its rows")
    return()
  }
  # test if sigma is a number and > 0
  if(length(sigma) != 1 || sigma <= 0)
  {
    print("sigma must be a number > 0 specifying the rbf-kernel width")
    return()
  }
  if(!is.null(Y))
  {
    # test if Y is a matrix
    if(!is.matrix(Y))
    {
      print("Y must be a matrix containing samples in its rows or NULL if it should not be used")
      return()
    }
    # test if vectors in X and Y have same dimension
    if(ncol(X) != ncol(Y))
    {
      print("The samples in the rows of X and Y must be of same dimension")
      return()
    }
  }
  
  n <- nrow(X) # number of samples in X
  
  if(is.null(Y))
  {
    # calculate distance matrix
    XtX <- tcrossprod(X)
    XX <- matrix(1, n) %*% diag(XtX)
    D <- XX - 2*XtX + t(XX) # distance matrix
  }
  else
  {
    m <- nrow(Y) # number of samples in Y
    # calculate distance matrix (between vectors of X and Y)
    XX <- matrix(apply(X ^ 2, 1, sum), n, m)
    YY <- matrix(apply(Y ^ 2, 1, sum), n, m, byrow = TRUE)
    XY <- tcrossprod(X, Y)
    D <- XX - 2*XY + YY
  }
  
  # calculate rbf-kernel matrix
  K <- exp(-D/(2*sigma))
  
  return(K)
}

### Function to computed initial bandwidth using the weighted Euclidean distances
wt_euclidean_bw <- function(X, w) {
  n <- nrow(X)
  dists <- numeric(n * (n - 1) / 2)
  weights <- numeric(n * (n - 1) / 2)
  
  idx <- 1
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      dists[idx] <- sum((X[i, ] - X[j, ])^2)
      weights[idx] <- w[i] * w[j]
      idx <- idx + 1
    }
  }
  
  # Normalize weights
  weights <- weights / sum(weights)
  
  # Compute weighted median
  ord <- order(dists)
  dists_sorted <- dists[ord]
  weights_sorted <- weights[ord]
  cum_weights <- cumsum(weights_sorted)
  
  # Find the first index where the cumulative weight is >= 0.5
  median_idx <- which(cum_weights >= 0.5)[1]
  bw <- dists_sorted[median_idx] / 2
  return(bw)
}

inv.mat <- function(matrix, eps = 1e-10, tol = 1e-12) {
  cond <- rcond(matrix)
  if (is.nan(cond) || cond < tol) {
    #message("Matrix is ill-conditioned. Adding regularization.")
    matrix <- matrix + diag(rep(eps, nrow(matrix)))
  }
  return(solve(matrix))
}

### Interior matrices
### Note: for Gamma_H provide K %*% qS instead of K
Gamm = function(K, lm, n, ident_n){
  0.25 * inv.mat(K + lm * log(n) * ident_n)
}

Gamm_wt = function(K, lm, n, m, ident_m){
  0.25 * inv.mat(K/n + lm * log(m)/m * ident_m)
}

### Bridge function computed using the closed-form expression using RKHSs
#### For outcome bridge function, A = K_G %*% Y
#### For treatment bridge function, Kint = KH %*% qS, Gamma = qS %*% Gamma_G_H,
####   and A = t(\tilde{KH}) %*% 1_n
bridge_fn_rkhs = function(Kext, Kint, ident_n, A, Gamm, lm, n){
  inv.mat(Gamm %*% Kint %*% Kext + lm* sqrt(log(n)*n) * ident_n) %*%
    Gamm %*% A
}

bridge_fn_rkhs_wt = function(Kext, Kint, ident_m, A, Gamm, lm, n, m){
  inv.mat(Gamm %*% Kint %*% Kext + n^2 * lm* sqrt(log(m)/m) * ident_m) %*%
    Gamm %*% A
}

cv_h = function(Theta, n, ident_n, n_te, arg_ext, arg_int, arg_ext_te,
                Y, Y_te,
                lm_int_list, lm_ext_list, bw0_int, bw0_ext,
                bw_ext_scale_list, bw_int_scale_list){
  
  fn_ext_norm_array = array(0, dim = c(length(bw_int_scale_list), length(bw_ext_scale_list), length(lm_ext_list), length(lm_int_list)))
  fn_int_norm_array = array(0, dim = c(length(bw_int_scale_list), length(bw_ext_scale_list), length(lm_ext_list), length(lm_int_list)))
  
  risk_array = array(0, dim = c(length(bw_int_scale_list), length(bw_ext_scale_list), length(lm_ext_list), length(lm_int_list)))
  alpha_array = array(0, dim = c(length(bw_int_scale_list), length(bw_ext_scale_list),length(lm_ext_list), length(lm_int_list), n))
  
  for(bw_int_ind in 1:length(bw_int_scale_list)){
    ### Kernel basis for internal function class at training data
    Kint = rbfkernel(arg_int, sigma = bw0_int*bw_int_scale_list[bw_int_ind])
    
    for (lm_int_index in 1:length(lm_int_list)){
      lm_int = lm_int_list[lm_int_index]
      
      Gamm_int = Gamm(K = Kint, lm = lm_int, n = n, ident_n = ident_n)
      
      for (lm_ext_index in 1:length(lm_ext_list)){
        lm_ext = lm_ext_list[lm_ext_index]
        
        for(bw_ext_ind in 1:length(bw_ext_scale_list)){
          Kext = rbfkernel(arg_ext, sigma = bw0_ext*bw_ext_scale_list[bw_ext_ind])
          Kext_te = rbfkernel(X = arg_ext_te, sigma = bw0_ext*bw_ext_scale_list[bw_ext_ind], Y = arg_ext)
          
          alpha =  bridge_fn_rkhs(Kext = Kext, Kint = Kint, ident_n = ident_n, A = Kint %*% Y, Gamm = Gamm_int, lm = lm_ext, n = n)
          
          ### Compute H_norm
          norm_ext = sqrt(as.numeric(t(alpha) %*% Kext  %*% alpha))
          
          xi0 = Y- Kext %*% alpha
          norm_int = sqrt(t(xi0) %*% Gamm_int %*% Kint %*% Gamm_int %*% xi0)/2
          
          xi_n = (Y_te- Kext_te %*% alpha) / n_te
          risk_array[bw_int_ind, bw_ext_ind, lm_ext_index, lm_int_index] = t(xi_n) %*% Theta %*% xi_n
          alpha_array[bw_int_ind, bw_ext_ind, lm_ext_index, lm_int_index,] = alpha
          fn_int_norm_array[bw_int_ind, bw_ext_ind, lm_ext_index, lm_int_index] = norm_int
          fn_ext_norm_array[bw_int_ind, bw_ext_ind, lm_ext_index, lm_int_index] = norm_ext
        }
      }
    }
  }  
  ind_bw_int = which(risk_array == min(risk_array), arr.ind = TRUE)[1, 1]
  ind_bw_ext = which(risk_array == min(risk_array), arr.ind = TRUE)[1, 2]
  ind_lm_ext = which(risk_array == min(risk_array), arr.ind = TRUE)[1, 3]
  ind_lm_int = which(risk_array == min(risk_array), arr.ind = TRUE)[1, 4]
  
  lm_ext = lm_ext_list[ind_lm_ext]
  lm_int = lm_int_list[ind_lm_int]
  
  bw_int = bw0_int*bw_int_scale_list[ind_bw_int]
  bw_ext = bw0_ext*bw_ext_scale_list[ind_bw_ext]
  
  ### Recovering h
  alpha = alpha_array[ind_bw_int, ind_bw_ext, ind_lm_ext, ind_lm_int, ]
  
  norm_int = fn_int_norm_array[ind_bw_int, ind_bw_ext, ind_lm_ext, ind_lm_int]
  norm_ext = fn_ext_norm_array[ind_bw_int, ind_bw_ext, ind_lm_ext, ind_lm_int]
  
  res = list(lm_ext, lm_int, bw_int, bw_ext, alpha, norm_int, norm_ext)
  names(res) = c("lm_ext", "lm_int", "bw_int", "bw_ext", "alpha", "norm_int", "norm_ext")
  return(res)
}

cv_g = function(Theta, n, ident_n, n_te, arg_ext, arg_int, arg_int_q, arg_ext_te,
                qS, qS_te, ind_S,
                Kint_te,  A_te,
                lm_int_list, lm_ext_list, bw0_int, bw0_ext,
                bw_ext_scale_list, bw_int_scale_list){
  
  fn_ext_norm_array = array(0, dim = c(length(bw_int_scale_list), length(bw_ext_scale_list), length(lm_ext_list), length(lm_int_list)))
  fn_int_norm_array = array(0, dim = c(length(bw_int_scale_list), length(bw_ext_scale_list), length(lm_ext_list), length(lm_int_list)))
  
  risk_array = array(0, dim = c(length(bw_int_scale_list), length(bw_ext_scale_list), length(lm_ext_list), length(lm_int_list)))
  alpha_array = array(0, dim = c(length(bw_int_scale_list), length(bw_ext_scale_list),length(lm_ext_list), length(lm_int_list), n))
  
  for(bw_int_ind in 1:length(bw_int_scale_list)){
    ### Kernel basis for internal function class at training data
    Kint_tilde = rbfkernel(arg_int_q, sigma = bw0_int*bw_int_scale_list[bw_int_ind], Y = arg_int)
    Kint = rbfkernel(arg_int, sigma = bw0_int*bw_int_scale_list[bw_int_ind])
    
    for (lm_int_index in 1:length(lm_int_list)){
      lm_int = lm_int_list[lm_int_index]
      
      Gamm_int = Gamm(K = Kint %*% qS, lm = lm_int, n = n, ident_n = ident_n)
      
      for (lm_ext_index in 1:length(lm_ext_list)){
        lm_ext = lm_ext_list[lm_ext_index]
        
        for(bw_ext_ind in 1:length(bw_ext_scale_list)){
          Kext = rbfkernel(arg_ext, sigma = bw0_ext*bw_ext_scale_list[bw_ext_ind])
          Kext_te = rbfkernel(X = arg_ext_te, sigma = bw0_ext*bw_ext_scale_list[bw_ext_ind], Y = arg_ext)
          
          A = colSums(ind_S %*% Kint_tilde)
          alpha =  bridge_fn_rkhs(Kext = Kext, Kint = Kint  %*%  qS, ident_n = ident_n, A = A, Gamm = qS %*%  Gamm_int, lm = lm_ext, n = n)
          
          ### Compute H_norm
          norm_ext = sqrt(as.numeric(t(alpha) %*% Kext  %*% alpha))
          
          xi0 =  A- Kint %*% qS %*% Kext %*% alpha
          norm_int = sqrt(t(xi0) %*% t(Gamm_int) %*% inv.mat(Kint) %*% Gamm_int %*% xi0)/2
          
          xi_n = (A_te- Kint_te %*% qS_te %*% Kext_te %*% alpha) / n_te
          risk_array[bw_int_ind, bw_ext_ind, lm_ext_index, lm_int_index] = t(xi_n) %*% Theta %*% xi_n
          alpha_array[bw_int_ind, bw_ext_ind, lm_ext_index, lm_int_index,] = alpha
          fn_int_norm_array[bw_int_ind, bw_ext_ind, lm_ext_index, lm_int_index] = norm_int
          fn_ext_norm_array[bw_int_ind, bw_ext_ind, lm_ext_index, lm_int_index] = norm_ext
        }
      }
    }
  }  
  ind_bw_int = which(risk_array == min(risk_array), arr.ind = TRUE)[1, 1]
  ind_bw_ext = which(risk_array == min(risk_array), arr.ind = TRUE)[1, 2]
  ind_lm_ext = which(risk_array == min(risk_array), arr.ind = TRUE)[1, 3]
  ind_lm_int = which(risk_array == min(risk_array), arr.ind = TRUE)[1, 4]
  
  lm_ext = lm_ext_list[ind_lm_ext]
  lm_int = lm_int_list[ind_lm_int]
  
  bw_int = bw0_int*bw_int_scale_list[ind_bw_int]
  bw_ext = bw0_ext*bw_ext_scale_list[ind_bw_ext]
  
  ### Recovering g
  alpha = alpha_array[ind_bw_int, ind_bw_ext, ind_lm_ext, ind_lm_int, ]
  
  norm_int = fn_int_norm_array[ind_bw_int, ind_bw_ext, ind_lm_ext, ind_lm_int]
  norm_ext = fn_ext_norm_array[ind_bw_int, ind_bw_ext, ind_lm_ext, ind_lm_int]
  
  res = list(lm_ext, lm_int, bw_int, bw_ext, alpha, norm_int, norm_ext)
  names(res) = c("lm_ext", "lm_int", "bw_int", "bw_ext", "alpha", "norm_int", "norm_ext")
  return(res)
}

cv_h_w = function(Theta, n, m, ident_m, n_te, arg_ext, arg_int, arg_ext_te,
                  Y, Y_te, W, W_te,
                  lm_int_list, lm_ext_list, bw0_int, bw0_ext,
                  bw_ext_scale_list, bw_int_scale_list){
  
  fn_ext_norm_array = array(0, dim = c(length(bw_int_scale_list), length(bw_ext_scale_list), length(lm_ext_list), length(lm_int_list)))
  fn_int_norm_array = array(0, dim = c(length(bw_int_scale_list), length(bw_ext_scale_list), length(lm_ext_list), length(lm_int_list)))
  
  risk_array = array(0, dim = c(length(bw_int_scale_list), length(bw_ext_scale_list), length(lm_ext_list), length(lm_int_list)))
  alpha_array = array(0, dim = c(length(bw_int_scale_list), length(bw_ext_scale_list),length(lm_ext_list), length(lm_int_list), m))
  
  for(bw_int_ind in 1:length(bw_int_scale_list)){
    ### Kernel basis for internal function class at training data
    Kint = rbfkernel(arg_int, sigma = bw0_int*bw_int_scale_list[bw_int_ind])
    
    for (lm_int_index in 1:length(lm_int_list)){
      lm_int = lm_int_list[lm_int_index]
      
      Gamm_int = Gamm_wt(K = Kint %*% W, lm = lm_int, n = n, m = m, ident_m = ident_m)
      
      for (lm_ext_index in 1:length(lm_ext_list)){
        lm_ext = lm_ext_list[lm_ext_index]
        
        for(bw_ext_ind in 1:length(bw_ext_scale_list)){
          Kext = rbfkernel(arg_ext, sigma = bw0_ext*bw_ext_scale_list[bw_ext_ind])
          Kext_te = rbfkernel(X = arg_ext_te, sigma = bw0_ext*bw_ext_scale_list[bw_ext_ind], Y = arg_ext)
          
          alpha =  bridge_fn_rkhs_wt(Kext = Kext, Kint = Kint %*% W, ident_m = ident_m, A = Kint %*% W %*% Y, Gamm = W %*% Gamm_int, lm = lm_ext, n = n, m = m)
          
          ### Compute H_norm
          norm_ext = sqrt(as.numeric(t(alpha) %*% Kext  %*% alpha))
          
          xi0 = W %*%(Y- Kext %*% alpha)
          norm_int = sqrt(t(xi0) %*% Gamm_int %*% Kint %*% t(Gamm_int) %*% xi0)/2
          
          xi_n = W_te %*% (Y_te- Kext_te %*% alpha) / n_te
          risk_array[bw_int_ind, bw_ext_ind, lm_ext_index, lm_int_index] = t(xi_n) %*% Theta %*% xi_n
          alpha_array[bw_int_ind, bw_ext_ind, lm_ext_index, lm_int_index,] = alpha
          fn_int_norm_array[bw_int_ind, bw_ext_ind, lm_ext_index, lm_int_index] = norm_int
          fn_ext_norm_array[bw_int_ind, bw_ext_ind, lm_ext_index, lm_int_index] = norm_ext
        }
      }
    }
  }  
  ind_bw_int = which(risk_array == min(risk_array), arr.ind = TRUE)[1, 1]
  ind_bw_ext = which(risk_array == min(risk_array), arr.ind = TRUE)[1, 2]
  ind_lm_ext = which(risk_array == min(risk_array), arr.ind = TRUE)[1, 3]
  ind_lm_int = which(risk_array == min(risk_array), arr.ind = TRUE)[1, 4]
  
  lm_ext = lm_ext_list[ind_lm_ext]
  lm_int = lm_int_list[ind_lm_int]
  
  bw_int = bw0_int*bw_int_scale_list[ind_bw_int]
  bw_ext = bw0_ext*bw_ext_scale_list[ind_bw_ext]
  
  ### Recovering h
  alpha = alpha_array[ind_bw_int, ind_bw_ext, ind_lm_ext, ind_lm_int, ]
  
  norm_int = fn_int_norm_array[ind_bw_int, ind_bw_ext, ind_lm_ext, ind_lm_int]
  norm_ext = fn_ext_norm_array[ind_bw_int, ind_bw_ext, ind_lm_ext, ind_lm_int]
  
  res = list(lm_ext, lm_int, bw_int, bw_ext, alpha, norm_int, norm_ext)
  names(res) = c("lm_ext", "lm_int", "bw_int", "bw_ext", "alpha", "norm_int", "norm_ext")
  return(res)
}

cv_g_w = function(Theta, n, m = m, ident_m, n_te, arg_ext, arg_int, arg_int_q, arg_ext_te,
                  qS, qS_te, ind_S, W = W_tr, W_te = W_te, 
                  Kint_te,  A_te,
                  lm_int_list, lm_ext_list, bw0_int, bw0_ext,
                  bw_ext_scale_list, bw_int_scale_list){
  
  fn_ext_norm_array = array(0, dim = c(length(bw_int_scale_list), length(bw_ext_scale_list), length(lm_ext_list), length(lm_int_list)))
  fn_int_norm_array = array(0, dim = c(length(bw_int_scale_list), length(bw_ext_scale_list), length(lm_ext_list), length(lm_int_list)))
  
  risk_array = array(0, dim = c(length(bw_int_scale_list), length(bw_ext_scale_list), length(lm_ext_list), length(lm_int_list)))
  alpha_array = array(0, dim = c(length(bw_int_scale_list), length(bw_ext_scale_list),length(lm_ext_list), length(lm_int_list), m))
  
  for(bw_int_ind in 1:length(bw_int_scale_list)){
    ### Kernel basis for internal function class at training data
    Kint_tilde = rbfkernel(arg_int_q, sigma = bw0_int*bw_int_scale_list[bw_int_ind], Y = arg_int)
    Kint = rbfkernel(arg_int, sigma = bw0_int*bw_int_scale_list[bw_int_ind])
    
    for (lm_int_index in 1:length(lm_int_list)){
      lm_int = lm_int_list[lm_int_index]
      
      Gamm_int = Gamm_wt(K = Kint %*% W %*% qS, lm = lm_int, n = n, m = m, ident_m = ident_m)
      
      for (lm_ext_index in 1:length(lm_ext_list)){
        lm_ext = lm_ext_list[lm_ext_index]
        
        for(bw_ext_ind in 1:length(bw_ext_scale_list)){
          Kext = rbfkernel(arg_ext, sigma = bw0_ext*bw_ext_scale_list[bw_ext_ind])
          Kext_te = rbfkernel(X = arg_ext_te, sigma = bw0_ext*bw_ext_scale_list[bw_ext_ind], Y = arg_ext)
          
          A = colSums(W %*% ind_S %*% Kint_tilde)
          alpha =  bridge_fn_rkhs_wt(Kext = Kext, Kint = Kint  %*% W %*%  qS, ident_m = ident_m, A = A, Gamm = qS %*% W %*%  Gamm_int, lm = lm_ext, n = n, m = m)
          
          ### Compute H_norm
          norm_ext = sqrt(as.numeric(t(alpha) %*% Kext  %*% alpha))
          
          xi0 =  A- Kint %*% W %*% qS %*% Kext %*% alpha
          norm_int = sqrt(t(xi0) %*% t(Gamm_int) %*% inv.mat(Kint) %*% Gamm_int %*% xi0)/2
          
          xi_n = (A_te- Kint_te %*% W_te %*%  qS_te %*% Kext_te %*% alpha) / n_te
          risk_array[bw_int_ind, bw_ext_ind, lm_ext_index, lm_int_index] = t(xi_n) %*% Theta %*% xi_n
          alpha_array[bw_int_ind, bw_ext_ind, lm_ext_index, lm_int_index,] = alpha
          fn_int_norm_array[bw_int_ind, bw_ext_ind, lm_ext_index, lm_int_index] = norm_int
          fn_ext_norm_array[bw_int_ind, bw_ext_ind, lm_ext_index, lm_int_index] = norm_ext
        }
      }
    }
  }  
  ind_bw_int = which(risk_array == min(risk_array), arr.ind = TRUE)[1, 1]
  ind_bw_ext = which(risk_array == min(risk_array), arr.ind = TRUE)[1, 2]
  ind_lm_ext = which(risk_array == min(risk_array), arr.ind = TRUE)[1, 3]
  ind_lm_int = which(risk_array == min(risk_array), arr.ind = TRUE)[1, 4]
  
  lm_ext = lm_ext_list[ind_lm_ext]
  lm_int = lm_int_list[ind_lm_int]
  
  bw_int = bw0_int*bw_int_scale_list[ind_bw_int]
  bw_ext = bw0_ext*bw_ext_scale_list[ind_bw_ext]
  
  ### Recovering g
  alpha = alpha_array[ind_bw_int, ind_bw_ext, ind_lm_ext, ind_lm_int, ]
  
  norm_int = fn_int_norm_array[ind_bw_int, ind_bw_ext, ind_lm_ext, ind_lm_int]
  norm_ext = fn_ext_norm_array[ind_bw_int, ind_bw_ext, ind_lm_ext, ind_lm_int]
  
  res = list(lm_ext, lm_int, bw_int, bw_ext, alpha, norm_int, norm_ext)
  names(res) = c("lm_ext", "lm_int", "bw_int", "bw_ext", "alpha", "norm_int", "norm_ext")
  return(res)
}


expand_names <- function(var, prefix) {
  if (is.null(dim(var))) {
    prefix
  } else {
    paste0(prefix, seq_len(ncol(var)))
  }
}


pmtp = function(data, trt = "X", outcome = "Y", covariates = c("L"),
                nct = "Z", nco = "W", weights = NULL,
                policy, ind_S = NULL,
                ind_qS = NULL,
                theta = 0, K_folds = 3,
                lm_H_list = 10^seq(-5, -1, 1),
                lm_Gh_list = 10^seq(-1, 2, 1),
                lm_G_list = 10^seq(-5, -1, 1),
                lm_Hg_list = 10^seq(-1, 2, 1),
                bw_int_scale_list = 1/4,
                bw_ext_scale_list = c(1/4, 1/2, 1, 2, 4),
                bw_int_fixed_scale = 1/4,
                bw0_H = NULL, bw0_G = NULL,
                inf.fun = FALSE, show.prog = FALSE,
                control.folds = FALSE,
                compute.TMLE.est = FALSE
){
  set.seed(1234)
  required_vars = c(trt, outcome, covariates, nct, nco)
  missing_vars = setdiff(required_vars, colnames(data))
  if (length(missing_vars) > 0) {
    stop(paste("Missing variables in data:", paste(missing_vars, collapse = ", ")))
  }
  
  # Extract variables
  X_val = data[[trt]]
  Y_val = data[[outcome]]
  L_val = as.matrix(data[, covariates, drop = FALSE])
  Z_val = as.matrix(data[, nct, drop = FALSE])
  W_val = as.matrix(data[, nco, drop = FALSE])
  
  if(is.null(weights)){
    wt_val = rep(1, nrow(data))
    no_weights = TRUE
  }else{
    missing_wt = setdiff(weights, colnames(data))
    if (length(missing_wt) > 0) {
      stop("Missing weight variable in data.")
    }
    wt_val = data[[weights]]
    no_weights = all(wt_val == 1)
  }
  
  relevant_data = data.frame(L_val, Z_val, W_val, X_val, Y_val, wt_val)
  complete_rows = complete.cases(relevant_data)
  
  if (!all(complete_rows)) {
    message(sprintf("Removing %d incomplete observations.", sum(!complete_rows)))
    
    # Subset all extracted variables to complete cases
    L_val = L_val[complete_rows, , drop = FALSE]
    Z_val = Z_val[complete_rows, , drop = FALSE]
    W_val = W_val[complete_rows, , drop = FALSE]
    X_val = X_val[complete_rows]
    Y_val = Y_val[complete_rows]
    wt_val = wt_val[complete_rows]
  }
  
  n_sample = length(X_val)
  n_total = sum(wt_val) ## It matches n_sample in the unweighted case
  
  if (is.null(policy)) {
    stop("A 'policy' must be provided.")
  }
  
  # If 'policy' is a single function, wrap it in a list
  if (is.function(policy)) {
    policy <- list(policy)
  }
  
  # Now it must be a list
  if (!is.list(policy)) {
    stop("'policy' must be a function or a list of functions.")
  }
  
  # Ensure all elements of the list are functions
  if (!all(sapply(policy, is.function))) {
    stop("All elements of 'policy' must be functions.")
  }
  
  if (!is.null(ind_S)) {
    if (!is.numeric(ind_S) || length(ind_S) != nrow(data) || any(!(ind_S %in% c(0, 1)))) {
      stop("When provided, 'ind_S' must be a numeric vector of 0s and 1s with length equal to the number of rows in the original dataset.")
    }
    ind_S_val = ind_S[complete_rows]
  }else{
    ind_S_val = rep(1, n_sample)
  }
  
  q_names = sapply(1:length(policy), FUN = function(x){
    paste0("Policy_q^", x)
  })
  
  mean_A = mean(X_val)
  sd_A = sd(X_val)
  A_val = (X_val-mean_A)/sd_A
  
  q_list = lapply(policy, function(f) (f(X_val) - mean_A) / sd_A)
  names(q_list) = q_names
  
  if (is.null(ind_qS)) {
    # Define support of X as its observed range
    x_min <- min(X_val)
    x_max <- max(X_val)
    
    X_target <- X_val[ind_S_val == 1]
    
    ind_qS <- lapply(seq_along(policy), function(m) {
      qS_vals <- policy[[m]](X_target)  # Evaluate policy on the target population S
      
      1 * (X_val >= min(qS_vals) & X_val <= max(qS_vals))
      
    })
  } else {
    # Ensure it's either a list or convert singleton to list
    if (!is.list(ind_qS)) {
      if (length(ind_qS) == nrow(data) && is.numeric(ind_qS[[1]])) {
        ind_qS = rep(list(ind_qS), length(policy))
      } else {
        stop("'ind_qS' must be a list of numeric indicator vectors or a single such vector of length equal to the number of rows in the original dataset.")
      }
    }
    
    if (length(ind_qS) != length(policy)) {
      stop("The length of 'ind_qS' must match the number of policies provided.")
    }
    
    # Validate structure of each element
    for (m in seq_along(ind_qS)) {
      v <- ind_qS[[m]]
      if (!is.numeric(v) || length(v) != nrow(data) || any(!(v %in% c(0, 1)))) {
        stop(paste0("Each 'ind_qS[[", m, "]]' must be a numeric vector of 0s and 1s of length equal to the number of rows in the original dataset."))
      }
      # Subset to complete rows
      ind_qS[[m]] <- v[complete_rows]
    }
  }
  
  ### Standardizing negative controls
  working_data = list(
    L = scale(L_val),
    Z = scale(Z_val),
    W = scale(W_val),
    X = X_val,
    Y = Y_val,
    wt = wt_val,
    ind_S = ind_S_val,
    A = A_val
  )
  
  working_data[q_names] = q_list
  
  K_folds = K_folds
  
  names_folds = sapply(1:K_folds, FUN = function(x){
    paste0("Fold", x)
  })
  
  if(control.folds){
    random_sample = createFolds(X_val, k = K_folds)
  }else{
    random_sample = createFolds(1:n_sample, k = K_folds)
  }
  
  OR_IF = matrix(0, nrow = n_sample, ncol = length(policy))
  IPW_IF = matrix(0, nrow = n_sample, ncol = length(policy))
  DR_IF = matrix(0, nrow = n_sample, ncol = length(policy))
  TMLE_IF = matrix(0, nrow = n_sample, ncol = length(policy))
  
  lm_H_cv_list  = rep(0, K_folds)
  lm_Gh_cv_list  = rep(0, K_folds)
  names(lm_H_cv_list) = names_folds
  names(lm_Gh_cv_list) = names_folds
  lm_G_cv_list  = matrix(0, nrow = K_folds, ncol = length(policy), dimnames = list(names_folds, q_names))
  lm_Hg_cv_list  = matrix(0, nrow = K_folds, ncol = length(policy), dimnames = list(names_folds, q_names))
  
  bw_H_cv_list  = rep(0, K_folds)
  bw_Gh_cv_list  = rep(0, K_folds)
  names(bw_H_cv_list) = names_folds
  names(bw_Gh_cv_list) = names_folds
  bw_G_cv_list  = matrix(0, nrow = K_folds, ncol = length(policy), dimnames = list(names_folds, q_names))
  bw_Hg_cv_list  = matrix(0, nrow = K_folds, ncol = length(policy), dimnames = list(names_folds, q_names))
  
  Gh_norm_cv_list  = rep(0, K_folds)
  H_norm_cv_list  = rep(0, K_folds)
  names(H_norm_cv_list) = names_folds
  names(Gh_norm_cv_list) = names_folds
  Hg_norm_cv_list  = matrix(0, nrow = K_folds, ncol = length(policy), dimnames = list(names_folds, q_names))
  G_norm_cv_list  = matrix(0, nrow = K_folds, ncol = length(policy), dimnames = list(names_folds, q_names))
  
  h.fn = rep(0, n_sample)
  g.fn = matrix(0, nrow = n_sample, ncol = length(policy))
  colnames(g.fn) = q_names
  
  names_L = expand_names(L_val, "L")
  names_Z = expand_names(Z_val, "Z")
  names_W = expand_names(W_val, "W")
  
  var_names = c(names_L, names_Z, names_W, "X", "Y", "wt", "ind_S", "A", q_names)
  
  for(cv_rep in 1:K_folds){
    if(show.prog){
      print(paste("Fold =", cv_rep))
    }
    
    te_inds = random_sample[[cv_rep]]
    eval_inds = random_sample[[cv_rep%%K_folds+1]]
    
    tr_inds = (1:n_sample)[-c(te_inds, eval_inds)]
    
    tr_data = as.data.frame(do.call(cbind, lapply(working_data, function(x) {
      if (is.data.frame(x) || is.matrix(x)) {
        x[tr_inds, , drop = FALSE]
      } else {
        x[tr_inds]
      }
    })))
    
    te_data = as.data.frame(do.call(cbind, lapply(working_data, function(x) {
      if (is.data.frame(x) || is.matrix(x)) {
        x[te_inds, , drop = FALSE]
      } else {
        x[te_inds]
      }
    })))
    
    eval_data = as.data.frame(do.call(cbind, lapply(working_data, function(x) {
      if (is.data.frame(x) || is.matrix(x)) {
        x[eval_inds, , drop = FALSE]
      } else {
        x[eval_inds]
      }
    })))
    
    names(tr_data) = var_names
    names(te_data) = var_names
    names(eval_data) = var_names
    
    n_tr = sum(tr_data$wt)
    n_te = sum(te_data$wt)
    n_eval = sum(eval_data$wt)
    
    if(show.prog){
      print(paste("(n_tr, n_te, n_eval) =", nrow(tr_data), nrow(te_data), nrow(eval_data)))
    }
    
    ident_n_tr = diag(rep(1, nrow(tr_data)))
    ident_n_te = diag(rep(1, nrow(te_data)))
    W_tr = diag(tr_data$wt)
    W_te = diag(te_data$wt)
    
    ### Kernel bandwidths based on data
    if(is.null(bw0_H)){
      if(no_weights){
        bw0_H = median(dist(as.matrix(tr_data[,c("A", names_L, names_W)])))^2/2
      }else{
        bw0_H = wt_euclidean_bw(X = as.matrix(tr_data[,c("A", names_L, names_W)]), w = tr_data$wt)
      }
    }
    if(is.null(bw0_G)){
      if(no_weights){
        bw0_G = median(dist(as.matrix(tr_data[,c("A", names_L, names_Z)])))^2/2
      }else{
        bw0_G = wt_euclidean_bw(X = as.matrix(tr_data[,c("A", names_L, names_Z)]), w = tr_data$wt)
      }
    }
    
    ## Computing Gaussian kernels using training data
    
    ## For estimation of h
    h_arg_tr = as.matrix(tr_data[,c("A", names_L, names_W)])
    g_arg_tr = as.matrix(tr_data[,c("A", names_L, names_Z)])
    
    h_arg_te = as.matrix(te_data[,c("A", names_L, names_W)])
    g_arg_te = as.matrix(te_data[,c("A", names_L, names_Z)])
    
    ## Common kernel basis for CV
    K_Gh_te = rbfkernel(g_arg_te, sigma = bw0_G*bw_int_fixed_scale)
    K_Hg_te = rbfkernel(h_arg_te, sigma = bw0_H*bw_int_fixed_scale)
    
    # CV to select h
    if(no_weights){
      Theta_Gh = 0.25 * n_te * K_Gh_te %*% inv.mat(K_Gh_te + theta * log(n_te) * ident_n_te)
      obj.cv.h = cv_h(Theta=Theta_Gh, n =n_tr, ident_n=ident_n_tr, n_te= n_te,
                      arg_ext = h_arg_tr,
                      arg_int = g_arg_tr,
                      arg_ext_te = h_arg_te, Y = tr_data$Y,
                      Y_te = te_data$Y,
                      lm_int_list = lm_Gh_list, lm_ext_list = lm_H_list, bw0_int =bw0_G,
                      bw0_ext = bw0_H,
                      bw_ext_scale_list = bw_ext_scale_list,
                      bw_int_scale_list = bw_int_scale_list)
    }else{
      m_te = nrow(te_data)
      Theta_Gh = 0.25 * K_Gh_te %*% inv.mat(W_te %*% K_Gh_te/n_te  + theta * log(m_te)/m_te * ident_n_te)
      obj.cv.h = cv_h_w(Theta=Theta_Gh, n =n_tr, m = nrow(tr_data),
                        ident_m=ident_n_tr, n_te= n_te,
                        arg_ext = h_arg_tr,
                        arg_int = g_arg_tr,
                        arg_ext_te = h_arg_te, Y = tr_data$Y,
                        Y_te = te_data$Y, W = W_tr, W_te = W_te,
                        lm_int_list = lm_Gh_list, lm_ext_list = lm_H_list, bw0_int =bw0_G,
                        bw0_ext = bw0_H,
                        bw_ext_scale_list = bw_ext_scale_list,
                        bw_int_scale_list = bw_int_scale_list)
    }
    
    lm_H_cv_list[cv_rep] = obj.cv.h[[1]]
    lm_Gh_cv_list[cv_rep] = obj.cv.h[[2]]
    bw_Gh_cv_list[cv_rep] = obj.cv.h[[3]]
    bw_H_cv_list[cv_rep] = obj.cv.h[[4]]
    Gh_norm_cv_list[cv_rep] = obj.cv.h[[6]]
    H_norm_cv_list[cv_rep] = obj.cv.h[[7]]
    
    ### Checkpoint
    if(show.prog){
      print("Estimation of the outcome bridge function complete.")
    }
    
    ## Evaluation
    h_arg_eval = as.matrix(eval_data[,c("A", names_L, names_W)])
    K_H_tr_eval = rbfkernel(X = h_arg_eval, sigma = obj.cv.h[[4]], Y = h_arg_tr)
    h.fn[eval_inds] = as.vector(K_H_tr_eval %*%obj.cv.h[[5]])
    
    for(q_index in 1:length(policy)){
      qS = diag(ind_qS[[q_index]][tr_inds])
      qS_te = diag(ind_qS[[q_index]][te_inds])
      
      h_mtp_arg_tr = as.matrix(tr_data[,c(q_names[q_index], names_L, names_W)])
      h_mtp_arg_te = as.matrix(te_data[,c(q_names[q_index], names_L, names_W)])
      
      ind_S_tr = diag(tr_data$ind_S)
      ind_S_te = diag(te_data$ind_S)
      ind_S_eval = diag(eval_data$ind_S)
      
      ### Common matrices for CV
      K_H_tilde_te = rbfkernel(h_mtp_arg_te, sigma = bw0_H*bw_int_fixed_scale, Y = h_arg_te)
      
      ## Common term in the loss for CV to select g 
      
      # CV to select g
      if(no_weights){
        A_te = colSums(ind_S_te %*% K_H_tilde_te)
        Theta_Hg = 0.25 * n_te * inv.mat(K_Hg_te) %*% inv.mat(K_Hg_te %*% qS_te + theta * log(n_te) * ident_n_te)
        obj.cv.g = cv_g(Theta=Theta_Hg, n = n_tr,
                        ident_n= ident_n_tr, n_te= n_te,
                        arg_ext = g_arg_tr,
                        arg_int = h_arg_tr,
                        arg_int_q = h_mtp_arg_tr,
                        arg_ext_te = g_arg_te, qS = qS,
                        qS_te = qS_te, ind_S = ind_S_tr,
                        Kint_te = K_Hg_te, A_te = A_te,
                        lm_int_list = lm_Hg_list, lm_ext_list = lm_G_list,
                        bw0_int =bw0_H, bw0_ext = bw0_G,
                        bw_ext_scale_list = bw_ext_scale_list,
                        bw_int_scale_list = bw_int_scale_list)
      }else{
        A_te = colSums(W_te %*% ind_S_te %*% K_H_tilde_te)
        Theta_Hg = 0.25 * inv.mat(K_Hg_te) %*% inv.mat(K_Hg_te %*% W_te %*% qS_te/n_te + theta * log(m_te)/m_te * ident_n_te)
        obj.cv.g = cv_g_w(Theta=Theta_Hg, n = n_tr, m = nrow(tr_data),
                          ident_m= ident_n_tr, n_te= n_te,
                          arg_ext = g_arg_tr,
                          arg_int = h_arg_tr,
                          arg_int_q = h_mtp_arg_tr,
                          arg_ext_te = g_arg_te, qS = qS,
                          qS_te = qS_te, ind_S = ind_S_tr,
                          W = W_tr, W_te = W_te,
                          Kint_te = K_Hg_te, A_te = A_te,
                          lm_int_list = lm_Hg_list, lm_ext_list = lm_G_list,
                          bw0_int =bw0_H, bw0_ext = bw0_G,
                          bw_ext_scale_list = bw_ext_scale_list,
                          bw_int_scale_list = bw_int_scale_list)
      }
      
      lm_G_cv_list[cv_rep, q_index] = obj.cv.g[[1]]
      lm_Hg_cv_list[cv_rep, q_index] = obj.cv.g[[2]]
      bw_Hg_cv_list[cv_rep, q_index] = obj.cv.g[[3]]
      bw_G_cv_list[cv_rep, q_index] = obj.cv.g[[4]]
      Hg_norm_cv_list[cv_rep, q_index] = obj.cv.g[[6]]
      G_norm_cv_list[cv_rep, q_index] = obj.cv.g[[7]]
      
      ### Checkpoint
      if(show.prog){
        print(paste0("Estimation of the treatment bridge function for policy ", q_index, " complete."))
      }
      
      ### For evaluation of h
      h_mtp_arg_eval = as.matrix(eval_data[,c(q_names[q_index], names_L, names_W)])
      K_H_tilde_tr_eval = rbfkernel(X = h_mtp_arg_eval, sigma = obj.cv.h[[4]], Y = h_arg_tr)
      
      ### For evaluation of g
      g_arg_eval = as.matrix(eval_data[,c("A", names_L, names_Z)])
      K_G_tr_eval = rbfkernel(X = g_arg_eval, sigma =obj.cv.g[[4]], Y = g_arg_tr)
      
      ## Evaluation OR ##################################################################
      OR_IF[eval_inds, q_index] = ind_S_eval%*%as.vector(K_H_tilde_tr_eval%*%obj.cv.h[[5]])
      
      ## Evaluation IPW #################################################################
      g0_aux = as.vector(K_G_tr_eval %*% obj.cv.g[[5]])
      
      IPW_IF[eval_inds, q_index] = g0_aux*as.vector(eval_data$Y)*ind_qS[[q_index]][eval_inds]
      
      g.fn[eval_inds, q_index] = g0_aux
      
      ## Evaluation DR ##################################################################
      DR_IF[eval_inds, q_index] = as.vector(OR_IF[eval_inds, q_index])+as.vector(g0_aux)*as.vector((eval_data$Y-h.fn[eval_inds]))*ind_qS[[q_index]][eval_inds]
    }
  }
  
  ### TMLE computation
  if(compute.TMLE.est){
    for(q_index in 1:length(policy)){
      data.aux = data.frame(
        Y = Y_val,
        wt = wt_val,
        h0_bounded = pmax(pmin(h.fn,1-1e-5),1e-5),
        hq_bounded = pmax(pmin(OR_IF[, q_index],1-1e-5),1e-5),
        g0 = g.fn[,q_index]*ind_qS[[q_index]]
      )
      if(no_weights){
        obj.tmle = glm(Y ~ 0+g0, offset = qlogis(h0_bounded), family = binomial(link = "logit"),
                       data = data.aux[data.aux$g0>0,])
      }else{
        obj.tmle = glm(Y ~ 0+g0, offset = qlogis(h0_bounded), weights = wt,
                       family = binomial(link = "logit"),   data = data.aux[data.aux$g0>0,])
      }
      
      h1 = plogis(obj.tmle$coefficients*data.aux$g0+qlogis(data.aux$h0_bounded))
      hq1 = plogis(obj.tmle$coefficients*data.aux$g0+qlogis(data.aux$hq_bounded))
      
      TMLE_IF[, q_index] = hq1 + data.aux$g0*(Y_val - h1)
    }
  }
  
  p_hat = mean(ind_S_val)
  
  if(no_weights){
    OR_est = sapply(data.frame(OR_IF), mean)/p_hat
    IPW_est = sapply(data.frame(IPW_IF), mean)/p_hat
    DR_est = sapply(data.frame(DR_IF), mean)/p_hat
    TMLE_est = sapply(data.frame(TMLE_IF), mean)/p_hat
  }else{
    OR_est = sapply(as.data.frame(OR_IF), FUN = function(x){
      weighted.mean(x, w = wt_val)
    })/p_hat
    IPW_est = sapply(as.data.frame(IPW_IF), FUN = function(x){
      weighted.mean(x, w = wt_val)
    })/p_hat
    DR_est = sapply(as.data.frame(DR_IF), FUN = function(x){
      weighted.mean(x, w = wt_val)
    })/p_hat
    TMLE_est = sapply(as.data.frame(TMLE_IF), FUN = function(x){
      weighted.mean(x, w = wt_val)
    })/p_hat
  }
  
  ### Preparing final object
  OR_IF = OR_IF/p_hat - t(OR_est %*% t(ind_S_val))/p_hat
  IPW_IF = IPW_IF/p_hat - t(IPW_est %*% t(ind_S_val))/p_hat
  DR_IF = DR_IF/p_hat - t(DR_est %*% t(ind_S_val))/p_hat
  TMLE_IF = TMLE_IF/p_hat - t(TMLE_est %*% t(ind_S_val))/p_hat
  
  #### Elements for the final object
  if(compute.TMLE.est){
    res1 = data.frame(rbind(DR_est, TMLE_est))
    colnames(res1) = q_names
    rownames(res1) = c("One-step estimator", "TMLE estimator")
  }else{
    res1 = as.data.frame(t(as.numeric(DR_est))) 
    colnames(res1) = q_names
    rownames(res1) = c("One-step estimator")
  }
  
  if(compute.TMLE.est){
    res2 = sapply(1:length(policy), FUN = function(x){
      c(sum(wt_val*DR_IF[,x]^2)/n_total, sum(wt_val*TMLE_IF[,x]^2)/n_total)
    })
    colnames(res2) = q_names
    rownames(res2) = c("One-step Est.Var", "TMLE Est.Var")
  }else{
    res2 = sapply(1:length(policy), FUN = function(x){
      sum(wt_val*DR_IF[,x]^2)/n_total
    })
    res2 = as.data.frame(t(res2))  # transpose to get a row
    colnames(res2) = q_names
    rownames(res2) = c("One-step Est.Var")
  }
  
  if(!inf.fun){
    res.out = list(res1, res2, n_total, n_sample, ind_S_val, ind_qS, lm_H_cv_list, lm_Gh_cv_list, lm_G_cv_list, lm_Hg_cv_list,
                   bw_H_cv_list, bw_Gh_cv_list, bw_G_cv_list, bw_Hg_cv_list,
                   Gh_norm_cv_list, H_norm_cv_list, Hg_norm_cv_list, G_norm_cv_list)
    names(res.out) = c("Estimators", "Var.est", "n", "n_val", "ind_S", "ind_qS", "lm_H_cv_list", "lm_Gh_cv_list", "lm_G_cv_list", "lm_Hg_cv_list",
                       "bw_H_cv_list", "bw_Gh_cv_list", "bw_G_cv_list", "bw_Hg_cv_list",
                       "Gh_norm_cv_list", "H_norm_cv_list", "Hg_norm_cv_list", "G_norm_cv_list")
  }else{
    res.out = list(res1, res2, n_total, n_sample, ind_S_val, ind_qS, lm_H_cv_list, lm_Gh_cv_list, lm_G_cv_list, lm_Hg_cv_list,
                   bw_H_cv_list, bw_Gh_cv_list, bw_G_cv_list, bw_Hg_cv_list,
                   DR_IF, h.fn, g.fn, Gh_norm_cv_list, H_norm_cv_list, Hg_norm_cv_list, G_norm_cv_list)
    names(res.out) = c("Estimators", "Var.est", "n", "n_val", "ind_S", "ind_qS", "lm_H_cv_list", "lm_Gh_cv_list", "lm_G_cv_list", "lm_Hg_cv_list",
                       "bw_H_cv_list", "bw_Gh_cv_list", "bw_G_cv_list", "bw_Hg_cv_list", 
                       "DR_IFs", "h.fn.est", "g.fn.est","Gh_norm_cv_list", "H_norm_cv_list", "Hg_norm_cv_list",
                       "G_norm_cv_list")
  }
  
  class(res.out) <- "proximal.mtp"
  print(res.out$Estimators)
  invisible(res.out)
}



###############################################

print.proximal.mtp <- function(x, sign.level = 0.05) {
  est <- x$Estimators
  se <- sqrt(x$Var.est / x$n)
  z <- qnorm(1-sign.level/2)
  est_names = c("Proximal DR cross-fitted", "Proximal TMLE cross-fitted")
  
  cat(paste0("Proximal Doubly Robust Estimators and ", round(100*(1-sign.level)),"% Confidence Intervals:\n\n"))
  
  for (i in 1:nrow(est)) {
    est_row <- unlist(est[i, ])
    se_row  <- unlist(se[i, ])
    
    lower <- est_row - z * se_row
    upper <- est_row + z * se_row
    
    ci_df <- data.frame(
      Estimate = est_row,
      `Std Error` = se_row,
      `CI Lower` = lower,
      `CI Upper` = upper,
      row.names = names(est)
    )
    
    cat(paste0("Estimator: ", est_names[i], "\n"))
    print(round(ci_df, 4))
    cat("\n")
  }
  
  invisible(x)
}

summary.proximal.mtp <- function(x, sign.level = 0.05) {
  est <- x$Estimators
  se <- sqrt(x$Var.est / x$n)
  z <- qnorm(1-sign.level/2)
  est_names = c("Proximal doubly-robust cross-fitted", "Proximal TMLE cross-fitted")
  
  cat("Summary of Proximal MTP Estimation\n")
  cat(strrep("-", 50), "\n")
  
  cat(paste0("Number of observations: ", x$n_val, "\n"))
  cat(paste0("Number in target population: ", sum(x$ind_S), " (", round(100*mean(x$ind_S)), "%)\n"))
  for (i in 1:ncol(est)){
    cat(paste0("-Proportion in the image of ", names(est)[i], ": ", round(100*mean(x$ind_qS[[i]][x$ind_S==1])),"%\n"))
  }
  cat("\n")
  cat(strrep("-", 50), "\n")
  
  cat(paste0("Proximal Estimators and ", round(100*(1-sign.level)),"% Confidence Intervals:\n"))
  cat(strrep("-", 50), "\n")
  for (i in 1:nrow(est)) {
    est_row <- unlist(est[i, ])
    se_row  <- unlist(se[i, ])
    
    lower <- est_row - z * se_row
    upper <- est_row + z * se_row
    
    ci_df <- data.frame(
      Estimate = est_row,
      `Std Error` = se_row,
      `CI Lower` = lower,
      `CI Upper` = upper,
      row.names = names(est)
    )
    
    cat(paste0("Estimator: ", est_names[i], "\n"))
    print(round(ci_df, 4))
    cat("\n")
  }
  
  
  invisible(x)
}



