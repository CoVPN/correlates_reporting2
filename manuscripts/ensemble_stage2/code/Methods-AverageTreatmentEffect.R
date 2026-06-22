require(hal9001)

EstimateCFMeans <- function(Y, A, W, weights = NULL) {
  
  miss <- is.na(Y) | is.na(A) | apply(W, 1, function(x) !all(!is.na(x)))
  Y <- Y[!miss]
  A <- A[!miss]
  W <- W[!miss,]
  
  outcome.regression <- fit_hal(X = cbind(A,W), Y = Y, weights = weights)
  propensity.score <- fit_hal(X = W, Y = A, family = "binomial",
                              weights = weights)
  
  Y.hat <- predict(outcome.regression, new_data = cbind(A,W))
  Y1.hat <- predict(outcome.regression, new_data = cbind(1,W))
  Y0.hat <- predict(outcome.regression, new_data = cbind(0,W))
  A.hat <- predict(propensity.score, new_data = W)
  N <- length(Y)
  
  cf.mean.1 <- sum(Y1.hat * weights)/sum(weights) + sum(A * ((Y - Y.hat)/A.hat) * weights)/sum(weights)
  cf.mean.0 <- sum(Y0.hat * weights)/sum(weights) + sum((1 - A) * ((Y - Y.hat)/(1 - A.hat)) * weights)/sum(weights)
  ate <- cf.mean.1 - cf.mean.0
  eif.1 <- (Y1.hat + A * ((Y - Y.hat)/A.hat)) * weights/mean(weights)
  eif.0 <- (Y0.hat + (1 - A) * ((Y - Y.hat)/(1 - A.hat))) * weights/mean(weights)
  eif.1 <- eif.1 - mean(eif.1)
  eif.0 <- eif.0 - mean(eif.0)
  eif.ate <- eif.1 - eif.0
  
  out <- list(cf.mean.0 = cf.mean.0,
              cf.mean.1 = cf.mean.1,
              ate = ate,
              eif.0 = eif.0,
              eif.1 = eif.1,
              eif.ate = eif.ate)
}
