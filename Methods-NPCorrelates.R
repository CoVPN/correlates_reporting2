# Estimate nuisance:
## 1) Conditional density
## 2) Conditional survival [HAL w/ Cox PH, since it allows sampling weights]
# Estimate pseudo outcome [Ted's paper]
# Kernel ridge regression estimator (can get sup confidence bands)
require(mvtnorm)
require(hal9001)
require(survival)
logit <- function(x) log(x/(1-x))
dlogit <- function(x) 1/(x * (1-x))
expit <- function(x) exp(x)/(1 + exp(x))

EstimateConditionalDensity <- function(a, w, no.folds = 10, no.eval = 20, bw = NULL, weights = NULL) {
  n <- length(a)
  a.eval <- seq(min(a), max(a), length.out = no.eval)
  
  # use bandwidth that is optimal for kernel estimation of margianl density
  # is sensible choice if con'l density is about as smooth in a as marginal density
  # could try other choices via cross-validation, but this is comp'l cheaper
  # con'l density estimation is generally hard, and needs further study
  
  if(is.null(bw)) {
    # select bandwidth via cross-validation
    folds <- sample(cut(1:n, breaks = no.folds, labels = FALSE))
    bw.seq <- exp(seq(log(.01 * sd(a)), log(25 * sd(a)), length.out = 100))
    
    pred.error <- matrix(NA, nrow = 100, ncol = no.folds) # prediction error
    for(i in 1:no.folds) {
      for(j in 1:100) {
        a.train <- a[folds != i]
        a.test <- a[folds == i]
        
        n.test <- length(a.test)
        fit.test <- numeric(n.test)
        
        for(k in 1:n.test) {
          kern <- dnorm(a[folds != i], mean = a.test[k], sd = bw.seq[j])
          fit.test[k] <- sum(kern * weights[folds != i])/sum(weights[folds != i])
        }
        
        pred.error[j,i] <- sum(-log(fit.test) * weights[folds == i])/sum(weights[folds == i])
      }
    }
    avg.pred.error <- apply(pred.error, 1, mean)
    se.pred.error <- apply(pred.error, 1, sd)/sqrt(no.folds)
    bw.opt <- max(bw.seq[avg.pred.error < min(avg.pred.error) + se.pred.error[which.min(avg.pred.error)]])
    # bw.opt <- bw.seq[which.min(apply(pred.error, 1, mean))]
  } else {
    # set bandwidth to a pre-specified value, if one is provided
    bw.opt <- bw
  }
  
  marg.dens <- numeric(n)
  for(i in 1:n) {
    marg.dens[i] <- mean(dnorm(a, mean = a[i], sd = 2 * bw.opt))
  }
  
  f.eval <- matrix(NA, nrow = n, ncol = no.eval)
  for(j in 1:no.eval) {
    kern.aj <- dnorm(a, mean = a.eval[j], sd = bw.opt)
    # hal.fit <- fit_hal(X = w, Y = kern.aj)
    hal.fit <- fit_hal(X = w, Y = kern.aj,
                       family = "poisson", return_x_basis = TRUE,
                       weights = weights)
    f.eval[,j] <- predict(hal.fit, new_data = w)
    # f.eval[,j] <- exp(hal.fit$X.basis %*% hal.fit$coefs)
  }
  
  # set cond.dens[i,j] = cond.dens(a[i], w[j])
  cond.dens <- matrix(NA, n, n)
  for(i in 1:n) {
    k <- min(which(a.eval - a[i] >= 0))
    for(j in 1:n) {
      if(k == 1) {
        cond.dens[i,j] <- f.eval[j,1]
      } else if(k == no.eval) {
        cond.dens[i,j] <- f.eval[j,no.eval]
      } else {
        t <- (a.eval[k] - a[i])/(a.eval[k] - a.eval[k-1])
        cond.dens[i,j] <- f.eval[j,k] * (1-t) + f.eval[j,k-1] * t
      }
    }
  }
  
  # marg.dens <- apply(cond.dens, 1, mean)
  marg.dens <- cond.dens %*% weights/sum(weights)
  prop.score <- c(marg.dens * diag(cond.dens^(-1)))
  
  out <- list(cond.dens = cond.dens,
              marg.dens = marg.dens,
              prop.score = prop.score)
  return(out)
}

EstimateConditionalSurvival <- function(time, event, a, w = NULL, weights = NULL, 
                                        t0 = quantile(time[event == 1], seq(.1, .9, .2)),
                                        nonparametric = TRUE) {
  # One option is to re-fit the cox-ph using the HAL predictions, then obtain the
  # baseline hazard using the basehaz function. 
  
  n <- length(time)
  
  # Estimate hazard ratio (relative to baseline hazard) using proportional hazards regression
  if(!is.null(w)) {
    if(nonparametric) {
      fit.event <- fit_hal(Y = Surv(time = time, event = event), X = cbind(a, w), family = "cox",
                           weights = weights)
      fit.cens <- fit_hal(Y = Surv(time = time, event = 1 - event), X = cbind(a, w), family = "cox",
                          weights = weights)
      
      hr.event <- exp(predict(fit.event, new_data = cbind(a,w), type = "link"))
      hr.cens <- exp(predict(fit.cens, new_data = cbind(a,w), type = "link"))
    } else {
      fit.event <- coxph(Surv(time = time, event = event) ~ a + w, 
                         weights = weights, robust = TRUE)
      fit.cens <- coxph(Surv(time = time, event = 1 - event) ~ a + w,
                        weights = weights, robust = TRUE)
      
      hr.event <- c(exp(cbind(a, w) %*% fit.event$coef))
      hr.cens <- c(exp(cbind(a, w) %*% fit.cens$coef))
    }
  } else {
    if(nonparametric) {
      fit.event <- fit_hal(Y = Surv(time = time, event = event), X = a, family = "cox",
                           weights = weights)
      fit.cens <- fit_hal(Y = Surv(time = time, event = 1 - event), X = a, family = "cox",
                          weights = weights)
      
      hr.event <- exp(predict(fit.event, new_data = a, type = "link"))
      hr.cens <- exp(predict(fit.cens, new_data = a, type = "link"))
    } else {
      fit.event <- coxph(Surv(time = time, event = event) ~ a, 
                         weights = weights, robust = TRUE)
      fit.cens <- coxph(Surv(time = time, event = 1 - event) ~ a,
                        weights = weights, robust = TRUE)
      
      hr.event <- c(exp(cbind(a) %*% fit.event$coef))
      hr.cens <- c(exp(cbind(a) %*% fit.cens$coef))
    }
  }
  
  # compute Breslow estimator for the baseline cumulative hazard
  base.cl.haz.event <- length(n)
  base.cl.haz.cens <- length(n)
  
  base.cl.haz.event.t0 <- length(n)
  base.cl.haz.cens.t0 <- length(n)
  
  base.haz.event <- length(n)
  base.haz.cens <- length(n)
  
  for(i in 1:n) {
    risk.set <- which(time >= time[i])
    base.haz.event[i] <- weights[i]/sum((hr.event * weights)[risk.set])
    base.haz.cens[i] <- weights[i]/sum((hr.cens * weights)[risk.set])
  }
  for(i in 1:n) {
    base.cl.haz.event[i] <- sum(base.haz.event[time <= time[i]])
    base.cl.haz.cens[i] <- sum(base.haz.cens[time <= time[i]])
  }
  for(i in 1:length(t0)) {
    base.cl.haz.event.t0[i] <- sum(base.haz.event[time <= t0[i]])
    base.cl.haz.cens.t0[i] <- sum(base.haz.cens[time <= t0[i]])
  }
  
  # estimate cumulative hazard
  # cond.surv <- array(NA, dim = c(n, n, n))
  # cond.cens <- array(NA, dim = c(n, n, n))
  cond.surv.t0 <- array(NA, dim = c(n, n, length(t0)))
  cond.cens.t0 <- array(NA, dim = c(n, n, length(t0)))
  fitted.surv <- matrix(NA, n, n)
  fitted.cens <- matrix(NA, n, n)
  fitted.surv.t0 <- matrix(NA, n, length(t0))
  fitted.cens.t0 <- matrix(NA, n, length(t0))
  fitted.haz.event <- matrix(NA, n, n)
  fitted.haz.cens <- matrix(NA, n, n)
  
  for(i in 1:n) {
    fitted.surv[i,] <- exp(-base.cl.haz.event * hr.event[i])
    fitted.cens[i,] <- exp(-base.cl.haz.cens * hr.event[i])
    
    fitted.surv.t0[i,] <- exp(-base.cl.haz.event.t0 * hr.event[i])
    fitted.cens.t0[i,] <- exp(-base.cl.haz.cens.t0 * hr.event[i])
    
    fitted.haz.cens[i,] <- hr.cens[i] * base.haz.cens
    fitted.haz.event[i,] <- hr.event[i] * base.haz.event
    if(!is.null(w)) {
      if(nonparametric) {
        hr.event.i <-  exp(predict(fit.event, new_data = cbind(a[i],w), type = "link"))
        hr.cens.i <- exp(predict(fit.cens, new_data = cbind(a[i],w), type = "link"))
      } else {
        hr.event.i <- exp(c(cbind(a[i], w) %*% fit.event$coef))
        hr.cens.i <- exp(c(cbind(a[i], w) %*% fit.cens$coef))
      }
    } else {
      if(nonparametric) {
        hr.event.i <-  exp(predict(fit.event, new_data = cbind(a[i]), type = "link"))
        hr.cens.i <- exp(predict(fit.cens, new_data = cbind(a[i]), type = "link"))
      } else {
        hr.event.i <- rep(exp(c(a[i] * fit.event$coef)), length(event))
        hr.cens.i <- rep(exp(c(a[i] * fit.cens$coef)), length(event))
      }
    }
    for(j in 1:n) {
      # cond.surv[i,j,] <- exp(-base.cl.haz.event * hr.event.i[j])
      # cond.cens[i,j,] <- exp(-base.cl.haz.cens * hr.cens.i[j])
      
      cond.surv.t0[i,j,] <- exp(-base.cl.haz.event.t0 * hr.event.i[j])
      cond.cens.t0[i,j,] <- exp(-base.cl.haz.cens.t0 * hr.cens.i[j])
    }
  }
  
  # out <- list(fitted.surv = fitted.surv,
  #             fitted.cens = fitted.cens,
  #             fitted.surv.t0 = fitted.surv.t0,
  #             fitted.cens.t0 = fitted.cens.t0,
  #             # cond.surv = cond.surv,
  #             # cond.cens = cond.cens,
  #             cond.surv.t0 = cond.surv.t0,
  #             cond.cens.t0 = cond.cens.t0,
  #             fitted.haz.event = fitted.haz.event,
  #             fitted.haz.cens = fitted.haz.cens)
  
  out <- list(cond.surv = cond.surv.t0,
              cond.cens = cond.cens.t0,
              fitted.surv = fitted.surv,
              fitted.cens = fitted.cens,
              fitted.haz.event = fitted.haz.event,
              fitted.haz.cens = fitted.haz.cens)
  
  return(out)
  
}

EstiamteLinearContrast <- function(time, event, a, w, t0, basis,
                                   weights = rep(1, length(time)),
                                   fitted.surv, fitted.cens,
                                   cond.surv, cond.cens,
                                   fitted.haz.event, fitted.haz.cens,
                                   prop.score, logit.scale = FALSE) {
  
  n <- length(time)
  basis.cent <- (diag(n) - matrix(weights/sum(weights), n, n, byrow = TRUE)) %*% basis
  
  plug.in <- matrix(NA, nrow = ncol(basis), ncol = length(t0))
  avg.plug.in <- numeric(length(t0))
  
  one.step <- matrix(NA, nrow = ncol(basis), ncol = length(t0))
  avg.one.step <- numeric(length(t0))
  one.step.cent <- matrix(NA, nrow = ncol(basis), ncol = length(t0))
  
  influence.function <- array(NA, dim = c(n, ncol(basis), length(t0)))
  avg.influence.function <- matrix(NA, nrow = n, ncol = length(t0))
  influence.function.cent <- array(NA, dim = c(n, ncol(basis), length(t0)))
  
  pseudo.outcome <- matrix(NA, nrow = n, ncol = length(t0))
  
  for(s in 1:length(t0)) {
    
    # get plug-in estimate
    plug.in[,s] <- t(basis) %*% diag(weights/sum(weights)) %*% cond.surv[,,s] %*% weights/sum(weights)
    avg.plug.in[s] <- t(weights/sum(weights)) %*% cond.surv[,,s] %*% weights/sum(weights)
    
    # obtain influence.function and pseudo outcome
    r <- numeric(n)
    for(i in 1:n) {
      r[i] <- sum((fitted.haz.event[i,]/(fitted.surv[i,] * fitted.cens[i,]) * 
                   ifelse(time <= min(c(t0[s], time[i])), 1, 0)))
    }
    # resid <- (ifelse(time <= t0[s] & event == 1, 1, 0)/diag(cond.surv[,,s] * cond.cens[,,s]) -
    #          (fitted.haz.event/(fitted.surv * fitted.cens)) %*% ifelse(time <= t0[s], 1, 0)) *
    #          (-diag(cond.surv[,,s]))
    resid <- (ifelse(time <= t0[s] & event == 1, 1, 0)/diag(fitted.surv * fitted.cens) - r) *
             (-diag(cond.surv[,,s]))
    
    pseudo.outcome[,s] <- cond.surv[,,s] %*% (weights/sum(weights)) +
                          prop.score * resid
    # pseudo.outcome[,s] <- (diag(cond.surv[,,s]) + prop.score * resid)
    
    influence.function[,,s] <- diag(c(cond.surv[,,s] %*% (weights/sum(weights)))) %*% basis +
                               t(cond.surv[,,s]) %*% diag(weights/sum(weights)) %*% basis +
                               diag(c((prop.score) * resid)) %*% basis -
                               matrix(2 * plug.in[,s], byrow = TRUE, nrow = n, ncol = ncol(basis))
    influence.function[,,s] <- diag(c(weights/mean(weights))) %*% influence.function[,,s]
                               
    
    avg.influence.function[,s] <- cond.surv[,,s] %*% (weights/sum(weights)) +
                                  t(cond.surv[,,s]) %*% (weights/sum(weights)) +
                                  prop.score * resid -
                                  2 * rep(avg.plug.in[s], n)
    # avg.influence.function[,s] <- avg.influence.function[,s] * weights/mean(weights) - 
    #                               (weights/mean(weights) - 1) %*% t(avg.plug.in[s])
    avg.influence.function[,s] <- avg.influence.function[,s] * weights/mean(weights)
    
    # obtain one-step estimators
    mean.basis <- t(basis) %*% weights/sum(weights)
    one.step[,s] <- plug.in[,s] + t(influence.function[,,s]) %*% rep(1/n, n) #(weights/sum(weights))
    avg.one.step[s] <- avg.plug.in[s] + mean(avg.influence.function[,s])
    one.step.cent[,s] <- one.step[,s] - t(basis) %*% weights/sum(weights) * avg.one.step[s]
    
    influence.function.cent[,,s] <- influence.function[,,s] - 
                                    avg.one.step[s] * (diag(weights/mean(weights)) %*% basis.cent) -
                                    matrix(avg.influence.function[,s], nrow = n, ncol = ncol(basis), byrow = FALSE) %*%
                                    diag(c(t(basis) %*% weights/sum(weights)))
  }
  
  if(logit.scale){
    logit.plug.in <- matrix(NA, nrow = ncol(basis), ncol = length(t0))
    logit.one.step <- matrix(NA, nrow = ncol(basis), ncol = length(t0))
    logit.influence.function <- array(NA, dim = c(n, ncol(basis), length(t0)))
    
    for(s in 1:length(t0)) {
      logit.plug.in[,s] <- logit(plug.in[,s])
      logit.influence.function[,,s] <- influence.function[,,s] %*% diag(dlogit(plug.in[,s]))
      logit.one.step[,s] <- logit.plug.in[,s] + t(logit.influence.function[,,s]) %*% rep(1/n, n)
    }
   } else{
      logit.plug.in <- NULL
      logit.one.step <- NULL
      logit.influence.function <- NULL 
  }
  
  # colnames(plug.in) <- paste0("time = ", t0)
  # colnames(one.step) <- paste0("time = ", t0)
  # colnames(one.step.cent) <- paste0("time = ", t0)
  # dimnames(influence.function[[3]]) <- paste0("time = ", t0)
  # dimnames(influence.function.cent) <- paste0("time = ", t0)
  
  out <- list(plug.in = plug.in,
              avg.plug.in = avg.plug.in,
              one.step = one.step,
              avg.one.step = avg.one.step,
              influence.function = influence.function,
              avg.influence.function = avg.influence.function,
              one.step.cent = one.step.cent,
              influence.function.cent = influence.function.cent,
              logit.plug.in = logit.plug.in,
              logit.one.step = logit.one.step,
              logit.influence.function = logit.influence.function,
              t0 = t0, pseudo.outcome = pseudo.outcome)
}

KernelSmoothFit <- function(time, event, a, w, 
                            no.mc = 5000, no.folds = 10,
                            a.eval = unique(quantile(a, seq(.01, .99, .01))),
                            weights = rep(1, length(a)), t0, 
                            fitted.surv, fitted.cens,
                            cond.surv, cond.cens,
                            fitted.haz.event, fitted.haz.cens,
                            prop.score, alpha = .05) {

  n <- length(time)
  
  # obtain pseudo outcome
  pseudo.outcome <- matrix(NA, nrow = n, ncol = length(t0))
  for(s in 1:length(t0)) {
    r <- numeric(n)
    for(i in 1:n) {
      r[i] <- sum((fitted.haz.event[i,]/(fitted.surv[i,] * fitted.cens[i,]) * 
                     ifelse(time <= min(c(t0[s], time[i])), 1, 0)))
    }
    resid <- (ifelse(time <= t0[s] & event == 1, 1, 0)/diag(fitted.surv * fitted.cens) - r) *
             (-diag(cond.surv[,,s]))
    
    pseudo.outcome[,s] <-   cond.surv[,,s] %*% (weights/sum(weights)) +
                            prop.score * resid
  }
  
  # regress pseudo outcomes on x using kernel ridge regression
  # repeat for each t0
  folds <- sample(1:no.folds, n, replace = TRUE)
  # bw.seq <- exp(seq(log(sd(a) * .1), log(sd(a) * 10), length.out = 20))
  bw.seq <- exp(seq(log(sd(a) * n^(-1/5) * .1), log(sd(a) * n^(-1/5) * 10), length.out = 20))
  test.error <- array(NA, dim = c(no.folds, length(bw.seq), length(t0)))
  
  
  for(k in 1:length(bw.seq)) {
    basis <- matrix(NA, nrow = n, ncol = n)
    for(i in 1:n) {
      basis[i,] <- exp(-(a - a[i])^2/(bw.seq[k]))
    }
    for(j in 1:no.folds) {
      train <- which(folds != j)
      test <- which(folds == j)
      mean.basis <- c(t(basis[train,]) %*% (weights[train]/sum(weights[train])))[test]
      est.lc <- EstiamteLinearContrast(time[train], event[train],
                                       basis = basis[train,test] %*% diag(1/mean.basis),
                                       a = a[train], w = w[train,], 
                                       weights = weights[train],
                                       t0 = t0, 
                                       fitted.surv = fitted.surv[train,train],
                                       fitted.cens = fitted.cens[train,train],
                                       cond.surv = cond.surv[train,train,],
                                       cond.cens = cond.cens[train,train,],
                                       fitted.haz.event = fitted.haz.event[train,train],
                                       fitted.haz.cens = fitted.haz.cens[train,train], 
                                       prop.score = prop.score[train],
                                       logit.scale = TRUE)
  
      for(s in 1:length(t0)) {
        drf.est <- expit(est.lc$logit.one.step[,s])
        test.error[j,k,s] <- sum((pseudo.outcome[test] - drf.est)^2*weights[test])/sum(weights[test])
      }
    }
  }
  
  drf.est <- matrix(NA, nrow = length(a.eval), ncol = length(t0))
  logit.drf.influence.function <- array(NA, dim = c(n, length(a.eval), length(t0))) 
  drf.cb.lower <- matrix(NA, nrow = length(a.eval), ncol = length(t0))
  drf.cb.upper <- matrix(NA, nrow = length(a.eval), ncol = length(t0))
  
  for(s in 1:length(t0)) {
    bw.opt <- max(bw.seq[which.min(apply(test.error[,,s], 2, mean))])
    
    basis <- matrix(NA, nrow = n, ncol = length(a.eval))
    for(i in 1:length(a.eval)) {
      basis[,i] <- exp(-(a - a.eval[i])^2/bw.opt)
    }
    basis.cent <- (diag(n) - matrix(weights/sum(weights), n, n, byrow = TRUE)) %*% basis
    
    mean.basis <- c(t(basis) %*% (weights/sum(weights)))
    est.lc <- EstiamteLinearContrast(time, event,
                                     basis = basis %*% diag(1/mean.basis),
                                     a = a, w = w, 
                                     weights = weights,
                                     t0 = t0, 
                                     fitted.surv = fitted.surv,
                                     fitted.cens = fitted.cens,
                                     cond.surv = cond.surv,
                                     cond.cens = cond.cens,
                                     fitted.haz.event = fitted.haz.event,
                                     fitted.haz.cens = fitted.haz.cens, 
                                     prop.score = prop.score,
                                     logit.scale = TRUE)
    
    drf.est[,s] <- expit(est.lc$logit.one.step[,s])
    # logit.drf.influence.function[,,s] <- (est.lc$influence.function[,,s] -  
    #                                       diag(weights/mean(weights)) %*% basis.cent %*% diag(c(est.lc$one.step[,s]/mean.basis)) +
    #                                       (weights/mean(weights) - 1) %*% t(mean.basis) %*% diag(c(est.lc$one.step[,s]/mean.basis)) ) %*%
    #                                      diag(dlogit(est.lc$plug.in[,s]))
    logit.drf.influence.function[,,s] <- (est.lc$influence.function[,,s] -  
                                          diag(weights/mean(weights)) %*% basis.cent %*% diag(c(est.lc$one.step[,s]/mean.basis))) %*%
                                          diag(dlogit(est.lc$plug.in[,s]))
    logit.drf.cov <- cov(logit.drf.influence.function[,,s])
    drf.cb.lower[,s] <- expit(est.lc$logit.one.step[,s] - qnorm(1-alpha/2) * sqrt(diag(logit.drf.cov)/n))
    drf.cb.upper[,s] <- expit(est.lc$logit.one.step[,s] + qnorm(1-alpha/2) * sqrt(diag(logit.drf.cov)/n))
  }
  
  out <- list(drf.est = drf.est,
              a.eval = a.eval, t0 = t0,
              drf.cb.lower = drf.cb.lower,
              drf.cb.upper = drf.cb.upper)
  return(out)
  
}

MonotoneTest <- function(time, event, a, w, 
                         thresholds = quantile(a, seq(.01, .99, .01)),
                         no.mc = 10000,
                         weights = rep(1, length(event)), t0, 
                         fitted.surv,
                         fitted.cens,
                         cond.surv,
                         cond.cens,
                         fitted.haz.event,
                         fitted.haz.cens,
                         prop.score) {
  
  n <- length(event)
  basis <- matrix(NA, n, length(thresholds))
  for(j in 1:length(thresholds)) {
    basis[,j] = ifelse(a >= thresholds[j], 1, 0)
  }
  
  est.lc <- EstiamteLinearContrast(time = time, event = event,
                                   basis = basis, a = a, w = w, 
                                   weights = weights, t0 = t0, 
                                   fitted.surv = fitted.surv,
                                   fitted.cens = fitted.cens,
                                   cond.surv = cond.surv,
                                   cond.cens = cond.cens,
                                   fitted.haz.event = fitted.haz.event,
                                   fitted.haz.cens = fitted.haz.cens,
                                   prop.score = prop.score)
  
  test.stats <- numeric(length(t0))
  p.vals <- numeric(length(t0))
  for(s in 1:length(t0)) {
    test.stats[s] <- max(abs(est.lc$one.step.cent[,s]))
    Sigma <- cov(est.lc$influence.function.cent[,,s])
    mc <- rmvnorm(no.mc, mean = rep(0, length(thresholds)), sigma = Sigma/length(time))
    mc.test.stats <- apply(mc, 1, function(x) max(abs(x)))
    p.vals[s] <- mean(mc.test.stats > test.stats[s])
  }
  
  out <- list(test.stats = test.stats,
              p.vals = p.vals,
              t0 = t0)
  return(out)
}



# set.seed(1122)
# n <- 1000
# x <- matrix(runif(2 * n, 0, 1), nrow = n, ncol = 2)
# a <- x[,1]
# w <- x[,2]
# s.time <- rexp(n, rate = exp(x[,1] + x[,2] - 2))
# c.time <- rexp(n, rate = 1/2 * exp(x[,1] + x[,2] - 2))
# time <- apply(cbind(s.time, c.time), 1, min)
# event <- ifelse(s.time < c.time, 1, 0)
# weights <-rep(1,n) #sample(c(1:3), 200, replace= TRUE)
# t0 <- c(.5, 1, 1.5, 2)
# 
# Y <- cbind(event, time)
# colnames(Y) <- c("event", "time")
# ph.hal <- fit_hal(Y = Surv(time = time, event = event), X = x, family = "cox")
# ph.cox <- coxph(Surv(time = time, event = event) ~ x)
# jkl <- predict(asdf, new_data = x, type = "link")
# plot(predict(ph.hal, new_data = x, type = "link"), predict(ph.cox))
