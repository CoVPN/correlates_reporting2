#!/usr/local/bin/Rscript

library(CFsurvival)
library(lubridate)
library(knitr)
library(sas7bdat)
library(survival)
library(MASS)
library(rms)
library(haven)
library(epitools)
library(ranger)
# 
# library(mgcv)
# detach("package:mgcv", unload = TRUE)

expit <- function(x) exp(x)/(1 + exp(x))

CumulativeIncidenceFit <- function(time, event, W, A, t0 = quantile(time[event == 1], .90), np = FALSE,
                                   fit.times = seq(0, t0, 7)) {
  # time -- time of endpoint or censoring (i.e., min(T, C))
  # event -- binary indicator of whether an event was observed or not
  # W -- n x p matrix of covariates
  # A -- n-dimensional vector of binary treatment variables
  # np -- boolean indicating whether to use nonparametric estimators for the
  #       nuisance parameters (uses Cox regression & logistic regresion when FALSE)
  # fit.times -- sequence of time points at which the cumulative risk should be estimated
  
  if(max(fit.times) != t0) fit.times <- c(fit.times, t0)
  
  fit <- NULL
  if(np == TRUE) {
    # Estimate the propensity score
    fit <- CFsurvival(time = time, event = event, treat = A, confounders = W,
                      fit.times = fit.times, contrasts = c("surv.diff", "risk.ratio"),
                      nuisance.options = list(event.SL.library = c("survSL.km", "survSL.coxph",
                                                                   "survSL.gam", "survSL.rfsrc"),
                                              cens.SL.library = c("survSL.km", "survSL.coxph",
                                                                  "survSL.gam", "survSL.rfsrc"),
                                              prop.SL.library = c("SL.mean", "SL.glm", 
                                                                  "SL.gam", "SL.ranger"),
                                              cross.fit = TRUE))
  } else {
    # if W == Null, perform no covariate adjustment
    if(is.null(W)) {
      W <- matrix(data = runif(2 * length(time)), nrow = length(time), ncol = 2)
      prop.fit <- rep(mean(A), length(time))
      fit <- CFsurvival(time = time, event = event, treat = A, confounders = W,
                        fit.times = fit.times, contrasts = c("surv.diff", "risk.ratio"),
                        nuisance.options = list(event.SL.library = "survSL.km",
                                                cens.SL.library = "survSL.km",
                                                prop.SL.library = NULL,
                                                prop.pred = prop.fit,
                                                cross.fit = FALSE))
    } else {
      glm.fit <- glm(A ~ W, family = binomial)
      prop.fit <- glm.fit$fitted.values
      fit <- CFsurvival(time = time, event = event, treat = A, confounders = W,
                        fit.times = fit.times, contrasts = c("surv.diff", "risk.ratio"),
                        nuisance.options = list(event.SL.library = "survSL.coxph",
                                                cens.SL.library = "survSL.coxph",
                                                prop.SL.library = NULL,
                                                prop.pred = prop.fit,
                                                cross.fit = FALSE))
    }
  }
  return(fit)
}

CumulativeIncidenceSummary <- function(time, event, W, A, t0, t.start = 1, fit,
                                       t.start.ci = 1,
                                       y.axis.bound.risk = NULL,
                                       y.axis.bound.riskdiff = NULL,
                                       y.axis.bound.riskratio = NULL,
                                       label = "",
                                       path.prefix = paste0(getwd(), "/"),
                                       path.suffix,
                                       time.axis.label = "time",
                                       tx0.label = "Vaccine",
                                       tx1.label = "Hybrid",
                                       tx0.label2 = "Vaccine",
                                       tx1.label2 = "Hybrid",
                                       override.t0 = FALSE,
                                       add.ticks = FALSE,
                                       time.axis.offset = 0) {
  
  ### Modify t0 if needed
  t0.index.rr <- which(fit$risk.ratio.df$time == t0)
  se.t0 <- fit$risk.ratio.df$se.log.ratio[t0.index.rr]
  if(is.na(se.t0) | is.nan(se.t0)) se.t0 <- Inf
  flag <- abs(log(3)/se.t0[1]) >= abs(qnorm(.025, mean = 0, sd = 1))
  if(!flag & !override.t0) {
    t0.new <- max(fit$risk.ratio.df$time[log(3)/fit$risk.ratio.df$se.log.ratio >= abs(qnorm(.025, mean = 0, sd = 1))],
                  na.rm = TRUE)
    if(!is.na(t0.new) & !is.nan(t0.new) & t0.new > 0) t0 <- t0.new
  }
  
  ### Summarize results for forest plots
  
  ## cumulative incidence in "treatment" group
  incidence.ctrl <- matrix(NA, nrow = 1, ncol = 3)
  colnames(incidence.ctrl) <- c("CI", "Low", "High")
  t0.index0 <- which(fit$surv.df$time == t0 & fit$surv.df$trt == 0)
  incidence.ctrl[1,1:3] <- c(1 - fit$surv.df[t0.index0, 3],
                             1 - fit$surv.df[t0.index0, 7],
                             1 - fit$surv.df[t0.index0, 6])
  
  ## cumulative incidence in "control" group
  incidence.tx <- matrix(NA, nrow = 1, ncol = 3)
  colnames(incidence.tx) <- c("CI", "Low", "High")
  t0.index1 <- which(fit$surv.df$time == t0 & fit$surv.df$trt == 1)
  incidence.tx[1,1:3] <- c(1 - fit$surv.df[t0.index1, 3],
                           1 - fit$surv.df[t0.index1, 7],
                           1 - fit$surv.df[t0.index1, 6])
  
  ## risk ratio
  rr.tab <- matrix(NA, nrow = 1, ncol = 4)
  colnames(rr.tab) <- c("CI", "Low", "High", "P-value")
  t0.index.rr <- which(fit$risk.ratio.df$time == t0)
  rr.tab[1,1:4] <- c(fit$risk.ratio.df[t0.index.rr, 3],
                     fit$risk.ratio.df[t0.index.rr, 7],
                     fit$risk.ratio.df[t0.index.rr,8],
                     pchisq(fit$risk.ratio.df$log.risk.ratio[t0.index.rr]^2/
                              fit$risk.ratio.df$se.log.ratio[t0.index.rr]^2,
                            df = 1, lower.tail = F))
  
  ## risk difference
  rd.tab <- matrix(NA, nrow = 1, ncol = 4)
  colnames(rd.tab) <- c("CI", "Low", "High", "P-value")
  t0.index.rd <- which(fit$surv.diff.df$time == t0)
  rd.tab[1,1:4] <- c(-fit$surv.diff.df[t0.index.rd, 2],
                     -fit$surv.diff.df[t0.index.rd, 6],
                     -fit$surv.diff.df[t0.index.rd, 5],
                     pchisq(fit$surv.diff.df$surv.diff[t0.index.rd]^2/
                              fit$surv.diff.df$se[t0.index.rd]^2, 
                            df = 1, lower.tail = F))
  
  
  t.keep <- seq(t.start, t0, 1)
  
  ### Plot cumulative incidence curves
  ## Produce survival curves
  fit$surv.df$unif.lower <- fit$surv.df$unif.logit.lower
  fit$surv.df$unif.upper <- fit$surv.df$unif.logit.upper
  if(is.null(y.axis.bound.risk)) {
    y.axis.bound <- max(abs(c(1 - min(fit$surv.df$unif.lower[fit$surv.df$time <= t0], na.rm = T), 
                              1 - max(fit$surv.df$unif.upper[fit$surv.df$time <= t0], na.rm = T))))
    y.axis.bound <- max(c(y.axis.bound, .05))
  } else {
    y.axis.bound <- y.axis.bound.risk
  }
  
  # par(las=1, cex.axis=1.1, cex.lab=1.35, tcl= -0.35, mar=c(1, 1, 1, 1), 
  #     oma = c(1, 1, 1, 1),  mgp = c(2.3, 0.5, 0)
  
  pdf(file = paste0(path.prefix, "incidence_", path.suffix, ".pdf"))
  par(las=1, cex.axis=.8, cex.lab=1.35, tcl= -0.35, mar=c(8, 6.25, 6, 0.5),
      oma = c(2, 1, 1, 1),  mgp = c(2.3, 0.5, 0))
  
  plot(c(0, t0), c(0, 1), 
       lty = 1, col = "orange", lwd = 2.5, type = "n", 
       ylim = c(0, y.axis.bound),
       main =  label, 
       xlab = time.axis.label, ylab = "Cumulative Incidence",
       xaxt = "n")
  
  
  abline(v = t0, lwd = 1, lty = 3, col = "grey")
  lines(subset(fit$surv.df, trt == 0)$time[subset(fit$surv.df, trt == 0)$time <= t0],
        1 - subset(fit$surv.df, trt == 0)$surv[subset(fit$surv.df, trt == 0)$time <= t0],
        type = "l", lty = 1, col = "#F2502A", lwd = 2.5)
  
  lines(subset(fit$surv.df, trt == 0)$time[subset(fit$surv.df, trt == 0)$time <= t0 & subset(fit$surv.df, trt == 0)$time >=  t.start.ci],
        1 - subset(fit$surv.df, trt == 0)$ptwise.logit.upper[subset(fit$surv.df, trt == 0)$time <= t0 & subset(fit$surv.df, trt == 0)$time >= t.start.ci],
        type = "l", lty = 2, col = "#F2502A", lwd = 1.5)
  lines(subset(fit$surv.df, trt == 0)$time[subset(fit$surv.df, trt == 0)$time <= t0 & subset(fit$surv.df, trt == 0)$time >=  t.start.ci],
        1 - subset(fit$surv.df, trt == 0)$ptwise.logit.lower[subset(fit$surv.df, trt == 0)$time <= t0 & subset(fit$surv.df, trt == 0)$time >= t.start.ci],
        type = "l", lty = 2, col = "#F2502A", lwd = 1.5)
  
  # lines(subset(fit$surv.df, trt == 0)$time[subset(fit$surv.df, trt == 0)$time <= t0 & subset(fit$surv.df, trt == 0)$time >= 14],
  #       1 - subset(fit$surv.df, trt == 0)$unif.upper[subset(fit$surv.df, trt == 0)$time <= t0 & subset(fit$surv.df, trt == 0)$time >= 14],
  #       type = "l", lty = 3, col = "#F2502A", lwd = 1.1)
  # lines(subset(fit$surv.df, trt == 0)$time[subset(fit$surv.df, trt == 0)$time <= t0 & subset(fit$surv.df, trt == 0)$time >= 14],
  #       1 - subset(fit$surv.df, trt == 0)$unif.lower[subset(fit$surv.df, trt == 0)$time <= t0 & subset(fit$surv.df, trt == 0)$time >= 14],
  #       type = "l", lty = 3, col = "#F2502A", lwd = 1.1)
  
  abline(v = t0, lwd = 1, lty = 3, col = "grey")
  lines(subset(fit$surv.df, trt == 1)$time[subset(fit$surv.df, trt == 1)$time <= t0], 
        1 - subset(fit$surv.df, trt == 1)$surv[subset(fit$surv.df, trt == 1)$time <= t0],
        type = "l", lty = 1, col = "blue", lwd = 2.5)
  
  lines(subset(fit$surv.df, trt == 1)$time[subset(fit$surv.df, trt == 1)$time <= t0 & subset(fit$surv.df, trt == 1)$time >=  t.start.ci],
        1 - subset(fit$surv.df, trt == 1)$ptwise.logit.upper[subset(fit$surv.df, trt == 1)$time <= t0 & subset(fit$surv.df, trt == 1)$time >= t.start.ci],
        type = "l", lty = 2, col = "blue", lwd = 1.5)
  lines(subset(fit$surv.df, trt == 1)$time[subset(fit$surv.df, trt == 1)$time <= t0 & subset(fit$surv.df, trt == 1)$time >=  t.start.ci],
        1 - subset(fit$surv.df, trt == 1)$ptwise.logit.lower[subset(fit$surv.df, trt == 1)$time <= t0 & subset(fit$surv.df, trt == 1)$time >= t.start.ci],
        type = "l", lty = 2, col = "blue", lwd = 1.5)
  
  # lines(subset(fit$surv.df, trt == 1)$time[subset(fit$surv.df, trt == 1)$time <= t0 & subset(fit$surv.df, trt == 1)$time >= 14],
  #       1 - subset(fit$surv.df, trt == 1)$unif.upper[subset(fit$surv.df, trt == 1)$time <= t0 & subset(fit$surv.df, trt == 1)$time >= 14],
  #       type = "l", lty = 3, col = "blue", lwd = 1.1)
  # lines(subset(fit$surv.df, trt == 1)$time[subset(fit$surv.df, trt == 1)$time <= t0 & subset(fit$surv.df, trt == 1)$time >= 14], 
  #       1 - subset(fit$surv.df, trt == 1)$unif.lower[subset(fit$surv.df, trt == 1)$time <= t0 & subset(fit$surv.df, trt == 1)$time >= 14],
  #       type = "l", lty = 3, col = "blue", lwd = 1.1)
  
  # Plot tick marks for censoring times
  if(add.ticks) {
    cens <- 1 - event
    cens0.fit <- survfit(Surv(time[A == 0], cens[A == 0]) ~ 1)
    cens1.fit <- survfit(Surv(time[A == 1], cens[A == 1]) ~ 1)
    
    cens0.time <- cens0.fit$time[cens0.fit$time <= t0]
    cens1.time <- cens1.fit$time[cens1.fit$time <= t0]
    
    cens0.risk <- numeric(length(cens0.time))
    cens1.risk <- numeric(length(cens1.time))
    
    for(i in 1:length(cens0.risk)) {
      j <- which(fit$surv.df$time == cens0.time[i] & fit$surv.df$trt == 0)
      cens0.risk[i] <- 1 - fit$surv.df$surv[j]
    }
    
    for(i in 1:length(cens1.risk)) {
      j <- which(fit$surv.df$time == cens1.time[i] & fit$surv.df$trt == 1)
      cens1.risk[i] <- 1 - fit$surv.df$surv[j]
    }
    
    # cens0.risk <- 1 - summary(fit, time = cens0.time)$surv[1:length(cens0.time)]
    # cens1.risk <- 1 - summary(fit, time = cens1.time)$surv[(length(cens1.time)  + 1):(2 * length(cens1.time))]
    
    segments(x0 = cens0.time, x1 = cens0.time, 
             y0 = cens0.risk - y.axis.bound/60,
             y1 = cens0.risk + y.axis.bound/60, lwd = 1, col = "#F2502A")
    
    segments(x0 = cens1.time, x1 = cens1.time, 
             y0 = cens1.risk - y.axis.bound/60,
             y1 = cens1.risk + y.axis.bound/60, lwd = 1, col = "blue")
  }
  
  legend("topleft", col = c("#F2502A", "blue"), lwd = c(2,2), lty = c(1,1),
         legend = c(tx0.label, tx1.label), cex = .8, bg = "white")
  
  ## Add event count and no. at risk below the margins
  x.time <- floor(seq(1, t0, length.out = 5))  #c(1, 30, 60, 90, 120, 150, 180)
  x.time.lab <- x.time + time.axis.offset
  # x.time <- seq(1, t0, 30)
  # axis(side = 1, at = c(0, x.time[-1]), cex = .25)
  axis(side = 1, at = x.time, cex = .25, labels = x.time.lab)
  survfit.A0 <- survfit(Surv(time[A == 0], event[A == 0]) ~ 1, type = "kaplan-meier")
  survfit.A1 <- survfit(Surv(time[A == 1], event[A == 1]) ~ 1, type = "kaplan-meier")
  
  n.risk.A0    <-  summary(survfit.A0, times=x.time, extend=TRUE)$n.risk
  n.risk.A1    <-  summary(survfit.A1, times=x.time, extend=TRUE)$n.risk
  
  n.event.A0    <-  cumsum(summary(survfit.A0, times=x.time, extend=TRUE)$n.event)
  n.event.A1    <-  cumsum(summary(survfit.A1, times=x.time, extend=TRUE)$n.event)
  
  usrX <- par("usr")[1:2]
  widthUser      <- diff(usrX)
  widthInches    <- par("pin")[1]
  leftMargInches <- par("mai")[2]
  
  ## Compute the 'user-coordinate' value that corresponds to the start of the
  ## *figure* region (par("usr") gives us user coordinates of the plot region).
  atLoc <- usrX[1] - (leftMargInches - 0) * (widthUser/widthInches) 
  
  mtext( expression(bold("No. at Risk")), side=1, line=3.4,
         at=atLoc, adj=0, font=2, cex=.85)
  mtext(tx0.label2, side=1, line=4.25, at=atLoc, adj=0, cex=.7)
  mtext(tx1.label2, side=1, line=4.95,    at=atLoc, adj=0, cex=.7)
  mtext("Total Events", side=1, line=6, at=atLoc, adj=0, font=2, cex=.85)
  mtext(tx0.label2, side=1, line=6.8, at=atLoc, adj=0, cex=.7)
  mtext(tx1.label2, side=1, line=7.5,    at=atLoc, adj=0, cex=.7)
  
  mtext(n.risk.A0, side=1, line=4.25, at=x.time, cex=.7)
  mtext(n.risk.A1, side=1, line=4.95,    at=x.time, cex=.7)
  
  mtext(n.event.A0, side=1, line=6.8, at=x.time, cex=.7)
  mtext(n.event.A1, side=1, line=7.5,    at=x.time, cex=.7)
  
  # dev.off()
  
  
  ### plot risk difference
  if(is.null(y.axis.bound.riskdiff)) {
    y.axis.bound <- max(abs(c(min(fit$surv.diff.df$unif.lower[fit$surv.diff.df$time <= t0], na.rm = T), 
                              max(fit$surv.diff.df$unif.upper[fit$surv.diff.df$time <= t0], na.rm = T))))
    y.axis.bound <- min(c(y.axis.bound, .25))
  } else {
    y.axis.bound <- y.axis.bound.riskdiff
  }
  
  # pdf(file = paste0(path.prefix, "riskdiff_", path.suffix, ".pdf"))
  plot(fit$surv.diff.df$time[fit$surv.diff.df$time <= t0], 
       -fit$surv.diff.df$surv.diff[fit$surv.diff.df$time <= t0], type = "l", 
       lty = 1, col = "orange", lwd = 2.5, 
       ylim = c(-y.axis.bound, y.axis.bound),
       main = label, 
       xlab = time.axis.label, ylab = "Risk Difference", xaxt = "n")
  x.time <- floor(seq(1, t0, length.out = 5))  #c(1, 30, 60, 90, 120, 150, 180)
  x.time.lab <- x.time + time.axis.offset
  # x.time <- seq(1, t0, 30)
  # axis(side = 1, at = c(0, x.time[-1]), cex = .25)
  axis(side = 1, at = x.time, cex = .25, labels = x.time.lab)
  
  abline(h = 0, lwd = 1, lty = 2, col = "grey")
  abline(v = t0, lwd = 1, lty = 3, col = "grey")
  lines(fit$surv.diff.df$time[fit$surv.diff.df$time <= t0 & fit$surv.diff.df$time >=  t.start.ci],
        -fit$surv.diff.df$ptwise.lower[fit$surv.diff.df$time <= t0 & fit$surv.diff.df$time >= t.start.ci], type = "l",
        lty = 2, col = "orange", lwd = 1.5)
  lines(fit$surv.diff.df$time[fit$surv.diff.df$time <= t0 & fit$surv.diff.df$time >=  t.start.ci],
        -fit$surv.diff.df$ptwise.upper[fit$surv.diff.df$time <= t0 & fit$surv.diff.df$time >= t.start.ci], type = "l",
        lty = 2, col = "orange", lwd = 1.5)
  # lines(fit$surv.diff.df$time[fit$surv.diff.df$time <= t0 & fit$surv.diff.df$time >= 14],
  #       -fit$surv.diff.df$unif.lower[fit$surv.diff.df$time <= t0 & fit$surv.diff.df$time >= 14], type = "l",
  #       lty = 3, col = "orange", lwd = 1.1)
  # lines(fit$surv.diff.df$time[fit$surv.diff.df$time <= t0 & fit$surv.diff.df$time >= 14],
  #       -fit$surv.diff.df$unif.upper[fit$surv.diff.df$time <= t0 & fit$surv.diff.df$time >= 14], type = "l",
  #       lty = 3, col = "orange", lwd = 1.1)
  
  # dev.off()
  
  ### plot risk ratio
  # pdf(file = paste0(path.prefix, "riskratio_", path.suffix, ".pdf"))
  
  if(is.null(y.axis.bound.riskratio)) {
    y.axis.bound <- max(abs(c(min(log2(fit$risk.ratio.df$unif.lower)[fit$risk.ratio.df$time <= t0], na.rm = T), 
                              max(log2(fit$risk.ratio.df$unif.upper)[fit$risk.ratio.df$time <= t0], na.rm = T))))
    y.axis.bound <- min(c(y.axis.bound, log2(16)))
  } else {
    y.axis.bound <- y.axis.bound.riskratio
  }
  
  plot(fit$risk.ratio.df$time[fit$risk.ratio.df$time <= t0],
       log2(fit$risk.ratio.df$risk.ratio)[fit$risk.ratio.df$time <= t0], type = "l", 
       lty = 1, col = "orange", lwd = 2.5, 
       ylim = c(-y.axis.bound, y.axis.bound),
       main = label, 
       xlab = time.axis.label, ylab = "Risk Ratio", yaxt = "n", xaxt = "n")
  x.time <- floor(seq(1, t0, length.out = 5))  #c(1, 30, 60, 90, 120, 150, 180)
  x.time.lab <- x.time + time.axis.offset
  # x.time <- seq(1, t0, 30)
  # axis(side = 1, at = c(0, x.time[-1]), cex = .25)
  axis(side = 1, at = x.time, cex = .25, labels = x.time.lab)
  
  step <- round(y.axis.bound/2)
  axis(side = 2,
       at = c(-step * 2, -step, 0, step, step * 2),
       labels = fractions(2^c(-step * 2, -step, 0, step, step * 2)),
       cex = .75)
  abline(h = 0, lwd = 1, lty = 2, col = "grey")
  abline(v = t0, lwd = 1, lty = 3, col = "grey")
  lines(fit$risk.ratio.df$time[fit$risk.ratio.df$time <= t0 & fit$risk.ratio.df$time >=  t.start.ci],
        log2(fit$risk.ratio.df$ptwise.lower)[fit$risk.ratio.df$time <= t0 & fit$risk.ratio.df$time >= t.start.ci], type = "l",
        lty = 2, col = "orange", lwd = 1.5)
  lines(fit$risk.ratio.df$time[fit$risk.ratio.df$time <= t0 & fit$risk.ratio.df$time >=  t.start.ci],
        log2(fit$risk.ratio.df$ptwise.upper)[fit$risk.ratio.df$time <= t0 & fit$risk.ratio.df$time >= t.start.ci], type = "l",
        lty = 2, col = "orange", lwd = 1.5)
  # lines(fit$risk.ratio.df$time[fit$risk.ratio.df$time <= t0 & fit$risk.ratio.df$time >= 14],
  #       log2(fit$risk.ratio.df$unif.lower[fit$risk.ratio.df$time <= t0 & fit$risk.ratio.df$time >= 14]), type = "l",
  #       lty = 3, col = "orange", lwd = 1.1)
  # lines(fit$risk.ratio.df$time[fit$risk.ratio.df$time <= t0 & fit$risk.ratio.df$time >= 14],
  #       log2(fit$risk.ratio.df$unif.upper[fit$risk.ratio.df$time <= t0 & fit$risk.ratio.df$time >= 14]), type = "l",
  #       lty = 3, col = "orange", lwd = 1.1)
  
  dev.off()
  
  out <- list(incidence.ctrl = incidence.ctrl,
              incidence.tx = incidence.tx,
              rr.tab = rr.tab,
              rd.tab = rd.tab,
              t0 = t0)
  return(out)
}
