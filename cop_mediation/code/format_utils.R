#' @param x numeric vector of point estimate, cil, ciu
#' @param digits number of digits to round to
format_ci <- function(x, digits = 3, ve = TRUE){
  paste0(sprintf(paste0("%.", digits, "f"), x[1]), " (",
         sprintf(paste0("%.", digits, "f"), x[ifelse(ve, 3, 2)]), ", ",
         sprintf(paste0("%.", digits, "f"), x[ifelse(ve, 2, 3)]), ")")
}

#' Format rows for a table as est (95% CI)
#' @param fit a natmed2 fit
format_row <- function(fit, digits = 3){
	# extract risk and effect estimates
  risk_est <- fit$risk[,1]
  total_eff <- fit$eff["Total", "one_step_est"]
  direct_eff <- fit$eff["Direct", "one_step_est"]
  indirect_eff <- fit$eff["Indirect", "one_step_est"]

  # compute proportion mediated
  prop_med <- 1 - log(direct_eff) / log(total_eff)
  g <- matrix(c(log(direct_eff)/(log(total_eff))^2*1/risk_est[1], log(indirect_eff)/(log(total_eff))^2*1/risk_est[2],
              -1/(risk_est[3]*log(total_eff)), 0))
  se_prop_med <- sqrt(t(g) %*% fit$cov %*% g)
  cil <- prop_med + qnorm(0.975) * se_prop_med
  ciu <- prop_med - qnorm(0.975) * se_prop_med

  this_row <- c(
    format_ci(1 - fit$eff["Direct", 2:4], digits = digits),
    format_ci(1 - fit$eff["Indirect", 2:4], digits = digits),
    format_ci(c(prop_med, cil, ciu), digits = digits)
  )
  return(this_row)
}



#' Output total, in/direct effects, CI's based on fitted results
#' @normal_survtmle_fit fitted object for estimating the total effect
#' @mediation_survtmle_fit fitted object for estimating mediation parameters
compute_mediation_params <- function(
  normal_survtmle_fit,
  mediation_survtmle_fit,
  ...
){
  # combine influence functions (ey00, ey11, ey10)
  all_ic <- cbind(normal_survtmle_fit$ic, mediation_survtmle_fit$ic)
  # estimates
  est <- c(normal_survtmle_fit$est[,1], mediation_survtmle_fit$est[,1])
  # covariance matrix
  cov_mat <- cov(all_ic) / dim(all_ic)[1]

  # confidence intervals of estimators
  est_cils <- est - 1.96*sqrt(diag(cov_mat))
  est_cius <- est + 1.96*sqrt(diag(cov_mat))
  out_est <- data.frame(est = est, cil = est_cils, ciu = est_cius)
  row.names(out_est) <- c("ey00", "ey11", "ey10")

  # delta method
  A <- matrix(c(
    -1 / est[1], 1 / est[2],     0      ,
    0       , 1 / est[2], -1 / est[3],
    -1 / est[1],      0    ,  1 / est[3]
  ), nrow = 3, ncol = 3, byrow = TRUE)

  log_total_eff <- log(est[2] / est[1])
  log_indirect_eff <- log(est[2] / est[3])
  log_direct_eff <- log(est[3] / est[1])
  log_effs <- c(log_total_eff, log_indirect_eff, log_direct_eff)
  ses_log_eff <- sqrt(diag(A %*% cov_mat %*% t(A)))

  cils <- exp(log_effs - 1.96 * ses_log_eff)
  cius <- exp(log_effs + 1.96 * ses_log_eff)
  out_eff <- data.frame(est = exp(log_effs), cil = cils, ciu = cius)

  # estimation of proportional mediated
  prop_med <- 1 - log_direct_eff/log_total_eff
  g <- matrix(c(log_indirect_eff/log_total_eff^2*1/est[1],
                log_direct_eff/log_total_eff^2*1/est[2],
               -1/(est[3]*log_total_eff)), ncol = 3, byrow = T)
  ses_prop_med <- sqrt(g %*% cov_mat %*% t(g))
  cil_prop <- prop_med - 1.96*ses_prop_med
  ciu_prop <- prop_med + 1.96*ses_prop_med

  out_eff <- rbind(out_eff, c(prop_med, cil_prop, ciu_prop))
  row.names(out_eff) <- c("Total", "Indirect", "Direct", "Prop_med")
  out <- list(risk = out_est,
              eff = out_eff)
  class(out) <- "survtmle_natmed"
  return(out)
}
