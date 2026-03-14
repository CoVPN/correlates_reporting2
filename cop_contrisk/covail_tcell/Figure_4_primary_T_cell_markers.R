dir.create("./Figure/naive", showWarnings = FALSE, recursive=TRUE) 
dir.create("./Figure/nonnaive", showWarnings = FALSE, recursive=TRUE) 
###############################################################################
# Load packages
library(dplyr)
library(vaccine)
library(ggplot2)
library(patchwork)
library(Hmisc)
###############################################################################
# Helper functions
###############################################################################
truncate_ce <- function(ests, low, high){
  ind_exclude = which(ests$cr$s < low | ests$cr$s > high)
  
  if(length(ind_exclude) == 0)
    return(ests)
  
  else {
    ests$cr$s = ests$cr$s[-ind_exclude]
    ests$cr$est = ests$cr$est[-ind_exclude]
    ests$cr$se = ests$cr$se[-ind_exclude]
    ests$cr$ci_lower = ests$cr$ci_lower[-ind_exclude]
    ests$cr$ci_upper = ests$cr$ci_upper[-ind_exclude]
    return(ests)
  }
}

# Helper function
# Conduct a basic CoR analysis
CoR_D1 <- function(dt, density,
                   arm_name, time, event, 
                   marker_name, adj_vars, t0,
                   title_name, xlab_name,
                   dir_name, x_range = c(0.1, 3),
                   x_break = c(0.1, 0.3, 1, 3, 10),
                   save_pdf = TRUE){
  
  dat <- load_data(
    time = time,
    event = event,
    vacc = arm_name,
    marker = marker_name,
    covariates = adj_vars,
    weights = "wt.D15.tcell",
    ph2 = "casecontrol",
    data = dt)
  
  ests_cox <- est_ce(dat = dat, type="Cox", t_0 = t0, return_p_value = TRUE)
  ests_np <- est_ce(dat = dat, type="NP", t_0 = t0, return_p_value = TRUE)
  
  low = wtd.quantile(density$s, weights = density$w, probs = c(0.025))
  high = wtd.quantile(density$s, weights = density$w, probs = c(0.975))
  
  Cox = truncate_ce(ests_cox, low, high)
  Nonparametric = truncate_ce(ests_np, low, high)
  
  p = plot_ce(Cox, Nonparametric,
              zoom_y = c(0, 1.1),zoom_x = c(log10(x_range[1])-0.3, log10(x_range[2])+0.3)) +
    geom_density(aes(x = x, y = after_stat(scaled), weight = w), 
                 data = data.frame(x = density$s, w = density$w),  
                 inherit.aes = F, bounds = c(min(density$s, na.rm = TRUE), 
                                             max(density$s, na.rm = TRUE)),
                 fill = "forestgreen", alpha = 0.15, color = NA) +
    geom_line(size = 1.5) +
    scale_x_continuous(breaks = log10(x_break), 
                       labels = function(x) paste0(10^x, "%")) +
    scale_y_continuous(labels = function(x) paste0(x*100, "%")) +
    xlab(xlab_name) + 
    ylab(paste0('Covariate-adjusted risk \n', t0, ' days post D15')) +
    labs(title = '') +
    theme_bw(base_size = 29) + 
    theme(legend.position = 'none')
  
  if (save_pdf)
    ggsave(dir_name, plot = p, width = 12, height = 8)
  
  return(list(plot = p,
              Cox_CoR = ests_cox,
              NP_CoR = ests_np))
}


# CoR analysis by baseline naive and non-naive status
CoR_D15_by_naive <- function(dt, arm_name, time, event, 
                             marker_name, mk_vec, 
                             adj_vars, t0,
                             title_name, xlab_name, dir_name,
                             x_range = c(0.1, 3),
                             x_break = c(0.1, 0.3, 1, 3, 10),
                             naive_only = FALSE,
                             vertical = TRUE){
  
  # Divide data into subgroups based on baseline marker level
  ind_naive = which(dt$naive == 1)
  ind_nonnaive = which(dt$naive == 0)
  
  dt_naive = dt[ind_naive,]
  dt_nonnaive = dt[ind_nonnaive, ]
  
  # Naive
  cat('Running naive', '\n')
  density_naive = list(s = mk_vec[ind_naive], w = dt$wt.D15.tcell[ind_naive])
  
  res_naive = CoR_D1(dt_naive, density_naive,
                   arm_name, time, event, 
                   marker_name, adj_vars, t0,
                   title_name = 'SARS-CoV-2 naive at baseline',
                   xlab_name,
                   dir_name = paste0('./Figure/naive/', dir_name), 
                   x_range = x_range,
                   x_break = x_break,
                   save_pdf = TRUE)
  
  p_naive = res_naive$plot
  
  
  # Non-naive
  if (!naive_only){
    cat('Running non-naive', '\n')
    density_nonnaive = list(s = mk_vec[ind_nonnaive], w = dt$wt.D15.tcell[ind_nonnaive])
    res_nonnaive = CoR_D1(dt_nonnaive, density_nonnaive,
                          arm_name, time, event, 
                          marker_name, adj_vars, t0,
                          title_name = 'SARS-CoV-2 non-naive at baseline',
                          xlab_name,
                          dir_name = paste0('./Figure/nonnaive/', dir_name),  
                          x_range = x_range,
                          x_break = x_break,
                          save_pdf = TRUE)
    
    p_nonnaive = res_nonnaive$plot
  }
}

# A wrap-up function for running analysis for all primary markers 

CoR_primary_analysis <- function(dt, mk_primary_list,
                                 mk_primary_names_list,
                                 arm_name = 'TrtonedosemRNA',
                                 title_name = 'Controlled risk among one-dose mRNA vaccinees',
                                 six_month_only = FALSE,
                                 naive_only = FALSE){
  for (i in 1:length(mk_primary_list)){
    var_name = mk_primary_list[i]
    mk_nm = mk_primary_names_list[i]
    cat('Running D22-D181', var_name, '\n')
    
    x_range = c(0.001, 10)
    x_break = c(0.001, 0.01, 0.1, 1, 10)
    
    # 6-Month
    cor_res = CoR_D15_by_naive(dt, arm_name = arm_name,
                               time = 'COVIDtimeD22toD181',
                               event = 'COVIDIndD22toD181',
                               marker_name = var_name,
                               mk_vec = pull(dt, var_name),
                               adj_vars = c("risk_score","FOIstandardized"), 
                               t0 = 188, 
                               title_name = title_name,
                               xlab_name = mk_nm,
                               dir_name = paste0('CoR_D22toD181_', var_name, '_', arm_name, '.pdf'),
                               x_range = x_range,
                               x_break = x_break,
                               naive_only = naive_only)
    
    if (!six_month_only){
      # 3-Month
      cor_res2 = CoR_D15_by_naive(dt, arm_name = arm_name,
                                  time = 'COVIDtimeD22toD181',
                                  event = 'COVIDIndD22toD181',
                                  marker_name = var_name,
                                  mk_vec = pull(dt, var_name),
                                  adj_vars = c("risk_score","FOIstandardized"), 
                                  t0 = 91, 
                                  title_name = title_name,
                                  xlab_name = mk_nm,
                                  dir_name = paste0('CoR_D22toD91_', var_name, '_', arm_name, '.pdf'),
                                  x_range = x_range,
                                  x_break = x_break,
                                  naive_only = naive_only)
    } 
  }
}

#########################################################################################
# Load the data
dt = read.csv('./covail_data_processed_20250612.csv')

# Define the PP correlates cohort
dt_PP = dt %>%
  filter(ph1.D15.tcell == 1) %>% # n = 932
  mutate(risk_score = ifelse(is.na(risk_score), median(risk_score), risk_score)) %>%
  mutate(risk_score = scale(risk_score)) %>%
  mutate(casecontrol = TwophasesampIndD15.tcell) # all are included in the phase 2 sampling


################################################################################
# Primary markers
mk_primary_list = c('Day15cd4_IFNg.IL2_BA.4.5.S', 
                    'Bcd4_IFNg.IL2_BA.4.5.S',
                    'Bcd4_IFNg.IL2_Wuhan.N')
mk_primary_names_list = c('D15 % CD4+ T cell expressing IFNg \nand/or IL-2 against Spike BA.4/5',
                          'D1 % CD4+ T cell expressing IFNg \nand/or IL-2 against Spike BA.4/5',
                          'D1 % CD4+ T cell expressing IFNg \nand/or IL-2 against Index N')

###############################################################################
###############################################################################
###############################################################################
# All one dose mRNA arms pooled
dt_PP_onedosemRNA = dt_PP %>%
  filter(TrtonedosemRNA == 1)

CoR_primary_analysis(dt = dt_PP_onedosemRNA,
                     mk_primary_list = mk_primary_list,
                     mk_primary_names_list = mk_primary_names_list,
                     title_name = 'Controlled risk among one-dose mRNA vaccinees',
                     arm_name = 'TrtonedosemRNA')

