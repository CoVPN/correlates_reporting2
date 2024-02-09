#Sys.setenv(TRIAL = "hvtn705second"); COR="D210"; Sys.setenv(VERBOSE = 1) 
#Sys.setenv(TRIAL = "moderna_boost"); COR="BD29naive"; Sys.setenv(VERBOSE = 1) 

renv::activate(project = here::here(".."))

source(here::here("..", "_common.R"))

source(here::here("code", "params.R"))

library(npthreshold)

library(sl3)
library(SuperLearner)
library(data.table)
library(mvtnorm)
library(uuid)
library(doMC)
library(earth)

#source(here::here("code", "tmleThresh.R"))
source(here::here("code", "learners.R"))
source(here::here("code", "survivalThresh", "Threshold_survivalCR.R"))
source(here::here("code", "survivalThresh", "fitting_likelihood.R"))
source(here::here("code", "survivalThresh", "Pooled_hazards.R"))
source(here::here("code", "survivalThresh", "sl3_learner_helpers.R"))
source(here::here("code", "survivalThresh", "survival_helper_functions.R"))
source(here::here("code", "survivalThresh", "targeting_functions.R"))
source(here::here("code", "survivalThresh", "task_generators.R"))

begin=Sys.time()


################################################################################
# code from clean_data.R

# Generate the outcome and censoring indicator variables

short_key <- COR

data <- dat.mock
data$Ttilde <- data$EventTimePrimary
data$Delta <- data$EventIndPrimary
data$TwophasesampInd <- data$ph2

##### Discretize time grid

max_t <- max(data[data$EventIndPrimary==1 & data$Trt == 1 & data$ph2 == 1, "EventTimePrimary" ])
size_time_grid <- 15 

time_grid <- unique(sort(c(max_t ,quantile(data$Ttilde[data$Ttilde <= max_t + 5 & data$TwophasesampInd ==1 & !is.na(data$Delta)& data$Delta==1 ], seq(0,1, length = size_time_grid)))))
time_grid[which.min(abs(time_grid -max_t))[1]] <- max_t
time_grid <- sort(unique(time_grid))

Ttilde_discrete <- findInterval(data$Ttilde, time_grid, all.inside = TRUE)
target_time <- findInterval(max_t, time_grid, all.inside = TRUE)
data$Ttilde <- Ttilde_discrete
data$target_time <- target_time
data$J <- 1


####
# Subset data

print(markers)
variables_to_keep <- c(covariates, markers, "TwophasesampInd", "wt", "Ttilde", "Delta", "target_time", "J")
variables_to_keep <- intersect(variables_to_keep, colnames(data))


if (TRIAL=="moderna_boost") {
  if (COR=='BD29naive') {
    keep = data[["ph1"]]==1 & data[["naive"]]==1  
  } else if (COR=='BD29nnaive') {
    keep = data[["ph1"]]==1 & data[["naive"]]==0
  } else stop('wrong COR: '%.% COR)
  
} else {
  keep = data[["ph1"]]==1 & data$Trt == 1 
}


data_firststage <- data[keep, variables_to_keep]

#data_firststage <- na.omit(data_firststage)
data_secondstage <- data_firststage[data_firststage$TwophasesampInd == 1, ]

write.csv(data_firststage,
          here::here('output', TRIAL, COR, "data_clean", paste0("data_firststage.csv")),
          row.names = F
)
write.csv(data_secondstage,
          here::here('output', TRIAL, COR, "data_clean", paste0("data_secondstage.csv")),
          row.names = F
)

markers.cpy=markers
markers=c()
for (marker in markers.cpy) {
  
  if (length(unique(data_secondstage[[marker]])) > threshold_grid_size) {
    # quantiles from cases
    thresh_grid <- report.assay.values(data_secondstage[[marker]][data_secondstage[["Delta"]]==1], marker, grid_size=threshold_grid_size)
    # minimum from both cases and controls
    min = min(data_secondstage[[marker]], na.rm=T)
    # no need to sort b/c report.assay.values sorts, but need unique
    thresh_grid = unique(c(min, thresh_grid))
    
  } else {
    thresh_grid <- sort(unique(data_secondstage[[marker]]))
  }
  
  # limit to ncases>=5
  ncases=sapply(thresh_grid, function(s) sum(data_secondstage[[marker]]>s & data_secondstage[["Delta"]]==1, na.rm=T))
  thresh_grid = thresh_grid[ncases>=5]
  
  if (length(thresh_grid)>1) {
    markers=c(markers, marker)
    write.csv(data.frame(thresh = thresh_grid), here::here('output', TRIAL, COR, "data_clean", "Thresholds_by_marker", paste0("thresholds_", marker, ".csv")), row.names = F)
  }
}



################################################################################
# Runs threshold analysis and saves raw results.

#lrnr <- get_learner(fast_analysis = fast_analysis, include_interactions = include_interactions, covariate_adjusted = covariate_adjusted)

#lrnr <- Lrnr_glm$new()
#lrnr_N <-lrnr
#lrnr_C <-lrnr
#lrnr_A <-lrnr

#' @param marker Marker to run threshold analysis for
run_threshold_analysis <- function(marker, direction = "above") {
  # marker=markers[1]; direction = "above"
  

  if("risk_score" %in% covariates) {
    form <- paste0("Y~h(.) + h(t,.) + h(A, .)")
  } else {
    form <- paste0("Y~h(.) + h(.,.)")
    
  }
  form <- paste0("Y~h(.) + h(t,.)  + h(A,.)")
  print(form)
  require(doMC)
  registerDoMC()
  ngrid_A <- 30
  lrnr_N <- Lrnr_hal9001_custom$new(max_degree = 2, smoothness_orders = 0, num_knots = c(10,3), reduce_basis =1e-3, formula_hal = form, fit_control = list(n_folds = 10, parallel = TRUE))
  lrnr_C <- Lrnr_hal9001_custom$new(max_degree = 2, smoothness_orders = 0, num_knots = c(5,1), reduce_basis=1e-3,  formula_hal = form, fit_control = list(n_folds = 10, parallel = TRUE))
  lrnr_A <- Lrnr_hal9001_custom$new(max_degree = 2, smoothness_orders = 0, num_knots = c(5,1), reduce_basis=1e-3, fit_control = list(n_folds = 10, parallel = TRUE))
  if(fast_analysis) {
    ngrid_A <- 15
    print("fast analysis")
    form <- paste0("Y~h(.) + h(t,.)")
    lrnr_N <- Lrnr_hal9001_custom$new(max_degree = 2, smoothness_orders = 0, num_knots = c(10,2), reduce_basis =1e-4, formula_hal = form, fit_control = list(n_folds = 10, parallel = TRUE))
    lrnr_C <- Lrnr_hal9001_custom$new(max_degree = 2, smoothness_orders = 0, num_knots = c(9,2), reduce_basis=1e-4,  formula_hal = form, fit_control = list(n_folds = 10, parallel = TRUE))
    lrnr_A <- Lrnr_hal9001_custom$new(max_degree = 2, smoothness_orders = 0, num_knots = c(10,3), reduce_basis=1e-4, fit_control = list(n_folds = 10, parallel = TRUE))
  }
  if(super_fast_analysis) {
    ngrid_A <- 15
    print("super fast analysis")
    form <- paste0("Y~h(.) + h(t,.)")
    lrnr_N <- Lrnr_hal9001_custom$new(max_degree = 2, smoothness_orders = 0, num_knots = c(3,1), reduce_basis =1e-4, formula_hal = form, fit_control = list(n_folds = 10, parallel = TRUE))
    lrnr_C <- Lrnr_hal9001_custom$new(max_degree = 2, smoothness_orders = 0, num_knots = c(3,1), reduce_basis=1e-4,  formula_hal = form, fit_control = list(n_folds = 10, parallel = TRUE))
    lrnr_A <- Lrnr_hal9001_custom$new(max_degree = 1, smoothness_orders = 0, num_knots = c(3,1), reduce_basis=1e-4, fit_control = list(n_folds = 10, parallel = TRUE))
  }
  
  thresholds <- read.csv(here::here('output', TRIAL, COR, "data_clean", "Thresholds_by_marker", paste0("thresholds_", marker, ".csv")))
  thresholds <- as.vector(unlist(thresholds[, 1])) 
  print(thresholds)
  #time <- marker_to_time[[marker]]
  data_full <- read.csv(here::here('output', TRIAL, COR, "data_clean", paste0("data_firststage", ".csv")))
  
  
  ####################################################
  # Run thresholdTMLE
  ####################################################
  print("Run thresholdTMLE")
  if(direction=="above") {
    #esttmle_full <- survivalThresh(as.data.table(data_full), trt = marker, Ttilde = "Ttilde",Delta = "Delta", J = "J", covariates = covariates, target_times = unique(data_full$target_time),cutoffs_A = thresholds, cutoffs_J = 1, type_J = "equal", lrnr =lrnr, lrnr_A = lrnr_A, lrnr_N = lrnr_N, lrnr_C = lrnr_C, biased_sampling_group= NULL, biased_sampling_indicator = "TwophasesampInd", weights_var = "wt", monotone_decreasing = T, ngrid_A = ngrid_A )
    node_list <- list("W" = covariates, "A" = marker, "Y" = "Delta", weights = "wt")
    lrnr <- Lrnr_sl$new(learners = Stack$new(
      Lrnr_glmnet$new(),
      Lrnr_gam$new(),
      Lrnr_earth$new(),
      Lrnr_mean$new()
    ), metalearner = Lrnr_cv_selector$new(loss_loglik_binomial))
    
    esttmle_full <- thresholdTMLE(data_full, node_list, thresholds = thresholds, biased_sampling_strata = NULL
                                  , biased_sampling_indicator = "TwophasesampInd"
                                  , lrnr_A = lrnr
                                  , lrnr_Y = lrnr
                                  , lrnr_Delta = Lrnr_glmnet$new()
                                  , monotone_decreasing = TRUE) 
    
    
  } else if(direction=="below") {
    stop("NO longer used")
    thresholds <- thresholds[-1]
    data_full[[marker]] <- -data_full[[marker]]
    esttmle_full <- survivalThresh(as.data.table(data_full), trt = marker, Ttilde = "Ttilde",Delta = "Delta", J = "J", covariates = covariates, target_times = unique(data_full$target_time),cutoffs_A = sort(-thresholds), cutoffs_J = 1, type_J = "equal", lrnr =lrnr, lrnr_A = lrnr_A, lrnr_N = lrnr_N, lrnr_C = lrnr_C, biased_sampling_group= NULL, biased_sampling_indicator = "TwophasesampInd", weights_var = "wt", monotone_decreasing = F, ngrid_A = min(2,ngrid_A) )
    # non-monotone
    if (!config$threshold_monotone_only) {
      esttmle_full[[1]][,1] <- - esttmle_full[[1]][,1]
      esttmle_full[[1]] <- esttmle_full[[1]][order(esttmle_full[[1]][,1]),]
    }
    esttmle_full[[2]][,1] <- - esttmle_full[[2]][,1]
    esttmle_full[[2]] <- esttmle_full[[2]][order(esttmle_full[[2]][,1]),]
    
    
  }
  
  direction_append <- ""
  if(direction=="below") {
    direction_append <- "_below"
  }
  print(direction_append)

  # Save estimates and CI of threshold-response function
  
  ## non-monotone
  if (!config$threshold_monotone_only) {
    esttmle <- esttmle_full[[1]]
    
    save("esttmle", file = here::here('output', TRIAL, COR,
                                      paste0("tmleThresh_",   marker, direction_append,".RData")
    ))
    write.csv(esttmle, file = here::here(
      "output", TRIAL, COR,
      paste0("tmleThresh_",  marker,direction_append, ".csv")
    ), row.names = F)
  }
  
  esttmle <- esttmle_full[[2]]
  save("esttmle", file = here::here(
    "output", TRIAL, COR,
    paste0("tmleThresh_monotone_", marker,direction_append, ".RData")
  ))
  write.csv(esttmle, file = here::here(
    "output", TRIAL, COR,
    paste0("tmleThresh_monotone_", marker,direction_append, ".csv")
  ), row.names = F)
  
  return(esttmle)
}

# print("WHATTTTTTTt")
#  print(markers)
# for (marker in markers) {
#   print(marker)
#     v <- run_threshold_analysis(marker, "below")
# }


for (marker in markers) {
  myprint(marker)
  run_threshold_analysis(marker, "above")
}


print("run time: "%.%format(Sys.time()-begin, digits=2))
