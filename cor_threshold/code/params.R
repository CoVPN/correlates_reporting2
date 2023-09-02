renv::activate(project = here::here(".."))
source(here::here("..", "_common.R"))



# Reference time to perform analysis. Y = 1(T <= tf) where T is event time of Covid.
# tf should be large enough that most events are observed but small enough so that not many people are right censored. For the practice dataset, tf = 170 works.
# Right-censoring is taken into account for  this analysis.
covariate_adjusted <-
  T #### Estimate threshold-response function with covariate adjustment
fast_analysis <-
  F ### Perform a fast analysis using glmnet at cost of accuracy
super_fast_analysis <- T
# hack
threshold_grid_size <- 2 ### Number of thresholds to estimate (equally spaced in quantiles). Should be 15 at least for the plots of the threshold-response and its inverse to be representative of the true functions.
plotting_assay_label_generator <- function(marker, above = T) {
  if (above) {
    add <- " (>=s)"
  } else {
    add <- " (<=s)"
  }
  day <- ""
  
  time <- paste0(DayPrefix, tpeak)
  assay <- marker_to_assay[[marker]]
  
  labx <- labels.axis[time, assay]
  labx <- paste0(labx, add)
  
  return(labx)
  
}

plotting_assay_title_generator <- function(marker) {
  day <- ""
  time <- paste0(DayPrefix, tpeak)
  assay <- marker_to_assay[[marker]]
  title <- labels.title[time, assay]
  
  return(title)
  
}

if (TRIAL=='moderna_boost') {
  assays=c("bindSpike_BA.1", "pseudoneutid50_BA.1", "bindSpike", "pseudoneutid50")
} else {
  assays <- config$assays
}


time <- paste0(DayPrefix, tpeak)
key <- COR
markers <- paste0(DayPrefix, tpeak, assays)

markers <- intersect(markers, colnames(dat.mock))
print(markers)
marker_to_assay <- sapply(markers, function(v) {
  unname(gsub(paste0(DayPrefix, tpeak),  "", v))
})




# max_t <- max(dat.mock[dat.mock$EventIndPrimary==1 & dat.mock$Trt == 1 & dat.mock$ph2 == 1, "EventTimePrimary" ])
max_t = tfinal.tpeak
# Covariates to adjust for. SHOULD BE AT LEAST TWO VARIABLES OR GLMNET WILL ERROR

data_name_updated <- sub(".csv", "_with_riskscore.csv", data_name)

covariates = strsplit(sub("~", "", config$covariates_riskscore), "\\+")[[1]][-1]
covariates = kyotil::trim(covariates)

if ("risk_score" %in% covariates) {
  append_data <- "_with_riskscore"
} else {
  append_data <- ""
}


####################
#### Internal variables
###################
# Threshold grid is generated equally spaced (on the quantile level)
# between the threshold at lower_quantile and threshold at upper_quantile
# Best to have these be away from 0 and 1.
lower_quantile <- 0
upper_quantile <- 1
risks_to_estimate_thresh_of_protection <- NULL ### Risk values at which to estimate threshold of protection...
###                Leaving at NULL (default) is recommended to ensure risks fall in range of estimate. Example: c(0, 0.0005, 0.001, 0.002, 0.003)
threshold_grid_size <- max(threshold_grid_size, 15)
