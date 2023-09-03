
# path for figures and tables etc
save.results.to = here::here("output", TRIAL, COR, 'figs', 'pointwise_CI'); if (!dir.exists(save.results.to))  dir.create(save.results.to, recursive = TRUE)
save.results.to = here::here("output", TRIAL, COR, 'figs', 'simultaneous_CI'); if (!dir.exists(save.results.to))  dir.create(save.results.to, recursive = TRUE)
save.results.to = here::here("output", TRIAL, COR, 'data_clean', 'Thresholds_by_marker'); if (!dir.exists(save.results.to))  dir.create(save.results.to, recursive = TRUE)


# Reference time to perform analysis. Y = 1(T <= tf) where T is event time of Covid.
# tf should be large enough that most events are observed but small enough so that not many people are right censored. For the practice dataset, tf = 170 works.
# Right-censoring is taken into account for  this analysis.
covariate_adjusted <- T #### Estimate threshold-response function with covariate adjustment

fast_analysis <- config$threshold_fast ### Perform a fast analysis using glmnet at cost of accuracy
super_fast_analysis <- config$threshold_superfast 

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
  assay <- marker_to_assay[[marker]]
  # time <- paste0(DayPrefix, tpeak)
  time = sub(assay,"",marker)
  title <- labels.title[time, assay]
  
  return(title)
  
}

if (TRIAL=='moderna_boost') {
  assays=c("bindSpike_BA.1")
  # assays=c("bindSpike_BA.1", "pseudoneutid50_BA.1", "bindSpike", "pseudoneutid50")
} else {
  assays <- config$assays
}


time <- paste0(DayPrefix, tpeak)
key <- COR
markers <- paste0(DayPrefix, tpeak, assays)
if (TRIAL=='moderna_boost') markers = c(markers, paste0("DeltaBD29overBD1", assays))

markers <- intersect(markers, colnames(dat.mock))
marker_to_assay <- sapply(markers, function(v) marker.name.to.assay(v))



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
