
#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
source(here::here("..", "_common.R"))
#-----------------------------------------------
print("ADE")
TRIAL <- Sys.getenv("TRIAL") 

### #####
#### NOTE Currently only supports Day29 and Day57 markers
#########
 

# Reference time to perform analysis. Y = 1(T <= tf) where T is event time of Covid.
# tf should be large enough that most events are observed but small enough so that not many people are right censored. For the practice dataset, tf = 170 works.
# Right-censoring is taken into account for  this analysis.
covariate_adjusted <- T #### Estimate threshold-response function with covariate adjustment
fast_analysis <- F ### Perform a fast analysis using glmnet at cost of accuracy
super_fast_analysis <- F
threshold_grid_size <- 30 ### Number of thresholds to estimate (equally spaced in quantiles). Should be 15 at least for the plots of the threshold-response and its inverse to be representative of the true functions.
plotting_assay_label_generator <- function(marker, above = T) {
  if(above) {
    add <- " (>=s)"
  } else {
    add <- " (<=s)"
  }
  day <- ""
  
  time <- paste0("Day", tpeak)
  assay <- marker_to_assay[[marker]]
 
  labx <- labels.axis[time, assay]
  labx <- paste0(labx, add)
  
  return(labx)
  
}

plotting_assay_title_generator <- function(marker) {
  day <- ""
  time <- paste0("Day", tpeak)
  assay <- marker_to_assay[[marker]]
  title <- labels.title[time, assay]
  
  return(title)
  
}

assays <- config$assays
time <- paste0("Day", tpeak)
key <- COR
markers <- paste0("Day", tpeak, assays)
 
 print("OK")
markers <- intersect(markers, colnames(dat.mock) ) 
print(markers)
marker_to_assay <- sapply(markers, function(v) {
 unname(gsub(paste0("Day", tpeak),  "", v))
})
 
 
 

 max_t <- max(dat.mock[dat.mock$EventIndPrimary==1 & dat.mock$Trt == 1 & dat.mock$ph2 == 1, "EventTimePrimary" ])
 
# Covariates to adjust for. SHOULD BE AT LEAST TWO VARIABLES OR GLMNET WILL ERROR
data_name_updated <- sub(".csv", "_with_riskscore.csv", data_name)
if (file.exists(here::here("..", "data_clean", data_name_updated))) {
  add_risk_score <- T
  covariates <- c("MinorityInd", "HighRiskInd", "risk_score") # , "BRiskScore") # Add "age"?
  
} else {
  add_risk_score <- F
  covariates <- c("MinorityInd", "HighRiskInd") # , "BRiskScore") # Add "age"?
  
}
if(Sys.getenv("TRIAL") =="hvtn705") {
  covariates <- c("Riskscore",  "BMI", "RSA")
}
 
if("risk_score" %in% covariates) {
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
