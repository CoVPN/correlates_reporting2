#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
source(here::here("..", "_common.R"))
#-----------------------------------------------

library(SuperLearner)
source(here::here("code", "sl_screen_fn.R"))
source(here::here("code", "format_utils.R"))
data <- dat.mock

tf_Day <- tfinal.tpeak

print(
  paste0("The follow-up day used to define primary binary endpoint is: ", tf_Day)
)

impute_placebo_to_lod_over_2 <- FALSE
if(study_name == "ENSEMBLE"){
  if(length(unique(data$Region)) == 1){
    covariates <- "risk_score"
  }else{
    covariates <- c("risk_score", "Region1", "Region2")
    data$Region1 <- as.numeric(data$Region == 1)
    data$Region2 <- as.numeric(data$Region == 2)
  }
  run_survtmle <- TRUE
  cut_size <- 6
  impute_placebo_to_lod_over_2 <- TRUE
}else if(study_name == "COVE"){
  covariates <- c("MinorityInd", "HighRiskInd", "risk_score")
  run_survtmle <- FALSE
  cut_size <- 6
}else if(study_name == "AZD1222"){
  covariates <- c("Age", "risk_score")
  run_survtmle <- TRUE
  cut_size <- 6
}else if(study_name == "HVTN705"){
  covariates <- c("RSA", "Age", "BMI", "Riskscore")
  run_survtmle <- FALSE
  cut_size <- 7
  impute_placebo_to_lod_over_2 <- TRUE
}


time_long <- paste0("Day", config.cor$tpeak)
time_short <- paste0("D", config.cor$tpeak)

include_assays <- NULL
threshold <- 0.05 # must have this % vax with no response

assay_names <- paste0(time_long, assays)
for(a in assay_names){
	placebo_assay_values <- data[[a]][data$Trt == 0]
	vaccine_assay_values <- data[[a]][data$Trt == 1]
	#max_placebo_assay_value <- max(placebo_assay_values, na.rm = TRUE)
  # prop_vaccine_less_than_max_placebo <- mean(vaccine_assay_values <= max_placebo_assay_value, na.rm = TRUE)
  placebo_pct90_value <- quantile(placebo_assay_values, p = 0.9, na.rm = TRUE)
  prop_vaccine_less_than_max_placebo <- mean(vaccine_assay_values <= placebo_pct90_value, na.rm = TRUE)
	print(paste0(a, " = ", round(prop_vaccine_less_than_max_placebo, 3)))
	if(prop_vaccine_less_than_max_placebo > threshold){
	  include_assays <- c(include_assays, a)
	}
}

# for binary endpoints
data$outcome <- as.numeric(
  data$EventIndPrimary == 1 & data$EventTimePrimary <= tf_Day                      
)

variables_to_keep <- c(
  covariates,
  include_assays,
  config.cor$ph2,
  config.cor$wt,
  "Trt",
  paste0("EventTimePrimary"),
  paste0("EventIndPrimary"),
  "outcome"
)

data_keep <- data[!is.na(data[[config.cor$wt]]), variables_to_keep]
row.names(data_keep) <- seq_len(dim(data_keep)[1])

# remove events after tf_Day from estimation
data_keep$EventTimePrimary[data_keep$EventTimePrimary >= tf_Day] <- tf_Day
data_keep$EventIndPrimary[data_keep$EventTimePrimary >= tf_Day] <- 0

W <- data_keep[, covariates, drop = FALSE]
A <- data_keep$Trt
R <- as.numeric(data_keep[[config.cor$ph2]])
Y <- data_keep$outcome

if(impute_placebo_to_lod_over_2){
  R[A == 0] <- 1
  for(a in include_assays){
    impute_idx <- which(
      is.na(data_keep[,a]) & (A == 0)
    )
    data_keep[impute_idx, a] <- min(data_keep[A == 0, a], na.rm = TRUE)
  }
  data_keep$wt[A == 0] <- 1
}

# summaries for report
tab1 <- data.frame(table(Y[R == 1], A[R==1]))
colnames(tab1) <- c("Endpoint", "Vaccine", "n")

tab2 <- data.frame(table(Y, A))
colnames(tab2) <- c("Endpoint", "Vaccine", "n")

summary_stats <- list(
  tab1 = tab1,
  tab2 = tab2,
  day_of_analysis = tf_Day
)
saveRDS(
  summary_stats, 
  file = here::here("output", 
    paste0("summary_stats_", 
           attr(config,"config"), "_",
           COR, ".rds")
  )
)

sl_library <- list(
  c("SL.mean", "screen_all"),
  c("SL.glm", "screen_all"),
  #c("SL.glmnet", "screen_all"), # no need for screens
  #c("SL.xgboost", "screen_all"), # no need for screens
  c("SL.ranger", "screen_all"),   # faster than cforest?
  c("SL.gam", "screen_all"),
  c("SL.earth", "screen_all")
)

sl_library <- "SL.mean"

# results without marker used to compute PM
# does not need to be updated with each marker
if(run_survtmle){
    fit1 <- survtmle::hazard_tmle(
      ftime = data_keep$EventTimePrimary,
      ftype = data_keep$EventIndPrimary,
      trt = data_keep$Trt,
      adjustVars = data_keep[ , covariates],
      t0 = tf_Day,
      SL.ctime = sl_library,
      SL.ftime = sl_library, 
      verbose = TRUE,
      returnModels = TRUE
    )
    print(fit1)
}

if(!is.null(Args[2])){
  assay_for_this_run <- Args[2]
}else{
  assay_for_this_run <- include_assays
}

quant_result <- day_col <- assay_col <- NULL
for (marker in include_assays[include_assays %in% assay_for_this_run]) {
  print(marker)

  if(!run_survtmle){
    fit <- natmed2::natmed2(
      W = W,
      A = A,
      R = R,
      S = data_keep[, marker, drop = TRUE],
      C = rep(1, nrow(data_keep)),
      Y = Y,
      gRn = 1 / data_keep$wt,
      glm_gA = ".",
      glm_gAS = NULL,
      SL_gAS = sl_library,
      glm_QY_WAS = NULL,
      SL_QY_WAS = sl_library,
      glm_QY_WACY = NULL,
      SL_QY_WACY = sl_library,
      glm_QD_WACY = NULL,
      SL_QD_WACY = sl_library,
      glm_QY_W = NULL,
      SL_QY_W = sl_library,
      glm_QY_WA = NULL,
      SL_QY_WA = sl_library
  	)
  }else{
    fit2 <- survtmle::hazard_tmle(
      ftime = data_keep$EventTimePrimary,
      ftype = data_keep$EventIndPrimary,
      trt = data_keep$Trt,
      adjustVars = data_keep[ , covariates, drop = FALSE],
      mediator = data_keep[ , marker, drop = FALSE],
      mediatorTrtVal = 0,
      trtOfInterest = 1,
      mediatorSampProb = 1 / data_keep$wt,
      mediatorInCensMod = FALSE,
      t0 = tf_Day,
      SL.ctime = fit1$ctimeMod,
      SL.ftime = sl_library,
      SL.mediator = sl_library,
      SL.trtMediator = sl_library,
      SL.eif = sl_library,
      verbose = TRUE
    )
    print(fit2)
    fit <- compute_mediation_params(fit1, fit2)
    print(fit)
  }
	this_row <- format_row(fit)
	quant_result <- rbind(quant_result, this_row)
  which_assay <- substr(marker, cut_size, nchar(marker)) ## this could be written better
	assay_col <- c(assay_col, gsub("\\%", "\\\\\\%", labels.assays[which_assay]))
}

full_result <- cbind(assay_col, quant_result)
row.names(full_result) <- NULL
colnames(full_result) <- c("Assay", "Direct VE", "Indirect VE", "Prop. mediated")

saveRDS(
  full_result, 
  file = here::here("output", 
    paste0("full_result_", 
           attr(config,"config"), "_",
           COR, "_", 
           paste0(include_assays[include_assays %in% assay_for_this_run], collapse = "_"),
           ".rds")
  )
)

# # re-labeling
# for(COR in c("D29IncludeNotMolecConfirmed", "D29IncludeNotMolecConfirmedstart1")){
#   bind_result <- readRDS(here::here("output",
#     paste0("full_result_janssen_pooled_real_", COR, ".rds")
#   ))
#   bind_result[,'Assay'] <- c("Binding Antibody to Spike", "Binding Antibody to RBD")
#   saveRDS(bind_result, file = here::here("output",
#     paste0("full_result_janssen_pooled_real_", COR, ".rds")
#   ))
#   adcp_result <- readRDS(here::here("output",
#     paste0("full_result_janssen_pooled_realADCP_", COR, ".rds")
#   ))
#   adcp_result[,'Assay'] <- "Phagocytic Score"
#   saveRDS(adcp_result, file = here::here("output",
#     paste0("full_result_janssen_pooled_realADCP_", COR, ".rds")
#   ))
# }