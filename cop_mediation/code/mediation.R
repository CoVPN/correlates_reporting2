#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
source(here::here("..", "_common.R"))
#-----------------------------------------------

library(SuperLearner)
source(here::here("code", "sl_screen_fn.R"))
source(here::here("code", "format_utils.R"))
data <- dat.mock

tf_Day <- max(data[data$EventIndPrimary==1 & data$Trt == 1 & data$ph2, "EventTimePrimary" ])
print(
  paste0("The follow-up day used to define primary binary endpoint is: ", tf_Day)
)
if(config$study_name_code == "ENSEMBLE"){
  times <- "Day29"
  covariates <- c("risk_score", "Region1", "Region2")
  data$Region1 <- as.numeric(data$Region == 1)
  data$Region2 <- as.numeric(data$Region == 2)
}else if(config$study_name_code == "COVE"){
  times <- c("Day29") # Day57 all have positivity issues
  covariates <- c("MinorityInd", "HighRiskInd", "risk_score")
}

include_assays <- NULL
for(a in assays){
	assay_name <- paste0(times, a)
	placebo_assay_values <- data[[assay_name]][data$Trt == 0]
	vaccine_assay_values <- data[[assay_name]][data$Trt == 1]
	#max_placebo_assay_value <- max(placebo_assay_values, na.rm = TRUE)
  # prop_vaccine_less_than_max_placebo <- mean(vaccine_assay_values <= max_placebo_assay_value, na.rm = TRUE)
  placebo_pct90_value <- quantile(placebo_assay_values, p = 0.9, na.rm = TRUE)
  prop_vaccine_less_than_max_placebo <- mean(vaccine_assay_values <= placebo_pct90_value, na.rm = TRUE)
	print(paste0(assay_name, " = ", round(prop_vaccine_less_than_max_placebo, 3)))
	if(prop_vaccine_less_than_max_placebo > 0.05){
	  include_assays <- c(include_assays, assay_name)
	}
}

data$outcome <- as.numeric(
  data$EventIndPrimary == 1 & data$EventTimePrimary <= tf_Day                      
)

variables_to_keep <- c(
  covariates,
  include_assays,
  config.cor$ph2,
  config.cor$wt,
  "Trt",
  "outcome"
)

data_keep <- data[!is.na(data[[config.cor$wt]]), variables_to_keep]
W <- data_keep[, covariates, drop = FALSE]
A <- data_keep$Trt
R <- as.numeric(data_keep[[config.cor$ph2]])
Y <- data_keep$outcome

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
  c("SL.glmnet", "screen_all"), # no need for screens
  c("SL.xgboost", "screen_all"), # no need for screens
  c("SL.ranger", "screen_all"),   # faster than cforest?
  c("SL.glm", "screen_glmnet"),
  c("SL.glm", "screen_univariate_logistic_pval"),
  c("SL.glm", "screen_highcor_random"),
  c("SL.glm.interaction", "screen_glmnet"),
  c("SL.glm.interaction", "screen_univariate_logistic_pval"),
  c("SL.glm.interaction", "screen_highcor_random"),
  c("SL.gam", "screen_glmnet"),
  c("SL.gam", "screen_univariate_logistic_pval"),
  c("SL.gam", "screen_highcor_random")
)

quant_result <- day_col <- assay_col <- NULL
for (marker in include_assays) {
  print(marker)

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
	this_row <- format_row(fit)
	quant_result <- rbind(quant_result, this_row)
  which_assay <- substr(marker, 6, nchar(marker)) ## this could be written better
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
           COR, ".rds")
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