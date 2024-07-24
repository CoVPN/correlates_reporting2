#Sys.setenv(TRIAL = "hvtn705second"); COR="D210"; Sys.setenv(VERBOSE = 1) 
#Sys.setenv(TRIAL = "moderna_boost"); COR="BD29naive"; Sys.setenv(VERBOSE = 1) 

renv::activate(project = here::here(".."))
source(here::here("..", "_common.R"))
source(here::here("code", "params.R"))

library(npthreshold)
# library(uuid)
# library(doMC)
# library(earth)

begin=Sys.time()


################################################################################
# data preparation

{
  
if(fast_analysis) {
  threshold_grid_size <- 15
} else if(super_fast_analysis) {
  threshold_grid_size <- 15
} # default is 30 in the config

  
# Generate the outcome and censoring indicator variables
short_key <- COR

data <- dat_proc
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

}

################################################################################
# Runs threshold analysis and saves raw results.


for (marker in markers) {
# marker=markers[1]
  
  form <- paste0("Y~h(.) + h(t,.)  + h(A,.)")
  
  thresholds <- read.csv(here::here('output', TRIAL, COR, "data_clean", "Thresholds_by_marker", paste0("thresholds_", marker, ".csv")))
  thresholds <- as.vector(unlist(thresholds[, 1])) 
  print(thresholds)
  #time <- marker_to_time[[marker]]
  
  data_full <- read.csv(here::here('output', TRIAL, COR, "data_clean", paste0("data_firststage", ".csv")))
  
  print("Run thresholdTMLE")
  
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

  # Save estimates and CI of threshold-response function
  
  direction_append <- ""
  
  # non-monotone
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
  
  # monotone
  esttmle <- esttmle_full[[2]]
  save("esttmle", file = here::here(
    "output", TRIAL, COR,
    paste0("tmleThresh_monotone_", marker,direction_append, ".RData")
  ))
  write.csv(esttmle, file = here::here(
    "output", TRIAL, COR,
    paste0("tmleThresh_monotone_", marker,direction_append, ".csv")
  ), row.names = F)
  
}



print("run time: "%.%format(Sys.time()-begin, digits=2))
