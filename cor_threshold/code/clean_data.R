#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))

#-----------------------------------------------
# load parameters

source(here::here("code", "params.R"))



# load data
#dat.mock <- read.csv(here::here("..", "data_clean", paste0(stringr::str_match(data_name,"(.+).csv")[,2],append_data,".csv")))

 
print("HERE")
for (a in assays) {
  for (t in "Day"%.%tpeak ) {
    
    print(t %.% a)
    print( log10(uloqs[a]))
    dat.mock[[t %.% a]] <- ifelse(dat.mock[[t %.% a]] > log10(uloqs[a]), log10(uloqs[a]), dat.mock[[t %.% a]])
  }
}    


 

 

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
  variables_to_keep <-
    c(
      covariates,
      markers,
      "TwophasesampInd",
      #"grp",
      "wt",
      "Ttilde",
      "Delta",
      "target_time",
      "J"
    )
  # TEMPORARY
  variables_to_keep <- intersect(variables_to_keep, colnames(data))
  
  
  
  #keep <- data[[Earlyendpoint]] ==0 & data$Trt == 1 & data$Bserostatus == 0 & data$Perprotocol==1 & !is.na(data$wt) & data[[Event_Time_variable[[time]]]] >=7 & !is.na(data$Wstratum)
  keep <- data[["ph1"]]==1 & data$Trt == 1  
  
  
  data_firststage <- data[keep, variables_to_keep]
  
  #data_firststage <- na.omit(data_firststage)
  data_secondstage <- data_firststage[data_firststage$TwophasesampInd == 1, ]
  
  write.csv(data_firststage,
            here::here("data_clean", paste0("data_firststage_", short_key, ".csv")),
            row.names = F
  )
  write.csv(data_secondstage,
            here::here("data_clean", paste0("data_secondstage_", short_key, ".csv")),
            row.names = F
  )
  
 
  thresholds_list <- list()
  for (marker in markers) {
    
    if (length(unique(data_secondstage[[marker]])) > threshold_grid_size) {
      # Choose grid of 15 thresholds
      # Make sure the thresholds are supported by data (so that A >=v has enough people)
      #  seq(0, 1, length.out = threshold_grid_size) will error code
      # Upper threshold must be less than the 0.99 quantile
      thresh_grid <-
        unique(quantile(
          data_secondstage[[marker]][data_secondstage[["Delta"]]==1],
          seq(0,1,length.out=30),
          na.rm = T
        ))
      
       
      thresh_mand <- report.assay.values(data_secondstage[[marker]][data_secondstage[["Delta"]]==1], marker)
      thresh_grid <- sort(union(thresh_grid, thresh_mand))
       
      #thresh_grid <- sort(union(thresh_mand, thresh_grid))
    } else {
      
      thresh_grid <- sort(unique(data_secondstage[[marker]]))
    }
    
    
    write.csv(data.frame(thresh = thresh_grid), here::here("data_clean", "Thresholds_by_marker", paste0("thresholds_", marker, ".csv")), row.names = F)
    
    #thresholds_list[[marker]] <- thresh_grid
  }
  






