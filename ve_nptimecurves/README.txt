This README is for implementing VE analysis for the CoVPN VAT00008 efficacy trial of bivalent vaccine. 

########################### 
# Software and package installation: 

## Install R 4.4.3, Rstudio 2024.12.1.563. Both softwares are free and available for most operating systems. 


## Install packages in Rstudio: open Rstudio and execute the following lines of R code:

install.packages("devtools")
library(devtools)
install_version("tidyverse", version = "2.0.0")
install_version("tidycmprsk, version = "1.1.0")

########################### 
# Step 00: Place the dataset file "vat08_combined_data_processed_20250321.csv" in the folder `data`

###########################

# Step 01: Generate cumulative incidence plot. 

## Source R file cumuincplot.R in R or Rstudio



###########################

# Step 02: Generate plot for instantaneous hazard-based VE. 

## Source R file instantaneousVE.R in R or Rstudio


###########################

# Step 03: Generate plot for cumulative incidence-based VE. 

## Source R file cumulativeVE.R in R or Rstudio

