#Sys.setenv(TRIAL = "moderna_boost"); COR="BD29"; Sys.setenv(VERBOSE = 1) 

renv::activate(project = here::here(".."))

library(cowplot)
library(scales)
library(knitr)
library(dplyr)
library(magrittr)
library(ggplot2)

source(here::here("code", "params.R"))

source(here::here("code", "learners.R"))
source(here::here("code", "tmleThresh.R"))
source(here::here("code", "plotting_helpers.R"))
ident <- function(x) x

# Plotting arguments


# for (key in keys) {
#     above <- F
#   get_plot(key, simultaneous_CI = F, monotone = F, above)
#   get_plot(key, simultaneous_CI = T, monotone = F, above)
#   get_plot(key, simultaneous_CI = F, monotone = T,above)
#   get_plot(key, simultaneous_CI = T, monotone = T,above)
#   generate_tables(key, num_show = 10, monotone = F,above)
#   generate_tables(key, num_show = 10, monotone = T,above)
#   #get_inverse_plot(marker, F)
#   #get_inverse_plot(marker, T)
# }


for (marker in markers) {
  get_plot(marker, simultaneous_CI = F, monotone = T, above=T)
  get_plot(marker, simultaneous_CI = T, monotone = T, above=T)
  generate_tables(marker, num_show = 10, monotone = T,above=T)
  
  # get_plot(marker, simultaneous_CI = F, monotone = F, above=T)
  # get_plot(marker, simultaneous_CI = T, monotone = F, above=T)
  # generate_tables(marker, num_show = 10, monotone = F,above=T)
  
  #get_inverse_plot(marker, F)
  #get_inverse_plot(marker, T)
}
