#Sys.setenv(TRIAL = "moderna_boost"); COR="BD29naive"; Sys.setenv(VERBOSE = 1) 

renv::activate(project = here::here(".."))

source(here::here("..", "_common.R"))

source(here::here("code", "params.R"))

library(cowplot)
library(scales)
library(knitr)
library(dplyr)
library(magrittr)
library(ggplot2)

source(here::here("code", "learners.R"))
source(here::here("code", "tmleThresh.R"))
source(here::here("code", "plotting_helpers.R"))
ident <- function(x) x


# need to redefine markers since not all markers have thresholds
files = dir(path=here::here('output', TRIAL, COR, "data_clean", "Thresholds_by_marker"), pattern="thresholds_.*\\.csv")
markers = sub(".csv", "", sub("thresholds_","",files))


for (marker in markers) {
  message("monotone")
  get_plot(marker, simultaneous_CI = F, monotone = T, above=T)
  get_plot(marker, simultaneous_CI = T, monotone = T, above=T)
  generate_tables(marker, num_show = 10, monotone = T,above=T)
  
  message("non-monotone")
  get_plot(marker, simultaneous_CI = F, monotone = F, above=T)
  get_plot(marker, simultaneous_CI = T, monotone = F, above=T)
  generate_tables(marker, num_show = 10, monotone = F,above=T)
  
  #get_inverse_plot(marker, F)
  #get_inverse_plot(marker, T)
}
