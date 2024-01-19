# this file should not source _common.R. ~ Youyi

library(here)
library(stringr)

#labels.title2 <- apply(labels.title, c(1, 2), function(st) {
#  str_replace(st, ":", "\n")
#})

if (study_name=="IARCHPV") {
  trt.labels <- c("Single-dose","Two-dose","Two doses default","Three-dose")
  trt.labels2 <- c("Single-dose","Two-dose (Days 1 and =180)","Two doses default (Days 1 and 60)","Three-dose Three_dose (Days 1, 60 and =180)")
  bstatus.labels <- "" # no baseline serostatus for the IARCHPV study
  bstatus.labels.2 <- ""
  
  } else if (study_name=="VAT08"){
    trt.labels <- c("Placebo", "Vaccine")
    bstatus.labels <-  c("Naive", "Non-naive")
    bstatus.labels.2 <- c("naive", "non-naive")
    
  } else {
    trt.labels <- c("Placebo", "Vaccine")
    bstatus.labels <- c("Baseline Neg", "Baseline Pos")
    bstatus.labels.2 <- c("BaselineNeg", "BaselinePos")
    
  }


