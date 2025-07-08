# this file should not source _common.R. ~ Youyi

library(here)
library(stringr)

labels.title2 <- apply(labels.title, c(1, 2), function(st) {
  str_replace(st, ":", "\n")
})

if (study_name=="IARCHPV") {
  trt.labels <- c("Single-dose","Two-dose","Two doses default","Three-dose")
  trt.labels2 <- c("Single-dose","Two-dose (Days 1 and =180)","Two doses default (Days 1 and 60)","Three-dose Three_dose (Days 1, 60 and =180)")
  bstatus.labels <- "" # no baseline serostatus for the IARCHPV study
  bstatus.labels.2 <- ""
  
  } else if (study_name=="VAT08" | grepl("stage2", COR)){
    trt.labels <- c("Placebo", "Vaccine")
    bstatus.labels <-  c("Naive", "Non-naive")
    bstatus.labels.2 <- c("naive", "non-naive")
    
  } else if (study_name == "NextGen_Mock") {
    trt.labels <- c("Comparator Vaccine", "Investigational Vaccine")
    bstatus.labels <-  c("Naive", "Non-naive")
    bstatus.labels.2 <- c("", "")
    
  } else {
    trt.labels <- c("Placebo", "Vaccine")
    bstatus.labels <- c("Baseline Neg", "Baseline Pos")
    bstatus.labels.2 <- c("BaselineNeg", "BaselinePos")
    
  }

# create a new variables times_ and set it to times unless the study needs plots for timepoint more than what's included in the times
if(attr(config,"config")=="janssen_pooled_partA") {
  times_ = c("B","Day29","Day71","Mon6")
  labels.time = c("Day 1","Day 29", "Day 71", "Month 6"); names(labels.time) = times_
} else if (attr(config,"config")=="prevent19_stage2") {
  times_ = c("Day35","C1","BD1","DD1")
  labels.time = c("Day 35", "Crossover Day 1", "Booster Day 1", "Disease Day 1"); names(labels.time) = times_
} else if (attr(config,"config")=="azd1222_stage2") {
  times_ = c("Day57","Day90","Day180","Day360")
  labels.time = c("Day 57","Day 90", "Day 180", "Day 360"); names(labels.time) = times_
} else if (attr(config,"config")=="nextgen_mock") {
  times_ = c("B", "Day31", "Delta31overB", "Day91", "Day181", "Day366")
  labels.time = c("Day 1","Day 31", "D31 fold-rise over D1", "Day 91", "Day 181", "Day 366"); names(labels.time) = times_
} else {times_ = times}




