# this file should not source _common.R. ~ Youyi

library(here)
library(stringr)

labels.title2 <- apply(labels.title, c(1, 2), function(st) {
  str_replace(st, ":", "\n")
})

if (study_name!="IARCHPV") {
    trt.labels <- c("Placebo", "Vaccine")
    bstatus.labels <- c("Baseline Neg", "Baseline Pos")
    bstatus.labels.2 <- c("BaselineNeg", "BaselinePos")
} else {
    trt.labels <- c("Single-dose","Two-dose","Two doses default","Three-dose")
    bstatus.labels <- "" # no baseline serostatus for the IARCHPV study
    bstatus.labels.2 <- ""
    }


