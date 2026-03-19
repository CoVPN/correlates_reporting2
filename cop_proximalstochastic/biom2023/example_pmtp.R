rm(list = ls())
### Load functions  
source("proxi_mtp_functions.R")

### Load data
sim_trial_data = read.csv("sim_trial_data.csv")

### Define policy functions

max(sim_trial_data$X)

policy_q1 <- function(x) {
    (x+0.4)*(x+0.4<=3.5-1)+(0.4*3.5+1*x)/(0.4+1)*(x+0.4>3.5-1)
  }
policy_q2 <- function(x) {
    (x+0.8)*(x+0.8<=3.5-1)+(0.8*3.5+1*x)/(0.8+1)*(x+0.8>3.5-1)
  }


### Recreate analysis showed in the supplement

### Estimation accounting for sampling weights

set.seed(123)
obj.pmtp = pmtp(data = sim_trial_data,
                trt = "X",
                outcome = "Y",
                covariates = c("L1", "L2", "L3"),
                nct = "Z",
                nco = "W",
                weights = "wt",
                policy = list(policy_q1, policy_q2),
                control.folds = TRUE)

summary(obj.pmtp)

### Estimation ignoring sampling weights

set.seed(123)
obj.pmtp.no.wt = pmtp(data = sim_trial_data,
                      trt = "X",
                      outcome = "Y",
                      covariates = c("L1", "L2", "L3"),
                      nct = "Z",
                      nco = "W",
                      weights = NULL,
                      policy = list(policy_q1, policy_q2),
                      control.folds = TRUE)

summary(obj.pmtp.no.wt)

