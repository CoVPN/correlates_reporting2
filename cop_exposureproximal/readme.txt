Exposure_proximal

Evaluating biomarkers as exposure-proximal correlates of disease risk and vaccine efficacy in baseline-naive population

This folder contains the following files:

1. Example.R  --- R code to read in the simulated data and conduct exposure-proximal correlates analysis,
it also includes the output of the R code implementation 

2. FunVacc.R  --- contains major functions for fitting Cox models with various components out of (x(t),x(0),t),
based on data from vaccine recipients only

3. FunVacPlac.R --- contains major functions for fitting Cox models with various components out of (x(t),x(0),t),
based on data from vaccine and placebo recipients together

4. example.Rdata --- simulated R data of an immune correlate study nested within a vaccine trial, it consists of the following two 
data object

#####
(1) data: the vaccine trial data, it includes the following variables:

### ID: participant id
### ZE: treatment indicator: 0, 1 for placebo and vaccine
### WE: covariate to adjust for in the linear mixed effects model (LME) of immune response trajectory
### WEc: confounder to adjust for in the Cox model
### RE: peak time point after vaccination, measured on calendar scale
### SE: minimum of infection and censoring time, measured on calendar scale
### TE: SE-RE, time post peak to SE
### XE0: peak immune response
### delta: case/control indicator, 1, 0 for case and control respectively
### incoh: binary indicator that a vaccine recipient is included in the immune correlates study
### weight: inverse of probability of sampling into the immune correlates study for vaccine recipients, 1 for placeob 

######

(2) datlong: the immune correlates study data among vaccine recipients in long format, it includes the following variables

### id: participant id, should be a subset of ID in "data"
### visit: visit number for measuring immune responses
### R: peak time point after vaccination, measured on calendar scale
### time: time post peak to immune response measuremnent
### X: immune response measure
### W: covariate to adjust for in the linear mixed effect model
### S: minimum of infection and censoring time, measured on calendar scale
### Z: treatment indicator: 0, 1 for placebo and vaccine
### weights: inverse of probability of sampling into the immune correlates study for vaccine recipients, 0 for placebo 
(i.e. LME is only fitted using vaccinee data)
