# to run this file:
# ~/correlates_reporting2/cor_coxph# Rscript code/explore_age.R

time.start = Sys.time()

renv::activate("..")
library(kyotil)
library(vaccine)

if (!dir.exists("output/age_risk/")) dir.create("output/age_risk/")


################################################################################
# AZ 


# stage 1
if (!file.exists("output/age_risk/azd1222.stage1.Rdata")) {
  print("AZ COVID")
  
  dat_azd1222_stage1=read.csv('/trials/covpn/p3002/analysis/correlates/Part_A_Blinded_Phase_Data/adata/azd1222_data_processed_with_riskscore.csv')
  
  
  for (trt in 0:1) {
    fits=list()
    
    dat=subset(dat_azd1222_stage1, ph1.D57==1 & Trt==trt & Bserostatus==0 & Country==2)
    dat$wt=1
    dat$Trt=1 # pretend to be vacc so that the code will work
    
    dat <- load_data(time="EventTimePrimaryD57", event="EventIndPrimaryD57", vacc="Trt",
                     marker="Age", covariates=c("Sex"),
                     weights="wt", ph2="wt", data=dat)
    
    tfinal.tpeak=92
    
    fits$ests_cox <- est_ce(dat=dat, type="Cox", t_0=tfinal.tpeak)
    
    fits$ests_cox_spline <- est_ce(dat=dat, type="Cox", t_0=tfinal.tpeak, params_cox = params_ce_cox(spline_df=3))
    
    # # non-monotone fit
    # fits$ests_np_nonmonotone <- est_ce(dat=dat, type="NP", t_0=tfinal.tpeak, params_np=params_ce_np(mono_cis=F))
    # monotone fit
    fits$ests_np <- est_ce(dat=dat, type="NP", t_0=tfinal.tpeak)
    
    if (trt==0) azd1222.stage1.plac = fits    
    if (trt==1) azd1222.stage1.vacc = fits    
  }

  save(azd1222.stage1.plac, azd1222.stage1.vacc, file="output/age_risk/azd1222.stage1.Rdata")
}



################################################################################
# NVX 

# prevent19 stage 1 
if (!file.exists("output/age_risk/prevent19.stage1.Rdata")) {
  print("NVX prevent19 stage 1")
  
  dat_prevent19_stage1=read.csv('/trials/covpn/p3004/analysis/correlates/Part_A_Blinded_Phase_Data/adata/prevent19_data_processed_20240325.csv')
  
  for (trt in 0:1) {
    fits=list()
  
    dat=subset(dat_prevent19_stage1, ph1.D35==1 & Trt==trt & Bserostatus==0 & Country==0)
    dat$wt=1
    dat$Trt=1 # pretend to be vacc so that the code will work
    
    dat <- load_data(time="EventTimePrimaryD35", event="EventIndPrimaryD35", vacc="Trt",
                     marker="Age", covariates=c("Sex"),
                     weights="wt", ph2="wt", data=dat)
    
    tfinal.tpeak=59
    
    fits$ests_cox <- est_ce(dat=dat, type="Cox", t_0=tfinal.tpeak)
    
    fits$ests_cox_spline <- est_ce(dat=dat, type="Cox", t_0=tfinal.tpeak, params_cox = params_ce_cox(spline_df=3))
    
    # monotone fit
    fits$ests_np <- est_ce(dat=dat, type="NP", t_0=tfinal.tpeak)
    
    if (trt==0) prevent19.stage1.plac = fits    
    if (trt==1) prevent19.stage1.vacc = fits    
  }
  
  save(prevent19.stage1.plac, prevent19.stage1.vacc, file="output/age_risk/prevent19.stage1.Rdata")
}



# prevent19 stage 2
if (!file.exists("output/age_risk/prevent19.stage2.Rdata")) {
  print("NVX prevent19 stage 2")
  
  dat_prevent19_stage2=read.csv('/trials/covpn/p3004/analysis/correlates/stage2/adata/prevent19_stage2_data_processed_20240401.csv')
  
  for (i in 1:2) { # 1: delta, 2: severe
    fits=list()
    
    dat=subset(dat_prevent19_stage2, ph1.D35_108==1 & Trt==1 & Bserostatus==0 & Country==0)
    dat$wt=1
    dat$Trt=1 # pretend to be vacc so that the code will work
    
    tfinal.tpeak=305
    
    if (i==1) {
      dat <- load_data(time="COVIDTimeD35to21Dec10", event="KnownOrImputedDeltaCOVIDIndD35_108to21Dec10", vacc="Trt",
                       marker="Age", covariates=c("Sex"),
                       weights="wt", ph2="wt", data=dat)
      
      fits$ests_cox <- est_ce(dat=dat, type="Cox", t_0=tfinal.tpeak)
      
      fits$ests_cox_spline <- est_ce(dat=dat, type="Cox", t_0=tfinal.tpeak, params_cox = params_ce_cox(spline_df=3))

      # monotone fit
      fits$ests_np <- est_ce(dat=dat, type="NP", t_0=tfinal.tpeak)
      
      prevent19.stage2.delta = fits    
      
    } else if (i==2) {
      dat <- load_data(time="SevereCOVIDTimeD35to21Dec10", event="SevereCOVIDIndD35_108to21Dec10", vacc="Trt",
                       marker="Age", covariates=c("Sex"),
                       weights="wt", ph2="wt", data=dat)
      
      fits$ests_cox <- est_ce(dat=dat, type="Cox", t_0=tfinal.tpeak)
      
      fits$ests_cox_spline <- est_ce(dat=dat, type="Cox", t_0=tfinal.tpeak, params_cox = params_ce_cox(spline_df=3))
      
      # monotone fit
      fits$ests_np <- est_ce(dat=dat, type="NP", t_0=tfinal.tpeak)
      
      prevent19.stage2.severe = fits
    }
    
  }

  save(prevent19.stage2.delta, prevent19.stage2.severe, file="output/age_risk/prevent19.stage2.Rdata")
}



# UK302 stage 1 
if (!file.exists("output/age_risk/uk302.Rdata")) {
  print("NVX UK302")
  
  dat_uk302=read.csv('/trials/covpn/p3004/analysis/correlates/UK302/adata/nvx_uk302_data_processed_20240320.csv')
  
  for (trt in 0:1) {
    fits=list()
    
    dat=subset(dat_uk302, ph1.D35==1 & Trt==trt & Bserostatus==0)
    dat$wt=1
    dat$Trt=1 # pretend to be vacc so that the code will work
    
    dat <- load_data(time="EventTimePrimaryD35", event="EventIndPrimaryD35", vacc="Trt",
                     marker="Age", covariates=c("Sex"),
                     weights="wt", ph2="wt", data=dat)
    
    tfinal.tpeak=134
    
    fits$ests_cox <- est_ce(dat=dat, type="Cox", t_0=tfinal.tpeak)
    
    fits$ests_cox_spline <- est_ce(dat=dat, type="Cox", t_0=tfinal.tpeak, params_cox = params_ce_cox(spline_df=3))
    
    # monotone fit
    fits$ests_np <- est_ce(dat=dat, type="NP", t_0=tfinal.tpeak)
    
    if (trt==0) uk302.plac = fits    
    if (trt==1) uk302.vacc = fits    
  }
  
  save(uk302.plac, uk302.vacc, file="output/age_risk/uk302.Rdata")
}




print("cor_coxph run time: " %.% format(Sys.time() - time.start, digits = 1))
