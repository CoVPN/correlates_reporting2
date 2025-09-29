

if(TRIAL == "covail_tcell"){
  
  non_naive = FALSE
  #non_naive = TRUE
  
  trt.arms.for.analysis <- c(1:2, 4:12)
  # trt.arms.for.analysis <- c(1:2, 4:15)
  
  # baseline risk factors
  bRiskFactors_includes_insert.stage.info = FALSE
  
  dat_proc_updated <- dat_proc %>% 
    mutate(Delta15overBpseudoneutid50_BA.1_2fold = ifelse(Day15pseudoneutid50_BA.1 > (Bpseudoneutid50_BA.1 + log10(2)), 1, 0),
           Delta15overBpseudoneutid50_BA.1_4fold = ifelse(Day15pseudoneutid50_BA.1 > (Bpseudoneutid50_BA.1 + log10(4)), 1, 0)) %>%
    mutate(Trt = 1,   # defined to sync with code below!
           stage = case_when(stage == 1 ~ "Moderna",
                             stage == 2 ~ "Pfizer.1",
                             stage == 3 ~ "Sanofi",
                             stage == 4 ~ "Pfizer.2"),
           Trtgrp = case_when(Trtgrp == "Beta with no Omicron" ~ "Beta.no.Omicron",
                              Trtgrp == "monovalent Prototype insert only" ~ "mono.Proto.insert.only",
                              Trtgrp == "Omicron" ~ "Omicron")) 
  
  # Identify the endpoint variable
  if(COR == "D15to91covail_tcell"){ 
    endpoint <- "COVIDIndD22toD91"
  } else if(COR == "D15to181covail_tcell"){
    endpoint <- "COVIDIndD22toD181"
  }
  
  
  if(non_naive == TRUE){
    briskfactors <- c("risk_score")
    briskfactors_correction <- "Y ~ x + X$risk_score"
  } else if(bRiskFactors_includes_insert.stage.info == FALSE){
    briskfactors <- c("risk_score", "FOIstandardized")
    briskfactors_correction <- "Y ~ x + X$risk_score + X$FOIstandardized"
  } else if (bRiskFactors_includes_insert.stage.info == TRUE){
    briskfactors <- c("risk_score", "FOIstandardized", "stage_Pfizer.1", "stage_Sanofi", "Trtgrp_mono.Proto.insert.only", "Trtgrp_Omicron")
    briskfactors_correction <- "Y ~ x + X$risk_score + X$FOIstandardized + X$stage_Pfizer.1 + X$stage_Sanofi +X$Trtgrp_mono.Proto.insert.only + X$Trtgrp_Omicron"
  } 
  
  wt <- "wt.D15.tcell"            # already set by COR
  ptidvar <- "Ptid"
  
  individualMarkers <- c("Bpseudoneutid50_BA.1",
                         "Day15pseudoneutid50_BA.1",
                         "Bfrnt50_BA.1",
                         "Day15frnt50_BA.1",
                         "Bfrnt80_BA.1",
                         "Day15frnt80_BA.1",
                         primary,
                         secondary, 
                         exploratory)
  
  # markers of interest
  markerVars <- c("Bpseudoneutid50_BA.1",               "Day15pseudoneutid50_BA.1",           "Delta15overBpseudoneutid50_BA.1_2fold",   "Delta15overBpseudoneutid50_BA.1_4fold", 
                  #"Bpseudoneutid50_BA.1cat",            "Day15pseudoneutid50_BA.1cat",        "Delta15overBpseudoneutid50_BA.1cat", 
                  "Bcd4_IFNg.IL2_Wuhan.N",              "Bcd4_FS_Wuhan.N",                   
                  "Bcd4_IFNg.IL2_BA.4.5.S",             "Day15cd4_IFNg.IL2_BA.4.5.S",         "Bcd4_FS_BA.4.5.S",                   "Day15cd4_FS_BA.4.5.S",              
                  "Bcd4_IFNg.IL2.154_Wuhan.N",          "Bcd8_IFNg.IL2_BA.4.5.S",             "Day15cd8_IFNg.IL2_BA.4.5.S",         "Bcd4_IFNg.IL2.154_BA.4.5.S",        
                  "Day15cd4_IFNg.IL2.154_BA.4.5.S",     "Day15cd4_IL21_BA.4.5.S",             "Bcd4_154_Wuhan.N",                   "Bcd4_IFNg_Wuhan.N",                 
                  "Bcd4_IL2_Wuhan.N",                   "Day15cd4_154_Wuhan.N",               "Day15cd4_IFNg_Wuhan.N",              "Day15cd4_IFNg.IL2_Wuhan.N",         
                  "Day15cd4_IFNg.IL2.154_Wuhan.N",      "Day15cd4_IL2_Wuhan.N",               "Bcd4_IFNg.IL2_COV2.CON.S",           "Day15cd4_IFNg.IL2_COV2.CON.S",      
                  "Bcd8_IFNg.IL2_COV2.CON.S",           "Day15cd8_IFNg.IL2_COV2.CON.S",       "Bcd4_IFNg.IL2.154_COV2.CON.S",       "Day15cd4_IFNg.IL2.154_COV2.CON.S",  
                  "Day15cd4_IL21_COV2.CON.S",           "Bcd4_IFNg_COV2.CON.S",               "Day15cd4_IFNg_COV2.CON.S",           "Bcd4_IL2_COV2.CON.S",               
                  "Day15cd4_IL2_COV2.CON.S",            "Bcd4_TNFa_COV2.CON.S",               "Day15cd4_TNFa_COV2.CON.S",           "Bcd8_IFNg_COV2.CON.S",              
                  "Day15cd8_IFNg_COV2.CON.S",           "Bcd8_IL2_COV2.CON.S",                "Day15cd8_IL2_COV2.CON.S",            "Bcd8_TNFa_COV2.CON.S",              
                  "Day15cd8_TNFa_COV2.CON.S",           "Bcd4_154_COV2.CON.S",                "Day15cd4_154_COV2.CON.S",            "Day15cd4_CXCR5.154_COV2.CON.S",     
                  "Bcd8_IFNg.IL2.TNFa_COV2.CON.S",      "Day15cd8_IFNg.IL2.TNFa_COV2.CON.S",  "Bcd4_IFNg_BA.4.5.S",                 "Day15cd4_IFNg_BA.4.5.S",            
                  "Bcd4_IL2_BA.4.5.S",                  "Day15cd4_IL2_BA.4.5.S",              "Bcd4_TNFa_BA.4.5.S",                 "Day15cd4_TNFa_BA.4.5.S",            
                  "Bcd8_IFNg_BA.4.5.S",                 "Day15cd8_IFNg_BA.4.5.S",             "Bcd8_IL2_BA.4.5.S",                  "Day15cd8_IL2_BA.4.5.S",             
                  "Bcd8_TNFa_BA.4.5.S",                 "Day15cd8_TNFa_BA.4.5.S",             "Bcd4_154_BA.4.5.S",                  "Day15cd4_154_BA.4.5.S",             
                  "Day15cd4_CXCR5.154_BA.4.5.S",        "Bcd8_IFNg.IL2.TNFa_BA.4.5.S",        "Day15cd8_IFNg.IL2.TNFa_BA.4.5.S")
  
  
  # Create combined new dataset which has imputed values of demographics (for phase 1 data) from dat.covar.imp AND
  # imputed values for markers (for phase 2 data) from dat.wide.v
  
  if(non_naive == FALSE){
    dat.ph1_init = dat_proc_updated %>%
      filter(ph1 == 1) %>%
      filter(naive == 1 & arm %in% trt.arms.for.analysis)
  } else if(non_naive == TRUE){
    dat.ph1_init = dat_proc_updated %>%
      filter(ph1 == 1) %>%
      filter(naive == 0 & arm %in% trt.arms.for.analysis)
  }
  
  ##############################
  if(bRiskFactors_includes_insert.stage.info){
    inputMod <- dat.ph1_init %>%
      mutate(stage = as.factor(stage),
             Trtgrp = as.factor(Trtgrp))
    
    rec <- recipe(~ stage + Trtgrp, data = inputMod)
    dummies <- rec %>%
      step_unknown(stage, Trtgrp, new_level = "missing") %>%
      step_dummy(stage, Trtgrp) %>%
      prep(training = inputMod)
    inputMod <- inputMod %>% bind_cols(bake(dummies, new_data = NULL))
    # %>%
    #   select(-c(Country, Region, CalDtEnrollIND))
    # names(inputMod)<-gsub("\\_",".",names(inputMod))
  } else {
    inputMod <- dat.ph1_init
  }
  
  #############################
  
  dat.ph1 <- inputMod %>%
    # Drop any observation with NA values in Ptid, Trt, briskfactors, endpoint and wt.D57
    drop_na(Ptid, all_of(briskfactors), all_of(endpoint), all_of(wt)) %>%
    arrange(desc(get(endpoint)))  
  
  # phase two data (all variables measured on all participants in phase 2 data)
  dat.ph2_init <- dat.ph1 %>%
    filter(ph2 == TRUE) %>%
    select(Ptid, Trt, all_of(briskfactors), all_of(endpoint), all_of(wt), any_of(markerVars)) %>%
    #drop_na(Day57bindSpike, Day57bindRBD, Day57pseudoneutid50, Day57pseudoneutid80) %>%
    arrange(desc(get(endpoint)))
  
} else if (TRIAL == "covail_xassays"){
  
  non_naive = FALSE
  #non_naive = TRUE
  
  trt.arms.for.analysis <- c(1:2, 4:12)
  # trt.arms.for.analysis <- c(1:2, 4:15)
  
  # baseline risk factors
  bRiskFactors_includes_insert.stage.info = FALSE
  
  dat_proc_updated <- dat_proc %>% 
    mutate(Delta15overBpseudoneutid50_BA.1_2fold = ifelse(Day15pseudoneutid50_BA.1 > (Bpseudoneutid50_BA.1 + log10(2)), 1, 0),
           Delta15overBpseudoneutid50_BA.1_4fold = ifelse(Day15pseudoneutid50_BA.1 > (Bpseudoneutid50_BA.1 + log10(4)), 1, 0)) %>%
    mutate(Trt = 1,   # defined to sync with code below!
           stage = case_when(stage == 1 ~ "Moderna",
                             stage == 2 ~ "Pfizer.1",
                             stage == 3 ~ "Sanofi",
                             stage == 4 ~ "Pfizer.2"),
           Trtgrp = case_when(Trtgrp == "Beta with no Omicron" ~ "Beta.no.Omicron",
                              Trtgrp == "monovalent Prototype insert only" ~ "mono.Proto.insert.only",
                              Trtgrp == "Omicron" ~ "Omicron")) 
  
  # Identify the endpoint variable
  if(COR == "D15to91covail_xassays"){ 
    endpoint <- "COVIDIndD22toD91"
  } else if(COR == "D15to181covail_xassays"){
    endpoint <- "COVIDIndD22toD181"
  }
  
  
  if(non_naive == TRUE){
    briskfactors <- c("risk_score")
    briskfactors_correction <- "Y ~ x + X$risk_score"
  } else if(bRiskFactors_includes_insert.stage.info == FALSE){
    briskfactors <- c("risk_score", "FOIstandardized")
    briskfactors_correction <- "Y ~ x + X$risk_score + X$FOIstandardized"
  } else if (bRiskFactors_includes_insert.stage.info == TRUE){
    briskfactors <- c("risk_score", "FOIstandardized", "stage_Pfizer.1", "stage_Sanofi", "Trtgrp_mono.Proto.insert.only", "Trtgrp_Omicron")
    briskfactors_correction <- "Y ~ x + X$risk_score + X$FOIstandardized + X$stage_Pfizer.1 + X$stage_Sanofi +X$Trtgrp_mono.Proto.insert.only + X$Trtgrp_Omicron"
  } 
  
  wt <- "wt.D15.tcell"            # already set by COR
  ptidvar <- "Ptid"
  
  individualMarkers <- c("Bpseudoneutid50_BA.1",
                         "Day15pseudoneutid50_BA.1",
                         "Bfrnt50_BA.1",
                         "Day15frnt50_BA.1",
                         "Bfrnt80_BA.1",
                         "Day15frnt80_BA.1",
                         primary,
                         secondary, 
                         exploratory)
  
  # markers of interest
  markerVars <- c("Bpseudoneutid50_BA.1",               "Day15pseudoneutid50_BA.1",           "Delta15overBpseudoneutid50_BA.1_2fold",   "Delta15overBpseudoneutid50_BA.1_4fold", 
                  #"Bpseudoneutid50_BA.1cat",            "Day15pseudoneutid50_BA.1cat",        "Delta15overBpseudoneutid50_BA.1cat", 
                  "Bcd4_IFNg.IL2_Wuhan.N",              "Bcd4_FS_Wuhan.N",                   
                  "Bcd4_IFNg.IL2_BA.4.5.S",             "Day15cd4_IFNg.IL2_BA.4.5.S",         "Bcd4_FS_BA.4.5.S",                   "Day15cd4_FS_BA.4.5.S",              
                  "Bcd4_IFNg.IL2.154_Wuhan.N",          "Bcd8_IFNg.IL2_BA.4.5.S",             "Day15cd8_IFNg.IL2_BA.4.5.S",         "Bcd4_IFNg.IL2.154_BA.4.5.S",        
                  "Day15cd4_IFNg.IL2.154_BA.4.5.S",     "Day15cd4_IL21_BA.4.5.S",             "Bcd4_154_Wuhan.N",                   "Bcd4_IFNg_Wuhan.N",                 
                  "Bcd4_IL2_Wuhan.N",                   "Day15cd4_154_Wuhan.N",               "Day15cd4_IFNg_Wuhan.N",              "Day15cd4_IFNg.IL2_Wuhan.N",         
                  "Day15cd4_IFNg.IL2.154_Wuhan.N",      "Day15cd4_IL2_Wuhan.N",               "Bcd4_IFNg.IL2_COV2.CON.S",           "Day15cd4_IFNg.IL2_COV2.CON.S",      
                  "Bcd8_IFNg.IL2_COV2.CON.S",           "Day15cd8_IFNg.IL2_COV2.CON.S",       "Bcd4_IFNg.IL2.154_COV2.CON.S",       "Day15cd4_IFNg.IL2.154_COV2.CON.S",  
                  "Day15cd4_IL21_COV2.CON.S",           "Bcd4_IFNg_COV2.CON.S",               "Day15cd4_IFNg_COV2.CON.S",           "Bcd4_IL2_COV2.CON.S",               
                  "Day15cd4_IL2_COV2.CON.S",            "Bcd4_TNFa_COV2.CON.S",               "Day15cd4_TNFa_COV2.CON.S",           "Bcd8_IFNg_COV2.CON.S",              
                  "Day15cd8_IFNg_COV2.CON.S",           "Bcd8_IL2_COV2.CON.S",                "Day15cd8_IL2_COV2.CON.S",            "Bcd8_TNFa_COV2.CON.S",              
                  "Day15cd8_TNFa_COV2.CON.S",           "Bcd4_154_COV2.CON.S",                "Day15cd4_154_COV2.CON.S",            "Day15cd4_CXCR5.154_COV2.CON.S",     
                  "Bcd8_IFNg.IL2.TNFa_COV2.CON.S",      "Day15cd8_IFNg.IL2.TNFa_COV2.CON.S",  "Bcd4_IFNg_BA.4.5.S",                 "Day15cd4_IFNg_BA.4.5.S",            
                  "Bcd4_IL2_BA.4.5.S",                  "Day15cd4_IL2_BA.4.5.S",              "Bcd4_TNFa_BA.4.5.S",                 "Day15cd4_TNFa_BA.4.5.S",            
                  "Bcd8_IFNg_BA.4.5.S",                 "Day15cd8_IFNg_BA.4.5.S",             "Bcd8_IL2_BA.4.5.S",                  "Day15cd8_IL2_BA.4.5.S",             
                  "Bcd8_TNFa_BA.4.5.S",                 "Day15cd8_TNFa_BA.4.5.S",             "Bcd4_154_BA.4.5.S",                  "Day15cd4_154_BA.4.5.S",             
                  "Day15cd4_CXCR5.154_BA.4.5.S",        "Bcd8_IFNg.IL2.TNFa_BA.4.5.S",        "Day15cd8_IFNg.IL2.TNFa_BA.4.5.S")
  
  
  # Create combined new dataset which has imputed values of demographics (for phase 1 data) from dat.covar.imp AND
  # imputed values for markers (for phase 2 data) from dat.wide.v
  
  if(non_naive == FALSE){
    dat.ph1_init = dat_proc_updated %>%
      filter(ph1 == 1) %>%
      filter(naive == 1 & arm %in% trt.arms.for.analysis)
  } else if(non_naive == TRUE){
    dat.ph1_init = dat_proc_updated %>%
      filter(ph1 == 1) %>%
      filter(naive == 0 & arm %in% trt.arms.for.analysis)
  }
  
  ##############################
  if(bRiskFactors_includes_insert.stage.info){
    inputMod <- dat.ph1_init %>%
      mutate(stage = as.factor(stage),
             Trtgrp = as.factor(Trtgrp))
    
    rec <- recipe(~ stage + Trtgrp, data = inputMod)
    dummies <- rec %>%
      step_unknown(stage, Trtgrp, new_level = "missing") %>%
      step_dummy(stage, Trtgrp) %>%
      prep(training = inputMod)
    inputMod <- inputMod %>% bind_cols(bake(dummies, new_data = NULL))
    # %>%
    #   select(-c(Country, Region, CalDtEnrollIND))
    # names(inputMod)<-gsub("\\_",".",names(inputMod))
  } else {
    inputMod <- dat.ph1_init
  }
  
  #############################
  
  dat.ph1 <- inputMod %>%
    # Drop any observation with NA values in Ptid, Trt, briskfactors, endpoint and wt.D57
    drop_na(Ptid, all_of(briskfactors), all_of(endpoint), all_of(wt)) %>%
    arrange(desc(get(endpoint)))  
  
  # phase two data (all variables measured on all participants in phase 2 data)
  dat.ph2_init <- dat.ph1 %>%
    filter(ph2 == TRUE) %>%
    select(Ptid, Trt, all_of(briskfactors), all_of(endpoint), all_of(wt), any_of(markerVars)) %>%
    #drop_na(Day57bindSpike, Day57bindRBD, Day57pseudoneutid50, Day57pseudoneutid80) %>%
    arrange(desc(get(endpoint)))
  
}