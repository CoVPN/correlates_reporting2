

if(TRIAL == "covail_tcell"){
  
  non_naive = FALSE
  #non_naive = TRUE
  
  trt.arms.for.analysis <- c(1:2, 4:12)
  # trt.arms.for.analysis <- c(1:2, 4:15)
  
  # baseline risk factors
  bRiskFactors_includes_insert.stage.info = TRUE
  
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
  
  non_naive = TRUE
  bRiskFactors_includes_insert.stage.info = FALSE
  trt_arms = "1dosemRNA"
  # trt_arms = "1dose"
  
  # Setting up file path
  if(isFALSE(non_naive)){ exposureStatus = "naive" } else {exposureStatus = "nonnaive"}
  if(COR == "D15to91covail_xassays"){ caseType = "proximal" } else {caseType = "allcases"}
  if(isFALSE(bRiskFactors_includes_insert.stage.info) & isFALSE(non_naive)){ referenceVars = "bRiskFactors" } else {referenceVars = "bRiskFactors.insert.stage"}
  if(isTRUE(non_naive)){ referenceVars = "briskscore" } 
  file_path = here("output", Sys.getenv("TRIAL"), paste0(exposureStatus, "_", trt_arms, "_", caseType, "_", referenceVars))
  Sys.setenv(file_path = file_path)
  
  # Create the directory only if it does not exist
  if (!dir.exists(here(file_path))) {
    dir.create(here(file_path), recursive = TRUE)
  }
  # Also create the figs directory within it
  if (!dir.exists(here(file_path, "figs"))) {
    dir.create(here(file_path, "figs"), recursive = TRUE)
  }
  
  if(trt_arms == "1dosemRNA"){
    trt.arms.for.analysis <- c(1:2, 4:12)
  } else if (trt_arms == "1dose"){
    trt.arms.for.analysis <- c(1:2, 4:15)
  }
  
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
  
  wt <- "wt.D15.xassays"            # already set by COR
  ptidvar <- "Ptid"
  
  cols <- colnames(dat_proc_updated) # This helps in parsing marker columns
  
  individualMarkers <- c(cols[grepl("pseudoneutid50", cols) & grepl("BA.1", cols) & !grepl("Day29|Day85|Day91|Day147|Day181|Delta29", cols)],
                         cols[grepl("frnt50", cols) & grepl("BA.1", cols) & !grepl("Day29|Day85|Day91|Day147|Day181|Delta29", cols)],
                         cols[grepl("frnt80", cols) & grepl("BA.1", cols) & !grepl("Day29|Day85|Day91|Day147|Day181|Delta29", cols)],
                         # primary
                         cols[grepl("IFNg.IL2_BA.4.5.S", cols)  & !grepl("BA.4.5.S1|BA.4.5.S2|Day91|Day181|vacresp|resp", cols)], #3612 3613 3616 3617    using  match(primary.ls$naive, colnames(dat_proc_updated))
                         # secondary
                         c(cols[grepl("FS_BA.4.5.S", cols) &  !grepl("BA.4.5.S1|BA.4.5.S2|Day91|Day181", cols)],     
                           cols[grepl("Day15cd4_IFNg.IL2.154_BA.4.5.S", cols) & !grepl("Day15cd4_IFNg.IL2.154_BA.4.5.S1|Day15cd4_IFNg.IL2.154_BA.4.5.S2|resp", cols)]),
                         # exploratory
                         cols[Reduce(`|`, lapply(exploratory$naive, grepl, x = cols)) & 
                                !grepl("\\.S1|\\.S2", cols) & 
                                !grepl("resp", cols)],
                         # Fc receptor and complement variables 
                         cols[grepl("G2AH|G2AH|GR2B|G3AV|G3AF|GR3B|FCAR|FCR7|FCR6|C1Q", cols) & grepl("BA1S|WA1N", cols) & !grepl("_One|_Two|Day91|Day181", cols)],
                         # IgG and IgA and IgM binding antibody variables
                         cols[grepl("TIGG|IGG1|IGG2|IGG3|IGG4|IGA1|IGA2|IGM", cols) & grepl("BA1S|WA1N", cols) & !grepl("_One|_Two|Day91|Day181", cols)],
                         # Combination markers variable set 
                         cols[grepl("PC1", cols) & !grepl("Day29|Day91|Day181|NLPC1", cols)],
                         cols[grepl("MDW", cols) & !grepl("Day29|Day91|Day181|Delta29", cols)])
  


  # Drop out all categorical markers
  individualMarkers <- individualMarkers[!grepl("cat$|Delta", individualMarkers)]
  
  # markers of interest
  markerVars <- individualMarkers
  

  if(isFALSE(non_naive)){
    dat.ph1_init = dat_proc_updated %>%
      filter(ph1 == 1) %>%
      filter(naive == 1 & arm %in% trt.arms.for.analysis)
  } else if(isTRUE(non_naive)){
    dat.ph1_init = dat_proc_updated %>%
      filter(ph1 == 1) %>%
      filter(naive == 0 & arm %in% trt.arms.for.analysis)
  }

  ##############################
  if(isTRUE(bRiskFactors_includes_insert.stage.info)){
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
  
  if(isTRUE(non_naive)){
    briskfactors <- c("risk_score")
  } else if(isFALSE(bRiskFactors_includes_insert.stage.info)){
    briskfactors <- c("risk_score", "FOIstandardized")
  } else if (isTRUE(bRiskFactors_includes_insert.stage.info)){
    if (all(c("stage_Pfizer.1", "stage_Pfizer.2", "stage_Sanofi") %in% colnames(inputMod))) {
      briskfactors <- c("risk_score", "FOIstandardized", "stage_Pfizer.1", "stage_Pfizer.2", "stage_Sanofi", 
                        "Trtgrp_mono.Proto.insert.only", "Trtgrp_Omicron")
    } else if (all(c("stage_Pfizer.1", "stage_Sanofi") %in% colnames(inputMod))){
      briskfactors <- c("risk_score", "FOIstandardized", "stage_Pfizer.1", "stage_Sanofi", 
                        "Trtgrp_mono.Proto.insert.only", "Trtgrp_Omicron")
    } else if (all(c("stage_Pfizer.1") %in% colnames(inputMod))){
      briskfactors <- c("risk_score", "FOIstandardized", "stage_Pfizer.1", 
                        "Trtgrp_mono.Proto.insert.only", "Trtgrp_Omicron")
    }
  }
  
  briskfactors_correction <- paste0("Y ~ x + ", paste0("X$", briskfactors, collapse = " + "))
  
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