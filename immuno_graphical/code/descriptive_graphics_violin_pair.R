#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
library(GGally)
library(stringr)
library(cowplot) # for function plot_grid
library(grid)
library(gridExtra)
library(wCorr) # weighted correlation
#library(ggnewscale) # for new_scale_color() 
library(survey)
# Suppress summarise info
options(dplyr.summarise.inform = FALSE)
# For stratum with 1 ppt
options(survey.lonely.psu="adjust")
library(tidyverse)

# There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
#if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))
source(here::here("code", "process_violin_pair_functions.R"))
Sys.setenv(DESCRIPTIVE = 1)
source(here::here("..", "_common.R"))
source(here::here("code", "params.R")) # load parameters
assay_lim <- readRDS(here("data_clean", "assay_lim.rds"))

# for the order of figure panels
if ("order_in_panel" %in% colnames(assay_metadata)){
    assay_order = assay_metadata %>% dplyr::arrange(panel, order_in_panel) %>% select(assay_label_short) %>% pull()
    assay_metadata = assay_metadata %>%
        mutate(assay_label_short = factor(assay_label_short,
                                          levels = assay_order
        ))
}

dat.longer.immuno.subset.plot1 <- readRDS(here::here("data_clean", "longer_immuno_data_plot1.rds"))
if(attr(config,"config") == "nextgen_mock"){
    dat.longer.immuno.subset.plot1.trackA = dat.longer.immuno.subset.plot1$dat_stats # for Track A RIS/RIS-PBMC (initial immuno)
    dat.longer.immuno.subset.plot1.whole = dat.longer.immuno.subset.plot1$dat_stats2 # for whole RIS/RIS-PBMC (final immuno)
}
if (attr(config,"config") == "nextgen_mock") {
    dat.longer.immuno.subset.plot1.3 <- readRDS(here::here("data_clean", "longer_immuno_data_plot1.3.rds"))
    dat.longer.immuno.subset.plot1.3.whole = dat.longer.immuno.subset.plot1.3$dat_stats2 # for whole RIS/RIS-PBMC (final immuno)
}
if(attr(config,"config") == "vat08_combined"){
    dat.longer.immuno.subset.plot1_stage1_stage2 <- readRDS(here::here("data_clean", "longer_immuno_data_plot1_stage1_stage2.rds"))}
if(attr(config,"config")=="janssen_partA_VL") {
    dat.longer.immuno.subset.plot1.2 <- readRDS(here::here("data_clean", "longer_immuno_data_plot1.2.rds"))}
dat.immuno.subset.plot3 <- readRDS(here::here("data_clean", "twophase_data.rds")); dat.immuno.subset.plot3$all_one <- 1 # as a placeholder for strata values
dat.immuno.subset.plot3 <- dat.immuno.subset.plot3 %>%
    mutate(Trt = factor(Trt, levels = c(1, 0), labels = trt.labels),
           nnaive = factor(Bserostatus, levels = c(0, 1), labels = bstatus.labels))

###### Set 1 plots: Ab distributions at main timepoints and delta by vaccine/placebo, naive/nnaive
if (attr(config, "config") == "nextgen_mock") {
    set1_times <- times_[c(2, 4, 6, 8)]
} else {
    set1_times <- times_[1:5]
}
# ID50, at B, D22, D43, D22-B, D43-B
# bindSpike, at B, D22, D43, D22-B, D43-B
for (panel in if (study_name == "NextGen_Mock") {c("IgG_sera", "pseudoneutid50_sera", "IgG_nasal", "pseudoneutid50_nasal", "IgG_saliva", "pseudoneutid50_saliva", "T4", "T8")
    } else {c("pseudoneutid50", "bindSpike")}){
    # by naive/non-naive, vaccine/placebo
    
    if (any(grepl(panel, assays)) == FALSE) next
    
    if (attr(config,"config") %in% c("janssen_partA_VL", "janssen_pooled_partA", "nextGen_mock")) next # janssen_partA_VL doesn't need these plots
    
    if (attr(config,"config")=="vat08_combined" & panel=="pseudoneutid50") {
        dat.longer.immuno.subset.plot1_ = dat.longer.immuno.subset.plot1 %>% filter(ph2.immuno.nAb==1)
    } else if (attr(config,"config")=="vat08_combined" & panel=="bindSpike") {
        dat.longer.immuno.subset.plot1_ = dat.longer.immuno.subset.plot1 %>% filter(ph2.immuno.bAb==1)
    } else if (study_name != "NextGen_Mock") {
        dat.longer.immuno.subset.plot1_ = dat.longer.immuno.subset.plot1
    }
    
    for (tm_subset in c("^B|Day", if(sum(grepl("Delta", set1_times))>0) "Delta")){
        
        for (tm in c(if(attr(config, "config") == "nextgen_mock") "Day initial", "Day whole")){
            
            if (tm == "Day initial") {
                set1_times_sub = set1_times[c(2, 4)]
                dat.longer.immuno.subset.plot1_  = dat.longer.immuno.subset.plot1.trackA %>% filter(Track == "A") # NextGen_Mock has different weights for track A and whole, for bAb/nAb and ICS
            } else if (tm == "Day whole" & attr(config, "config") == "nextgen_mock") {
                set1_times_sub = set1_times[c(1, 3)]
                dat.longer.immuno.subset.plot1_  = dat.longer.immuno.subset.plot1.whole # NextGen_Mock has different weights for track A and whole, for bAb/nAb and ICS
            } else {set1_times_sub = set1_times[grepl(tm_subset, set1_times)]}
            
            # set time to be a factor
            dat.longer.immuno.subset.plot1_ <- dat.longer.immuno.subset.plot1_ %>% 
                mutate(time = factor(time, levels = times_))
            
            f_1 <- f_by_time_assay(
                dat = dat.longer.immuno.subset.plot1_ %>% mutate(x="1"),
                assays = assays[grepl(panel, assays)],
                times = set1_times_sub,
                ylim = c(floor(min(assay_lim[assays[grepl(panel, assays)], set1_times_sub, "lb"], na.rm = T)),
                         ceiling(max(assay_lim[assays[grepl(panel, assays)], set1_times_sub, "ub"], na.rm = T))),#if (grepl("^B|Day", tm_subset) & panel=="bindSpike") {c(2, 7)
                    #} else if (grepl("^B|Day", tm_subset) & panel=="pseudoneutid50") {c(1, 6.5)
                     #   } else if (grepl("Delta", tm_subset)) {c(-3, 4.2)} else if (study_name == "NextGen_Mock") {c(-3, 7)} else {c(-3, 4.5)},
                axis.x.text.size = 20,
                strip.x.text.size = ifelse(study_name == "NextGen_Mock", 25, 10),
                strip.y.text.size = ifelse(grepl("pseudoneutid50$", panel), 25, ifelse(grepl("pseudoneutid50", panel), 15, ifelse(grepl("T4|T8", panel), 18, ifelse(grepl("IgG|IgA", panel), 8, 10)))),
                panel.text.size = ifelse(grepl("pseudoneutid50$", panel), 7, ifelse(grepl("T4|T8|pseudoneutid50", panel), 6, ifelse(grepl("IgG|IgA", panel), 10, 4.5))),
                facet.y.var = if (study_name == "NextGen_Mock") vars(assay_label2) else vars(Trt),
                facet.x.var = if (study_name == "NextGen_Mock") vars(Trt) else vars(assay_label2),
                label_format = ifelse(panel %in% c("T4", "T8"), "percent", "log10"), 
                color.map = c("Investigational Vaccine" = "#1749FF", "Comparator Vaccine" = "#009E73"))
            
            for (i in 1:length(set1_times_sub)){
                
                if (grepl("IgG|IgA", panel) & study_name == "NextGen_Mock") {p_wrapped <- ggdraw(f_1[[i]]) + theme(plot.margin = unit(c(0, 20, 0, 20), "cm"))
                } else if (grepl("pseudoneutid50", panel) & study_name == "NextGen_Mock") {p_wrapped <- ggdraw(f_1[[i]]) + theme(plot.margin = unit(c(0, 5, 0, 5), "cm"))
                } else {p_wrapped <- ggdraw(f_1[[i]])}
                
                file_name <- paste0("/", panel, "_at_", set1_times_sub[i], 
                                    ifelse(study_name == "NextGen_Mock" & tm == "Day whole", "_final", 
                                           ifelse(study_name == "NextGen_Mock" & tm == "Day initial", "_initial", "")), ".pdf")
                ggsave(plot = p_wrapped, filename = paste0(save.results.to, file_name), 
                       width = ifelse(grepl("IgG|IgA", panel) & study_name == "NextGen_Mock", 40, 
                                      ifelse(grepl("pseudoneutid50", panel) & study_name == "NextGen_Mock", 22, 16)), 
                       height = ifelse(grepl("IgG|IgA", panel) & study_name == "NextGen_Mock", 40, 
                                       ifelse(grepl("pseudoneutid50", panel) & study_name == "NextGen_Mock", 22, 16)))
            }
        }
    }
}

# pooling two arms at D01 for NextGen_Mock
if (study_name == "NextGen_Mock" & study_name != "NextGen_Mock") {
    for (panel in c("IgG_sera", "pseudoneutid50_sera", "IgG_nasal", "pseudoneutid50_nasal", "IgG_saliva", "pseudoneutid50_saliva", "T4", "T8")) {
        
        if (any(grepl(panel, assays)) == FALSE) next
            
        f_1 <- f_by_time_assay(
            dat = dat.longer.immuno.subset.plot1.3.whole %>% mutate(x="1"),
            assays = assays[grepl(panel, assays)],
            times = "B",
            ylim = c(floor(min(assay_lim[assays[grepl(panel, assays)], "B", "lb"], na.rm = T)),
                     ceiling(max(assay_lim[assays[grepl(panel, assays)], "B", "ub"], na.rm = T))),
            axis.x.text.size = 20,
            strip.x.text.size = ifelse(study_name == "NextGen_Mock", 25, 10),
            strip.y.text.size = ifelse(grepl("pseudoneutid50$", panel), 25, ifelse(grepl("pseudoneutid50", panel), 15, ifelse(grepl("T4|T8", panel), 18, ifelse(grepl("IgG|IgA", panel), 8, 10)))),
            panel.text.size = ifelse(grepl("pseudoneutid50$", panel), 7, ifelse(grepl("T4|T8|pseudoneutid50", panel), 6, ifelse(grepl("IgG|IgA", panel), 10, 4.5))),
            facet.y.var = vars(assay_label2),
            facet.x.var = vars(Trt),
            label_format = ifelse(panel %in% c("T4", "T8"), "percent", "log10"),
            color.map = c("Pooled Arm" = "#FF6F1B"))
        
        for (i in 1){
            
            if (grepl("IgG|IgA", panel) & study_name == "NextGen_Mock") {p_wrapped <- ggdraw(f_1[[i]]) + theme(plot.margin = unit(c(0, 30, 0, 30), "cm"))
            } else if (grepl("pseudoneutid50", panel) & study_name == "NextGen_Mock") {p_wrapped <- ggdraw(f_1[[i]]) + theme(plot.margin = unit(c(0, 10, 0, 10), "cm"))
            } else {p_wrapped <- ggdraw(f_1[[i]])}
            
            file_name <- paste0("/", panel, "_at_B_final.pdf")
            ggsave(plot = p_wrapped, filename = paste0(save.results.to, file_name), 
                   width = ifelse(grepl("IgG|IgA", panel) & study_name == "NextGen_Mock", 40, 
                                  ifelse(grepl("pseudoneutid50", panel) & study_name == "NextGen_Mock", 22, 16)), 
                   height = ifelse(grepl("IgG|IgA", panel) & study_name == "NextGen_Mock", 40, 
                                   ifelse(grepl("pseudoneutid50", panel) & study_name == "NextGen_Mock", 22, 16)))
        }
    }
}

# just for janssen_partA_VL, loop through panel and region for vaccine at day 29
if (attr(config,"config")=="janssen_partA_VL") {
    for (panel in c("bind", "pseudo")){
        
        panel_lb = case_when(panel=="bind" ~ "bAb",
                             panel=="pseudo" ~ "nAb")
        
        for (region in c(1, 2)){
            
            region_lb = case_when(region==0 ~ "NAM",
                                  region==1 ~ "LATAM",
                                  region==2 ~ "ZA")
            region_lb_long = case_when(region==0 ~ "Northern America",
                                       region==1 ~ "Latin America",
                                       region==2 ~ "Southern Africa")
            
            if (panel=="bind") {assays_ = subset(assays, grepl(panel, assays))
            } else if (panel %in% c("pseudo") & region == 1) { # Latin America = 1
                selected = assays[grepl("Reference|Zeta|Mu|Gamma|Lambda", assay_metadata$assay_label_short)]
                assays_ = subset(selected, grepl(panel, selected))
            } else if (panel %in% c("pseudo") & region == 2) { # Southern Africa = 2
                selected = assays[grepl("Reference|Beta|Delta", assay_metadata$assay_label_short)]
                assays_ = subset(selected, grepl(panel, selected))
            }
            
            f_1 <- f_by_time_assay(
                dat = dat.longer.immuno.subset.plot1.2 %>% filter(Trt=="Vaccine" & Bserostatus=="0" & Region==region) %>% mutate(x="1"),
                assays = assays_,
                times = "Day29",
                ylim = c(0, 4),
                rate.pos = 3.5, 
                axis.x.text.size = 20,
                strip.x.text.size = ifelse(panel=="bind" | region==1, 8, 16),
                panel.text.size = 7,
                facet.y.var = vars(Trt_nnaive),
                facet.x.var = vars(assay_label_short))
                
                file_name <- paste0("/", panel_lb, "_at_Day29_vaccine_", region_lb, ".pdf")
                ggsave(plot = f_1[[1]], filename = paste0(save.results.to, file_name), width = 16, height = 16)
        }
    }
}

###### Set 2 plots: Longitudinal plots D1, D22, D43
set2_assays = assays
# bindSpike
# ID50
# bindSpike mdw
for (panel in if (study_name == "NextGen_Mock") {
     c("bindSpike_IgG_sera$|bindSpike_IgG_sera_delta_AY.4$",
       "bindSpike_IgG_sera_BA.5$|bindSpike_IgG_sera_BA.2.86$",
       "bindSpike_IgG_sera_XBB.1.5$|bindSpike_IgG_sera_JN.1$",
       "bindSpike_IgG_sera_KP.2$|bindSpike_IgG_sera_KP.3$",
       "bindSpike_IgG_sera_LB.1$|bindN_IgG_sera$",
       
       "bindSpike_IgA_sera$|bindSpike_IgA_sera_delta_AY.4$",
       "bindSpike_IgA_sera_BA.5$|bindSpike_IgA_sera_BA.2.86$",
       "bindSpike_IgA_sera_XBB.1.5$|bindSpike_IgA_sera_JN.1$",
       "bindSpike_IgA_sera_KP.2$|bindSpike_IgA_sera_KP.3$",
       "bindSpike_IgA_sera_LB.1$",
       
       "pseudoneutid50_sera_KP.2$|pseudoneutid50_sera_XBB.1.5$", 
       
       "bindSpike_IgG_nasal$|bindSpike_IgG_nasal_delta_AY.4$",
       "bindSpike_IgG_nasal_BA.5$|bindSpike_IgG_nasal_BA.2.86$",
       "bindSpike_IgG_nasal_XBB.1.5$|bindSpike_IgG_nasal_JN.1$",
       "bindSpike_IgG_nasal_KP.2$|bindSpike_IgG_nasal_KP.3$",
       "bindSpike_IgG_nasal_LB.1$|bindN_IgG_nasal$",
       
       "bindSpike_IgA_nasal$|bindSpike_IgA_nasal_delta_AY.4$",
       "bindSpike_IgA_nasal_BA.5$|bindSpike_IgA_nasal_BA.2.86$",
       "bindSpike_IgA_nasal_XBB.1.5$|bindSpike_IgA_nasal_JN.1$",
       "bindSpike_IgA_nasal_KP.2$|bindSpike_IgA_nasal_KP.3$",
       "bindSpike_IgA_nasal_LB.1$",
       
       "pseudoneutid50_nasal_KP.2$|pseudoneutid50_nasal_XBB.1.5$",
       
       "bindSpike_IgG_saliva$|bindSpike_IgG_saliva_delta_AY.4$",
       "bindSpike_IgG_saliva_BA.5$|bindSpike_IgG_saliva_BA.2.86$",
       "bindSpike_IgG_saliva_XBB.1.5$|bindSpike_IgG_saliva_JN.1$",
       "bindSpike_IgG_saliva_KP.2$|bindSpike_IgG_saliva_KP.3$",
       "bindSpike_IgG_saliva_LB.1$|bindN_IgG_saliva$",
       
       "bindSpike_IgA_saliva$|bindSpike_IgA_saliva_delta_AY.4$",
       "bindSpike_IgA_saliva_BA.5$|bindSpike_IgA_saliva_BA.2.86$",
       "bindSpike_IgA_saliva_XBB.1.5$|bindSpike_IgA_saliva_JN.1$",
       "bindSpike_IgA_saliva_KP.2$|bindSpike_IgA_saliva_KP.3$",
       "bindSpike_IgA_saliva_LB.1$",
       
       "pseudoneutid50_saliva_KP.2$|pseudoneutid50_saliva_XBB.1.5$",
       
       "T4_IFNg_OR_IL2_N_Index$",
       
       "T4_IFNg_OR_IL2_Spike_KP.2$|T4_IFNg_OR_IL2_Spike_Index.D614$",
       
       "T8_IFNg_OR_IL2_N_Index$",
       
       "T8_IFNg_OR_IL2_Spike_KP.2$|T8_IFNg_OR_IL2_Spike_Index.D614$"
       )
    } else {c("pseudoneutid50", "bindSpike")}){
    # by naive/non-naive, vaccine/placebo
    
    strs <- gsub("\\$$", "", unlist(strsplit(panel, "\\|")))
    if (any(strs %in% assays) == FALSE) next
    
    if (attr(config,"config") %in% c("janssen_partA_VL", "janssen_pooled_partA")) next # janssen_partA_VL doesn't need these plots
    
    if (attr(config,"config")=="vat08_combined" & panel=="pseudoneutid50") {
        dat.longer.immuno.subset.plot1_ = dat.longer.immuno.subset.plot1 %>% filter(ph2.immuno.nAb==1)
    } else if (attr(config,"config")=="vat08_combined" & panel=="bindSpike") {
        dat.longer.immuno.subset.plot1_ = dat.longer.immuno.subset.plot1 %>% filter(ph2.immuno.bAb==1)
    } else if (attr(config,"config") == "nextgen_mock") {
        dat.longer.immuno.subset.plot1_  = dat.longer.immuno.subset.plot1.trackA %>%  # NextGen_Mock has different weights for track A and whole, for bAb/nAb and ICS
            mutate(time = factor(time, levels = c("B", "Day31", "Day91", "Day181", "Day366")))
    } else {
        dat.longer.immuno.subset.plot1_ = dat.longer.immuno.subset.plot1
    }
    
    
    f_2 <- f_longitude_by_assay(
        dat = dat.longer.immuno.subset.plot1_,
        x.var = "time",
        x.lb = if (study_name == "NextGen_Mock") {c("D01","D31","D91","D181","D366")} else if (study_name == "VAT08"){c("D1","D22","D43")},
        assays = if (study_name == "NextGen_Mock") {set2_assays[grepl(panel, set2_assays)]} else {set2_assays[grepl(panel, set2_assays) & !grepl("mdw", set2_assays)]},
        ylim = c(floor(min(assay_lim[assays[grepl(panel, set2_assays)], , "lb"], na.rm = T)),
                 ceiling(max(assay_lim[assays[grepl(panel, set2_assays)], , "ub"], na.rm = T))),#if (grepl("T4|T8", panel)) {c(-2.2, 2.5)} else {c(1, 6.5)},
        ybreaks = if (study_name == "NextGen_Mock") {seq(-2, 7, 1)} else {c(0,2,4,6)},
        times = if (study_name == "NextGen_Mock") {c("B","Day31","Day91","Day181","Day366")} else if (study_name == "VAT08"){c("B","Day22","Day43")},
        strip.text.x.size = ifelse(panel=="pseudoneutid50", 25, ifelse(study_name == "NextGen_Mock", 25, 12)),
        strip.text.y.size = ifelse(study_name == "NextGen_Mock", 13, 25), 
        panel.text.size = ifelse(panel=="pseudoneutid50", 6, 4),
        facet.y.var = if (study_name == "NextGen_Mock") vars(assay_label_short) else vars(Trt), 
        facet.x.var = if (study_name == "NextGen_Mock") vars(Trt) else vars(assay_label_short),
        y.axis.lb = ifelse(study_name == "NextGen_Mock", " ", ""),
        label_format = ifelse(grepl("T4|T8", panel), "percent", "log10"),
        color.map = c("Investigational Vaccine" = "#1749FF", "Comparator Vaccine" = "#009E73")
    )
    
    file_name <- paste0("/", ifelse(panel=="pseudoneutid50", "nAb", ifelse(panel=="bindSpike", "bAb", gsub("\\$", "", gsub("\\|", "__", panel)))), "_longitudinal", if (study_name == "NextGen_Mock") "_initial", ".pdf")
    ggsave(plot = f_2, filename = paste0(save.results.to, file_name), width = 16, height = 11)
    
}

if (attr(config,"config") == "vat08_combined"){
    for (panel in c("pseudoneutid50", "bindSpike")){
        # by naive/non-naive, vaccine/placebo
        
        if (attr(config,"config") %in% c("janssen_partA_VL", "janssen_pooled_partA")) next # janssen_partA_VL doesn't need these plots
        
        if (attr(config,"config")=="vat08_combined" & panel=="pseudoneutid50") {
            dat.longer.immuno.subset.plot1_ = dat.longer.immuno.subset.plot1 %>% filter(ph2.immuno.nAb==1)
        } else if (attr(config,"config")=="vat08_combined" & panel=="bindSpike") {
            dat.longer.immuno.subset.plot1_ = dat.longer.immuno.subset.plot1 %>% filter(ph2.immuno.bAb==1)
        } else {
            dat.longer.immuno.subset.plot1_ = dat.longer.immuno.subset.plot1
        }
        
        f_2 <- f_longitude_by_assay(
            dat = dat.longer.immuno.subset.plot1_,
            x.var = "time",
            x.lb = c("D1","D22","D43"),
            assays = set2_assays[grepl(panel, set2_assays) & grepl("mdw", set2_assays)],
            ylim = c(1, 6), #c(0, 5.2),
            times = c("B","Day22","Day43"),
            strip.text.x.size = ifelse(panel=="pseudoneutid50", 25, 12),
            panel.text.size = ifelse(panel=="pseudoneutid50", 6, 4),
            facet.y.var = vars(Trt_nnaive), 
            facet.x.var = vars(assay_label_short)
        )
        
        file_name <- paste0("/", ifelse(panel=="pseudoneutid50", "nAb", ifelse(panel=="bindSpike", "bAb", "")), "_mdw_longitudinal.pdf")
        ggsave(plot = f_2, filename = paste0(save.results.to, file_name), width = 11, height = 11)
    }
}

# print this adhoc figures only in sanofi stage 1 report, which includes both stage 1 and stage 2 data
if (attr(config,"config")=="vat08_combined") {
    if (dat.longer.immuno.subset.plot1_$Trialstage[1]==1) {
        # longitudinal plots for stage 1 and stage 2, non-naive ppt
        
        for (asy in set2_assays){
            
            if (attr(config,"config")=="vat08_combined" & grepl("bind", asy)) {
                subdat_adhoc = dat.longer.immuno.subset.plot1_stage1_stage2 %>% filter(ph2.immuno.bAb==1)
            } else if (attr(config,"config")=="vat08_combined" & grepl("pseudo", asy)){
                subdat_adhoc = dat.longer.immuno.subset.plot1_stage1_stage2 %>% filter(ph2.immuno.nAb==1)
            }
            
            f_2 <- f_longitude_by_assay(
                dat = subdat_adhoc %>% 
                    filter(Bserostatus == 1) %>%
                    mutate(Trt_nnaive = factor(paste(Trt, nnaive),
                                               levels = c("Vaccine Non-naive", "Placebo Non-naive"),
                                               labels = c("Vaccine\nNon-naive", "Placebo\nNon-naive")),
                           Trialstage = ifelse(Trialstage == 1, "Stage 1", "Stage 2"),
                           time = factor(time, levels = c("B","Day22","Day43","Day78","Day134","Day202","Day292","Day387"))
                    ),
                x.var = "time",
                x.lb = c("D1","D22","D43","D78","D134","D202","D292","D387"),
                assays = asy,
                ylim = if (grepl("bindSpike", asy)) {c(2, 7)} else if (grepl("pseudoneutid50", asy)) {c(1, 6.5)},
                times = c("B","Day22","Day43","Day78","Day134","Day202","Day292","Day387"),
                strip.text.x.size = 14,
                panel.text.size = 6,
                axis.text.x.size = 14,
                facet.y.var = vars(Trt_nnaive), 
                facet.x.var = vars(Trialstage),
                y.axis.lb = gsub(" \\(AU/ml\\)", "", labels.assays.short[asy])
            )
            
            file_name <- paste0("/", asy, "_longitudinal_nonnaive_stage1stage2.pdf")
            ggsave(plot = f_2, filename = paste0(save.results.to, file_name), width = 16, height = 11)
        }
        
        
        # Generate a full table with all the estimates: GMTR
        # non-naive ppts only
        # Vaccine vs. Placebo, by Stage 1 and 2
        
        sub.by <- c("all")
        ds.i.bAb <- filter(dat_proc %>%
                           mutate(all = 1) %>% # fake column to satisfy sub.by
                           filter(Bserostatus == 1) %>% # only do this for non-naive ppts
                           mutate(Trt=ifelse(Trt==1, "Vaccine", "Placebo")), ph1.immuno.bAb==1)
        ds.i.nAb <- filter(dat_proc %>%
                               mutate(all = 1) %>% # fake column to satisfy sub.by
                               filter(Bserostatus == 1) %>% # only do this for non-naive ppts
                               mutate(Trt=ifelse(Trt==1, "Vaccine", "Placebo")), ph1.immuno.nAb==1)
        gm.v <- apply(expand.grid(times_[!grepl("Delta",times_)], assays), 1, paste0, collapse="")
        
        subs <- "Trt"
        comp.i <- c("Vaccine", "Placebo")
        
        rgmt_stage1_bAb_vac <- get_rgmt(ds.i.bAb %>% filter(Trialstage == 1), 
                                        gm.v[grepl("bindSpike", gm.v)], subs, comp_lev=comp.i, sub.by, strata="Wstratum", 
                                        weights="wt.immuno.bAb", subset="ph2.immuno.bAb") %>% mutate(Trialstage = "Stage 1")
        rgmt_stage1_nAb_vac <- get_rgmt(ds.i.nAb %>% filter(Trialstage == 1), 
                                        gm.v[grepl("pseudoneutid50", gm.v)][c(1:5,9:11,17:19,25:27,33:35,41:43)], subs, comp_lev=comp.i, sub.by, strata="Wstratum", 
                                        weights="wt.immuno.nAb", subset="ph2.immuno.nAb") %>% mutate(Trialstage = "Stage 1") 
        rgmt_stage2_bAb_vac <- get_rgmt(ds.i.bAb %>% filter(Trialstage == 2), 
                                        gm.v[grepl("bindSpike", gm.v)], subs, comp_lev=comp.i, sub.by, strata="Wstratum", 
                                        weights="wt.immuno.bAb", subset="ph2.immuno.bAb") %>% mutate(Trialstage = "Stage 2") 
        rgmt_stage2_nAb_vac <- get_rgmt(ds.i.nAb %>% filter(Trialstage == 2), 
                                        gm.v[grepl("pseudoneutid50", gm.v)], subs, comp_lev=comp.i, sub.by, strata="Wstratum", 
                                        weights="wt.immuno.nAb", subset="ph2.immuno.nAb") %>% mutate(Trialstage = "Stage 2")
        for (asy in set2_assays){
            f_2 <- f_longitude_by_assay(
                dat = rgmt_stage1_bAb_vac %>% 
                    bind_rows(rgmt_stage1_nAb_vac) %>%
                    bind_rows(rgmt_stage2_bAb_vac) %>%
                    bind_rows(rgmt_stage2_nAb_vac) %>%
                    filter(!grepl("Delta", mag_cat)) %>%
                    mutate(value = geomean, 
                           Ptid = Trialstage,
                           RespRate = "",
                           time = gsub(paste0(assays, collapse="|"), "", mag_cat),
                           time = factor(time, levels = c("B","Day22","Day43","Day78","Day134","Day202","Day292","Day387")),
                           assay = gsub(paste0("^", times_, collapse="|"), "", mag_cat),
                           nnaive = "Non-naive",
                           lb = "",
                           lbval = -99,
                           lb2 = "",
                           lbval2 = -99),
                x.var = "time",
                x.lb = c("D1","D22","D43","D78","D134","D202","D292","D387"),
                assays = asy,
                ylim = if (grepl("bindSpike", asy)) {c(0, 60)} else if (grepl("pseudoneutid50", asy)) {c(0, 155)},
                ybreaks = if (grepl("bindSpike", asy)) {seq(0, 60, 10)} else if (grepl("pseudoneutid50", asy)) {seq(0, 150, 20)},
                times = c("B","Day22","Day43","Day78","Day134","Day202","Day292","Day387"),
                strip.text.x.size = 14,
                panel.text.size = 6,
                axis.text.x.size = 14,
                facet.y.var = vars(nnaive), 
                facet.x.var = vars(Trialstage),
                y.axis.lb = paste0("Geometric Mean Ratio (Vaccine/Placebo) of\n", gsub(" \\(AU/ml\\)", "", labels.assays.short[asy])),
                y.lb.scale = "original"
            )
            
            file_name <- paste0("/", asy, "_longitudinal_gmr_nonnaive_stage1stage2.pdf")
            ggsave(plot = f_2, filename = paste0(save.results.to, file_name), width = 16, height = 11)
        }
    }
}

# adhoc figure needed for the severe manuscript reviewer
if (attr(config,"config") == "janssen_pooled_partA"){
    assay_metadata$assay_label = labels.assays[match(assay_metadata$assay, names(labels.assays))]
    assay_metadata$assay_label = case_when(assay_metadata$assay_label=="Binding Antibody to Spike" ~ "Binding Antibody to Spike\nBAU/mL",
                                           assay_metadata$assay_label=="Binding Antibody to RBD" ~ "Binding Antibody to RBD\nBAU/mL",
                                           assay_metadata$assay_label=="PsV Neutralization 50% Titer" ~ "PsV Neutralization 50% Titer\nIU/mL",
                                           TRUE ~ assay_metadata$assay_label)
    panel=""
    
    dat_plot_ = subset(dat.longer.immuno.subset.plot1, 
                      Trt=="Vaccine" & nnaive=="Baseline Neg" & time %in% c("Day29","Day71") & assay %in% c("bindSpike","bindRBD","pseudoneutid50")) %>% 
        filter(!is.na(value))
    # only 104 ptids have id50 data at Day 71
    ids_id50_day71 <- subset(dat_plot_, time=="Day71" & assay=="pseudoneutid50")$Ptid %>% unique()
    stopifnot(length(ids_id50_day71) == 104)
    
    # 2 of 104 do not have day 71 bindRBD and bindSpike
    ids_day71 <- ids_id50_day71[ids_id50_day71 %in% subset(dat_plot_, time=="Day71" & assay=="bindRBD")$Ptid]
    stopifnot(length(ids_day71) == 102)
    
    all(ids_day71 %in% subset(dat_plot_, time=="Day29" & assay=="pseudoneutid50")$Ptid)
    all(ids_day71 %in% subset(dat_plot_, time=="Day29" & assay=="bindRBD")$Ptid)
    all(ids_day71 %in% subset(dat_plot_, time=="Day29" & assay=="bindSpike")$Ptid)
    all(ids_day71 %in% subset(dat_plot_, time=="Day71" & assay=="bindRBD")$Ptid)
    all(ids_day71 %in% subset(dat_plot_, time=="Day71" & assay=="bindSpike")$Ptid)

    set.seed(20240320)
    random25 <- sample(ids_day71, 25)
    dat_plot = dat_plot_ %>% filter(Ptid %in% random25) %>% mutate(RespRate=" ")
    
    f_2 <- f_longitude_by_assay(
        dat = dat_plot,
        x.var = "time",
        x.lb = c("D29","D71"),
        assays = c("bindSpike","bindRBD","pseudoneutid50"),
        ylim = c(0, 3),
        times = c("Day29","Day71"),
        strip.text.x.size = 25,
        panel.text.size = 6,
        facet.y.var = NULL, 
        facet.x.var = vars(assay_label)
    )
    
    file_name <- paste0("/adhoc_longitudinal_vaccine_baselineneg_random25.pdf")
    ggsave(plot = f_2, filename = paste0(save.results.to, file_name), width = 16, height = 11)
}
#f_2 <- f_longitude_by_assay(
#    dat = dat.longer.immuno.subset.plot1,
#    x.var = "time",
#    x.lb = c("D1","D22","D43"),
#    assays = "bindSpike_mdw",
#    ylim = c(-6, 11),
#    times = c("B","Day22","Day43"),
#    strip.text.x.size = 25,
#    strip.text.y.size = 22,
#    panel.text.size = 6,
#    facet.y.var = vars(Trt_nnaive), 
#    facet.x.var = vars(assay_label_short)
#)

#file_name <- paste0("/crossreactivity_longitudinal.pdf")
#ggsave(plot = f_2, filename = paste0(save.results.to, file_name), width = 6, height = 11)


###### Set 3 plots: Correlation plots across markers at a given time point
# 15 markers, non-naive pooled vaccine & placebo, three timepoints
# 15 markers, naive vaccine, three timepoints
for (grp in c("non_naive_vac_pla", 
              if (study_name == "NextGen_Mock") "non_naive_vac", 
              if (study_name == "NextGen_Mock") "non_naive_pla", 
              if (study_name != "NextGen_Mock") "naive_vac")){
    
    if (attr(config,"config") %in% c("janssen_partA_VL", "janssen_pooled_partA")) next # janssen_partA_VL doesn't need these plots
    
    for (t in if (study_name == "NextGen_Mock") {c("B", "Day31", "Day181")} else {c("B","Day22","Day43")}) {
        
        if (grp == "naive_vac" && t=="B") next # this is not needed for VAT08
        if (grp == "non_naive_vac_pla" & study_name == "NextGen_Mock" & t != "B") next # pooled results are only needed at baseline for NextGen_Mock
        if (grp != "non_naive_vac_pla" & study_name == "NextGen_Mock" & t == "B") next # pooled results are only needed at baseline for NextGen_Mock
        
        if (grp == "non_naive_vac_pla") {
            dat.plot = subset(dat.immuno.subset.plot3, Bserostatus==1)
            grp_lb = if (study_name == "NextGen_Mock") {"pooled group participants"
                } else {paste0(gsub("-", "", tolower(bstatus.labels.3[2])), " pooled group participants")}
            assays_sub_ = assays
        } else if (grp == "naive_vac"){
            dat.plot = subset(dat.immuno.subset.plot3, Bserostatus==0 & as.character(Trt)==trt.labels[2])
            grp_lb = if (study_name == "NextGen_Mock") {paste0(tolower(trt.labels[2]), " group participants")
            } else {paste0(bstatus.labels.3[1], " ", tolower(trt.labels[2]), " group participants")}
            assays_sub_ = assays
        } else if (study_name == "NextGen_Mock" & t %in% c("B", "Day31", "Day181") &
                   grp %in% c("non_naive_vac", "non_naive_pla")) {
            
            trt_val <- ifelse(grp == "non_naive_vac", trt.labels[2], trt.labels[1])
            #label_val <- bstatus.labels.3[2]
            
            dat.plot <- subset(dat.immuno.subset.plot3, as.character(Trt) == trt_val)
            grp_lb <- paste0(tolower(trt_val), " group participants")
            assays_sub_ <- assays
        } #else if (study_name == "NextGen_Mock" & t %in% c("Day91", "Day366") &
           #        grp %in% c("non_naive_vac", "non_naive_pla")) {
            
            #trt_val <- ifelse(grp == "non_naive_vac", trt.labels[2], trt.labels[1])
            #label_val <- bstatus.labels.3[2]
            
            #dat.plot <- subset(dat.immuno.subset.plot3, as.character(Trt) == trt_val & Track == "A")
            #grp_lb <- paste0(tolower(trt_val), " group participants")
            #assays_sub_ <- assays
        #}
        
        if (attr(config,"config")=="vat08_combined") {
            dat.plot = dat.plot %>% filter(ph2.immuno.bAb==1)
        }
        
        for (asy in if (study_name != "NextGen_Mock") {""} else {
            c("IgG_sera", "bindSpike_IgA_sera", "pseudoneutid50_sera", 
              "IgG_nasal", "bindSpike_IgA_nasal", "pseudoneutid50_nasal",
              "IgG_saliva", "bindSpike_IgA_saliva", "pseudoneutid50_saliva",
              "T4", "T8", 
              
              "IgG_sera_KP.2|bindSpike_IgA_sera_KP.2", "IgG_sera_KP.2|pseudoneutid50_sera", "IgG_sera_KP.2|T4", "IgG_sera_KP.2|T8",
              "bindSpike_IgA_sera_KP.2|pseudoneutid50_sera", "bindSpike_IgA_sera_KP.2|T4", "bindSpike_IgA_sera_KP.2|T8",
              "pseudoneutid50_sera|T4", "pseudoneutid50_sera|T8", 
              
              "IgG_nasal_KP.2|bindSpike_IgA_nasal_KP.2", "IgG_nasal_KP.2|pseudoneutid50_nasal", "IgG_nasal_KP.2|T4", "IgG_nasal_KP.2|T8",
              "bindSpike_IgA_nasal_KP.2|pseudoneutid50_nasal", "bindSpike_IgA_nasal_KP.2|T4", "bindSpike_IgA_nasal_KP.2|T8",
              "pseudoneutid50_nasal|T4", "pseudoneutid50_nasal|T8", 
              
              "IgG_saliva_KP.2|bindSpike_IgA_saliva_KP.2", "IgG_saliva_KP.2|pseudoneutid50_saliva", "IgG_saliva_KP.2|T4", "IgG_saliva_KP.2|T8",
              "bindSpike_IgA_saliva_KP.2|pseudoneutid50_saliva", "bindSpike_IgA_saliva_KP.2|T4", "bindSpike_IgA_saliva_KP.2|T8",
              "pseudoneutid50_saliva|T4", "pseudoneutid50_saliva|T8", 
              
              "T4|T8")}) { 
            
            parts <- strsplit(asy, "\\|")[[1]]
            if (!  (any(grepl(parts[1], assays)) & any(grepl(parts[length(parts)], assays)))  ) next
            
            if (asy == "") {assays_sub = assays_sub_
            } else {assays_sub = assays_sub_[grepl(asy, assays_sub_)]}
        
            covid_corr_pairplots(
                plot_dat = dat.plot,
                time = t,
                assays = assays_sub,
                strata = "all_one",
                weight = ifelse(attr(config,"config")=="vat08_combined", "wt.immuno.bAb", 
                                ifelse(study_name == "NextGen_Mock" & t %in% c("B", "Day31", "Day181"), "wt.immuno", 
                                       ifelse(study_name == "NextGen_Mock" & t %in% c("Day91", "Day366"), "wt.AB.immuno", 
                                              "wt.subcohort"))),
                plot_title = paste0(
                    "Correlations of ", length(assays_sub), " ", ifelse(t=="B", "D01", t), " antibody markers in ", grp_lb, "\nCorr = Weighted Spearman Rank Correlation."
                ),
                column_labels = assay_metadata$assay_label_short[match(assays_sub, assay_metadata$assay)],
                height = max(1.3 * length(assays_sub) + 0.1, 5.5),
                width = max(1.3 * length(assays_sub), 5.5),
                column_label_size = ifelse(max(nchar(paste(assay_metadata$assay_label_short))) > 28 & length(assays_sub) >= 4, 4, ifelse(length(assays_sub) == 2, 6.5, 6)),
                label_format = ifelse((grepl("T4", asy) | grepl("T8", asy)) & !grepl("pseudo|bind", asy), "percent", "log10"),
                filename = paste0(
                    save.results.to, "/pairs_by_time_", t,
                    "_", ifelse(asy == "", length(assays_sub), gsub("\\|", "_and_", asy)), "_markers_", grp, 
                    ifelse(study_name == "NextGen_Mock" & t %in% c("B", "Day31", "Day181"), "_final", 
                           ifelse(study_name == "NextGen_Mock" & t %in% c("Day91", "Day366"), "_initial", "")), 
                    ".pdf"
                )
            )
        }
    }
}

###### Set 4 plots: Correlation plots for a given marker across time points
# all markers, by naive/non-naive, vaccine/placebo, pooling cases and non-cases
for (a in assays){
    
    if (attr(config,"config") %in% c("janssen_partA_VL", "janssen_pooled_partA")) next # janssen_partA_VL doesn't need these plots
    
    panels_set <- list()

    
    if (study_name == "NextGen_Mock") {
        times_list <- list(
            c("B", "Day31", "Day181"),
            c("B", "Day31", "Day91", "Day181", "Day366")
        )
    } else {
        times_list <- list(c("B", paste0("Day", timepoints)))
    }
    

    for (times_sub in times_list) {
    
        i <- 1
        for (trt in trt.labels){
            for (bsero in c(0, 1)){
            
                dat.plot4_ = dat.immuno.subset.plot3 %>% filter(Trt == trt & Bserostatus == bsero)
                if (nrow(dat.plot4_) == 0) next
                
                if (attr(config,"config")=="vat08_combined") {
                    dat.plot4 = dat.plot4_ %>% filter(ph2.immuno.bAb==1)
                }
                
                if (study_name == "NextGen_Mock") {
                    if (length(times_sub) == 3) {
                        dat.plot4 = dat.plot4_ %>% mutate(wt = wt.immuno) # RIS    
                    } else if (length(times_sub) == 5) {
                        dat.plot4 = dat.plot4_ %>% filter(Track == "A") %>% mutate(wt = wt.AB.immuno) # Track A
                    }
                }
                panels_set[[i]] = covid_corr_pairplots(
                    plot_dat = dat.plot4,
                    time = times_sub,
                    assays = a,
                    strata = "all_one",
                    weight = ifelse(attr(config,"config")=="vat08_combined", "wt.immuno.bAb", 
                                    ifelse(study_name == "NextGen_Mock", "wt", "wt.subcohort")),
                    plot_title = "",
                    column_labels = gsub("B", "Day01", times_sub),
                    axis_label_size = ifelse(study_name == "NextGen_Mock", 7, 9),
                    height = 5.5,
                    width = 5.5,
                    column_label_size = 10,
                    label_format = ifelse((grepl("T4", a) | grepl("T8", a)) & !grepl("pseudo|bind", a), "percent", "log10"),
                    write_to_file = F
                ) 
                i = i + 1
            }
        }

    
        y.grob.1 <- textGrob(paste0("(A) ", trt.labels[2]), gp=gpar(fontface="bold", col="black", fontsize=9)) # Vaccine
        y.grob.2 <- textGrob(paste0("(A) ", trt.labels[1]), gp=gpar(fontface="bold", col="black", fontsize=9)) # Placebo
        y.grob.3 <- textGrob(bstatus.labels.3[1], gp=gpar(fontface="bold", col="black", fontsize=9)) # Naive
        y.grob.4 <- textGrob(bstatus.labels.3[2], gp=gpar(fontface="bold", col="black", fontsize=9)) # Non-naive
        
        #add to plot
        if (study_name == "NextGen_Mock") {
            row_plot <- arrangeGrob(
                arrangeGrob(ggmatrix_gtable(panels_set[[1]]), top = y.grob.1),
                arrangeGrob(ggmatrix_gtable(panels_set[[2]]), top = y.grob.2),
                nrow = 1,
                top = textGrob(labels.assays.short[a], gp = gpar(fontsize = 10, fontface = "bold"))
            )
            # Add bottom space with an empty grob
            bottom_space <- textGrob(" ", gp = gpar(fontsize = 10)) 
            combined_p <- grid.arrange(row_plot, bottom_space, ncol = 1, heights = c(1, 0.7))
            
        } else {combined_p <- grid.arrange(
            grid.arrange(arrangeGrob(plot_grid(
                arrangeGrob(ggmatrix_gtable(panels_set[[1]]), top = y.grob.3), arrangeGrob(ggmatrix_gtable(panels_set[[2]]), top = y.grob.4)), left = y.grob.1), nrow=1),
            grid.arrange(arrangeGrob(plot_grid(
                ggmatrix_gtable(panels_set[[3]]), ggmatrix_gtable(panels_set[[4]])), left = y.grob.2), nrow=1),
            #bottom = x.grob,
            ncol = 1,
            top = textGrob(labels.assays.short[a], gp = gpar(fontsize = 10, fontface = "bold"))
        )}
        
        ggsave(filename = paste0(
            save.results.to, "/pairs_across_timepoints_", a, 
            if (study_name == "NextGen_Mock" & length(times_sub) == 5) "_initial", 
            if (study_name == "NextGen_Mock" & length(times_sub) == 3) "_final", ".pdf"), plot = combined_p, width = 8, height = 8, units="in")
    }
}

