#-----------------------------------------------
# obligatory to append to the top of each script
# renv::activate(project = here::here(".."))
library(GGally) # ggpairs
library(stringr)
library(ggplot2)
library(cowplot) # plot_grid
library(grid) # textGrob
library(gridExtra)
library(wCorr) # weighted correlation
library(ggnewscale) # for new_scale_color() 
library(spatstat.geom) # for ewcdf
library(here)
#install.packages("reldist", repos="http://cran.us.r-project.org")
#library(reldist) # for wtd.quantile
#install.packages("lemon", repos="http://cran.us.r-project.org")
#library(lemon)

# There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
#if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))

source(here::here("..", "_common.R"))
source(here::here("code", "cor_graphics_functions.R"))
source(here::here("code", "cor_process_function.R"))
source(here::here("code", "params.R"))
#colnames(assay_metadata) = gsub("X[.]+","", colnames(assay_metadata))
# add panel for bindN in assay_metadata
assay_metadata[which(assay_metadata$assay=="bindN"), "panel"] = "bindN"


# for the order of figure panels
assay_order = assay_metadata %>% dplyr::arrange(panel, order_in_panel) %>% select(assay_label_short) %>% pull()
assay_metadata = assay_metadata %>%
    mutate(assay_label_short = factor(assay_label_short,
                                      levels = assay_order
))

dat.longer.cor.subset.plot1 <- readRDS(here("data_clean", "longer_cor_data_plot1.rds")) # at level of trt and assay
dat.cor.subset.plot3 <- readRDS(here("data_clean", "cor_data.rds"));dat.cor.subset.plot3$all_one <- 1 # as a placeholder for strata values

stage <- 1
day <- ifelse(study_name=="VAT08" & stage==1, 180, 150)

cases_lb <- if (study_name=="VAT08"){ 
     
    if (nrow(subset(dat.longer.cor.subset.plot1, cohort_event=="C1"))!=0) {c("C1", "C2", "C3")
    } else {c("C2")}
    
    } else if (attr(config,"config") %in% c("prevent19_stage2","azd1222_stage2")) {paste0(config.cor$txt.endpoint, " Cases")
    } else {"Post-Peak Cases"}

cases_lb2 <- if (study_name=="VAT08"){ 
    
    if (nrow(subset(dat.longer.cor.subset.plot1, cohort_event=="C1"))!=0) {c("C1"="C1: 7-27 days PD2 cases", "C2"=sprintf("C2: 28-%s days PD2 cases", day), "C3"=sprintf("C3: 7-%s days PD2 cases", day))
    } else {c("C2"=sprintf("C2: 28-%s days PD2 cases", day))}
    
    } else {cases_lb}

# path for figures and tables etc
save.results.to = here::here("output")
if (!dir.exists(save.results.to))  dir.create(save.results.to)
save.results.to = paste0(here::here("output"), "/", attr(config,"config"),"/");
if (!dir.exists(save.results.to))  dir.create(save.results.to)
if (study_name=="VAT08") {
    save.results.to = paste0(save.results.to, "/", COR, "_stage", unique(dat.longer.cor.subset.plot1$Trialstage), "/");
    if (!dir.exists(save.results.to))  dir.create(save.results.to)
} else {
    save.results.to = paste0(save.results.to, "/", COR,"/");
    if (!dir.exists(save.results.to))  dir.create(save.results.to)
}
print(paste0("save.results.to equals ", save.results.to))

###### Set 1 plots: Ab distributions for assays of one panel, at set1_times, by case/non-case (by naive/non-naive, vaccine/placebo)
set1_times <- labels.time[!grepl(paste0("over D", tinterm), labels.time)] # "Day 1" "Day 22" "Day 43" "D22 fold-rise over D1"  "D43 fold-rise over D1"

assay_metadata <- assay_metadata %>% filter(assay %in% c("bindSpike", "pseudoneutid50"))

# adhoc for vat08_combined: for all timepoints and delta, one plot for naive, the other for non-naive

dat_long_ <- dat.longer.cor.subset.plot1 %>%
    filter(assay %in% c("bindSpike", "pseudoneutid50")) %>% 
    filter(Trt=="Vaccine", Bserostatus=="Non-naive") %>% 
    mutate(subgrp = case_when(assay=="bindSpike" & value < LLoQ ~ "\n< LLOQ",
                              assay=="bindSpike" & value >= LLoQ ~ "\n\u2265 LLOQ",
                              assay=="pseudoneutid50" & value < LLoD ~ "\n< LOD",
                              assay=="pseudoneutid50" & value >= LLoD ~ "\n\u2265 LOD"),
           subgrp_time=gsub(" fold-rise over D1", "", time),
           subgrp_time=gsub("ay ", "", subgrp_time)) %>%
    filter(Bserostatus == "Non-naive") %>% 
    ungroup() 

dat_long <- left_join(dat_long_, 
                      dat_long_ %>% 
                         filter(time %in% c("Day 22", "Day 43")) %>% 
                         mutate(subgrp_time=gsub("ay ", "", time)) %>% 
                         rename(subgrp_=subgrp) %>% 
                         distinct(Ptid, assay, subgrp_time, subgrp_)) %>% 
    mutate(subgrp=ifelse(grepl("fold-rise", time), subgrp_, subgrp)) %>% 
    mutate(Trt_nnaive = factor(paste(Trt, Bserostatus, subgrp), 
                               levels =
                                   paste(rep(c("Vaccine Naive", "Vaccine Non-naive", "Placebo Naive", "Placebo Non-naive"), each=4), 
                                         rep(c("\n< LLOQ", "\n≥ LLOQ","\n< LOD", "\n≥ LOD"), 4))),
           subgrp_time = ifelse(grepl("fold-rise", time), gsub(" fold-rise over D1", "", time), NA)) %>% 
    group_by(Trt, Bserostatus, subgrp, assay, time, cohort_event) %>% 
    mutate(rr=sum(response, na.rm = T)/sum(!is.na(response)), 
           N_RespRate=ifelse(grepl("fold-rise", time), "", sprintf("%s\n%.1f%%", sum(!is.na(response)), rr*100)))


if (attr(config,"config") == "vat08_combined"){

    for (panel in c("pseudoneutid50", "bindSpike")){
    
        for (na in c("Non-naive")){
            
            if (sum(grepl(substr(panel, 1, 4), assay_metadata$assay))==0) next
            
            assay_num = length(assays[grepl(substr(panel,1,4), assays)])
            
            for (tm_subset in set1_times){
                
                # by naive/non-naive, vaccine/placebo
                f_1 <- f_case_non_case_by_time_assay_wrap(
                    dat = dat_long,
                    
                    facet.x.var = "assay_label_short",
                    facet.y.var = "Trt_nnaive",
                    assays = assays[grepl(substr(panel, 1, 4), assays)],
                    times = tm_subset,
                    ylim = if (grepl("Day", tm_subset) & panel=="bindSpike") {c(2, 7)} else if (grepl("Day", tm_subset) & panel=="pseudoneutid50") {c(1, 6.5)} else if (grepl("fold", tm_subset)) {c(-3, 4.2)}, 
                    ybreaks = if (grepl("Day", tm_subset) & panel=="bindSpike") {c(2, 3, 4, 5, 6)} else if (grepl("Day", tm_subset) & panel=="pseudoneutid50") {c(1, 2, 3, 4, 5, 6)} else if (grepl("fold", tm_subset)) {c(-3, -2, -1, 0, 1, 2, 3, 4)},
                    axis.x.text.size = 30,
                    strip.x.text.size = 30,
                    panel.text.size = 14,
                    scale.x.discrete.lb = c(cases_lb, "Non-Cases"),
                    lgdbreaks = c(cases_lb, "Non-Cases", "Non-Responders"),
                    lgdlabels = if (study_name=="VAT08") {c(cases_lb2, "Non-Cases"="Non-Cases", "Non-Responders"="Non-Responders")} else {c(cases_lb, "Non-Cases", "Non-Responders")},
                    chtcols = setNames(c(if(length(cases_lb)==3) "#1749FF", "#FF6F1B", if(length(cases_lb)==3) "#D92321", "#0AB7C9", "#8F8F8F"), c(cases_lb, "Non-Cases", "Non-Responders")), # BLUE, ORANGE, RED, LIGHT BLUE, GRAY
                    chtpchs = setNames(c(if(length(cases_lb)==3) 19, 19, if(length(cases_lb)==3) 19, 19, 2), c(cases_lb, "Non-Cases", "Non-Responders")))
                
                for (i in 1:length(tm_subset)){
                    
                    file_name <- paste0(panel, "_by_case_non_case_at_", tm_subset, ifelse(na=="Naive", "_n", "_nn"), ".pdf")
                    # ggsave(plot = f_1[[i]], filename = paste0(save.results.to, file_name), width = 30, height = 28)
                    cairo_pdf(filename = paste0(save.results.to, file_name), width = 30, height = 28)
                        print(f_1[[i]])
                    dev.off()
                }
            }
        }
    }
}

# adhoc for prevent19_stage2: show y-axis variable in percentile
