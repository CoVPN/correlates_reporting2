#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
library(GGally) # ggpairs
library(stringr)
library(cowplot) # plot_grid
library(grid) # textGrob
library(gridExtra)
library(wCorr) # weighted correlation
library(ggnewscale) # for new_scale_color() 
library(spatstat.geom) # for ewcdf
#install.packages("reldist", repos="http://cran.us.r-project.org")
#library(reldist) # for wtd.quantile
#install.packages("lemon", repos="http://cran.us.r-project.org")
#library(lemon)

# There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
#if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))

source(here::here("..", "_common.R"))
source(here::here("code", "cor_graphics_functions.R"))
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

cases_lb <- if (study_name=="VAT08"){ 
     
    if (nrow(subset(dat.longer.cor.subset.plot1, cohort_event=="C1"))!=0) {c("C1", "C2", "C3")
    } else {c("C2")}
    
    } else if (attr(config,"config") %in% c("prevent19_stage2","azd1222_stage2")) {paste0(config.cor$txt.endpoint, " Cases")
    } else {"Post-Peak Cases"}

cases_lb2 <- if (study_name=="VAT08"){ 
    
    if (nrow(subset(dat.longer.cor.subset.plot1, cohort_event=="C1"))!=0) {c("C1"="C1: 7-27 days PD2 cases", "C2"="C2: 28-180 days PD2 cases", "C3"="C3: 7-180 days PD2 cases")
    } else {c("C2"="C2: 28-180 days PD2 cases")}
    
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
if (attr(config,"config") == "prevent19_stage2"){
    set1_times <- set1_times[!set1_times %in% c("Booster Day 1", "Disease Day 1")]
} else if (attr(config,"config") == "azd1222_stage2") {set1_times <- set1_times[set1_times!="Day 360"]
} else if (attr(config,"config") == "prevent19nvx") {set1_times <- set1_times[set1_times!="Day 1"]}

for (panel in c("pseudoneutid50", if(attr(config,"config")!="prevent19_stage2") "bindSpike", if(attr(config,"config") %in% c("prevent19_stage2","azd1222_stage2")) "bindSpike_sub_stage2")){
    
    if (sum(grepl(substr(panel, 1, 4), assay_metadata$assay))==0) next
    
    assay_num = length(assays[grepl(substr(panel,1,4), assays)])
    if (panel=="bindSpike_sub_stage2") assay_num = 2
    
    for (tm_subset in c("Day", if(sum(grepl("fold", set1_times))>0) "fold")){
        
        set1_times_sub = set1_times[grepl(tm_subset, set1_times)]
        
        # by naive/non-naive, vaccine/placebo
        f_1 <- f_case_non_case_by_time_assay(
            dat = dat.longer.cor.subset.plot1 %>%
                mutate(Trt_nnaive = factor(paste(Trt, Bserostatus), 
                                           levels = paste(rep(c("Vaccine","Placebo"),each=2), bstatus.labels),
                                           labels = paste0(rep(c("Vaccine","Placebo"),each=2), "\n", bstatus.labels.2))),
            
            facet.y.var = vars(Trt_nnaive),
            assays = if(panel=="bindSpike_sub_stage2") {c("bindSpike_D614","bindSpike_Delta1")
                } else if (attr(config,"config")=="prevent19nvx") {assays} else {assays[grepl(substr(panel, 1, 4), assays)]},
            times = set1_times_sub,
            ylim = if (attr(config,"config") == "nvx_uk302") {c(1, 6)} else if (attr(config,"config") == "prevent19nvx") {c(0,6.6)} else if (study_name == "VAT08" & tm_subset == "Day") {c(0, 4.2)} else if (study_name == "VAT08" & tm_subset == "fold") {c(-3, 4.2)} else {c(0, 5.5)}, 
            ybreaks = if (attr(config,"config") == "nvx_uk302") {c(1,2,3,4,5)} else if (attr(config,"config") == "prevent19nvx") {c(0,1,2,3,4,5,6)} else if (study_name == "VAT08" & tm_subset == "Day") {c(0, 1, 2, 3, 4)} else if (study_name == "VAT08" & tm_subset == "fold") {c(-3, -2, -1, 0, 1, 2, 3, 4)} else {c(0,1,2,3,4,5)},
            axis.x.text.size = ifelse(assay_num > 7 & length(cases_lb)==3, 13, ifelse(assay_num > 5, 20, ifelse(assay_num > 3, 25, 32))),
            strip.x.text.size = ifelse(assay_num > 7, 10, ifelse(assay_num > 5, 18, ifelse(assay_num > 3, 25, 32))),
            panel.text.size = ifelse(assay_num > 7 && length(cases_lb)==3, 3, ifelse(assay_num > 5, 5, ifelse(assay_num > 3, 6, 12))),
            scale.x.discrete.lb = c(cases_lb, "Non-Cases"),
            lgdbreaks = c(cases_lb, "Non-Cases", "Non-Responders"),
            lgdlabels = if (study_name=="VAT08") {c(cases_lb2, "Non-Cases"="Non-Cases", "Non-Responders"="Non-Responders")} else {c(cases_lb, "Non-Cases", "Non-Responders")},
            chtcols = setNames(c(if(length(cases_lb)==3) "#1749FF", "#FF6F1B", if(length(cases_lb)==3) "#D92321", "#0AB7C9", "#8F8F8F"), c(cases_lb, "Non-Cases", "Non-Responders")), # BLUE, ORANGE, RED, LIGHT BLUE, GRAY
            chtpchs = setNames(c(if(length(cases_lb)==3) 19, 19, if(length(cases_lb)==3) 19, 19, 2), c(cases_lb, "Non-Cases", "Non-Responders")))
        
        for (i in 1:length(set1_times_sub)){
            
            file_name <- paste0(panel, "_by_case_non_case_at_", set1_times_sub[i], ".pdf")
            ggsave(plot = f_1[[i]], filename = paste0(save.results.to, file_name), width = 30, height = 16)
        }
    }
}

# adhoc for prevent19_stage2: show y-axis variable in percentile
if(attr(config,"config") == "prevent19_stage2"){
    
    for (panel in c("pseudoneutid50","bindSpike_sub_stage2")){
        
        # by naive/non-naive, vaccine/placebo
        f_1 <- f_case_non_case_by_time_assay_adhoc(
            dat = dat.longer.cor.subset.plot1 %>%
                mutate(Trt_nnaive = factor(paste(Trt, Bserostatus), 
                                           levels = paste(rep(c("Vaccine","Placebo"),each=2), bstatus.labels),
                                           labels = paste0(rep(c("Vaccine","Placebo"),each=2), "\n", bstatus.labels.2))) %>%
                group_by(Trt, Bserostatus, time, assay) %>%
                mutate(lb = "",
                       lb2 = "",
                       lbval = 1/3,
                       lbval2 = 2/3,
                       value = ewcdf(value, wt)(value)
                       ) %>%
                group_by(Trt, Bserostatus, time, assay, cohort_event) %>%
                mutate(lower = weighted_percentile(x=value, percentile=0.25, weights=wt),
                       middle = weighted_percentile(x=value, percentile=0.5, weights=wt),
                       upper = weighted_percentile(x=value, percentile=0.75, weights=wt),
                       IQR = upper - lower, 
                       ymax = min(weighted_percentile(x=value, percentile=1, weights=wt), upper + 1.5*IQR),
                       ymin = max(weighted_percentile(x=value, percentile=0, weights=wt), lower - 1.5*IQR)) %>%
                ungroup() %>%
                mutate(facet.y.var = Trt_nnaive, 
                       facet.x.var = assay),
            
            #facet.y.var = vars(Trt_nnaive),
            assays = if(panel=="bindSpike_sub_stage2") {c("bindSpike_D614","bindSpike_Delta1")
            } else {assays[grepl(substr(panel, 1, 4), assays)]},
            times = set1_times,
            ylim = c(0, 1.2), 
            ybreaks = c(0,0.2,0.4,0.6,0.8,1),
            axis.x.text.size = 32,
            strip.x.text.size = 32,
            panel.text.size = 12,
            scale.x.discrete.lb = c(cases_lb, "Non-Cases"),
            lgdbreaks = c(cases_lb, "Non-Cases", "Non-Responders"),
            chtcols = setNames(c(if(length(cases_lb)==2) "#1749FF","#D92321","#0AB7C9", "#8F8F8F"), c(cases_lb, "Non-Cases", "Non-Responders")),
            chtpchs = setNames(c(if(length(cases_lb)==2) 19, 19, 19, 2), c(cases_lb, "Non-Cases", "Non-Responders")))
        
        for (i in 1:length(set1_times)){
            
            file_name <- paste0(panel, "_by_case_non_case_at_", set1_times[i], "_percentile.pdf")
            ggsave(plot = f_1[[i]], filename = paste0(save.results.to, file_name), width = 30, height = 16)
        }
    }
}

###### Set 2 plots: Longitudinal violin plots, by cases and non-cases (by naive/non-naive, vaccine/placebo)
set2.1_assays = assays[!assays %in% c("bindSpike_mdw")]
if(attr(config,"config") == "prevent19_stage2"){set2.1_assays <- set2.1_assays[grepl("Delta$|Delta1$|D614", set2.1_assays)]}

if(study_name=="VAT08"){
    
    stopifnot(nrow(subset(dat.longer.cor.subset.plot1, time=="Day 43" & cohort_event=="C1"))==0)
    
    if (nrow(subset(dat.longer.cor.subset.plot1, cohort_event=="C1"))!=0){
        time_cohort.lb <- c(paste0(labels.time[1:3], "\n", "Non-Cases"), paste0(labels.time[1:2], "\n", cases_lb[1]), paste0(labels.time[1:3], "\n", cases_lb[2]), paste0(labels.time[1:2], "\n", cases_lb[3]))
    } else {time_cohort.lb <- c(paste0(labels.time[1:3], "\n", "Non-Cases"), paste0(labels.time[1:3], "\n", cases_lb))}

} else if(attr(config,"config") == "prevent19_stage2"){
    
    time_cohort.lb = c(paste0(labels.time, "\n", "Non-Cases"), paste0(labels.time, "\n", cases_lb[length(cases_lb)]))
    time_cohort.lb <- time_cohort.lb[!time_cohort.lb %in% c("Booster Day 1\nDelta COVID Cases", "Booster Day 1\nSevere COVID Cases", "Disease Day 1\nNon-Cases")]

} else {time_cohort.lb = c(paste0(labels.time, "\n", "Non-Cases"), paste0(labels.time, "\n", cases_lb[length(cases_lb)]))
}

# two assays per plot
for (i in 1:length(set2.1_assays)) {
    
    if (attr(config,"config") %in% c("nvx_uk302","prevent19nvx")) next # no need for nvx_uk302, prevent19nvx
    
    if (i%%2==0 & attr(config,"config") != "azd1222_stage2") next     # skip even i for all studies but AZ stage 2
    
    # modify ptid for VAT08 
    if (study_name == "VAT08"){
        dat.longer.cor.subset.plot1_ = dat.longer.cor.subset.plot1 %>%
            mutate(Ptid = paste0(Ptid, cohort_event))
            
    } else {dat.longer.cor.subset.plot1_ = dat.longer.cor.subset.plot1}
    
    f_2 <- f_longitude_by_assay(
        dat = dat.longer.cor.subset.plot1_ %>%
            filter(paste0(time, "\n", cohort_event) %in% time_cohort.lb) %>%
            mutate(time_cohort = factor(paste0(time, "\n", cohort_event), 
                                        levels = time_cohort.lb,
                                        labels = time_cohort.lb),
                   Trt_nnaive = factor(paste(Trt, Bserostatus), 
                                       levels = paste(rep(c("Vaccine","Placebo"),each=2), bstatus.labels),
                                       labels = paste0(rep(c("Vaccine","Placebo"),each=2), "\n", bstatus.labels.2))),
        x.var = "time_cohort",
        x.lb = time_cohort.lb,
        facet.y.var = vars(Trt_nnaive),
        
        assays = if(attr(config,"config") == "azd1222_stage2"){set2.1_assays[i]} else {set2.1_assays[c(i,i+1)]},
        panel.text.size = ifelse(study_name=="VAT08" & length(cases_lb)==3, 2, ifelse(study_name=="VAT08" & length(cases_lb)==1, 4, 5.8)),
        ylim = c(0,4.5), 
        ybreaks = c(0,1,2,3,4),
        axis.text.x.size = ifelse(attr(config,"config") == "prevent19_stage2" | (study_name=="VAT08" & length(cases_lb)==3), 8.4, 9.5),
        lgdbreaks = c(cases_lb, "Non-Cases", "Non-Responders"),
        lgdlabels = if (study_name=="VAT08") {c(cases_lb2, "Non-Cases"="Non-Cases", "Non-Responders"="Non-Responders")} else {c(cases_lb, "Non-Cases", "Non-Responders")},
        chtcols = setNames(c(if(length(cases_lb)==3) "#1749FF", "#FF6F1B", if(length(cases_lb)==3) "#D92321", "#0AB7C9", "#8F8F8F"), c(cases_lb, "Non-Cases", "Non-Responders")), # BLUE, ORANGE, RED, LIGHT BLUE, GRAY
        chtpchs = setNames(c(if(length(cases_lb)==3) 19, 19, if(length(cases_lb)==3) 19, 19, 2), c(cases_lb, "Non-Cases", "Non-Responders")))
    
    file_name <- paste0(paste0(if(attr(config,"config") == "azd1222_stage2"){set2.1_assays[i]} else {set2.1_assays[c(i,i+1)]}, 
                               collapse="_"), 
                        "_longitudinal_by_case_non_case.pdf")
    ggsave(plot = f_2[[1]], filename = paste0(save.results.to, file_name), width = 16, height = 11)
}

if (study_name=="VAT08"){
    # one assay per plot
    set2.2_assays = c("bindSpike_mdw")
    time_cohort.lb = c(paste0(labels.time[1:3], "\n", "Non-Cases"), paste0(labels.time[1:2], "\n", cases_lb[1]), paste0(labels.time[1:3], "\n", cases_lb[2]), paste0(labels.time[1:2], "\n", cases_lb[3]))
    
    dat.longer.cor.subset.plot1_ = dat.longer.cor.subset.plot1 %>%
            mutate(Ptid = paste0(Ptid, cohort_event))

    for (a in set2.2_assays) {
        
        f_2 <- f_longitude_by_assay(
            dat = dat.longer.cor.subset.plot1_ %>%
                filter(paste0(time, "\n", cohort_event) %in% time_cohort.lb) %>%
                mutate(time_cohort = factor(paste0(time, "\n", cohort_event), 
                                            levels = time_cohort.lb,
                                            labels = time_cohort.lb),
                       Trt_nnaive = factor(paste(Trt, Bserostatus), 
                                           levels = paste(rep(c("Vaccine","Placebo"),each=2), bstatus.labels),
                                           labels = paste0(rep(c("Vaccine","Placebo"),each=2), "\n", bstatus.labels.2))),
            x.var = "time_cohort",
            x.lb = time_cohort.lb,
            facet.y.var = vars(Trt_nnaive),
            
            assays = a,
            panel.text.size = ifelse(study_name=="VAT08" & length(cases_lb)==3, 2, ifelse(study_name=="VAT08" & length(cases_lb)==1, 4, 5.8)),
            ylim = c(0,4.5), 
            ybreaks = c(0,1,2,3,4),
            axis.text.x.size = ifelse(attr(config,"config") == "prevent19_stage2" | (study_name=="VAT08" & length(cases_lb)==3), 8.4, 9.5),
            lgdbreaks = c(cases_lb, "Non-Cases", "Non-Responders"),
            lgdlabels = if (study_name=="VAT08") {c(cases_lb2, "Non-Cases"="Non-Cases", "Non-Responders"="Non-Responders")} else {c(cases_lb, "Non-Cases", "Non-Responders")},
            chtcols = setNames(c(if(length(cases_lb)==3) "#1749FF", "#FF6F1B", if(length(cases_lb)==3) "#D92321", "#0AB7C9", "#8F8F8F"), c(cases_lb, "Non-Cases", "Non-Responders")), # BLUE, ORANGE, RED, LIGHT BLUE, GRAY
            chtpchs = setNames(c(if(length(cases_lb)==3) 19, 19, if(length(cases_lb)==3) 19, 19, 2), c(cases_lb, "Non-Cases", "Non-Responders")))
        
        file_name <- paste0(a, "_longitudinal_by_case_non_case.pdf")
        ggsave(plot = f_2[[1]], filename = paste0(save.results.to, file_name), width = 8, height = 11)
    }
}

###### Set 3 plots: Correlation plots across markers at a given time point
set3_times = if (attr(config,"config") == "vat08_combined") {times_[!grepl("Delta", times_)] # B, Day22, Day43
} else if (attr(config,"config") == "prevent19_stage2") {times_[!grepl("DD1|C1|BD1", times_)] # Day35
} else if (attr(config,"config") == "prevent19nvx") {"Day35"
} else {times_}

for (grp in c("non_naive_vac_pla", "naive_vac")){
    for (t in set3_times) {
        
        if (attr(config,"config") %in% c("nvx_uk302")) next # no need for nvx_uk302
        
        if (grp == "naive_vac" && t=="B") next # this is not needed for VAT08
        
        if (grp == "non_naive_vac_pla") {
            dat.plot = subset(dat.cor.subset.plot3, Bserostatus==1)
            grp_lb = paste0(gsub("-","",bstatus.labels.2[2]), " participants")
        } else if (grp == "naive_vac"){
            dat.plot = subset(dat.cor.subset.plot3, Bserostatus==0 & Trt==1)
            grp_lb = paste0(bstatus.labels.2[1], " vaccine group participants")
        }
        
        if (nrow(dat.plot)==0) next
        
        for (assay_list in c(1, if (attr(config,"config") == "prevent19_stage2") 2)){
        
            if (attr(config,"config") == "prevent19_stage2" & assay_list==2) {
                assay_metadata_ = assay_metadata %>% filter(assay %in% c("pseudoneutid50_D614G", "pseudoneutid50_Delta", "bindSpike_D614", "bindSpike_Delta1"))
            } else {assay_metadata_ = assay_metadata}
            
            if (study_name=="VAT08" & t %in% c("B","Day22")){
                dat.plot_ = dat.plot %>% filter(ph2.D22.both == 1) %>% mutate(wt = wt.D22.both)
            } else if (study_name=="VAT08" & t %in% c("Day43")){
                dat.plot_ = dat.plot %>% filter(ph2.D43.both == 1 & EarlyinfectionD43==0) %>% mutate(wt = wt.D43.both)
            }else {dat.plot_ = dat.plot}
            
            covid_corr_pairplots(
                plot_dat = dat.plot_,
                time = t,
                assays = assay_metadata_$assay,
                strata = "all_one",
                weight = "wt",#ifelse(grepl(tpeak, t), paste0("wt.D", tpeak), paste0("wt.D", tinterm)),
                plot_title = paste0(
                    "Correlations of ", length(assay_metadata_$assay), " ", t, " antibody markers in ", grp_lb, ",\nCorr = Weighted Spearman Rank Correlation."
                ),
                column_labels = paste(t, assay_metadata_$assay_label_short),
                height = max(1.3 * length(assay_metadata_$assay) + 0.1, 5.5),
                width = max(1.3 * length(assay_metadata_$assay), 5.5),
                column_label_size = ifelse(max(nchar(paste(t, assay_metadata_$assay_label_short)))>40, 3.4, 
                                           ifelse(length(assay_metadata_$assay)<=3, 6, 3.8)),
                filename = paste0(
                    save.results.to, "/pairs_by_time_", t,
                    "_", length(assay_metadata_$assay), "_markers_", grp, ".pdf"
                )
            )
            
            # adhoc for azd1222_stage2: day 57, a pairwise scatterplot of D57 binding Abs (Reference and Delta) vs. neutralizing antibodies (Reference and Delta)
            #                           we use the bAb weights and we can note in the caption that this is the case
            if (attr(config,"config") == "azd1222_stage2" & COR %in% c("D57azd1222_stage2_delta_bAb") & t=="Day57"){
                
                assay_metadata_adhoc = data.frame(assays = c("bindSpike_D614",
                                                             "bindSpike_Delta1",
                                                             "pseudoneutid50_D614G",
                                                             "pseudoneutid50_Delta"),
                                                  assay_label_short = c("Day57 Anti Spike IgG - Reference (AU/ml)",
                                                                        "Day57 Anti Spike IgG - Delta AY.4.2 (AU/ml)",
                                                                        "Day57 nAb ID50 D614G (AU/ml)",
                                                                        "Day57 nAb DI50 Delta (AU/ml)"))
                
                covid_corr_pairplots(
                    plot_dat = dat.plot,
                    time = t,
                    assays = assay_metadata_adhoc$assays,
                    strata = "all_one",
                    weight = "wt",#ifelse(grepl(tpeak, t), paste0("wt.D", tpeak), paste0("wt.D", tinterm)),
                    plot_title = paste0(
                        "Correlations of 4 ", t, " antibody markers in ", grp_lb, ",\nCorr = Weighted Spearman Rank Correlation. (bAb weights are used)"
                    ),
                    column_labels = assay_metadata_adhoc$assay_label_short,
                    height = max(1.3 * 4 + 0.1, 5.5),
                    width = max(1.3 * 4, 5.5),
                    column_label_size = 3.4,
                    filename = paste0(
                        save.results.to, "/pairs_by_time_", t,
                        "_4_markers_", grp, "_adhoc.pdf"
                    )
                )
                
            }
        }
    }
}

###### Set 4 plots: Correlation plots for a given marker across time points
# all markers, by naive/non-naive, vaccine/placebo, (pooling cases and non-cases)
if (study_name == "VAT08") {
    for (a in assays){
        panels_set <- list()
        i <- 1
        
        for (trt in c(1, 0)){
            for (bsero in c(0, 1)){
                times_sub = c("B",paste0("Day", tinterm), paste0("Day", tpeak))
                
                if(nrow(dat.cor.subset.plot3 %>% filter(Trt == trt & Bserostatus == bsero))==0) next
                
                if (study_name=="VAT08" & t %in% c("B","Day22")){
                    dat.cor.subset.plot3_ = dat.cor.subset.plot3 %>% filter(ph2.D22.both == 1) %>% mutate(wt = wt.D22.both)
                } else if (study_name=="VAT08" & t %in% c("Day43")){
                    dat.cor.subset.plot3_ = dat.cor.subset.plot3 %>% filter(ph2.D43.both == 1 & EarlyinfectionD43==0) %>% mutate(wt = wt.D43.both)
                }else {dat.cor.subset.plot3_ = dat.cor.subset.plot3}
                
                panels_set[[i]] = covid_corr_pairplots(
                    plot_dat = dat.cor.subset.plot3_ %>% filter(Trt == trt & Bserostatus == bsero),
                    time = times_sub,
                    assays = a,
                    strata = "all_one",
                    weight = "wt",#paste0("wt.D", tpeak),
                    plot_title = "",
                    column_labels = times_sub,
                    height = 5.5,
                    width = 5.5,
                    column_label_size = 10,
                    write_to_file = F
                ) 
                i = i + 1
            }
        }
        
        y.grob.1 <- textGrob("Vaccine", gp=gpar(fontface="bold", col="black", fontsize=9))
        y.grob.2 <- textGrob("Placebo", gp=gpar(fontface="bold", col="black", fontsize=9))
        y.grob.3 <- textGrob("Naive", gp=gpar(fontface="bold", col="black", fontsize=9))
        y.grob.4 <- textGrob("Non-naive", gp=gpar(fontface="bold", col="black", fontsize=9))
        
        #add to plot
        combined_p <- grid.arrange(
            grid.arrange(arrangeGrob(plot_grid(
                arrangeGrob(ggmatrix_gtable(panels_set[[1]]), top = y.grob.3), arrangeGrob(ggmatrix_gtable(panels_set[[2]]), top = y.grob.4)), left = y.grob.1), nrow=1),
            grid.arrange(arrangeGrob(plot_grid(
                ggmatrix_gtable(panels_set[[3]]), ggmatrix_gtable(panels_set[[4]])), left = y.grob.2), nrow=1),
            #bottom = x.grob,
            ncol = 1
        )
        ggsave(filename = paste0(
            save.results.to, "/pairs_across_timepoints_", a, ".pdf"), plot = combined_p, width = 8, height = 10, units="in")
    }
}

if (attr(config,"config") %in% c("prevent19_stage2","azd1222_stage2")) {
    for (a in assays){
        panels_set <- list()
        i <- 1
        dat.cor.subset.plot3$Trt_nnaive = "Vaccine Naive"
        
        for (tn in c("Vaccine Naive")){
            for (ce in c("Non-Cases", cases_lb)){
                times_sub = if (attr(config,"config") == "prevent19_stage2") {c("Day35","C1",if(ce=="Non-Cases") "BD1", if(ce==cases_lb) "DD1")
                    } else if (attr(config,"config") == "azd1222_stage2") {c("Day57","Day90","Day180", if(ce=="Non-Cases") "Day360")}
                
                panels_set[[i]] = covid_corr_pairplots(
                    plot_dat = dat.cor.subset.plot3 %>% filter(Trt_nnaive == tn & cohort_event == ce),
                    time = times_sub,
                    assays = a,
                    strata = "all_one",
                    weight = "wt",
                    plot_title = "",
                    column_labels = times_sub,
                    height = 5.5,
                    width = 5.5,
                    column_label_size = 10,
                    write_to_file = F
                ) 
                i = i + 1
            }
        }
        
        y.grob.1 <- textGrob("Vaccine     \nNaive", gp=gpar(fontface="bold", col="black", fontsize=9))
        #y.grob.2 <- textGrob("Vaccine   \nNon-naive", gp=gpar(fontface="bold", col="black", fontsize=9))
        #y.grob.3 <- textGrob("Placebo     \nNaive", gp=gpar(fontface="bold", col="black", fontsize=9))
        #y.grob.4 <- textGrob("Placebo   \nNon-naive", gp=gpar(fontface="bold", col="black", fontsize=9))
        y.grob.5 <- textGrob("Non-cases", gp=gpar(fontface="bold", col="black", fontsize=9))
        y.grob.6 <- textGrob(cases_lb, gp=gpar(fontface="bold", col="black", fontsize=9))
        
        #add to plot
        combined_p <- grid.arrange(
            grid.arrange(arrangeGrob(plot_grid(
                arrangeGrob(ggmatrix_gtable(panels_set[[1]]), top = y.grob.5), arrangeGrob(ggmatrix_gtable(panels_set[[2]]), top = y.grob.6)), left = y.grob.1), nrow=1),
            #grid.arrange(arrangeGrob(plot_grid(
            #    ggmatrix_gtable(panels_set[[3]]), ggmatrix_gtable(panels_set[[4]])), left = y.grob.2), nrow=1),
            #grid.arrange(arrangeGrob(plot_grid(
            #    ggmatrix_gtable(panels_set[[5]]), ggmatrix_gtable(panels_set[[6]])), left = y.grob.3), nrow=1),
            #grid.arrange(arrangeGrob(plot_grid(
            #    ggmatrix_gtable(panels_set[[7]]), ggmatrix_gtable(panels_set[[8]])), left = y.grob.4), nrow=1),
            #bottom = x.grob,
            ncol = 1
        )
        
        ggsave(filename = paste0(
            save.results.to, "/pairs_across_timepoints_", a, ".pdf"), plot = combined_p, width = 8, height = 5, units="in")
    }
}
