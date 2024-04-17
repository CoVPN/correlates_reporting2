#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
library(GGally)
library(stringr)
require(devtools)
install_version("dummies", version = "1.5.6", repos = "http://cran.us.r-project.org")
library(cowplot) # for function plot_grid
library(grid)
library(gridExtra)
install.packages("wCorr", repos = "http://cran.us.r-project.org") # weighted correlation
library(wCorr)
install.packages("ggnewscale", repos = "http://cran.us.r-project.org")
library(ggnewscale) # for new_scale_color() 

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

cases_lb <- if (attr(config,"config")=="vat08_combined"){c("7-27 days PD2 cases", "28-180 days PD2 cases")
    } else {paste0(config.cor$txt.endpoint, " Cases")}

# path for figures and tables etc
save.results.to = here::here("output")
if (!dir.exists(save.results.to))  dir.create(save.results.to)
save.results.to = paste0(here::here("output"), "/", attr(config,"config"),"/");
if (!dir.exists(save.results.to))  dir.create(save.results.to)
save.results.to = paste0(save.results.to, "/", COR,"/");
if (!dir.exists(save.results.to))  dir.create(save.results.to)
print(paste0("save.results.to equals ", save.results.to))

###### Set 1 plots: Ab distributions for assays of one panel, at set1_times, by case/non-case (by naive/non-naive, vaccine/placebo)
set1_times <- labels.time[!grepl(paste0("over D", tinterm), labels.time)] # "Day 1" "Day 22" "Day 43" "D22 fold-rise over D1"  "D43 fold-rise over D1"
if(attr(config,"config") == "prevent19_stage2"){set1_times <- set1_times[set1_times!="Disease Day 1"]}

for (panel in c("pseudoneutid50", "bindSpike", if(attr(config,"config") == "prevent19_stage2") "bindSpike_sub_nvx_stage2")){
    # by naive/non-naive, vaccine/placebo
    f_1 <- f_case_non_case_by_time_assay(
        dat = dat.longer.cor.subset.plot1,
        assays = if(panel=="bindSpike_sub_nvx_stage2") {c("bindSpike_D614","bindSpike_Delta1")
            } else {assays[grepl(panel, assays)]},
        times = set1_times,
        ylim = c(0,4.5), 
        ybreaks = c(0,1,2,3,4),
        axis.x.text.size = ifelse(length(assays[grepl(panel, assays)]) > 7, 13, 25),
        strip.x.text.size = ifelse(length(assays[grepl(panel, assays)]) > 7, 13, 25),
        panel.text.size = ifelse(length(assays[grepl(panel, assays)]) > 7, 4.5, 8),
        scale.x.discrete.lb = c(cases_lb, "Non-Cases"),
        lgdbreaks = c(cases_lb, "Non-Cases", "Non-Responders"),
        chtcols = setNames(c(if(length(cases_lb)==2) "#1749FF","#D92321","#0AB7C9", "#8F8F8F"), c(cases_lb, "Non-Cases", "Non-Responders")),
        chtpchs = setNames(c(if(length(cases_lb)==2) 19, 19, 19, 2), c(cases_lb, "Non-Cases", "Non-Responders")))
    
    for (i in 1:length(set1_times)){
        
        file_name <- paste0(panel, "_by_case_non_case_at_", set1_times[i], ".pdf")
        ggsave(plot = f_1[[i]], filename = paste0(save.results.to, file_name), width = 30, height = 16)
    }

}

###### Set 2 plots: Longitudinal violin plots, by cases and non-cases (by naive/non-naive, vaccine/placebo)
set2.1_assays = assays[!assays %in% c("bindSpike_mdw")]
if(attr(config,"config") == "prevent19_stage2"){set2.1_assays <- set2.1_assays[grepl("Delta$|Delta1$|D614", set2.1_assays)]}
time_cohort.lb = c(paste0(labels.time[!grepl("over", labels.time)], "\n", "Non-Cases"), paste0(labels.time[!grepl("over", labels.time)], "\n", cases_lb[length(cases_lb)]))
if(attr(config,"config") == "prevent19_stage2"){time_cohort.lb <- time_cohort.lb[time_cohort.lb!="Disease Day 1\nNon-Cases"]}

# two assays per plot
for (i in 1:length(set2.1_assays)) {
    
    if (i%%2==0) next     # skip even i
    
    f_2 <- f_longitude_by_assay(
        dat = dat.longer.cor.subset.plot1 %>%
            filter(paste0(time, "\n", cohort_event) %in% time_cohort.lb) %>%
            mutate(time_cohort = factor(paste0(time, "\n", cohort_event), 
                                        levels = time_cohort.lb,
                                        labels = time_cohort.lb)),
        x.lb = time_cohort.lb,
        assays = set2.1_assays[c(i,i+1)],
        panel.text.size = 6,
        ylim = c(0,4.5), 
        ybreaks = c(0,2,4),
        lgdbreaks = c(cases_lb, "Non-Cases", "Non-Responders"),
        chtcols = setNames(c(if(length(cases_lb)==2) "#1749FF","#D92321","#0AB7C9", "#8F8F8F"), c(cases_lb, "Non-Cases", "Non-Responders")),
        chtpchs = setNames(c(if(length(cases_lb)==2) 19, 19, 19, 2), c(cases_lb, "Non-Cases", "Non-Responders")))
    
    file_name <- paste0(paste0(set2.1_assays[c(i,i+1)], collapse="_"), "_longitudinal_by_case_non_case.pdf")
    ggsave(plot = f_2[[1]], filename = paste0(save.results.to, file_name), width = 16, height = 11)
}

if (attr(config,"config") == "vat08_combined"){
    # one assay per plot
    set2.2_assays = c("bindSpike_mdw")
    time_cohort.lb = c(paste(labels.time[!grepl("over", labels.time)], "Non-Cases"), paste(labels.time[!grepl("over", labels.time)], cases_lb[length(cases_lb)]))
    
    for (a in set2.2_assays) {
        
        f_2 <- f_longitude_by_assay(
            dat = dat.longer.cor.subset.plot1 %>%
                filter(time %in% labels.time[!grepl("over", labels.time)]) %>%
                mutate(time_cohort = factor(paste0(time, cohort_event), 
                                            levels = time_cohort.lb,
                                            labels = time_cohort.lb)),
            x.lb = time_cohort.lb,
            assays = a,
            panel.text.size = 6,
            lgdbreaks = c(cases_lb, "Non-Cases", "Non-Responders"),
            chtcols = setNames(c(c("#1749FF","#D92321")[length(cases_lb)],"#0AB7C9", "#8F8F8F"), c(cases_lb, "Non-Cases", "Non-Responders")),
            chtpchs = setNames(c(c(19, 19)[length(cases_lb)], 19, 2), c(cases_lb, "Non-Cases", "Non-Responders")))
        
        file_name <- paste0(a, "_longitudinal_by_case_non_case.pdf")
        ggsave(plot = f_2[[1]], filename = paste0(save.results.to, file_name), width = 8, height = 11)
    }
}

###### Set 3 plots: Correlation plots across markers at a given time point
set3_times = if (attr(config,"config") == "vat08_combined") {times_[!grepl("Delta", times_)] # B, Day22, Day43
} else if (attr(config,"config") == "prevent19_stage2") {times_[!grepl("DD1", times_)]} # BD1, Day35, Delta35overBD1

# add delta values for prevent19_stage2
if (attr(config,"config") == "prevent19_stage2"){
    for (a in assays) {
        dat.cor.subset.plot3[, paste0("Delta35overBD1",a)] = dat.cor.subset.plot3[, paste0("Day35", a)] - dat.cor.subset.plot3[, paste0("BD1", a)]
    }
}

for (grp in c("non_naive_vac_pla", "naive_vac")){
    for (t in set3_times) {
        
        if (grp == "naive_vac" && t=="B") next # this is not needed for VAT08
        
        if (grp == "non_naive_vac_pla") {
            dat.plot = subset(dat.cor.subset.plot3, Bserostatus==1)
            grp_lb = paste0(gsub("-","",bstatus.labels.2[2]), " participants")
        } else if (grp == "naive_vac"){
            dat.plot = subset(dat.cor.subset.plot3, Bserostatus==0 & Trt==1)
            grp_lb = paste0(bstatus.labels.2[1], " vaccine group participants")
        }
        
        if (nrow(dat.plot)==0) next
        
        covid_corr_pairplots(
            plot_dat = dat.plot,
            time = t,
            assays = assays,
            strata = "all_one",
            weight = "wt",#ifelse(grepl(tpeak, t), paste0("wt.D", tpeak), paste0("wt.D", tinterm)),
            plot_title = paste0(
                "Correlations of ", length(assays), " ", t, " antibody markers in ", grp_lb, ", Corr = Weighted Spearman Rank Correlation."
            ),
            column_labels = paste(t, assay_metadata$assay_label_short),
            height = max(1.3 * length(assays) + 0.1, 5.5),
            width = max(1.3 * length(assays), 5.5),
            column_label_size = ifelse(max(nchar(paste(t, assay_metadata$assay_label_short)))>40, 3.3, 3.8),
            filename = paste0(
                save.results.to, "/pairs_by_time_", t,
                "_", length(assays), "_markers_", grp, ".pdf"
            )
        )
    }
}

###### Set 4 plots: Correlation plots for a given marker across time points
# all markers, by naive/non-naive, vaccine/placebo, pooling cases and non-cases
for (a in assays){
    panels_set <- list()
    i <- 1
    
    for (trt in c(1, 0)){
        for (bsero in c(0, 1)){
            times_sub = if (attr(config,"config") == "prevent19_stage2") {c("BD1","Day35","DD1")
                } else {c("B",paste0("Day", tinterm), paste0("Day", tpeak))}
            
            if(nrow(dat.cor.subset.plot3 %>% filter(Trt == trt & Bserostatus == bsero))==0) next
            
            panels_set[[i]] = covid_corr_pairplots(
                plot_dat = dat.cor.subset.plot3 %>% filter(Trt == trt & Bserostatus == bsero),
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
    if (i == 4){
        combined_p <- grid.arrange(
            grid.arrange(arrangeGrob(plot_grid(
                arrangeGrob(ggmatrix_gtable(panels_set[[1]]), top = y.grob.3), arrangeGrob(ggmatrix_gtable(panels_set[[2]]), top = y.grob.4)), left = y.grob.1), nrow=1),
            grid.arrange(arrangeGrob(plot_grid(
                ggmatrix_gtable(panels_set[[3]]), ggmatrix_gtable(panels_set[[4]])), left = y.grob.2), nrow=1),
            #bottom = x.grob,
            ncol = 1
        )
    } else {
        combined_p = arrangeGrob(ggmatrix_gtable(panels_set[[1]]), top = y.grob.3, left = y.grob.1)
    }
    ggsave(filename = paste0(
        save.results.to, "/pairs_across_timepoints_", a, ".pdf"), plot = combined_p, width = 8, height = 10, units="in")
}
