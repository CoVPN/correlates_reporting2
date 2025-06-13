#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
library(GGally)
library(stringr)
require(devtools)
#install_version("dummies", version = "1.5.6", repos = "http://cran.us.r-project.org")
library("dummies")
library(cowplot) # for function plot_grid
library(grid)
library(gridExtra)
#install.packages("wCorr", repos = "http://cran.us.r-project.org") # weighted correlation
library(wCorr)
library(ggnewscale) # for new_scale_color() 

# There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
#if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))

source(here::here("code", "cor_graphics_functions_covail.R"))
source(here::here("..", "_common.R"))

assay_metadata <- assay_metadata %>% filter(panel=="id50")
assays <- assay_metadata$assay

# for the order of figure panels
assay_order = assay_metadata %>%
    mutate(order_in_panel=c(1,3,2,4,5,6)) %>%
    dplyr::arrange(panel, order_in_panel) %>%
    select(assay_label_short) %>%
    pull()

assay_metadata = assay_metadata %>%
    mutate(assay_label_short = factor(assay_label_short,levels = assay_order))


dat.longer.cor.subset.plot1.1 <- readRDS(here::here("data_clean", "longer_cor_data_plot1.1.rds"))
dat.longer.cor.subset.plot1.1.2 <- readRDS(here::here("data_clean", "longer_cor_data_plot1.1.2.rds"))
dat.longer.cor.subset.plot1.2 <- readRDS(here::here("data_clean", "longer_cor_data_plot1.2.rds"))
dat.longer.cor.subset.plot1.3 <- readRDS(here::here("data_clean", "longer_cor_data_plot1.3.rds"))
dat.longer.cor.subset.plot1.4 <- readRDS(here::here("data_clean", "longer_cor_data_plot1.4.rds"))
dat.longer.cor.subset.plot1.4.1 <- readRDS(here::here("data_clean", "longer_cor_data_plot1.4.1.rds"))
dat.longer.cor.subset.plot1.4.2 <- readRDS(here::here("data_clean", "longer_cor_data_plot1.4.2.rds"))
dat.longer.cor.subset.plot1.5 <- readRDS(here::here("data_clean", "longer_cor_data_plot1.5.rds"))
dat.cor.subset.plot3 <- readRDS(here::here("data_clean", "cor_data_plot3.rds")); dat.cor.subset.plot3$all_one <- 1 # as a placeholder for strata values

dat.longer.cor.subset.plot2.1 <- readRDS(here::here("data_clean", "longer_cor_data_plot2.1.rds"))
dat.longer.cor.subset.plot2.2 <- readRDS(here::here("data_clean", "longer_cor_data_plot2.2.rds"))

dat.longer.cor.subset.plot1.1.3 <- readRDS(here::here("data_clean", "longer_cor_data_plot1.1.3.rds"))

# path for figures and tables etc
save.results.to = here::here("output")
if (!dir.exists(save.results.to))  dir.create(save.results.to)
save.results.to = paste0(here::here("output"), "/", attr(config,"config"),"/");
if (!dir.exists(save.results.to))  dir.create(save.results.to)
#save.results.to = paste0(save.results.to, "/", COR,"/");
#if (!dir.exists(save.results.to))  dir.create(save.results.to)
print(paste0("save.results.to equals ", save.results.to))

###### Set 1 plots: BD1 and BD29 Ab distributions by case/non-case
set1_times <- c("B","Day15","Delta15overB")
set1_times_label <- c("Baseline","Day15","Delta15overB")

# pseudoneutid50 D614G, Beta, Delta, BA.1, BA.4.BA.5, MDW at BD1, BD29, BD29-BD1
for (panel in c("id50")){
    
    # 13 one-dose mRNA booster arms pooled
    f_1.1 <- f_case_non_case_by_time_assay(
        dat = dat.longer.cor.subset.plot1.1,
        assays = paste0("pseudoneutid50",c("_D614G", "_Delta", "_Beta", "_BA.1", "_BA.4.BA.5","_MDW")),
        times = set1_times,
        axis.x.text.size = 20,
        strip.x.text.size = 20,
        panel.text.size = 5,
        facet.y.var = NULL, 
        facet.x.var = vars(assay_label_short),
        pointby = "cohort_col",
        lgdbreaks = c("Non-Cases","D22-D91 Cases","D92-D181 Cases","D22-D181 Cases", "Non-Responders"),
        chtcols = setNames(c("#FF6F1B", "#0AB7C9","#0AB7C9", "#0AB7C9", "#8F8F8F"), c("Non-Cases","D22-D91 Cases","D92-D181 Cases","D22-D181 Cases", "Non-Responders")),
        chtpchs = setNames(c(19, 19,19, 19, 2), c("Non-Cases","D22-D91 Cases","D92-D181 Cases","D22-D181 Cases", "Non-Responders")),
        pool_label="(13 one-dose mRNA booster arms pooled)",
        label_angle = -60,
        label_hjust = 0.5)
    
    for (i in 1:length(set1_times)){
        
        file_name <- paste0(panel, "_6_strain_by_case_non_case_at_", set1_times[i]," (13 one-dose mRNA booster arms pooled)", ".pdf")
        ggsave(plot = f_1.1[[i]], filename = paste0(save.results.to, file_name), width = 32, height = 16)
    }
    
    # 13 one-dose mRNA booster arms pooled, naive vs. non-naive
    f_1.1.2 <- f_case_non_case_by_time_assay(
        dat = dat.longer.cor.subset.plot1.1.2,
        assays = paste0("pseudoneutid50",c("_D614G", "_Delta", "_Beta", "_BA.1", "_BA.4.BA.5","_MDW")),
        times = set1_times,
        axis.x.text.size = 20,
        strip.x.text.size = 20,
        strip.y.text.size = 14,
        panel.text.size = 5,
        facet.y.var = vars(assay_label_short), 
        facet.x.var = vars(naive),
        pointby = "cohort_col",
        lgdbreaks = c("Non-Cases","D22-D91 Cases","D92-D181 Cases","D22-D181 Cases", "Non-Responders"),
        chtcols = setNames(c("#FF6F1B", "#0AB7C9","#0AB7C9", "#0AB7C9", "#8F8F8F"), c("Non-Cases","D22-D91 Cases","D92-D181 Cases","D22-D181 Cases", "Non-Responders")),
        chtpchs = setNames(c(19, 19,19, 19, 2), c("Non-Cases","D22-D91 Cases","D92-D181 Cases","D22-D181 Cases", "Non-Responders")),
        pool_label="\n(13 one-dose mRNA booster arms pooled: naive vs. non-naive)",
        label_angle = -10,
        label_hjust = 0.5,
        ncol=5,
        nrow=1)
    
    for (i in 1:length(set1_times)){
        
        file_name <- paste0(panel, "_6_strain_by_case_non_case_at_", set1_times[i]," (13 one-dose mRNA booster arms pooled: naive vs. non-naive)", ".pdf")
        ggsave(plot = f_1.1.2[[i]], filename = paste0(save.results.to, file_name), width = 16, height = 25)
    }
    
    # 3 Sanofi booster arms pooled
    f_1.2 <- f_case_non_case_by_time_assay(
        dat = dat.longer.cor.subset.plot1.2,
        assays = paste0("pseudoneutid50",c("_D614G", "_Delta", "_Beta", "_BA.1", "_BA.4.BA.5","_MDW")),
        times = set1_times,
        axis.x.text.size = 20,
        strip.x.text.size = 20,
        panel.text.size = 5,
        facet.y.var = NULL, 
        facet.x.var = vars(assay_label_short),
        pointby = "cohort_col",
        lgdbreaks = c("Non-Cases","D22-D91 Cases","D92-D181 Cases","D22-D181 Cases", "Non-Responders"),
        chtcols = setNames(c("#FF6F1B", "#0AB7C9","#0AB7C9", "#0AB7C9", "#8F8F8F"), c("Non-Cases","D22-D91 Cases","D92-D181 Cases","D22-D181 Cases", "Non-Responders")),
        chtpchs = setNames(c(19, 19,19, 19, 2), c("Non-Cases","D22-D91 Cases","D92-D181 Cases","D22-D181 Cases", "Non-Responders")),
        pool_label="(3 Sanofi booster arms pooled)",
        label_angle = -60,
        label_hjust = 0.5)
    
    for (i in 1:length(set1_times)){
        
        file_name <- paste0(panel, "_6_strain_by_case_non_case_at_", set1_times[i]," (3 Sanofi booster arms pooled)", ".pdf")
        ggsave(plot = f_1.2[[i]], filename = paste0(save.results.to, file_name), width = 32, height = 16)
    }
    
    # mRNA Moderna vs. mRNA Pfizer
    f_1.3 <- f_case_non_case_by_time_assay(
        dat = dat.longer.cor.subset.plot1.3,
        assays = paste0("pseudoneutid50",c("_D614G", "_Delta", "_Beta", "_BA.1", "_BA.4.BA.5","_MDW")),
        times = set1_times,
        axis.x.text.size = 20,
        strip.x.text.size = 20,
        panel.text.size = 5,
        facet.y.var = vars(TrtA), 
        facet.x.var = vars(assay_label_short),
        pointby = "cohort_col",
        lgdbreaks = c("Non-Cases","D22-D91 Cases","D92-D181 Cases","D22-D181 Cases", "Non-Responders"),
        chtcols = setNames(c("#FF6F1B", "#0AB7C9","#0AB7C9", "#0AB7C9", "#8F8F8F"), c("Non-Cases","D22-D91 Cases","D92-D181 Cases","D22-D181 Cases", "Non-Responders")),
        chtpchs = setNames(c(19, 19,19, 19, 2), c("Non-Cases","D22-D91 Cases","D92-D181 Cases","D22-D181 Cases", "Non-Responders")),
        pool_label="(mRNA Moderna vs. mRNA Pfizer)",
        label_angle = -60,
        label_hjust = 0.5
    )
    
    for (i in 1:length(set1_times)){
        
        file_name <- paste0(panel, "_6_strain_by_case_non_case_at_", set1_times[i], " (mRNA Moderna vs. mRNA Pfizer)", ".pdf")
        ggsave(plot = f_1.3[[i]], filename = paste0(save.results.to, file_name), width = 32, height = 16)
    }
    
    # mRNA Prototype vs. mRNA Omicron-Containing
    f_1.4 <- f_case_non_case_by_time_assay(
        dat = dat.longer.cor.subset.plot1.4,
        assays = paste0("pseudoneutid50",c("_D614G", "_Delta", "_Beta", "_BA.1", "_BA.4.BA.5","_MDW")),
        times = set1_times,
        axis.x.text.size = 20,
        strip.x.text.size = 20,
        panel.text.size = 5,
        facet.y.var = vars(TrtB), 
        facet.x.var = vars(assay_label_short),
        pointby = "cohort_col",
        lgdbreaks = c("Non-Cases","D22-D91 Cases","D92-D181 Cases","D22-D181 Cases", "Non-Responders"),
        chtcols = setNames(c("#FF6F1B", "#0AB7C9","#0AB7C9", "#0AB7C9", "#8F8F8F"), c("Non-Cases","D22-D91 Cases","D92-D181 Cases","D22-D181 Cases", "Non-Responders")),
        chtpchs = setNames(c(19, 19,19, 19, 2), c("Non-Cases","D22-D91 Cases","D92-D181 Cases","D22-D181 Cases", "Non-Responders")),
        pool_label="(mRNA Prototype vs. mRNA Omicron-Containing)",
        label_angle = -60,
        label_hjust = 0.5
    )
    
    for (i in 1:length(set1_times)){
        
        file_name <- paste0(panel, "_6_strain_by_case_non_case_at_", set1_times[i]," (mRNA Prototype vs. mRNA Omicron-Containing)", ".pdf")
        ggsave(plot = f_1.4[[i]], filename = paste0(save.results.to, file_name), width = 32, height = 16)
    }
    
    # Moderna: mRNA Prototype vs. mRNA Omicron-Containing
    f_1.4.1 <- f_case_non_case_by_time_assay(
        dat = dat.longer.cor.subset.plot1.4.1,
        assays = paste0("pseudoneutid50",c("_D614G", "_Delta", "_Beta", "_BA.1", "_BA.4.BA.5","_MDW")),
        times = set1_times,
        axis.x.text.size = 20,
        strip.x.text.size = 20,
        panel.text.size = 5,
        facet.y.var = vars(TrtB_moderna), 
        facet.x.var = vars(assay_label_short),
        pointby = "cohort_col",
        lgdbreaks = c("Non-Cases","D22-D91 Cases","D92-D181 Cases","D22-D181 Cases", "Non-Responders"),
        chtcols = setNames(c("#FF6F1B", "#0AB7C9","#0AB7C9", "#0AB7C9", "#8F8F8F"), c("Non-Cases","D22-D91 Cases","D92-D181 Cases","D22-D181 Cases", "Non-Responders")),
        chtpchs = setNames(c(19, 19,19, 19, 2), c("Non-Cases","D22-D91 Cases","D92-D181 Cases","D22-D181 Cases", "Non-Responders")),
        pool_label="(Moderna: mRNA Prototype vs. mRNA Omicron-Containing)",
        label_angle = -60,
        label_hjust = 0.5
    )
    
    for (i in 1:length(set1_times)){
        
        file_name <- paste0(panel, "_6_strain_by_case_non_case_at_", set1_times[i]," (Moderna: mRNA Prototype vs. mRNA Omicron-Containing)", ".pdf")
        ggsave(plot = f_1.4.1[[i]], filename = paste0(save.results.to, file_name), width = 32, height = 16)
    }
    
    # Pfizer: mRNA Prototype vs. mRNA Omicron-Containing
    f_1.4.2 <- f_case_non_case_by_time_assay(
        dat = dat.longer.cor.subset.plot1.4.2,
        assays = paste0("pseudoneutid50",c("_D614G", "_Delta", "_Beta", "_BA.1", "_BA.4.BA.5","_MDW")),
        times = set1_times,
        axis.x.text.size = 20,
        strip.x.text.size = 20,
        panel.text.size = 5,
        facet.y.var = vars(TrtB_pfizer), 
        facet.x.var = vars(assay_label_short),
        pointby = "cohort_col",
        lgdbreaks = c("Non-Cases","D22-D91 Cases","D92-D181 Cases","D22-D181 Cases", "Non-Responders"),
        chtcols = setNames(c("#FF6F1B", "#0AB7C9","#0AB7C9", "#0AB7C9", "#8F8F8F"), c("Non-Cases","D22-D91 Cases","D92-D181 Cases","D22-D181 Cases", "Non-Responders")),
        chtpchs = setNames(c(19, 19,19, 19, 2), c("Non-Cases","D22-D91 Cases","D92-D181 Cases","D22-D181 Cases", "Non-Responders")),
        pool_label="(Pfizer: mRNA Prototype vs. mRNA Omicron-Containing)",
        label_angle = -60,
        label_hjust = 0.5
    )
    
    for (i in 1:length(set1_times)){
        
        file_name <- paste0(panel, "_6_strain_by_case_non_case_at_", set1_times[i]," (Pfizer: mRNA Prototype vs. mRNA Omicron-Containing)", ".pdf")
        ggsave(plot = f_1.4.2[[i]], filename = paste0(save.results.to, file_name), width = 32, height = 16)
    }
    
    # mRNA Bivalent vs. mRNA Monovalent
    f_1.5 <- f_case_non_case_by_time_assay(
        dat = dat.longer.cor.subset.plot1.5,
        assays = paste0("pseudoneutid50",c("_D614G", "_Delta", "_Beta", "_BA.1", "_BA.4.BA.5","_MDW")),
        times = set1_times,
        axis.x.text.size = 20,
        strip.x.text.size = 20,
        panel.text.size = 5,
        facet.y.var = vars(TrtC), 
        facet.x.var = vars(assay_label_short),
        pointby = "cohort_col",
        lgdbreaks = c("Non-Cases","D22-D91 Cases","D92-D181 Cases","D22-D181 Cases", "Non-Responders"),
        chtcols = setNames(c("#FF6F1B", "#0AB7C9","#0AB7C9", "#0AB7C9", "#8F8F8F"), c("Non-Cases","D22-D91 Cases","D92-D181 Cases","D22-D181 Cases", "Non-Responders")),
        chtpchs = setNames(c(19, 19,19, 19, 2), c("Non-Cases","D22-D91 Cases","D92-D181 Cases","D22-D181 Cases", "Non-Responders")),
        pool_label="(mRNA Bivalent vs. mRNA Monovalent)",
        label_angle = -60,
        label_hjust = 0.5
    )
    
    for (i in 1:length(set1_times)){
        
        file_name <- paste0(panel, "_6_strain_by_case_non_case_at_", set1_times[i]," (mRNA Bivalent vs. mRNA Monovalent)" , ".pdf")
        ggsave(plot = f_1.5[[i]], filename = paste0(save.results.to, file_name), width = 32, height = 16)
    }
    
}




###### Set 3 plots: Correlation plots across markers at a given time point
# 9 markers, three timepoints
for (t in c("B","Day15","Day29","Day91")) {
    covid_corr_pairplots(
        plot_dat = subset(dat.cor.subset.plot3,dat.cor.subset.plot3$TrtonedosemRNA==1),
        time = t,
        assays = assays,
        strata = "all_one",
        weight = "wt.D15",
        plot_title = ifelse(t=="B", paste0(
            "Correlations of nAb ID50 markers at ", sub("B", "Baseline",t), " (13 one-dose mRNA booster arms pooled), Corr = Weighted Spearman Rank"
        ),
        paste0(
            "Correlations of nAb ID50 markers at ", sub("Day([0-9]+)", "Day \\1", t), " (13 one-dose mRNA booster arms pooled), Corr = Weighted Spearman Rank"
        )),
        # plot_title = paste0(
        #     "Correlations of nAb ID50 markers at", t, "(13 one-dose mRNA booster arms pooled), Corr = Weighted Spearman Rank"
        # ),
        column_labels = paste(t, assay_metadata$assay_label_short),
        height = max(1.3 * length(assays) + 0.1, 5.5),
        width = max(1.3 * length(assays), 5.5),
        column_label_size = ifelse(max(nchar(paste(t, assay_metadata$assay_label_short)))>40, 3.3, 3.8),
        filename = paste0(
            save.results.to, "/pairs_by_time_", t,
            "_6_markers (13 one-dose mRNA booster arms pooled).pdf"
        )
    )
    
    covid_corr_pairplots(
        plot_dat = subset(dat.cor.subset.plot3,dat.cor.subset.plot3$stage==3),
        time = t,
        assays = assays,
        strata = "all_one",
        weight = "wt.D15",
        plot_title = ifelse(t=="B", paste0(
            "Correlations of nAb ID50 markers at ", sub("B", "Baseline",t), " (3 Sanofi booster arms pooled), Corr = Weighted Spearman Rank"
        ),
        paste0(
            "Correlations of nAb ID50 markers at ", sub("Day([0-9]+)", "Day \\1", t), " (3 Sanofi booster arms pooled), Corr = Weighted Spearman Rank"
        )),
        # plot_title = paste0(
        #     "Correlations of nAb ID50 markers at ", t, " (3 Sanofi booster arms pooled), Corr = Weighted Spearman Rank"
        # ),
        column_labels = paste(t, assay_metadata$assay_label_short),
        height = max(1.3 * length(assays) + 0.1, 5.5),
        width = max(1.3 * length(assays), 5.5),
        column_label_size = ifelse(max(nchar(paste(t, assay_metadata$assay_label_short)))>40, 3.3, 3.8),
        filename = paste0(
            save.results.to, "/pairs_by_time_", t,
            "_6_markers (3 Sanofi booster arms pooled).pdf"
        )
    )
}


###### Set 4 plots: Correlation plots for a given marker across time points
# 6 markers (13 one-dose mRNA booster arms pooled)
for (a in assays){
    
    panels_set <- list()
    i <- 1
    
    for (tn in c(1)){
        #for (ce in c(1)){
        times_sub = c("B","Day15","Day29","Day91")
        
        panels_set[[i]] = covid_corr_pairplots(
            plot_dat = dat.cor.subset.plot3 %>% filter(TrtonedosemRNA== tn),
            time = times_sub,
            assays = a,
            strata = "all_one",
            weight = "wt.D15",
            plot_title = paste0(
                "Correlations of nAb ID50 against ", sub("pseudoneutid50_", "", a), " across time points (13 one-dose mRNA booster arms pooled),\nCorr = Weighted Spearman Rank"
            ),
            column_labels = times_sub,
            height = 5.5,
            width = 5.5,
            column_label_size = 10,
            write_to_file = F
        ) 
        i = i + 1
        #}
    }
    
    
    #add to plot
    combined_p <- grid.arrange(
        grid.arrange(arrangeGrob(plot_grid(
            ggmatrix_gtable(panels_set[[1]]))))
    )
    
    ggsave(filename = paste0(
        save.results.to, "/pairs_by_timepoints_", a, " (13 one-dose mRNA booster arms pooled).pdf"), plot = combined_p, width = 8, height = 10, units="in")
}



# 6 markers (3 Sanofi booster arms pooled)
for (a in assays){
    
    panels_set <- list()
    i <- 1
    
    for (tn in c(3)){
        
        times_sub = c("B","Day15","Day29","Day91")
        
        panels_set[[i]] = covid_corr_pairplots(
            plot_dat = dat.cor.subset.plot3 %>% filter(stage== tn),
            time = times_sub,
            assays = a,
            strata = "all_one",
            weight = "wt.D15",
            plot_title = paste0(
                "Correlations of nAb ID50 against ", sub("pseudoneutid50_", "", a), " across time points (3 Sanofi booster arms pooled),\nCorr = Weighted Spearman Rank"
            ),
            column_labels = times_sub,
            height = 5.5,
            width = 5.5,
            column_label_size = 10,
            write_to_file = F
        ) 
        i = i + 1
        
    }
    
    
    #add to plot
    combined_p <- grid.arrange(
        grid.arrange(arrangeGrob(plot_grid(
            ggmatrix_gtable(panels_set[[1]]))))
    )
    
    ggsave(filename = paste0(
        save.results.to, "/pairs_by_timepoints_", a, " (3 Sanofi booster arms pooled).pdf"), plot = combined_p, width = 8, height = 10, units="in")
}



####### for manuscript 1 ###################
############################################

# 13 one-dose mRNA booster arms pooled, naive vs. non-naive, BA.1

set1_times <- c("B","Day15","Delta15overB")
set1_times_label <- c("Baseline","Day15","Delta15overB")

# pseudoneutid50 BA.1
for (panel in c("id50")){
    dat.longer.cor.subset.plot1.1.2_ms<- dat.longer.cor.subset.plot1.1.2 %>%
        mutate(naive=factor(naive,levels = c("Naive","non-Naive")))%>%
        filter(assay=="pseudoneutid50_BA.1")
    
    f_ms <- f_case_non_case_by_time_assay(
        dat = dat.longer.cor.subset.plot1.1.2_ms,
        assays = paste0("pseudoneutid50",c("_BA.1")),
        times = set1_times,
        axis.x.text.size = 20,
        strip.x.text.size = 20,
        strip.y.text.size = 14,
        panel.text.size = 5,
        facet.y.var = vars(assay_label_short), 
        facet.x.var = vars(naive),
        pointby = "cohort_col",
        lgdbreaks = c("Non-Cases","D22-D91 Cases","D92-D181 Cases","D22-D181 Cases", "Non-Responders"),
        chtcols = setNames(c("red3", "#0AB7C9","#0AB7C9", "#0AB7C9", "#8F8F8F"), c("Non-Cases","D22-D91 Cases","D92-D181 Cases","D22-D181 Cases", "Non-Responders")),
        chtpchs = setNames(c(19, 19,19, 19, 2), c("Non-Cases","D22-D91 Cases","D92-D181 Cases","D22-D181 Cases", "Non-Responders")),
        pool_label="\n(13 one-dose mRNA booster arms pooled_naive vs. non-naive)",
        label_angle = -10,
        label_hjust = 0.5,
        ncol=5,
        nrow=1)
    
    for (i in 1:length(set1_times)){
        
        file_name <- paste0("ms_",panel, "_BA.1_by_case_non_case_at_", set1_times[i]," (13 one-dose mRNA booster arms pooled_naive vs. non-naive)", ".pdf")
        ggsave(plot = f_ms[[i]], filename = paste0(save.results.to, file_name), width = 12, height = 8)
    }
    
}

# YL: for the purpose of manuscript revision on 12/11/2024, requested by Lindsay:
sum_stat_fig2cd = dat.longer.cor.subset.plot1.1.2_ms %>% filter(assay == "pseudoneutid50_BA.1" & time %in% c("B", "Day15")) %>%
    select(cohort_event, time, assay, naive, RespRate, counts, min, q1, median, q3, max) %>%
    unique()
write.csv(sum_stat_fig2cd, "./output/covail/sum_stat_fig2cd_20241211.csv", row.names = F)



# 13 one-dose mRNA booster arms pooled, naive vs. non-naive, female vs male, BA.1
set1_times <- c("B","Day15","Delta15overB")
set1_times_label <- c("Baseline","Day15","Delta15overB")

# pseudoneutid50 BA.1
for (panel in c("id50")){
    dat.longer.cor.subset.plot1.1.3_ms<- dat.longer.cor.subset.plot1.1.3 %>%
        mutate(naive=factor(naive,levels = c("Naive","non-Naive")))%>%
        filter(assay=="pseudoneutid50_BA.1")
    
    f_ms <- f_case_non_case_by_time_assay(
        dat = dat.longer.cor.subset.plot1.1.3_ms,
        assays = paste0("pseudoneutid50",c("_BA.1")),
        times = set1_times,
        axis.x.text.size = 20,
        strip.x.text.size = 20,
        strip.y.text.size = 14,
        panel.text.size = 5,
        facet.y.var = vars(Sex), 
        facet.x.var = vars(naive),
        pointby = "cohort_col",
        lgdbreaks = c("Non-Cases","D22-D91 Cases","D92-D181 Cases","D22-D181 Cases", "Non-Responders"),
        chtcols = setNames(c("red3", "#0AB7C9","#0AB7C9", "#0AB7C9", "#8F8F8F"), c("Non-Cases","D22-D91 Cases","D92-D181 Cases","D22-D181 Cases", "Non-Responders")),
        chtpchs = setNames(c(19, 19,19, 19, 2), c("Non-Cases","D22-D91 Cases","D92-D181 Cases","D22-D181 Cases", "Non-Responders")),
        pool_label="\n(13 one-dose mRNA booster arms pooled_naive vs. non-naive)",
        label_angle = -10,
        label_hjust = 0.5,
        ncol=5,
        nrow=1)
    
    for (i in 1:length(set1_times)){
        
        file_name <- paste0("ms_",panel, "_BA.1_by_case_non_case_sex_at_", set1_times[i]," (13 one-dose mRNA booster arms pooled_naive vs. non-naive)", ".pdf")
        ggsave(plot = f_ms[[i]], filename = paste0(save.results.to, file_name), width = 12, height = 8)
    }
    
}



# update Supplementary Figure 10 y-axis range from 10^-1 to 10^-4
for (panel in c("id50")){
    # 13 one-dose mRNA booster arms pooled, naive vs. non-naive
    f_1.1.2 <- f_case_non_case_by_time_assay(
        dat = dat.longer.cor.subset.plot1.1.2,
        assays = paste0("pseudoneutid50",c("_D614G", "_Delta", "_Beta", "_BA.1", "_BA.4.BA.5","_MDW")),
        times = set1_times,
        axis.x.text.size = 20,
        strip.x.text.size = 20,
        strip.y.text.size = 14,
        panel.text.size = 5,
        facet.y.var = vars(assay_label_short), 
        facet.x.var = vars(naive),
        pointby = "cohort_col",
        lgdbreaks = c("Non-Cases","D22-D91 Cases","D92-D181 Cases","D22-D181 Cases", "Non-Responders"),
        chtcols = setNames(c("red3", "#0AB7C9","#0AB7C9", "#0AB7C9", "#8F8F8F"), c("Non-Cases","D22-D91 Cases","D92-D181 Cases","D22-D181 Cases", "Non-Responders")),
        chtpchs = setNames(c(19, 19,19, 19, 2), c("Non-Cases","D22-D91 Cases","D92-D181 Cases","D22-D181 Cases", "Non-Responders")),
        pool_label="\n(13 one-dose mRNA booster arms pooled: naive vs. non-naive)",
        label_angle = -10,
        label_hjust = 0.5,
        ncol=5,
        nrow=1)
    
    for (i in 1:length(set1_times)){
        
        file_name <- paste0("ms_",panel, "_6_strain_by_case_non_case_at_", set1_times[i]," (13 one-dose mRNA booster arms pooled_naive vs. non-naive)", ".pdf")
        ggsave(plot = f_1.1.2[[i]], filename = paste0(save.results.to, file_name), width = 16, height = 25)
    }
}






######## Manuscript #2 - Sanofi brief communication to Clin Inf Dis ###################
#######################################################################################

#### supp figures 1, Violin box plots of the distribution of fore of infection (FOI) score of 
#vaccine recipients by COVID-19 outcome status by individual vaccine and pooled over vaccine arms. 
#Do not need boost-proximal/boost-distal cases just need all cases vs. non-cases.


plot_theme <- theme_bw(base_size = 25) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size = 18,angle = -30,hjust = 0.5),
          axis.text.y = element_text(size = 25, face="plain"),
          axis.title = element_text(size = 24, face="bold"),
          strip.text.x = element_text(size = 18), # facet label size
          strip.text.y = element_text(size = 25),
          strip.background = element_rect(fill=NA,colour=NA),
          strip.placement = "outside",
          legend.position = "bottom", 
          legend.text = element_text(size = 16, face="plain"),
          legend.key = element_blank(), # remove square outside legend key
          plot.caption = element_text(size = 26, hjust=0, face="plain"), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          plot.margin = margin(5.5, 12, 5.5, 5.5, "pt")) 

#by individual vaccine arms
p2.1.1 <- ggplot(data=dat.longer.cor.subset.plot2.1, aes(x = cohort_event, y = FOIoriginal)) +
    facet_grid(col = vars(arm)) +
    geom_violin(aes(color = cohort_event), scale = "width", na.rm = TRUE, show.legend = FALSE) +
    geom_boxplot(aes(color = cohort_event), width = 0.25, lwd = 1.5, alpha = 0.3, stat = "boxplot", outlier.shape = NA, show.legend = FALSE) +
    scale_color_manual(name = "",
                       values = setNames(c("#FF6F1B", "#0AB7C9"), c("Non-Cases","D22-D181 Cases")),
                       guide = "none") + 
    geom_point(alpha=0.2,
               aes(color = cohort_event),
               position=position_jitter(width = 0.25, height = 0),
               size = 2,
               show.legend = TRUE)+
    scale_shape_manual(name = "",
                       values = setNames(c(19, 19), c("Non-Cases","D22-D181 Cases")),
                       breaks = c("Non-Cases","D22-D181 Cases"),
                       drop=FALSE) +
    scale_x_discrete(labels = c("Non-Cases","D22-D181 Cases"), drop=FALSE) +
    scale_y_continuous() +
    labs(x = "Cohort", 
         y = "FOI score",
         title = "FOI score by non-case/case (individual vaccine arms)") +
    plot_theme +
    guides(color = guide_legend(nrow = 1), shape = guide_legend(ncol = 1))

p2.1.1
ggsave(plot = p2.1.1, filename = paste0(save.results.to, "ms2_supp1_foi_score_ind.pdf"), width = 14, height = 8)

#by pooled over vaccine arms
p2.1.2 <- ggplot(data=dat.longer.cor.subset.plot2.1, aes(x = cohort_event, y = FOIoriginal)) +
    geom_violin(aes(color = cohort_event), scale = "width", na.rm = TRUE, show.legend = FALSE) +
    geom_boxplot(aes(color = cohort_event), width = 0.25, lwd = 1.5, alpha = 0.3, stat = "boxplot", outlier.shape = NA, show.legend = FALSE) +
    scale_color_manual(name = "",
                       values = setNames(c("#FF6F1B", "#0AB7C9"), c("Non-Cases","D22-D181 Cases")),
                       guide = "none") + 
    geom_point(alpha=0.2,
               aes(color = cohort_event),
               position=position_jitter(width = 0.25, height = 0),
               size = 2,
               show.legend = TRUE)+
    scale_shape_manual(name = "",
                       values = setNames(c(19, 19), c("Non-Cases","D22-D181 Cases")),
                       breaks = c("Non-Cases","D22-D181 Cases"),
                       drop=FALSE) +
    scale_x_discrete(labels = c("Non-Cases","D22-D181 Cases"), drop=FALSE) +
    scale_y_continuous() +
    labs(x = "Cohort", 
         y = "FOI score",
         title = "FOI score by non-case/case (pooled vaccine arms)") +
    plot_theme +
    guides(color = guide_legend(nrow = 1), shape = guide_legend(ncol = 1))

p2.1.2
ggsave(plot = p2.1.2, filename = paste0(save.results.to, "ms2_supp1_foi_score_pooled.pdf"), width = 11, height = 8)


#### supp figures 9, 10, 11, 12, 13, 14, Violin box plots of D1, D15, fold-rise (D1 to D15) levels of the six nAb titer markers 
#(D614G, Delta, Beta, BA.1, BA.4/BA.5, weighted average) by
#non-cases, booster-proximal COVID-19 endpoint cases, booster-distal COVID-19 endpoint cases,
#and COVID-19 endpoint cases (proximal + distal) for participants SARS-CoV-2 naïve and non-naïve at enrollment.
#Calculations are the same as in Supplementary Figure 3 (f_1.2). 
# 3 Sanofi booster arms pooled

set1_times <- c("B","Day15","Delta15overB")
set1_times_label <- c("Baseline","Day15","Delta15overB")

# pseudoneutid50 D614G, Beta, Delta, BA.1, BA.4.BA.5, MDW at BD1, BD29, BD29-BD1
panel <- c("id50")

#naïve
f_2.2.1 <- f_case_non_case_by_time_assay(
    dat = dat.longer.cor.subset.plot2.2 %>% filter(naive=="Naive"),
    assays = paste0("pseudoneutid50",c("_D614G", "_Delta", "_Beta", "_BA.1", "_BA.4.BA.5","_MDW")),
    times = set1_times,
    axis.x.text.size = 20,
    strip.x.text.size = 20,
    panel.text.size = 5,
    facet.y.var = NULL, 
    facet.x.var = vars(assay_label_short),
    pointby = "cohort_col",
    lgdbreaks = c("Non-Cases","D22-D91 Cases","D92-D181 Cases","D22-D181 Cases", "Non-Responders"),
    chtcols = setNames(c("#FF6F1B", "#0AB7C9","#0AB7C9", "#0AB7C9", "#8F8F8F"), c("Non-Cases","D22-D91 Cases","D92-D181 Cases","D22-D181 Cases", "Non-Responders")),
    chtpchs = setNames(c(19, 19,19, 19, 2), c("Non-Cases","D22-D91 Cases","D92-D181 Cases","D22-D181 Cases", "Non-Responders")),
    pool_label="(3 Sanofi booster arms pooled_naive)",
    label_angle = -60,
    label_hjust = 0.5)

for (i in 1:length(set1_times)){
    
    file_name <- paste0("ms2_",panel, "_6_strain_by_case_non_case_at_", set1_times[i]," (3 Sanofi booster arms pooled_naive)", ".pdf")
    ggsave(plot = f_2.2.1[[i]], filename = paste0(save.results.to, file_name), width = 32, height = 16)
}

#non-naïve
f_2.2.2 <- f_case_non_case_by_time_assay(
    dat = dat.longer.cor.subset.plot2.2 %>% filter(naive=="non-Naive"),
    assays = paste0("pseudoneutid50",c("_D614G", "_Delta", "_Beta", "_BA.1", "_BA.4.BA.5","_MDW")),
    times = set1_times,
    axis.x.text.size = 20,
    strip.x.text.size = 20,
    panel.text.size = 5,
    facet.y.var = NULL, 
    facet.x.var = vars(assay_label_short),
    pointby = "cohort_col",
    lgdbreaks = c("Non-Cases","D22-D91 Cases","D92-D181 Cases","D22-D181 Cases", "Non-Responders"),
    chtcols = setNames(c("#FF6F1B", "#0AB7C9","#0AB7C9", "#0AB7C9", "#8F8F8F"), c("Non-Cases","D22-D91 Cases","D92-D181 Cases","D22-D181 Cases", "Non-Responders")),
    chtpchs = setNames(c(19, 19,19, 19, 2), c("Non-Cases","D22-D91 Cases","D92-D181 Cases","D22-D181 Cases", "Non-Responders")),
    pool_label="(3 Sanofi booster arms pooled_non-naive)",
    label_angle = -60,
    label_hjust = 0.5)

for (i in 1:length(set1_times)){
    
    file_name <- paste0("ms2_",panel, "_6_strain_by_case_non_case_at_", set1_times[i]," (3 Sanofi booster arms pooled_non_naive)", ".pdf")
    ggsave(plot = f_2.2.2[[i]], filename = paste0(save.results.to, file_name), width = 32, height = 16)
}




#Supplementary Figure 15. Correlations of Wt. Avg. nAb titer at D1 and fold-rise from D1 to D15, D1 to D29, and D1 to D91
#among the three Sanofi booster arms pooled. 
for (a in assays){
    
    panels_set <- list()
    i <- 1
    
    for (tn in c(3)){
        
        times_sub = c("B","Delta15overB","Delta29overB","Delta91overB")
        
        panels_set[[i]] = covid_corr_pairplots(
            plot_dat = dat.cor.subset.plot3 %>% filter(stage== tn),
            time = times_sub,
            assays = a,
            strata = "all_one",
            weight = "wt.D15",
            plot_title = paste0(
                "Correlations of nAb ID50 against ", sub("pseudoneutid50_", "", a), " across time points (3 Sanofi booster arms pooled),\nCorr = Weighted Spearman Rank"
            ),
            column_labels = times_sub,
            height = 5.5,
            width = 5.5,
            column_label_size = 10,
            write_to_file = F
        ) 
        i = i + 1
        
    }
    
    
    #add to plot
    combined_p <- grid.arrange(
        grid.arrange(arrangeGrob(plot_grid(
            ggmatrix_gtable(panels_set[[1]]))))
    )
    
    ggsave(filename = paste0(
        save.results.to, "/ms2_pairs_by_timepoints_fr_", a, " (3 Sanofi booster arms pooled).pdf"), plot = combined_p, width = 8, height = 10, units="in")
}
