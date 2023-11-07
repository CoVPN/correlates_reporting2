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
library(ggnewscale) # for new_scale_color() 

# There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
#if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))
source(here::here("code", "process_violin_pair_functions.R"))
source(here::here("..", "_common.R"))

# for the order of figure panels
assay_order = assay_metadata %>% dplyr::arrange(panel, order_in_panel) %>% select(assay_label_short) %>% pull()
assay_metadata = assay_metadata %>%
    mutate(assay_label_short = factor(assay_label_short,
                                      levels = assay_order
    ))

dat.longer.immuno.subset.plot1 <- readRDS(here::here("data_clean", "longer_immuno_data_plot1.rds")) %>%
    mutate(Trt = factor(Trt, levels = c(0, 1), labels = c("Placebo", "Vaccine")),
           nnaive = factor(Bserostatus,
                            levels = c(0, 1),
                            labels = c("Naive", "Non-naive")))
dat.immuno.subset.plot3 <- readRDS(here::here("data_clean", "twophase_data.rds")); dat.immuno.subset.plot3$all_one <- 1 # as a placeholder for strata values
dat.immuno.subset.plot3 <- dat.immuno.subset.plot3 %>%
    mutate(Trt = factor(Trt, levels = c(0, 1), labels = c("Placebo", "Vaccine")),
           nnaive = factor(Bserostatus,
                            levels = c(0, 1),
                            labels = c("Naive", "Non-naive")))

###### Set 1 plots: Ab distributions at main timepoints and delta by vaccine/placebo, naive/nnaive
set1_times <- times[1:5]
# ID50, at B, D22, D43, D22-B, D43-B
# bindSpike, at B, D22, D43, D22-B, D43-B
for (panel in c("pseudoneutid50", "bindSpike")){
    # by naive/non-naive, vaccine/placebo
    f_1 <- f_by_time_assay(
        dat = dat.longer.immuno.subset.plot1 %>% mutate(x="1"),
        assays = assays[grepl(panel, assays)],
        times = set1_times,
        ylim = c(ifelse(panel=="pseudoneutid50", -5, -6), ifelse(panel=="pseudoneutid50", 5, 11)),
        axis.x.text.size = 20,
        strip.x.text.size = ifelse(panel=="pseudoneutid50", 25, 10),
        panel.text.size = ifelse(panel=="pseudoneutid50", 7, 4.5),
        facet.y.var = vars(Trt_nnaive),
        facet.x.var = vars(assay_label_short))
    
    for (i in 1:length(set1_times)){
        
        file_name <- paste0("/", panel, "_at_", set1_times[i], ".pdf")
        ggsave(plot = f_1[[i]], filename = paste0(save.results.to, file_name), width = 16, height = 16)
    }
}


###### Set 2 plots: Longitudinal plots D1, D22, D43
set2_assays = assays
# bindSpike
# ID50
# bindSpike mdw
for (panel in c("pseudoneutid50", "bindSpike")){
    # by naive/non-naive, vaccine/placebo
    f_2 <- f_longitude_by_assay(
        dat = dat.longer.immuno.subset.plot1,
        x.var = "time",
        x.lb = c("D1","D22","D43"),
        assays = set2_assays[grepl(panel, set2_assays) & !grepl("mdw", set2_assays)],
        ylim = c(0, 4.5),
        times = c("B","Day22","Day43"),
        strip.text.x.size = ifelse(panel=="pseudoneutid50", 25, 12),
        panel.text.size = ifelse(panel=="pseudoneutid50", 6, 4),
        facet.y.var = vars(Trt_nnaive), 
        facet.x.var = vars(assay_label_short)
    )
    
    file_name <- paste0("/", panel, "_longitudinal.pdf")
    ggsave(plot = f_2, filename = paste0(save.results.to, file_name), width = 16, height = 11)
}

f_2 <- f_longitude_by_assay(
    dat = dat.longer.immuno.subset.plot1,
    x.var = "time",
    x.lb = c("D1","D22","D43"),
    assays = "bindSpike_mdw",
    ylim = c(-6, 11),
    times = c("B","Day22","Day43"),
    strip.text.x.size = 25,
    strip.text.y.size = 22,
    panel.text.size = 6,
    facet.y.var = vars(Trt_nnaive), 
    facet.x.var = vars(assay_label_short)
)

file_name <- paste0("/crossreactivity_longitudinal.pdf")
ggsave(plot = f_2, filename = paste0(save.results.to, file_name), width = 6, height = 11)


###### Set 3 plots: Correlation plots across markers at a given time point
# 15 markers, non-naive pooled vaccine & placebo, three timepoints
# 15 markers, naive vaccine, three timepoints
for (grp in c("non_naive_vac_pla", "naive_vac")){
    for (t in c("B","Day22","Day43")) {
        
        if (grp == "non_naive_vac_pla") {dat.plot = subset(dat.immuno.subset.plot3, Bserostatus==1)
        } else {dat.plot = subset(dat.immuno.subset.plot3, Bserostatus==0 & Trt=="Vaccine")}
        
        if (grp == "naive_vac" && t=="B") next
        
        covid_corr_pairplots(
            plot_dat = dat.plot,
            time = t,
            assays = assays,
            strata = "all_one",
            weight = "wt.subcohort",
            plot_title = paste0(
                "Correlations of 15 ", ifelse(t=="B", "D1", t), " antibody markers in ", grp, " grp, Corr = Weighted Spearman Rank Correlation."
            ),
            column_labels = paste(t, assay_metadata$assay_label_short),
            height = max(1.3 * length(assays) + 0.1, 5.5),
            width = max(1.3 * length(assays), 5.5),
            column_label_size = ifelse(max(nchar(paste(t, assay_metadata$assay_label_short)))>40, 3.3, 3.8),
            filename = paste0(
                save.results.to, "/pairs_by_time_", t,
                "_15_markers_", grp, ".pdf"
            )
        )
        
    }
}

###### Set 4 plots: Correlation plots for a given marker across time points
# 15 markers, by naive/non-naive, vaccine/placebo
for (a in assays){
    panels_set <- list()
    i <- 1
    dat.cor.subset.plot3$Trt_nnaive = with(dat.cor.subset.plot3, 
                                           case_when(Trt == 1 & nnaive == 0 ~ "Vaccine Naive", 
                                                     Trt == 1 & nnaive == 1 ~ "Vaccine Non-naive", 
                                                     Trt == 0 & nnaive == 0 ~ "Placebo Naive", 
                                                     Trt == 0 & nnaive == 1 ~ "Placebo Non-naive"))
    for (tn in c("Vaccine Naive", "Vaccine Non-naive", "Placebo Naive", "Placebo Non-naive")){
        for (ce in c("Non-Cases", "Omicron Cases")){
            times_sub = c("BD1","BD29",if(ce=="Omicron Cases") "DD1")
            
            panels_set[[i]] = covid_corr_pairplots(
                plot_dat = dat.cor.subset.plot3 %>% filter(Trt_nnaive == tn & cohort_event == ce),
                time = times_sub,
                assays = a,
                strata = "all_one",
                weight = "wt.BD29",
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
    y.grob.2 <- textGrob("Vaccine   \nNon-naive", gp=gpar(fontface="bold", col="black", fontsize=9))
    y.grob.3 <- textGrob("Placebo     \nNaive", gp=gpar(fontface="bold", col="black", fontsize=9))
    y.grob.4 <- textGrob("Placebo   \nNon-naive", gp=gpar(fontface="bold", col="black", fontsize=9))
    
    #x.grob <- textGrob(paste0(
    #    "Correlations of ", a, " Levels at BD1, BD29 (and DD1) by Cases/Non-cases, Corr = Weighted Spearman Rank Correlation."
    #), gp=gpar(fontface="bold", col="black", fontsize=9))
    
    #add to plot
    combined_p <- grid.arrange(
        grid.arrange(arrangeGrob(plot_grid(
            ggmatrix_gtable(panels_set[[1]]), ggmatrix_gtable(panels_set[[2]])), left = y.grob.1), nrow=1),
        grid.arrange(arrangeGrob(plot_grid(
            ggmatrix_gtable(panels_set[[3]]), ggmatrix_gtable(panels_set[[4]])), left = y.grob.2), nrow=1),
        grid.arrange(arrangeGrob(plot_grid(
            ggmatrix_gtable(panels_set[[5]]), ggmatrix_gtable(panels_set[[6]])), left = y.grob.3), nrow=1),
        grid.arrange(arrangeGrob(plot_grid(
            ggmatrix_gtable(panels_set[[7]]), ggmatrix_gtable(panels_set[[8]])), left = y.grob.4), nrow=1),
        #bottom = x.grob,
        ncol = 1
    )
    
    ggsave(filename = paste0(
        save.results.to, "/pairs_by_timepoints_", a, ".pdf"), plot = combined_p, width = 8, height = 10, units="in")
}


# 9 markers, by naive, non-naive
for (a in assays){
    panels_set <- list()
    i <- 1
    dat.cor.subset.plot3$nnaive_lb = with(dat.cor.subset.plot3, 
                                          case_when(nnaive == 0 ~ "Naive", 
                                                    nnaive == 1 ~ "Non-naive"))
    for (tn in c("Naive", "Non-naive")){
        for (ce in c("Non-Cases", "Omicron Cases")){
            times_sub = c("BD1","BD29",if(ce=="Omicron Cases") "DD1")
            
            panels_set[[i]] = covid_corr_pairplots(
                plot_dat = dat.cor.subset.plot3 %>% filter(nnaive_lb == tn & cohort_event == ce),
                time = times_sub,
                assays = a,
                strata = "all_one",
                weight = "wt.BD29",
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
    
    y.grob.1 <- textGrob("Naive", gp=gpar(fontface="bold", col="black", fontsize=9))
    y.grob.2 <- textGrob("Non-\nnaive", gp=gpar(fontface="bold", col="black", fontsize=9))
    
    #x.grob <- textGrob(paste0(
    #    "Correlations of ", a, " Levels at BD1, BD29 (and DD1) by Cases/Non-cases, Corr = Weighted Spearman Rank Correlation."
    #), gp=gpar(fontface="bold", col="black", fontsize=9))
    
    #add to plot
    combined_p <- grid.arrange(
        grid.arrange(arrangeGrob(plot_grid(
            ggmatrix_gtable(panels_set[[1]]), ggmatrix_gtable(panels_set[[2]])), left = y.grob.1), nrow=1),
        grid.arrange(arrangeGrob(plot_grid(
            ggmatrix_gtable(panels_set[[3]]), ggmatrix_gtable(panels_set[[4]])), left = y.grob.2), nrow=1),
        #bottom = x.grob,
        ncol = 1
    )
    
    ggsave(filename = paste0(
        save.results.to, "/pairs_by_timepoints_", a, "_pooled.pdf"), plot = combined_p, width = 8, height = 10, units="in")
}


## adhoc figure 5
assays_adhoc <- c("bindSpike","pseudoneutid50","bindSpike_BA.1","pseudoneutid50_BA.1")
for (i in 1:length(assays_adhoc)){
    x.var <- paste0("BD1", assays_adhoc)[i]
    y.var <- paste0("DeltaBD29overBD1", assays_adhoc)[i]
    x.lb <- paste("BD1", assay_labels_short[match(assays_adhoc, names(assay_labels_short))])[i]
    y.lb <- paste("DeltaBD29overBD1\n", assay_labels_short[match(assays_adhoc, names(assay_labels_short))])[i]
    
    dat_plot <- dat.long.cor.subset.plot5 %>%
        mutate(Trt_nnaive2 = factor(paste(nnaive, cohort_event), 
                                    levels = c("Naive Omicron Cases", "Naive Non-Cases", "Non-naive Omicron Cases", "Non-naive Non-Cases"),
                                    labels = c("Naive\nOmicron Cases", "Naive\nNon-Cases", "Non-naive\nOmicron Cases", "Non-naive\nNon-Cases")),
               Trt = factor(Trt, levels = c("Vaccine", "Placebo")))
    
    plot_theme <- theme_bw() +
        theme(plot.title = element_text(hjust = 0.5),
              axis.title = element_text(size = 24, face="bold"),
              axis.text = element_text(size = 20),
              strip.text.x = element_text(size = 25), # facet label size
              strip.text.y = element_text(size = 13),
              strip.background = element_rect(fill=NA,colour=NA),
              strip.placement = "outside",
              legend.position = "bottom", 
              legend.text = element_text(size = 16, face="plain"),
              legend.key = element_blank(), # remove square outside legend key
              plot.caption = element_text(size = 26, hjust=0, face="plain"), 
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              plot.margin = margin(5.5, 12, 5.5, 5.5, "pt")) 
    
    p2 <- 
        ggplot(dat_plot %>% filter(assay == assays_adhoc[i]), 
               aes_string(x.var, y.var, color = "Trt")) +
        facet_grid(cols = vars(Trt_nnaive2)) +
        geom_point(size = 2) +
        geom_smooth(method = "loess", se=FALSE, color="red") +
        geom_rug(alpha = 0.6, position = "jitter") +
        geom_hline(aes(yintercept = 0), linetype = "dashed", color = "gray", na.rm = TRUE) +
        scale_color_manual(name = "", values = c("goldenrod2","#378252"), drop=FALSE) +
        scale_x_continuous(
            #limits = c(1.5, 5),
            labels = scales::label_math(10^.x)
        ) +
        scale_y_continuous(
            #limits = c(-1, 3),
            labels = scales::label_math(10^.x)
        ) +
        xlab(x.lb) + 
        ylab(y.lb) + 
        theme_bw() +
        coord_fixed(ratio = 1) +
        plot_theme
    
    file_name <- paste0(assays_adhoc[i], "_scatter_BD1_DeltaBD29overBD1_adhoc2.pdf")
    ggsave(plot = p2, filename = paste0(save.results.to, file_name), width = 16, height = 11)
    
}
