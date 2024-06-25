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

# There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
#if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))
source(here::here("code", "process_violin_pair_functions.R"))
Sys.setenv(DESCRIPTIVE = 1)
source(here::here("..", "_common.R"))
source(here::here("code", "params.R")) # load parameters

# for the order of figure panels
if ("order_in_panel" %in% colnames(assay_metadata)){
    assay_order = assay_metadata %>% dplyr::arrange(panel, order_in_panel) %>% select(assay_label_short) %>% pull()
    assay_metadata = assay_metadata %>%
        mutate(assay_label_short = factor(assay_label_short,
                                          levels = assay_order
        ))
}

dat.longer.immuno.subset.plot1 <- readRDS(here::here("data_clean", "longer_immuno_data_plot1.rds"))
if(attr(config,"config")=="janssen_partA_VL") {
    dat.longer.immuno.subset.plot1.2 <- readRDS(here::here("data_clean", "longer_immuno_data_plot1.2.rds"))}
dat.immuno.subset.plot3 <- readRDS(here::here("data_clean", "twophase_data.rds")); dat.immuno.subset.plot3$all_one <- 1 # as a placeholder for strata values
dat.immuno.subset.plot3 <- dat.immuno.subset.plot3 %>%
    mutate(Trt = factor(Trt, levels = c(0, 1), labels = c("Placebo", "Vaccine")),
           nnaive = factor(Bserostatus,
                            levels = c(0, 1),
                            labels = bstatus.labels))

###### Set 1 plots: Ab distributions at main timepoints and delta by vaccine/placebo, naive/nnaive
set1_times <- times[1:5]
# ID50, at B, D22, D43, D22-B, D43-B
# bindSpike, at B, D22, D43, D22-B, D43-B
for (panel in c("pseudoneutid50", "bindSpike")){
    # by naive/non-naive, vaccine/placebo
    
    if (attr(config,"config") %in% c("janssen_partA_VL", "janssen_pooled_partA")) next # janssen_partA_VL doesn't need these plots
    
    f_1 <- f_by_time_assay(
        dat = dat.longer.immuno.subset.plot1 %>% mutate(x="1"),
        assays = assays[grepl(panel, assays)],
        times = set1_times,
        ylim = c(-3, 4.5), #c(ifelse(panel=="pseudoneutid50", -5, -6), ifelse(panel=="pseudoneutid50", 6, 11)),
        axis.x.text.size = 20,
        strip.x.text.size = ifelse(panel=="pseudoneutid50", 25, 10),
        panel.text.size = ifelse(panel=="pseudoneutid50", 7, 4.5),
        facet.y.var = vars(Trt_nnaive),
        facet.x.var = vars(assay_label2))
    
    for (i in 1:length(set1_times)){
        
        file_name <- paste0("/", panel, "_at_", set1_times[i], ".pdf")
        ggsave(plot = f_1[[i]], filename = paste0(save.results.to, file_name), width = 16, height = 16)
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
for (panel in c("pseudoneutid50", "bindSpike")){
    # by naive/non-naive, vaccine/placebo
    
    if (attr(config,"config") %in% c("janssen_partA_VL", "janssen_pooled_partA")) next # janssen_partA_VL doesn't need these plots
    
    f_2 <- f_longitude_by_assay(
        dat = dat.longer.immuno.subset.plot1,
        x.var = "time",
        x.lb = c("D1","D22","D43"),
        assays = set2_assays[grepl(panel, set2_assays) & !grepl("mdw", set2_assays)],
        ylim = c(0, 4.5), #c(0, 5.2),
        times = c("B","Day22","Day43"),
        strip.text.x.size = ifelse(panel=="pseudoneutid50", 25, 12),
        panel.text.size = ifelse(panel=="pseudoneutid50", 6, 4),
        facet.y.var = vars(Trt_nnaive), 
        facet.x.var = vars(assay_label_short)
    )
    
    file_name <- paste0("/", ifelse(panel=="pseudoneutid50", "nAb", ifelse(panel=="bindSpike", "bAb", "")), "_longitudinal.pdf")
    ggsave(plot = f_2, filename = paste0(save.results.to, file_name), width = 16, height = 11)
}

if (attr(config,"config") == "vat08_combined"){
    for (panel in c("pseudoneutid50", "bindSpike")){
        # by naive/non-naive, vaccine/placebo
        
        if (attr(config,"config") %in% c("janssen_partA_VL", "janssen_pooled_partA")) next # janssen_partA_VL doesn't need these plots
        
        f_2 <- f_longitude_by_assay(
            dat = dat.longer.immuno.subset.plot1,
            x.var = "time",
            x.lb = c("D1","D22","D43"),
            assays = set2_assays[grepl(panel, set2_assays) & grepl("mdw", set2_assays)],
            ylim = c(0, 4.5), #c(0, 5.2),
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
for (grp in c("non_naive_vac_pla", "naive_vac")){
    
    if (attr(config,"config") %in% c("janssen_partA_VL", "janssen_pooled_partA")) next # janssen_partA_VL doesn't need these plots
    
    for (t in c("B","Day22","Day43")) {
        
        if (grp == "naive_vac" && t=="B") next # this is not needed for VAT08
        
        if (grp == "non_naive_vac_pla") {
            dat.plot = subset(dat.immuno.subset.plot3, Bserostatus==1)
            grp_lb = paste0(gsub("-","",bstatus.labels.3[2]), " participants")
            assays_sub = assays
        } else if (grp == "naive_vac"){
            dat.plot = subset(dat.immuno.subset.plot3, Bserostatus==0 & as.character(Trt)=="Vaccine")
            grp_lb = paste0(bstatus.labels.3[1], " vaccine group participants")
            assays_sub = assays
        }
        
        covid_corr_pairplots(
            plot_dat = dat.plot,
            time = t,
            assays = assays_sub,
            strata = "all_one",
            weight = "wt.subcohort",
            plot_title = paste0(
                "Correlations of 15 ", ifelse(t=="B", "D1", t), " antibody markers in ", grp_lb, ", Corr = Weighted Spearman Rank Correlation."
            ),
            column_labels = assay_metadata$assay_label_short[match(assays_sub, assay_metadata$assay)],
            height = max(1.3 * length(assays_sub) + 0.1, 5.5),
            width = max(1.3 * length(assays_sub), 5.5),
            column_label_size = ifelse(max(nchar(paste(assay_metadata$assay_label_short)))>28, 4, 6.5),
            filename = paste0(
                save.results.to, "/pairs_by_time_", t,
                "_", length(assays_sub), "_markers_", grp, ".pdf"
            )
        )
        
    }
}

###### Set 4 plots: Correlation plots for a given marker across time points
# all markers, by naive/non-naive, vaccine/placebo, pooling cases and non-cases
for (a in assays){
    
    if (attr(config,"config") %in% c("janssen_partA_VL", "janssen_pooled_partA")) next # janssen_partA_VL doesn't need these plots
    
    panels_set <- list()
    i <- 1
    
    for (trt in c("Vaccine", "Placebo")){
        for (bsero in c(0, 1)){
            times_sub = c("B", paste0("Day", timepoints))
            
            panels_set[[i]] = covid_corr_pairplots(
                plot_dat = dat.immuno.subset.plot3 %>% filter(Trt == trt & Bserostatus == bsero),
                time = times_sub,
                assays = a,
                strata = "all_one",
                weight = "wt.subcohort",
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
