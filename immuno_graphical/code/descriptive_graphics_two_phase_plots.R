#Sys.setenv(TRIAL = "vat08m")
#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
# There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))
source(here::here("..", "_common.R"))
#-----------------------------------------------

# install.packages(c("ggpubr", "GGally", "SWIM", "scales", "dummies",
# "gridExtra", "PResiduals"))
library(here)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(GGally)
library(spatstat.geom)
library(scales)
#library(dummies) # this package got archived on 2022-04-29
require(devtools)
install_version("dummies", version = "1.5.6", repos = "http://cran.us.r-project.org")
library(gridExtra)
library(PResiduals)

# produce geom_statistics w/ resampling-based covariate-adjusted Spearman
source(here("code", "params.R"))
source(here("code", "ggally_cor_resample.R"))
source(here("code", "covid_corr_plot_functions.R"))

set.seed(12345)
# load cleaned data
dat.long.twophase.sample <- readRDS(here(
  "data_clean",
  "long_twophase_data.rds"
))

dat.twophase.sample <- readRDS(here("data_clean", "twophase_data.rds"))

assay_lim <- readRDS(here("data_clean", "assay_lim.rds"))

tps_no_delta_over_tinterm <-  times[!times %in% c(paste0("Delta",timepoints[length(timepoints)],"over",timepoints[1]))] #c("B", "Day29", "Delta29overB", "Day57", "Delta57overB")
tps_no_B_and_delta_over_tinterm <-  times[!times %in% c("B",paste0("Delta",timepoints[length(timepoints)],"over",timepoints[1]))] #c("Day29", "Delta29overB", "Day57", "Delta57overB")
tps_no_fold_change <- times[!grepl("Delta", times)]
tps_no_B_and_fold_change <- times[!grepl("Delta", times) & times!="B"]
tps_delta_over_B <- times[grepl("overB",times)]

#-----------------------------------------------
# PAIR PLOTS
#-----------------------------------------------
# - The correlation of each pair of Day 1, Day tinterm, Day tpeak and Fold-change antibody marker readouts,
#   stratified by treatment group and baseline serostatus
# - Pairs plots/scatterplots and baseline strata-adjusted Spearman rank correlations are used.
#-----------------------------------------------
for (country in c("Nvx_US_Mex", if(study_name=="PREVENT19") "Nvx_US")) { # this loop is only for prevent19, prevent19 needs one set of US+MEX, and one set of US only
    
  if (!length(assay_immuno)==1){ # AZ two datasets only have one marker in each as of 5/13/2022, can't do pair 
    
    print("Pair plots 1:")
    
    for (tp in tps_no_delta_over_tinterm) { # "B", "Day29", "Day57", "Day29overB", "Day57overB"
      for (trt in 0:1) {
        # Don't produce figures for placebo baseline negative to improve build time
        if(trt==0) {bstatus.range <- 1} else {bstatus.range <- unique(dat.twophase.sample$Bserostatus)}
    
        for (bserostatus in bstatus.range) {
          if (!bserostatus %in% unique(dat.twophase.sample$Bserostatus)) next
          
          tt=match(tp, times)
          
          subdat <- dat.twophase.sample %>%
            dplyr::filter(Bserostatus == bserostatus & Trt == trt)
          
          if(study_name=="PREVENT19" & country=="Nvx_US") {subdat=subset(subdat, Country==0)} # Nvx Country: (USA = 0, MEX  = 1)
          
          covid_corr_pairplots(
            plot_dat = subdat,
            time = tp,
            assays = assay_immuno,
            strata = "Bstratum",
            weight = "wt.subcohort",
            plot_title = paste0(
              gsub("ay ","", labels.time)[tt],
              " Ab markers: baseline ",
              ifelse(bserostatus, "positive", "negative"), ", ",
              c("placebo", "vaccine")[trt + 1], " arm"
            ),
            column_labels = labels.axis[tp, seq_along(assay_immuno)] %>% unlist(),
            height = max(1.3 * length(assay_immuno) + 0.1, 5.5),
            width = max(1.3 * length(assay_immuno), 5.5),
            column_label_size = ifelse(max(str_length(assay_labels_short)) > 28, 5, 6.5),
            filename = paste0(
              save.results.to, "/pairs_", tp,
              "_Markers_", bstatus.labels.2[bserostatus + 1],
              c("_placebo_arm", "_vaccine_arm")[trt + 1], "_", ifelse(country=="Nvx_US", "US_only_",""),
              study_name, ".pdf"
            )
          )
        }
      }
    }
  }
    
  ## pairplots of assay readouts for multiple timepoints
  ## pairplots by baseline serostatus
  print("Pair plots 2:")
  for (trt in 0:1) {
    # Don't produce figures for placebo baseline negative to improve build time
    if(trt==0) {bstatus.range <- 1} else {bstatus.range <- unique(dat.twophase.sample$Bserostatus)}
  
    for (bserostatus in bstatus.range) {
      if (!bserostatus %in% unique(dat.twophase.sample$Bserostatus)) next
      
      subdat <- dat.twophase.sample %>%
        dplyr::filter(Bserostatus == bserostatus & Trt == trt)
      
      times_selected <- if(study_name=="VAT08m") {tps_no_delta_over_tinterm 
        # "B", "Day29", "Day57", "Day29overB", "Day57overB", only show fold_change for Sanofi study
        } else {tps_no_fold_change} # "B", "Day29", "Day57"
      
      if(study_name=="PREVENT19" & country=="Nvx_US") {subdat=subset(subdat, Country==0)} # Nvx Country: (USA = 0, MEX  = 1)
      
      for (aa in assay_immuno) {
        covid_corr_pairplots_by_time(
          plot_dat = subdat,
          times = times_selected,
          assay = aa,
          strata = "Bstratum",
          weight = "wt.subcohort",
          plot_title = paste0(
            labels.assays[aa], ": baseline ",
            ifelse(bserostatus, "positive ", "negative "),
            c("placebo", "vaccine")[trt + 1], " arm"
          ),
          column_labels = paste(gsub("ay ","", labels.time[times_selected]),
                                "\n", labels.axis[, aa][1]),
          column_label_size = ifelse(study_name=="VAT08m", 4.5, 
                                     ifelse(max(str_length(assay_labels_short)) > 28, 5, 6.5)),
          axis_label_size = ifelse(study_name=="VAT08m", 7, 9),
          filename = paste0(
            save.results.to, "/pairs_", aa, "_by_times_",
            bstatus.labels.2[bserostatus + 1], "_", c("placebo_", "vaccine_")[trt + 1], ifelse(country=="Nvx_US", "US_only_",""),
            study_name, ".pdf"
          )
        )
      }
    }
  }
}
  
#-----------------------------------------------
# - Reverse empirical cdf (rcdf) plots for the Day tinterm, Day tpeak, and Fold-change antibody marker readouts,
#   stratified by treatment group and baseline serostatus
# - We made multiple ggplot objects, each for one assay, and combine them with ggarrange()
#-----------------------------------------------
print("RCDF 1:")
for (tp in tps_no_B_and_delta_over_tinterm) { # "Day29", "Day57", "Day29overB", "Day57overB"
  covid_corr_rcdf_facets(
    plot_dat = dat.long.twophase.sample,
    x = tp,
    facet_by = "assay",
    color = "trt_bstatus_label",
    weight = "wt.subcohort",
    xlim = assay_lim[rep(assay_immuno, ifelse(length(assay_immuno)==1, 2, 1)), tp, ], # call the same marker twice if only one marker exists
    arrange_ncol = 3,
    arrange_nrow = ceiling(length(assay_immuno) / 3),
    panel_titles = labels.title2[tp, ] %>% unlist(),
    axis_titles = labels.axis[tp, ] %>% unlist(),
    filename = paste0(
      save.results.to, "/Marker_Rcdf_", tp,
      "_trt_both_bstatus_both_", study_name, ".pdf"
    )
  )
}

#-----------------------------------------------
# RCDF plot of multiple bAb and nAb assay readouts in one plot, 
# one plot for each Day tinterm and Day tpeak
# line-types distinguishing the baseline serostatus
# one type with both baseline serostatus in the same plot, one type in the different plot
#-----------------------------------------------
dat.long.twophase.sample$assay_labels <-
  factor(dat.long.twophase.sample$assay,
         levels = assay_immuno,
         labels = labels.assays.short)
print("RCDF 2:")
# plot bAb, PsV and ADCP assays separately
for (Ab in c(1, 2, 3)) {
  
  if (Ab == 1) {
    rcdf_assays <- assay_immuno[grepl("bind", assay_immuno)]
  } else if (Ab == 2) {
    rcdf_assays <- assay_immuno[grepl("neut", assay_immuno)]
  } else  {
    rcdf_assays <- assay_immuno[grepl("ADCP", assay_immuno)]
  }
  
  if (length(rcdf_assays) > 0) {
    
    for (tp in tps_no_B_and_delta_over_tinterm) { # "Day29", "Day57", "Day29overB", "Day57overB"
      covid_corr_rcdf(
        plot_dat = subset(dat.long.twophase.sample, Trt == "Vaccine" & assay %in% rcdf_assays),
        x = tp,
        color = "assay_labels",
        lty = "Bserostatus",
        weight = "wt.subcohort",
        xlab = paste0(gsub("ay ", "", labels.time[tp]), " Ab Markers"
        ),
        xlim = c(min(assay_lim[rcdf_assays, tp, 1]), 
                 max(assay_lim[rcdf_assays, tp, 2])),
        xbreaks = seq(min(assay_lim[rcdf_assays, tp, 1]), 
                      max(assay_lim[rcdf_assays, tp, 2]), 
                      1),
        plot_title = paste0(labels.time[tp], " Ab Markers"),
        filename = paste0(
          save.results.to, "/Marker_Rcdf_", c("bAb", "nAb", "other")[Ab], "_", tp,
          "_trt_vaccine_bstatus_both_", study_name, ".pdf"
        )
      )
    }
    
    #-----------------------------------------------
    # RCDF plot 
    # different baseline serostatus in different plot
    #-----------------------------------------------
    print("RCDF 3:")
    for (bstatus in 1:2) {
      if (nrow(subset(dat.long.twophase.sample, Bserostatus==bstatus.labels[bstatus]))==0) next 
      
      for (tp in tps_no_B_and_delta_over_tinterm) { # "Day29", "Day57", "Delta29overB", "Delta57overB"
        covid_corr_rcdf(
          plot_dat = filter(dat.long.twophase.sample, Trt == "Vaccine", 
                            Bserostatus == bstatus.labels[bstatus],
                            assay %in% rcdf_assays),
          x = tp,
          color = "assay_labels",
          lty = NULL,
          weight = "wt.subcohort",
          xlab = paste0(gsub("ay ", "", labels.time[tp]), " Ab Markers"
          ),
          xlim = c(min(assay_lim[rcdf_assays, tp, 1]), 
                   max(assay_lim[rcdf_assays, tp, 2])),
          xbreaks = seq(min(assay_lim[rcdf_assays, tp, 1]), 
                        max(assay_lim[rcdf_assays, tp, 2]), 
                        1),
          plot_title = paste0(labels.time[tp], " Ab Markers"
          ),
          filename = paste0(
            save.results.to, "/Marker_Rcdf_", c("bAb", "nAb", "other")[Ab], "_", tp,
            "_trt_vaccine_bstatus_", c("Neg", "Pos")[bstatus], "_", study_name, ".pdf"
          )
        )
      }
    }
  }
}

#-----------------------------------------------
# BOX PLOTS
#-----------------------------------------------
# - Box plots across treatment groups.
# - made a ggplot object for every assay and use the ggarrange() function to combine the resulted plots.
# - boxplots of assay readouts at Day 1, Day tinterm and Day tpeak, 
# - by treatment groups or by baseline serostatus
#-----------------------------------------------
print("Boxplots 1:")
for (bstatus in 1:2) {
  if (nrow(subset(dat.long.twophase.sample, Bserostatus==bstatus.labels[bstatus]))==0) next 
  
  for (tp in tps_no_delta_over_tinterm) { # "B", "Day29", "Day57", "Delta57overB", "Delta29overB"

    covid_corr_boxplot_facets(
      plot_dat = subset(
        dat.long.twophase.sample,
        Bserostatus == bstatus.labels[bstatus]
      ),
      x = "Trt",
      y = tp,
      color = "Trt",
      facet_by = "assay",
      ylim = assay_lim[, tp, ],
      plot_LLOX = !grepl("Delta", tp), # "B", "Day29", "Day57"
      POS.CUTOFFS = log10(pos.cutoffs[assay_immuno]),
      LLOX = log10(lloxs[assay_immuno]),
      ULOQ = log10(uloqs[assay_immuno]),
      arrange_ncol = 3,
      arrange_nrow = ceiling(length(assay_immuno) / 3),
      legend = c("Placebo", "Vaccine"),
      axis_titles_y = labels.axis[tp, ] %>% unlist(),
      panel_titles = labels.title2[tp, ] %>% unlist(),
      filename = paste0(
        save.results.to, "/boxplots_", tp, "_x_trt_", bstatus.labels.2[bstatus],
        "_", study_name, ".pdf"
      )
    )
  }
}

#-----------------------------------------------
# - Box plots of the assay readouts versus baseline sero-status, stratified by treatment groups
# - Make seperate plots for Placebo and Vaccine arms
#-----------------------------------------------
for (trt in 1:2) {
  for (tp in tps_no_delta_over_tinterm) {
    
    covid_corr_boxplot_facets(
      plot_dat = subset(dat.long.twophase.sample, as.numeric(Trt) == trt),
      x = "Bserostatus",
      y = tp,
      color = "Bserostatus",
      facet_by = "assay",
      ylim = assay_lim[, tp,],
      plot_LLOX = !grepl("Delta", tp), # "B", "Day29", "Day57"
      POS.CUTOFFS = log10(pos.cutoffs[assay_immuno]),
      LLOX = log10(lloxs[assay_immuno]),
      ULOQ = log10(uloqs[assay_immuno]),
      arrange_ncol = 3,
      arrange_nrow = ceiling(length(assay_immuno) / 3),
      legend = c("Baseline Negative", "Baseline Positive"),
      axis_titles_y = labels.axis[tp, ] %>% unlist(),
      panel_titles = labels.title2[tp, ] %>% unlist(),
      filename = paste0(
        save.results.to, "/boxplots_", tp,
        "_x_bstatus_", c("placebo_arm_", "vaccine_arm_")[trt],
        study_name, ".pdf"
      )
    )
  }
}

#-----------------------------------------------
# Spaghetti PLOTS
#-----------------------------------------------
# - Spaghetti plots of antibody marker change over time
#-----------------------------------------------
print("Spaghetti plots:")
## in each baseline serostatus group, randomly select 10 placebo recipients and 20 vaccine recipients
set.seed(12345)
var_names <- expand.grid(times = times[!grepl("Delta",times)], # "B", "Day29", "Day57" 
                         assays = assay_immuno) %>%
  mutate(var_names = paste0(times, assay_immuno)) %>%
  .[, "var_names"]

spaghetti_ptid <- dat.twophase.sample[, c("Ptid", "Bserostatus", "Trt", var_names)] %>%
  filter(., complete.cases(.)) %>%
  transmute(BT = paste0(as.character(Bserostatus), as.character(Trt)),
            Ptid = Ptid) %>%
  split(., .$BT) %>%
  lapply(function(xx) {
    if (xx$BT[1] %in% c("10", "00")) {
      sample(xx$Ptid, 20, ifelse(length(xx$Ptid)<20, T, F))  ## sample 10 placebo recipients
      # add ifelse(length(xx$Ptid)<20, F, T) because some subset has small sample e.g. janssen_sa_partA
    } else {
      sample(xx$Ptid, 20, ifelse(length(xx$Ptid)<20, T, F))  ## sample 20 vaccine recipients
    }
  }) %>% unlist %>% as.character

spaghetti_dat <- dat.long.twophase.sample[, c("Ptid", "Bserostatus", "Trt", "assay",
                                              times[!grepl("Delta",times)] # "B", "Day29", "Day57"
                                              )] %>%
  filter(Ptid %in% spaghetti_ptid) %>%
  pivot_longer(cols = times[!grepl("Delta",times)], # "B", "Day29", "Day57"
               names_to = "time") %>%
  mutate(assay = factor(assay, levels = assay_immuno, labels = assay_immuno),
         time_label = factor(time, levels = times[!grepl("Delta",times)], # "B", "Day29", "Day57"
                             labels = gsub("B","D1",gsub("ay","", times[!grepl("Delta",times)]))# "D1", "D29", "D57"
                             )) %>%
  as.data.frame

for (bstatus in 1:2) {
  subdat <- subset(spaghetti_dat, Bserostatus == bstatus.labels[bstatus])
  if(nrow(subdat)==0) next
  covid_corr_spaghetti_facets(plot_dat = subdat,
                              x = "time_label",
                              y = "value",
                              id = "Ptid",
                              color = "Trt",
                              facet_by = "assay",
                              ylim = assay_lim[, times[!grepl("B|Delta",times)][1] # "Day29", "Day57"
                                               ,],
                              panel_titles = labels.assays.short,
                              plot_title = paste0(
                                "Baseline ",
                                c("Negative", "Positive")[bstatus],
                                " PP Placebo + Vaccine group"
                              ),
                              arrange_ncol = 3,
                              arrange_nrow = ceiling(length(assay_immuno) / 3),
                              filename = paste0(
                                save.results.to, "/spaghetti_plot_",
                                bstatus.labels.2[bstatus], "_",
                                study_name, ".pdf"
                              ))
}

#-----------------------------------------------
# Scatter PLOTS
#-----------------------------------------------
# - Scatter plots assay vs. age in years, (Day 1) Day tinterm, Day tpeak
#-----------------------------------------------
print("Scatter plots:")
for (tp in tps_no_fold_change) {
  for (trt in 1:2) {
    for (bstatus in 1:2) {
      
      subdat <- dat.long.twophase.sample %>%
        filter(Bserostatus == bstatus.labels[bstatus], Trt == trt.labels[trt])
      if(nrow(subdat)==0) next
      
      ## setting the range of the axes
      xrange <- range(dat.long.twophase.sample$Age)
      
      
      scatter_plot_list <- vector("list", length = length(assay_immuno))
      
      for (aa in 1:length(assay_immuno)) {
        scatter_plot_list[[aa]] <- ggplot(data = subset(subdat, assay == assay_immuno[aa]),
                                          mapping = aes_string("Age", tp)) +
          geom_point() +
          xlab("Age") +
          ylab(labels.axis[tp, aa]) +
          ggtitle(labels.title[tp, aa]) +
          stat_smooth(method = "loess", color = "red", se = TRUE, lwd = 1) +
          scale_x_continuous(limits = xrange) +
          scale_y_continuous(
            labels = label_math(10^.x), limits = assay_lim[aa, tp,],
            breaks = seq(assay_lim[aa, tp, 1], assay_lim[aa, tp, 2], by = 1)
          ) +
          theme_pubr() +
          theme(
            plot.title = element_text(hjust = 0.5, size = ifelse(max(str_length(assays)) > 14, 9.5, 10)),
            panel.border = element_rect(fill = NA),
            panel.grid.minor.y = element_line(),
            panel.grid.major.y = element_line(),
            axis.title = element_text(size = 10),
            axis.text = element_text(size = 10),
            legend.title = element_blank()
          )
      }
      
      output_plot <- ggarrange(
        plotlist = scatter_plot_list, ncol = 3,
        nrow = ceiling(length(assay_immuno) / 3), legend = "none", align = "h"
      )
      
      ggsave(
        filename = paste0(
          save.results.to, "/scatter_", tp, "_vs_age_",
          "trt_", trt.labels[trt], "_", bstatus.labels.2[bstatus],
          "_", study_name, ".pdf"
        ), 
        plot = output_plot, 
        width = 9,
        height = 0.5 + 3 * ceiling(length(assay_immuno) / 3), 
        units = "in"
      )
      
    }
  }
}
