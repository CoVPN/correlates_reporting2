#Sys.setenv(TRIAL = "vat08_combined")
#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
# There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))
Sys.setenv(DESCRIPTIVE = 1)
source(here::here("..", "_common.R"))
uloqs=assay_metadata$uloq; names(uloqs)=assays
pos.cutoffs=assay_metadata$pos.cutoff; names(pos.cutoffs)=assays
lloxs=ifelse(assay_metadata$llox_label=="lloq", assay_metadata$lloq, assay_metadata$lod); names(lloxs)=assays
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
library(grid) # textGrob
#library(dummies) # this package got archived on 2022-04-29
require(devtools)
install_version("dummies", version = "1.5.6", repos = "http://cran.us.r-project.org")
library(gridExtra)
library(PResiduals)
install.packages("fmsb", repos = "http://cran.us.r-project.org") # radar plot
library(fmsb) # radarchart()
install.packages("wCorr", repos = "http://cran.us.r-project.org") # weighted correlation
library(wCorr)

# produce geom_statistics w/ resampling-based covariate-adjusted Spearman
source(here("code", "params.R"))
if (study_name=="VAT08"){
  source(here("code", "covid_corr_plot_functions.R"))
  source(here("code", "process_violin_pair_functions.R")) # pair functions in the first program are overwritten by those in the second program
  # for VAT08, pairplots are non-bstratum-adjusted, no resampling, IPS-weighted spearman correlation
} else {
  source(here("code", "ggally_cor_resample.R"))
  source(here("code", "covid_corr_plot_functions.R"))
}

set.seed(12345)
# load cleaned data
dat.long.twophase.sample <- readRDS(here(
  "data_clean",
  "long_twophase_data.rds"
))

dat.twophase.sample <- readRDS(here("data_clean", "twophase_data.rds")); dat.twophase.sample$all_one <- 1 # as a placeholder for strata values

assay_lim <- readRDS(here("data_clean", "assay_lim.rds"))

tps_no_delta_over_tinterm <-  times[!times %in% c(paste0("Delta",timepoints[length(timepoints)],"over",timepoints[1]))] #c("B", "Day29", "Delta29overB", "Day57", "Delta57overB")
tps_no_B_and_delta_over_tinterm <-  times[!times %in% c("B",paste0("Delta",timepoints[length(timepoints)],"over",timepoints[1]))] #c("Day29", "Delta29overB", "Day57", "Delta57overB")
tps_no_fold_change <- times[!grepl("Delta", times)]
tps_no_B_and_fold_change <- times[!grepl("Delta", times) & times!="B"]
tps_delta_over_B <- times[grepl("overB",times)]

# adhoc request for profiscov
if (F){
  limits = c(-1.5, 5)
  breaks = seq(-1, 4)
  
  p1 <- ggplot(subset(dat.twophase.sample, Trt==1 & Bserostatus==0), aes(Day43bindSpike_P.1, Day91bindSpike_P.1)) +
    geom_point(size = 0.5) +
    geom_smooth(method = "loess", se=FALSE, color="red") +
    scale_x_continuous(
      limits = limits , breaks = breaks,
      labels = label_math(10^.x)
    ) +
    scale_y_continuous(
      limits = limits , breaks = breaks,
      labels = label_math(10^.x)
    ) +
    xlab("Day 43\nAnti Spike P.1 IgG (BAU/ml)") + 
    ylab("Day 91\nAnti Spike P.1 IgG (BAU/ml)") + 
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    coord_fixed(ratio = 1)
  
  p2 <- ggplot(subset(dat.twophase.sample, Trt==1 & Bserostatus==0), aes(Day43bindRBD_P.1, Day91bindRBD_P.1)) +
    geom_point(size = 0.5) +
    geom_smooth(method = "loess", se=FALSE, color="red") +
    scale_x_continuous(
      limits = limits , breaks = breaks,
      labels = label_math(10^.x)
    ) +
    scale_y_continuous(
      limits = limits , breaks = breaks,
      labels = label_math(10^.x)
    ) +
    xlab("Day 43\nAnti RBD P.1 IgG (BAU/ml)") + 
    ylab("Day 91\nAnti RBD P.1 IgG (BAU/ml)") + 
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    coord_fixed(ratio = 1)
  
  p3 <- ggplot(subset(dat.twophase.sample, Trt==1 & Bserostatus==0), aes(Day43bindN, Day91bindN)) +
    geom_point(size = 0.5) +
    geom_smooth(method = "loess", se=FALSE, color="red") +
    scale_x_continuous(
      limits = limits , breaks = breaks,
      labels = label_math(10^.x)
    ) +
    scale_y_continuous(
      limits = limits , breaks = breaks,
      labels = label_math(10^.x)
    ) +
    xlab("Day 43\nAnti N IgG (BAU/ml)") + 
    ylab("Day 91\nAnti N IgG (BAU/ml)") + 
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    coord_fixed(ratio = 1)
  
  scatters <- grid.arrange(p1, p2, p3, nrow = 1, widths = c(1, 1, 1),
                           top = textGrob("D43 vs D91: baseline negative vaccine arm"))
  ggsave(
    filename = paste0(
      save.results.to, "/scatterplot_D43vsD91_bAbSpikeP1_bAbRBDP1_N.pdf"), plot = scatters, width = 10, height = 3,
    units = "in"
  )
}

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
    
    for (tp in if(study_name!="VAT08") {tps_no_B_and_delta_over_tinterm} else {tps_no_fold_change}) { # "B", "Day29", "Day57", "Day29overB", "Day57overB"
      for (trt in 0:1) {
        # Don't produce figures for placebo baseline negative to improve build time
        if(trt==0 & study_name!="VAT08") {bstatus.range <- 1} else {bstatus.range <- unique(dat.twophase.sample$Bserostatus)}
    
        for (bserostatus in bstatus.range) {
          if (!bserostatus %in% unique(dat.twophase.sample$Bserostatus)) next
          
          tt=match(tp, times)
          
          subdat <- dat.twophase.sample %>%
            dplyr::filter(Bserostatus == bserostatus & Trt == trt)
          
          if(study_name=="PREVENT19" & country=="Nvx_US") {subdat=subset(subdat, Country==0)} # Nvx Country: (USA = 0, MEX  = 1)
          
          if (study_name=="VAT08" && bserostatus==0 && tp=="B") { # psv_mdw doesn't have any value for naive at baseline
            assay_immuno_ = assay_immuno[assay_immuno!="pseudoneutid50_mdw"]
          } else {assay_immuno_ = assay_immuno}
          
          covid_corr_pairplots(
            plot_dat = subdat,
            time = tp,
            assays = assay_immuno_, # adhoc request by David: assay_immuno = c("bindSpike", "bindSpike_P.1", "bindRBD", "bindRBD_P.1", "bindN")
                                   # adhoc request 2 by David: assay_immuno = c("liveneutmn50", "bindSpike_P.1", "bindRBD_P.1", "bindN")
            strata = ifelse(study_name=="VAT08", "all_one", "Bstratum"),
            weight = "wt.subcohort",
            plot_title = paste0(
              gsub("ay ","", labels.time)[tt],
              " Ab markers: ",
              bstatus.labels.3[bserostatus + 1], ", ",
              c("placebo", "vaccine")[trt + 1], " arm"
            ),
            column_labels = labels.axis[tp, seq_along(assay_immuno_)] %>% unlist(), # adhoc request by David: labels.axis[tp, match(assay_immuno, colnames(labels.axis))]
            height = max(1.3 * length(assay_immuno_) + 0.1, 5.5),
            width = max(1.3 * length(assay_immuno_), 5.5),
            column_label_size = ifelse(max(str_length(labels.axis[1,])) > 28, 4, 6.5),
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
    if(trt==0 & study_name!="VAT08") {bstatus.range <- 1} else {bstatus.range <- unique(dat.twophase.sample$Bserostatus)}
  
    for (bserostatus in bstatus.range) {
      if (!bserostatus %in% unique(dat.twophase.sample$Bserostatus)) next
      
      subdat <- dat.twophase.sample %>%
        dplyr::filter(Bserostatus == bserostatus & Trt == trt)
      
      times_selected <- if(study_name=="VAT08") {tps_no_delta_over_tinterm[c(1,4,5)]
        # "B", "Day29", "Day57", "Day29overB", "Day57overB", only show B and fold_change for Sanofi study
        } else {tps_no_fold_change} # "B", "Day29", "Day57"
      
      if(study_name=="PREVENT19" & country=="Nvx_US") {subdat=subset(subdat, Country==0)} # Nvx Country: (USA = 0, MEX  = 1)
      
      for (aa in assay_immuno) {
        if (study_name=="VAT08" && aa=="pseudoneutid50_mdw" && bserostatus==0) next # psv_mdw doesn't have any value for naive at baseline
        
        covid_corr_pairplots_by_time(
          plot_dat = subdat,
          times = times_selected,
          assay = aa,
          strata = ifelse(study_name=="VAT08", "all_one", "Bstratum"),
          weight = "wt.subcohort",
          plot_title = paste0(
            labels.assays[aa], ": ",
            bstatus.labels.3[bserostatus + 1], " ",
            c("placebo", "vaccine")[trt + 1], " arm"
          ), 
          column_labels = paste(gsub("ay ","", labels.time[times_selected]),
                                "\n", labels.axis[, aa][1]),
          column_label_size = ifelse(study_name=="VAT08", 4.5, 
                                     ifelse(max(str_length(labels.axis[1,])) > 28, 4.3, 6.5)),
          axis_label_size = ifelse(study_name=="VAT08", 7, 9),
          filename = paste0(
            save.results.to, "/pairs_", aa, "_by_times_",
            bstatus.labels.2[bserostatus + 1], "_", c("placebo_", "vaccine_")[trt + 1], ifelse(country=="Nvx_US", "US_only_",""),
            study_name, ".pdf"
          )
        )
      }
    }
  }
  
  print("Pair plots 3:")
  
  if (study_name=="VAT08") { # request only for this study, at day 1, pool over vaccine and placebo
    tp = "B"
    for (bserostatus in bstatus.range) {
      
        if (bserostatus==0) { # VAT08 psv_mdw doesn't have any value for naive at baseline
          assay_immuno_ = assay_immuno[assay_immuno!="pseudoneutid50_mdw"]
        } else {assay_immuno_ = assay_immuno}
      
        if (!bserostatus %in% unique(dat.twophase.sample$Bserostatus)) next
        
        tt=match(tp, times)
        
        subdat <- dat.twophase.sample %>%
          dplyr::filter(Bserostatus == bserostatus)
        
        covid_corr_pairplots(
          plot_dat = subdat,
          time = tp,
          assays = assay_immuno_,
          strata = ifelse(study_name=="VAT08", "all_one", "Bstratum"),
          weight = "wt.subcohort",
          plot_title = paste0(
            gsub("ay ","", labels.time)[tt],
            " Ab markers: ",
            bstatus.labels.3[bserostatus + 1], ", pooled arm"
          ),
          column_labels = labels.axis[tp, seq_along(assay_immuno_)] %>% unlist(), # adhoc request by David: labels.axis[tp, match(assay_immuno, colnames(labels.axis))]
          height = max(1.3 * length(assay_immuno_) + 0.1, 5.5),
          width = max(1.3 * length(assay_immuno_), 5.5),
          column_label_size = ifelse(max(str_length(labels.axis[1,])) > 28, 4, 6.5),
          filename = paste0(
            save.results.to, "/pairs_", tp,
            "_Markers_", bstatus.labels.2[bserostatus + 1],
            "_pooled_arm_", ifelse(country=="Nvx_US", "US_only_",""),
            study_name, ".pdf"
          )
        )
      }
  }
}
  
#-----------------------------------------------
# - Reverse empirical cdf (rcdf) plots for the Day tinterm, Day tpeak, and Fold-change antibody marker readouts,
#   stratified by treatment group and baseline serostatus
# - We made multiple ggplot objects, each for one assay, and combine them with ggarrange()
#-----------------------------------------------
print("RCDF 1:")
for (tp in if(study_name!="VAT08") {tps_no_B_and_delta_over_tinterm} else {tps_no_fold_change}) { # "Day29", "Day57", "Day29overB", "Day57overB" for most studies; if VAT08, "Day1", "Day22", "Day43"
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
    xbreaks = ifelse(study_name=="VAT08", 2, 1),
    axis_title_size = 10,
    axis_size = 10,
    panel_title_size = ifelse(study_name=="VAT08", 8, 10),
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
for (Ab in c(1, 2, if(study_name!="VAT08") 3)) {
  
  if (Ab == 1) {
    rcdf_assays <- assay_immuno[grepl("bind", assay_immuno)]
  } else if (Ab == 2) {
    rcdf_assays <- assay_immuno[grepl("neut", assay_immuno)]
  } else  {
    rcdf_assays <- assay_immuno[grepl("ADCP", assay_immuno)]
  }
  
  if (length(rcdf_assays) > 0) {
    
    for (tp in if(study_name!="VAT08") {tps_no_B_and_delta_over_tinterm} else {tps_no_B_and_fold_change}) { # "Day29", "Day57", "Day29overB", "Day57overB" for most studies; if VAT08, "Day22", "Day43"
        
      for (trt in c("Vaccine", if(study_name=="VAT08") "Placebo")){
        covid_corr_rcdf(
          plot_dat = subset(dat.long.twophase.sample, Trt == trt & assay %in% rcdf_assays),
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
                        ifelse(study_name=="VAT08", 3, 1)),
          plot_title = paste0(labels.time[tp], " Ab Markers"),
          filename = paste0(
            save.results.to, "/Marker_Rcdf_", c("bAb", "nAb", "other")[Ab], "_", tp,
            "_trt_", tolower(trt), "_bstatus_both_", study_name, ".pdf"
          )
        )
      }
    }
    
    #-----------------------------------------------
    # RCDF plot 
    # different baseline serostatus in different plot, one treatment arm per plot
    #-----------------------------------------------
    if (study_name!="VAT08"){# VAT08 doesn't need plots for one arm * one baseline status
      print("RCDF 3:")
      for (bstatus in 1:2) {
        if (nrow(subset(dat.long.twophase.sample, Bserostatus==bstatus.labels[bstatus]))==0) next 
        
        for (tp in if(study_name!="VAT08") {tps_no_B_and_delta_over_tinterm} else {tps_no_B_and_fold_change}) { # "Day29", "Day57", "Day29overB", "Day57overB" for most studies; if VAT08, "Day22", "Day43"
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
    
    #-----------------------------------------------
    # RCDF plot 
    # different treatment arm in different plot, one baseline status per plot
    #-----------------------------------------------
    if (study_name=="VAT08"){
      print("RCDF 4:")
      for (bstatus in 1:2) {
        if (nrow(subset(dat.long.twophase.sample, Bserostatus==bstatus.labels[bstatus]))==0) next
        
        for (tp in if(study_name!="VAT08") {tps_no_B_and_delta_over_tinterm} else {tps_no_B_and_fold_change}) { # "Day29", "Day57", "Day29overB", "Day57overB" for most studies; if VAT08, "Day22", "Day43"
          covid_corr_rcdf(
            plot_dat = subset(dat.long.twophase.sample, Bserostatus == bstatus.labels[bstatus] & assay %in% rcdf_assays),
            x = tp,
            color = "assay_labels",
            lty = "Trt",
            weight = "wt.subcohort",
            xlab = paste0(gsub("ay ", "", labels.time[tp]), " Ab Markers"
            ),
            xlim = c(min(assay_lim[rcdf_assays, tp, 1]), 
                     max(assay_lim[rcdf_assays, tp, 2])),
            xbreaks = seq(min(assay_lim[rcdf_assays, tp, 1]), 
                          max(assay_lim[rcdf_assays, tp, 2]), 
                          ifelse(study_name=="VAT08", 3, 1)),
            plot_title = paste0(labels.time[tp], " Ab Markers"),
            filename = paste0(
              save.results.to, "/Marker_Rcdf_", c("bAb", "nAb", "other")[Ab], "_", tp,
              "_trt_both_bstatus_", c("Neg", "Pos")[bstatus], "_", study_name, ".pdf"
            )
          )
        }
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
  
  for (tp in if(study_name!="VAT08") {tps_no_B_and_delta_over_tinterm} else {tps_no_fold_change}) { # "Day29", "Day57", "Day29overB", "Day57overB" for most studies; if VAT08, "Day1", "Day22", "Day43"
    
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
      panel_title_size = ifelse(study_name=="VAT08", 8, 10),
      height = ifelse(study_name=="VAT08", 11, 3 * arrange_nrow + 0.5),
      filename = paste0(
        save.results.to, "/boxplots_", tp, "_x_trt_", bstatus.labels.2[bstatus],
        "_", study_name, ".pdf"
      )
    )
  }
}

#-----------------------------------------------
# - Box plots of the assay readouts versus baseline sero-status, stratified by treatment groups
# - Make separate plots for Placebo and Vaccine arms
#-----------------------------------------------
for (trt in 1:2) {
  for (tp in if (study_name=="VAT08") {tps_no_fold_change} else {tps_no_delta_over_tinterm}) {
    
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
      legend = stringr::str_to_title(bstatus.labels.3),
      axis_titles_y = labels.axis[tp, ] %>% unlist(),
      panel_titles = labels.title2[tp, ] %>% unlist(),
      panel_title_size = ifelse(study_name=="VAT08", 8, 10),
      height = ifelse(study_name=="VAT08", 11, 3 * arrange_nrow + 0.5),
      filename = paste0(
        save.results.to, "/boxplots_", tp,
        "_x_bstatus_", c("placebo_arm_", "vaccine_arm_")[trt],
        study_name, ".pdf"
      )
    )
  }
}


#-----------------------------------------------
# - Box plots of the assay readouts, stratified by baseline sero-status and treatment groups
# - One plot for Placebo and Vaccine arms, Baseline Neg and Pos
#-----------------------------------------------
if (study_name=="VAT08") {# this is only reported for VAT08
  for (tp in tps_no_fold_change) {
    
    covid_corr_boxplot_facets(
      plot_dat = dat.long.twophase.sample %>% mutate(BseroTrt = factor(paste0(Bserostatus,"\n",Trt),
                                                     levels = paste0(rep(bstatus.labels, 2), "\n", rep(c("Vaccine", "Placebo"), each=2))
                                                     )),
      x = "BseroTrt",
      y = tp,
      color = "BseroTrt",
      facet_by = "assay",
      ylim = assay_lim[, tp,],
      plot_LLOX = !grepl("Delta", tp), # "B", "Day29", "Day57"
      POS.CUTOFFS = log10(pos.cutoffs[assay_immuno]),
      LLOX = log10(lloxs[assay_immuno]),
      ULOQ = log10(uloqs[assay_immuno]),
      arrange_ncol = 3,
      arrange_nrow = ceiling(length(assay_immuno) / 3),
      legend = paste0(rep(stringr::str_to_title(bstatus.labels.3), 2), ", ", rep(c("Vaccine", "Placebo"), each=2)),
      axis_titles_y = labels.axis[tp, ] %>% unlist(),
      panel_titles = labels.title2[tp, ] %>% unlist(),
      panel_title_size = ifelse(study_name=="VAT08", 8, 10),
      height = ifelse(study_name=="VAT08", 11, 3 * arrange_nrow + 0.5),
      filename = paste0(
        save.results.to, "/boxplots_", tp,
        "_x_trt_bstatus_",
        study_name, ".pdf"
      )
    )
  }
}


#-----------------------------------------------
# - Box plots of the assay readouts, stratified by baseline sero-status, treatment groups and sex at birth
# - One plot for Placebo and Vaccine arms, Baseline Neg and Pos, Females and Males
#-----------------------------------------------
if (study_name=="VAT08") {# this is only reported for VAT08
  for (tp in tps_no_fold_change) {
    
    covid_corr_boxplot_facets(
      plot_dat = dat.long.twophase.sample %>% mutate(BseroTrtGender = factor(paste0(Bserostatus,"\n",Trt,"\n",sex_label),
                                                                       levels = paste0(rep(bstatus.labels, 2), "\n", rep(c("Vaccine", "Placebo"), each=2) , "\n", rep(c("Female", "Male"), each=4))
      )),
      x = "BseroTrtGender",
      y = tp,
      color = "BseroTrtGender",
      facet_by = "assay",
      ylim = assay_lim[, tp,],
      plot_LLOX = !grepl("Delta", tp), # "B", "Day29", "Day57"
      POS.CUTOFFS = log10(pos.cutoffs[assay_immuno]),
      LLOX = log10(lloxs[assay_immuno]),
      ULOQ = log10(uloqs[assay_immuno]),
      arrange_ncol = 3,
      arrange_nrow = ceiling(length(assay_immuno) / 3),
      legend = paste0(rep(stringr::str_to_title(bstatus.labels.3), 2), ", ", rep(c("Vaccine", "Placebo"), each=2) , ", ", rep(c("Female", "Male"), each=4)),
      axis_titles_y = labels.axis[tp, ] %>% unlist(),
      panel_titles = labels.title2[tp, ] %>% unlist(),
      panel_title_size = ifelse(study_name=="VAT08", 8, 10),
      height = ifelse(study_name=="VAT08", 11, 3 * arrange_nrow + 0.5),
      filename = paste0(
        save.results.to, "/boxplots_", tp,
        "_x_trt_bstatus_gender_",
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
if (study_name!="VAT08"){ # no spaghetti plots for VAT08

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
  if(study_name!="VAT08"){
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
  }
}


print("Spider plots:")
if(study_name=="VAT08"){
  
  ## load data 
  dat.spider <- readRDS(here::here("data_clean", "twophase_data.rds"))
  
  # spider plot showing geometric means calculated using IPS weighting, by trt and baseline status
  
  # calculate geometric mean of IPS weighted readouts
  dat.spider.by.time <- dat.spider %>%
    select(one_of(paste0("B", assays), "Bserostatus", "Trt", "wt.subcohort")) %>%
    rename_with(~str_remove(., "^B")) %>%
    mutate(time="B") %>%
    rename(Bserostatus=serostatus) %>%
    bind_rows(
      dat.spider %>%
        select(one_of(paste0("Day22", assays), "Bserostatus", "Trt", "wt.subcohort")) %>%
        rename_with(~str_remove(., "^Day22")) %>%
        mutate(time="Day22")
    ) %>%
    bind_rows(
      dat.spider %>%
        select(one_of(paste0("Day43", assays), "Bserostatus", "Trt", "wt.subcohort")) %>%
        rename_with(~str_remove(., "^Day43")) %>%
        mutate(time="Day43")
    ) %>%
    bind_rows(
      dat.spider %>%
        select(one_of(paste0("Delta22overB", assays), "Bserostatus", "Trt", "wt.subcohort")) %>%
        rename_with(~str_remove(., "^Delta22overB")) %>%
        mutate(time="Delta22overB")
    ) %>%
    bind_rows(
      dat.spider %>%
        select(one_of(paste0("Delta43overB", assays), "Bserostatus", "Trt", "wt.subcohort")) %>%
        rename_with(~str_remove(., "^Delta43overB")) %>%
        mutate(time="Delta43overB")
    ) %>%
    mutate(Bserostatus = ifelse(Bserostatus == 1, "Pos", "Neg"),
           Trt = ifelse(Trt == 1, "vaccine", "placebo")) %>%
    group_by(time, Bserostatus, Trt) %>%
    summarise(across(assays, ~ exp(sum(log(.x * wt.subcohort), na.rm=T) / sum(wt.subcohort)))) %>%
    unique() %>%
    as.data.frame()
  
  # those without any data will have a weighted geomean equal to 1, set these to NA
  dat.spider.by.time[dat.spider.by.time == 1] <- NA
  
  rownames(dat.spider.by.time) <- paste0(dat.spider.by.time$time, dat.spider.by.time$Bserostatus, dat.spider.by.time$Trt)
  dat.spider.by.time$time <- NULL
  dat.spider.by.time$Bserostatus <- NULL
  dat.spider.by.time$Trt <- NULL
  
  # stack with max and min values
  max_min <- rbind(rep(1.8,ncol(dat.spider.by.time)), 
                   rep(0,ncol(dat.spider.by.time)))
  colnames(max_min) <- colnames(dat.spider.by.time)
  rownames(max_min) <- c("max", "min")
  
  dat.spider.by.time <- rbind(max_min, 
                              dat.spider.by.time)
  
  # setup pdf file
  for (ab in c("bAb", "nAb")) {
    
    for (time in c("day1day22day43", "delta")) {
      
      filename = paste0(save.results.to, "/radar_plot_weighted_geomean_", time, "_", ab, ".pdf")
      pdf(filename, width=5.5, height=6.5)
      par(mfrow=c(2,2), mar=c(0.1,0.1,1,0.1))
      
      for (bsero in c("Neg", "Pos")){
        for(trt in c("placebo", "vaccine")){
          
          #filename = paste0(save.results.to, "/radar_plot_weighted_geomean_", ab, "_trt_", trt, "_bstatus_", bsero, ".pdf")
          #pdf(filename, width=5.5, height=6)
          
          dat.plot <- dat.spider.by.time[c(1,2),] %>%
            bind_rows(dat.spider.by.time[grepl(paste0(bsero, trt), rownames(dat.spider.by.time)),]) %>%
            select(starts_with(ifelse(ab=="bAb", "bindSpike", "pseudoneutid50")))
          
          colnames(dat.plot) <- assay_metadata$assay_label[match( colnames(dat.plot) , assay_metadata$assay)]
          
          colnames(dat.plot) <- gsub("PsV Neutralization to |PsV Neutralization |Binding Antibody to Spike ", "", colnames(dat.plot))
          
          if (time == "day1day22day43") {
            dat.plot.sub = dat.plot[1:5, ]
            color = c("#0AB7C9","#FF6F1B","#FF5EBF")
            legend_lb = c("B","Day22","Day43")
          } else {
            dat.plot.sub = dat.plot[c(1,2,6,7), ]
            color = c("dodgerblue","chartreuse3")
            legend_lb = c("Delta22overB","Delta43overB")}
          
          radarchart(dat.plot.sub, 
                     axistype=1 , 
                     # Customize the polygon
                     pcol = scales::alpha(color, 0.7), plwd=1.5, pty=c(15), plty=1,
                     pfcol = scales::alpha(color, 0.2),
                     #custom the grid
                     cglcol="grey", cglty=1, axislabcol="grey", cglwd=0.8, caxislabels=paste0("10^",seq(0.1,1.7,0.4)), 
                     #label size
                     vlcex=0.4,
                     #title
                     title=paste0("GeoMean ", ifelse(ab=="bAb", "of bAb Markers, ", "of nAb Markers, "), 
                                  ifelse(bsero=="Neg", "naive ", "non-naive "),
                                  trt),
                     #title size
                     cex.main=0.7)
          
          #par(xpd=NA)
          
          #legend
          legend("bottom", legend=legend_lb, lty=5, pch=c(15),
                 col=color, bty="n", ncol=3, cex=0.7,
                 inset=c(-0.25,0))
          
          #dev.off()
        }
      }
      par(xpd=NA)
      dev.off()
    }
  }
  
}  
