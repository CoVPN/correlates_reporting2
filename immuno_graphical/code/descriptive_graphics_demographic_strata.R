#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
source(here::here("..", "_common.R"))
#-----------------------------------------------

library(here)
library(dplyr)
library(stringr)
library(GGally)
library(ggpubr)
#renv::install("SWIM")
library(SWIM)
library(ggplot2)
library(scales)

set.seed(12345)
# load helper scripts and parameters
source(here("code", "ggally_cor_resample.R"))
source(here("code", "covid_corr_plot_functions.R"))
source(here("code", "params.R"))



# load cleaned data
dat.long.twophase.sample <- readRDS(here(
  "data_clean",
  "long_twophase_data.rds"
))
dat.twophase.sample <- readRDS(here("data_clean", "twophase_data.rds"))

assay_lim <- readRDS(here("data_clean", "assay_lim.rds"))
## ============================================================================
## boxplots and weighted rcdf plots of assay readouts at time points versus
##  (1)  age >= 65 or age < 65
##  (2)  at risk / not at risk
##  (3)  age * high risk
##  (4)  sex at birth
##  (5)  age * sex at birth
##  (6)  ethnicity
##  (7)  race
##  (8)  minority status
##  (9)  age * minority status
##  (10) HIV positivity
##  (11) Country
## plot for each treatment group by baseline status
## ============================================================================

for (tp in times[!times %in% c("B",paste0("Delta",timepoints[length(timepoints)],"over",timepoints[1]))]) {  #c("Day29", "Day57", "Delta29overB", "Delta57overB")
  for (trt in 1:2) {
    # Don't produce figures for placebo baseline negative to improve build time
    if(trt==1) {bstatus.range <- 2} else {bstatus.range <- unique(dat.twophase.sample$Bserostatus)+1}
    for (bstatus in bstatus.range) {
      
      subdat <- subset(
        dat.long.twophase.sample,
        Bserostatus == bstatus.labels[bstatus] &
          Trt == trt.labels[trt]
      )
      
      if (nrow(subdat)==0) next

      ##  (1) age >= 65 or age < 65
      covid_corr_boxplot_facets(
        plot_dat = subdat,
        x = "age_geq_65_label",
        y = tp,
        facet_by = "assay",
        plot_LLOX = !grepl("Delta", tp),
        POS.CUTOFFS = log10(pos.cutoffs[assay_immuno]),
        LLOX = log10(lloxs[assay_immuno]),
        ULOQ = log10(uloqs[assay_immuno]),
        ylim = assay_lim[rep(assay_immuno, ifelse(length(assay_immuno)==1, 2, 1)), tp, ], # call the same marker twice if only one marker exists
        axis_titles_y = labels.axis[tp, ] %>% unlist(),
        panel_titles = labels.title2[tp, ] %>% unlist(),
        arrange_ncol = 3,
        arrange_nrow = ceiling(length(assay_immuno) / 3),
        filename = paste0(
          save.results.to,
          "/demographics/boxplots_",
          tp, "_",
          bstatus.labels.2[bstatus],
          "_trt_", trt.labels[trt],
          "_by_age_group_",
          study_name, ".pdf"
        )
      )

      covid_corr_rcdf_facets(
        plot_dat = subdat,
        x = tp,
        facet_by = "assay",
        xlim = assay_lim[rep(assay_immuno, ifelse(length(assay_immuno)==1, 2, 1)), tp, ], # call the same marker twice if only one marker exists
        color = "age_geq_65_label",
        weight = "wt.subcohort",
        panel_titles = labels.title2[tp, ] %>% unlist(),
        axis_titles = labels.axis[tp, ] %>% unlist(),
        arrange_nrow = ceiling(length(assay_immuno) / 3),
        arrange_ncol = 3,
        filename = paste0(
          save.results.to,
          "/demographics/Marker_Rcdf_",
          tp, "_trt_",
          c("placebo_", "vaccine_")[trt],
          bstatus.labels.2[bstatus],
          "_by_age_group_", study_name,
          ".pdf"
        )
      )

      ##  (2) at risk / not at risk
      covid_corr_boxplot_facets(
        plot_dat = subdat,
        x = "highrisk_label",
        y = tp,
        facet_by = "assay",
        ylim = assay_lim[rep(assay_immuno, ifelse(length(assay_immuno)==1, 2, 1)), tp, ], # call the same marker twice if only one marker exists
        plot_LLOX = !grepl("Delta", tp),
        POS.CUTOFFS = log10(pos.cutoffs[assay_immuno]),
        LLOX = log10(lloxs[assay_immuno]),
        ULOQ = log10(uloqs[assay_immuno]),
        axis_titles_y = labels.axis[tp, ] %>% unlist(),
        panel_titles = labels.title2[tp, ] %>% unlist(),
        arrange_nrow = ceiling(length(assay_immuno) / 3),
        arrange_ncol = 3,
        filename = paste0(
          save.results.to,
          "/demographics/boxplots_",
          tp, "_",
          bstatus.labels.2[bstatus],
          "_trt_", trt.labels[trt],
          "_by_risk_group_",
          study_name, ".pdf"
        )
      )

      covid_corr_rcdf_facets(
        plot_dat = subdat,
        x = tp,
        facet_by = "assay",
        xlim = assay_lim[rep(assay_immuno, ifelse(length(assay_immuno)==1, 2, 1)), tp, ], # call the same marker twice if only one marker exists
        color = "highrisk_label",
        weight = "wt.subcohort",
        panel_titles = labels.title2[tp, ] %>% unlist(),
        axis_titles = labels.axis[tp, ] %>% unlist(),
        arrange_nrow = ceiling(length(assay_immuno) / 3),
        arrange_ncol = 3,
        filename = paste0(
          save.results.to,
          "/demographics/Marker_Rcdf_",
          tp, "_trt_",
          c("placebo_", "vaccine_")[trt],
          bstatus.labels.2[bstatus],
          "_by_risk_group_", study_name,
          ".pdf"
        )
      )

      ##  (3) age * high risk
      covid_corr_boxplot_facets(
        plot_dat = subdat,
        x = "age_risk_label",
        y = tp,
        facet_by = "assay",
        ylim = assay_lim[rep(assay_immuno, ifelse(length(assay_immuno)==1, 2, 1)), tp, ], # call the same marker twice if only one marker exists
        plot_LLOX = !grepl("Delta", tp),
        POS.CUTOFFS = log10(pos.cutoffs[assay_immuno]),
        LLOX = log10(lloxs[assay_immuno]),
        ULOQ = log10(uloqs[assay_immuno]),
        axis_titles_y = labels.axis[tp, ] %>% unlist(),
        panel_titles = labels.title2[tp, ] %>% unlist(),
        arrange_nrow = ceiling(length(assay_immuno) / 3),
        arrange_ncol = 3,
        filename = paste0(
          save.results.to,
          "/demographics/boxplots_",
          tp, "_",
          bstatus.labels.2[bstatus],
          "_trt_", trt.labels[trt],
          "_by_age_x_risk_group_",
          study_name, ".pdf"
        )
      )

      covid_corr_rcdf_facets(
        plot_dat = subdat,
        x = tp,
        facet_by = "assay",
        xlim = assay_lim[rep(assay_immuno, ifelse(length(assay_immuno)==1, 2, 1)), tp, ], # call the same marker twice if only one marker exists
        color = "age_risk_label",
        weight = "wt.subcohort",
        panel_titles = labels.title2[tp, ] %>% unlist(),
        axis_titles = labels.axis[tp, ] %>% unlist(),
        arrange_nrow = ceiling(length(assay_immuno) / 3),
        arrange_ncol = 3,
        filename = paste0(
          save.results.to,
          "/demographics/Marker_Rcdf_",
          tp, "_trt_",
          c("placebo_", "vaccine_")[trt],
          bstatus.labels.2[bstatus],
          "_by_age_risk_group_",
          study_name, ".pdf"
        )
      )

      ##  (4) sex at birth
      covid_corr_boxplot_facets(
        plot_dat = subdat,
        x = "sex_label",
        y = tp,
        facet_by = "assay",
        ylim = assay_lim[rep(assay_immuno, ifelse(length(assay_immuno)==1, 2, 1)), tp, ], # call the same marker twice if only one marker exists
        plot_LLOX = !grepl("Delta", tp),
        POS.CUTOFFS = log10(pos.cutoffs[assay_immuno]),
        LLOX = log10(lloxs[assay_immuno]),
        ULOQ = log10(uloqs[assay_immuno]),
        axis_titles_y = labels.axis[tp, ] %>% unlist(),
        panel_titles = labels.title2[tp, ] %>% unlist(),
        arrange_nrow = ceiling(length(assay_immuno) / 3),
        arrange_ncol = 3,
        filename = paste0(
          save.results.to,
          "/demographics/boxplots_",
          tp, "_",
          bstatus.labels.2[bstatus],
          "_trt_", trt.labels[trt],
          "_by_sex_", study_name,
          ".pdf"
        )
      )

      covid_corr_rcdf_facets(
        plot_dat = subdat,
        x = tp,
        facet_by = "assay",
        xlim = assay_lim[rep(assay_immuno, ifelse(length(assay_immuno)==1, 2, 1)), tp, ], # call the same marker twice if only one marker exists
        color = "sex_label",
        weight = "wt.subcohort",
        panel_titles = labels.title2[tp, ] %>% unlist(),
        axis_titles = labels.axis[tp, ] %>% unlist(),
        arrange_nrow = ceiling(length(assay_immuno) / 3),
        arrange_ncol = 3,
        filename = paste0(
          save.results.to,
          "/demographics/Marker_Rcdf_",
          tp, "_trt_",
          c("placebo_", "vaccine_")[trt],
          bstatus.labels.2[bstatus],
          "_by_sex_group_", study_name,
          ".pdf"
        )
      )

      ##  (5) age * sex at birth
      covid_corr_boxplot_facets(
        plot_dat = subdat,
        x = "age_sex_label",
        y = tp,
        facet_by = "assay",
        ylim = assay_lim[rep(assay_immuno, ifelse(length(assay_immuno)==1, 2, 1)), tp, ], # call the same marker twice if only one marker exists
        plot_LLOX = !grepl("Delta", tp),
        POS.CUTOFFS = log10(pos.cutoffs[assay_immuno]),
        LLOX = log10(lloxs[assay_immuno]),
        ULOQ = log10(uloqs[assay_immuno]),
        axis_titles_y = labels.axis[tp, ] %>% unlist(),
        panel_titles = labels.title2[tp, ] %>% unlist(),
        arrange_nrow = ceiling(length(assay_immuno) / 3),
        arrange_ncol = 3,
        filename = paste0(
          save.results.to,
          "/demographics/boxplots_",
          tp, "_",
          bstatus.labels.2[bstatus],
          "_trt_", trt.labels[trt],
          "_by_age_x_sex_group_",
          study_name, ".pdf"
        )
      )

      covid_corr_rcdf_facets(
        plot_dat = subdat,
        x = tp,
        facet_by = "assay",
        xlim = assay_lim[rep(assay_immuno, ifelse(length(assay_immuno)==1, 2, 1)), tp, ], # call the same marker twice if only one marker exists
        color = "age_sex_label",
        weight = "wt.subcohort",
        panel_titles = labels.title2[tp, ] %>% unlist(),
        axis_titles = labels.axis[tp, ] %>% unlist(),
        arrange_nrow = ceiling(length(assay_immuno) / 3),
        arrange_ncol = 3,
        filename = paste0(
          save.results.to,
          "/demographics/Marker_Rcdf_",
          tp, "_trt_",
          c("placebo_", "vaccine_")[trt],
          bstatus.labels.2[bstatus],
          "_by_age_sex_group_",
          study_name, ".pdf"
        )
      )

      # (6) ethnicity
      covid_corr_boxplot_facets(
        plot_dat = subdat,
        x = "ethnicity_label",
        y = tp,
        facet_by = "assay",
        ylim = assay_lim[rep(assay_immuno, ifelse(length(assay_immuno)==1, 2, 1)), tp, ], # call the same marker twice if only one marker exists
        plot_LLOX = !grepl("Delta", tp),
        POS.CUTOFFS = log10(pos.cutoffs[assay_immuno]),
        LLOX = log10(lloxs[assay_immuno]),
        ULOQ = log10(uloqs[assay_immuno]),
        axis_titles_y = labels.axis[tp, ] %>% unlist(),
        panel_titles = labels.title2[tp, ] %>% unlist(),
        arrange_nrow = ceiling(length(assay_immuno) / 3),
        arrange_ncol = 3,
        filename = paste0(
          save.results.to,
          "/demographics/boxplots_",
          tp, "_",
          bstatus.labels.2[bstatus],
          "_trt_", trt.labels[trt],
          "_by_ethnicity_", study_name,
          ".pdf"
        )
      )

      covid_corr_rcdf_facets(
        plot_dat = subdat,
        x = tp,
        facet_by = "assay",
        xlim = assay_lim[rep(assay_immuno, ifelse(length(assay_immuno)==1, 2, 1)), tp, ], # call the same marker twice if only one marker exists
        color = "ethnicity_label",
        weight = "wt.subcohort",
        panel_titles = labels.title2[tp, ] %>% unlist(),
        axis_titles = labels.axis[tp, ] %>% unlist(),
        arrange_nrow = ceiling(length(assay_immuno) / 3),
        arrange_ncol = 3,
        filename = paste0(
          save.results.to,
          "/demographics/Marker_Rcdf_",
          tp, "_trt_",
          c("placebo_", "vaccine_")[trt],
          bstatus.labels.2[bstatus],
          "_by_ethnicity_", study_name,
          ".pdf"
        )
      )

      # (7) race
      covid_corr_boxplot_facets(
        plot_dat =
          subset(subdat, !(race == "White" &
            WhiteNonHispanic == 0)),
        x = "race",
        y = tp,
        facet_by = "assay",
        ylim = assay_lim[rep(assay_immuno, ifelse(length(assay_immuno)==1, 2, 1)), tp, ], # call the same marker twice if only one marker exists
        plot_LLOX = !grepl("Delta", tp),
        POS.CUTOFFS = log10(pos.cutoffs[assay_immuno]),
        LLOX = log10(lloxs[assay_immuno]),
        ULOQ = log10(uloqs[assay_immuno]),
        axis_titles_y = labels.axis[tp, ] %>% unlist(),
        panel_titles = labels.title2[tp, ] %>% unlist(),
        arrange_nrow = ceiling(length(assay_immuno) / 3),
        arrange_ncol = 3,
        filename = paste0(
          save.results.to,
          "/demographics/boxplots_",
          tp, "_",
          bstatus.labels.2[bstatus],
          "_trt_", trt.labels[trt],
          "_by_race_", study_name,
          ".pdf"
        )
      )

      covid_corr_rcdf_facets(
        plot_dat =
          subset(subdat, !(race == "White" &
            WhiteNonHispanic == 0)),
        x = tp,
        facet_by = "assay",
        xlim = assay_lim[rep(assay_immuno, ifelse(length(assay_immuno)==1, 2, 1)), tp, ], # call the same marker twice if only one marker exists
        color = "race",
        weight = "wt.subcohort",
        panel_titles = labels.title2[tp, ] %>% unlist(),
        axis_titles = labels.axis[tp, ] %>% unlist(),
        arrange_nrow = ceiling(length(assay_immuno) / 3),
        arrange_ncol = 3,
        filename = paste0(
          save.results.to,
          "/demographics/Marker_Rcdf_",
          tp, "_trt_",
          c("placebo_", "vaccine_")[trt],
          bstatus.labels.2[bstatus],
          "_by_race_", study_name,
          ".pdf"
        )
      )

      if(study_name=="ENSEMBLE" | study_name=="MockENSEMBLE") {
        minority.data <- subset(subdat, Country==0)
      } else {
        minority.data <- subdat
      }
      
      ##  (8) minority status
      if(!attr(config,"config") %in% c("janssen_la_partA","janssen_sa_partA",
                                       "janssen_la_partAsenior","janssen_la_partAnonsenior",
                                       "janssen_sa_partAnonsenior")){
        covid_corr_boxplot_facets(
          plot_dat = minority.data,
          x = "minority_label",
          y = tp,
          facet_by = "assay",
          ylim = assay_lim[rep(assay_immuno, ifelse(length(assay_immuno)==1, 2, 1)), tp, ], # call the same marker twice if only one marker exists
          plot_LLOX = !grepl("Delta", tp),
          POS.CUTOFFS = log10(pos.cutoffs[assay_immuno]),
          LLOX = log10(lloxs[assay_immuno]),
          ULOQ = log10(uloqs[assay_immuno]),
          axis_titles_y = labels.axis[tp, ] %>% unlist(),
          panel_titles = labels.title2[tp, ] %>% unlist(),
          arrange_nrow = ceiling(length(assay_immuno) / 3),
          arrange_ncol = 3,
          filename = paste0(
            save.results.to,
            "/demographics/boxplots_",
            tp, "_",
            bstatus.labels.2[bstatus],
            "_trt_", trt.labels[trt],
            "_by_minority_group_",
            study_name, ".pdf"
          )
        )
  
        covid_corr_rcdf_facets(
          plot_dat = minority.data,
          x = tp,
          facet_by = "assay",
          xlim = assay_lim[rep(assay_immuno, ifelse(length(assay_immuno)==1, 2, 1)), tp, ], # call the same marker twice if only one marker exists
          color = "minority_label",
          weight = "wt.subcohort",
          panel_titles = labels.title2[tp, ] %>% unlist(),
          axis_titles = labels.axis[tp, ] %>% unlist(),
          arrange_nrow = ceiling(length(assay_immuno) / 3),
          arrange_ncol = 3,
          filename = paste0(
            save.results.to,
            "/demographics/Marker_Rcdf_",
            tp, "_trt_",
            c("placebo_", "vaccine_")[trt],
            bstatus.labels.2[bstatus],
            "_by_minority_group_",
            study_name, ".pdf"
          )
        )
  
        ##  (9) age * minority status
        covid_corr_boxplot_facets(
          plot_dat = minority.data,
          x = "age_minority_label",
          y = tp,
          facet_by = "assay",
          ylim = assay_lim[rep(assay_immuno, ifelse(length(assay_immuno)==1, 2, 1)), tp, ], # call the same marker twice if only one marker exists
          plot_LLOX = !grepl("Delta", tp),
          POS.CUTOFFS = log10(pos.cutoffs[assay_immuno]),
          LLOX = log10(lloxs[assay_immuno]),
          ULOQ = log10(uloqs[assay_immuno]),
          axis_titles_y = labels.axis[tp, ] %>% unlist(),
          panel_titles = labels.title2[tp, ] %>% unlist(),
          arrange_nrow = ceiling(length(assay_immuno) / 3),
          arrange_ncol = 3,
          filename = paste0(
            save.results.to,
            "/demographics/boxplots_",
            tp, "_",
            bstatus.labels.2[bstatus],
            "_trt_", trt.labels[trt],
            "_by_age_x_minority_",
            study_name, ".pdf"
          )
        )
  
        covid_corr_rcdf_facets(
          plot_dat = minority.data,
          x = tp,
          facet_by = "assay",
          xlim = assay_lim[rep(assay_immuno, ifelse(length(assay_immuno)==1, 2, 1)), tp, ], # call the same marker twice if only one marker exists
          color = "age_minority_label",
          weight = "wt.subcohort",
          panel_titles = labels.title2[tp, ] %>% unlist(),
          axis_titles = labels.axis[tp, ] %>% unlist(),
          arrange_nrow = ceiling(length(assay_immuno) / 3),
          arrange_ncol = 3,
          filename = paste0(
            save.results.to,
            "/demographics/Marker_Rcdf_",
            tp, "_trt_",
            c("placebo_", "vaccine_")[trt],
            bstatus.labels.2[bstatus],
            "_by_age_minority_group_",
            study_name, ".pdf"
          )
        )
      }
      
      if(study_name=="ENSEMBLE" | study_name=="MockENSEMBLE") {
        
        ##  (10) country
        covid_corr_boxplot_facets(
          plot_dat = subdat,
          x = "country_label",
          y = tp,
          facet_by = "assay",
          ylim = assay_lim[assay_immuno, tp, ],
          plot_LLOX = !grepl("Delta", tp),
          POS.CUTOFFS = log10(pos.cutoffs[assay_immuno]),
          LLOX = log10(lloxs[assay_immuno]),
          ULOQ = log10(uloqs[assay_immuno]),
          axis_titles_y = labels.axis[tp, ] %>% unlist(),
          panel_titles = labels.title2[tp, ] %>% unlist(),
          arrange_nrow = ceiling(length(assay_immuno) / 3),
          arrange_ncol = 3,
          filename = paste0(
            save.results.to,
            "/demographics/boxplots_",
            tp, "_",
            bstatus.labels.2[bstatus],
            "_trt_", trt.labels[trt],
            "_by_country_",
            study_name, ".pdf"
          )
        )
        
        covid_corr_rcdf_facets(
          plot_dat = subdat,
          x = tp,
          facet_by = "assay",
          xlim = assay_lim[assay_immuno, tp, ],
          color = "country_label",
          weight = "wt.subcohort",
          panel_titles = labels.title2[tp, ] %>% unlist(),
          axis_titles = labels.axis[tp, ] %>% unlist(),
          arrange_nrow = ceiling(length(assay_immuno) / 3),
          arrange_ncol = 3,
          filename = paste0(
            save.results.to,
            "/demographics/Marker_Rcdf_",
            tp, "_trt_",
            c("placebo_", "vaccine_")[trt],
            bstatus.labels.2[bstatus],
            "_by_country_",
            study_name, ".pdf"
          )
        )
        
        ##  (10) HIV positivity
        covid_corr_boxplot_facets(
          plot_dat = subdat,
          x = "hiv_label",
          y = tp,
          facet_by = "assay",
          ylim = assay_lim[assay_immuno, tp, ],
          plot_LLOX = !grepl("Delta", tp),
          POS.CUTOFFS = log10(pos.cutoffs[assay_immuno]),
          LLOX = log10(lloxs[assay_immuno]),
          ULOQ = log10(uloqs[assay_immuno]),
          axis_titles_y = labels.axis[tp, ] %>% unlist(),
          panel_titles = labels.title2[tp, ] %>% unlist(),
          arrange_nrow = ceiling(length(assay_immuno) / 3),
          arrange_ncol = 3,
          filename = paste0(
            save.results.to,
            "/demographics/boxplots_",
            tp, "_",
            bstatus.labels.2[bstatus],
            "_trt_", trt.labels[trt],
            "_by_hiv_group_",
            study_name, ".pdf"
          )
        )
        
        covid_corr_rcdf_facets(
          plot_dat = subdat,
          x = tp,
          facet_by = "assay",
          xlim = assay_lim[assay_immuno, tp, ],
          color = "hiv_label",
          weight = "wt.subcohort",
          panel_titles = labels.title2[tp, ] %>% unlist(),
          axis_titles = labels.axis[tp, ] %>% unlist(),
          arrange_nrow = ceiling(length(assay_immuno) / 3),
          arrange_ncol = 3,
          filename = paste0(
            save.results.to,
            "/demographics/Marker_Rcdf_",
            tp, "_trt_",
            c("placebo_", "vaccine_")[trt],
            bstatus.labels.2[bstatus],
            "_by_hiv_group_",
            study_name, ".pdf"
          )
        )
      }

      
      print(paste0("Done with loop of ", tp, ", trt=",
                   trt," and bstatus=", bstatus))
    }
  }
}