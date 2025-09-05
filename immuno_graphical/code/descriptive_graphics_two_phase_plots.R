#Sys.setenv(TRIAL = "vat08_combined")
#Sys.setenv(TRIAL = "nextgen_mock")
#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
# There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))
Sys.setenv(DESCRIPTIVE = 1)
source(here::here("..", "_common.R"))
if (!is.null(config$assay_metadata)) {pos.cutoffs = assay_metadata$pos.cutoff}
#-----------------------------------------------

library(here)
library(tidyverse)
library(stringr)
library(ggplot2)
library(ggpubr)
library(GGally)
library(spatstat.geom)
library(scales)
library(grid) # textGrob
library(gridExtra)
library(PResiduals)
library(fmsb) # radarchart()
library(wCorr) # weighted correlation


source(here("code", "params.R"))
if (attr(config,"config")=="vat08_combined"){
  timepoints_ = timepoints
  times_ = times
}

assay_lim <- readRDS(here("data_clean", "assay_lim.rds"))
if (study_name %in% c("VAT08", "ENSEMBLE", "NextGen_Mock") | attr(config,"config")=="prevent19_stage2"){
  source(here("code", "covid_corr_plot_functions.R"))
  source(here("code", "process_violin_pair_functions.R")) # pair plot functions in this program are overwritten by those in the second program
  # pairplots are non-stratum-adjusted, no resampling, IPS-weighted spearman correlation
} #else {
#  source(here("code", "ggally_cor_resample.R"))
#  source(here("code", "covid_corr_plot_functions.R")) # pair plot functions in this program produces geom_statistics w/ resampling-based covariate-adjusted Spearman
#}

set.seed(12345)
# load cleaned data
dat.long.twophase.sample <- readRDS(here(
  "data_clean",
  "long_twophase_data.rds"
))

dat.twophase.sample <- readRDS(here::here("data_clean", "twophase_data.rds")); dat.twophase.sample$all_one <- 1 # as a placeholder for strata values
dat.spider <- readRDS(here::here("data_clean", "twophase_data.rds"))

tps_no_delta_over_tinterm <-  times_[!times_ %in% c(paste0("Delta",timepoints_[length(timepoints_)],"over",timepoints_[1]))] #c("B", "Day29", "Delta29overB", "Day57", "Delta57overB")
tps_no_B_and_delta_over_tinterm <-  times_[!times_ %in% c("B",paste0("Delta",timepoints_[length(timepoints_)],"over",timepoints_[1]))] #c("Day29", "Delta29overB", "Day57", "Delta57overB")
tps_no_fold_change <- times_[!grepl("Delta", times_)]
tps_no_B_and_fold_change <- times_[!grepl("Delta", times_) & times_!="B"]
tps_delta_over_B <- times_[grepl("overB",times_)]

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
# - Pairs plots/scatterplots and simple spearman rank correlations (or resampling-based strata-adjusted partial spearman rank correlation) are used.
#-----------------------------------------------
for (country in if(attr(config,"config")=="prevent19") {c("Nvx_US_Mex","Nvx_US")} else if (attr(config,"config")=="janssen_partA_VL") {c(1,2)} else {c("all")}) { 
  # this loop is for prevent19 and janssen_partA_VL, prevent19 needs to be looped through all people (US+Mex) and US only, janssen_partA_VL needs to be looped through regions
    
  if (length(assay_immuno)==1) next # AZ two datasets only have one marker in each as of 5/13/2022, can't do pair 
  
    print("Pair plots 1:")
  
    country_lb = case_when(country=="Nvx_US" ~ "US_only_", 
                           country==0 ~ "NAM_",
                           country==1 ~ "LATAM_",
                           country==2 ~ "ZA_",
                           TRUE ~ "")
    country_lb_long = case_when(country=="Nvx_US" ~ "US only", 
                                country==0 ~ "Northern America",
                                country==1 ~ "Latin America",
                                country==2 ~ "Southern Africa",
                                TRUE ~ "")
    
    for (tp in if(study_name=="VAT08") {tps_no_fold_change} else if (attr(config,"config")=="janssen_partA_VL") {paste0("Day", timepoints_)} else {tps_no_B_and_delta_over_tinterm}) { # "B", "Day29", "Day57", "Day29overB", "Day57overB"
      if (study_name == "NextGen_Mock") next # NextGen doesn't need pairplot at single timepoint
      
      for (trt in 0:1) {
        
        # Don't produce figures for placebo baseline negative except for VAT08 to improve build time
        if(trt==0 & study_name!="VAT08" & attr(config,"config")!="janssen_partA_VL") {bstatus.range <- 1} else {bstatus.range <- unique(dat.twophase.sample$Bserostatus)}
    
        for (bserostatus in bstatus.range) {
          if (!bserostatus %in% unique(dat.twophase.sample$Bserostatus)) next
          if (attr(config,"config") %in% c("janssen_partA_VL","prevent19_stage2") & (bserostatus==1|trt==0)) next # skip baseline positive and placebo for janssen_partA_VL, prevent19_stage2
          
          tt=match(tp, times_)
          
          subdat_ <- dat.twophase.sample %>%
            dplyr::filter(Bserostatus == bserostatus & Trt == trt)
          
          if(attr(config,"config")=="prevent19" & country=="Nvx_US") {subdat=subset(subdat, Country==0)} # Nvx Country: (USA = 0, MEX  = 1)
          if(attr(config,"config")=="janssen_partA_VL") {subdat=subset(subdat, Region==country)} # janssen_partA_VL Region: (Northern America = 0, Latin America = 1, Southern Africa = 2)
          
          for (asy in c("*", 
                        if (attr(config,"config") %in% c("janssen_partA_VL","vat08_combined")) "bind", 
                        if (attr(config,"config") %in% c("janssen_partA_VL","vat08_combined")) "pseudo",
                        
                        ## janssen_partA_VL needs some Spike-ID50 paired results
                        
                        # for Latin America, Reference|Mu|Gamma|Lambda
                        if (attr(config,"config")=="janssen_partA_VL") "Reference",#"bindSpike$|pseudoneutid50$", 
                        if (attr(config,"config")=="janssen_partA_VL") "Mu",#"bindSpike_B.1.621$|pseudoneutid50_Mu$",
                        if (attr(config,"config")=="janssen_partA_VL") "Gamma",#"bindSpike_P.1$|pseudoneutid50_Gamma$",
                        if (attr(config,"config")=="janssen_partA_VL") "Lambda",#"bindSpike_C.37$|pseudoneutid50_Lambda$",
                        
                        # for South Africa, Index/Reference, Delta-Score/Delta, Beta
                        #if (attr(config,"config")=="janssen_partA_VL") "bindSpike$|pseudoneutid50$",
                        if (attr(config,"config")=="janssen_partA_VL") "Delta",#"bindSpike_DeltaMDW$|pseudoneutid50_Delta$",
                        if (attr(config,"config")=="janssen_partA_VL") "Beta"#"bindSpike_B.1.351$|pseudoneutid50_Beta$"
                        
                        )){
            # this loop is only for janssen_partA_VL becasue it needs to be looped through all, bAb and nAb
            assay_lb = case_when(asy=="*" ~ "",
                                 asy=="bind" ~ "bAb_",
                                 asy=="pseudo" ~ "nAb_",
                                 TRUE ~ paste0(asy, "_"))
            
            #if (study_name=="VAT08" && bserostatus==0 && tp=="B") { # psv_mdw doesn't have any value for naive at baseline
              #assay_immuno_ = assay_immuno[assay_immuno!="pseudoneutid50_mdw"]
              #} else 
            if (asy=="bind") {assay_immuno_ = subset(assay_immuno, grepl(asy, assay_immuno))
              } else if (asy=="pseudo" & country == "all") {assay_immuno_ = subset(assay_immuno, grepl(asy, assay_immuno))
              } else if (asy %in% c("*", "pseudo") & country == 1) { # Latin America = 1
                selected_latam = assay_immuno[grepl("Reference|Zeta|Mu|Gamma|Lambda", assay_metadata$assay_label_short)]
                assay_immuno_ = subset(selected_latam, grepl(asy, selected_latam))
              } else if (asy %in% c("*", "pseudo") & country == 2) { # Southern Africa = 2
                selected_za = assay_immuno[grepl("Reference|Beta|Delta", assay_metadata$assay_label_short)]
                assay_immuno_ = subset(selected_za, grepl(asy, selected_za))
              } else if (asy %in% c("Reference","Mu","Gamma","Lambda","Delta","Beta")) {
                assay_immuno_ = assay_immuno[grepl(asy, assay_metadata$assay_label_short)]
              } else {assay_immuno_ = assay_immuno}
            
            if (sum(complete.cases(subdat[, paste0(tp, assay_immuno_)]))==0) next # skip if no assay data available
            
            if (asy %in% c("Mu","Gamma","Lambda") & country == 2) next # don't need these figures for janssen_partA_VL
            if (asy %in% c("Delta","Beta") & country == 1) next # don't need these figures for janssen_partA_VL
            
            # subset for prevent19_stage2
            if(attr(config,"config")=="prevent19_stage2") {
              subdat = subdat_
              subdat = subdat[which(subdat[, paste0("ph2.immuno.",gsub("ay","",tp))]==1), ]
            } else if (attr(config,"config")=="vat08_combined" & asy %in% c("bind", "*")) {
              subdat = subdat_ %>% filter(ph2.immuno.bAb == 1)
            } else if (attr(config,"config")=="vat08_combined" & asy=="pseudo") {
              subdat = subdat_ %>% filter(ph2.immuno.nAb == 1)
            } else (subdat = subdat_)
            
            covid_corr_pairplots(
              plot_dat = subdat,
              time = tp,
              assays = assay_immuno_, # adhoc request by David: assay_immuno = c("bindSpike", "bindSpike_P.1", "bindRBD", "bindRBD_P.1", "bindN")
                                     # adhoc request 2 by David: assay_immuno = c("liveneutmn50", "bindSpike_P.1", "bindRBD_P.1", "bindN")
              strata = "all_one",
              weight = ifelse(attr(config,"config")=="prevent19_stage2" & tp=="Day35", "wt.immuno.D35",
                              ifelse(attr(config,"config")=="prevent19_stage2" & tp=="C1", "wt.immuno.C1",
                                     ifelse(attr(config,"config")=="prevent19_stage2" & tp=="BD1", "wt.immuno.BD1",
                                          ifelse(attr(config,"config")=="vat08_combined" & asy!="pseudo", "wt.immuno.bAb",
                                                ifelse(attr(config,"config")=="vat08_combined" & asy=="pseudo", "wt.immuno.nAb",
                              "wt.subcohort"))))),
              plot_title = paste0(
                gsub("ay ","", labels.time)[tt],
                ifelse(assay_lb=="*"," Ab", paste0(" ", gsub("_", "", assay_lb))), " markers: ",
                bstatus.labels.3[bserostatus + 1], ", ",
                c("placebo", "vaccine")[trt + 1], " arm", if (attr(config,"config")=="janssen_partA_VL") paste0(", ", country_lb_long)
              ),
              column_labels = labels.axis[tp, assay_immuno_] %>% unlist(), # adhoc request by David: labels.axis[tp, match(assay_immuno, colnames(labels.axis))]
              height = max(1.3 * length(assay_immuno_) + 0.1, 5.5),
              width = max(1.3 * length(assay_immuno_), 5.5),
              column_label_size = ifelse(length(assay_immuno_) <= 2, 8,
                                         ifelse(max(str_length(labels.axis[1,])) > 28, 4, 6.5)),
              filename = paste0(
                save.results.to, "/pairs_", tp,
                "_Markers_", bstatus.labels.2[bserostatus + 1],
                c("_placebo_arm", "_vaccine_arm")[trt + 1], "_", country_lb, assay_lb,
                study_name, ".pdf"
              )
            )
          }
        }
      }
    }
  
    
  ## pairplots of assay readouts for multiple timepoints
  print("Pair plots 2:")
  if (study_name != "ENSEMBLE"){
    for (trt in 0:1) {
      # Don't produce figures for placebo baseline negative for studies other than VAT08 to improve build time
      if(trt==0 & study_name!="VAT08") {bstatus.range <- 1} else {bstatus.range <- unique(dat.twophase.sample$Bserostatus)}
    
      for (bserostatus in bstatus.range) {
        if (!bserostatus %in% unique(dat.twophase.sample$Bserostatus)) next
        
        subdat_ <- dat.twophase.sample %>%
          dplyr::filter(Bserostatus == bserostatus & Trt == trt)
        
        times_selected <- if(study_name=="VAT08") {list(tps_no_delta_over_tinterm[c(1,4,5)])
          # "B", "Day22", "Day43", "Day22overB", "Day43overB", only show B and fold_change for Sanofi study
          } else if (study_name == "NextGen_Mock") {list(
              tps_no_delta_over_tinterm[c(1,2,6,3,7)]#, # "B", "Day31", "Day181", "Delta31overB", "Delta181overB"
              #tps_no_delta_over_tinterm[c()] # 
              ) 
          } else {list(tps_no_fold_change)} # "B", "Day29", "Day57"
        
        if(attr(config,"config")=="prevent19" & country=="Nvx_US") {subdat = subset(subdat_, Country==0)} # Nvx Country: (USA = 0, MEX  = 1)
        
        for (tm in seq(length(times_selected))){
            for (aa in assay_immuno) {
              #if (study_name=="VAT08" && aa=="pseudoneutid50_mdw" && bserostatus==0) next # psv_mdw doesn't have any value for naive at baseline
              
              # subset for prevent19_stage2
              if(attr(config,"config")=="prevent19_stage2") {
                subdat = subdat_ %>% filter(ph2.immuno.BD1 == 1)
              } else if (attr(config,"config")=="vat08_combined" & grepl("bind", aa)) {
                subdat = subdat_ %>% filter(ph2.immuno.bAb == 1) 
              } else if (attr(config,"config")=="vat08_combined" & grepl("pseudo", aa)) {
                subdat = subdat_ %>% filter(ph2.immuno.nAb == 1)
              } else {subdat = subdat_}
              
              covid_corr_pairplots_by_time(
                plot_dat = subdat,
                times = times_selected[[tm]],
                assay = aa,
                strata = "all_one",
                weight = ifelse(attr(config,"config")=="prevent19_stage2", "wt.immuno.BD1",
                                ifelse(attr(config,"config")=="vat08_combined" & grepl("bind", aa), "wt.immuno.bAb",
                                       ifelse(attr(config,"config")=="vat08_combined" & grepl("pseudo", aa), "wt.immuno.nAb",
                                              ifelse(attr(config,"config")=="nextgen_mock" & grepl("bind|pseudo", aa), "wt.immuno",
                                                     ifelse(attr(config,"config")=="nextgen_mock" & grepl("T4|T8", aa), "wt.AB.immuno",
                                "wt.subcohort"))))),
                plot_title = paste0(ifelse(study_name == "NextGen_Mock" & trt == 0, "(B) ",
                                           ifelse(study_name == "NextGen_Mock" & trt == 1, "(A) ", "")),
                  labels.assays[aa], ": ",
                  if (study_name == "NextGen_Mock") {""} else {paste0(bstatus.labels.3[bserostatus + 1], " ")},
                  tolower(trt.labels)[trt + 1], " arm"
                ), 
                plot_title_size = ifelse(study_name == "NextGen_Mock", 7, 10), 
                column_labels = paste(gsub("ay ","", labels.time[times_selected[[tm]]]),
                                      "\n", labels.axis[, aa][1]),
                column_label_size = ifelse(study_name=="VAT08", 4.5, 
                                           max(5.4, min(8, 8 - (str_length(labels.axis[1, aa]) - 14) * (8 - 5.4) / (34 - 14))) * 0.59
                                           ),
                axis_label_size = ifelse(grepl("T4|T8", aa), 4.7, ifelse(study_name == "NextGen_Mock", 6, ifelse(study_name=="VAT08", 7, 9))),
                label_format = ifelse(grepl("T4|T8", aa), "percent", "log10"),
                filename = paste0(
                  save.results.to, "/pairs_", aa, "_by_times_", ifelse(tm!=1, paste0(tm, "_"), ""), 
                  bstatus.labels.2[bserostatus + 1], "_", paste0(gsub(" ", "_", tolower(trt.labels)), "_")[trt + 1], country_lb,
                  study_name, ifelse(study_name == "NextGen_Mock", "_final", ""), ".pdf"
                )
              )
            }
          }
      }
    }
  }
  
  print("Pair plots 3:")
  if (study_name=="VAT08") { # request only for this study, at day 1, pooling over vaccine and placebo
    tp = "B"
    for (bserostatus in bstatus.range) {
      
        #if (bserostatus==0) { # VAT08 psv_mdw doesn't have any value for naive at baseline
        #  assay_immuno_ = assay_immuno[assay_immuno!="pseudoneutid50_mdw"]
        #} else {
      
      for (asy in c("*", "bind", "pseudo")){
        assay_immuno_ = if (asy=="*") {assay_immuno} else {assay_immuno[grepl(asy, assay_immuno)]}
        
        assay_lb = case_when(asy=="*" ~ "",
                             asy=="bind" ~ "bAb_",
                             asy=="pseudo" ~ "nAb_")
      
        if (!bserostatus %in% unique(dat.twophase.sample$Bserostatus)) next
        
        tt=match(tp, times_)
        
        subdat_ <- dat.twophase.sample %>%
          dplyr::filter(Bserostatus == bserostatus)
        
        if (attr(config,"config")=="vat08_combined" & asy %in% c("bind", "*")) {
          subdat = subdat_ %>% filter(ph2.immuno.bAb == 1)
        } else if (attr(config,"config")=="vat08_combined" & asy=="pseudo") {
          subdat = subdat_ %>% filter(ph2.immuno.nAb == 1)
        } else {subdat = subdat_}
        
        covid_corr_pairplots(
          plot_dat = subdat,
          time = tp,
          assays = assay_immuno_,
          strata = "all_one",
          weight = ifelse(attr(config,"config")=="vat08_combined" & grepl("bind", aa), "wt.immuno.bAb",
                                 ifelse(attr(config,"config")=="vat08_combined" & grepl("pseudo", aa), "wt.immuno.nAb")),
          plot_title = paste0(
            gsub("ay ","", labels.time)[tt],
            " Ab markers: ",
            bstatus.labels.3[bserostatus + 1], ", pooled arm"
          ),
          column_labels = labels.axis[tp, assay_immuno_] %>% unlist(),
          height = max(1.3 * length(assay_immuno_) + 0.1, 5.5),
          width = max(1.3 * length(assay_immuno_), 5.5),
          column_label_size = ifelse(max(str_length(labels.axis[1, aa])) > 28, 4, 6.5),
          filename = paste0(
            save.results.to, "/pairs_", tp,
            "_Markers_", bstatus.labels.2[bserostatus + 1],
            "_pooled_arm_", country_lb, assay_lb, 
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
for (tp in if(!study_name %in% c("VAT08")) {tps_no_B_and_delta_over_tinterm} else if (study_name == "NextGen_Mock") {
  tps_no_fold_change # "B", "Day31", "Day91", "Day181", "Day366"
} else {tps_no_fold_change}) { # "Day29", "Day57", "Day29overB", "Day57overB" for most studies; if VAT08, "Day1", "Day22", "Day43"
  
  if (attr(config,"config") %in% c("janssen_partA_VL", "nextgen_mock")) next # NextGen_Mock and janssen_partA_VL don't need these plots
  
  # subset for prevent19_stage2
  if(attr(config,"config")=="prevent19_stage2") {
    subdat_rcdf1 = dat.long.twophase.sample[which(dat.long.twophase.sample[, paste0("ph2.immuno.",gsub("ay","",tp))]==1), ]
  } else if (attr(config,"config")=="vat08_combined") {
    subdat_rcdf1 = dat.long.twophase.sample[which(dat.long.twophase.sample$ph2.immuno.bAb==1), ]
  } else if (study_name == "NextGen_Mock" & tp %in% c("B", "Day31", "Day181")) {
    subdat_rcdf1 = dat.long.twophase.sample
    subdat_rcdf1$wt = with(subdat_rcdf1, ifelse(grepl("T4|T8", assay), wt.AB.immuno, wt.immuno)) # ICS assay use wt.AB.immuno as weight for whole RIS/RIS-PBMC
  } else if (study_name == "NextGen_Mock" & tp %in% c("Day91", "Day366")) {
    subdat_rcdf1 = dat.long.twophase.sample[which(dat.long.twophase.sample$Track == "A"), ]
    subdat_rcdf1$wt = subdat_rcdf1$wt.AB.immuno
  } else {subdat_rcdf1 = dat.long.twophase.sample}
  
  categories = c(paste0(trt.labels[1], ", ", bstatus.labels[1]), 
                 paste0(trt.labels[1], ", ", bstatus.labels[2]), 
                 paste0(trt.labels[2], ", ", bstatus.labels[1]),
                 paste0(trt.labels[2], ", ", bstatus.labels[2]))
  colors = c("#1749FF", "#D92321", "#0AB7C9", "#FF6F1B")
  
  covid_corr_rcdf_facets(
    plot_dat = subdat_rcdf1,
    x = tp,
    facet_by = "assay",
    color = "trt_bstatus_label",
    palette = setNames(colors, categories),
    legend = setNames(categories, categories),
    weight = ifelse(attr(config,"config")=="prevent19_stage2" & tp=="Day35", "wt.immuno.D35",
                    ifelse(attr(config,"config")=="prevent19_stage2" & tp=="C1", "wt.immuno.C1",
                           ifelse(attr(config,"config")=="prevent19_stage2" & tp=="BD1", "wt.immuno.BD1",
                                  ifelse(attr(config,"config")=="vat08_combined", "wt.immuno.bAb",
                                         ifelse(attr(config,"config")=="nextgen_mock", "wt",
                                                "wt.subcohort"))))),
    xlim = assay_lim[rep(assay_immuno, ifelse(length(assay_immuno)==1, 2, 1)), tp, ], # call the same marker twice if only one marker exists
    arrange_ncol = ifelse(study_name %in% c("VAT08", "NextGen_Mock"), 4, 3),
    arrange_nrow = ceiling(length(assay_immuno) / 3),
    panel_titles = paste0(ifelse(study_name == "NextGen_Mock" & trt == 0, "(B) ",
                                 ifelse(study_name == "NextGen_Mock" & trt == 1, "(A) ", "")), 
                          labels.title2[tp, ] %>% unlist()),
    axis_titles = labels.axis[tp, ] %>% unlist(),
    xbreaks = ifelse(study_name=="VAT08", 2, 1),
    axis_title_size = ifelse(study_name == "NextGen_Mock", 8, 10),
    axis_size = ifelse(study_name == "NextGen_Mock", 8, 10),
    panel_title_size = ifelse(study_name == "VAT08", 8, ifelse(study_name == "NextGen_Mock", 7, 10)),
    height = ifelse(attr(config,"config")=="prevent19_stage2", 10, 
                    ifelse(study_name=="VAT08", 3 * ceiling(length(assay_immuno) / 4) + 0.5,
                           ifelse(study_name=="NextGen_Mock", 3 * ceiling(length(assay_immuno) / 4) + 0.5,
                                 3 * ceiling(length(assay_immuno) / 3) + 0.5))),
    filename = paste0(
      save.results.to, "/Marker_Rcdf_", tp,
      "_trt_both_bstatus_both_", study_name, 
      ifelse(study_name == "NextGen_Mock" & tp %in% c("B", "Day31", "Day181"), "_final", 
             ifelse(study_name == "NextGen_Mock" & tp %in% c("Day91", "Day366"), "_initial", "")), ".pdf"
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

# plot bAb, PsV and ADCP assays separately
for (Ab in c(if (study_name != "NextGen_Mock") "bind", 
             if (study_name != "NextGen_Mock") "pseudo", 
             "bind.*IgG_sera", "bind.*IgG_nasal", "bind.*IgG_saliva", 
             "bind.*IgA_sera", "bind.*IgA_nasal", "bind.*IgA_saliva",
             "pseudo.*sera", "pseudo.*nasal", "pseudo.*saliva", 
             "ADCP", "T4|T8"#,
             #"pseudoneutid50_sera_XBB.1.5"
             )) {
  
  Ab_lb = case_when(Ab=="ADCP" ~ "other_",
                    Ab=="bind" ~ "bAb_",
                    Ab=="pseudo" ~ "nAb_", 
                    
                    Ab=="bind.*IgG_sera" ~ "bAb_IgG_sera_",
                    Ab=="bind.*IgG_nasal" ~ "bAb_IgG_nasal_",
                    Ab=="bind.*IgG_saliva" ~ "bAb_IgG_saliva_",
                    
                    Ab=="bind.*IgA_sera" ~ "bAb_IgA_sera_",
                    Ab=="bind.*IgA_nasal" ~ "bAb_IgA_nasal_",
                    Ab=="bind.*IgA_saliva" ~ "bAb_IgA_saliva_",
                    
                    Ab=="pseudo.*sera" ~ "nAb_sera_",
                    Ab=="pseudo.*nasal" ~ "nAb_nasal_",
                    Ab=="pseudo.*saliva" ~ "nAb_saliva_",
                    
                    Ab=="T4|T8" ~ "ics_"#,
                    
                    #Ab=="pseudoneutid50_sera_XBB.1.5" ~ "pseudoneutid50_sera_XBB.1.5_"
                    )
  
  rcdf_assays <- assay_immuno[grepl(Ab, assay_immuno)]
  
  if (length(rcdf_assays) == 0) next
    
  #-----------------------------------------------
  # RCDF plot 
  # one treatment arm and two baseline status per plot
  #-----------------------------------------------
  print("RCDF 2:")
  for (tp in if(!study_name %in% c("VAT08", "NextGen_Mock")) {tps_no_B_and_delta_over_tinterm
    } else if (study_name == "NextGen_Mock") {
      c(paste0(tps_no_fold_change, "_initial"), tps_no_fold_change[c(1,2,4)]) # "B"      "Day31"  "Day91"  "Day181" "Day366"
  } else {tps_no_B_and_fold_change}) { # "Day29", "Day57", "Day29overB", "Day57overB" for most studies; if VAT08, "Day22", "Day43"
      
    if (attr(config,"config") %in% c("janssen_partA_VL", "prevent19_stage2", "nextgen_mock")) next # janssen_partA_VL, prevent19_stage2 doesn't need these plots
    
    for (trt in c(trt.labels[2], if(study_name %in% c("VAT08", "NextGen_Mock")) trt.labels[1])){
      
      subdat_rcdf2_ = subset(dat.long.twophase.sample, Trt == trt & assay %in% rcdf_assays)
      
      if (attr(config,"config")=="vat08_combined" & Ab=="bind") {
        subdat_rcdf2 = subdat_rcdf2_ %>% filter(ph2.immuno.bAb==1)
      } else if (attr(config,"config")=="vat08_combined" & Ab=="pseudo") {
        subdat_rcdf2 = subdat_rcdf2_ %>% filter(ph2.immuno.nAb==1)
      } else if (study_name == "NextGen_Mock" & tp %in% c("B", "Day31", "Day181")) {
        subdat_rcdf2 = subdat_rcdf2_ %>%
          mutate(wt = ifelse(grepl("T4|T8", assay), wt.AB.immuno, wt.immuno))  # ICS assay use wt.AB.immuno as weight for whole RIS/RIS-PBMC
      } else if (study_name == "NextGen_Mock" & grepl("_initial", tp)) {
        subdat_rcdf2 = subdat_rcdf2_ %>% 
          filter(Track == "A") %>%
          mutate(wt = wt.AB.immuno)
      }
      
      covid_corr_rcdf(
        plot_dat = subdat_rcdf2,
        x = gsub("_initial", "", tp),
        color = "assay_labels",
        lty = if (study_name == "NextGen_Mock") NULL else "Bserostatus",
        weight = ifelse(attr(config,"config")=="vat08_combined" & Ab=="bind", "wt.immuno.bAb",
                        ifelse(attr(config,"config")=="vat08_combined" & Ab=="pseudo", "wt.immuno.nAb",
                               ifelse(attr(config,"config")=="nextgen_mock", "wt",
                                      "wt.subcohort"))),
        xlab = if (Ab == "T4|T8") {"Percent of T cells expressing indicated function"
        } else if (grepl("bind", Ab)) {"Concentration of binding antibodies (AU/ml)"
        } else if (grepl("pseudo", Ab)) {"nAb ID50 titer (AU/ml)"} else {paste0(gsub("ay ", "", labels.time[gsub("_initial", "", tp)]), " Ab Markers")},
        xlim = c(min(assay_lim[rcdf_assays, gsub("_initial", "", tp), 1]), 
                 max(assay_lim[rcdf_assays, gsub("_initial", "", tp), 2])),
        xbreaks = seq(floor(min(assay_lim[rcdf_assays, gsub("_initial", "", tp), 1])), 
                      ceiling(max(assay_lim[rcdf_assays, gsub("_initial", "", tp), 2])), 
                      ifelse(study_name=="VAT08", 3, 1)),
        plot_title = paste0(ifelse(study_name == "NextGen_Mock" & trt == 0, "(B) ",
                                   ifelse(study_name == "NextGen_Mock" & trt == 1, "(A) ", "")), 
                            paste0(labels.time[gsub("_initial", "", tp)], " Ab Markers")),
        legend_size = ifelse(length(rcdf_assays) > 15, 5, ifelse(length(rcdf_assays) >= 4, 8, 14)), 
        axis_size = ifelse(attr(config,"config")=="nextgen_mock", 10, 16), 
        label_format = ifelse(Ab == "T4|T8", "percent", "log10"),
        legend_nrow = ifelse(length(rcdf_assays) < 15, length(rcdf_assays), ceiling(length(rcdf_assays)/2)),
        filename = paste0(
          save.results.to, "/Marker_Rcdf_", Ab_lb, 
          gsub("_initial", "", tp),
          "_trt_", tolower(gsub(" ", "_", trt)), "_bstatus_both_", study_name, 
          ifelse(study_name == "NextGen_Mock" & tp %in% c("B", "Day31", "Day181"), "_final", 
                 ifelse(study_name == "NextGen_Mock" & grepl("_initial", tp), "_initial", "")), ".pdf"
        )
      )
    }
  }
    
  #-----------------------------------------------
  # RCDF plot 
  # one treatment arm and one baseline status per plot, in vaccine/placebo arm
  #-----------------------------------------------
  print("RCDF 3:")
  for (bstatus in 1:2) {
    if (attr(config,"config") %in% c("vat08_combined","prevent19_stage2","nextgen_mock")) next # vat08_combined, prevent19_stage2 doesn't need these plots
    
    if (nrow(subset(dat.long.twophase.sample, Bserostatus==bstatus.labels[bstatus]))==0) next
    if (attr(config,"config")=="janssen_partA_VL" && bstatus==2) next # do not plot baseline positive for janssen_partA_VL
    
    for (tp in if (attr(config,"config")=="janssen_partA_VL") {paste0("Day", timepoints_)} else {tps_no_B_and_delta_over_tinterm}) { # "Day29", "Day57", "Day29overB", "Day57overB" for most studies; if VAT08, "Day22", "Day43"
      
      for (country in if (attr(config,"config")=="janssen_partA_VL") {c(1,2)} else {"all"}) { # loop through regions for janssen_partA_VL
      
        country_lb = case_when(country==0 ~ "NAM_",
                               country==1 ~ "LATAM_",
                               country==2 ~ "ZA_",
                               TRUE ~ "")
        country_lb_long = case_when(country==0 ~ "Northern America",
                                    country==1 ~ "Latin America",
                                    country==2 ~ "Southern Africa",
                                    TRUE ~ "")
        
        # this loop is just for janssen_partA_VL to further subsetting neut assays for different regions
        if (attr(config,"config")!="janssen_partA_VL") {rcdf_assays_ = rcdf_assays 
        } else if (Ab %in% c("bind")) {
          rcdf_assays_ = rcdf_assays
        } else if (Ab %in% c("pseudo") & country == 1) { # Latin America = 1
          selected = assays[grepl("Reference|Zeta|Mu|Gamma|Lambda", assay_metadata$assay_label_short)]
          rcdf_assays_ = subset(selected, grepl(Ab, selected))
        } else if (Ab %in% c("pseudo") & country == 2) { # Southern Africa = 2
          selected = assays[grepl("Reference|Beta|Delta", assay_metadata$assay_label_short)]
          rcdf_assays_ = subset(selected, grepl(Ab, selected))
        }
        
        
        if(attr(config,"config")=="janssen_partA_VL") {
          subdat_rcdf3 = subset(dat.long.twophase.sample, Region == country & Trt == trt & assay %in% rcdf_assays)
          } else {subdat_rcdf3 = subset(dat.long.twophase.sample, Trt == trt & assay %in% rcdf_assays)
          } # janssen_partA_VL Region: (Northern America = 0, Latin America = 1, Southern Africa = 2)
        
        if (attr(config,"config")=="vat08_combined" & Ab=="bind") {
          subdat_rcdf3 = subdat_rcdf3 %>% filter(ph2.immuno.bAb==1)
        } else if (attr(config,"config")=="vat08_combined" & Ab=="pseudo") {
          subdat_rcdf3 = subdat_rcdf3 %>% filter(ph2.immuno.nAb==1)
        }
        
        covid_corr_rcdf(
          plot_dat = subdat_rcdf3,
          x = tp,
          color = "assay_labels",
          lty = NULL,
          weight = ifelse(attr(config,"config")=="vat08_combined" & Ab=="bind", "wt.immuno.bAb",
                          ifelse(attr(config,"config")=="vat08_combined" & Ab=="pseudo", "wt.immuno.nAb",
                                 "wt.subcohort")),
          xlab = if (Ab == "T4|T8") {"Percent of T cells expressing indicated function"
        } else if (grepl("bind", Ab)) {"Concentration of binding antibodies (AU/ml)"
        } else if (grepl("pseudo", Ab)) {"nAb ID50 titer (AU/ml)"} else {paste0(gsub("ay ", "", labels.time[gsub("_initial", "", tp)]), " Ab Markers")},
        xlim = c(min(assay_lim[rcdf_assays_, tp, 1]), 
                   max(assay_lim[rcdf_assays_, tp, 2])),
          xbreaks = seq(min(assay_lim[rcdf_assays_, tp, 1]), 
                        max(assay_lim[rcdf_assays_, tp, 2]), 
                        1),
          plot_title = paste0(labels.time[tp], " Ab Markers: ", bstatus.labels.3[bstatus], ", vaccine arm", if (country %in% c(0, 1, 2)) paste0("\n", country_lb_long)
          ),
          filename = paste0(
            save.results.to, "/Marker_Rcdf_", Ab_lb, tp,
            "_trt_vaccine_bstatus_", c("Neg", "Pos")[bstatus], "_", country_lb, study_name, ".pdf"
          )
        )
      }
    }
  }
    
  #-----------------------------------------------
  # RCDF plot 
  # two treatment arms, one baseline status per plot
  #-----------------------------------------------
  if (study_name %in% c("VAT08"#, "NextGen_Mock"
                        )){
    print("RCDF 4:")
    for (bstatus in 1:2) {
      if (nrow(subset(dat.long.twophase.sample, Bserostatus==bstatus.labels[bstatus]))==0) next
      
      for (tp in if (study_name == "NextGen_Mock") {
        c(paste0(tps_no_fold_change, "_initial"), tps_no_fold_change[c(1,2,4)]) # "B"      "Day31"  "Day91"  "Day181" "Day366"
      } else {tps_no_B_and_fold_change}) { # "Day29", "Day57", "Day29overB", "Day57overB" for most studies; if VAT08, "Day22", "Day43"
        
        subdat_rcdf4_ = subset(dat.long.twophase.sample, Bserostatus == bstatus.labels[bstatus] & assay %in% rcdf_assays)
        
        if (attr(config,"config")=="vat08_combined" & Ab=="bind") {subdat_rcdf4 = subdat_rcdf4_ %>% filter(ph2.immuno.bAb == 1)
        } else if (attr(config,"config")=="vat08_combined" & Ab=="pseudo") {subdat_rcdf4 = subdat_rcdf4_ %>% filter(ph2.immuno.nAb == 1)
        } else if (study_name == "NextGen_Mock" & tp %in% c("B", "Day31", "Day181")) {
          subdat_rcdf4 = subdat_rcdf4_ %>%
            mutate(wt = ifelse(grepl("T4|T8", assay), wt.AB.immuno, wt.immuno))  # ICS assay use wt.AB.immuno as weight for whole RIS/RIS-PBMC
        } else if (study_name == "NextGen_Mock" & grepl("_initial", tp)) {
          subdat_rcdf4 = subdat_rcdf4_ %>% 
            filter(Track == "A") %>%
            mutate(wt = wt.AB.immuno)
        } else {subdat_rcdf4 = subdat_rcdf4_}
        
        covid_corr_rcdf(
          plot_dat = subdat_rcdf4,
          x = gsub("_initial", "", tp),
          color = "assay_labels",
          lty = "Trt",
          weight =  ifelse(attr(config,"config")=="vat08_combined" & Ab=="bind", "wt.immuno.bAb",
                           ifelse(attr(config,"config")=="vat08_combined" & Ab=="pseudo", "wt.immuno.nAb", 
                                  ifelse(attr(config,"config")=="nextgen_mock", "wt", ""))),
          xlab = if (Ab == "T4|T8") {"Percent of T cells expressing indicated function"
          } else if (grepl("bind", Ab)) {"Concentration of binding antibodies (AU/ml)"
          } else if (grepl("pseudo", Ab)) {"nAb ID50 titer (AU/ml)"} else {paste0(gsub("ay ", "", labels.time[gsub("_initial", "", tp)]), " Ab Markers")},
          xlim = c(min(assay_lim[rcdf_assays, gsub("_initial", "", tp), 1]), 
                   max(assay_lim[rcdf_assays, gsub("_initial", "", tp), 2])),
          xbreaks = seq(floor(min(assay_lim[rcdf_assays, gsub("_initial", "", tp), 1])), 
                        ceiling(max(assay_lim[rcdf_assays, gsub("_initial", "", tp), 2])), 
                        ifelse(study_name=="VAT08", 3, 1)),
          plot_title = paste0(ifelse(study_name == "NextGen_Mock" & trt == 0, "(B) ",
                                     ifelse(study_name == "NextGen_Mock" & trt == 1, "(A) ", "")), 
                              paste0(labels.time[gsub("_initial", "", tp)], " Ab Markers")),
          legend_size = ifelse(length(rcdf_assays) > 15, 5, ifelse(length(rcdf_assays) >= 4, 8, 14)), 
          axis_size = ifelse(attr(config,"config")=="nextgen_mock", 10, 16), 
          label_format = ifelse(Ab == "T4|T8", "percent", "log10"),
          legend_nrow = ifelse(length(rcdf_assays) < 15, length(rcdf_assays), ceiling(length(rcdf_assays)/2)),
          filename = paste0(
            save.results.to, "/Marker_Rcdf_", Ab_lb, gsub("_initial", "", tp),
            "_trt_both_bstatus_", c("Neg", "Pos")[bstatus], "_", study_name, 
            ifelse(study_name == "NextGen_Mock" & tp %in% c("B", "Day31", "Day181"), "_final", 
                   ifelse(study_name == "NextGen_Mock" & grepl("_initial", tp), "_initial", "")), ".pdf"
          )
        )
      }
    }
  }
}

#-----------------------------------------------
# RCDF plot 
# two treatment arms, one baseline status side by side, all timepoints per plot
#-----------------------------------------------

if (study_name %in% c("NextGen_Mock")){
  print("RCDF 5:")
  
  for (Ab in assays) {

  Ab_lb = paste0(Ab, "_")

    for (bstatus in 2) {
      
      for (tp in c("_initial", "_final")){
        
        for (trt in trt.labels){
      
        subdat_rcdf5_ = subset(dat.long.twophase.sample, Bserostatus == bstatus.labels[bstatus] & assay %in% Ab & trt == trt)
        
        if (study_name == "NextGen_Mock" & grepl("_final", tp)) {
          subdat_rcdf5 = subdat_rcdf5_ %>%
            mutate(wt = ifelse(grepl("T4|T8", Ab), wt.AB.immuno, wt.immuno))  # ICS assay use wt.AB.immuno as weight for whole RIS/RIS-PBMC
        } else if (study_name == "NextGen_Mock" & grepl("_initial", tp)) {
          subdat_rcdf5 = subdat_rcdf5_ %>% 
            filter(Track == "A") %>%
            mutate(wt = wt.AB.immuno)
        } else {subdat_rcdf5 = subdat_rcdf5_}
        
        if (tp == "_initial") {
          subdat_rcdf5_long = subdat_rcdf5 %>% 
            pivot_longer(
              cols = c("B", "Day31", "Day91", "Day181", "Day366"),
              names_to = "time",
              values_to = "value"
            )
        } else if (tp == "_final"){
          subdat_rcdf5_long = subdat_rcdf5 %>% 
            pivot_longer(
              cols = c("B", "Day31", "Day181"),
              names_to = "time",
              values_to = "value"
            )
        }
        
        subdat_rcdf5_long$time_labels <-
          factor(subdat_rcdf5_long$time,
                 levels = tps_no_fold_change,
                 labels = labels.time[tps_no_fold_change])
        
        covid_corr_rcdf_facet_adhoc(
          plot_dat = subdat_rcdf5_long,
          x = "value", # at each time
          color = "time_labels",
          lty = NULL,
          facet_by = "Trt",
          weight =  "wt",
          xlab = if (grepl("T4|T8", Ab)) {"Percent of T cells expressing indicated function"
          } else if (grepl("bind", Ab)) {"Concentration of binding antibodies (AU/ml)"
          } else if (grepl("pseudo", Ab)) {"nAb ID50 titer (AU/ml)"},
          xlim = c(min(assay_lim[Ab, , 1]), 
                   max(assay_lim[Ab, , 2])),
          xbreaks = 1,
          panel_titles = c(paste0("(A) ", trt.labels[2]), paste0("(B) ", trt.labels[1])),
          legend_size = 8,
          legend = setNames(trt.labels, trt.labels), 
          axis_size = 5, 
          overall_title = labels.assays.short[Ab],
          label_format = ifelse(grepl("T4|T8", Ab), "percent", "log10"),
          legend_nrow = 1,
          arrange_ncol = 2,
          arrange_nrow = 1,
          width = 5,
          filename = paste0(
            save.results.to, "/Marker_Rcdf_", Ab_lb,
            "_trt_both_bstatus_", c("Neg", "Pos")[bstatus], "_", study_name, tp, ".pdf"
          )
        )
        } # end of trt
      } # end of initial vs final
  } # end of bsero
  } # end of assay
}

#-----------------------------------------------
# BOX PLOTS
#-----------------------------------------------
# - Box plots across treatment groups.
# - made a ggplot object for every assay and use the ggarrange() function to combine the resulted plots.
# - boxplots of assay readouts at Day 1, Day tinterm and Day tpeak, 
# - by treatment groups or by baseline serostatus
#-----------------------------------------------

#-----------------------------------------------
# box plot 
# two treatment arms, one baseline status per plot
#-----------------------------------------------
print("Boxplots 1:")
for (bstatus in 1:2) {
  
  if (attr(config,"config") %in% c("janssen_partA_VL")) next # janssen_partA_VL doesn't need these plots
  
  if (nrow(subset(dat.long.twophase.sample, Bserostatus==bstatus.labels[bstatus]))==0) next 
  
  for (tp in if(!study_name %in% c("VAT08", "NextGen_Mock")) {tps_no_B_and_delta_over_tinterm} else {tps_no_fold_change}) { # "Day29", "Day57", "Day29overB", "Day57overB" for most studies; if VAT08, "Day1", "Day22", "Day43"
    
    for (pn in if (study_name != "NextGen_Mock") {""} else {c("IgG_sera", "IgA_sera", "pseudoneutid50_sera", 
                                                              "IgG_nasal", "IgA_nasal", "pseudoneutid50_nasal", 
                                                              "IgG_saliva", "IgA_saliva", "pseudoneutid50_saliva", 
                                                              "T4_T8")}) {
      if (study_name == "NextGen_Mock" & any(grepl(pn, assays)) == FALSE) next
      
      # subset for prevent19_stage2
      if(attr(config,"config")=="prevent19_stage2") {
        subdat_box1_ = dat.long.twophase.sample[which(dat.long.twophase.sample[, paste0("ph2.immuno.", gsub("ay","",tp))]==1), ]
      } else if (attr(config,"config")=="vat08_combined") {
        subdat_box1_ = dat.long.twophase.sample %>% filter((grepl("bind", assay) & ph2.immuno.bAb==1) | (grepl("pseudo", assay) & ph2.immuno.nAb==1) )
      } else if (attr(config,"config") == "nextgen_mock" & tp %in% c("Day91", "Day366")){
        subdat_box1_ = dat.long.twophase.sample %>% filter(Track == "A")
      } else {subdat_box1_ = dat.long.twophase.sample}
      
      if (pn == "") {subdat_box1 = subdat_box1_
      } else {subdat_box1 = subdat_box1_ %>% filter(grepl(gsub("T4_T8", "T4|T8", pn), assay)) %>% mutate(assay = droplevels(assay))}
      
      # reorder per Peter's request
      if (study_name == "NextGen_Mock" & grepl("IgG", pn)) {
        subdat_box1$assay <- factor(subdat_box1$assay,
                                     levels = levels(subdat_box1$assay)[c(2:11, 1)])
      }
      
      if (study_name == "NextGen_Mock" & tp == "B") {subdat_box1$Trt = "Pooled Arm"}
      
      assay_sub = levels(subdat_box1$assay)
        
      covid_corr_boxplot_facets(
        plot_dat = subset(subdat_box1,
          Bserostatus == bstatus.labels[bstatus]
        ),
        x = "Trt",
        y = tp,
        color = "Trt",
        facet_by = "assay",
        palette = if (study_name == "NextGen_Mock" & tp == "B") {c("#FF6F1B")
          } else if (study_name == "NextGen_Mock" & tp != "B") {c("#1749FF", "#378252")
          } else {c(
          "#1749FF", "#D92321",
          "#0AB7C9", "#FF6F1B",
          "#810094", "#378252",
          "#FF5EBF", "#3700A5",
          "#8F8F8F", "#787873"
        )},
        ylim = assay_lim[assay_sub, tp, ],
        plot_LLOX = !grepl("Delta", tp), # "B", "Day29", "Day57"
        POS.CUTOFFS = log10(pos.cutoffs[assay_sub]),
        LLOX = log10(lloxs[assay_sub]),
        ULOQ = log10(uloqs[assay_sub]),
        arrange_ncol = ifelse(study_name == "VAT08", 4, ifelse(study_name == "NextGen_Mock" & grepl("Ig", pn), 4, 2)),
        arrange_nrow = ifelse(study_name %in% c("VAT08"), 4, ifelse(study_name == "NextGen_Mock" & grepl("Ig", pn), 3, ifelse(study_name == "NextGen_Mock" & grepl("bindN", pn), 1, 2))),
        #legend = setNames(trt.labels, trt.labels),
        axis_titles_y = labels.axis[tp, assay_sub] %>% unlist(),
        label_format = ifelse(all(grepl("T4|T8", assay_sub)==1), "percent", "log10"),
        panel_titles = labels.title2[tp, assay_sub] %>% unlist(),
        panel_title_size = ifelse(study_name=="VAT08", 8, ifelse(study_name == "NextGen_Mock", 6, 10)),
        height = ifelse(study_name %in% c("VAT08"), 11, 
                        ifelse(attr(config,"config")=="prevent19_stage2", 10, 
                               ifelse(study_name == "NextGen_Mock" & grepl("Ig", pn), 10.5,
                                      ifelse(study_name == "NextGen_Mock" & grepl("bindN", pn), 3.5,
                                             6.8)))),
        add_violin = ifelse(study_name == "NextGen_Mock", T, F),
        filename = paste0(
          save.results.to, "/boxplots_", tp, "_x_trt_", 
          ifelse(study_name == "NextGen_Mock", paste0(pn, "_"), ""),
          bstatus.labels.2[bstatus], "_", study_name, 
          ifelse(study_name == "NextGen_Mock" & tp %in% c("B", "Day31", "Day181"), "_final", 
                 ifelse(study_name == "NextGen_Mock" & tp %in% c("Day91", "Day366"), "_initial", "")), ".pdf"
        )
      )
    }
  }
}

# box plot 
# one treatment arm, two baseline status per plot
#-----------------------------------------------
print("Boxplots 2:")
for (trt in 1:2) {
  
  if (attr(config,"config") %in% c("janssen_partA_VL", "prevent19_stage2", "nextgen_mock")) next # janssen_partA_VL, prevent19_stage2 doesn't need these plots
  
  for (tp in if (study_name!="VAT08") {tps_no_delta_over_tinterm} else {tps_no_fold_change}) {
    
    if (attr(config,"config")=="vat08_combined") {
      subdat_box2 = dat.long.twophase.sample %>% filter((grepl("bind", assay) & ph2.immuno.bAb==1) | (grepl("pseudo", assay) & ph2.immuno.nAb==1) )
    } else {subdat_box2 = dat.long.twophase.sample}
    
    covid_corr_boxplot_facets(
      plot_dat = subset(subdat_box2, as.numeric(Trt) == trt),
      x = "Bserostatus",
      y = tp,
      color = "Bserostatus",
      facet_by = "assay",
      ylim = assay_lim[, tp,],
      plot_LLOX = !grepl("Delta", tp), # "B", "Day29", "Day57"
      POS.CUTOFFS = log10(pos.cutoffs[assay_immuno]),
      LLOX = log10(lloxs[assay_immuno]),
      ULOQ = log10(uloqs[assay_immuno]),
      arrange_ncol = ifelse(study_name=="VAT08", 4, 3),
      arrange_nrow = ifelse(study_name=="VAT08", 4, ceiling(length(assay_immuno) / 3)),
      legend = stringr::str_to_title(bstatus.labels.3),
      axis_titles_y = labels.axis[tp, ] %>% unlist(),
      panel_titles = labels.title2[tp, ] %>% unlist(),
      panel_title_size = ifelse(study_name=="VAT08", 8, 10),
      height = ifelse(study_name=="VAT08", 11, 3 * ceiling(length(assay_immuno) / 3) + 0.5),
      filename = paste0(
        save.results.to, "/boxplots_", tp,
        "_x_bstatus_", c("placebo_arm_", "vaccine_arm_")[trt],
        study_name, ".pdf"
      )
    )
  }
}


#-----------------------------------------------
# box plot 
# one treatment arm, one baseline status per plot
#-----------------------------------------------
print("Boxplots 3:")
if (study_name %in% c("VAT08")) {# this is only reported for VAT08
  for (tp in tps_no_fold_change) {
    
    subdat_box3 = dat.long.twophase.sample %>% filter((grepl("bind", assay) & ph2.immuno.bAb==1) | (grepl("pseudo", assay) & ph2.immuno.nAb==1) )
      
    covid_corr_boxplot_facets(
      plot_dat = subdat_box3 %>% 
        mutate(BseroTrt = factor(paste0(Bserostatus,"\n",Trt),
                                 levels = paste0(rep(bstatus.labels, 2), "\n", rep(trt.labels[2:1], each=2))
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
      arrange_ncol = 4,
      arrange_nrow = ifelse(study_name=="VAT08", 4, ceiling(length(assay_immuno) / 3)),
      legend = paste0(rep(stringr::str_to_title(bstatus.labels.3), 2), ", ", rep(trt.labels[2:1], each=2)),
      axis_titles_y = labels.axis[tp, ] %>% unlist(),
      panel_titles = labels.title2[tp, ] %>% unlist(),
      panel_title_size = ifelse(study_name=="VAT08", 8, 10),
      height = ifelse(study_name=="VAT08", 11, 3 * ceiling(length(assay_immuno) / 3) + 0.5),
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
print("Boxplots 4:")
if (study_name=="VAT08" & F) {# this is only reported for VAT08
  for (tp in tps_no_fold_change) {
    
    subdat_box4 = dat.long.twophase.sample %>% filter(ph2.immuno.bAb==1)
    
    covid_corr_boxplot_facets(
      plot_dat = subdat_box4 %>% mutate(BseroTrtGender = factor(paste0(Bserostatus,"\n",Trt,"\n",sex_label),
                                                                levels = paste0(rep(bstatus.labels, 2), "\n", rep(trt.labels[2:1], each=2) , "\n", rep(c("Female", "Male"), each=4))
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
      arrange_ncol = ifelse(study_name=="VAT08", 4, 3),
      arrange_nrow = ifelse(study_name=="VAT08", 4, ceiling(length(assay_immuno) / 3)),
      legend = paste0(rep(stringr::str_to_title(bstatus.labels.3), 2), ", ", rep(trt.labels[2:1], each=2) , ", ", rep(c("Female", "Male"), each=4)),
      axis_titles_y = labels.axis[tp, ] %>% unlist(),
      panel_titles = labels.title2[tp, ] %>% unlist(),
      panel_title_size = ifelse(study_name=="VAT08", 8, 10),
      height = ifelse(study_name=="VAT08", 11, 3 * ceiling(length(assay_immuno) / 3) + 0.5),
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
if (!attr(config,"config") %in% c("vat08_combined", "prevent19_stage2", "nextgen_mock")){ # no spaghetti plots for VAT08, prevent19_stage2, NextGen_Mock

  print("Spaghetti plots:")
  ## in each baseline serostatus group, randomly select 10 placebo recipients and 20 vaccine recipients
  set.seed(12345)
  
  for (plot in if(attr(config,"config")=="janssen_partA_VL"){c(1,2,3,4)} else {NA}) { # janssen_partA_VL needs four plots but other just needs one
  
    if (attr(config,"config")=="janssen_partA_VL") {
      
      # this is ad-hoc request for janssen_partA_VL
      # at day 29, day 71, month 6, only for vaccine and baseline negative, only for reference assays, only for Latin America 
      times_ = c("Day29","Day71","Mon6") 
      
      assay_immuno_ = if(plot==1) {c("bindSpike")
        } else if(plot==2) {c("pseudoneutid50")
        } else if(plot==3) {assays[grepl("bind", assays)]
            } else if(plot==4) {assays[grepl("pseudo", assays) & !grepl("Delta|Beta", assays)]}
  
      spaghetti_ptid <- dat.twophase.sample %>% filter(Trt==1 & Bserostatus==0 & Region==1) %>% dplyr::select(Ptid) %>% pull()
      
      } else {
      
      times_ = times_
      
      assay_immuno_ = assay_immuno
      
      var_names <- do.call(paste0, expand.grid(times[!grepl("Delta", times)], assay_immuno))
      
      spaghetti_ptid <- dat.twophase.sample[, c("Ptid", "Bserostatus", "Trt", var_names)] %>%
        filter(., complete.cases(.)) %>%
        transmute(BT = paste0(as.character(Bserostatus), as.character(Trt)),
                  Ptid = Ptid) %>%
        split(., .$BT) %>%
        lapply(function(xx) {
          if (xx$BT[1] %in% c("10", "00")) {
            sample(xx$Ptid, 20, ifelse(length(xx$Ptid)<20, T, F))  ## sample 20 placebo recipients
            # add ifelse here because some subset has small sample e.g. janssen_sa_partA
          } else {
            sample(xx$Ptid, 20, ifelse(length(xx$Ptid)<20, T, F))  ## sample 20 vaccine recipients
          }
        }) %>% unlist %>% as.character
    }
    
    spaghetti_dat <- dat.long.twophase.sample[, c("Ptid", "Bserostatus", "Trt", "assay",
                                                  times_[!grepl("Delta",times_)] # "B", "Day29", "Day57"
                                                  )] %>%
      filter(Ptid %in% spaghetti_ptid & assay %in% assay_immuno_) %>%
      pivot_longer(cols = times_[!grepl("Delta",times_)], # "B", "Day29", "Day57"
                   names_to = "time") %>%
      mutate(assay = factor(assay, levels = assay_immuno_, labels = assay_immuno_),
             time_label = factor(time, levels = times_[!grepl("Delta",times_)], # "B", "Day29", "Day57"
                                 labels = gsub("B","D1",gsub("ay","", times_[!grepl("Delta",times_)]))# "D1", "D29", "D57"
                                 )) %>%
      as.data.frame
    
    for (bstatus in 1:2) {
      
      if (attr(config,"config")=="janssen_partA_VL" && bstatus==2) next # skip baseline positive plots for janssen_partA_VL
      
      subdat <- subset(spaghetti_dat, Bserostatus == bstatus.labels[bstatus])
      
      if(nrow(subdat)==0) next
      
      covid_corr_spaghetti_facets(plot_dat = subdat,
                                  x = "time_label",
                                  y = "value",
                                  id = "Ptid",
                                  color = "Trt",
                                  facet_by = "assay",
                                  ylim = assay_lim[, times_[!grepl("B|Delta",times_)][1] # "Day29", "Day57"
                                                   ,],
                                  panel_titles = labels.assays.short[assay_immuno_],
                                  plot_title = ifelse(attr(config,"config")=="janssen_partA_VL", "Baseline Negative Vaccine Group\nLatin America", paste0(
                                    "Baseline ",
                                    c("Negative", "Positive")[bstatus],
                                    " PP Placebo + Vaccine Group"
                                  )),
                                  panel_title_size = ifelse(length(assay_immuno_)==1, 6.5, 8),
                                  arrange_ncol = ifelse(length(assay_immuno_)==1, 1, 3),
                                  arrange_nrow = ceiling(length(assay_immuno_) / 3),
                                  plot_title_size = ifelse(length(assay_immuno_)==1, 8, 12),
                                  axis_size = ifelse(length(assay_immuno_)==1, 8, 12),
                                  axis_title_size = ifelse(length(assay_immuno_)==1, 8, 12),
                                  filename = paste0(
                                    save.results.to, "/spaghetti_plot", ifelse(!is.na(plot), plot, ""), "_",
                                    bstatus.labels.2[bstatus], "_",
                                    study_name, ".pdf"
                                  ))
    }
  }
  
  #-----------------------------------------------
  # Scatter PLOTS
  #-----------------------------------------------
  # - Scatter plots assay vs. age in years, (Day 1) Day tinterm, Day tpeak
  #-----------------------------------------------
  print("Scatter plots:")
  if(study_name!="VAT08" && attr(config,"config")!="janssen_partA_VL"){
    
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
if(attr(config,"config") %in% c("vat08_combined", "janssen_partA_VL", "nextgen_mock")){
  
  # setup pdf file
  for (ab in if(study_name == "NextGen_Mock") {
    c("bind.*sera", "bind.*nasal", "bind.*saliva", "pseudo.*sera", "pseudo.*nasal", "pseudo.*saliva", "T4", "T8")
    } else {c("bind", "pseudo")}) {
    
    for (tm in c(if(attr(config,"config") != "nextgen_mock") "Day", 
                 if(attr(config,"config") == "nextgen_mock") "Day initial", # "B", "Day31", "Day181", "Day91", "Day366" Track A
                 if(attr(config,"config") == "nextgen_mock") "Day whole" # "B", "Day31", "Day181" whole
                 #, if (attr(config,"config")!="janssen_partA_VL") "Delta"
                 )) {
      
      for (bsero in c(if(attr(config,"config") != "nextgen_mock") 0, 1)) {
        
        for (trt_pair in list(c(0, 1))) {
          
          plot_list <- list()
          
          for (trt in trt_pair) {
          
            for (reg in if (attr(config,"config")=="janssen_partA_VL") {c(1,2)} else {"all"}) {
              
              reg_lb = case_when(reg==0 ~ "NAM_",
                                 reg==1 ~ "LATAM_",
                                 reg==2 ~ "ZA_",
                                 TRUE ~ "")
              
              reg_lb_long = case_when(reg==0 ~ "Northern America",
                                      reg==1 ~ "Latin America",
                                      reg==2 ~ "Southern Africa",
                                      TRUE ~ "")
              
              if (attr(config,"config")=="janssen_partA_VL" & (trt==trt.labels[1] | bsero=="Pos")) next
              #if (attr(config,"config") == "nextgen_mock" & tm %in% c("Day whole", "Day initial") & ab == "ics") next # include negative values
              
              # calculate geometric mean of IPS weighted readouts
              times_spider = if (study_name=="VAT08") {times_[c(1,2,3)]
              } else if (tm=="Day") {times_[!grepl("Delta", times_)]
              } else if (tm == "Day initial") {times_[c(1, 2, 4, 6, 8)]
              } else if (tm == "Day whole") {times_[c(1, 2, 6)]
              } else {times_[grepl("Delta", times_)]}
              
              if (!"Region" %in% colnames(dat.spider)) {dat.spider$Region=reg}
              
              assays_ = if (ab %in% c("bind", "bind.*sera", "bind.*nasal", "bind.*saliva",
                                      "pseudo", "pseudo.*sera", "pseudo.*nasal", "pseudo.*saliva")) {assays[grepl(ab, assays) & !grepl("mdw", assays)]
              } else if (ab %in% c("T4","T8")) {assays[grepl(ab, assays) & !grepl("mdw", assays)]
              } else {assays}
              
              if (length(assays_)==0) next
              
              # define cohort and create weight variable
              if (attr(config,"config")=="janssen_partA_VL") {
                dat.spider.by.time_ = dat.spider
                dat.spider.by.time_$wt = dat.spider.by.time_$wt.subcohort
              } else if (attr(config,"config")=="vat08_combined" & ab=="bind") {
                dat.spider.by.time_ = dat.spider %>% filter(ph2.immuno.bAb==1)
                dat.spider.by.time_$wt = dat.spider.by.time_$wt.immuno.bAb
              } else if (attr(config,"config")=="vat08_combined" & ab=="pseudo") {
                dat.spider.by.time_ = dat.spider %>% filter(ph2.immuno.nAb==1)
                dat.spider.by.time_$wt = dat.spider.by.time_$wt.immuno.nAb
              } else if (attr(config,"config")=="nextgen_mock" & tm == "Day initial") {
                dat.spider.by.time_ = dat.spider %>% filter(Track == "A")
                dat.spider.by.time_$wt = dat.spider.by.time_$wt.AB.immuno
              } else if (attr(config,"config")=="nextgen_mock" & tm == "Day whole" & ab != "ics") {
                dat.spider.by.time_ = dat.spider
                dat.spider.by.time_$wt = dat.spider.by.time_$wt.immuno
              } else if (attr(config,"config")=="nextgen_mock" & tm == "Day whole" & ab == "ics") {
                dat.spider.by.time_ = dat.spider
                dat.spider.by.time_$wt = dat.spider.by.time_$wt.AB.immuno
              }
              
              dat.spider.by.time <- dat.spider.by.time_ %>%
                select(one_of("Ptid", "Bserostatus", "Region", "wt", "Trt", 
                              do.call(paste0, expand.grid(times_spider, assays_)))) %>%
                pivot_longer(!Ptid:Trt, names_to = "time_assay", values_to = "value") %>%
                mutate(assay = gsub(paste0(paste0("^",times_spider), collapse="|"), "", time_assay),
                       time = gsub(paste0(assays_, collapse="|"), "", time_assay),
                       time_assay = NULL,
                       value = 10^value
                       ) %>% 
                # in order to make this work for the negative ICS value, use 10^.x to calculate geomean 
                pivot_wider(names_from = assay, values_from = value) %>%
                group_by(time, Bserostatus, Region, Trt) %>%
                summarise(across(all_of(assays_), ~ exp(sum(log(.x * wt), na.rm=T) / sum(wt, na.rm=T)))) %>%
                unique() %>%
                ungroup() %>%
                as.data.frame()
              
              # stack with max and min values
              find_max = round(max(dat.spider.by.time %>% summarise(across(where(is.numeric), ~ max(.x, na.rm = TRUE)))), 1)
              find_min = round(min(dat.spider.by.time %>% summarise(across(where(is.numeric), ~ min(.x, na.rm = TRUE)))), 1)
              
              max_min <- rbind(rep(find_max,
                                   ncol(dat.spider.by.time)), 
                               rep(min(0, floor(find_min)),
                                 #0, 
                                 ncol(dat.spider.by.time)))
              colnames(max_min) <- colnames(dat.spider.by.time)
              rownames(max_min) <- c("max", "min")
              
              dat.spider.by.time <- rbind(max_min, 
                                          dat.spider.by.time)
              
              dat.plot <- dat.spider.by.time[c(1,2), ] %>%
                bind_rows(
                  dat.spider.by.time[2:nrow(dat.spider.by.time), ] %>%
                    filter(grepl(gsub("Day whole", "B|Day31|Day181", gsub("Day initial", "B|Day31|Day181|Day91|Day366", tm)), time) & Bserostatus %in% bsero & Trt %in% trt & Region %in% reg)
                  ) %>%
                mutate(time = NULL, Bserostatus=NULL, Trt=NULL) %>%
                select(if (grepl("bind", ab)) {matches(ab)
                } else if (ab %in% c("T4","T8")) {matches(ab)
                } else if (ab=="pseudo" && reg==1) {matches("pseudoneutid50$|pseudoneutid50_Zeta|pseudoneutid50_Mu|pseudoneutid50_Gamma|pseudoneutid50_Lambda")
                } else if (ab=="pseudo" && reg==2) {matches("pseudoneutid50$|pseudoneutid50_Delta|pseudoneutid50_Beta")
                } else {contains("pseudoneutid50")})
              
              # those without any data will have a weighted geomean equal to 1 because exp(0)=1, set these to NA
              idx <- dat.plot == exp(10^0)#1
              idx[1:2, ] <- FALSE
              dat.plot[idx] <- NA
              
              if (nrow(dat.plot)==2) next
              
              filename = paste0(save.results.to, "/radar_plot_weighted_geomean_", tolower(gsub(" ", "_", tm)), "_", ifelse(reg!="all", reg_lb, ""), gsub("\\.\\*", "_", ab), "_", tolower(bstatus.labels.2[bsero + 1]), "_", "trt_comparison", #trt.labels.2[trt + 1], 
                                ifelse(study_name == "NextGen_Mock" & tm == "Day whole", "_final", 
                                       ifelse(study_name == "NextGen_Mock" & tm == "Day initial", "_initial", "")), ".pdf")
              try(pdf(filename, width = ifelse(study_name == "NextGen_Mock", 15, 5.5), height = 6.5, onefile = TRUE), silent = FALSE)
              cat("Saving to:", filename, "\n")
              par(mfrow=#if (study_name=="VAT08") {c(2,2)} else {
                    c(1, 2) #c(1,1)#}
                  , mar=c(0.1,0.1,1,0.1))
              
              legend_lb = labels.time[times_spider]
              
              spider_range = if(attr(config,"config")=="janssen_partA_VL") {10^seq(1, 1.2, (1.2-1)/4) # hard code for the range here
                #seq(1, 1.2, (1.2-1)/4)} else {seq(0, ceiling(find_max), (ceiling(find_max))/4)}
              } else if (study_name == "NextGen_Mock" & !ab %in% c("T4","T8")) {10^seq(2, 6, 1)#seq(min(ceiling(find_min), 10^0), ceiling(find_max), (ceiling(find_max))/4)
              } else if (study_name == "NextGen_Mock" & ab %in% c("T4","T8")) {10^seq(-2, 2, 1)#10^seq(floor(log10(find_min)), ceiling(log10(find_max)), by = 1)
              }
              
              color = c(if(study_name=="VAT08") "#0AB7C9", "#FF6F1B", "#FF5EBF", "dodgerblue", "chartreuse3", "#009E73")[1:length(times_spider)]
              
              # Save plot data for use in side-by-side plotting
              plot_list[[as.character(trt)]] <- list(
                dat.plot = dat.plot,
                color = color,
                legend_lb = legend_lb,
                title = paste0(ifelse(trt==0, "(B) ", "(A) "), "Geometric Means ", ifelse(grepl("bind", ab), "of bAb Markers, ", ifelse(grepl("pseudo", ab), "of nAb Markers, ", ifelse(grepl("T4|T8", ab), paste0("of CD", gsub("T", "", ab), "+ Markers, "), ""))), 
                               if (study_name == "NextGen_Mock") {""} else {paste0(bstatus.labels.2[bsero + 1], " ")}, trt.labels[trt + 1],
                               ifelse(reg!="all", paste0(", ", reg_lb_long), "")),
                spider_range = spider_range
              )
          } # end of region
          } # end of trt
            
            ############# figure start here
            # Proceed only if both trt=0 and trt=1 plots are available
            if (length(plot_list) == 2) {
              
              for (trt in trt_pair[c(2,1)]) {
                
                p <- plot_list[[as.character(trt)]]
                
                colnames(p$dat.plot) <- assay_metadata$assay_label[match( colnames(dat.plot) , assay_metadata$assay)]
                
                colnames(p$dat.plot) <- gsub("PsV Neutralization to |PsV Neutralization |Binding Antibody to Spike |Binding Antibody to |Binding Antibody |T cells expressing", "", 
                                           gsub("Binding IgG Antibody", "bAb IgG", 
                                                gsub("Binding IgA Antibody", "bAb IgA", 
                                                     gsub("neutralization to", "nAb", colnames(p$dat.plot)))))
                
                radarchart(p$dat.plot, 
                           axistype=1, 
                           # Customize the polygon
                           pcol = scales::alpha(p$color, 0.7), plwd=1.5, pty=c(15), plty=1,
                           pfcol = scales::alpha(p$color, 0.2),
                           #custom the grid
                           cglcol="grey", cglty=1, axislabcol="grey", cglwd=0.8, 
                           caxislabels=if (study_name == "NextGen_Mock" & ab %in% c("T4","T8")) {paste0(p$spider_range, "%")} else {paste0("10^", round(log10(p$spider_range), 2))}, 
                           #label size
                           vlcex=ifelse(study_name=="VAT08", 0.4, ifelse(length(assays_) > 12 | max(nchar(assays_)) > 25, 0.45, 1)),
                           #title
                           title=p$title,
                           #title size
                           cex.main=0.7)
                
                legend("bottomleft", legend=p$legend_lb, lty=5, pch=c(15), col=p$color, bty="n", ncol=3, cex=0.7, inset=c(0.01, 0))
              }
              
              dev.off()
            } # end of plot_list
          } # end of trt_pair
        } # end of bsero
    } # end of tm
  } # end of ab
  
}  




if (F){
  radarchart(dat.plot, 
             axistype=1 , 
             # Customize the polygon
             pcol = scales::alpha(color, 0.7), plwd=1.5, pty=c(15), plty=1,
             pfcol = scales::alpha(color, 0.2),
             #custom the grid
             cglcol="grey", cglty=1, axislabcol="grey", cglwd=0.8, caxislabels=#paste0("10^", spider_range),#
               if (study_name == "NextGen_Mock" & ab == "ics") {paste0(spider_range, "%")} else {paste0("10^", round(log10(spider_range), 2))}, 
             #label size
             vlcex=ifelse(study_name=="VAT08", 0.4, ifelse(length(assays_) > 12 | max(nchar(assays_)) > 25, 0.45, 1)),
             #title
             title=paste0("Geometric Means ", ifelse(grepl("bind", ab), "of bAb Markers, ", ifelse(grepl("pseudo", ab), "of nAb Markers, ", ifelse(grepl("ics", ab), "of CD4+ and CD8+ Markers, ", ""))), 
                          if (study_name == "NextGen_Mock") {""} else {paste0(bstatus.labels.2[bsero + 1], " ")}, trt.labels[trt + 1],
                          ifelse(reg!="all", paste0(", ", reg_lb_long), "")),
             #title size
             cex.main=0.7)
  
  par(xpd=NA)
  
  #legend
  legend("bottomleft", legend=legend_lb, lty=5, pch=c(15),
         col=color, bty="n", ncol=3, cex=0.7,
         inset=c(0.01, 0))
  # print(data.frame(Time = times_spider, Label = legend_lb, Color = color))
}