#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
library(GGally)
library(stringr)
require(devtools)
install_version("dummies", version = "1.5.6", repos = "http://cran.us.r-project.org")
library(grid)
library(gridExtra)
install.packages("wCorr", repos = "http://cran.us.r-project.org") # for the weightedCorr() in pairplot, weighted correlation
library(wCorr)
    
# There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))
    
Args <- commandArgs(trailingOnly=TRUE)
COR=Args[1]
if (grepl("IncludeNotMolecConfirmed", COR)) {incNotMol <- "IncludeNotMolecConfirmed"
} else {incNotMol <- ""}
#-----------------------------------------------

source(here::here("code", "covid_corr_plot_functions.R"))
source(here::here("..", "_common.R"))
# load parameters
source(here::here("code", "params.R"))

## load data 
dat.cor.data.pair <- readRDS(here::here("data_clean", "cor_data.rds")); dat.cor.data.pair$all_one <- 1 # as a placeholder for strata values
config.cor <- config::get(config = COR)

# path for figures and tables etc
save.results.to = here::here("output")
if (!dir.exists(save.results.to))  dir.create(save.results.to)
save.results.to = paste0(here::here("output"), "/", attr(config,"config"));
if (!dir.exists(save.results.to))  dir.create(save.results.to)
save.results.to = paste0(save.results.to, "/", COR,"/");
if (!dir.exists(save.results.to))  dir.create(save.results.to)
print(paste0("save.results.to equals ", save.results.to))


###### Correlation plots across markers at a given time point
# 3 markers (Anc, Delta, Beta), SA, Day 29
if (study_name == "janssen_partA_VL" & COR == "D29VLvariant") {

  for (t in "Day29"){
    for (trt in c(1)){
      assay_metadata_sub_sa <- subset(assay_metadata, assay %in% c("pseudoneutid50", "pseudoneutid50_Delta",
                                                                "pseudoneutid50_Beta"))
      dat.cor.data.pair.SA <- subset(dat.cor.data.pair, Region == 2 & Trt==1) 
      # 111 rows, 13 Post-Peak Cases, 11 of 13 is Beta Cases (Post-Peak Cases-Beta), 2 no variant call (Post-Peak Cases), these two doesn't have pseudoneutid50_Delta, pseudoneutid50_Beta values 
      #           98 Non-Cases, all with pseudoneutid50, pseudoneutid50_Delta, pseudoneutid50_Beta
      
      covid_corr_pairplots(
        plot_dat = dat.cor.data.pair.SA,
        time = t,
        assays = assay_metadata_sub_sa$assay,
        strata = "all_one",
        weight = config.cor$wt,
        plot_title = paste0(
          "Correlations of 3 ", t, " Antibody Markers in South Africa,\nCorr = Weighted Spearman Rank Correlation."
        ),
        column_labels = paste(t, assay_metadata_sub_sa$assay_label_short),
        height = max(1.3 * length(assay_metadata_sub_sa$assay) + 0.1, 5.5),
        width = max(1.3 * length(assay_metadata_sub_sa$assay), 5.5),
        column_label_size = ifelse(max(nchar(paste(t, assay_metadata_sub_sa$assay_label_short)))>40, 4.2, 4.3),
        filename = paste0(
          save.results.to, "/pairs_by_time_", t,
          "_", ifelse(trt==1, "vaccine", "placebo"), "_bseroneg_NAb_SA.pdf"
        )
      )
      
      assay_metadata_sub_la <- subset(assay_metadata, assay %in% c("pseudoneutid50", "pseudoneutid50_Zeta",
                                                                "pseudoneutid50_Mu", "pseudoneutid50_Gamma",
                                                                "pseudoneutid50_Lambda"))
      
      dat.cor.data.pair.LA <- subset(dat.cor.data.pair, Region == 1 & Trt==1) 
      # 406 rows, 260 Post-Peak Cases, 53 Reference, 71 Gamma, 43 Lambda, 37 Mu, 18 Zeta, 38 no variant call 
      #           146 Non-Cases, all with pseudoneutid50, pseudoneutid50_Gamma, pseudoneutid50_Zeta, pseudoneutid50_Mu, pseudoneutid50_Lambda
      # (just use 146 non-cases to draw the plot)
      
      covid_corr_pairplots(
        plot_dat = dat.cor.data.pair.LA,
        time = t,
        assays = assay_metadata_sub_la$assay,
        strata = "all_one",
        weight = config.cor$wt,
        plot_title = paste0(
          "Correlations of 5 ", t, " Antibody Markers in Latin America,\nCorr = Weighted Spearman Rank Correlation."
        ),
        column_labels = paste(t, assay_metadata_sub_la$assay_label_short),
        height = max(1.3 * length(assay_metadata_sub_la$assay) + 0.1, 5.5),
        width = max(1.3 * length(assay_metadata_sub_la$assay), 5.5),
        column_label_size = ifelse(max(nchar(paste(t, assay_metadata_sub_la$assay_label_short)))>40, 3.45, 4.3),
        filename = paste0(
          save.results.to, "/pairs_by_time_", t,
          "_", ifelse(trt==1, "vaccine", "placebo"), "_bseroneg_NAb_LA.pdf"
        )
      )
    
      }
  }
} else if (study_name=="IARCHPV"){
  
  for (t in "M18"){
    
    for (asy in c("bind", "pseudoneutid50", "some_bind", "some_pseudoneutid50")) {
      
      if (!grepl("pseudoneutid50", assays) && grepl("pseudoneutid50", asy)) next
      if (!grepl("bind", assays) && grepl("bind", asy)) next
      
      # all markers but the marker score
      if(asy %in% c("bind", "pseudoneutid50")) {
        assay_metadata_sub = subset(assay_metadata, grepl(asy, assay) & !grepl("mdw", assay))
      } else if (asy=="some_bind") {
        assay_metadata_sub = subset(assay_metadata, assay %in% c("bind_HPV6","bind_HPV11","bind_HPV16","bind_HPV18","bind_HPV31","bind_mdw"))
      } else if (asy=="some_pseudoneutid50") {
        assay_metadata_sub = subset(assay_metadata, assay %in% c("pseudoneutid50_HPV6","pseudoneutid50_HPV11","pseudoneutid50_HPV16",
                                                                 "pseudoneutid50_HPV18","pseudoneutid50_HPV31","pseudoneutid50_mdw"))
      }
      
      trt_lb = ""
      
      covid_corr_pairplots(
        plot_dat = dat.cor.data.pair,
        time = t,
        assays = assay_metadata_sub$assay,
        strata = "all_one",
        weight = config.cor$wt,
        plot_title = paste0(
          "Pairwise Correlations of ", paste0(t, if(COR=="M18sus") "sus"), " Antibody Markers\nCorr = Weighted Spearman Rank Correlation."
        ),
        column_labels = paste(t, assay_metadata_sub$assay_label_short),
        height = max(1.3 * length(assay_metadata_sub$assay) + 0.1, 5.5),
        width = max(1.3 * length(assay_metadata_sub$assay), 5.5),
        column_label_size = ifelse(max(nchar(paste(t, assay_metadata_sub$assay_label_short)))>40, 4.2, 4.3),
        filename = paste0(
          save.results.to, "/pairs_by_time_", paste0(t, if(COR=="M18sus") "sus"), # COR: M18, M18sus
          "_pooled_", gsub("bind", "BAb", gsub("pseudoneutid50", "NAb", asy)), ".pdf"
        )
      )
    }
  }
  
}
