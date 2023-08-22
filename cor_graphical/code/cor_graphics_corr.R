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

## load data 
dat.cor.data.pair <- readRDS(here("data_clean", "cor_data_pair.rds")); dat.cor.data.pair$all_one <- 1 # as a placeholder for strata values

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
if (attr(config,"config") == "janssen_partA_VL" & COR == "D29VLvariant") {

  for (t in "Day29"){
    for (trt in c(1)){
      assay_metadata_sub_sa <- subset(assay_metadata, assay %in% c("pseudoneutid50", "pseudoneutid50_Delta",
                                                                "pseudoneutid50_Beta"))
      dat.cor.data.pair.SA <- subset(dat.cor.data.pair, Region == 2 & Trt==1)
      dat.cor.data.pair.SA$Day29pseudoneutid50_Delta = sample(dat.cor.data.pair.SA$Day29pseudoneutid50, dim(dat.cor.data.pair.SA)[1])
      dat.cor.data.pair.SA$Day29pseudoneutid50_Beta = sample(dat.cor.data.pair.SA$Day29pseudoneutid50, dim(dat.cor.data.pair.SA)[1])
        
      covid_corr_pairplots(
        plot_dat = dat.cor.data.pair.SA,
        time = t,
        assays = assay_metadata_sub_sa$assay,
        strata = "all_one",
        weight = "wt.D29",
        plot_title = paste0(
          "Correlations of 3 ", t, " antibody markers in Southern Africa,\nCorr = Weighted Spearman Rank Correlation."
        ),
        column_labels = paste(t, assay_metadata_sub_sa$assay_label_short),
        height = max(1.3 * length(assay_metadata_sub_sa$assay) + 0.1, 5.5),
        width = max(1.3 * length(assay_metadata_sub_sa$assay), 5.5),
        column_label_size = ifelse(max(nchar(paste(t, assay_metadata_sub_sa$assay_label_short)))>40, 3.8, 4.3),
        filename = paste0(
          save.results.to, "/pairs_by_time_", t,
          "_markers_",ifelse(trt==1, "vaccine", "placebo"), "_NAb_SA.pdf"
        )
      )
      
      assay_metadata_sub_la <- subset(assay_metadata, assay %in% c("pseudoneutid50", "pseudoneutid50_Zeta",
                                                                "pseudoneutid50_Mu", "pseudoneutid50_Gamma",
                                                                "pseudoneutid50_Lambda"))
      
      dat.cor.data.pair.LA <- subset(dat.cor.data.pair, Region == 1 & Trt==1)
        
      covid_corr_pairplots(
        plot_dat = dat.cor.data.pair.LA,
        time = t,
        assays = assay_metadata_sub_la$assay,
        strata = "all_one",
        weight = "wt.D29",
        plot_title = paste0(
          "Correlations of 5 ", t, " antibody markers in Latin America,\nCorr = Weighted Spearman Rank Correlation."
        ),
        column_labels = paste(t, assay_metadata_sub_la$assay_label_short),
        height = max(1.3 * length(assay_metadata_sub_la$assay) + 0.1, 5.5),
        width = max(1.3 * length(assay_metadata_sub_la$assay), 5.5),
        column_label_size = ifelse(max(nchar(paste(t, assay_metadata_sub_la$assay_label_short)))>40, 3.3, 4.3),
        filename = paste0(
          save.results.to, "/pairs_by_time_", t,
          "_markers_", ifelse(trt==1, "vaccine", "placebo"), "_NAb_LA.pdf"
        )
      )
    
      }
  }
}
