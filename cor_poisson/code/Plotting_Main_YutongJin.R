library(pheatmap)
library(RColorBrewer)
library(stringr)
library(grid)
library(dplyr)
library(tidyr) 
# for gather/spread
# library(pdf_bookmark)

source("code/Rfunc/figure_correlation.R")
source("code/Rfunc/figure_clustering.R")

config.reporting <- config::get(config = "iliad_ib202p", file="../config.yml") 
dat_202_<-read.csv(config.reporting$data_cleaned,na.strings=c("n/a","NA","N/A","","."))

config.reporting <- config::get(config = "iliad_ib201p", file="../config.yml") 
dat_201_<-read.csv(config.reporting$data_cleaned,na.strings=c("n/a","NA","N/A","","."))


# Trial 201: filter features and markers of interests (339 * 46)
dat_201 <- dat_201_ %>% 
  mutate(
    Trial = "201",
    Gender = ifelse(Sex == 0, "Male", "Female"),
    Trt_Nasal = case_when(Trt == "BPZE1_PBO" ~ 1,
                          Trt == "BPZE1_Tdap" ~ 1,
                          Trt == "PBO_Tdap" ~ 0),
    Trt_tdap = case_when(Trt == "BPZE1_PBO" ~ 0,
                          Trt == "BPZE1_Tdap" ~ 1,
                          Trt == "PBO_Tdap" ~ 1)
    ) %>% 
  select(Ptid, Trial, Gender, 
         Trt, Trt_Nasal, Trt_tdap,
         BFHA_Nasal_IgA: Day28WCE_Serum_IgG)
  

# Trial 202: filter features and markers of interests (45 * 46)
dat_202 <- dat_202_ %>%
  mutate(
    Trial = "202",
    Gender = ifelse(Sex == 0, "Male", "Female"),
    Trt_Nasal = case_when(Trt == "BPZE1" ~ 1,
                          Trt == "PBO" ~ 0),
    Trt_tdap = 0
    ) %>%
  select(Ptid, Trial, Gender,
         Trt, Trt_Nasal, Trt_tdap,
         BFHA_Nasal_IgA: Day28WCE_Serum_IgG)


# output
dir="output/"

################################
#
#  correlation heatmap
#
################################
dat_plot_201 <- dat_201[complete.cases(dat_201), -c(1:6)] # 319 * 40
# 20 participants are excluded due to incomplete data
dat_plot_202 <- dat_202[complete.cases(dat_202), -c(1:6)] # 45 * 40
# no one is exluded due to imcomplete data

pdf(file = paste0(dir, "heatmap_correlation_201_", format(Sys.time(), "%m%d"),".pdf"), height = 12, width = 15)
grid.newpage()
# Baseline markers only
cor_pheatmap(dat_plot = cor(dat_plot_201 %>% select(starts_with("B"))),
             size=12,
             main = "Trial 201 Correlation Heatmap - Baseline Markers")

grid.newpage()
# Peak markers only
cor_pheatmap(dat_plot = cor(dat_plot_201 %>% select(starts_with("Day28"))),
             size=12,
             main = "Trial 201 Correlation Heatmap - Peak Markers")

grid.newpage()
# All markers
cor_pheatmap(dat_plot = cor(dat_plot_201),
             size=12,
             main = "Trial 201 Correlation Heatmap - All Markers")
dev.off()



pdf(file = paste0(dir, "heatmap_correlation_202_", format(Sys.time(), "%m%d"),".pdf"), height = 12, width = 15)
grid.newpage()
# Baseline markers only
cor_pheatmap(dat_plot = cor(dat_plot_202 %>% select(starts_with("B"))),
             size=12,
             main = "Trial 202 Correlation Heatmap - Baseline Markers")

grid.newpage()
# Peak markers only
cor_pheatmap(dat_plot = cor(dat_plot_202 %>% select(starts_with("Day28"))),
             size=12,
             main = "Trial 202 Correlation Heatmap - Peak Markers")

grid.newpage()
# All markers
cor_pheatmap(dat_plot = cor(dat_plot_202),
             size=12,
             main = "Trial 202 Correlation Heatmap - All Markers")
dev.off()







################################
#
#  clustering heatmap
#
################################
# stack dataset from two study
dat_all <- rbind(dat_201, dat_202) # 384 * 46
# subset complete data
dat_all <- dat_all[complete.cases(dat_all),] # 364 * 46
# add rownames using ptid
rownames(dat_all) <- dat_all$Ptid
# remove redundant columns (Ptid, trt)
dat_all <- dat_all %>% select(-Trt, -Ptid) # 364 * 44

# prepare dataset without Nasal IgA
# dat_all_noNasalIgA <- dat_all %>% select(-c(5:14))
dat_all_noNasalIgA <- dat_all %>% select(-(matches("Nasal_IgA$") & !matches("Norm_Nasal_IgA$")))


#--- old version
pdf(file = paste0(dir, "heatmap_clustering_", format(Sys.time(), "%m%d"),".pdf"), height = 12, width = 15)
grid.newpage()
# Baseline markers only
fig_clustering(dat_plot = dat_all_noNasalIgA %>% select(-starts_with("Day28")),
               main = "Clustering Heatmap - Baseline Markers",
               size=12,
               show_assay = TRUE,
               show_antigen = TRUE,
               cluster_method = "euclidean")

grid.newpage()
# Peak markers only
fig_clustering(dat_plot = dat_all_noNasalIgA %>% select(-starts_with("B")),
               main = "Clustering Heatmap - Peak Markers",
               size=12,
               show_assay = TRUE,
               show_antigen = TRUE,
               cluster_method = "euclidean")

grid.newpage()
# All markers
fig_clustering(dat_plot = dat_all_noNasalIgA,
               main = "Clustering Heatmap - All Markers",
               size=12,
               show_assay = TRUE,
               show_antigen = TRUE,
               cluster_method = "euclidean")
dev.off()


#--- new version
pdf(file = paste0(dir, "heatmap_clustering_", format(Sys.time(), "%m%d"),".pdf"), height = 12, width = 15)
grid.newpage()
# Trial 201 Baseline markers only
fig_clustering(dat_plot = dat_all_noNasalIgA %>% select(-starts_with("Day28")) %>% filter(Trial == "201"),
               main = "Trial 201 Clustering Heatmap - Baseline Markers",
               size=12,
               show_assay = TRUE,
               show_antigen = TRUE,
               cluster_method = "euclidean")

grid.newpage()
# Trial 202 Baseline markers only
fig_clustering(dat_plot = dat_all_noNasalIgA %>% select(-starts_with("Day28")) %>% filter(Trial == "202"),
               main = "Trial 202 Clustering Heatmap - Baseline Markers",
               size=12,
               show_assay = TRUE,
               show_antigen = TRUE,
               cluster_method = "euclidean")

grid.newpage()
# Baseline markers only
fig_clustering(dat_plot = dat_all_noNasalIgA %>% select(-starts_with("Day28")),
               main = "Clustering Heatmap - Baseline Markers",
               size=12,
               show_assay = TRUE,
               show_antigen = TRUE,
               cluster_method = "euclidean")

grid.newpage()
# Trial 201 Peak markers only
fig_clustering(dat_plot = dat_all_noNasalIgA %>% select(-starts_with("B")) %>% filter(Trial == "201"),
               main = "Trial 201 Clustering Heatmap - Peak Markers",
               size=12,
               show_assay = TRUE,
               show_antigen = TRUE,
               cluster_method = "euclidean")

grid.newpage()
# Trial 202 Peak markers only
fig_clustering(dat_plot = dat_all_noNasalIgA %>% select(-starts_with("B")) %>% filter(Trial == "202"),
               main = "Trial 202 Clustering Heatmap - Peak Markers",
               size=12,
               show_assay = TRUE,
               show_antigen = TRUE,
               cluster_method = "euclidean")

grid.newpage()
# Peak markers only
fig_clustering(dat_plot = dat_all_noNasalIgA %>% select(-starts_with("B")),
               main = "Clustering Heatmap - Peak Markers",
               size=12,
               show_assay = TRUE,
               show_antigen = TRUE,
               cluster_method = "euclidean")

grid.newpage()
# All markers
fig_clustering(dat_plot = dat_all_noNasalIgA,
               main = "Clustering Heatmap - All Markers",
               size=12,
               show_assay = TRUE,
               show_antigen = TRUE,
               cluster_method = "euclidean")
dev.off()



