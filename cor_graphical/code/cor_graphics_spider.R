#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
library(GGally)
library(stringr)
require(devtools)
library(tidyverse)
install_version("dummies", version = "1.5.6", repos = "http://cran.us.r-project.org")
library(grid)
library(gridExtra)
install.packages("wCorr", repos = "http://cran.us.r-project.org") # for the weightedCorr() in pairplot, weighted correlation
library(wCorr)
install.packages("fmsb", repos = "http://cran.us.r-project.org") # radar plot
library(fmsb) # radarchart()

# There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))

Args <- commandArgs(trailingOnly=TRUE)
COR=Args[1]
if (grepl("IncludeNotMolecConfirmed", COR)) {incNotMol <- "IncludeNotMolecConfirmed"
} else {incNotMol <- ""}
#-----------------------------------------------

source(here::here("..", "_common.R"))

## load data 
dat.cor.data.spider <- readRDS(here::here("data_clean", "cor_data.rds"))

# path for figures and tables etc
save.results.to = here::here("output")
if (!dir.exists(save.results.to))  dir.create(save.results.to)
save.results.to = paste0(here::here("output"), "/", attr(config,"config"));
if (!dir.exists(save.results.to))  dir.create(save.results.to)
save.results.to = paste0(save.results.to, "/", COR,"/");
if (!dir.exists(save.results.to))  dir.create(save.results.to)
print(paste0("save.results.to equals ", save.results.to))


# spider plot showing geometric means calculated using IPS weighting in vaccine arm for 2 region: SA and LA

# calculate geometric mean of IPS weighted readouts
assays_variant <- assays[!assays %in% c("bindSpike", "pseudoneutid50")]
assays_ancestral <- c("bindSpike", "pseudoneutid50")

spider_plot_ <- dat.cor.data.spider %>%
    filter(Region %in% c(0, 1, 2) & Trt==1 & cohort_event2!="Post-Peak Cases") %>%
    mutate(Day29pseudoneutid50 = ifelse(cohort_event2 %in% c("Post-Peak Cases-Reference", "Non-Cases"), Day29pseudoneutid50, NA),
           Day29pseudoneutid50_Beta = ifelse(cohort_event2 %in% c("Post-Peak Cases-Beta", "Non-Cases"), Day29pseudoneutid50_Beta, NA),
           Day29pseudoneutid50_Delta = ifelse(cohort_event2 %in% c("Post-Peak Cases-Delta", "Non-Cases"), Day29pseudoneutid50_Delta, NA),
           Day29pseudoneutid50_Gamma = ifelse(cohort_event2 %in% c("Post-Peak Cases-Gamma", "Non-Cases"), Day29pseudoneutid50_Gamma, NA),
           Day29pseudoneutid50_Lambda = ifelse(cohort_event2 %in% c("Post-Peak Cases-Lambda", "Non-Cases"), Day29pseudoneutid50_Lambda, NA),
           Day29pseudoneutid50_Mu = ifelse(cohort_event2 %in% c("Post-Peak Cases-Mu", "Non-Cases"), Day29pseudoneutid50_Mu, NA),
           Day29pseudoneutid50_Zeta = ifelse(cohort_event2 %in% c("Post-Peak Cases-Zeta", "Non-Cases"), Day29pseudoneutid50_Zeta, NA)
           ) %>%
    # only keep variant data for the corresponding case strain, i.e, delete wild and delta values for beta cases
    ungroup() %>%
    select(one_of(paste0("Day29", assays_variant), paste0("Day29", assays_ancestral), "Region", "cohort_event", "wt.D29variant", "wt.D29")) %>%
    group_by(Region, cohort_event) %>%
    summarise(across(c(paste0("Day29",assays_variant)), ~ exp(sum(log(.x * wt.D29variant), na.rm=T) / sum(wt.D29variant))),
              across(c(paste0("Day29",assays)), ~ exp(sum(log(.x * wt.D29), na.rm=T) / sum(wt.D29)))) %>% # ancestral assays use wt.D29
    unique() %>%
    rename_with(~ gsub("Day29", "", .), everything()) %>%
    as.data.frame()

# stack with max and min values
max_min <- as.data.frame(rbind(rep(1.6,ncol(spider_plot_)), 
                 rep(0.1,ncol(spider_plot_)))) %>%
    setNames(colnames(spider_plot_)) %>%
    `rownames<-`(c("max","min"))

spider_plot_ <- rbind(max_min, 
                         spider_plot_)

# setup pdf file
for (ab in c("bAb", "nAb")) {
    
    for (reg in c(0, 1, 2)) {
        
        if (reg==0 && ab=="nAb") next # doesn't need nAb for US
        
        reg_lb = case_when(reg==0 ~ "NAM",
                           reg==1 ~ "LATAM",
                           reg==2 ~ "ZA",
                           TRUE ~ "")
        
        reg_lb_long = case_when(reg==0 ~ "Northern America",
                                reg==1 ~ "Latin America",
                                reg==2 ~ "Southern Africa",
                                TRUE ~ "")
        
        spider_plot <- spider_plot_[c(1,2), ] %>%
            bind_rows(
                spider_plot_[2:nrow(spider_plot_), ] %>%
                    filter(Region %in% reg)
            ) %>%
            mutate(Region = NULL) %>%
            select(if(ab=="bAb" && reg!=0) {starts_with("bindSpike")
            } else if (ab=="bAb" && reg==0) {matches("bindSpike$|bindSpike_P.1")
            } else if (ab=="nAb" && reg==1) {matches("pseudoneutid50$|pseudoneutid50_Zeta|pseudoneutid50_Mu|pseudoneutid50_Gamma|pseudoneutid50_Lambda")
            } else if (ab=="nAb" && reg==2) {matches("pseudoneutid50$|pseudoneutid50_Delta|pseudoneutid50_Beta")
            } else {contains("pseudoneutid50")})
        
        if (ncol(spider_plot)<=2) next
        
        # those without any data will have a weighted geomean equal to 1 because exp(0)=1, set these to NA
        spider_plot[spider_plot == 1] <- NA
        
        colnames(spider_plot) <- assay_metadata$assay_label[match(colnames(spider_plot), assay_metadata$assay)]
        colnames(spider_plot) <- gsub("PsV Neutralization to |PsV Neutralization |Binding Antibody to Spike |Binding Antibody to | Spike|Binding Antibody ", "", colnames(spider_plot))
    
        # figure starts here
        filename = paste0(save.results.to, "radar_weighted_geomean_Day29_vaccine_bseroneg_", ab, "_", reg_lb, " .pdf")
        pdf(filename, width=5.5, height=6.5)
        par(mfrow=c(1,1), mar=c(0.1,0.1,1,0.1))
        
        radarchart(spider_plot, 
                   axistype=1 , pcol=c("#FF6F1B","#0AB7C9"),
                   plwd=1.5, pty=c(15), plty=2,
                   #custom the grid
                   cglcol="grey", cglty=1, axislabcol="grey", cglwd=0.8, caxislabels=paste0("10^",seq(0.1,1.6,(1.6-0.1)/4)), 
                   #label size
                   vlcex=0.8,
                   #title
                   title=paste0("Geometric Mean ", ifelse(ab=="bAb", "of bAb Markers", "of nAb Markers"), " at Day 29", ", ", reg_lb_long),
                   #title size
                   cex.main=0.8)
        
        par(xpd=NA)
        
        #legend
        legend("bottom", legend=c("Peak-Peak Cases","Non-Cases"), lty=5, pch=c(15),
               col=c("#FF6F1B","#0AB7C9"), bty="n", ncol=1, cex=0.7,
               inset=c(-0.25,0))
    
        dev.off()
    }
}
if (length(dev.list()!=0)) {dev.off()}
