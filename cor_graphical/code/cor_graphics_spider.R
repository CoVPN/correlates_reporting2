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
dat.cor.data.spider <- readRDS(here::here("data_clean", "longer_cor_data_plot1.rds"))

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
assays_id50 <- c("pseudoneutid50",
                 "pseudoneutid50_Beta","pseudoneutid50_Delta",
                 "pseudoneutid50_Gamma","pseudoneutid50_Lambda","pseudoneutid50_Mu","pseudoneutid50_Zeta")

dat.cor.data.spider.id50 <- dat.cor.data.spider %>% 
    filter(time == "Day 29" & Region %in% c(1, 2) & Trt=="Vaccine") %>%
    ungroup() %>%
    select(one_of(paste0("Day29", assays_id50), "Region", "cohort_event", "wt.D29variant")) %>%
    group_by(Region, cohort_event) %>%
    summarise(across(c(paste0("Day29",assays_id50)), ~ exp(sum(log(.x * wt.D29variant), na.rm=T) / sum(wt.D29variant)))) %>%
    unique() %>%
    as.data.frame()

rownames(dat.cor.data.spider.id50) <- paste("Region",dat.cor.data.spider.id50$Region, dat.cor.data.spider.id50$cohort_event)
dat.cor.data.spider.id50$Region <- NULL
dat.cor.data.spider.id50$cohort_event <- NULL
colnames(dat.cor.data.spider.id50) <- c("Reference","Beta","Delta","Gamma","Lambda","Mu","Zeta")

# stack with max and min values
max_min <- rbind(rep(1.1,ncol(dat.cor.data.spider.id50)), 
                 rep(0.1,ncol(dat.cor.data.spider.id50)))
colnames(max_min) <- c("Reference","Beta","Delta","Gamma","Lambda","Mu","Zeta")
rownames(max_min) <- c("max", "min")

dat.cor.data.spider.id50 <- rbind(max_min, 
                                  dat.cor.data.spider.id50)

# setup pdf file
for (outtype in c("PDF")) {
    
    # Latin America
    if(outtype=="PDF"){
        filename = paste0(save.results.to, "radar_plot_weighted_geomean_Day29_vaccine_bseroneg_NAb_LA.pdf")
        #pdf(filename, width=9.1, height=8.5)
        pdf(filename, width=5.5, height=6)
        par(mfrow=c(1,1), mar=c(0.1,0.1,1,0.1))
    }
    
    radarchart(dat.cor.data.spider.id50[c(1,2,3,4), c("Reference","Zeta","Mu","Gamma","Lambda")], 
               axistype=1 , pcol=c("#FF6F1B","#0AB7C9"),
               plwd=1.5, pty=c(15), plty=2,
               #custom the grid
               cglcol="grey", cglty=1, axislabcol="grey", cglwd=0.8, caxislabels=paste0("10^",seq(0.1,1.1,0.25)), 
               #label size
               vlcex=0.8,
               #title
               title="Geometric Mean of PsV-nAb ID50 at Day 29, in Latin America",
               #title size
               cex.main=0.8)
    
    par(xpd=NA)
    
    #legend
    legend("bottom", legend=c("Peak-Peak Cases","Non-Cases"), lty=5, pch=c(15),
           col=c("#FF6F1B","#0AB7C9"), bty="n", ncol=1, cex=0.7,
           inset=c(-0.25,0))

    if(outtype=="PDF") dev.off()
    
    
    # Southern America
    if(outtype=="PDF"){
        filename = paste0(save.results.to, "radar_plot_weighted_geomean_Day29_vaccine_bseroneg_NAb_SA.pdf")
        #pdf(filename, width=9.1, height=8.5)
        pdf(filename, width=5.5, height=6)
        par(mfrow=c(1,1), mar=c(0.1,0.1,1,0.1))
    }
    
    radarchart(dat.cor.data.spider.id50[c(1,2,5,6), c("Reference","Delta","Beta")], 
               axistype=1 , pcol=c("#FF6F1B","#0AB7C9"),
               plwd=1.5, pty=c(15, 17), plty=2,
               #custom the grid
               cglcol="grey", cglty=1, axislabcol="grey", cglwd=0.8, caxislabels=paste0("10^",seq(0.1,1.1,0.25)), 
               #label size
               vlcex=0.8,
               #title
               title="Geometric Mean of PsV-nAb ID50 at Day 29, in South Africa",
               #title size
               cex.main=0.8)
    
    par(xpd=NA)
    
    #legend
    legend("bottom", legend=c("Peak-Peak Cases","Non-Cases"), lty=5, pch=c(17),
           col=c("#FF6F1B","#0AB7C9"), bty="n", ncol=1, cex=0.7,
           inset=c(-0.25,0))
    
    if(outtype=="PDF") dev.off()
}

