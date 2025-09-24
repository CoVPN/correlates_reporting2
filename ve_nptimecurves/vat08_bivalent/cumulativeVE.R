rm(list=ls(all=TRUE))
here::i_am("README.txt")
repoDir <- here::here()
datDir <- file.path(repoDir,"data")
outputDir <- file.path(repoDir, "output")
library(tidyverse)
source(file.path(repoDir, "code/confbandsurv_2020Jan.R"))

datPP <- read.csv(file.path(datDir, "vat08_combined_data_processed_20250321.csv")) 
datPP_D22 <- dplyr::select(datPP, all_of(c("Trt", "Bserostatus", "Trialstage", 
                                       "EventTimeOmicronD22M12hotdeck10", 
                                       "EventIndOmicronD22M12hotdeck10")))
#restrict to endpoints 13 days post dose 2
datPP_D22 <- filter(datPP_D22, EventTimeOmicronD22M12hotdeck10 >=14)
datPP_D22$Trt <- datPP_D22$Trt + 1
colnames(datPP_D22) <- c("tx", "Bserostatus", "Trialstage", "AVAL", "EVENT")


cumVE <- function(data, load.ind = TRUE, saveFile=NULL, saveDir=NULL, ngrid = 300){
  # set input parameters
  N <- 10000
  nsims <- 10000

  
  # set u1=t1 to the maximum of the smallest treatment-specific failure times
  t1 <- max(tapply(data$AVAL, data$tx, cgwtools::minn, nth = 1))
  
  # set u2=t2 to the minimum of the 2nd largest treatment-specific failure times
  t2 <- min(tapply(data$AVAL, data$tx, cgwtools::maxn, nth = 2))
  
  # save the output from confidenceband.loghazardratio() which is the bottleneck of this code
  if(load.ind){
    load(file.path(saveDir, saveFile))
  }else{
    logHRfit <- fit <- Confbandsurv(data$AVAL, data$EVENT, data$tx, nsims, t1, t2, ngrid, 
                                    "VE", bootstrap=FALSE)  
    save(logHRfit, file=file.path(saveDir, saveFile)) 
  }
  
  return(logHRfit)
}

cumVE_NonNaive_Stage2 <- cumVE(filter(datPP_D22, Bserostatus == 1 & Trialstage ==2),  load = FALSE,
                                      saveFile=paste0("cumVE_nonnaive.rds"), saveDir=outputDir, ngrid = 300)


readRdata <- function(logHRfit){
  fit <- data.frame(logHRfit[,1:6])
  colnames(fit) <- c("time", "pointEst", "ub", "lb","ub2","lb2")
  fit[,2:6] <- 100*fit[,2:6]
  return(fit)
}

fit <- readRdata(cumVE_NonNaive_Stage2)
p <- ggplot(data = fit)+
  geom_line(aes(x = time, y = pointEst), color = "orange", linewidth = 1.7)+
  geom_ribbon(aes( x = time, ymin = lb, ymax = ub, alpha = 0.2), fill = "orange")+ #turn off show.legend is a key to modify alpha in the legend later
  coord_cartesian(ylim = c(0,100), xlim = c(0, 365), expand = FALSE)+
  scale_y_continuous(breaks = c(0, 20, 40, 60, 80, 100))+
  scale_x_continuous(breaks = c(0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360))+
  guides(alpha = "none")+
  geom_hline(yintercept = 0)+
  xlab ("Time (Days) from Dose 2 Visit")+
  ylab( "Cumulative Incidence-based VE (%)")+
  theme_bw()+
  theme(legend.key.width = unit(2, "cm"),
        legend.key.height = unit(1, "cm"),
        plot.margin = unit(c(1,0.5,0.5,0.5), "cm"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.title=element_text(size = 15),
        legend.margin=margin(grid::unit(0,"cm")),
        legend.text=element_text(size=14),
        legend.position = c(0.9,0.9),
        legend.justification = c(1,0),
        legend.key = element_blank(),
        legend.background = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5, vjust = 3, size = 22),
        axis.title.x = element_text(size = 19, vjust = -1),
        axis.title.y = element_text(size = 19, hjust = 0.5, vjust = 1),
        axis.text.x = element_text(size = 17, colour = "black"),
        axis.text.y = element_text(size = 17, colour = "black"))


ggplot2::ggsave(filename = paste0("cumVE_nonnaive.pdf"), 
                plot = p, 
                path = outputDir,
                dpi = 320,
                width = 9, height = 9,units = "in") 
