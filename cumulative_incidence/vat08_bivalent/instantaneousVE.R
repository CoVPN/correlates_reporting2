rm(list=ls(all=TRUE))
here::i_am("README.txt")
repoDir <- here::here()
datDir <- file.path(repoDir,"data")
outputDir <- file.path(repoDir, "output")
library(tidyverse)
source(file.path(repoDir, "confidenceband.loghazardratio.R"))

datPP <- read.csv(config::get(config = "vat08_combined", file="../../correlates_reporting2/config.yml")$data_cleaned)


datPP <- dplyr::select(datPP, all_of(c("Ptid","Trt", "Bserostatus", "Trialstage", 
                                       "EventTimeOmicronD22M12hotdeck10", "EventIndOmicronD22M12hotdeck10")))
#remove observations with negative event time
datPP <- filter(datPP, EventTimeOmicronD22M12hotdeck10 >=14)
datPP$Trt <- datPP$Trt + 1
colnames(datPP) <- c("Ptid","tx", "Bserostatus", "Trialstage", "AVAL", "EVENT")


smoothHazVE <- function(data, nBoot=NULL, load.ind = TRUE, saveFile=NULL, saveDir=NULL){
  # set input parameters
  N <- 10000
  
  # set u1=t1 to the maximum of the smallest treatment-specific failure times
  t1 <- max(tapply(data$AVAL, data$tx, cgwtools::minn, nth = 1))
  
  # set u2=t2 to the minimum of the 2nd largest treatment-specific failure times
  t2 <- min(tapply(data$AVAL, data$tx, cgwtools::maxn, nth = 2))
  
  if (is.null(nBoot)){
    band12 <- band11 <- (t2 - t1) / 2  
  } else {
    band12 <- band11 <- 0  
  }
  
  # default settings for the rest of the parameters
  band1 <- 0
  band2 <- 0
  ngrid <- 100
  biasadjust <- TRUE
  tailsl <- TRUE
  tailsu <- TRUE
  Nv <- 0
  
  # save the output from confidenceband.loghazardratio() which is the bottleneck of this code
  if(load.ind){
    load(file.path(saveDir, saveFile))
  }else{
    logHRfit <- confidenceband.loghazardratio(data$AVAL, data$EVENT, data$tx, t1, t2, t1, t2, N, band11, band1, band12, band2, ngrid, biasadjust, tailsl, tailsu, nBoot, Nv)  
    save(logHRfit, file=file.path(saveDir, saveFile)) 
  }
  
  return(logHRfit)
}

smoothHazVE_NonNaive_Stage2 <- smoothHazVE(filter(datPP, Bserostatus == 1 & Trialstage ==2), nBoot=500, load = FALSE,
                                      saveFile=paste0("smoothHazVE_nonnaive.rds"), saveDir=outputDir)



readRdata <- function(logHRfit){
  fit <- data.frame(logHRfit[,1:6])
  colnames(fit) <- c("time", "pointEst", "ub", "lb","ub2","lb2")
  fit[,2:6] <- 100*(1-exp(fit[,2:6]))
  return(fit)
}

allFits <- readRdata(smoothHazVE_NonNaive_Stage2)
fit <- drop_na(allFits)
p <- ggplot(data = fit)+
  geom_line(aes(x = time, y = pointEst), color = "orange", linewidth = 1.7)+
  geom_ribbon(aes( x = time, ymin = lb, ymax = ub, alpha = 0.2) , fill = "orange")+ 
  coord_cartesian( ylim = c(-100, 100), expand = FALSE)+
  scale_y_continuous(breaks = c(-100,-80,-60,-40, -20, 0, 20, 40, 60, 80, 100))+
  scale_x_continuous(lim = c(0,245), breaks = c(0, 30, 60, 90, 120, 150, 180, 210, 240))+
  guides(alpha = "none")+
  geom_hline(yintercept = 0)+
  xlab ("Time (Days) from Dose 2 Visit")+
  ylab( "Instantaneous Hazard-based VE (%)")+
  theme_bw()+
  theme(legend.key.width = unit(2, "cm"),
        legend.key.height = unit(1, "cm"),
        plot.margin = unit(c(1,0.5,0.5,0.5), "cm"),
        legend.title=element_text(size = 15),
        legend.margin=margin(grid::unit(0,"cm")),
        legend.text=element_text(size=14),
        legend.position = c(0.8,0.85),
        legend.justification = c(1,0),
        legend.key = element_blank(),
        legend.background = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5, vjust = 3, size = 22),
        axis.title.x = element_text(size = 19, vjust = -1),
        axis.title.y = element_text(size = 19, hjust = 0.5, vjust = 1),
        axis.text.x = element_text(size = 17, colour = "black"),
        axis.text.y = element_text(size = 17, colour = "black"))

ggplot2::ggsave(filename = paste0("smoothHazVE_nonnaive.pdf"), 
                plot = p, 
                path = outputDir,
                dpi = 320,
                width = 9, height = 9,units = "in") 

