rm(list=ls(all=TRUE))
here::i_am("README.txt")
repoDir <- here::here()
datDir <- file.path(repoDir,"data")
outputDir <- file.path(repoDir, "output")
library(tidyverse)
library(tidycmprsk)
datPP <- read.csv(file.path(datDir, "vat08_combined_data_processed_20250321.csv")) 
datPP_D22 <- dplyr::select(datPP, all_of(c("Trt", "Bserostatus", "Trialstage", 
                                           "EventTimeOmicronD22M12hotdeck10", 
                                           "EventIndOmicronD22M12hotdeck10", "prev_inf", 
                                           "D01_S_pos_only_in_non_naive_group")))
#remove observations with negative event time
datPP_D22 <- filter(datPP_D22, EventTimeOmicronD22M12hotdeck10 >=14)
colnames(datPP_D22) <- c("tx", "Bserostatus", "Trialstage", "AVAL", "EVENT", "prev_inf", "D01_S_pos_only_in_non_naive_group")
datPP_D22 <- filter(datPP_D22, Trialstage == 2)

#obtain cumulative incidence
#formula with Surv() on LHS and covariates on RHS. The event status variable must be a factor, 
#with the first level indicating 'censor' and subsequent levels the competing risks
cuminc_naive_placebo <- cuminc(Surv(AVAL, as.factor(EVENT)) ~ 1, filter(datPP_D22, Bserostatus == 0 & tx == 0))
cuminc_naive_vaccine<- cuminc(Surv(AVAL, as.factor(EVENT)) ~ 1, filter(datPP_D22, Bserostatus == 0 & tx == 1))
cuminc_nonnaive_placebo <- cuminc(Surv(AVAL, as.factor(EVENT)) ~ 1, filter(datPP_D22, Bserostatus == 1 & tx == 0))
cuminc_nonnaive_vaccine <- cuminc(Surv(AVAL, as.factor(EVENT)) ~ 1, filter(datPP_D22, Bserostatus == 1 & tx == 1))

#create cumulative incidence for a grid
cumIncF <- function(cumincObj, timeGrid){
    time <- cumincObj$time
    est <- cumincObj$estimate
    CL <- cumincObj$conf.low
    CU <- cumincObj$conf.high
    estStepF<- stepfun(time, c(0, est)) #right = FALSE, i.e., closed on the left; f(x_i) = y_i; x_0 = 0
    CLStepF<- stepfun(time, c(0, CL))
    CUStepF<- stepfun(time, c(0, CU))
    
    return(data.frame(timeGrid = timeGrid, 
                      est = estStepF(timeGrid),
                      CL = CLStepF(timeGrid),
                      CU = CUStepF(timeGrid)))
}

#plot cumulative incidence
timemax <- max(datPP_D22$AVAL)
timeGrid <- seq(0.1, timemax, 0.1)
nmax <- length(timeGrid)

df_naive_placebo <- cumIncF(cuminc_naive_placebo$tidy, timeGrid)
df_naive_vaccine <- cumIncF(cuminc_naive_vaccine$tidy, timeGrid)
df_nonnaive_placebo <- cumIncF(cuminc_nonnaive_placebo$tidy, timeGrid)
df_nonnaive_vaccine <- cumIncF(cuminc_nonnaive_vaccine$tidy, timeGrid)

df_naive <- data.frame(rbind(df_naive_placebo, df_naive_vaccine), tx = c(rep("Placebo", nmax), rep("Vaccine", nmax)))
df_nonnaive <- data.frame(rbind(df_nonnaive_placebo, df_nonnaive_vaccine), tx = c(rep("Placebo", nmax), rep("Vaccine", nmax)))


#naive 
rangeData <- data.frame(cbind(c(min(df_naive$timeGrid), max(df_naive$timeGrid)), c(0,1)))
colnames(rangeData) <- c("timeGrid","cumInc")
xBreaksAt <- c(0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360)
xlimits <- c(0, 365)
df_naive$tx <- factor(df_naive$tx, levels = c("Vaccine", "Placebo"))
p <- ggplot(data = rangeData) +
  geom_line( aes(x = timeGrid, y = est, colour = tx), lwd = 1.2, data = df_naive)+
  geom_ribbon( aes(x = timeGrid, ymin = CL, ymax = CU, fill = tx), alpha = 0.3, data = df_naive) + 
  scale_x_continuous(limits = xlimits, breaks = xBreaksAt, minor_breaks = NULL , expand = c(0, 0))+
  scale_y_continuous(limits = c(0, 0.18), breaks = c(0, 0.025, 0.050, 0.075, 0.100, 0.125, 0.150, 0.175), 
                     minor_breaks = NULL, expand = c(0, 0))+
  coord_cartesian( clip = "off")+
  theme_classic()+
  theme(plot.margin = margin(t = 1, r = 1, b = 8, l = 2, unit = "cm"),
        axis.text.x = element_text(size = 17),
        axis.text.y = element_text(size = 17),
        axis.title.x = element_text(size = 19, margin = margin(t = 10)),
        axis.title.y = element_text(size = 19, margin = margin(r = 10)),
        legend.key.size = unit(1, "cm"),
        legend.position=c(0.25,0.8),legend.direction="vertical",
        legend.margin=margin(grid::unit(0.1,"cm")),
        legend.key = element_blank(),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.text=element_text(size=16,face="plain"),
        legend.key.height=grid::unit(0.8,"cm"),
        legend.key.width=grid::unit(0.8,"cm"))+
  ylab("Cumulative Incidence of Omicron COVID-19")+
  xlab("Time (Days) from Dose 2")
p <-  p + scale_color_manual(name="",
                     breaks = c("Vaccine", "Placebo"),
                     values = c("Vaccine" = "steelblue1", "Placebo" = "gray40"),
                     labels = c("Vaccine", "Placebo")) +
      scale_fill_manual(name="",
                        breaks = c("Vaccine", "Placebo"),
                        values = c("Vaccine" = "steelblue1", "Placebo" = "gray40"),
                        labels = c("Vaccine", "Placebo"))

#create cumulative incidence table
cumIncTable <-  tibble("time" = numeric(), "tx" = character(), "No. at Risk" = numeric(), 
                       "cum No. of Cases" = numeric())

groupsLabel <- c("Vaccine" , "Placebo")
names(groupsLabel) <- c("1" , "0")
tableData <- filter(datPP_D22, Bserostatus == 0)
for(x in xBreaksAt){
  for(group in names(groupsLabel)){
    cumIncTable  <- add_row(.data = cumIncTable , "time" = x, "tx" = groupsLabel[group] , 
                            "No. at Risk" = length(filter(tableData, tx == group  & AVAL >= x)$EVENT),
                            "cum No. of Cases" =  sum(filter(tableData, tx == group  & AVAL <= x)$EVENT))
  }
  
}

NoAtRisk <- pivot_wider(cumIncTable[,-4], names_from = "tx", values_from = "No. at Risk")
cumCases <- pivot_wider(cumIncTable[,-3], names_from = "tx", values_from = "cum No. of Cases")
size.text = 5.2
#Add the number at risk
leftPos = -55
rowPos = -0.025
p2 <- p + geom_text(x = leftPos, y = rowPos, fontface = "bold", size = size.text , label = "Number at Risk", hjust = 0, check_overlap = TRUE)
for(i in 1:2){
  rowPos = rowPos-0.012
  p2 <- p2 + geom_text(x = leftPos, y = rowPos , size = size.text , label = groupsLabel[i], hjust = 0,check_overlap = TRUE)
  
  for(j in 1:length(NoAtRisk$time)){
    time = NoAtRisk$time[j]
    
    label = as.character(NoAtRisk[j, i+1])
    p2 <- p2 + geom_text(x = time, y = rowPos, size = size.text , label = label, hjust = 0.5,check_overlap = TRUE)
  }
}
#Add the cumulative incidence
rowPos = rowPos-0.012
p3 <- p2 + geom_text(x = leftPos, y = rowPos, fontface = "bold", size = size.text , label = "Cumulative Omicron Events", lineheight = 0.7, hjust = 0, check_overlap = TRUE)

for(i in 1:2){
  rowPos = rowPos-0.012
  p3 <- p3 + geom_text(x = leftPos, y = rowPos , size = size.text , label = groupsLabel[i], hjust = 0,check_overlap = TRUE)
  
  for(j in 1:length(cumCases$time)){
    time = cumCases$time[j]
    label = as.character(cumCases[j, i+1])
    p3 <- p3 + geom_text(x = time, y = rowPos, size = size.text , label = label, hjust = 0.5,check_overlap = TRUE)
  }
}

ggplot2::ggsave(filename = paste0("stage2_pp_naivecumInc.pdf"), 
                plot = p3, 
                path = outputDir,
                dpi = 320,
                width = 10, height = 10, units = "in") 



#nonnaive group
rangeData <- data.frame(cbind(c(min(df_nonnaive$timeGrid), max(df_nonnaive$timeGrid)), c(0,1)))
colnames(rangeData) <- c("timeGrid","cumInc")
xBreaksAt <- c(0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360)
xlimits <- c(0, 365)
df_nonnaive$tx <- factor(df_nonnaive$tx, levels = c("Vaccine", "Placebo"))
#restrict to time to 366 days

breaksAt <- c(0, 0.005, 0.010, 0.015, 0.020, 0.025, 0.030)


p <- ggplot(data = rangeData) +
  geom_line( aes(x = timeGrid, y = est, colour = tx), lwd = 1.2, data = filter(df_nonnaive, timeGrid <=365))+
  geom_ribbon( aes(x = timeGrid, ymin = CL, ymax = CU, fill = tx), alpha = 0.3, data = filter(df_nonnaive, timeGrid <=365)) + 
  scale_x_continuous(limits = xlimits, breaks = xBreaksAt, minor_breaks = NULL , expand = c(0, 0))+
  scale_y_continuous(limits = c(0, 0.03), breaks = breaksAt, 
                     minor_breaks = NULL, expand = c(0, 0))+
  coord_cartesian( clip = "off")+
  theme_classic()+
  theme(plot.margin = margin(t = 1, r = 1, b = 8, l = 2, unit = "cm"),
        axis.text.x = element_text(size = 17),
        axis.text.y = element_text(size = 17),
        axis.title.x = element_text(size = 19, margin = margin(t = 12)),
        axis.title.y = element_text(size = 19, margin = margin(r = 12)),
        legend.key.size = unit(1, "cm"),
        legend.position=c(0.25,0.8),legend.direction="vertical",
        legend.margin=margin(grid::unit(0.1,"cm")),
        legend.key = element_blank(),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.text=element_text(size=16,face="plain"),
        legend.key.height=grid::unit(0.8,"cm"),
        legend.key.width=grid::unit(0.8,"cm"))+
  ylab("Cumulative Incidence of Omicron COVID-19")+
  xlab("Time (Days) from Dose 2")
p <-  p + scale_color_manual(name="",
                             breaks = c("Vaccine", "Placebo"),
                             values = c("Vaccine" = "steelblue1", "Placebo" = "gray40"),
                             labels = c("Vaccine", "Placebo")) +
  scale_fill_manual(name="",
                    breaks = c("Vaccine", "Placebo"),
                    values = c("Vaccine" = "steelblue1", "Placebo" = "gray40"),
                    labels = c("Vaccine", "Placebo"))

#create cumulative incidence table
cumIncTable <-  tibble("time" = numeric(), "tx" = character(), "No. at Risk" = numeric(), 
                       "cum No. of Cases" = numeric())

groupsLabel <- c("Vaccine" , "Placebo")
names(groupsLabel) <- c("1" , "0")
tableData <- filter(datPP_D22, Bserostatus == 1)
for(x in xBreaksAt){
  for(group in names(groupsLabel)){
    cumIncTable  <- add_row(.data = cumIncTable , "time" = x, "tx" = groupsLabel[group] , 
                            "No. at Risk" = length(filter(tableData, tx == group  & AVAL >= x)$EVENT),
                            "cum No. of Cases" =  sum(filter(tableData, tx == group  & AVAL <= x)$EVENT))
  }
  
}

NoAtRisk <- pivot_wider(cumIncTable[,-4], names_from = "tx", values_from = "No. at Risk")
cumCases <- pivot_wider(cumIncTable[,-3], names_from = "tx", values_from = "cum No. of Cases")
size.text = 5.2
#Add the number at risk
leftPos = -55
rowPos = -0.0045
p2 <- p + geom_text(x = leftPos, y = rowPos, fontface = "bold", size = size.text , label = "Number at Risk", hjust = 0, check_overlap = TRUE)
for(i in 1:2){
  rowPos = rowPos- 0.002
  p2 <- p2 + geom_text(x = leftPos, y = rowPos , size = size.text , label = groupsLabel[i], hjust = 0,check_overlap = TRUE)
  
  for(j in 1:length(NoAtRisk$time)){
    time = NoAtRisk$time[j]
    
    label = as.character(NoAtRisk[j, i+1])
    p2 <- p2 + geom_text(x = time, y = rowPos, size = size.text , label = label, hjust = 0.5,check_overlap = TRUE)
  }
}
#Add the cumulative incidence
rowPos = rowPos- 0.002
p3 <- p2 + geom_text(x = leftPos, y = rowPos, fontface = "bold", size = size.text , label = "Cumulative Omicron Events", lineheight = 0.7, hjust = 0, check_overlap = TRUE)

for(i in 1:2){
  rowPos = rowPos- 0.002
  p3 <- p3 + geom_text(x = leftPos, y = rowPos , size = size.text , label = groupsLabel[i], hjust = 0,check_overlap = TRUE)
  
  for(j in 1:length(cumCases$time)){
    time = cumCases$time[j]
    label = as.character(cumCases[j, i+1])
    p3 <- p3 + geom_text(x = time, y = rowPos, size = size.text , label = label, hjust = 0.5,check_overlap = TRUE)
  }
}

ggplot2::ggsave(filename = paste0("stage2_pp_nonnaivecumInc.pdf"), 
                plot = p3, 
                path = outputDir,
                dpi = 320,
                width = 10, height = 10, units = "in") 



