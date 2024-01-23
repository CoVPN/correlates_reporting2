library(survey)
library(kyotil)

dat_mapped=read.csv('/trials/covpn/COVAILcorrelates/analysis/mapping_immune_correlates/adata/covail_mapped_data_20240112.csv')
dat_proc=read.csv('/trials/covpn/COVAILcorrelates/analysis/correlates/adata/covail_data_processed_20240122.csv')

assay_metadata=read.csv('~/correlates_reporting2/assay_metadata/covail_assay_metadata.csv')
assays=assay_metadata$assay
assays


################################################################################
# coxph showing all coefficients

dat=subset(dat.mock, ph1.D15 & TrtonedosemRNA==1) 
coxph(Surv(COVIDtimeD22toD181, COVIDIndD22toD181) ~ FOIstandardized + standardized_risk_score + naive + Day15pseudoneutid50_MDW, dat)



################################################################################
# correlation between m15 and baseline markers

dat=subset(dat.mock, ph1.D15 & TrtonedosemRNA==1) 

mypairs(dat[,paste0("B",assays)])

mypairs(dat[,paste0("Day15",assays)])

mypairs(dat[,paste0(c("B","Day15", "Delta15overB"),assays[1])])



################################################################################
# marker missingness

table(dat_proc$AsympInfectIndD15to29, dat_proc$ph1.D15)



table(dat_proc$arm, dat_proc$ph1.D15)

table(dat_mapped$ph1.D15, dat_mapped$Immunemarkerset, dat_mapped$arm==3)


# across assays, all or none
summary(subset(dat_proc, ph1.D15 & is.na(Day29pseudoneutid50_D614G))["Day29"%.%assays[1:5]])
summary(subset(dat_proc, ph1.D15 & is.na(Day91pseudoneutid50_D614G))["Day91"%.%assays[1:5]])
summary(subset(dat_proc, ph1.D15 & is.na(Day181pseudoneutid50_D614G))["Day181"%.%assays[1:5]])

dat_proc$kp = dat_proc$ph1.D15==1 & dat_proc$COVIDIndD22toend!=1 & dat_proc$AsympInfectIndD15to271!=1 
dat_proc$kp = dat_proc$ph1.D15==1 
for (i in 1:1) {
  with(dat_proc[dat_proc$kp,], print(table(!is.na(get("Day29"%.%assays[i])), !is.na(get("Day15"%.%assays[i])))))
  with(dat_proc[dat_proc$kp,], print(table(!is.na(get("Day91"%.%assays[i])), !is.na(get("Day15"%.%assays[i])))))
  with(dat_proc[dat_proc$kp,], print(table(!is.na(get("Day181"%.%assays[i])), !is.na(get("Day15"%.%assays[i])))))
}


myboxplot(
list(dat_proc$Day15pseudoneutid50_MDW[dat_proc$treatment_actual=="Omicron + Wildtype/Prototype (Pfizer 1)"], 
     dat_proc$Day15pseudoneutid50_MDW[dat_proc$treatment_actual=="1 Dose Omicron + Prototype (Moderna)"])
)

myboxplot(
  list(dat_proc$Day15pseudoneutid50_D614G[dat_proc$treatment_actual=="Omicron + Wildtype/Prototype (Pfizer 1)"], 
       dat_proc$Day15pseudoneutid50_D614G[dat_proc$treatment_actual=="1 Dose Omicron + Prototype (Moderna)"])
)



# my.interaction.plot(subset(dat_proc, ph1==1, 
#                            c(Bpseudoneutid50_BA.1, Day15pseudoneutid50_BA.1)), 
#                     x.ori = 0, xaxislabels = c("B", "D15"), cex.axis = 1, add = T, xlab = "", ylab = "", pcol = NULL, lcol = NULL)  
# 
# my.interaction.plot(subset(dat_proc, ph1==1, 
#                            c(Day15pseudoneutid50_BA.1, Day29pseudoneutid50_BA.1)), 
#                     x.ori = 1, xaxislabels = c("D15", "D29"), cex.axis = 1, add = T, xlab = "", ylab = "", pcol = NULL, lcol = NULL)  
# 
# my.interaction.plot(subset(dat_proc, ph1==1, 
#                            c(Day29pseudoneutid50_BA.1, Day91pseudoneutid50_BA.1)), 
#                     x.ori = 2, xaxislabels = c("D29", "D91"), cex.axis = 1, add = T, xlab = "", ylab = "", pcol = NULL, lcol = NULL)  
# 
# my.interaction.plot(subset(dat_proc, ph1==1, 
#                            c(Day91pseudoneutid50_BA.1, Day181pseudoneutid50_BA.1)), 
#                     x.ori = 3, xaxislabels = c("D91", "D181"), cex.axis = 1, add = T, xlab = "", ylab = "", pcol = NULL, lcol = NULL)  
# 


{
# marker level over time
par(mfrow=c(4,1), mar=c(2,4,1,1))

dat_proc$ph1 = dat_proc$ph1.D15==1 & dat_proc$COVIDIndD22toend!=1 & dat_proc$AsympInfectIndD15to181!=1 & dat_proc$AsympInfectIndD182to271!=1

plot(0,0,type='n', xlim=c(1,5), ylim=c(1,5), xaxt="n", xlab="", ylab="ID50_BA.1")

my.interaction.plot(subset(dat_proc, ph1==1 & Bpseudoneutid50_BA.1<1.5, 
                           c(Bpseudoneutid50_BA.1, Day15pseudoneutid50_BA.1)), 
                    x.ori = 0, xaxislabels = c("B", "D15"), cex.axis = 1, add = T, xlab = "", ylab = "", pcol = NULL, lcol = NULL)  

my.interaction.plot(subset(dat_proc, ph1==1 & Bpseudoneutid50_BA.1<1.5, 
                           c(Day15pseudoneutid50_BA.1, Day29pseudoneutid50_BA.1)), 
                    x.ori = 1, xaxislabels = c("D15", "D29"), cex.axis = 1, add = T, xlab = "", ylab = "", pcol = NULL, lcol = NULL)  

my.interaction.plot(subset(dat_proc, ph1==1 & Bpseudoneutid50_BA.1<1.5, 
                           c(Day29pseudoneutid50_BA.1, Day91pseudoneutid50_BA.1)), 
                    x.ori = 2, xaxislabels = c("D29", "D91"), cex.axis = 1, add = T, xlab = "", ylab = "", pcol = NULL, lcol = NULL)  

my.interaction.plot(subset(dat_proc, ph1==1 & Bpseudoneutid50_BA.1<1.5, 
                           c(Day91pseudoneutid50_BA.1, Day181pseudoneutid50_BA.1)), 
                    x.ori = 3, xaxislabels = c("D91", "D181"), cex.axis = 1, add = T, xlab = "", ylab = "", pcol = NULL, lcol = NULL)  


plot(0,0,type='n', xlim=c(1,5), ylim=c(1,5), xaxt="n", xlab="", ylab="ID50_BA.1")

my.interaction.plot(subset(dat_proc, ph1==1 & Bpseudoneutid50_BA.1>1.5, 
                           c(Bpseudoneutid50_BA.1, Day15pseudoneutid50_BA.1)), 
                    x.ori = 0, xaxislabels = c("B", "D15"), cex.axis = 1, add = T, xlab = "", ylab = "", pcol = NULL, lcol = NULL)  

my.interaction.plot(subset(dat_proc, ph1==1 & Bpseudoneutid50_BA.1>1.5, 
                           c(Day15pseudoneutid50_BA.1, Day29pseudoneutid50_BA.1)), 
                    x.ori = 1, xaxislabels = c("D15", "D29"), cex.axis = 1, add = T, xlab = "", ylab = "", pcol = NULL, lcol = NULL)  

my.interaction.plot(subset(dat_proc, ph1==1 & Bpseudoneutid50_BA.1>1.5, 
                           c(Day29pseudoneutid50_BA.1, Day91pseudoneutid50_BA.1)), 
                    x.ori = 2, xaxislabels = c("D29", "D91"), cex.axis = 1, add = T, xlab = "", ylab = "", pcol = NULL, lcol = NULL)  

my.interaction.plot(subset(dat_proc, ph1==1 & Bpseudoneutid50_BA.1>1.5, 
                           c(Day91pseudoneutid50_BA.1, Day181pseudoneutid50_BA.1)), 
                    x.ori = 3, xaxislabels = c("D91", "D181"), cex.axis = 1, add = T, xlab = "", ylab = "", pcol = NULL, lcol = NULL)  



dat_proc$ph1 = dat_proc$ph1.D92 & dat_proc$ph1.D15==1 & dat_proc$COVIDIndD22toend!=1 & dat_proc$AsympInfectIndD15to181!=1 & dat_proc$AsympInfectIndD182to271!=1

plot(0,0,type='n', xlim=c(1,5), ylim=c(1,5), xaxt="n", xlab="", ylab="ID50_BA.1")

my.interaction.plot(subset(dat_proc, ph1==1 & Bpseudoneutid50_BA.1<1.5, 
                           c(Bpseudoneutid50_BA.1, Day15pseudoneutid50_BA.1)), 
                    x.ori = 0, xaxislabels = c("B", "D15"), cex.axis = 1, add = T, xlab = "", ylab = "", pcol = NULL, lcol = NULL)  

my.interaction.plot(subset(dat_proc, ph1==1 & Bpseudoneutid50_BA.1<1.5, 
                           c(Day15pseudoneutid50_BA.1, Day29pseudoneutid50_BA.1)), 
                    x.ori = 1, xaxislabels = c("D15", "D29"), cex.axis = 1, add = T, xlab = "", ylab = "", pcol = NULL, lcol = NULL)  

my.interaction.plot(subset(dat_proc, ph1==1 & Bpseudoneutid50_BA.1<1.5, 
                           c(Day29pseudoneutid50_BA.1, Day91pseudoneutid50_BA.1)), 
                    x.ori = 2, xaxislabels = c("D29", "D91"), cex.axis = 1, add = T, xlab = "", ylab = "", pcol = NULL, lcol = NULL)  

my.interaction.plot(subset(dat_proc, ph1==1 & Bpseudoneutid50_BA.1<1.5, 
                           c(Day91pseudoneutid50_BA.1, Day181pseudoneutid50_BA.1)), 
                    x.ori = 3, xaxislabels = c("D91", "D181"), cex.axis = 1, add = T, xlab = "", ylab = "", pcol = NULL, lcol = NULL)  


plot(0,0,type='n', xlim=c(1,5), ylim=c(1,5), xaxt="n", xlab="", ylab="ID50_BA.1")

my.interaction.plot(subset(dat_proc, ph1==1 & Bpseudoneutid50_BA.1>1.5, 
                           c(Bpseudoneutid50_BA.1, Day15pseudoneutid50_BA.1)), 
                    x.ori = 0, xaxislabels = c("B", "D15"), cex.axis = 1, add = T, xlab = "", ylab = "", pcol = NULL, lcol = NULL)  

my.interaction.plot(subset(dat_proc, ph1==1 & Bpseudoneutid50_BA.1>1.5, 
                           c(Day15pseudoneutid50_BA.1, Day29pseudoneutid50_BA.1)), 
                    x.ori = 1, xaxislabels = c("D15", "D29"), cex.axis = 1, add = T, xlab = "", ylab = "", pcol = NULL, lcol = NULL)  

my.interaction.plot(subset(dat_proc, ph1==1 & Bpseudoneutid50_BA.1>1.5, 
                           c(Day29pseudoneutid50_BA.1, Day91pseudoneutid50_BA.1)), 
                    x.ori = 2, xaxislabels = c("D29", "D91"), cex.axis = 1, add = T, xlab = "", ylab = "", pcol = NULL, lcol = NULL)  

my.interaction.plot(subset(dat_proc, ph1==1 & Bpseudoneutid50_BA.1>1.5, 
                           c(Day91pseudoneutid50_BA.1, Day181pseudoneutid50_BA.1)), 
                    x.ori = 3, xaxislabels = c("D91", "D181"), cex.axis = 1, add = T, xlab = "", ylab = "", pcol = NULL, lcol = NULL)  

}


sum(dat_proc$ph1.D15==1)

dat_proc$kp = dat_proc$ph1.D15==1 & dat_proc$COVIDIndD22toend!=1 & dat_proc$AsympInfectIndD15to271!=1 &
  !is.na(dat_proc$Day29pseudoneutid50_D614G) & !is.na(dat_proc$Day91pseudoneutid50_D614G) & !is.na(dat_proc$Day181pseudoneutid50_D614G) 

nrow(subset(dat_proc, kp))

nrow(subset(dat_proc, kp & Bpseudoneutid50_BA.1>1.5))

nrow(subset(dat_proc, kp & Bpseudoneutid50_BA.1<1.5))






################################################################################
myboxplot(dat_mapped[, c("B"%.%assays[1:5], "Day15"%.%assays[1:5])], names=sub("pseudoneutid50_", "", rep(assays[1:5],2)))



summary(dat_mapped)
sort(names(dat_mapped))

table(dat_mapped$treatment_actual, dat_mapped$stage)


table(dat_mapped$Immunemarkerset, dat_mapped$Perprotocol, useNA='ifany')
table(dat_mapped$Immunemarkerset, dat_mapped$eligibility_deviation, useNA='ifany')


# Perprotocol is more stringent than D1_D15_flag
table(dat_mapped$Perprotocol, dat_mapped$D1_D15_flag, useNA='ifany')


tmp = subset(dat_mapped, Perprotocol==1, c(D15neutBA1, D15neutBeta, D15neutD614G, D15neutDelta))
table(complete.cases(tmp))

tmp = subset(dat_mapped, D1_D15_flag==1, c(D29neutBA4BA5, D15neutBA1, D15neutBeta, D15neutD614G, D15neutDelta, D1neutBA1, D1neutBeta, D1neutD614G, D1neutDelta))
table(complete.cases(tmp))

tmp = subset(dat_mapped, Perprotocol==1, c(D29neutBA4BA5, D15neutBA1, D15neutBeta, D15neutD614G, D15neutDelta, D1neutBA1, D1neutBeta, D1neutD614G, D1neutDelta))
table(complete.cases(tmp))

