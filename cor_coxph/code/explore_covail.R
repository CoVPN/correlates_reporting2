library(survey)
library(kyotil)

dat_mapped=read.csv('/trials/covpn/COVAILcorrelates/analysis/mapping_immune_correlates/adata/covail_mapped_data_20240102.csv')
dat_proc=read.csv('/trials/covpn/COVAILcorrelates/analysis/correlates/adata/covail_data_processed_20240103.csv')

assay_metadata=read.csv('~/correlates_reporting2/assay_metadata/covail_assay_metadata.csv')
assays=assay_metadata$assay
assays


################################################################################
# marker missingness

for (i in 1:5) {
  print(table(!is.na(dat_proc[["B"%.%assays[i]]][dat_proc$ph1.D15==1 & dat_proc$arm!=3]), !is.na(dat_proc[["Day15"%.%assays[i]]][dat_proc$ph1.D15==1 & dat_proc$arm!=3])))
}

for (i in 1:5) {
  print(table(!is.na(dat_proc[["Day29"%.%assays[i]]][dat_proc$ph1.D15==1 & dat_proc$arm!=3]), !is.na(dat_proc[["Day15"%.%assays[i]]][dat_proc$ph1.D15==1 & dat_proc$arm!=3])))
}

for (i in 1:5) {
  print(table(!is.na(dat_proc[["Day91"%.%assays[i]]][dat_proc$ph1.D15==1 & dat_proc$arm!=3]), !is.na(dat_proc[["Day15"%.%assays[i]]][dat_proc$ph1.D15==1 & dat_proc$arm!=3])))
}


par(mfrow=c(3,1))

plot(0,0,type='n', xlim=c(1,4), ylim=c(1,5), xaxt="n", xlab="")

my.interaction.plot(subset(dat_proc, ph1.D15 & arm!=3 & !is.na(Bpseudoneutid50_BA.1) & !is.na(Day15pseudoneutid50_BA.1) & !is.na(Day29pseudoneutid50_BA.1) & !is.na(Day91pseudoneutid50_BA.1), 
                           c(Bpseudoneutid50_BA.1, Day15pseudoneutid50_BA.1)), 
                    x.ori = 0, xaxislabels = c("B", "D15"), cex.axis = 1, add = T, xlab = "", ylab = "", pcol = NULL, lcol = NULL)  

my.interaction.plot(subset(dat_proc, ph1.D15 & arm!=3 & !is.na(Bpseudoneutid50_BA.1) & !is.na(Day15pseudoneutid50_BA.1) & !is.na(Day29pseudoneutid50_BA.1) & !is.na(Day91pseudoneutid50_BA.1), 
                           c(Day15pseudoneutid50_BA.1, Day29pseudoneutid50_BA.1)), 
                    x.ori = 1, xaxislabels = c("D15", "D29"), cex.axis = 1, add = T, xlab = "", ylab = "", pcol = NULL, lcol = NULL)  

my.interaction.plot(subset(dat_proc, ph1.D15 & arm!=3 & !is.na(Bpseudoneutid50_BA.1) & !is.na(Day15pseudoneutid50_BA.1) & !is.na(Day29pseudoneutid50_BA.1) & !is.na(Day91pseudoneutid50_BA.1), 
                           c(Day29pseudoneutid50_BA.1, Day91pseudoneutid50_BA.1)), 
                    x.ori = 2, xaxislabels = c("D29", "D91"), cex.axis = 1, add = T, xlab = "", ylab = "", pcol = NULL, lcol = NULL)  


plot(0,0,type='n', xlim=c(1,4), ylim=c(1,5), xaxt="n", xlab="")

my.interaction.plot(subset(dat_proc, ph1.D15 & arm!=3 & Bpseudoneutid50_BA.1<1.5 & !is.na(Bpseudoneutid50_BA.1) & !is.na(Day15pseudoneutid50_BA.1) & !is.na(Day29pseudoneutid50_BA.1) & !is.na(Day91pseudoneutid50_BA.1), 
                           c(Bpseudoneutid50_BA.1, Day15pseudoneutid50_BA.1)), 
                    x.ori = 0, xaxislabels = c("B", "D15"), cex.axis = 1, add = T, xlab = "", ylab = "", pcol = NULL, lcol = NULL)  

my.interaction.plot(subset(dat_proc, ph1.D15 & arm!=3 & Bpseudoneutid50_BA.1<1.5 & !is.na(Bpseudoneutid50_BA.1) & !is.na(Day15pseudoneutid50_BA.1) & !is.na(Day29pseudoneutid50_BA.1) & !is.na(Day91pseudoneutid50_BA.1), 
                           c(Day15pseudoneutid50_BA.1, Day29pseudoneutid50_BA.1)), 
                    x.ori = 1, xaxislabels = c("D15", "D29"), cex.axis = 1, add = T, xlab = "", ylab = "", pcol = NULL, lcol = NULL)  

my.interaction.plot(subset(dat_proc, ph1.D15 & arm!=3 & Bpseudoneutid50_BA.1<1.5 & !is.na(Bpseudoneutid50_BA.1) & !is.na(Day15pseudoneutid50_BA.1) & !is.na(Day29pseudoneutid50_BA.1) & !is.na(Day91pseudoneutid50_BA.1), 
                           c(Day29pseudoneutid50_BA.1, Day91pseudoneutid50_BA.1)), 
                    x.ori = 2, xaxislabels = c("D29", "D91"), cex.axis = 1, add = T, xlab = "", ylab = "", pcol = NULL, lcol = NULL)  


plot(0,0,type='n', xlim=c(1,4), ylim=c(1,5), xaxt="n", xlab="")

my.interaction.plot(subset(dat_proc, ph1.D15 & arm!=3 & Bpseudoneutid50_BA.1>1.5 & !is.na(Bpseudoneutid50_BA.1) & !is.na(Day15pseudoneutid50_BA.1) & !is.na(Day29pseudoneutid50_BA.1) & !is.na(Day91pseudoneutid50_BA.1), 
                           c(Bpseudoneutid50_BA.1, Day15pseudoneutid50_BA.1)), 
                    x.ori = 0, xaxislabels = c("B", "D15"), cex.axis = 1, add = T, xlab = "", ylab = "", pcol = NULL, lcol = NULL)  

my.interaction.plot(subset(dat_proc, ph1.D15 & arm!=3 & Bpseudoneutid50_BA.1>1.5 & !is.na(Bpseudoneutid50_BA.1) & !is.na(Day15pseudoneutid50_BA.1) & !is.na(Day29pseudoneutid50_BA.1) & !is.na(Day91pseudoneutid50_BA.1), 
                           c(Day15pseudoneutid50_BA.1, Day29pseudoneutid50_BA.1)), 
                    x.ori = 1, xaxislabels = c("D15", "D29"), cex.axis = 1, add = T, xlab = "", ylab = "", pcol = NULL, lcol = NULL)  

my.interaction.plot(subset(dat_proc, ph1.D15 & arm!=3 & Bpseudoneutid50_BA.1>1.5 & !is.na(Bpseudoneutid50_BA.1) & !is.na(Day15pseudoneutid50_BA.1) & !is.na(Day29pseudoneutid50_BA.1) & !is.na(Day91pseudoneutid50_BA.1), 
                           c(Day29pseudoneutid50_BA.1, Day91pseudoneutid50_BA.1)), 
                    x.ori = 2, xaxislabels = c("D29", "D91"), cex.axis = 1, add = T, xlab = "", ylab = "", pcol = NULL, lcol = NULL)  



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

