library(survey)
library(kyotil)

dat_mapped=read.csv('/networks/cavd/Objective 4/GH-VAP/ID27-Sankaranarayanan/analysis/mapping_immune_correlates/adata/ID27_IARC_HPV_mapped_20240104.csv')
dat_proc=read.csv('/networks/cavd/Objective 4/GH-VAP/ID27-Sankaranarayanan/analysis/correlates/adata/id27hpv_data_processed_20240110.csv')

assay_metadata=read.csv('~/correlates_reporting2/assay_metadata/id27hpv_assay_metadata.csv')
assays=assay_metadata$assay
assays



################################################################################

# is single.dose confounder?
par(mfrow=c(1,3))
myboxplot(M18pseudoneutid50_HPV52~single.dose, dat_proc, test='w', main='ID50 HPV52', names=c('">1 dose"','1 dose'), col=ifelse(dat_proc$EventIndPrimaryAnyHPV==1,2,1), cex=1.2)
myboxplot(M18pseudoneutid50_HPV58~single.dose, dat_proc, test='w', main='ID50 HPV58', names=c('">1 dose"','1 dose'), col=ifelse(dat_proc$EventIndPrimaryAnyHPV==1,2,1), cex=1.2)
myboxplot(M18pseudoneutid50_mdw~single.dose, dat_proc, test='w', main='ID50 mdw', names=c('">1 dose"','1 dose'), col=ifelse(dat_proc$EventIndPrimaryAnyHPV==1,2,1), cex=1.2)

par(mfrow=c(1,1))
myboxplot(M18pseudoneutid50_mdw~single.dose, dat_proc, test='w', main='ID50 mdw', names=c('">1 dose"','1 dose'), col=ifelse(dat_proc$EventIndPrimaryAnyHPV==1,6,1), cex=1.3)

# some descriptives
with(subset(dat.ph1, ph2.M18==1), table(M18bind_HPV18cat, AgeGroup, EventIndPrimaryAnyHPV))
with(subset(dat.ph1, ph2.M18.sus==1), table(M18bind_HPV18cat, AgeGroup, EventIndPrimaryAnyHPV))
with(subset(dat.ph1, ph2.M18.sus==1), table(M18pseudoneutid50_HPV52cat, AgeGroup, EventIndPrimaryAnyHPV))



################################################################################
# marker missingness

summary(dat_mapped)

dat_mapped0=read.csv('/networks/cavd/Objective 4/GH-VAP/ID27-Sankaranarayanan/analysis/mapping_immune_correlates/adata/ID27_IARC_HPV_mapped_20231227.csv')
dat_mapped=read.csv('/networks/cavd/Objective 4/GH-VAP/ID27-Sankaranarayanan/analysis/mapping_immune_correlates/adata/ID27_IARC_HPV_mapped_20240103.csv')

subset(dat_mapped0, Ptid==25)
subset(dat_mapped, Ptid==25)

# correlation
mypairs(dat_mapped[, "M18"%.%assays[11:19]], ylim=c(1,5), xlim=c(1,5))

mypairs(dat_mapped[, "M18"%.%assays[1:9]], ylim=c(-1,2.5), xlim=c(-1,2.5))

mypairs(dat_mapped[, "M18"%.%assays[c(1,11,2,12,3,13,4,14,5,15)]])


for (i in 1:9) {
  print(table(!is.na(dat_mapped[["M18"%.%assays[i+10]]]), 
              !is.na(dat_mapped[["M18"%.%assays[i]]])))
}

dat=subset(dat_mapped, Perprotocol == 1 & EligibilityorinitialsamplingTimeM18>0)

table(!is.na(dat$M18pseudoneutid50_HPV11), dat$EventIndPrimaryAnyHPV)
table(!is.na(dat$M18bind_HPV11), dat$EventIndPrimaryAnyHPV)

corplot(M18pseudoneutid50_HPV11~M18bind_HPV11, dat, add.diagonal.line=F)
corplot(M18pseudoneutid50_HPV31~M18bind_HPV31, dat, add.diagonal.line=F)




# M18bind_HPV31 has 40% response rate and 0.5 correlation with other markers
# with(dat_proc, mypairs(cbind(M18bind_HPV6, M18bind_HPV11, M18bind_HPV16, M18bind_HPV18, M18bind_HPV31)))

with(dat_proc, summary(cbind(M18bind_HPV6, M18bind_HPV11, M18bind_HPV16, M18bind_HPV18, M18bind_HPV31)))
with(dat_proc, table(!is.na(M18bind_HPV6), !is.na(M18bind_HPV11)))


with(dat_proc, myboxplot(data.frame(M18pseudoneutid50_HPV6, M18pseudoneutid50_HPV11, M18pseudoneutid50_HPV16, M18pseudoneutid50_HPV18, M18pseudoneutid50_HPV31)))

with(dat_proc, sapply(data.frame(M18pseudoneutid50_HPV6, M18pseudoneutid50_HPV11, M18pseudoneutid50_HPV16, M18pseudoneutid50_HPV18, M18pseudoneutid50_HPV31), sd, na.rm=T))



myboxplot(dat_proc[,"M18"%.%assays[1:10]])
myboxplot(dat_proc[,"M18"%.%assays[11:19]])

min(dat_proc[,"M18"%.%assays[19]], na.rm=T)
with(dat_proc, sapply(dat_proc[,"M18"%.%assays[11:19]], function(x) mean(x>1.31, na.rm=T)))


with(subset(dat_proc, ph1.M18==1), table(ph2.M18, EventIndPrimaryAnyHPV))
with(subset(dat_proc, Perprotocol == 1 & EligibilityorinitialsamplingTimeM18>0), table(ph2.M18, EventIndPrimaryAnyHPV))


