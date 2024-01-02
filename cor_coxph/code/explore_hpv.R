library(survey)
library(kyotil)

dat_mapped=read.csv('/networks/cavd/Objective 4/GH-VAP/ID27-Sankaranarayanan/analysis/mapping_immune_correlates/adata/ID27_IARC_HPV_mapped_20231227.csv')
dat_proc=read.csv('/networks/cavd/Objective 4/GH-VAP/ID27-Sankaranarayanan/analysis/correlates/adata/id27hpv_data_processed_20231228.csv')

assay_metadata=read.csv('~/correlates_reporting2/assay_metadata/id27hpv_assay_metadata.csv')
assays=assay_metadata$assay
assays


################################################################################
# marker missingness
summary(dat_mapped)


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
