library(survey)
library(kyotil)

dat_mapped=read.csv('/trials/covpn/COVAILcorrelates/analysis/mapping_immune_correlates/adata/covail_mapped_data_12222023.csv')

assay_metadata=read.csv('~/correlates_reporting2/assay_metadata/covail_assay_metadata.csv')
assays=assay_metadata$assay
assays

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

