library(survey)
library(kyotil)

dat_mapped=read.csv('/trials/covpn/p3003/analysis/mapping_immune_correlates/adata/COVID_ENSEMBLE_PartAComplete_variant_mapped_20240108.csv')
dat_proc = read.csv('/trials/covpn/p3003/analysis/correlates/Part_A_Blinded_Phase_Data/adata/janssen_partA_VL_data_processed_20240110.csv')
assay_metadata=read.csv('~/correlates_reporting2/assay_metadata/janssen_partA_VL_assay_metadata.csv')
assays=assay_metadata$assay

must_have_assays <- c("bindSpike", "bindRBD")


################################################################################
# missingness pattern

# an indicator for non-cases in the subcohort
select = with(dat_mapped, SubcohortInd & (is.na(EventIndPrimaryIncludeNotMolecConfirmedD29) | EventIndPrimaryIncludeNotMolecConfirmedD29 == 0))
# an indicator for non-cases in the whole cohort
select.ctrl = with(dat_mapped, (is.na(EventIndPrimaryIncludeNotMolecConfirmedD29) | EventIndPrimaryIncludeNotMolecConfirmedD29 == 0))

# bindRBD and bindSpike is all or none
table(!is.na(dat_mapped[select,"Day29bindSpike"]), !is.na(dat_mapped[select,"Day29bindRBD"]))


tmp = complete.cases(dat_mapped[,c("Day29"%.%must_have_assays)]) 

with(dat_mapped[select & dat_mapped$Region==0, ], table(!is.na(Day29bindSpike), !is.na(Day29bindSpike_D614)))
with(dat_mapped[select & dat_mapped$Region==2, ], table(!is.na(Day29bindSpike), !is.na(Day29bindSpike_D614)))
with(dat_mapped[select & dat_mapped$Region==1, ], table(!is.na(Day29bindSpike), !is.na(Day29bindSpike_D614)))




with(dat_mapped[select & dat_mapped$Region==0, ], table(!is.na(Day29bindSpike), !is.na(Day29pseudoneutid50)))
with(dat_mapped[select & dat_mapped$Region==2, ], table(!is.na(Day29bindSpike), !is.na(Day29pseudoneutid50_Beta)))
with(dat_mapped[select & dat_mapped$Region==1, ], table(!is.na(Day29bindSpike), (!is.na(Day29pseudoneutid50_Gamma) | !is.na(Day29pseudoneutid50_Lambda) | !is.na(Day29pseudoneutid50_Mu) | !is.na(Day29pseudoneutid50_Zeta))))

with(dat_mapped[select & dat_mapped$Region==1, ], table(!is.na(Day29pseudoneutid50_Gamma), !is.na(Day29pseudoneutid50_Lambda), !is.na(Day29pseudoneutid50_Mu), !is.na(Day29pseudoneutid50_Zeta)))

with(dat_mapped[select & dat_mapped$Region==1, ], table(!is.na(Day29pseudoneutid50_Gamma), !is.na(Day29pseudoneutid50_Zeta)))
with(dat_mapped[select & dat_mapped$Region==1, ], table(!is.na(Day29pseudoneutid50_Lambda), !is.na(Day29pseudoneutid50_Zeta)))
with(dat_mapped[select & dat_mapped$Region==1, ], table(!is.na(Day29pseudoneutid50_Mu), !is.na(Day29pseudoneutid50_Zeta)))

# no delta in region 1, but there is in 0 and 2


# missingness across MSD variants bAb in all
dat=dat_mapped
for (i in 4:13) {
  dat[['tmp'%.%i]]=ifelse(is.na(dat[["Day29"%.%assays[i]]]), 0, 1)
}
dat$tmp=with(dat,paste0(tmp4,tmp5,tmp6,tmp7,tmp8,tmp9,tmp10,tmp11,tmp12,tmp13))
table(dat$tmp)


# missingness across variants nAb
# among controls, all or nothing separately withinn {zeta, mu, gamma, lambda} and {Delta, Beta}
dat=dat_mapped[select.ctrl,]
for (i in 15:20) {
  dat[['tmp'%.%i]]=ifelse(is.na(dat[["Day29"%.%assays[i]]]), 0, 1)
}
dat$tmp=with(dat,paste0(tmp15,tmp16,tmp17,tmp18,tmp19,tmp20))
table(dat$tmp)


with(dat_mapped, table(!is.na(Day29bindSpike), !is.na(Day29bindSpike_D614), Region))
with(dat_mapped[dat_mapped$SubcohortInd==1,], table(!is.na(Day29bindSpike), !is.na(Day29bindSpike_D614), Region))

with(dat_mapped,                              table(!is.na(Day29pseudoneutid50), !is.na(Day29pseudoneutid50_Beta), Region))
with(dat_mapped[dat_mapped$SubcohortInd==1,], table(!is.na(Day29pseudoneutid50), !is.na(Day29pseudoneutid50_Beta), Region))


with(dat_mapped[select.ctrl,], table(!is.na(Day29pseudoneutid50), !is.na(Day29pseudoneutid50_Beta), Region))
with(dat_mapped[select ,], table(!is.na(Day29pseudoneutid50), !is.na(Day29pseudoneutid50_Beta), Region))


with(dat_mapped[select.ctrl,], table(!is.na(Day29bindSpike), !is.na(Day29bindSpike_D614), Region))
with(dat_mapped[select ,], table(!is.na(Day29bindSpike), !is.na(Day29bindSpike_D614), Region))


# distribution of variant bAb in controls by Region and subcohort
with(dat_mapped[select.ctrl & !is.na(dat_mapped$Day29bindSpike_D614),], table(Region, SubcohortInd))
with(dat_mapped, table(Region, SubcohortInd))

# distribution of variant bAb in cases by Region and subcohort
with(dat_mapped[!select.ctrl & !is.na(dat_mapped$Day29bindSpike), ,], table(Region, !is.na(Day29bindSpike_D614)))
with(dat_mapped[!select.ctrl , ], table(Region, !is.na(Day29bindSpike_D614)))


# distribution of Beta ID50 in controls by Region and subcohort
with(dat_mapped[select.ctrl & !is.na(dat_mapped$Day29pseudoneutid50_Beta),], table(Region, SubcohortInd))
# distribution of variant nAb in cases by Region and subcohort
with(dat_mapped[!select.ctrl & !is.na(dat_mapped$Day29pseudoneutid50), ,], table(Region, !is.na(Day29pseudoneutid50_Beta)))
with(dat_mapped[!select.ctrl , ], table(Region, !is.na(Day29pseudoneutid50_Beta)))


# distribution of Lambda ID50 in controls by Region and subcohort
with(dat_mapped[select.ctrl & !is.na(dat_mapped$Day29pseudoneutid50_Lambda),], table(Region, SubcohortInd))
# distribution of variant nAb in cases by Region and subcohort
with(dat_mapped[!select.ctrl & !is.na(dat_mapped$Day29pseudoneutid50), ,], table(Region, !is.na(Day29pseudoneutid50_Lambda)))
with(dat_mapped[!select.ctrl , ], table(Region, !is.na(Day29pseudoneutid50_Lambda)))
