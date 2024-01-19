{
library(survey)
library(kyotil)

dat_mapped=read.csv('/trials/covpn/p3003/analysis/mapping_immune_correlates/adata/COVID_ENSEMBLE_PartAComplete_variant_mapped_20240108.csv')
dat_proc = read.csv('/trials/covpn/p3003/analysis/correlates/Part_A_Blinded_Phase_Data/adata/janssen_partA_VL_data_processed_20240110.csv')
assay_metadata=read.csv('~/correlates_reporting2/assay_metadata/janssen_partA_VL_assay_metadata.csv')
assays=assay_metadata$assay

country.codes=c("USA", "ARG", "BRA", "CHL", "COL", "MEX", "PER", "ZAF")
dat_mapped$cc=country.codes[dat_mapped$Country+1]
dat_proc$cc=country.codes[dat_proc$Country+1]
}



################################################################################
# correlation between two ancestral ID50 markers

corplot(Day29bindSpike~Day29bindSpike_D614, dat_mapped)

corplot(Day29bindSpike~Day29bindSpike_D614, dat_mapped[dat_mapped$Day29bindSpike>1,])

mytable(!is.na(dat_mapped$Day29bindSpike), !is.na(dat_mapped$Day29bindSpike_D614))
mytable(!is.na(dat_mapped$Day29pseudoneutid50), !is.na(dat_mapped$Day29pseudoneutid50_Lambda))

dat=dat_mapped[,"Day29"%.%assays[c(1,4:13)]]
names(dat)=sub("Day29bindSpike","",names(dat))
mypairs(dat)

################################################################################
# missingness pattern

# an indicator for non-cases in the whole cohort
kp = with(dat_mapped, (is.na(EventIndPrimaryIncludeNotMolecConfirmedD29) | EventIndPrimaryIncludeNotMolecConfirmedD29 == 0))
# an indicator for cases in the whole cohort
kpc= with(dat_mapped, !(is.na(EventIndPrimaryIncludeNotMolecConfirmedD29) | EventIndPrimaryIncludeNotMolecConfirmedD29 == 0))
# an alias for subcohort
subc = dat_mapped$SubcohortInd==1


mytable(!is.na(dat_mapped$Day29bindSpike), !is.na(dat_mapped$Day29bindSpike_D614), kpc)


# bindRBD and bindSpike is all or none
table(!is.na(dat_mapped[,"Day29bindSpike"]), !is.na(dat_mapped[,"Day29bindRBD"]))

# 10% bAb ancestral from outside subhochort
table(!is.na(dat_mapped[kp,"Day29bindSpike"]), dat_mapped$SubcohortInd[kp])
table(!is.na(dat_mapped[kp,"Day29bindSpike"]), dat_mapped$SubcohortInd[kp], dat_mapped$Region[kp])

# half nAb ancestral from outside subcohort
table(!is.na(dat_mapped[kp,"Day29pseudoneutid50"]), dat_mapped$SubcohortInd[kp])
table(!is.na(dat_mapped[kp,"Day29pseudoneutid50"]), dat_mapped$SubcohortInd[kp], dat_mapped$Region[kp])



# are there any cases with variant Ab data but not D614/D614G Ab data?  
table(!is.na(dat_mapped[,"Day29bindSpike_D614"]), !is.na(dat_mapped[,"Day29bindSpike"]), subc)
table(!is.na(dat_mapped[,"Day29bindSpike_D614"]), !is.na(dat_mapped[,"Day29bindSpike"]), kpc)
table(!is.na(dat_mapped[kp & subc,"Day29bindSpike_D614"]), !is.na(dat_mapped[kp & subc,"Day29bindSpike"]))

table(!is.na(dat_mapped[,"Day29pseudoneutid50_Lambda"]), !is.na(dat_mapped[,"Day29pseudoneutid50"]), subc)
table(!is.na(dat_mapped[,"Day29pseudoneutid50_Lambda"]), !is.na(dat_mapped[,"Day29pseudoneutid50"]), kpc)
table(!is.na(dat_mapped[kp & subc,"Day29pseudoneutid50_Lambda"]), !is.na(dat_mapped[kp & subc,"Day29pseudoneutid50"]))


# variants are measured on the same set of ptids within region of interest
table(!is.na(dat_mapped[,"Day29pseudoneutid50_Lambda"])[kp & dat_mapped$Region==1], !is.na(dat_mapped[,"Day29bindSpike_D614"])[kp & dat_mapped$Region==1])
table(!is.na(dat_mapped[,"Day29pseudoneutid50_Beta"])[kp & dat_mapped$Region==2], !is.na(dat_mapped[,"Day29bindSpike_D614"])[kp & dat_mapped$Region==2])

with(dat_mapped[kp & dat_mapped$Region==1 & !is.na(dat_mapped$Day29pseudoneutid50_Lambda),], table(SubcohortInd,cc))

with(dat_mapped[kp & dat_mapped$Region==2 & !is.na(dat_mapped$Day29pseudoneutid50_Beta),], table(SubcohortInd,cc))



# in the subcohort, bindSpike and ID50 D614G are mostly all or none, spike has just 4 more
# but there are many more id50 outside subcohort
table(!is.na(dat_mapped[,"Day29bindSpike"]), !is.na(dat_mapped[,"Day29pseudoneutid50"]))
table(!is.na(dat_mapped[kp,"Day29bindSpike"]), !is.na(dat_mapped[kp,"Day29pseudoneutid50"]))
table(!is.na(dat_mapped[kp & SubcohortInd,"Day29bindSpike"]), !is.na(dat_mapped[kp & SubcohortInd,"Day29pseudoneutid50"]))


table(!is.na(dat_mapped[,"Day29pseudoneutid50_Lambda"]), !is.na(dat_mapped[,"Day29pseudoneutid50"]))
table(!is.na(dat_mapped[kp,"Day29pseudoneutid50_Lambda"]), !is.na(dat_mapped[kp,"Day29pseudoneutid50"]))
table(!is.na(dat_mapped[kp & dat_mapped$SubcohortInd,"Day29pseudoneutid50_Lambda"]), !is.na(dat_mapped[kp & dat_mapped$SubcohortInd,"Day29pseudoneutid50"]))


# 98 COL ptids sampled for variants study
dat=subset(dat_mapped, cc=="COL" & !subc & !is.na(Day29bindSpike_D614) & !kpc); nrow(dat)
summary(dat)


table(!is.na(dat_mapped[kp,"Day29bindSpike_D614"]), !is.na(dat_mapped[kp,"Day29bindSpike"]), subc[kp])

table(!is.na(dat_mapped[kp,"Day29pseudoneutid50_Lambda"]), !is.na(dat_mapped[kp,"Day29pseudoneutid50"]), subc[kp])


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


# missingness is all or none for the MSD panel bAb 
dat=dat_mapped
for (i in c(1,4:13)) {
  dat[['tmp'%.%i]]=ifelse(is.na(dat[["Day29"%.%assays[i]]]), 0, 1)
}
dat$tmp=with(dat,paste0(tmp1, tmp4,tmp5,tmp6,tmp7,tmp8,tmp9,tmp10,tmp11,tmp12,tmp13))
table(dat$tmp[dat$Region==1], kpc[dat$Region==1])


# missingness across variants nAb
# among controls, all or nothing separately withinn {zeta, mu, gamma, lambda} and {Delta, Beta}
dat=dat_mapped#[kp,]
for (i in 14:20) {
  dat[['tmp'%.%i]]=ifelse(is.na(dat[["Day29"%.%assays[i]]]), 0, 1)
}
dat$tmp=with(dat,paste0(tmp14,tmp15,tmp16,tmp17,tmp18,tmp19,tmp20))
table(dat$tmp[dat$Region==1], kpc[dat$Region==1])


with(dat_mapped, table(!is.na(Day29bindSpike), !is.na(Day29bindSpike_D614), Region))
with(dat_mapped[dat_mapped$SubcohortInd==1,], table(!is.na(Day29bindSpike), !is.na(Day29bindSpike_D614), Region))

with(dat_mapped,                              table(!is.na(Day29pseudoneutid50), !is.na(Day29pseudoneutid50_Beta), Region))
with(dat_mapped[dat_mapped$SubcohortInd==1,], table(!is.na(Day29pseudoneutid50), !is.na(Day29pseudoneutid50_Beta), Region))


with(dat_mapped[kp,], table(!is.na(Day29pseudoneutid50), !is.na(Day29pseudoneutid50_Beta), Region))
with(dat_mapped[select ,], table(!is.na(Day29pseudoneutid50), !is.na(Day29pseudoneutid50_Beta), Region))


with(dat_mapped[kp,], table(!is.na(Day29bindSpike), !is.na(Day29bindSpike_D614), Region))
with(dat_mapped[select ,], table(!is.na(Day29bindSpike), !is.na(Day29bindSpike_D614), Region))


# distribution of variant bAb in controls by Region and subcohort
with(dat_mapped[kp & !is.na(dat_mapped$Day29bindSpike_D614),], table(Region, SubcohortInd))
with(dat_mapped, table(Region, SubcohortInd))

# distribution of variant bAb in cases by Region and subcohort
with(dat_mapped[!kp & !is.na(dat_mapped$Day29bindSpike), ,], table(Region, !is.na(Day29bindSpike_D614)))
with(dat_mapped[!kp , ], table(Region, !is.na(Day29bindSpike_D614)))


# distribution of Beta ID50 in controls by Region and subcohort
with(dat_mapped[kp & !is.na(dat_mapped$Day29pseudoneutid50_Beta),], table(Region, SubcohortInd))
# distribution of variant nAb in cases by Region and subcohort
with(dat_mapped[!kp & !is.na(dat_mapped$Day29pseudoneutid50), ,], table(Region, !is.na(Day29pseudoneutid50_Beta)))
with(dat_mapped[!kp , ], table(Region, !is.na(Day29pseudoneutid50_Beta)))


# distribution of Lambda ID50 in controls by Region and subcohort
with(dat_mapped[kp & !is.na(dat_mapped$Day29pseudoneutid50_Lambda),], table(Region, SubcohortInd))
# distribution of variant nAb in cases by Region and subcohort
with(dat_mapped[!kp & !is.na(dat_mapped$Day29pseudoneutid50), ,], table(Region, !is.na(Day29pseudoneutid50_Lambda)))
with(dat_mapped[!kp , ], table(Region, !is.na(Day29pseudoneutid50_Lambda)))
