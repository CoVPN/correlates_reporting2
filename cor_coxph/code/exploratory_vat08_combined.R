library(survey)
library(kyotil)

# "Colombia" = 1, "Ghana" = 2, "Honduras" = 3, "India" = 4, "Japan" = 5, "Kenya" = 6, "Nepal" = 7, "United States" = 8, "Mexico" = 9, "Uganda" = 10, "Ukraine" = 11

dat_mapped=read.csv('/trials/covpn/p3005/analysis/mapping_immune_correlates/combined/adata/COVID_Sanofi_stage1and2_mapped_20231228.csv')
dat_proc = read.csv('/trials/covpn/p3005/analysis/correlates/Part_A_Blinded_Phase_Data/adata/vat08_combined_data_processed_20231228.csv')
assay_metadata=read.csv('~/correlates_reporting2/assay_metadata/vat08_combined_assay_metadata.csv')
assays=assay_metadata$assay

dat_mapped$EventIndPrimaryD1 = dat_mapped$EventIndFirstInfectionD1
dat_mapped$EventTimePrimaryD1 = dat_mapped$EventTimeFirstInfectionD1
dat_mapped$EventTimePrimaryD43= dat_mapped$EventTimeFirstInfectionD43
dat_mapped$EarlyendpointD43 <- with(dat_mapped, EarlyinfectionD43==1 | (EventIndPrimaryD1==1 & EventTimePrimaryD1 < NumberdaysD1toD43 + 7),1,0)
dat_mapped$ph1.D43= with(dat_mapped, EarlyendpointD43==0 & Perprotocol==1 & EventTimePrimaryD43 >= 7)



################################################################################
# select variables to adjust using ph1 data 

dat=subset(dat_proc, Trialstage==2 & ph1.D43) 

f=Surv(EventTimeOmicronD43M6hotdeck1, EventIndOmicronD43M6hotdeck1) ~ Bserostatus+Trt+strata(Country)

fit.1=coxph(update(f, ~.+ standardized_risk_score + FOI + Sex), dat)
fit.2=coxph(update(f, ~.+ FOI + Sex), dat)
fit.3=coxph(update(f, ~.+ standardized_risk_score + Sex), dat)
fit.4=coxph(update(f, ~.+ standardized_risk_score + FOI), dat)

summary(fit.1)
anova(fit.1, fit.2)
anova(fit.1, fit.3)
anova(fit.1, fit.4)

fit.2.1=coxph(update(f, ~.+ Sex), dat)
fit.2.2=coxph(update(f, ~.+ FOI), dat)

anova(fit.2, fit.2.1)
anova(fit.2, fit.2.2)


dat=subset(dat_proc, Trialstage==2 & Trt==0 & ph1.D43) 

fit.1=coxph(Surv(EventTimeOmicronD43M6hotdeck1, EventIndOmicronD43M6hotdeck1) ~ standardized_risk_score+FOI+Sex+Bserostatus+strata(Country), dat)
fit.4=coxph(Surv(EventTimeOmicronD43M6hotdeck1, EventIndOmicronD43M6hotdeck1) ~ Sex+Bserostatus+strata(Country), dat)
fit.5=coxph(Surv(EventTimeOmicronD43M6hotdeck1, EventIndOmicronD43M6hotdeck1) ~ Bserostatus+strata(Country), dat)

summary(fit.1)
anova(fit.1, fit.4)
anova(fit.1, fit.5)


dat=subset(dat_proc, Trialstage==2 & Trt==1 & ph1.D43) 

fit.1=coxph(Surv(EventTimeOmicronD43M6hotdeck1, EventIndOmicronD43M6hotdeck1) ~ standardized_risk_score+FOI+Sex+Bserostatus+strata(Country), dat)
fit.4=coxph(Surv(EventTimeOmicronD43M6hotdeck1, EventIndOmicronD43M6hotdeck1) ~ Sex+Bserostatus+strata(Country), dat)
fit.5=coxph(Surv(EventTimeOmicronD43M6hotdeck1, EventIndOmicronD43M6hotdeck1) ~ Bserostatus+strata(Country), dat)

summary(fit.1)
anova(fit.1, fit.4)
anova(fit.1, fit.5)



dat=subset(dat_proc, Trialstage==1 & Trt==1 & ph1.D43 & Bserostatus==1) 

fit.1=coxph(Surv(EventTimeOmicronD43M6hotdeck1, EventIndOmicronD43M6hotdeck1) ~ standardized_risk_score+FOI+Sex+strata(Country), dat)
fit.4=coxph(Surv(EventTimeOmicronD43M6hotdeck1, EventIndOmicronD43M6hotdeck1) ~ Sex+strata(Country), dat)
fit.5=coxph(Surv(EventTimeOmicronD43M6hotdeck1, EventIndOmicronD43M6hotdeck1) ~ strata(Country), dat)

summary(fit.1)
anova(fit.1, fit.4)
anova(fit.1, fit.5)


summary(coxph(Surv(EventTimeOmicronD43M6hotdeck1, EventIndOmicronD43M6hotdeck1) ~ standardized_risk_score+FOI+Sex, dat))


################################################################################
# Senior and nAbBatch
# suggests that batch 2 is complementary, almost no senior in batch 2

dat=subset(dat_mapped, Trt==1 & Trialstage==2 & Bserostatus==1)  
dat$ph2=with(dat, !is.na(Bpseudoneutid50) & !is.na(Day43pseudoneutid50))
with(dat[dat$ph2,], table(Age>=60, nAbBatch, EventIndPrimaryD1))

# using bAb gives same conclusion
dat$ph2=with(dat, !is.na(BbindSpike) & !is.na(Day43bindSpike))
with(dat[dat$ph2,], table(Age>=60, nAbBatch, EventIndPrimaryD1))


with(subset(dat_mapped, Trt==1 & Trialstage==2 & Bserostatus==1), 
                 table(!is.na(Day43pseudoneutid50),
                       !is.na(Day43bindSpike),
                       nAbBatch))

with(subset(dat_mapped, Trt==1 & Trialstage==2 & Bserostatus==1), table(!is.na(Day43pseudoneutid50)))
with(subset(dat_mapped, Trt==1 & Trialstage==2 & Bserostatus==1), table(!is.na(Day43bindSpike)))



with(subset(dat_proc, Trialstage==2 & Bserostatus==1 & Trt==1 & ph1.D43 & !is.na(Bpseudoneutid50) & !is.na(Day43pseudoneutid50)), table(nAbBatch, Country))

################################################################################
# check missingness pattern 

# ID50

with(subset(dat_proc, Trialstage==2 & Bserostatus==1 & Trt==1 & ph1.D43), 
     table(!is.na(Bpseudoneutid50), !is.na(Day43pseudoneutid50)))

with(subset(dat_proc, Trialstage==2 & Bserostatus==1 & Trt==1 & ph1.D43 & !EventIndOmicronD43M12hotdeck1), 
     table(!is.na(Bpseudoneutid50), !is.na(Day43pseudoneutid50)))

with(subset(dat_proc, Trialstage==2 & Bserostatus==1 & Trt==1 & ph1.D43 & EventIndOmicronD43M12hotdeck1), 
     table(!is.na(Bpseudoneutid50), !is.na(Day43pseudoneutid50)))

with(subset(dat_proc, Trialstage==2 & Bserostatus==1 & Trt==1 & ph1.D43), 
     table(!is.na(Bpseudoneutid50), !is.na(Day43pseudoneutid50), nAbBatch))

# plot fold change
par(mfrow=c(2,4))
  corplot(Day43pseudoneutid50~Bpseudoneutid50, subset(dat_proc, Trialstage==2 & Bserostatus==1 & Trt==1 & ph1.D43), ylim=c(0,4))
  myboxplot(Day43pseudoneutid50~I(Bpseudoneutid50>0.2), subset(dat_proc, Trialstage==2 & Bserostatus==1 & Trt==1 & ph1.D43), names=c("Undetectable at B","Detectable at B"), ylab="D43", ylim=c(0,4))
  my.interaction.plot(subset(dat_proc, Trialstage==2 & Bserostatus==1 & Trt==1 & ph1.D43 & !is.na(Bpseudoneutid50) & !is.na(Day43pseudoneutid50), c(Bpseudoneutid50, Day43pseudoneutid50)), 
                      x.ori = 0, xaxislabels = c("B", "D43"), cex.axis = 1, add = FALSE, xlab = "", ylab = "", pcol = NULL, lcol = NULL)  
  corplot(Delta43overBpseudoneutid50~Bpseudoneutid50, subset(dat_proc, Trialstage==2 & Bserostatus==1 & Trt==1 & ph1.D43))
  
  # BA.1
  corplot(Day43pseudoneutid50_BA.1~Bpseudoneutid50_BA.1, subset(dat_proc, Trialstage==2 & Bserostatus==1 & Trt==1 & ph1.D43), ylim=c(0,4))
  myboxplot(Day43pseudoneutid50_BA.1~I(Bpseudoneutid50_BA.1>0.2), subset(dat_proc, Trialstage==2 & Bserostatus==1 & Trt==1 & ph1.D43), names=c("Undetectable at B","Detectable at B"), ylab="D43", ylim=c(0,4))
  my.interaction.plot(subset(dat_proc, Trialstage==2 & Bserostatus==1 & Trt==1 & ph1.D43 & !is.na(Bpseudoneutid50_BA.1) & !is.na(Day43pseudoneutid50_BA.1), c(Bpseudoneutid50_BA.1, Day43pseudoneutid50_BA.1)), 
                      x.ori = 0, xaxislabels = c("B", "D43"), cex.axis = 1, add = FALSE, xlab = "", ylab = "", pcol = NULL, lcol = NULL)  
  corplot(Delta43overBpseudoneutid50_BA.1~Bpseudoneutid50_BA.1, subset(dat_proc, Trialstage==2 & Bserostatus==1 & Trt==1 & ph1.D43))
  

# the conclusion from the above is that we need to restrict to samples with both B and D43 markers
  
  
################################################################################
# coxph M6

dat_proc$ph2=with(dat_proc, !is.na(Bpseudoneutid50) & !is.na(Day43pseudoneutid50) )#& nAbBatch==2
dat_proc$Wstratum=dat_proc$Wstratum.nAb # Wstratum.original
dat=subset(dat_proc, Trialstage==2 & Trt==1 & Bserostatus==1 & ph1.D43)   # !Senior
with(dat, table(ph2, Wstratum))

# select variables to adjust using ph1 data 
coxph(Surv(EventTimeOmicronD43M6hotdeck1, EventIndOmicronD43M6hotdeck1) ~ standardized_risk_score, dat)
coxph(Surv(EventTimeOmicronD43M6hotdeck1, EventIndOmicronD43M6hotdeck1) ~ FOI, dat)
coxph(Surv(EventTimeOmicronD43M6hotdeck1, EventIndOmicronD43M6hotdeck1) ~ standardized_risk_score+Sex, dat)
coxph(Surv(EventTimeOmicronD43M6hotdeck1, EventIndOmicronD43M6hotdeck1) ~ FOI+Sex, dat)
coxph(Surv(EventTimeOmicronD43M6hotdeck1, EventIndOmicronD43M6hotdeck1) ~ standardized_risk_score+Sex+strata(Country), dat)
coxph(Surv(EventTimeOmicronD43M6hotdeck1, EventIndOmicronD43M6hotdeck1) ~ FOI+Sex+strata(Country), dat)

fit.1=coxph(Surv(EventTimeOmicronD43M6hotdeck1, EventIndOmicronD43M6hotdeck1) ~ standardized_risk_score+FOI+Sex+strata(Country), dat)
fit.4=coxph(Surv(EventTimeOmicronD43M6hotdeck1, EventIndOmicronD43M6hotdeck1) ~ Sex+strata(Country), dat)
fit.5=coxph(Surv(EventTimeOmicronD43M6hotdeck1, EventIndOmicronD43M6hotdeck1) ~ strata(Country), dat)
fit.6=coxph(Surv(EventTimeOmicronD43M6hotdeck1, EventIndOmicronD43M6hotdeck1) ~ Sex, dat)

summary(fit.1)
anova(fit.1, fit.4)
anova(fit.1, fit.5)
anova(fit.1, fit.6)

# based on the above, we choose
f.base = Surv(EventTimeOmicronD43M6hotdeck1, EventIndOmicronD43M6hotdeck1) ~ Sex + strata(Country)

# Baseline detectable or not 
svycoxph(update(f.base, ~. + I(Bpseudoneutid50>0.116)), twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dat) )
# Baseline
svycoxph(update(f.base, ~. + Bpseudoneutid50), twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dat) )
# D43
svycoxph(update(f.base, ~. + Day43pseudoneutid50), twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dat) )
# Baseline detectable or not * fold change
svycoxph(update(f.base, ~. + I(Bpseudoneutid50>0.116) * scale(Delta43overBpseudoneutid50)), twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dat) )
# Baseline * fold change
svycoxph(update(f.base, ~. + scale(Bpseudoneutid50)*scale(Delta43overBpseudoneutid50)), twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dat) )


# BA.1
# Baseline detectable or not 
svycoxph(update(f.base, ~. + I(Bpseudoneutid50_BA.1>0.116)), twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dat) )
# Baseline
svycoxph(update(f.base, ~. + Bpseudoneutid50_BA.1), twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dat) )
# D43
svycoxph(update(f.base, ~. + Day43pseudoneutid50_BA.1), twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dat) )
# Baseline detectable or not * fold change
svycoxph(update(f.base, ~. + I(Bpseudoneutid50_BA.1>0.116) * scale(Delta43overBpseudoneutid50_BA.1)), twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dat) )
# Baseline * fold change
svycoxph(update(f.base, ~. + scale(Bpseudoneutid50_BA.1)*scale(Delta43overBpseudoneutid50_BA.1)), twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dat) )


# other markers
svycoxph(update(f.base, ~. + Bpseudoneutid50_B.1.351*Delta43overBpseudoneutid50_B.1.351),
         twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dat) )
svycoxph(update(f.base, ~. + Bpseudoneutid50_BA.1*Delta43overBpseudoneutid50_BA.1),
         twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dat) )
svycoxph(update(f.base, ~. + Bpseudoneutid50_BA.2*Delta43overBpseudoneutid50_BA.2),
         twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dat) )
svycoxph(update(f.base, ~. + Bpseudoneutid50_BA.4.5*Delta43overBpseudoneutid50_BA.4.5),
         twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dat) )
svycoxph(update(f.base, ~. + Bpseudoneutid50_mdw*Delta43overBpseudoneutid50_mdw),
         twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dat) )

svycoxph(update(f.base, ~. + I(Bpseudoneutid50_B.1.351>0.116)*Delta43overBpseudoneutid50_B.1.351),
         twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dat) )
svycoxph(update(f.base, ~. + I(Bpseudoneutid50_BA.1>0.116)*Delta43overBpseudoneutid50_BA.1),
         twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dat) )
svycoxph(update(f.base, ~. + I(Bpseudoneutid50_BA.2>0.116)*Delta43overBpseudoneutid50_BA.2),
         twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dat) )
svycoxph(update(f.base, ~. + I(Bpseudoneutid50_BA.4.5>0.116)*Delta43overBpseudoneutid50_BA.4.5),
         twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dat) )


################################################################################
# batch 2

dat_proc$ph2=with(dat_proc, !is.na(Bpseudoneutid50) & !is.na(Day43pseudoneutid50) & nAbBatch==2)
dat_proc$Wstratum=dat_proc$Wstratum.original
dat=subset(dat_proc, Trialstage==2 & Trt==1 & Bserostatus==1 & ph1.D43 & !Senior)
with(dat, table(ph2, Wstratum))

# based on the above, we choose
f.base = Surv(EventTimeOmicronD43M6hotdeck1, EventIndOmicronD43M6hotdeck1) ~ Sex + strata(Country)

# Baseline detectable or not 
svycoxph(update(f.base, ~. + I(Bpseudoneutid50>0.116)), twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dat) )
# Baseline
svycoxph(update(f.base, ~. + Bpseudoneutid50), twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dat) )
# D43
svycoxph(update(f.base, ~. + Day43pseudoneutid50), twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dat) )
# Baseline detectable or not * fold change
svycoxph(update(f.base, ~. + I(Bpseudoneutid50>0.116) * scale(Delta43overBpseudoneutid50)), twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dat) )
# Baseline * fold change
svycoxph(update(f.base, ~. + scale(Bpseudoneutid50)*scale(Delta43overBpseudoneutid50)), twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dat) )



################################################################################
# bAb

with(subset(dat_proc, Trialstage==2 & Bserostatus==1 & Trt==1 & ph1.D43), 
     table(!is.na(BbindSpike), !is.na(Day43bindSpike)))

with(subset(dat_proc, Trialstage==2 & Bserostatus==1 & Trt==1 & ph1.D43 & !EventIndOmicronD43M12hotdeck1), 
     table(!is.na(BbindSpike), !is.na(Day43bindSpike)))

with(subset(dat_proc, Trialstage==2 & Bserostatus==1 & Trt==1 & ph1.D43 & EventIndOmicronD43M12hotdeck1), 
     table(!is.na(BbindSpike), !is.na(Day43bindSpike)))

with(subset(dat_proc, Trialstage==2 & Bserostatus==1 & Trt==1 & ph1.D43), 
     table(!is.na(BbindSpike), !is.na(Day43bindSpike), nAbBatch))

par(mfrow=c(2,2))
corplot(Day43bindSpike~BbindSpike, subset(dat_proc, Trialstage==2 & Bserostatus==1 & Trt==1 & ph1.D43))
myboxplot(Day43bindSpike~I(BbindSpike>0.5), subset(dat_proc, Trialstage==2 & Bserostatus==1 & Trt==1 & ph1.D43), names=c("Undetectable at B","Detectable at B"))
my.interaction.plot(subset(dat_proc, Trialstage==2 & Bserostatus==1 & Trt==1 & ph1.D43 & !is.na(BbindSpike) & !is.na(Day43bindSpike), c(BbindSpike, Day43bindSpike)), 
                    x.ori = 0, xaxislabels = c("B", "D43"), cex.axis = 1, add = FALSE, xlab = "", ylab = "", pcol = NULL, lcol = NULL)  
corplot(Delta43overBbindSpike~BbindSpike, subset(dat_proc, Trialstage==2 & Bserostatus==1 & Trt==1 & ph1.D43))


# the conclusion from the above is that we need to restrict to samples with both B and D43 markers



dat_proc$ph2=with(dat_proc, !is.na(BbindSpike) & !is.na(Day43bindSpike))
dat_proc$Wstratum=dat_proc$Wstratum.original
dat=subset(dat_proc, Trt==1 & Trialstage==2 & Bserostatus==1 & ph1.D43)  

svycoxph(Surv(EventTimeOmicronD43M12hotdeck1, EventIndOmicronD43M12hotdeck1) ~ Sex + Day43bindSpike + strata(Country),
         twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dat) )

svycoxph(Surv(EventTimeOmicronD43M12hotdeck1, EventIndOmicronD43M12hotdeck1) ~ Sex + BbindSpike + strata(Country),
         twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dat) )

svycoxph(Surv(EventTimeOmicronD43M12hotdeck1, EventIndOmicronD43M12hotdeck1) ~ Sex + BbindSpike + Delta43overBbindSpike + strata(Country),
         twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dat) )

svycoxph(Surv(EventTimeOmicronD43M12hotdeck1, EventIndOmicronD43M12hotdeck1) ~ Sex + BbindSpike*Delta43overBbindSpike + strata(Country),
         twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dat) )




with(subset(dat_mapped, Trialstage==2 & Bserostatus==1 & Trt==1 & ph1.D43 & !EventIndOmicronD43M12hotdeck1), 
     table(!is.na(Day43bindSpike) & !is.na(BbindSpike), !is.na(Day43pseudoneutid50) & !is.na(Bpseudoneutid50), nAbBatch))

with(subset(dat_mapped, Trialstage==2 & Bserostatus==1 & Trt==1 & ph1.D43 & !EventIndOmicronD43M12hotdeck1), 
     table(ph2.D43.bAb, ph2.D43.nAb, nAbBatch))

with(subset(dat_mapped, Trialstage==2 & Bserostatus==1 & Trt==1 & ph1.D43 & !EventIndOmicronD43M12hotdeck1), 
     table(!is.na(Day43pseudoneutid50) & !is.na(Bpseudoneutid50), ph2.D43.nAb, nAbBatch))

subset(dat_mapped, Trialstage==2 & Bserostatus==1 & Trt==1 & ph1.D43 & !EventIndOmicronD43M12hotdeck1 & 
         !is.na(Day43pseudoneutid50) & !is.na(Bpseudoneutid50) & !ph2.D43.nAb)[1:2,]

with(subset(dat_mapped, Trialstage==2 & Bserostatus==1 & Trt==1 & ph1.D43 & EventIndOmicronD43M12hotdeck1), 
     table(!is.na(Day43bindSpike), !is.na(Day43pseudoneutid50)))


with(subset(dat_mapped, Trialstage==2 & Bserostatus==1 & Trt==1 & ph1.D43 & EventIndOmicronD43M12hotdeck1), 
     table(ph2.D43.nAb, !is.na(Day43pseudoneutid50)))



################################################################################

# some numbers for Peter's slide
with(subset(dat_proc, ph2.D22.nAb==1), table(1-Trt, EventIndOmicronD43M6hotdeck1, Bserostatus, Trialstage))

with(subset(dat_proc, ph2.D22.bAb==1), table(1-Trt, EventIndOmicronD43M6hotdeck1, Bserostatus, Trialstage))

with(subset(dat_proc, ph2.D22.nAb==1 & (EventIndOmicronD43M6hotdeck1 | nAbBatch==1)), table(Trt, EventIndOmicronD43M6hotdeck1, Bserostatus, Trialstage))




################################################################################
# batch

# both batches refer to stage 2 BV trial
with(subset(dat_mapped), table(nAbBatch, Trialstage))

# batch 2 has a lot more naive vaccinees
with(subset(dat_mapped, Trialstage==2), table(Trt, Bserostatus, nAbBatch))

with(subset(dat_mapped, Trialstage==2 & nAbBatch==1 & EventIndFirstInfectionD1!=1), table(Age>=60, Bserostatus, Country))

with(subset(dat_mapped, Trialstage==2 & nAbBatch==2 & EventIndFirstInfectionD1!=1), table(Age>=60, Bserostatus, Country))
with(subset(dat_mapped, Trialstage==2 & nAbBatch==2 & EventIndFirstInfectionD1!=1), table(Trt, Age>=60, Bserostatus))


with(subset(dat_mapped, Trialstage==2 & EventIndFirstInfectionD1!=1), table(nAbBatch, SubcohortInd, useNA='ifany'))
with(subset(dat_mapped, Trialstage==2 & EventIndFirstInfectionD1==1), table(nAbBatch, useNA='ifany'))

with(subset(dat_proc, Trialstage==2), table(nAbBatch, SubcohortInd, useNA='ifany'))

subset(dat_proc, Ptid%in%c('VAT00008-170-0002-20006','VAT00008-170-0002-20040','VAT00008-170-0002-21001'), SubcohortInd)

# batch 1 is subcohortind==1

with(subset(dat_proc, Trialstage==2 & EventIndOmicronD43M12hotdeck1==1), table(nAbBatch, useNA='ifany'))
with(subset(dat_proc, Trialstage==2 & EventIndOmicronD43M6hotdeck1==1), table(nAbBatch, useNA='ifany'))
with(subset(dat_proc, Trialstage==2 & EventIndOmicronD43M12hotdeck1==1 & Bserostatus==1 & ph1.D43 & Trt==1), table(nAbBatch, useNA='ifany'))

with(subset(dat_proc, Trialstage==2 & EventIndOmicronD1M12hotdeck1==1 ), table(nAbBatch, useNA='ifany'))



# # verified: no new nAb or bAb data from 1 to 2
# dat1=read.csv('/trials/covpn/p3005/analysis/mapping_immune_correlates/combined/adata/COVID_Sanofi_stage1and2_mapped_20231030.csv')
# dat2=read.csv('/trials/covpn/p3005/analysis/mapping_immune_correlates/combined/adata/COVID_Sanofi_stage1and2_mapped_20231107.csv')


country.codes=c("Colombia", "Ghana", "Honduras", "India", "Japan", "Kenya", "Nepal", "United States", "Mexico", "Uganda", "Ukraine")



################################################################################
# compare ID50 between step 1 and 2

dat=subset(dat_proc, Trialstage==2 & Trt==1 & Bserostatus==1 & ph1.D43 & !EventIndOmicronD43M6hotdeck1)
dat.1=subset(dat_proc, Trialstage==2 & Trt==1 & Bserostatus==1 & ph1.D43 & EventIndOmicronD43M6hotdeck1)

myboxplot(list(dat$Day43pseudoneutid50[dat$nAbBatch==1], 
               dat$Day43pseudoneutid50[dat$nAbBatch==2], 
               dat.1$Day43pseudoneutid50),
          names=c("Batch 1 non-cases", "Batch 2 non-cases", 'cases'))

summary(subset(dat, SubcohortInd==1, Day43pseudoneutid50, drop=T))
summary(subset(dat, SubcohortInd==0, Day43pseudoneutid50, drop=T))

summary(subset(dat, nAbBatch==1, Day43pseudoneutid50, drop=T))
summary(subset(dat, nAbBatch==2, Day43pseudoneutid50, drop=T))

with(dat, table(nAbBatch, !is.na(Day43pseudoneutid50), useNA='ifany'))

# by country
par(mfrow=c(1,1))
myboxplot(Day43pseudoneutid50~nAbBatch+Senior+Country, dat)
abline(v=c(4.5,8.5,11.5,14.5,18.5),lty=2)
# myboxplot(Day43pseudoneutid50~nAbBatch, dat, test='w', names=c('batch 1','batch 2'))

for (c in c(1,2,6,7,9,10)) {
  cat(c)
  print(wilcox.test(Day43pseudoneutid50~SubcohortInd, subset(dat,Country==c))$p.value)
}


################################################################################
# For Stage 2 NN, a scatterplot, with one point for each country, where the x-axis is the D43 GM ID50 and the y-axis is the country’s COVID-19 rate.
# to investigate confounding by country.  And each point with text for ncase:ncontrols.  

# create a table with one row for each country
tab = cbind(
  aggregate(EventIndOmicronD43M6hotdeck1~Country, subset(dat.mock, Trialstage==2 & Trt==1 & Bserostatus==1 & ph1.D43), mean),
  id50 = 
  aggregate(Day43pseudoneutid50*wt.D43.nAb~Country, subset(dat.mock, Trialstage==2 & Trt==1 & Bserostatus==1 & ph1.D43), sum)[,2]/
  aggregate(wt.D43.nAb~Country, subset(dat.mock, Trialstage==2 & Trt==1 & Bserostatus==1 & ph1.D43), sum)[,2],
  ncases=aggregate(EventIndOmicronD43M6hotdeck1~Country, subset(dat.mock, Trialstage==2 & Trt==1 & Bserostatus==1 & ph1.D43), sum)[,2],
  ncases.ph2=aggregate(EventIndOmicronD43M6hotdeck1~Country, subset(dat.mock, Trialstage==2 & Trt==1 & Bserostatus==1 & ph2.D43.nAb), sum)[,2],
  ncontrols=aggregate(!EventIndOmicronD43M6hotdeck1~Country, subset(dat.mock, Trialstage==2 & Trt==1 & Bserostatus==1 & ph1.D43), sum)[,2],
  ncontrols.ph2=aggregate(!EventIndOmicronD43M6hotdeck1~Country, subset(dat.mock, Trialstage==2 & Trt==1 & Bserostatus==1 & ph2.D43.nAb), sum)[,2]
)
names(tab)[2]='covid.rate'
let=ifelse(tab$Country==10,rep('A',nrow(tab)),tab$Country%.%"")
with(tab, plot(covid.rate, id50, pch=let))
mylegend(x=8, legend=country.codes[tab[,1]], pch=let%.%"")
tab



################################################################################
# auc for risk score after adding country

dat=subset(dat.mock, Trialstage==2 & Trt==1 & Bserostatus==1 & ph1.D43)

trainauc.glm(glm(EventIndOmicronD43M6hotdeck1~as.factor(Country)+Sex, dat, family=binomial))

trainauc.glm(glm(EventIndOmicronD43M6hotdeck1~as.factor(Country)+Sex+standardized_risk_score, dat, family=binomial))

trainauc.glm(glm(EventIndOmicronD43M6hotdeck1~Sex+standardized_risk_score, dat, family=binomial))

trainauc.glm(glm(EventIndOmicronD43M6hotdeck1~standardized_risk_score, dat, family=binomial))
trainauc.glm(glm(EventIndOmicronD43M6hotdeck1~as.factor(Country), dat, family=binomial))
trainauc.glm(glm(EventIndOmicronD43M6hotdeck1~as.factor(Country)+standardized_risk_score, dat, family=binomial))

trainauc.glm(glm(EventIndOmicronD43M6hotdeck1~as.factor(Country)+Sex+FOI, dat, family=binomial))


################################################################################
# examine which countries have cases

# E.g., For the non-naïve cohort in a given trial, Exclude if either 
# (1) zero eligible per-protocol Omicron COVID-19 endpoints >= 7 days post D43 through 180 days post dose 2 (pooling over vaccine, placebo) 
#     with D01 and D43 nAb, using the simplified definition of Omicron COVID-19 that assumes Omicron if an unknown lineage and event after January 17, 2022; 

# filter out country 2, 4 for stage 2 naive, filter out 4 for stage 2 nnaive
with(subset(dat_proc, Trialstage==2 & ph1.D43), table(EventIndOmicronD43M12hotdeck1, Country, Bserostatus))
with(subset(dat_proc, Trialstage==2 & ph2.D43.nAb), table(EventIndOmicronD43M6hotdeck1, Country, Bserostatus))

with(subset(dat_mapped, Trialstage==2 & ph1.D43), table(EventIndPrimaryD1, Country, Bserostatus)) # 2 in Ind, both Delta


# filter out country 4, 5, 6 for stage 1 naive, filter out 5, 6 for stage 1 nnaive
with(subset(dat_proc, Trialstage==1 & ph1.D43), table(EventIndOmicronD43M6hotdeck1, Country, Bserostatus))
with(subset(dat_proc, Trialstage==1 & ph2.D43.nAb), table(EventIndOmicronD43M6hotdeck1, Country, Bserostatus))

# if we use M12 followup, 5 is not filtered out for stage 1 naive, but we don't study stage 1 naive in correlates
with(subset(dat_proc, Trialstage==1 & ph1.D43), table(EventIndOmicronD43M12hotdeck1, Country, Bserostatus))
with(subset(dat_proc, Trialstage==1 & ph2.D43.nAb), table(EventIndOmicronD43M12hotdeck1, Country, Bserostatus))


# (2) min(np,nv) < N where np is the number of eligible per-protocol placebo non-cases with D01 and D43 nAb 
#     and nv is the number of eligible vaccine recipient non-cases with D01 and D43 nAb. 

# suppose N=1

# for stage 2 naive filter out country 1,2,4,7,10 and keep 6 and 9 only
# for stage 2 nnaive filter out none
with(subset(dat_proc, Trialstage==2 & ph2.D43.nAb & !EventIndOmicronD43M6hotdeck1), table(Trt, Country, Bserostatus))

# for stage 1 naive filter out country 4,6,7 
# for stage 1 nnaive filter out 6,7
with(subset(dat_proc, Trialstage==1 & ph2.D43.nAb & !EventIndOmicronD43M6hotdeck1), table(Trt, Country, Bserostatus))


nrow(subset(dat_proc, Wstratum.original2==60))
nrow(subset(dat_proc, Wstratum.original2==60 & SubcohortInd ))

nrow(subset(dat_proc, Wstratum.original2==51))
nrow(subset(dat_proc, Wstratum.original2==51 & SubcohortInd ))

with(subset(dat_proc, Trialstage==2 & ph1.D43 & SubcohortInd), table(Senior, Country, Trt, Bserostatus))



with(subset(dat_proc, Trialstage==2), table(!is.na(Day43pseudoneutid50), SubcohortInd))

with(subset(dat_proc, Trialstage==2), table(!is.na(Day43bindSpike), SubcohortInd))

with(subset(dat_proc, Trialstage==2), table(!is.na(Day43pseudoneutid50), !is.na(Day43bindSpike), SubcohortInd, Trt, Bserostatus))

dat_proc$ph1=dat_proc$ph1.D43
dat_proc$ph2=dat_proc$ph2.D43

dat_proc$ph2.all.nAb = !is.na(dat_proc$Day43pseudoneutid50)
dat_proc$batch2 = dat_proc$ph2.all.nAb & !dat_proc$ph2

with(subset(dat_proc, Trialstage==2), table(!is.na(Day43pseudoneutid50)))
with(subset(dat_proc, Trialstage==2), table(!is.na(Day43pseudoneutid50), SubcohortInd))

with(subset(dat_proc, Trialstage==2), table(!is.na(Day43bindSpike)))
with(subset(dat_proc, Trialstage==2), table(!is.na(Day43bindSpike), SubcohortInd))

with(subset(dat_proc, Trialstage==2 & ph1.D43 & Trt==1), table(EventIndOmicronD43M12hotdeck1, EventIndFirstInfectionD43))

with(subset(dat_proc, Trialstage==2 & ph1.D43 & Trt==1 & Bserostatus==1 & !is.na(Day43pseudoneutid50)), table(SubcohortInd, EventIndFirstInfectionD43))
with(subset(dat_proc, Trialstage==2 & ph1.D43 & Trt==1 & Bserostatus==0 & !is.na(Day43pseudoneutid50)), table(SubcohortInd, EventIndFirstInfectionD43))
with(subset(dat_proc, Trialstage==2 & ph1.D43 & Trt==0 & Bserostatus==1 & !is.na(Day43pseudoneutid50)), table(SubcohortInd, EventIndFirstInfectionD43))
with(subset(dat_proc, Trialstage==2 & ph1.D43 & Trt==0 & Bserostatus==0 & !is.na(Day43pseudoneutid50)), table(SubcohortInd, EventIndFirstInfectionD43))


################################################################################
# explore data

# SubcohortInd 
# in stage 1 data,
# 1 if IMMSUB1=="Y" or Subjectid in the list of p3005_stage1_step2_sampling_list_28MAR2022 or ADSL.COUNTRY %in% c("JPN","USA"), else 0 
# in stage 2 data, 
# 1 if Subjectid in immunosubset$USUBJID, else 0 



# country 4 (India) and 11 (Ukraine) have no cases
with(subset(dat_proc, Trt==1 & Trialstage==2 & Bserostatus==1 & ph1.D43), table(EventIndOmicronD43M12hotdeck1, Country))
with(subset(dat_proc, Trt==1 & Trialstage==2 & Bserostatus==1 & ph1.D22), table(EventIndOmicronD22M12hotdeck1, Country))

with(subset(dat_proc, Trt==1 & Trialstage==2 & Bserostatus==0 & ph1.D43), table(EventIndOmicronD43M12hotdeck1, Country))
with(subset(dat_proc, Trt==1 & Trialstage==2 & Bserostatus==0 & ph1.D22), table(EventIndOmicronD22M12hotdeck1, Country))

with(subset(dat_proc, Trt==1 & Trialstage==1 & Bserostatus==1 & ph1.D43), table(EventIndOmicronD43M12hotdeck1, Country))
with(subset(dat_proc, Trt==1 & Trialstage==1 & Bserostatus==1 & ph1.D22), table(EventIndOmicronD22M12hotdeck1, Country))




## do Sex improve prediction?
## not in placeo

with(subset(dat_proc, Trt==0 & Bserostatus==1), fastauc(standardized_risk_score, EventIndOmicronD22M12hotdeck1, quiet=FALSE))
glm.fit=glm(EventIndOmicronD22M12hotdeck1~Sex+standardized_risk_score, subset(dat_proc, Trt==0 & Bserostatus==1), family=binomial())
trainauc.glm(glm.fit)

with(subset(dat_proc, Trt==0 & Bserostatus==1 & Trialstage==2), fastauc(standardized_risk_score, EventIndOmicronD22M12hotdeck1, quiet=FALSE))
glm.fit=glm(EventIndOmicronD22M12hotdeck1~Sex+standardized_risk_score, subset(dat_proc, Trt==0 & Bserostatus==1 & Trialstage==2), family=binomial())
trainauc.glm(glm.fit)

with(subset(dat_proc, Trt==1 & Bserostatus==1 & Trialstage==2), fastauc(standardized_risk_score, EventIndOmicronD22M12hotdeck1, quiet=FALSE))
glm.fit=glm(EventIndOmicronD22M12hotdeck1~Sex+standardized_risk_score, subset(dat_proc, Trt==1 & Bserostatus==1 & Trialstage==2), family=binomial())
trainauc.glm(glm.fit)

with(subset(dat_proc, Trt==0 & Bserostatus==1 & Trialstage==2), table(Sex, EventIndOmicronD22M12hotdeck1))
with(subset(dat_proc, Trt==1 & Bserostatus==1 & Trialstage==2), table(Sex, EventIndOmicronD22M12hotdeck1))

with(subset(dat_proc, Bserostatus==1 & Trialstage==2 & ph1.D22), table(Trt, EventIndOmicronD22M12hotdeck1, Sex))

with(subset(dat_proc, Bserostatus==1 & Trialstage==1 & ph1.D22), table(Trt, EventIndOmicronD22M12hotdeck1, Sex))



##

with(subset(dat_proc,ph1),  table(EventIndFirstInfectionD1, batch2, Bserostatus, Trt, Trialstage))
with(subset(dat_proc,!ph1), table(EventIndFirstInfectionD1, batch2, Bserostatus, Trt, Trialstage))

with(subset(dat_proc,EventIndFirstInfectionD1==1), table(ph1.D43, batch2))
with(subset(dat_proc,EventIndFirstInfectionD1==1), table(ph1.D22, batch2))

with(subset(dat_proc,EventIndOmicronD1M12hotdeck1==1), table(ph1.D43, batch2))
with(subset(dat_proc,EventIndOmicronD1M12hotdeck1==1), table(ph1.D22, batch2))


with(subset(dat_proc, Trt==1 & Trialstage==2 & Bserostatus==1 & ph1), table(ph2, ph2.D43.nAb, is.na(Bpseudoneutid50), SubcohortInd))

# subset(dat_proc, Trt==1 & Trialstage==2 & Bserostatus==1 & ph1 & !ph2 & !ph2.D43.nAb & !is.na(Day43pseudoneutid50))[1:2,]

library(kyotil)

par(mfrow=c(1,9), mar=c(4,1,3,1), oma=c(0,0,0,0))

# ph2 
myboxplot(Day43pseudoneutid50~EventIndOmicronD43M12hotdeck1, 
          subset(dat_proc, Trt==1 & Trialstage==2 & Bserostatus==1 & ph2 &
                   Region=="MEX_Stage2"), 
          main='batch1', ylim=c(2,4), test='w', col=4)
myboxplot(Day43pseudoneutid50~EventIndOmicronD43M12hotdeck1, 
          subset(dat_proc, Trt==1 & Trialstage==2 & Bserostatus==1 & ph2.D43.nAb &
                   Region=="MEX_Stage2"), 
          main='+batch1.5', ylim=c(2,4), test='w', col=2)
title(main="MEX", line=2)
myboxplot(Day43pseudoneutid50~EventIndOmicronD43M12hotdeck1, 
          subset(dat_proc, Trt==1 & Trialstage==2 & Bserostatus==1 &
                   Region=="MEX_Stage2"), 
          main='+batch2', ylim=c(2,4), test='w')

myboxplot(Day43pseudoneutid50~EventIndOmicronD43M12hotdeck1, 
          subset(dat_proc, Trt==1 & Trialstage==2 & Bserostatus==1 & ph2 &
                   Region=="NotMEXIND_Stage2"), 
          main='batch1', ylim=c(2,4), test='w', col=4)
myboxplot(Day43pseudoneutid50~EventIndOmicronD43M12hotdeck1, 
          subset(dat_proc, Trt==1 & Trialstage==2 & Bserostatus==1 & ph2.D43.nAb &
                   Region=="NotMEXIND_Stage2"), 
          main='+batch1.5', ylim=c(2,4), test='w', col=2)
title(main="NotMEXIND", line=2)
myboxplot(Day43pseudoneutid50~EventIndOmicronD43M12hotdeck1, 
          subset(dat_proc, Trt==1 & Trialstage==2 & Bserostatus==1 &
                   Region=="NotMEXIND_Stage2"), 
          main='+batch2', ylim=c(2,4), test='w')

myboxplot(Day43pseudoneutid50~EventIndOmicronD43M12hotdeck1, 
          subset(dat_proc, Trt==1 & Trialstage==2 & Bserostatus==1 & ph2 &
                   Region!="IND_Stage2"), 
          main='batch1', ylim=c(2,4), test='w', col=4)
myboxplot(Day43pseudoneutid50~EventIndOmicronD43M12hotdeck1, 
          subset(dat_proc, Trt==1 & Trialstage==2 & Bserostatus==1 & ph2.D43.nAb &
                   Region!="IND_Stage2"), 
          main='+batch1.5', ylim=c(2,4), test='w', col=2)
title(main="Not Ind", line=2)
myboxplot(Day43pseudoneutid50~EventIndOmicronD43M12hotdeck1, 
          subset(dat_proc, Trt==1 & Trialstage==2 & Bserostatus==1 &
                   Region!="IND_Stage2"), 
          main='+batch2', ylim=c(2,4), test='w')




# batch effect is significant in MEX
par(mfrow=c(1,3), mar=c(4,4,3,2))

myboxplot(Day43pseudoneutid50~batch2, 
          subset(dat_proc, Trt==1 & Trialstage==2 & Bserostatus==1 &
                   !EventIndOmicronD43M12hotdeck1 & Country==9), 
          main='MEX, Non-cases', ylim=c(2,4), test='w', names=c("Batch 1", "Batch 2"), ylab='ID50, Day43')

myboxplot(Day43pseudoneutid50~batch2, 
          subset(dat_proc, Trt==1 & Trialstage==2 & Bserostatus==1 &
                   !EventIndOmicronD43M12hotdeck1 & Country==4), 
          main='IND, Non-cases', ylim=c(2,4), test='w', names=c("Batch 1", "Batch 2"), ylab='ID50, Day43')

myboxplot(Day43pseudoneutid50~batch2, 
          subset(dat_proc, Trt==1 & Trialstage==2 & Bserostatus==1 &
                   !EventIndOmicronD43M12hotdeck1 & Region=='NotMEXIND_Stage2'), 
          main='NonINDMEX, Non-cases', ylim=c(2,4), test='w', names=c("Batch 1", "Batch 2"), ylab='ID50, Day43')


myboxplot(FOI~batch2, 
          subset(dat_proc, Trt==1 & Trialstage==2 & Bserostatus==1 &
                   !EventIndOmicronD43M12hotdeck1 & Country==9 & ph2.all.nAb), 
          main='MEX, Non-cases', ylim=NULL, test='w', names=c("Batch 1", "Batch 2"), ylab='FOI')

myboxplot(standardized_risk_score~batch2, 
          subset(dat_proc, Trt==1 & Trialstage==2 & Bserostatus==1 &
                   !EventIndOmicronD43M12hotdeck1 & Country==9 & ph2.all.nAb), 
          main='MEX, Non-cases', ylim=NULL, test='w', names=c("Batch 1", "Batch 2"), ylab='standardized_risk_score')



# bAb
# region confounding may be more serious for bAb than for nAb, which could explain why nAb is a better correlate

par(mfrow=c(1,12), mar=c(4,1,1,1))

# all
myboxplot(Day43bindSpike~EventIndOmicronD43M12hotdeck1, 
          subset(dat_proc, Trt==1 & Trialstage==2 & Bserostatus==1 &
                   Region=="IND_Stage2"), 
          main='Ind', ylim=c(2,4))

myboxplot(Day43bindSpike~EventIndOmicronD43M12hotdeck1, 
          subset(dat_proc, Trt==1 & Trialstage==2 & Bserostatus==1 &
                   Region=="MEX_Stage2"), 
          main='Mex', ylim=c(2,4), test='w')

myboxplot(Day43bindSpike~EventIndOmicronD43M12hotdeck1, 
          subset(dat_proc, Trt==1 & Trialstage==2 & Bserostatus==1 &
                   Region=="NotMEXIND_Stage2"), 
          main='NotMEXIND', ylim=c(2,4), test='w')

myboxplot(Day43bindSpike~EventIndOmicronD43M12hotdeck1, 
          subset(dat_proc, Trt==1 & Trialstage==2 & Bserostatus==1 &
                   Region!="IND_Stage2"), 
          main='Not Ind', ylim=c(2,4), test='w')


# ph2.D43.nAb 
myboxplot(Day43bindSpike~EventIndOmicronD43M12hotdeck1, 
          subset(dat_proc, Trt==1 & Trialstage==2 & Bserostatus==1 & ph2.D43.nAb &
                   Region=="IND_Stage2"), 
          main='Ind', ylim=c(2,4), col=2)

myboxplot(Day43bindSpike~EventIndOmicronD43M12hotdeck1, 
          subset(dat_proc, Trt==1 & Trialstage==2 & Bserostatus==1 & ph2.D43.nAb &
                   Region=="MEX_Stage2"), 
          main='Mex', ylim=c(2,4), test='w', col=2)

myboxplot(Day43bindSpike~EventIndOmicronD43M12hotdeck1, 
          subset(dat_proc, Trt==1 & Trialstage==2 & Bserostatus==1 & ph2.D43.nAb &
                   Region=="NotMEXIND_Stage2"), 
          main='NotMEXIND', ylim=c(2,4), test='w', col=2)

myboxplot(Day43bindSpike~EventIndOmicronD43M12hotdeck1, 
          subset(dat_proc, Trt==1 & Trialstage==2 & Bserostatus==1 & ph2.D43.nAb &
                   Region!="IND_Stage2"), 
          main='Not Ind', ylim=c(2,4), test='w')


# moreph2 
myboxplot(Day43bindSpike~EventIndOmicronD43M12hotdeck1, 
          subset(dat_proc, Trt==1 & Trialstage==2 & Bserostatus==1 & ph2 &
                   Region=="IND_Stage2"), 
          main='Ind', ylim=c(2,4), col=4)

myboxplot(Day43bindSpike~EventIndOmicronD43M12hotdeck1, 
          subset(dat_proc, Trt==1 & Trialstage==2 & Bserostatus==1 & ph2 &
                   Region=="MEX_Stage2"), 
          main='Mex', ylim=c(2,4), test='w', col=4)

myboxplot(Day43bindSpike~EventIndOmicronD43M12hotdeck1, 
          subset(dat_proc, Trt==1 & Trialstage==2 & Bserostatus==1 & ph2 &
                   Region=="NotMEXIND_Stage2"), 
          main='NotMEXIND', ylim=c(2,4), test='w', col=4)

myboxplot(Day43bindSpike~EventIndOmicronD43M12hotdeck1, 
          subset(dat_proc, Trt==1 & Trialstage==2 & Bserostatus==1 & ph2 &
                   Region!="IND_Stage2"), 
          main='Not Ind', ylim=c(2,4), test='w')



# number of countries
with(subset(dat_proc, Trt==1 & Trialstage==2 & Bserostatus==1 & Region!='IND_Stage2' & ph1), table(Country) )



################################################################################
# Cox models

# the original
dat_proc$ph2=dat_proc$ph2.D43.original
dat_proc$Wstratum=dat_proc$Wstratum.original
svycoxph(Surv(EventTimeOmicronD43M12hotdeck1, EventIndOmicronD43M12hotdeck1) ~ Sex + Day43pseudoneutid50 + strata(Country),
         twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, 
                  data=subset(dat_proc, Trt==1 & Trialstage==2 & Bserostatus==1 & ph1.D43)) )

# all nAb data
dat_proc$ph2=dat_proc$ph2.D43.nAb
dat_proc$Wstratum=dat_proc$Wstratum.original
svycoxph(Surv(EventTimeOmicronD43M12hotdeck1, EventIndOmicronD43M12hotdeck1) ~ Sex + Day43pseudoneutid50 + strata(Country),
         twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, 
                  data=subset(dat_proc, Trt==1 & Trialstage==2 & Bserostatus==1 & ph1.D43)) )


dat_proc$ph2=dat_proc$ph2.D43.nAb & dat_proc$ph2.D43.original
dat_proc$Wstratum=dat_proc$Wstratum.original
svycoxph(Surv(EventTimeOmicronD43M12hotdeck1, EventIndOmicronD43M12hotdeck1) ~ Sex + Day43pseudoneutid50 + strata(Country),
         twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, 
                  data=subset(dat_proc, Trt==1 & Trialstage==2 & Bserostatus==1 & ph1.D43)) )


dat_proc$ph2=dat_proc$ph2.D43.nAb & (dat_proc$nAbBatch==1 | dat_proc$EventIndOmicronD43M12hotdeck1)
dat_proc$Wstratum=dat_proc$Wstratum.original
svycoxph(Surv(EventTimeOmicronD43M12hotdeck1, EventIndOmicronD43M12hotdeck1) ~ Sex + Day43pseudoneutid50 + strata(Country),
         twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, 
                  data=subset(dat_proc, Trt==1 & Trialstage==2 & Bserostatus==1 & ph1.D43)) )


dat_proc$ph2=dat_proc$ph2.D43.nAb & (dat_proc$ph2.D43.bAb==1 | dat_proc$EventIndOmicronD43M12hotdeck1)
dat_proc$Wstratum=dat_proc$Wstratum.original
svycoxph(Surv(EventTimeOmicronD43M12hotdeck1, EventIndOmicronD43M12hotdeck1) ~ Sex + Day43pseudoneutid50 + strata(Country),
         twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, 
                  data=subset(dat_proc, Trt==1 & Trialstage==2 & Bserostatus==1 & ph1.D43)) )


dat_proc$ph2=dat_proc$ph2.D43.nAb & (dat_proc$ph2.D43.bAb==1 & dat_proc$nAbBatch==1 | dat_proc$EventIndOmicronD43M12hotdeck1)
dat_proc$Wstratum=dat_proc$Wstratum.original
svycoxph(Surv(EventTimeOmicronD43M12hotdeck1, EventIndOmicronD43M12hotdeck1) ~ Sex + Day43pseudoneutid50 + strata(Country),
         twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, 
                  data=subset(dat_proc, Trt==1 & Trialstage==2 & Bserostatus==1 & ph1.D43)) )


dat_proc$ph2=dat_proc$ph2.D43.nAb & (dat_proc$ph2.D43.bAb==1 | dat_proc$EventIndOmicronD43M12hotdeck1)
dat_proc$Wstratum=dat_proc$Wstratum.nAb
svycoxph(Surv(EventTimeOmicronD43M12hotdeck1, EventIndOmicronD43M12hotdeck1) ~ Sex + Day43pseudoneutid50 + strata(Country),
         twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, 
                  data=subset(dat_proc, Trt==1 & Trialstage==2 & Bserostatus==1 & ph1.D43)) )




with (dat_proc, table(ph2.D43.nAb & (nAbBatch==1 | EventIndOmicronD43M12hotdeck1), ph2.D43.nAb & ph2.D43.original))
with (subset(dat_proc, Trialstage==2 & Bserostatus==1 & Trt==1 & ph1.D43==1), table(ph2.D43.nAb & (nAbBatch==1 | EventIndOmicronD43M12hotdeck1), ph2.D43.nAb & ph2.D43.original))
with (dat_proc, table(ph2.D43.nAb & (nAbBatch==1 | EventIndOmicronD43M12hotdeck1), ph2.D43.nAb & SubcohortInd))

with(subset(dat_proc, Trialstage==2 & Bserostatus==1 & Trt==1 & ph1.D43==1 & SubcohortInd==1),
     table(complete.cases(dat_proc[,c("B", "Day43")%.%'bindSpike']), complete.cases(dat_proc[,c("B", "Day43")%.%'pseudoneutid50'])))


# is it possible to adjust both as.factor(Country) and strata(Country)?
# no
svycoxph(Surv(EventTimeOmicronD43M12hotdeck1, EventIndOmicronD43M12hotdeck1) ~ 
           Sex + FOI + standardized_risk_score + Day43pseudoneutid50 + strata(Country),
         twophase(id=list(~1,~1), strata=list(NULL,~Wstratum.nAb), subset=~ph2.D43.nAb, data=subset(dat_proc, Trt==1 & Trialstage==2 & Bserostatus==1 & Region!='IND_Stage2' & ph1.D43) ))

svycoxph(Surv(EventTimeOmicronD43M12hotdeck1, EventIndOmicronD43M12hotdeck1) ~ 
           Sex + FOI + standardized_risk_score + Day43pseudoneutid50 + as.factor(Country),
         twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=subset(dat_proc, Trt==1 & Trialstage==2 & Bserostatus==1 & Region!='IND_Stage2' & ph1) ))

svycoxph(Surv(EventTimeOmicronD43M12hotdeck1, EventIndOmicronD43M12hotdeck1) ~ 
           Sex + FOI + standardized_risk_score + Day43pseudoneutid50 + strata(Region) + as.factor(Country),
         twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=subset(dat_proc, Trt==1 & Trialstage==2 & Bserostatus==1 & Region!='IND_Stage2' & ph1) ))



# using batch 1 data and pre-specified model
svycoxph(Surv(EventTimeOmicronD43M12hotdeck1, EventIndOmicronD43M12hotdeck1) ~ 
           # Sex + FOI + standardized_risk_score + Day43pseudoneutid50 + as.factor(Country),
           # Sex + FOI + standardized_risk_score + Day43pseudoneutid50 + strata(Country),
         twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=subset(dat_proc, Trt==1 & Trialstage==2 & Bserostatus==1 & Region!='IND_Stage2' & ph1) ))

# sex is a precision variable
svycoxph(Surv(EventTimeOmicronD43M12hotdeck1, EventIndOmicronD43M12hotdeck1) ~ 
           FOI + standardized_risk_score + Day43pseudoneutid50 + strata(Region) + Sex, 
         twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=subset(dat_proc, Trt==1 & Trialstage==2 & Bserostatus==1 & Region!='IND_Stage2' & ph1) ))

# restrict to NotMEXIND, not significant
svycoxph(Surv(EventTimeOmicronD43M12hotdeck1, EventIndOmicronD43M12hotdeck1) ~ 
           FOI + standardized_risk_score + Day43pseudoneutid50 + strata(Region), 
         twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=subset(dat_proc, Trt==1 & Trialstage==2 & Bserostatus==1 & Region=='NotMEXIND_Stage2' & ph1) ))

# ph2.D43.nAb
svycoxph(Surv(EventTimeOmicronD43M12hotdeck1, EventIndOmicronD43M12hotdeck1) ~ 
           FOI + standardized_risk_score + Day43pseudoneutid50 + strata(Region), 
         twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2.D43.nAb, data=subset(dat_proc, Trt==1 & Trialstage==2 & Bserostatus==1 & Region!='IND_Stage2' & ph1) ))

# ph2.D43.nAb, strata(Country)
svycoxph(Surv(EventTimeOmicronD43M12hotdeck1, EventIndOmicronD43M12hotdeck1) ~ 
           FOI + standardized_risk_score + Sex + Day43pseudoneutid50 + strata(Region), 
         twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2.all.nAb.D43, data=subset(dat_proc, Trt==1 & Trialstage==2 & Bserostatus==1 & Region!='IND_Stage2' & ph1) ))



with(subset(dat_proc, Trialstage==2 & Bserostatus==1 & Region=='IND_Stage2' ), table(Trt, EventIndOmicronD43M12hotdeck1))
with(subset(dat_proc, Trialstage==2 & Bserostatus==1 & Region=='IND_Stage2' ), table(Trt, EventIndFirstInfectionD1))



with(subset(dat_proc, Trialstage==2 & Bserostatus==1), table(ph2, ph2.D43.nAb, EventIndFirstInfectionD43, Trt))


