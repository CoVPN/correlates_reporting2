library(survey)

dat.mock = read.csv('/trials/covpn/p3005/analysis/correlates/Part_A_Blinded_Phase_Data/adata/vat08_combined_data_processed_20231118.csv')

dat.mock$ph1=dat.mock$ph1.D43
dat.mock$ph2=dat.mock$ph2.D43

dat.mock$ph2.all.nAb = !is.na(dat.mock$Day43pseudoneutid50)
dat.mock$batch2 = dat.mock$ph2.all.nAb & !dat.mock$ph2


with(subset(dat.mock,ph1),  table(EventIndFirstInfectionD1, batch2, Bserostatus, Trt, Trialstage))
with(subset(dat.mock,!ph1), table(EventIndFirstInfectionD1, batch2, Bserostatus, Trt, Trialstage))

with(subset(dat.mock,EventIndFirstInfectionD1==1), table(ph1.D43, batch2))
with(subset(dat.mock,EventIndFirstInfectionD1==1), table(ph1.D22, batch2))


with(subset(dat.mock, Trt==1 & Trialstage==2 & Bserostatus==1 & ph1), table(ph2, ph2.D43.nAb, is.na(Bpseudoneutid50), SubcohortInd))

# subset(dat.mock, Trt==1 & Trialstage==2 & Bserostatus==1 & ph1 & !ph2 & !ph2.D43.nAb & !is.na(Day43pseudoneutid50))[1:2,]

library(kyotil)

par(mfrow=c(1,9), mar=c(4,1,3,1), oma=c(0,0,0,0))

# ph2 
myboxplot(Day43pseudoneutid50~EventIndOmicronD43M12hotdeck1, 
          subset(dat.mock, Trt==1 & Trialstage==2 & Bserostatus==1 & ph2 &
                   Region=="MEX_Stage2"), 
          main='batch1', ylim=c(2,4), test='w', col=4)
myboxplot(Day43pseudoneutid50~EventIndOmicronD43M12hotdeck1, 
          subset(dat.mock, Trt==1 & Trialstage==2 & Bserostatus==1 & ph2.D43.nAb &
                   Region=="MEX_Stage2"), 
          main='+batch1.5', ylim=c(2,4), test='w', col=2)
title(main="MEX", line=2)
myboxplot(Day43pseudoneutid50~EventIndOmicronD43M12hotdeck1, 
          subset(dat.mock, Trt==1 & Trialstage==2 & Bserostatus==1 &
                   Region=="MEX_Stage2"), 
          main='+batch2', ylim=c(2,4), test='w')

myboxplot(Day43pseudoneutid50~EventIndOmicronD43M12hotdeck1, 
          subset(dat.mock, Trt==1 & Trialstage==2 & Bserostatus==1 & ph2 &
                   Region=="NotMEXIND_Stage2"), 
          main='batch1', ylim=c(2,4), test='w', col=4)
myboxplot(Day43pseudoneutid50~EventIndOmicronD43M12hotdeck1, 
          subset(dat.mock, Trt==1 & Trialstage==2 & Bserostatus==1 & ph2.D43.nAb &
                   Region=="NotMEXIND_Stage2"), 
          main='+batch1.5', ylim=c(2,4), test='w', col=2)
title(main="NotMEXIND", line=2)
myboxplot(Day43pseudoneutid50~EventIndOmicronD43M12hotdeck1, 
          subset(dat.mock, Trt==1 & Trialstage==2 & Bserostatus==1 &
                   Region=="NotMEXIND_Stage2"), 
          main='+batch2', ylim=c(2,4), test='w')

myboxplot(Day43pseudoneutid50~EventIndOmicronD43M12hotdeck1, 
          subset(dat.mock, Trt==1 & Trialstage==2 & Bserostatus==1 & ph2 &
                   Region!="IND_Stage2"), 
          main='batch1', ylim=c(2,4), test='w', col=4)
myboxplot(Day43pseudoneutid50~EventIndOmicronD43M12hotdeck1, 
          subset(dat.mock, Trt==1 & Trialstage==2 & Bserostatus==1 & ph2.D43.nAb &
                   Region!="IND_Stage2"), 
          main='+batch1.5', ylim=c(2,4), test='w', col=2)
title(main="Not Ind", line=2)
myboxplot(Day43pseudoneutid50~EventIndOmicronD43M12hotdeck1, 
          subset(dat.mock, Trt==1 & Trialstage==2 & Bserostatus==1 &
                   Region!="IND_Stage2"), 
          main='+batch2', ylim=c(2,4), test='w')




# batch effect is significant in MEX
par(mfrow=c(1,3), mar=c(4,4,3,2))

myboxplot(Day43pseudoneutid50~batch2, 
          subset(dat.mock, Trt==1 & Trialstage==2 & Bserostatus==1 &
                   !EventIndOmicronD43M12hotdeck1 & Country==9), 
          main='MEX, Non-cases', ylim=c(2,4), test='w', names=c("Batch 1", "Batch 2"), ylab='ID50, Day43')

myboxplot(Day43pseudoneutid50~batch2, 
          subset(dat.mock, Trt==1 & Trialstage==2 & Bserostatus==1 &
                   !EventIndOmicronD43M12hotdeck1 & Country==4), 
          main='IND, Non-cases', ylim=c(2,4), test='w', names=c("Batch 1", "Batch 2"), ylab='ID50, Day43')

myboxplot(Day43pseudoneutid50~batch2, 
          subset(dat.mock, Trt==1 & Trialstage==2 & Bserostatus==1 &
                   !EventIndOmicronD43M12hotdeck1 & Region=='NotMEXIND_Stage2'), 
          main='NonINDMEX, Non-cases', ylim=c(2,4), test='w', names=c("Batch 1", "Batch 2"), ylab='ID50, Day43')


myboxplot(FOI~batch2, 
          subset(dat.mock, Trt==1 & Trialstage==2 & Bserostatus==1 &
                   !EventIndOmicronD43M12hotdeck1 & Country==9 & ph2.all.nAb), 
          main='MEX, Non-cases', ylim=NULL, test='w', names=c("Batch 1", "Batch 2"), ylab='FOI')

myboxplot(risk_score~batch2, 
          subset(dat.mock, Trt==1 & Trialstage==2 & Bserostatus==1 &
                   !EventIndOmicronD43M12hotdeck1 & Country==9 & ph2.all.nAb), 
          main='MEX, Non-cases', ylim=NULL, test='w', names=c("Batch 1", "Batch 2"), ylab='risk_score')



# bAb
# region confounding may be more serious for bAb than for nAb, which could explain why nAb is a better correlate

par(mfrow=c(1,12), mar=c(4,1,1,1))

# all
myboxplot(Day43bindSpike~EventIndOmicronD43M12hotdeck1, 
          subset(dat.mock, Trt==1 & Trialstage==2 & Bserostatus==1 &
                   Region=="IND_Stage2"), 
          main='Ind', ylim=c(2,4))

myboxplot(Day43bindSpike~EventIndOmicronD43M12hotdeck1, 
          subset(dat.mock, Trt==1 & Trialstage==2 & Bserostatus==1 &
                   Region=="MEX_Stage2"), 
          main='Mex', ylim=c(2,4), test='w')

myboxplot(Day43bindSpike~EventIndOmicronD43M12hotdeck1, 
          subset(dat.mock, Trt==1 & Trialstage==2 & Bserostatus==1 &
                   Region=="NotMEXIND_Stage2"), 
          main='NotMEXIND', ylim=c(2,4), test='w')

myboxplot(Day43bindSpike~EventIndOmicronD43M12hotdeck1, 
          subset(dat.mock, Trt==1 & Trialstage==2 & Bserostatus==1 &
                   Region!="IND_Stage2"), 
          main='Not Ind', ylim=c(2,4), test='w')


# ph2.D43.nAb 
myboxplot(Day43bindSpike~EventIndOmicronD43M12hotdeck1, 
          subset(dat.mock, Trt==1 & Trialstage==2 & Bserostatus==1 & ph2.D43.nAb &
                   Region=="IND_Stage2"), 
          main='Ind', ylim=c(2,4), col=2)

myboxplot(Day43bindSpike~EventIndOmicronD43M12hotdeck1, 
          subset(dat.mock, Trt==1 & Trialstage==2 & Bserostatus==1 & ph2.D43.nAb &
                   Region=="MEX_Stage2"), 
          main='Mex', ylim=c(2,4), test='w', col=2)

myboxplot(Day43bindSpike~EventIndOmicronD43M12hotdeck1, 
          subset(dat.mock, Trt==1 & Trialstage==2 & Bserostatus==1 & ph2.D43.nAb &
                   Region=="NotMEXIND_Stage2"), 
          main='NotMEXIND', ylim=c(2,4), test='w', col=2)

myboxplot(Day43bindSpike~EventIndOmicronD43M12hotdeck1, 
          subset(dat.mock, Trt==1 & Trialstage==2 & Bserostatus==1 & ph2.D43.nAb &
                   Region!="IND_Stage2"), 
          main='Not Ind', ylim=c(2,4), test='w')


# moreph2 
myboxplot(Day43bindSpike~EventIndOmicronD43M12hotdeck1, 
          subset(dat.mock, Trt==1 & Trialstage==2 & Bserostatus==1 & ph2 &
                   Region=="IND_Stage2"), 
          main='Ind', ylim=c(2,4), col=4)

myboxplot(Day43bindSpike~EventIndOmicronD43M12hotdeck1, 
          subset(dat.mock, Trt==1 & Trialstage==2 & Bserostatus==1 & ph2 &
                   Region=="MEX_Stage2"), 
          main='Mex', ylim=c(2,4), test='w', col=4)

myboxplot(Day43bindSpike~EventIndOmicronD43M12hotdeck1, 
          subset(dat.mock, Trt==1 & Trialstage==2 & Bserostatus==1 & ph2 &
                   Region=="NotMEXIND_Stage2"), 
          main='NotMEXIND', ylim=c(2,4), test='w', col=4)

myboxplot(Day43bindSpike~EventIndOmicronD43M12hotdeck1, 
          subset(dat.mock, Trt==1 & Trialstage==2 & Bserostatus==1 & ph2 &
                   Region!="IND_Stage2"), 
          main='Not Ind', ylim=c(2,4), test='w')





# Cox models

# using batch 1 data and pre-specified model
svycoxph(Surv(EventTimeOmicronD43M12hotdeck1, EventIndOmicronD43M12hotdeck1) ~ 
           FOI + risk_score + Day43pseudoneutid50 + strata(Region), 
         twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=subset(dat.mock, Trt==1 & Trialstage==2 & Bserostatus==1 & Region!='IND_Stage2' & ph1) ))

# sex is a precision variable
svycoxph(Surv(EventTimeOmicronD43M12hotdeck1, EventIndOmicronD43M12hotdeck1) ~ 
           FOI + risk_score + Day43pseudoneutid50 + strata(Region) + Sex, 
         twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=subset(dat.mock, Trt==1 & Trialstage==2 & Bserostatus==1 & Region!='IND_Stage2' & ph1) ))

# restrict to NotMEXIND, not significant
svycoxph(Surv(EventTimeOmicronD43M12hotdeck1, EventIndOmicronD43M12hotdeck1) ~ 
           FOI + risk_score + Day43pseudoneutid50 + strata(Region), 
         twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=subset(dat.mock, Trt==1 & Trialstage==2 & Bserostatus==1 & Region=='NotMEXIND_Stage2' & ph1) ))

# ph2.D43.nAb
svycoxph(Surv(EventTimeOmicronD43M12hotdeck1, EventIndOmicronD43M12hotdeck1) ~ 
           FOI + risk_score + Day43pseudoneutid50 + strata(Region), 
         twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2.D43.nAb, data=subset(dat.mock, Trt==1 & Trialstage==2 & Bserostatus==1 & Region!='IND_Stage2' & ph1) ))

# ph2.D43.nAb, strata(Country)
svycoxph(Surv(EventTimeOmicronD43M12hotdeck1, EventIndOmicronD43M12hotdeck1) ~ 
           FOI + risk_score + Sex + Day43pseudoneutid50 + strata(Region), 
         twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2.all.nAb.D43, data=subset(dat.mock, Trt==1 & Trialstage==2 & Bserostatus==1 & Region!='IND_Stage2' & ph1) ))



with(subset(dat.mock, Trialstage==2 & Bserostatus==1 & Region=='IND_Stage2' ), table(Trt, EventIndOmicronD43M12hotdeck1))
with(subset(dat.mock, Trialstage==2 & Bserostatus==1 & Region=='IND_Stage2' ), table(Trt, EventIndFirstInfectionD1))



with(subset(dat.mock, Trialstage==2 & Bserostatus==1), table(ph2, ph2.D43.nAb, EventIndFirstInfectionD43, Trt))
