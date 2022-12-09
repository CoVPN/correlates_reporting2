Sys.setenv(TRIAL = "janssen_pooled_EUA"); COR="D29IncludeNotMolecConfirmedstart1"; Sys.setenv(VERBOSE = 1) 
#Sys.setenv(TRIAL = "prevent19"); COR="D35"; Sys.setenv(VERBOSE = 1)
renv::activate(project = here::here(".."))     
    # There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
    if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))    
source(here::here("..", "_common.R"))
library(kyotil) # p.adj.perm, getFormattedSummary
library(marginalizedRisk)
library(survey)
source(here::here("code", "params.R"))
# uloq censoring, done here b/c should not be done for immunogenicity reports
# note that if delta are used, delta needs to be recomputed
for (a in assays) {
  for (t in "Day"%.%tpeak ) {
    dat.mock[[t %.% a]] <- ifelse(dat.mock[[t %.% a]] > log10(uloqs[a]), log10(uloqs[a]), dat.mock[[t %.% a]])
  }
}    
dat.vac.seroneg=subset(dat.mock, Trt==1 & ph1)
dat.pla.seroneg=subset(dat.mock, Trt==0 & ph1)



# dat.vac.seroneg.id50.na etc are defined in cor_coxph_collate_meta_VE.R

#ID50.old=read.csv("/trials/covpn/p3003/analysis/mapping_immune_correlates/adata/COVID_ENSEMBLE_realdata_20211220.csv")
#ID50.new=read.csv("/trials/covpn/p3003/analysis/mapping_immune_correlates/adata/COVID_ENSEMBLE_realdata_20220125.csv")
#
#rbind(subset(ID50.old, Subjectid=="VAC31518COV3001-1155419"), subset(ID50.new, Subjectid=="VAC31518COV3001-1155419"))




fit.risk.1=survival::coxph(Surv(EventTimePrimary, EventIndPrimary) ~ as.factor(Region), subset(dat.vac.seroneg, ph2==1), weights=wt, model=TRUE)         
tmp=subset(dat.vac.seroneg, ph2==1) 
tmp[[time.var]]=66
pred=predict(fit.risk.1, newdata=tmp, type="expected")
weighted.mean(1 - exp(-pred), tmp$wt)
# 0.0293782

fit.risk.1=survival::coxph(Surv(EventTimePrimary, EventIndPrimary) ~ as.factor(Region), subset(dat.vac.seroneg, ph2==1), weights=wt, model=TRUE)         
tmp=subset(dat.vac.seroneg, ph2==1) 
tmp[[time.var]]=55
pred=predict(fit.risk.1, newdata=tmp, type="expected")
weighted.mean(1 - exp(-pred), tmp$wt)
# 0.0142193



#fit.risk.1=survival::coxph(Surv(EventTimePrimary, EventIndPrimary) ~ as.factor(Region), subset(dat.b, ph2==1), weights=wt, model=TRUE)         
#tmp=subset(dat.b, ph2==1) 
#tmp[[time.var]]=66
#pred=predict(fit.risk.1, newdata=tmp, type="expected")
#weighted.mean(1 - exp(-pred), tmp$wt)
## 0.0845114
#
#fit.risk.1=survival::coxph(Surv(EventTimePrimary, EventIndPrimary) ~ as.factor(Region), subset(dat.b, ph2==1), weights=wt, model=TRUE)         
#tmp=subset(dat.b, ph2==1) 
#tmp[[time.var]]=55
#pred=predict(fit.risk.1, newdata=tmp, type="expected")
#weighted.mean(1 - exp(-pred), tmp$wt)
## 0.0131479



fit.risk.1=survival::coxph(Surv(EventTimePrimary, EventIndPrimary) ~ risk_score, subset(dat.vac.seroneg, ph2==1), weights=wt, model=TRUE)         
tmp=subset(dat.vac.seroneg, ph2==1) 
tmp[[time.var]]=66
pred=predict(fit.risk.1, newdata=tmp, type="expected")
weighted.mean(1 - exp(-pred), tmp$wt)
# 0.0188766

fit.risk.1=survival::coxph(Surv(EventTimePrimary, EventIndPrimary) ~ risk_score, subset(dat.vac.seroneg, ph2==1), weights=wt, model=TRUE)         
tmp=subset(dat.vac.seroneg, ph2==1) 
tmp[[time.var]]=55
pred=predict(fit.risk.1, newdata=tmp, type="expected")
weighted.mean(1 - exp(-pred), tmp$wt)
# 0.012488


#fit.risk.1=survival::coxph(Surv(EventTimePrimary, EventIndPrimary) ~ risk_score, subset(dat.b, ph2==1), weights=wt, model=TRUE)         
#tmp=subset(dat.b, ph2==1) 
#tmp[[time.var]]=66
#pred=predict(fit.risk.1, newdata=tmp, type="expected")
#weighted.mean(1 - exp(-pred), tmp$wt)
## 0.0456686
#
#fit.risk.1=survival::coxph(Surv(EventTimePrimary, EventIndPrimary) ~ risk_score, subset(dat.b, ph2==1), weights=wt, model=TRUE)         
#tmp=subset(dat.b, ph2==1) 
#tmp[[time.var]]=55
#pred=predict(fit.risk.1, newdata=tmp, type="expected")
#weighted.mean(1 - exp(-pred), tmp$wt)
## 0.0126222


#
#
#tab=cbind(with(dat.b, table(tps.stratum, EventIndPrimary)), with(dat.b.ph2, table(tps.stratum, EventIndPrimary)), with(dat.b.ph2, table(tps.stratum, EventIndPrimary))[,2]*2.05263157894737)
#
#colSums(tab[1:8,])
#colSums(tab[9:12,])
#colSums(tab[13:16,])
#
#km <- survfit(Surv(EventTimePrimary, EventIndPrimary) ~ as.factor(Region), dat.vac.seroneg)
#km.2 <- survfit(Surv(EventTimePrimary, EventIndPrimary) ~ as.factor(Region), subset(dat.vac.seroneg, ph2==1), weights=wt)
#km.3 <- survfit(Surv(EventTimePrimary, EventIndPrimary) ~ risk_score, subset(dat.vac.seroneg, ph2==1), weights=wt)
#
#km.1 <- survfit(Surv(EventTimePrimary, EventIndPrimary) ~ as.factor(Region), dat.b)
#km.1.2 <- survfit(Surv(EventTimePrimary, EventIndPrimary) ~ as.factor(Region), dat.b.ph2, weights=wt)
#km.1.3 <- survfit(Surv(EventTimePrimary, EventIndPrimary) ~ 1, dat.b.ph2, weights=wt)
#
#par(mfrow=c(2,3))
#plot(km,   col=1:3, main="ph1 real", ylim=c(.9,1), xlim=c(0,66))
#plot(km.2, col=1:3, main="ph2 real", ylim=c(.9,1), xlim=c(0,66))
#plot(km.3, main="ph2 real", ylim=c(.9,1), xlim=c(0,66))
#plot(km.1,   col=1:3, main="ph1 seed=3", ylim=c(.9,1), xlim=c(0,66))
#plot(km.1.2, col=1:3, main="ph2 seed=3", ylim=c(.9,1), xlim=c(0,66))
#plot(km.1.3, main="ph2 seed=3", ylim=c(.9,1), xlim=c(0,66))


#data=dat.b #dat.vac.seroneg
#days=c(37:66); names(days)=days
#out=sapply (days, function(day) {
#    tmp=data
#    fit.risk = coxph(Surv(EventTimePrimary, EventIndPrimary) ~ risk_score + as.factor(Region), tmp)
#    tmp$EventTimePrimary=day
#    res=c(ph1=mean(1 - exp(-predict(fit.risk, newdata=tmp, type="expected"))))
#    
#    tmp=subset(data, ph2==1)     
#    fit.risk.1=coxph(Surv(EventTimePrimary, EventIndPrimary) ~ risk_score + as.factor(Region), tmp, weights=tmp$wt)         
#    tmp$EventTimePrimary=day
#    res=c(res, ph2=weighted.mean(1 - exp(-predict(fit.risk.1, newdata=tmp, type="expected")), tmp$wt))
#    
#    tmp=data
#    tmp.design=twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=tmp)
#    fit.risk=svycoxph(Surv(EventTimePrimary, EventIndPrimary) ~ risk_score + as.factor(Region), design=tmp.design)
#    tmp$EventTimePrimary=day
#    res=c(res, svy=mean(1 - exp(-predict(fit.risk, newdata=tmp, type="expected"))))    
#})
#mymatplot(days, t(out), xlab="Days", col=1:3)




10**with(dat.vac.seroneg.id50.na, wtd.quantile(Day29pseudoneutid50, wt, c(0.025, 0.975)))
10**with(dat.cove.1, wtd.quantile(Day57pseudoneutid50, wt, c(0.025, 0.975)))


with(dat.vac.seroneg.id50.la, summary(EventTimePrimary))
with(subset(dat.vac.seroneg.id50.la, ph2==1), summary(EventTimePrimary))
with(subset(dat.vac.seroneg.id50.la, ph2==1 & EventIndPrimary==0), summary(EventTimePrimary))

tmp=subset(dat.vac.seroneg.id50, ph2==1, select=c(EventIndPrimary, EventTimePrimary))
tail(tmp[order(tmp[,2]),], 50)




tmp=subset(dat.vac.seroneg.id50.na, ph2==1, select=c(EventIndPrimary, EventTimePrimary))
tmp[order(tmp[,2]),]
# case 66 42

tmp=subset(dat.vac.seroneg.id50.la, ph2==1, select=c(EventIndPrimary, EventTimePrimary))
tail(tmp[order(tmp[,2]),], 50)
# case 58 48

tmp=subset(dat.vac.seroneg.id50.sa, ph2==1, select=c(EventIndPrimary, EventTimePrimary))
tmp[order(tmp[,2]),]
# case 26

#tmp=subset(dat.b.ph2, Region==1, select=c(EventIndPrimary, EventTimePrimary))
#tmp[order(tmp[,2]),]
#
#
#subset(dat.b.ph2, Region==1 & EventTimePrimary==58)


with(dat.vac.seroneg, table(EventIndPrimary, !is.na(Day29pseudoneutid50) & !is.na(Bpseudoneutid50) ))
with(dat.vac.seroneg, table(SevereEventIndPrimaryIncludeNotMolecConfirmedD29, !is.na(Day29pseudoneutid50) & !is.na(Bpseudoneutid50) ))


with(dat.vac.seroneg, table(EventIndPrimary, !is.na(Day29pseudoneutid50) & !is.na(Bpseudoneutid50) ))
with(dat.vac.seroneg, table(SevereEventIndPrimaryIncludeNotMolecConfirmedD29, !is.na(Day29pseudoneutid50) & !is.na(Bpseudoneutid50) ))


with(dat.vac.seroneg, summary(Day35bindSpike))



###################################################################################################
