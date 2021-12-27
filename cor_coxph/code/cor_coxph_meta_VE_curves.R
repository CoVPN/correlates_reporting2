# putting controlled VE curves in one plot
renv::activate(project = here::here(".."))    
    # There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
    if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))
    
library(kyotil)
library(tools) # toTitleCase
library(xtable) # this is a dependency of kyotil
library(survey)
library(Hmisc)
get.trial=function(x, assay) {
    if (startsWith(x,"janssen")) {
        if(startsWith(assay,"pseudo")) x=paste0(x, "PsV")
        if(startsWith(assay,"ADCP")) x=paste0(x, "ADCP")
    }
    x
}
verbose=1
if(verbose) print("meta ...")

ve.az=read.csv("../data_clean/AZChAd26UKphase3FengetalCorrelates.csv")


#### COVE + ENSEMBLE + AZ

assays=c("bindSpike","bindRBD","pseudoneutid50","ADCP")
ylim=c(0, 1)    
hist.shrink=c(bindSpike=3,bindRBD=3,pseudoneutid50=3,ADCP=3)

TRIALS=c("moderna_real", "janssen_pooled_real")
studies=c("COVE","ENSEMBLE","AZ-COV002")
cols=c("blue","green","orange")
hist.col.ls=list() # orange, darkgoldenrod2, olivedrab3, lightblue
for (i in 1:length(cols)) {hist.col <- c(col2rgb(cols[i])); hist.col <- rgb(hist.col[1], hist.col[2], hist.col[3], alpha=255*0.3, maxColorValue=255); hist.col.ls[[i]]=hist.col}



for (a in assays) {
#a="pseudoneutid50"
    myprint(a)
    ## combine xlim from COVE and ENSEMBLE. Do ENSEMBLE first because we want to run _common.R for moderna last so that we get the proper lloxs
    ## source _commom.R
    ## save markers
    xlim.ls=list()
    marker=list()
    for (i in length(TRIALS):1) {    
        TRIAL=get.trial(TRIALS[i], a)
        Sys.setenv("TRIAL"=TRIAL)
        COR=ifelse (i==1,"D57","D29IncludeNotMolecConfirmedstart1")
        if(exists("Args")) rm("Args"); source(here::here("..", "_common.R"))
        
        # uloq censoring    
        tmp="Day"%.%config.cor$tpeak %.% a
        dat.mock[[tmp]] <- ifelse(dat.mock[[tmp]] > log10(uloqs[a]), log10(uloqs[a]), dat.mock[[tmp]])
        
        dat.vac.seroneg=subset(dat.mock, Trt==1 & ph1)
        xlim=range(dat.vac.seroneg[[tmp]], log10(llods[a]/2), na.rm=T)
        delta=(xlim[2]-xlim[1])/20     
        xlim.ls[[i]]=c(xlim[1]-delta, xlim[2]+delta)
        
        marker[[i]]=dat.vac.seroneg[[tmp]]
    }
    xlim=c(min(xlim.ls[[1]][1], xlim.ls[[2]][1]), max(xlim.ls[[1]][2], xlim.ls[[2]][2]))
    myprint(xlim)

    mypdf(file=paste0("output/meta_controlled_ve_curves_",a))
        par(las=1, cex.axis=0.9, cex.lab=1)# axis label orientation
        
        # depends on several variables from sourcing _common.R: lloxs, labels.assays, draw.x.axis.cor
        overall.ve.ls=list()
        # draw ve for i=1 COVE and i=2 ENSEMBLE
        for (i in 1:length(TRIALS)) {    
            TRIAL=get.trial(TRIALS[i], a)
            COR=ifelse (i==1,"D57","D29IncludeNotMolecConfirmedstart1")
            
            load(here::here("output", TRIAL, COR, "marginalized.risk.no.marker.Rdata"))
            load(here::here("output", TRIAL, COR, "marginalized.risk.Rdata"))
            risks=get("risks.all.1")[[a]]        
            overall.ve.ls[[i]]=overall.ve
            
            est = 1 - risks$prob/res.plac.cont["est"]
            boot = 1 - t( t(risks$boot)/res.plac.cont[2:(1+ncol(risks$boot))] )                         
            ci.band=apply(boot, 1, function (x) quantile(x, c(.025,.975)))                
        
            if (i!=1) shown=risks$marker>=quantile(marker[[i]], 2.5/100, na.rm=T) & risks$marker<=quantile(marker[[i]], 1-2.5/100, na.rm=T) else shown = risks$marker>=log10(10) & risks$marker<=quantile(marker[[i]], 1-2.5/100, na.rm=T) 
            mymatplot(risks$marker[shown], t(rbind(est, ci.band))[shown,], type="l", lty=c(1,3,3), lwd=2.5, make.legend=F, col=cols[i], ylab=paste0("Controlled VE"), xlab=labels.assays.short[a]%.%" (=s)", 
                #main=paste0(labels.assays.long["Day"%.%tpeak,a]),
                ylim=ylim, xlim=xlim, yaxt="n", xaxt="n", draw.x.axis=F, add=i==2)
            draw.x.axis.cor(xlim, NA)
            yat=seq(-1,1,by=.1)
            axis(side=2,at=yat,labels=(yat*100)%.%"%")            
        
            # add histogram
    #        par(new=TRUE) #this changes ylim, so we cannot use it in this loop
            tmp=hist(marker[[i]],breaks=ifelse(i==1,25,15),plot=F) # 15 is treated as a suggestion and the actual number of breaks is determined by pretty()
            tmp$density=tmp$density/hist.shrink[a] # so that it will fit vertically
            #tmp=hist(dat.vac.seroneg[["Day"%.%tpeak%.%a]],breaks=seq(min(dat.vac.seroneg[["Day"%.%tpeak%.%a]],na.rm=T), max(dat.vac.seroneg[["Day"%.%tpeak%.%a]],na.rm=T), len = 15),plot=F)
            plot(tmp,col=hist.col.ls[[i]],axes=F,labels=F,main="",xlab="",ylab="",border=0,freq=F,xlim=xlim, ylim=c(0,max(tmp$density*1.25)), add=T) 
        }
        
        # add az curve
        lines(log10(ve.az[[a]]), ve.az$VE/100, col=cols[3], lwd=2.5)
        lines(log10(ve.az[[a%.%"LL"]]), ve.az$VE/100, col=cols[3], lwd=2.5, lty=3)
        lines(log10(ve.az[[a%.%"UL"]]), ve.az$VE/100, col=cols[3], lwd=2.5, lty=3)
    
        # legend
        tmp.1=formatDouble(overall.ve.ls[[1]]*100,1)%.%"%"        
        tmp.2=formatDouble(overall.ve.ls[[2]]*100,1)%.%"%"        
        mylegend(x=6, col=cols, legend=c(
                paste0(studies[1], " (overall ",tmp.1[1], ")"), 
                paste0(studies[2], " (overall ",tmp.2[1], ")"),
                paste0(studies[3], " (overall 66.7%)") # based on Feng et al
            ), lty=1, lwd=2, cex=.7)
    
    dev.off()

} # end assays
    


#### COVE + ENSEMBLE regions

ylim=c(0, 1)    
assays=c("pseudoneutid50","bindSpike","bindRBD","ADCP")
hist.shrink=c(ADCP=3,pseudoneutid50=3,bindSpike=3,bindRBD=3)

TRIALS=c("moderna_real", "janssen_na_real", "janssen_la_real", "janssen_sa_real")
studies=c("COVE","ENSEMBLE NA","ENSEMBLE LA","ENSEMBLE SA")
cols=  c("blue","green","olivedrab3","purple")
hist.col.ls=list() # orange, darkgoldenrod2, olivedrab3, lightblue
for (i in 1:length(cols)) {hist.col <- c(col2rgb(cols[i])); hist.col <- rgb(hist.col[1], hist.col[2], hist.col[3], alpha=255*0.3, maxColorValue=255); hist.col.ls[[i]]=hist.col}

for (a in assays) {
#a="pseudoneutid50"
    myprint(a)
    .subset=ifelse(a=="ADCP",2,1):length(TRIALS)
    
    ## get markers data
    ## get xlim by combining trials. Do ENSEMBLE first because we want to source _common.R for moderna last so that we get the proper lloxs
    ## source _commom.R
    xlim.ls=list()
    marker=list()
    for (i in rev(.subset)) {    
        TRIAL=get.trial(TRIALS[i], a)
        Sys.setenv("TRIAL"=TRIAL)
        COR=ifelse (i==1,"D57","D29IncludeNotMolecConfirmedstart1")
        if(exists("Args")) rm("Args"); source(here::here("..", "_common.R"))
        
        # uloq censoring    
        tmp="Day"%.%config.cor$tpeak %.% a
        dat.mock[[tmp]] <- ifelse(dat.mock[[tmp]] > log10(uloqs[a]), log10(uloqs[a]), dat.mock[[tmp]])
        
        dat.vac.seroneg=subset(dat.mock, Trt==1 & ph1)
        xlim=range(dat.vac.seroneg[[tmp]], log10(llods[a]/2), na.rm=T)
        delta=(xlim[2]-xlim[1])/20     
        xlim.ls[[i]]=c(xlim[1]-delta, xlim[2]+delta)
        
        marker[[i]]=dat.vac.seroneg[[tmp]]
    }
    xlim=c(min(xlim.ls[[1]][1], xlim.ls[[2]][1], xlim.ls[[3]][1]), max(xlim.ls[[1]][2], xlim.ls[[2]][2], xlim.ls[[3]][2]))
    myprint(xlim)

    mypdf(file=paste0("output/meta_controlled_ve_curves2_",a))
        par(las=1, cex.axis=0.9, cex.lab=1)# axis label orientation
        
        # depends on several variables from sourcing _common.R: lloxs, labels.assays, draw.x.axis.cor
        overall.ve.ls=list()
        # draw ve for i=1 COVE and i>1 ENSEMBLE
        for (i in .subset) {    
            TRIAL=get.trial(TRIALS[i], a)
            COR=ifelse (i==1,"D57","D29IncludeNotMolecConfirmedstart1")
            
            load(here::here("output", TRIAL, COR, "marginalized.risk.no.marker.Rdata"))
            load(here::here("output", TRIAL, COR, "marginalized.risk.Rdata"))
            risks=get("risks.all.1")[[a]]        
            overall.ve.ls[[i]]=overall.ve
            
            est = 1 - risks$prob/res.plac.cont["est"]
            boot = 1 - t( t(risks$boot)/res.plac.cont[2:(1+ncol(risks$boot))] )                         
            ci.band=apply(boot, 1, function (x) quantile(x, c(.025,.975)))                
        
            if (i!=1) shown=risks$marker>=quantile(marker[[i]], 2.5/100, na.rm=T) & risks$marker<=quantile(marker[[i]], 1-2.5/100, na.rm=T) else shown = risks$marker>=log10(10) & risks$marker<=quantile(marker[[i]], 1-2.5/100, na.rm=T) 
            mymatplot(risks$marker[shown], t(rbind(est, ci.band))[shown,], type="l", lty=c(1,3,3), lwd=2.5, make.legend=F, col=cols[i], ylab=paste0("Controlled VE"), xlab=labels.assays.short[a]%.%" (=s)", 
                #main=paste0(labels.assays.long["Day"%.%tpeak,a]),
                ylim=ylim, xlim=xlim, yaxt="n", xaxt="n", draw.x.axis=F, add=i>min(.subset))
            draw.x.axis.cor(xlim, NA)
            yat=seq(-1,1,by=.1)
            axis(side=2,at=yat,labels=(yat*100)%.%"%")            
        
            # add histogram
    #        par(new=TRUE) #this changes ylim, so we cannot use it in this loop
            tmp=hist(marker[[i]],breaks=ifelse(i==1,25,15),plot=F) # 15 is treated as a suggestion and the actual number of breaks is determined by pretty()
            tmp$density=tmp$density/hist.shrink[a] # so that it will fit vertically
            #tmp=hist(dat.vac.seroneg[["Day"%.%tpeak%.%a]],breaks=seq(min(dat.vac.seroneg[["Day"%.%tpeak%.%a]],na.rm=T), max(dat.vac.seroneg[["Day"%.%tpeak%.%a]],na.rm=T), len = 15),plot=F)
            plot(tmp,col=hist.col.ls[[i]],axes=F,labels=F,main="",xlab="",ylab="",border=0,freq=F,xlim=xlim, ylim=c(0,max(tmp$density*1.25)), add=T) 
        }
        
        mylegend(x=6, col=cols[.subset], legend=paste0(studies[.subset], ", ",sapply(overall.ve.ls[.subset], function(x) formatDouble(x*100,1)%.%"%")[1,]), lty=1, lwd=2, cex=.7)
    
    dev.off()

} # end assays
    




####################################################################
# compare controlled VE for ID50 <30, 100> between cove and ensemble

# pick a day for cove such that the overall risk in placebo match that of ensemble: 89
# placebo risk 0.053
get.marginalized.risk.no.marker(Surv(EventTimePrimary, EventIndPrimary) ~ risk_score + MinorityInd + HighRiskInd, dat.cove.0, day=89)
# 0.0530523
get.marginalized.risk.no.marker(Surv(EventTimePrimary, EventIndPrimary) ~ risk_score + MinorityInd + HighRiskInd, dat.cove.1, day=89)
# 0.00328746

get.marginalized.risk.no.marker(Surv(EventTimePrimary, EventIndPrimary) ~ risk_score + as.factor(Region), dat.ense.0, day=66)
# 0.0525176
get.marginalized.risk.no.marker(Surv(EventTimePrimary, EventIndPrimary) ~ risk_score + as.factor(Region), dat.ense.1, day=66)
# 0.0196666


# ensemble
dat.vac.seroneg.id50$Day29pseudoneutid50cat2=factor(cut(dat.vac.seroneg.id50$Day29pseudoneutid50, breaks=c(-Inf,log10(c(30,100)),Inf) ))
f1=Surv(EventTimePrimary, EventIndPrimary) ~ risk_score + as.factor(Region) + Day29pseudoneutid50cat2
tmp.design=twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dat.vac.seroneg.id50)
fit.risk=try(svycoxph(f1, design=tmp.design)) # since we don't need se, we could use coxph, but the weights computed by svycoxph are a little different from the coxph due to fpc
prob=marginalized.risk(fit.risk, "Day29pseudoneutid50cat2", data=subset(dat.vac.seroneg.id50, ph2=1), ss=NULL, weights=subset(dat.vac.seroneg.id50, ph2=1)$wt, t=66, categorical.s=T, verbose=F)        
prob
#(-Inf,1.48]    (1.48,2]    (2, Inf] 
# 0.02776172  0.01533599  0.00675338 

# moderna
dat.cove.1$Day57pseudoneutid50cat2=factor(cut(dat.cove.1$Day57pseudoneutid50, breaks=c(-Inf,log10(c(30,100)),Inf) ))
f1=Surv(EventTimePrimary, EventIndPrimary) ~ risk_score + MinorityInd + HighRiskInd + Day57pseudoneutid50cat2
tmp.design=twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dat.cove.1)
fit.risk=try(svycoxph(f1, design=tmp.design)) # since we don't need se, we could use coxph, but the weights computed by svycoxph are a little different from the coxph due to fpc
prob=marginalized.risk(fit.risk, "Day57pseudoneutid50cat2", data=subset(dat.cove.1, ph2=1), ss=NULL, weights=subset(dat.cove.1, ph2=1)$wt, t=89, categorical.s=T, verbose=F)        
prob
#(-Inf,1.48]    (1.48,2]    (2, Inf] 
# 0.00832744  0.00780136  0.00234533 

# Cannot do this analysis for north america/ensemble due to sparsit of data
#> table(dat.ense.1.na$Day29pseudoneutid50cat2, dat.ense.1.na$EventIndPrimary)
#                0   1
#  (-Inf,1.48] 299  22
#  (1.48,2]     62   2
#  (2, Inf]     22   0
