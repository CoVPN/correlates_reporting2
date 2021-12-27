# putting controlled VE curves in one plot
renv::activate(project = here::here(".."))    
    # There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
    if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))
    
library(kyotil)
library(tools) # toTitleCase
library(xtable) # this is a dependency of kyotil
library(survey)
library(Hmisc)
source(here::here("code", "params.R"))
Sys.setenv(VERBOSE=1)
print("meta ...")
ve.az=read.csv("../data_clean/AZChAd26UKphase3FengetalCorrelates.csv")


# TRIALS is a subset of all.trials
# a is an assay
draw.ve.curves=function(a, TRIALS, file.name, include.az=FALSE) {
#a="bindSpike"; TRIALS=c("moderna_real", "janssen_pooled_real"); file.name=1
    myprint(a)
    ylim=c(0, 1)    
    hist.shrink=c(ADCP=3,pseudoneutid50=3,bindSpike=3,bindRBD=3)
    
    all.trials=c("moderna_real", "janssen_pooled_real", "janssen_na_real", "janssen_la_real", "janssen_sa_real", "AZ-COV002")
    studies=c("COVE","ENSEMBLE","ENSEMBLE NA","ENSEMBLE LA","ENSEMBLE SA","AZ-COV002"); names(studies)=all.trials
    cols=  c("blue","green","green","olivedrab3","darkseagreen4","orange"); names(cols)=all.trials
    hist.col.ls=lapply(cols, function(col) {hist.col <- c(col2rgb(col)); rgb(hist.col[1], hist.col[2], hist.col[3], alpha=255*0.3, maxColorValue=255)})
    
    .subset=match(TRIALS, all.trials)
    
    ## get markers data
    ## get xlim by combining trials. Do ENSEMBLE first because we want to source _common.R for moderna last so that we get the proper lloxs
    ## source _commom.R
    xlim.ls=list()
    marker=list()
    for (x in TRIALS) {    
        TRIAL=get.trial(x, a)
        myprint(TRIAL)
        Sys.setenv("TRIAL"=TRIAL)
        COR = ifelse (x=="moderna_real","D57","D29IncludeNotMolecConfirmedstart1")
        # keep to have local = T
        source(here::here("..", "_common.R"), local=T)
        
        # uloq censoring    
        tmp="Day"%.%config.cor$tpeak %.% a
        dat.mock[[tmp]] <- ifelse(dat.mock[[tmp]] > log10(uloqs[a]), log10(uloqs[a]), dat.mock[[tmp]])
        
        dat.vac.seroneg=subset(dat.mock, Trt==1 & ph1)
        xlim=range(dat.vac.seroneg[[tmp]], log10(llods[a]/2), na.rm=T)
        delta=(xlim[2]-xlim[1])/20     
        xlim.ls[[x]]=c(xlim[1]-delta, xlim[2]+delta)
        
        marker[[x]]=dat.vac.seroneg[[tmp]]
    }    
    xlim=c(min(sapply(xlim.ls, function(x) x[1])), max(sapply(xlim.ls, function(x) x[2])))
    myprint(xlim)
    
    mypdf(file=paste0("output/meta_controlled_ve_curves_",file.name,"_",a))
        par(las=1, cex.axis=0.9, cex.lab=1)# axis label orientation
        
        # need several variables from sourcing _common.R: lloxs, labels.assays, draw.x.axis.cor
        overall.ve.ls=list()
        for (x in TRIALS) {    
            TRIAL=get.trial(x, a)
            COR=ifelse (x=="moderna_real","D57","D29IncludeNotMolecConfirmedstart1")            
            load(here::here("output", TRIAL, COR, "marginalized.risk.no.marker.Rdata"))
            load(here::here("output", TRIAL, COR, "marginalized.risk.Rdata"))
            risks=get("risks.all.1")[[a]]        
            
            est = 1 - risks$prob/res.plac.cont["est"]
            boot = 1 - t( t(risks$boot)/res.plac.cont[2:(1+ncol(risks$boot))] )                         
            ci.band=apply(boot, 1, function (x) quantile(x, c(.025,.975)))                
        
            shown=risks$marker>=ifelse(x=="moderna_real",log10(10),quantile(marker[[x]], 2.5/100, na.rm=T)) & risks$marker<=quantile(marker[[x]], 1-2.5/100, na.rm=T)
            mymatplot(risks$marker[shown], t(rbind(est, ci.band))[shown,], type="l", lty=c(1,3,3), lwd=2.5, make.legend=F, col=cols[x], ylab=paste0("Controlled VE"), xlab=labels.assays.short[a]%.%" (=s)", 
                #main=paste0(labels.assays.long["Day"%.%tpeak,a]),
                ylim=ylim, xlim=xlim, yaxt="n", xaxt="n", draw.x.axis=F, add=x!=TRIALS[1])
            draw.x.axis.cor(xlim, NA)
            yat=seq(-1,1,by=.1)
            axis(side=2,at=yat,labels=(yat*100)%.%"%")            
        
            # add histogram
    #        par(new=TRUE) #this changes ylim, so we cannot use it in this loop
            tmp=hist(marker[[x]],breaks=ifelse(x=="moderna_real",25,15),plot=F) # 15 is treated as a suggestion and the actual number of breaks is determined by pretty()
            tmp$density=tmp$density/hist.shrink[a] # so that it will fit vertically
            #tmp=hist(dat.vac.seroneg[["Day"%.%tpeak%.%a]],breaks=seq(min(dat.vac.seroneg[["Day"%.%tpeak%.%a]],na.rm=T), max(dat.vac.seroneg[["Day"%.%tpeak%.%a]],na.rm=T), len = 15),plot=F)
            plot(tmp,col=hist.col.ls[[x]],axes=F,labels=F,main="",xlab="",ylab="",border=0,freq=F,xlim=xlim, ylim=c(0,max(tmp$density*1.25)), add=T) 
            
            overall.ve.ls[[x]]=overall.ve
        }        
    
        # add az curve
        if(include.az) {
            lines(log10(ve.az[[a]]), ve.az$VE/100, col=cols["AZ-COV002"], lwd=2.5)
            lines(log10(ve.az[[a%.%"LL"]]), ve.az$VE/100, col=cols["AZ-COV002"], lwd=2.5, lty=3)
            lines(log10(ve.az[[a%.%"UL"]]), ve.az$VE/100, col=cols["AZ-COV002"], lwd=2.5, lty=3)
        }
    
        # legend
        legend=paste0(studies[TRIALS], ", ",sapply(overall.ve.ls, function(x) formatDouble(x*100,1)%.%"%")[1,])
        if (include.az) legend=c(legend, "AZ-COV002, 66.7%")
        mylegend(x=6, col=cols[c(TRIALS, if(include.az) "AZ-COV002")], legend=legend, lty=1, lwd=2, cex=.7)
    
    dev.off()    
    
}


# COVE + ENSEMBLE + AZ
for (a in c("pseudoneutid50","bindSpike","bindRBD")) {
    draw.ve.curves(a, TRIALS=c("moderna_real", "janssen_pooled_real"), file.name="1", include.az=T)
}


# COVE + ENSEMBLE regions
for (a in c("pseudoneutid50","bindSpike","bindRBD","ADCP")) {
    draw.ve.curves(a, TRIALS=c(if(a!="ADCP") "moderna_real", "janssen_na_real", "janssen_la_real", "janssen_sa_real"), file.name="2", include.az=F)
}
