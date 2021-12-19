renv::activate(project = here::here(".."))    
    # There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
    if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))
    
library(kyotil)
library(tools) # toTitleCase
library(xtable) # this is a dependency of kyotil


# 
ve.az=read.csv("../data_clean/AZChAd26UKphase3FengetalCorrelates.csv")

assays=c("bindSpike","pseudoneutid50","bindRBD")
hist.shrink=c(bindSpike=3,bindRBD=3,pseudoneutid50=3)
studies=c("COVE","ENSEMBLE","AZ-COV002")
cols=c("blue","green","orange")
hist.col.ls=list() # orange, darkgoldenrod2, olivedrab3, lightblue
hist.col <- c(col2rgb("blue")) 
hist.col <- rgb(hist.col[1], hist.col[2], hist.col[3], alpha=255*0.3, maxColorValue=255)
hist.col.ls[[1]]=hist.col
hist.col <- c(col2rgb("green")) 
hist.col <- rgb(hist.col[1], hist.col[2], hist.col[3], alpha=255*0.3, maxColorValue=255)
hist.col.ls[[2]]=hist.col
ylim=c(0, 1)    


for (a in assays) {
#a="bindSpike"
    myprint(a)
    ## combine xlim from COVE and ENSEMBLE. Do ENSEMBLE first because we want to source _common.R for moderna last so that we get the proper lloxs
    ## source _commom.R
    ## save markers
    xlim.ls=list()
    marker=list()
    for (i in 2:1) {        
        TRIAL=ifelse (i==1, "moderna_real", ifelse(startsWith(a,"pseudo"), "janssen_pooled_realPsV", "janssen_pooled_real"))
        Sys.setenv("TRIAL"=TRIAL)
        COR=ifelse (i==1,"D57","D29IncludeNotMolecConfirmedstart1")
        source(here::here("..", "_common.R"))
        
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

    myfigure()
        par(las=1, cex.axis=0.9, cex.lab=1)# axis label orientation
        
        # depends on several variables from sourcing _common.R: lloxs, labels.assays, draw.x.axis.cor
        overall.ve.ls=list()
        # draw ve for i=1 COVE and i=2 ENSEMBLE
        for (i in 1:2) {
            TRIAL=ifelse (i==1, "moderna_real", ifelse(startsWith(a,"pseudo"), "janssen_pooled_realPsV", "janssen_pooled_real"))
            COR=ifelse (i==1,"D57","D29IncludeNotMolecConfirmedstart1")
            study_name=studies[i]
    #        config <- config::get(config = Sys.getenv("TRIAL"))
    #        config.cor <- config::get(config = COR)
            
            load(here::here("output", TRIAL, COR, "marginalized.risk.no.marker."%.%study_name%.%".Rdata"))
            load(here::here("output", TRIAL, COR, "marginalized.risk."%.%study_name%.%".Rdata"))
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
    
    mydev.off(file=paste0("output/meta_controlled_ve_curves_",a))
} # end assays
    
