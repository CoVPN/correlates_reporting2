renv::activate(project = here::here(".."))    
    # There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
    if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))
    
library(kyotil)
library(tools) # toTitleCase
library(xtable) # this is a dependency of kyotil


#eq.geq=3
#for (a in c("bindSpike")) {
a="bindSpike"
 myfigure(file=paste0("output/", a, "_controlled_ve_curves_cove_ensemble"))
    par(las=1, cex.axis=0.9, cex.lab=1)# axis label orientation
    ylim=c(0, 1)
    
    
    ## combine xlim from different trials. do janssen first because we want to source _common.R for moderna last so that we get the proper lloxs
    ## source _commom.R
    xlim.ls=list()
    for (i in 2:1) {        
        Sys.setenv("TRIAL"=ifelse (i==1, "moderna_real", "janssen_pooled_real"))
        COR=ifelse (i==1,"D57","D29IncludeNotMolecConfirmed")
        source(here::here("..", "_common.R"))
        
        # uloq censoring    
        tmp="Day"%.%config.cor$tpeak %.% a
        dat.mock[[tmp]] <- ifelse(dat.mock[[tmp]] > log10(uloqs[a]), log10(uloqs[a]), dat.mock[[tmp]])
        
        dat.vac.seroneg=subset(dat.mock, Trt==1 & ph1)
        xlim=get.range.cor(dat.vac.seroneg, a, config.cor$tpeak)
        xlim.ls[[i]]=xlim
    }
    xlim=c(min(xlim.ls[[1]][1], xlim.ls[[2]][1]), max(xlim.ls[[1]][2], xlim.ls[[2]][2]))
    
    # depends on several variables from sourcing _common.R: lloxs, labels.assays, draw.x.axis.cor
    for (i in 1:2) {
        TRIAL=ifelse (i==1, "moderna_real", "janssen_pooled_real")
        COR=ifelse (i==1,"D57","D29IncludeNotMolecConfirmed")
        study_name=ifelse (i==1,"COVE","ENSEMBLE")
#        config <- config::get(config = Sys.getenv("TRIAL"))
#        config.cor <- config::get(config = COR)
        
        load(here::here("output", TRIAL, COR, "marginalized.risk.no.marker."%.%study_name%.%".Rdata"))
        load(here::here("output", TRIAL, COR, "marginalized.risk."%.%study_name%.%".Rdata"))
        risks=get("risks.all.1")[[a]]        
        
        est = 1 - risks$prob/res.plac.cont["est"]
        boot = 1 - t( t(risks$boot)/res.plac.cont[2:(1+ncol(risks$boot))] )                         
        ci.band=apply(boot, 1, function (x) quantile(x, c(.025,.975)))                
    
        mymatplot(risks$marker, t(rbind(est, ci.band)), type="l", lty=c(1,2,2), lwd=2.5, make.legend=F, col=ifelse(i==1,"blue","green"), ylab=paste0("Controlled VE"), xlab=labels.assays.short[a]%.%" (=s)", 
            #main=paste0(labels.assays.long["Day"%.%tpeak,a]),
            ylim=ylim, xlim=xlim, yaxt="n", xaxt="n", draw.x.axis=F, add=i==2)
        draw.x.axis.cor(xlim, lloxs[a])
        yat=seq(-1,1,by=.1); axis(side=2,at=yat,labels=(yat*100)%.%"%")        
    
    
#        # add histogram
#        par(new=TRUE) 
#        col <- c(col2rgb("olivedrab3")) # orange, darkgoldenrod2
#        col <- rgb(col[1], col[2], col[3], alpha=255*0.4, maxColorValue=255)
#        tmp=hist(dat.vac.seroneg[["Day"%.%tpeak%.%a]],breaks=15,plot=F) # 15 is treated as a suggestion and the actual number of breaks is determined by pretty()
#        #tmp=hist(dat.vac.seroneg[["Day"%.%tpeak%.%a]],breaks=seq(min(dat.vac.seroneg[["Day"%.%tpeak%.%a]],na.rm=T), max(dat.vac.seroneg[["Day"%.%tpeak%.%a]],na.rm=T), len = 15),plot=F)
#        plot(tmp,col=col,axes=F,labels=F,main="",xlab="",ylab="",border=0,freq=F,xlim=xlim, ylim=c(0,max(tmp$density*1.25))) 
    }
    
    # legend
    tmp=formatDouble(overall.ve*100,1)%.%"%"        
    legend.x=9
    mylegend(x=legend.x,legend=c(
            paste0("Overall VE ",tmp[1]," (",tmp[2],", ",tmp[3],")"), 
            "Controlled VE",
            if(eq.geq==1) "Controlled VE Sens. Analysis"), 
        col=c("white", if(eq.geq==3 | eq.geq==2) "black" else "pink", if(eq.geq==1) "red"), 
        lty=1, lwd=2, cex=.8)

  mydev.off()
#} # end for
    
