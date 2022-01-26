# putting controlled VE curves in one plot
renv::activate(project = here::here(".."))    
    # There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
    if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))
    
library(kyotil)
library(tools) # toTitleCase
library(xtable) # this is a dependency of kyotil
library(survey)
library(Hmisc)
library(plotrix) # wtd.hist
source(here::here("code", "params.R"))
Sys.setenv(VERBOSE=1)
print("meta ...")
ve.az=read.csv("../data_clean/AZChAd26UKphase3FengetalCorrelates.csv")

# path for figures and tables etc
save.results.to = here::here("output")
if (!dir.exists(save.results.to))  dir.create(save.results.to)
save.results.to = paste0(here::here("output"), "/meta");
if (!dir.exists(save.results.to))  dir.create(save.results.to)



# TRIALS is a subset of all.trials
# a is an assay
draw.ve.curves=function(a, TRIALS, file.name, include.az=FALSE, log="") {
#a="pseudoneutid50"; TRIALS=c("moderna_real", "janssen_pooled_real"); file.name=1; include.az=T; log="y"
    myprint(a)
    if(log=="") ylim=c(0, 1) else ylim=-log(1-c(0,.98))
    hist.shrink=1/c(ADCP=2,pseudoneutid50=1.2,bindSpike=1.3,bindRBD=1.3)
    
    all.trials=c("moderna_real", "janssen_pooled_real", "janssen_na_real", "janssen_la_real", "janssen_sa_real", "AZ-COV002")
    studies=c("COVE","ENSEMBLE","ENSEMBLE US","ENSEMBLE LA","ENSEMBLE SA","AZ-COV002"); names(studies)=all.trials
    cols=  c("blue","green","green","olivedrab3","darkseagreen4","orange"); names(cols)=all.trials
    hist.col.ls=lapply(cols, function(col) {hist.col <- c(col2rgb(col)); rgb(hist.col[1], hist.col[2], hist.col[3], alpha=255*0.3, maxColorValue=255)})
    
    .subset=match(TRIALS, all.trials)
    
    ## get markers data
    ## get xlim by combining trials. Do ENSEMBLE first because we want to source _common.R for moderna last so that we get the proper lloxs
    ## source _commom.R
    xlim.ls=list()
    markers.x=list()
    weight=list()
    for (x in TRIALS) {    
        TRIAL=get.trial(x, a)
        Sys.setenv("TRIAL"=TRIAL)
        COR = ifelse (x=="moderna_real","D57","D29IncludeNotMolecConfirmedstart1")
        # key to have local = T
        source(here::here("..", "_common.R"), local=T)
        
        # uloq censoring    
        tmp="Day"%.%config.cor$tpeak %.% a
        dat.mock[[tmp]] <- ifelse(dat.mock[[tmp]] > log10(uloqs[a]), log10(uloqs[a]), dat.mock[[tmp]])
        
        dat.vac.seroneg=subset(dat.mock, Trt==1 & ph1)
        xlim=range(dat.vac.seroneg[[tmp]], log10(llods[a]/2), na.rm=T)
        delta=(xlim[2]-xlim[1])/20     
        xlim.ls[[x]]=c(xlim[1]-delta, xlim[2]+delta)
        
        markers.x[[x]]=dat.vac.seroneg[[tmp]][dat.vac.seroneg$ph2]
        weight[[x]]=dat.vac.seroneg[["wt"]][dat.vac.seroneg$ph2]
    }    
    xlim=c(min(sapply(xlim.ls, function(x) x[1])), max(sapply(xlim.ls, function(x) x[2])))
    myprint(xlim)
    
    mypdf(file=paste0("output/meta/meta_controlled_ve_curves",ifelse(log=="","","log"),"_",file.name,"_",a), width=5.2, height=5.2)
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
        
            shown=risks$marker>=ifelse(x=="moderna_real",log10(10),quantile(markers.x[[x]], 2.5/100, na.rm=T)) & risks$marker<=quantile(markers.x[[x]], 1-2.5/100, na.rm=T)
            y=t(rbind(est, ci.band))[shown,]
            if(log=="y") y=-log(1-y)
            mymatplot(risks$marker[shown], y, type="l", lty=c(1,3,3), lwd=2.5, make.legend=F, col=cols[x], ylab=paste0("Controlled VE"), xlab=labels.assays.short[a]%.%" (=s)", 
                #main=paste0(labels.assays.long["Day"%.%tpeak,a]),
                ylim=ylim, xlim=xlim, yaxt="n", xaxt="n", draw.x.axis=F, add=x!=TRIALS[1])
            draw.x.axis.cor(xlim, NA)
            if (log=="") {
                yat=seq(-1,1,by=.1)
                axis(side=2,at=yat,labels=(yat*100)%.%"%")            
            } else {
                yat=c(seq(0,.90,by=.1),.95)
                axis(side=2,at=-log(1-yat),labels=(yat*100)%.%"%")            
            }
        
            # add histogram
            #  par(new=TRUE) #this changes ylim, so we cannot use it in this loop
            # first call hist to get breaks, then call weighted.hist
            tmp.1=hist(markers.x[[x]],breaks=ifelse(x=="moderna_real",25,15),plot=F)  # 15 is treated as a suggestion and the actual number of breaks is determined by pretty()
            tmp=weighted.hist(markers.x[[x]],weight[[x]], breaks=tmp.1$breaks, plot=F)
            attr(tmp,"class")="histogram" 
            if (log=="") {
                tmp$density=tmp$density/hist.shrink[a] # so that it will fit vertically
            } else{
                tmp$density=tmp$density/hist.shrink[a]*3 # so that it will fit vertically
            }
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
        mylegend(x=ifelse(log=="",6,1), col=cols[c(TRIALS, if(include.az) "AZ-COV002")], legend=legend, lty=1, lwd=2, cex=.7)
    
    dev.off()    
    
}


# COVE + ENSEMBLE + AZ
for (a in c("pseudoneutid50","bindSpike","bindRBD")) {
    draw.ve.curves(a, TRIALS=c("moderna_real", "janssen_pooled_real"), file.name="1", include.az=T)
    draw.ve.curves(a, TRIALS=c("moderna_real", "janssen_pooled_real"), file.name="1", include.az=T, log="y")
}


# COVE + ENSEMBLE regions
for (a in c("pseudoneutid50","bindSpike","bindRBD")) {
    draw.ve.curves(a, TRIALS=c("moderna_real", "janssen_na_real", "janssen_la_real", "janssen_sa_real"), file.name="2", include.az=F)
    draw.ve.curves(a, TRIALS=c("moderna_real", "janssen_na_real", "janssen_la_real", "janssen_sa_real"), file.name="2", include.az=F, log="y")
}


# ENSEMBLE regions
for (a in c("pseudoneutid50","bindSpike","bindRBD","ADCP")) {
    draw.ve.curves(a, TRIALS=c("janssen_na_real", "janssen_la_real", "janssen_sa_real"), file.name="3", include.az=F)
    draw.ve.curves(a, TRIALS=c("janssen_na_real", "janssen_la_real", "janssen_sa_real"), file.name="3", include.az=F, log="y")
}

# COVE + ENSEMBLE/US + AZ
for (a in c("pseudoneutid50","bindSpike","bindRBD")) {
    draw.ve.curves(a, TRIALS=c("moderna_real", "janssen_na_real"), file.name="4", include.az=T)
    draw.ve.curves(a, TRIALS=c("moderna_real", "janssen_na_real"), file.name="4", include.az=T, log="y")
}

# COVE + ENSEMBLE/US
for (a in c("pseudoneutid50","bindSpike","bindRBD")) {
    draw.ve.curves(a, TRIALS=c("moderna_real", "janssen_na_real"), file.name="5", include.az=F)
    draw.ve.curves(a, TRIALS=c("moderna_real", "janssen_na_real"), file.name="5", include.az=F, log="y")
}
