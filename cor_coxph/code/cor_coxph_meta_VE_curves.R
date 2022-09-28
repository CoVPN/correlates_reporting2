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
source(here::here("code", "params.R")) # get.trial
Sys.setenv(VERBOSE=1)
print("meta ...")
ve.az=read.csv("../data_clean/AZChAd26UKphase3FengetalCorrelates.csv")

# path for figures and tables etc
save.results.to = here::here("output")
if (!dir.exists(save.results.to))  dir.create(save.results.to)
save.results.to = paste0(here::here("output"), "/meta");
if (!dir.exists(save.results.to))  dir.create(save.results.to)


# draw VE curves for several trials, the same marker
# TRIALS is a subset of all.trials
# a is an assay
draw.ve.curves=function(a, TRIALS, file.name, include.az=FALSE, log="", save2file=T, add.hist=T, show.cb=T) {
#a="bindSpike"; TRIALS=c("moderna_real", "janssen_pooled_real", "prevent19"); file.name=1; include.az=T; log=""
    myprint(a)
    
    transf=if(log=="y") function(y) -log(1-y) else identity 
    
    ylim=if(log=="y") transf(c(0,.98)) else c(0, 1) 
    hist.shrink=1/c(ADCP=2,pseudoneutid50=1.2,bindSpike=1.3,bindRBD=1.3)
    
    all.trials=c("moderna_real", "janssen_pooled_real", "janssen_na_real",    "janssen_la_real",    "janssen_sa_real",    "AZ-COV002", "prevent19",      "azd1222")
    studies   =c("Moderna COVE", "Janssen ENSEMBLE",    "Janssen ENSEMBLE US","Janssen ENSEMBLE LA","Janssen ENSEMBLE SA","AZCOV002",  "NVX PREVENT-19", "AZD1222")
    cols      =c("purple",       "green",               "green",              "olivedrab3",         "darkseagreen4",      "orange",    "cyan",           "tan")
    names(studies)=all.trials
    names(cols)=all.trials
    hist.col.ls=lapply(cols, function(col) {hist.col <- c(col2rgb(col)); rgb(hist.col[1], hist.col[2], hist.col[3], alpha=255*0.5, maxColorValue=255)})
    
    .subset=match(TRIALS, all.trials)
    
    
    ## get markers data
    ## get xlim by combining trials. Do ENSEMBLE first because we want to source _common.R for moderna last so that we get the proper lloxs
    ## source _commom.R
    xlim.ls=list()
    markers.x=list()
    weight=list()
    for (x in TRIALS) {    
        TRIAL=get.trial(x, a)
        COR = switch(x, moderna_real="D57",
            janssen_pooled_real="D29IncludeNotMolecConfirmedstart1", janssen_na_real="D29IncludeNotMolecConfirmedstart1", 
            janssen_la_real="D29IncludeNotMolecConfirmedstart1", janssen_sa_real="D29IncludeNotMolecConfirmedstart1",
            prevent19="D35", 
            azd1222="D57",
            stop("wrong trial label for COR")
        )
        # key to have local = T
        Sys.setenv("TRIAL"=TRIAL)
        source(here::here("..", "_common.R"), local=T)
        
        # uloq censoring    
        tmp="Day"%.%config.cor$tpeak %.% a
        dat.mock[[tmp]] <- ifelse(dat.mock[[tmp]] > log10(uloqs[a]), log10(uloqs[a]), dat.mock[[tmp]])
        
        dat.vac.seroneg=subset(dat.mock, Trt==1 & ph1) # should not need  & Bserostatus==0
        xlim=range(dat.vac.seroneg[[tmp]], log10(lloxs[a]/2), na.rm=T)
        delta=(xlim[2]-xlim[1])/20     
        xlim.ls[[x]]=c(xlim[1]-delta, xlim[2]+delta)
        
        markers.x[[x]]=dat.vac.seroneg[[tmp]][dat.vac.seroneg$ph2]
        weight[[x]]=dat.vac.seroneg[["wt"]][dat.vac.seroneg$ph2]
    }    
    xlim=c(min(sapply(xlim.ls, function(x) x[1])), max(sapply(xlim.ls, function(x) x[2])))
    myprint(xlim)
    
    
    if(save2file) mypdf(file=paste0("output/meta/meta_controlled_ve_curves",ifelse(log=="","","log"),"_",file.name,"_",a), width=5.2, height=5.2)
        par(las=1, cex.axis=0.9, cex.lab=1)# axis label orientation
        
        # need several variables from sourcing _common.R: lloxs, labels.assays, draw.x.axis.cor
        overall.ve.ls=list()
        for (x in TRIALS) {    
            myprint(x, a)
            TRIAL=get.trial(x, a)
            COR = switch(x, moderna_real="D57",
                janssen_pooled_real="D29IncludeNotMolecConfirmedstart1", janssen_na_real="D29IncludeNotMolecConfirmedstart1", 
                janssen_la_real="D29IncludeNotMolecConfirmedstart1", janssen_sa_real="D29IncludeNotMolecConfirmedstart1",
                prevent19="D35", 
                azd1222="D57",
                stop("wrong trial label for COR")
            )
            load(here::here("output", TRIAL, COR, "marginalized.risk.no.marker.Rdata"))
            load(here::here("output", TRIAL, COR, "marginalized.risk.Rdata"))
            # key to have local = T
            Sys.setenv("TRIAL"=TRIAL)
            source(here::here("..", "_common.R"), local=T)
            tmp="Day"%.%config.cor$tpeak %.% a            
            risks=get("risks.all.1")[[tmp]] 
            if (is.null(risks)) {
                myprint(x, a)
                stop("risks is NULL")
            }
            #risks=get(ifelse(x=="prevent19","risks.all.3","risks.all.1"))[[a]] # .3 is trichotomized
            
            est = 1 - risks$prob/res.plac.cont["est"]
            boot = 1 - t( t(risks$boot)/res.plac.cont[2:(1+ncol(risks$boot))] )                         
            ci.band=apply(boot, 1, function (x) quantile(x, c(.025,.975)))    
            
#            # if use .3 earlier
#            if(x=="prevent19") {
#                # to make it plottable by mymatplot
#                cutpoints=as.numeric(strsplit(risks$marker[2],"\\(|\\,|\\]")[[1]][-1])
#                risks$marker=c(quantile(markers.x[[x]], 2.5/100, na.rm=T), cutpoints[1], cutpoints[1], cutpoints[2], cutpoints[2], quantile(markers.x[[x]], 1-2.5/100, na.rm=T))
#                est=rep(est, each=2)
#                ci.band=rep.matrix(ci.band, each=2, by.row=F)
#            }            
        
            shown=risks$marker>=wtd.quantile(markers.x[[x]], weight[[x]], 2.5/100) & risks$marker<=wtd.quantile(markers.x[[x]], weight[[x]], 1-2.5/100)
            #shown=risks$marker>=ifelse(x=="moderna_real",log10(10),quantile(markers.x[[x]], 2.5/100, na.rm=T)) & risks$marker<=quantile(markers.x[[x]], 1-2.5/100, na.rm=T)
            mymatplot(risks$marker[shown], transf(t(rbind(est, if(show.cb) ci.band))[shown,]), type="l", lty=c(1,3,3), lwd=2.5, make.legend=F, col=cols[x], ylab=paste0("Controlled VE against COVID-19"), xlab=labels.assays.short[a]%.%" (=s)", 
                #main=paste0(labels.assays.long["Day"%.%tpeak,a]),
                ylim=ylim, xlim=xlim, yaxt="n", xaxt="n", draw.x.axis=F, add=x!=TRIALS[1])
            draw.x.axis.cor(xlim, NA)
            if (log=="") {
                yat=seq(-1,1,by=.1)
                axis(side=2,at=yat,labels=(yat*100)%.%"%")            
            } else {
                yat=c(seq(0,.90,by=.1),.95)
                axis(side=2,at=transf(yat),labels=(yat*100)%.%"%")            
            }
        
            img.dat=cbind(risks$marker[shown], transf(t(rbind(est, ci.band))[shown,]))
            mywrite.csv(img.dat, file=paste0("output/meta/meta_controlled_ve_curves",ifelse(log=="","","log"),"_",file.name, "_",a, "_",x))
            
            # add histogram
            #  par(new=TRUE) #this changes ylim, so we cannot use it in this loop
            if(add.hist) {
                tmp=get.marker.histogram(markers.x[[x]], weight[[x]], x)
                if (log=="") {
                    tmp$density=tmp$density/hist.shrink[a] # so that it will fit vertically
                } else{
                    tmp$density=tmp$density/hist.shrink[a]*3 # so that it will fit vertically
                }            
                plot(tmp,col=hist.col.ls[[x]],axes=F,labels=F,border=0,freq=F,add=T) 
            }
                        
            overall.ve.ls[[x]]=overall.ve
        }        
    
        # add az curve
        if(include.az) {
            lines(log10(ve.az[[a]]),        transf(ve.az$VE/100), col=cols["AZ-COV002"], lwd=2.5)
            if(show.cb) lines(log10(ve.az[[a%.%"LL"]]), transf(ve.az$VE/100), col=cols["AZ-COV002"], lwd=2.5, lty=3)
            if(show.cb) lines(log10(ve.az[[a%.%"UL"]]), transf(ve.az$VE/100), col=cols["AZ-COV002"], lwd=2.5, lty=3)
        }
    
        # legend
#        legend=paste0(studies[TRIALS], ", ", sapply(overall.ve.ls, function(x) formatDouble(x*100,1)%.%"%"))
#        if (include.az) legend=c(legend, "AZCOV002, 71.8%")
        # remove VE from legend
        legend=studies[TRIALS]
        if (include.az) legend=c(legend, "AZCOV002")
        mylegend(x=ifelse(log=="",6,1), col=cols[c(TRIALS, if(include.az) "AZ-COV002")], legend=legend, lty=1, lwd=2, cex=.7)
    
    if(save2file) dev.off()    
    
    return (list(markers.x, weight))
    
}


# draw VE curves for several markers
draw.ve.curves.aa=function(aa, TRIALS, file.name, log="") {
#aa=c("pseudoneutid50","pseudoneutid50la"); TRIALS=c("janssen_la_real"); file.name=tmp; log=""
    
    all.trials=c("moderna_real", "janssen_pooled_real", "janssen_na_real",    "janssen_la_real",    "janssen_sa_real",    "AZ-COV002", "prevent19",      "azd1222")
    studies   =c("Moderna COVE", "Janssen ENSEMBLE",    "Janssen ENSEMBLE US","Janssen ENSEMBLE LA","Janssen ENSEMBLE SA","AZCOV002",  "NVX PREVENT-19", "AZD1222")
    names(studies)=all.trials
    
    if (length(TRIALS)==1) {
        same.trial=T
        # use assay to index colors
        cols      =c("purple",       "green",               "green",              "olivedrab3",         "darkseagreen4",      "orange",    "cyan",           "tan")
        cols = cols[1:length(aa)]
        names(cols)=aa
        hist.col.ls=lapply(cols, function(col) {hist.col <- c(col2rgb(col)); rgb(hist.col[1], hist.col[2], hist.col[3], alpha=255*0.5, maxColorValue=255)})
        # make TRIALS the same length as aa
        TRIALS=rep(TRIALS, length(aa))
    } else if(length(TRIALS)==length(aa)) {
        same.trial=F
        # use trials to index colors
        cols      =c("purple",       "green",               "green",              "olivedrab3",         "darkseagreen4",      "orange",    "cyan",           "tan")
        names(cols)=all.trials
        hist.col.ls=lapply(cols, function(col) {hist.col <- c(col2rgb(col)); rgb(hist.col[1], hist.col[2], hist.col[3], alpha=255*0.5, maxColorValue=255)})
    } else stop("TRIALS has to be either of length 1 or same length as aa")
    
    x=TRIALS
    
    transf=if(log=="y") function(y) -log(1-y) else identity 
    
    ylim=if(log=="y") transf(c(0,.98)) else c(0, 1) 
    hist.shrink=1/c(ADCP=2,pseudoneutid50=1.2,bindSpike=1.3,bindRBD=1.3)
    
    .subset=match(TRIALS, all.trials)    
    
    ## get markers data
    ## get xlim by combining trials. Do ENSEMBLE first because we want to source _common.R for moderna last so that we get the proper lloxs
    ## source _commom.R
    xlim.ls=list()
    markers.x=list()
    weight=list()
    for (i in 1:length(aa)) {  
        a=aa[i]  
        TRIAL=get.trial(x[i], a)
        COR = switch(x[i], 
            moderna_real="D57",
            janssen_pooled_real="D29IncludeNotMolecConfirmedstart1", janssen_na_real="D29IncludeNotMolecConfirmedstart1", janssen_la_real="D29IncludeNotMolecConfirmedstart1", janssen_sa_real="D29IncludeNotMolecConfirmedstart1",
            prevent19="D35", 
            azd1222="D57",
            stop("wrong trial label for COR")
        )
        # key to have local = T
        Sys.setenv("TRIAL"=TRIAL)
        source(here::here("..", "_common.R"), local=T)
        
        # uloq censoring    
        tmp="Day"%.%config.cor$tpeak %.% a
        dat.mock[[tmp]] <- ifelse(dat.mock[[tmp]] > log10(uloqs[a]), log10(uloqs[a]), dat.mock[[tmp]])
        
        dat.vac.seroneg=subset(dat.mock, Trt==1 & ph1) # should not need  & Bserostatus==0
        xlim=range(dat.vac.seroneg[[tmp]], log10(lloxs[a]/2), na.rm=T)
        delta=(xlim[2]-xlim[1])/20     
        xlim.ls[[a]]=c(xlim[1]-delta, xlim[2]+delta)
        
        markers.x[[a]]=dat.vac.seroneg[[tmp]][dat.vac.seroneg$ph2]
        weight[[a]]=dat.vac.seroneg[["wt"]][dat.vac.seroneg$ph2]
    }    
    xlim=c(min(sapply(xlim.ls, function(x) x[1])), max(sapply(xlim.ls, function(x) x[2])))
    myprint(xlim)
    
    mypdf(file=paste0("output/meta/meta_controlled_ve_curves",ifelse(log=="","","log"),"_",file.name), width=5.2, height=5.2)
        par(las=1, cex.axis=0.9, cex.lab=1)# axis label orientation
        
        # need several variables from sourcing _common.R: lloxs, labels.assays, draw.x.axis.cor
        overall.ve.ls=list()
        for (i in 1:length(aa)) {  
            a=aa[i]  
            TRIAL=get.trial(x[i], a)
            COR = switch(x[i], moderna_real="D57",
                janssen_pooled_real="D29IncludeNotMolecConfirmedstart1", janssen_na_real="D29IncludeNotMolecConfirmedstart1", 
                janssen_la_real="D29IncludeNotMolecConfirmedstart1", janssen_sa_real="D29IncludeNotMolecConfirmedstart1",
                prevent19="D35", 
                azd1222="D57",
                stop("wrong trial label for COR")
            )
            load(here::here("output", TRIAL, COR, "marginalized.risk.no.marker.Rdata"))
            load(here::here("output", TRIAL, COR, "marginalized.risk.Rdata"))
    
            # key to have local = T
            Sys.setenv("TRIAL"=TRIAL)
            source(here::here("..", "_common.R"), local=T)
            tmp="Day"%.%config.cor$tpeak %.% a            
            risks=get("risks.all.1")[[tmp]] 
            #risks=get(ifelse(x=="prevent19","risks.all.3","risks.all.1"))[[a]] # .3 is trichotomized
            
            est = 1 - risks$prob/res.plac.cont["est"]
            boot = 1 - t( t(risks$boot)/res.plac.cont[2:(1+ncol(risks$boot))] )                         
            ci.band=apply(boot, 1, function (x) quantile(x, c(.025,.975)))    
            
#            # if use .3 earlier
#            if(x=="prevent19") {
#                # to make it plottable by mymatplot
#                cutpoints=as.numeric(strsplit(risks$marker[2],"\\(|\\,|\\]")[[1]][-1])
#                risks$marker=c(quantile(markers.x[[a]], 2.5/100, na.rm=T), cutpoints[1], cutpoints[1], cutpoints[2], cutpoints[2], quantile(markers.x[[a]], 1-2.5/100, na.rm=T))
#                est=rep(est, each=2)
#                ci.band=rep.matrix(ci.band, each=2, by.row=F)
#            }            
        
            shown=risks$marker>=wtd.quantile(markers.x[[a]], weight[[a]], 2.5/100) & risks$marker<=wtd.quantile(markers.x[[a]], weight[[a]], 1-2.5/100)
            #shown=risks$marker>=ifelse(x=="moderna_real",log10(10),quantile(markers.x[[a]], 2.5/100, na.rm=T)) & risks$marker<=quantile(markers.x[[a]], 1-2.5/100, na.rm=T)
            mymatplot(risks$marker[shown], transf(t(rbind(est, ci.band))[shown,]), type="l", lty=c(1,3,3), lwd=2.5, make.legend=F, col=if(same.trial) cols[a] else cols[TRIAL], ylab=paste0("Controlled VE against COVID-19"), 
                xlab=labels.assays.short[a]%.%" (=s)", 
                #main=paste0(labels.assays.long["Day"%.%tpeak,a]),
                ylim=ylim, xlim=xlim, yaxt="n", xaxt="n", draw.x.axis=F, 
                add=a!=aa[1])
            draw.x.axis.cor(xlim, NA)
            if (log=="") {
                yat=seq(-1,1,by=.1)
                axis(side=2,at=yat,labels=(yat*100)%.%"%")            
            } else {
                yat=c(seq(0,.90,by=.1),.95)
                axis(side=2,at=transf(yat),labels=(yat*100)%.%"%")            
            }
            
            img.dat=cbind(risks$marker[shown], transf(t(rbind(est, ci.band))[shown,]))
            mywrite.csv(img.dat, file=paste0("output/meta/meta_controlled_ve_curves",ifelse(log=="","","log"),"_",file.name, "_",i))
        
            # add histogram
            #  par(new=TRUE) #this changes ylim, so we cannot use it in this loop
            tmp=get.marker.histogram(markers.x[[a]], weight[[a]], x[i], markers.x[[aa[1]]])
            if (log=="") {
                tmp$density=tmp$density/hist.shrink[aa[1]] # so that it will fit vertically
            } else{
                tmp$density=tmp$density/hist.shrink[aa[1]]*3 # so that it will fit vertically
            }            
            plot(tmp,col=if(same.trial) hist.col.ls[[a]] else hist.col.ls[[TRIAL]],axes=F,labels=F,border=0,freq=F,add=T) 
            
            overall.ve.ls[[a]]=overall.ve
        }        
    
        # legend
        legend=aa
        mylegend(x=ifelse(log=="",6,1), col=if(same.trial) cols[aa] else cols[TRIALS], legend=legend, lty=1, lwd=2, cex=.7)
    
    dev.off()    
    
}


# for NEJM perspective
a="pseudoneutid50"
res=draw.ve.curves(a, TRIALS=c("moderna_real", "janssen_na_real", "prevent19", "azd1222"), file.name="9", include.az=T, log="y", save2file=F, add.hist=F, show.cb=F)





# for ENSEMBLE manuscript 1 
draw.ve.curves.aa(aa=c("pseudoneutid50","pseudoneutid50sa","pseudoneutid50la"), TRIALS=c("janssen_na_real", "janssen_sa_real", "janssen_la_real"), file.name="circ")
draw.ve.curves.aa(aa=c("pseudoneutid50","pseudoneutid50sa","pseudoneutid50la"), TRIALS=c("janssen_na_real", "janssen_sa_real", "janssen_la_real"), file.name="circ", log="y")

# 
draw.ve.curves.aa(aa=c("pseudoneutid50","pseudoneutid50la"), TRIALS="janssen_la_real", file.name="janssen_la_real")
draw.ve.curves.aa(aa=c("pseudoneutid50","pseudoneutid50sa"), TRIALS="janssen_sa_real", file.name="janssen_sa_real")
draw.ve.curves.aa(aa=c("pseudoneutid50","pseudoneutid50la"), TRIALS="janssen_la_real", file.name="janssen_la_real", log="y")
draw.ve.curves.aa(aa=c("pseudoneutid50","pseudoneutid50sa"), TRIALS="janssen_sa_real", file.name="janssen_sa_real", log="y")


# COVE + ENSEMBLE/US + AZ + PREVENT19 + AZD1222
for (a in c("pseudoneutid50","bindSpike")) {
    draw.ve.curves(a, TRIALS=c("moderna_real", "janssen_na_real", "prevent19", "azd1222"), file.name="9", include.az=T)
    draw.ve.curves(a, TRIALS=c("moderna_real", "janssen_na_real", "prevent19", "azd1222"), file.name="9", include.az=T, log="y")
}

# COVE + ENSEMBLE/US + AZ + AZD1222
for (a in c("pseudoneutid50","bindSpike")) {
    draw.ve.curves(a, TRIALS=c("moderna_real", "janssen_na_real", "azd1222"), file.name="10", include.az=T)
    draw.ve.curves(a, TRIALS=c("moderna_real", "janssen_na_real", "azd1222"), file.name="10", include.az=T, log="y")
}

# AZ + AZD1222
for (a in c("pseudoneutid50","bindSpike")) {
    draw.ve.curves(a, TRIALS=c("azd1222"), file.name="11", include.az=T)
    draw.ve.curves(a, TRIALS=c("azd1222"), file.name="11", include.az=T, log="y")
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
for (a in c("pseudoneutid50","bindSpike","bindRBD")) {
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

# COVE + ENSEMBLE regions + AZ
for (a in c("pseudoneutid50","bindSpike","bindRBD")) {
    draw.ve.curves(a, TRIALS=c("moderna_real", "janssen_na_real", "janssen_la_real", "janssen_sa_real"), file.name="6", include.az=T)
    draw.ve.curves(a, TRIALS=c("moderna_real", "janssen_na_real", "janssen_la_real", "janssen_sa_real"), file.name="6", include.az=T, log="y")
}

# COVE + ENSEMBLE/US + AZ + PREVENT19
for (a in c("bindSpike", "pseudoneutid50")) {
    draw.ve.curves(a, TRIALS=c("moderna_real", "janssen_na_real", "prevent19"), file.name="7", include.az=T)
    draw.ve.curves(a, TRIALS=c("moderna_real", "janssen_na_real", "prevent19"), file.name="7", include.az=T, log="y")
}

# COVE + ENSEMBLE/US + PREVENT19
for (a in c("bindSpike", "pseudoneutid50")) {
    draw.ve.curves(a, TRIALS=c("moderna_real", "janssen_na_real", "prevent19"), file.name="8", include.az=F)
    draw.ve.curves(a, TRIALS=c("moderna_real", "janssen_na_real", "prevent19"), file.name="8", include.az=F, log="y")
}
