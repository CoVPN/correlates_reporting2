renv::activate(project = here::here(".."))    
    # There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
    if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))
    
library(kyotil)
library(Hmisc)
library(marginalizedRisk)
library(tools) # toTitleCase
library(xtable) # this is a dependency of kyotil
source(here::here("..", "_common.R"))
source(here::here("code", "params.R"))
Sys.setenv(VERBOSE=1)
print("collate ...")



# reading in data for COVE and ENSEMBLE
dat.vac.seroneg.bAb=load.data("janssen_pooled_real", "D29IncludeNotMolecConfirmedstart1")
dat.vac.seroneg.id50=load.data("janssen_pooled_realPsV", "D29IncludeNotMolecConfirmedstart1")
dat.vac.seroneg.adcp=load.data("janssen_pooled_realADCP", "D29IncludeNotMolecConfirmedstart1")
stopifnot(all(dat.vac.seroneg.id50$ptid==dat.vac.seroneg.bAb$ptid))
stopifnot(all(dat.vac.seroneg.id50$ptid==dat.vac.seroneg.adcp$ptid))
# combine markers into one data frame
dat.ense.1=cbind(dat.vac.seroneg.bAb, 
    Day29pseudoneutid50=dat.vac.seroneg.id50$Day29pseudoneutid50, 
    Day29pseudoneutid50cat=dat.vac.seroneg.id50$Day29pseudoneutid50cat, 
    Day29ADCP=dat.vac.seroneg.adcp$Day29ADCP, 
    Day29ADCPcat=dat.vac.seroneg.adcp$Day29ADCPcat)

dat.ense.0=load.data("janssen_pooled_real", "D29IncludeNotMolecConfirmedstart1", trt=0)
    
dat.cove.1=load.data("moderna_real", "D57")
dat.cove.0=load.data("moderna_real", "D57", trt=0)



###################################################################################################
# make tables of HR that contain all markers

for (region in c("pooled","na","la","sa")) {
    if(verbose) myprint(region)

    trials=c("janssen_"%.%region%.%"_real", "janssen_"%.%region%.%"_realPsV", "janssen_"%.%region%.%"_realADCP")
    
    for (COR in c("D29IncludeNotMolecConfirmed", "D29IncludeNotMolecConfirmedstart1")) {
        if(verbose) myprint(COR)
        out=lapply(trials, function (trial) {
            print(trial)
            load (file=paste0(here::here("output", trial, COR), "/coxph_slopes.Rdata")) 
            list(tab.cont, tab.cat, save.s.1, save.s.2, pvals.adj)
        })
    
        pvals.adj=do.call(rbind, sapply(out, function (res) res[[5]]))
        pvals.adj.fdr=p.adjust(pvals.adj[,1], method="fdr")
        pvals.adj.hol=p.adjust(pvals.adj[,1], method="holm")
        
        tab.cont=do.call(rbind, sapply(out, function (res) res[[1]]))
        p.1=formatDouble(pvals.adj.hol[startsWith(names(pvals.adj.hol), "cont.")], 3, remove.leading0=F); p.1=sub("0.000","<0.001",p.1)
        p.2=formatDouble(pvals.adj.fdr[startsWith(names(pvals.adj.fdr), "cont.")], 3, remove.leading0=F); p.2=sub("0.000","<0.001",p.2)
        tab.cont[,"p.1"]=p.1
        tab.cont[,"p.2"]=p.2
        
        tab.cat=do.call(rbind, sapply(out, function (res) res[[2]]))
        overall.p.1=formatDouble(pvals.adj.hol[startsWith(names(pvals.adj.hol), "tri.")], 3, remove.leading0=F); p.1=sub("0.000","<0.001",p.1)
        overall.p.2=formatDouble(pvals.adj.fdr[startsWith(names(pvals.adj.fdr), "tri.")], 3, remove.leading0=F); p.2=sub("0.000","<0.001",p.2)
        tab.cat[,"overall.p.1"]=c(rbind(overall.p.1, NA,NA))
        tab.cat[,"overall.p.2"]=c(rbind(overall.p.2, NA,NA))
                
        mytex(tab.cont, file.name="CoR_univariable_svycoxph_pretty_ENSEMBLE_"%.%region%.%"_"%.%COR, align="c", include.colnames = F, save2input.only=T, input.foldername=here::here("output"),
            col.headers=paste0("\\hline\n 
                 \\multicolumn{1}{l}{} & \\multicolumn{1}{c}{No. cases /}   & \\multicolumn{2}{c}{HR per 10-fold incr.}                     & \\multicolumn{1}{c}{P-value}   & \\multicolumn{1}{c}{q-value}   & \\multicolumn{1}{c}{FWER} \\\\ 
                 \\multicolumn{1}{l}{Immunologic Marker}            & \\multicolumn{1}{c}{No. at-risk**} & \\multicolumn{1}{c}{Pt. Est.} & \\multicolumn{1}{c}{95\\% CI} & \\multicolumn{1}{c}{(2-sided)} & \\multicolumn{1}{c}{***} & \\multicolumn{1}{c}{} \\\\ 
                 \\hline\n 
            ")
        )    
        
        mytex(tab.cat, file.name="CoR_univariable_svycoxph_cat_pretty_ENSEMBLE_"%.%region%.%"_"%.%COR, align="c", include.colnames = F, save2input.only=T, input.foldername=here::here("output"),
            col.headers=paste0("\\hline\n 
                 \\multicolumn{1}{l}{} & \\multicolumn{1}{c}{Tertile}   & \\multicolumn{1}{c}{No. cases /}   & \\multicolumn{1}{c}{Attack}   & \\multicolumn{2}{c}{Haz. Ratio}                     & \\multicolumn{1}{c}{P-value}   & \\multicolumn{1}{c}{Overall P-}      & \\multicolumn{1}{c}{Overall q-}   & \\multicolumn{1}{c}{Overall} \\\\ 
                 \\multicolumn{1}{l}{Immunologic Marker}            & \\multicolumn{1}{c}{}          & \\multicolumn{1}{c}{No. at-risk**} & \\multicolumn{1}{c}{rate}   & \\multicolumn{1}{c}{Pt. Est.} & \\multicolumn{1}{c}{95\\% CI} & \\multicolumn{1}{c}{(2-sided)} & \\multicolumn{1}{c}{value***} & \\multicolumn{1}{c}{value $\\dagger$} & \\multicolumn{1}{c}{FWER} \\\\ 
                 \\hline\n 
            "),        
            add.to.row=list(list(nrow(tab.cat)), # insert at the beginning of table, and at the end of, say, the first table
                c(paste0(" \n \\multicolumn{8}{l}{} \\\\ \n", 
                          "\n \\multicolumn{2}{l}{Placebo} & ", 
                         out[[1]][[3]], "&",  
                         out[[1]][[4]], "&",  
                         "\\multicolumn{4}{l}{}  \\\\ \n")
                 )
            )
        )    
        
    }
    
}





###################################################################################################
# correlations between markers

# ENSEMBLE
mypdf(mfrow=c(2,2), file="output/ensemble_corplot")
with(dat.ense.1, {
    corplot(Day29pseudoneutid50, Day29ADCP, xlab="ID50", ylab="ADCP")
    corplot(Day29pseudoneutid50, Day29bindRBD, xlab="ID50", ylab="bAb RBD")
    corplot(Day29bindSpike, Day29bindRBD, xlab="bAb Spike", ylab="bAb RBD")
    corplot(Day29bindRBD, Day29ADCP, xlab="bAb RBD", ylab="ADCP")
})
dev.off()

with(dat.ense.1, table(Day29pseudoneutid50cat, Day29ADCPcat))
with(dat.ense.1, table(Day29pseudoneutid50cat, Day29bindRBDcat))
with(dat.ense.1, table(Day29bindSpikecat, Day29bindRBDcat))
with(dat.ense.1, table(Day29bindRBDcat, Day29ADCPcat))

# COVE
mypdf(mfrow=c(2,2), file="output/cove_corplot")
with(dat.cove.1, {
    corplot(Day57pseudoneutid50, Day57bindSpike, xlab="ID50", ylab="bAb Spike")
    corplot(Day57pseudoneutid50, Day57bindRBD, xlab="ID50", ylab="bAb RBD")
    corplot(Day57bindSpike, Day57bindRBD, xlab="bAb Spike", ylab="bAb RBD")
    corplot(Day57pseudoneutid50, Day57pseudoneutid80, xlab="ID50", ylab="ID80")
})
dev.off()

with(dat.cove.1, table(Day57pseudoneutid50cat, Day57bindSpikecat))
with(dat.cove.1, table(Day57pseudoneutid50cat, Day57bindRBDcat))
with(dat.cove.1, table(Day57bindSpikecat, Day57bindRBDcat))
with(dat.cove.1, table(Day57pseudoneutid50cat, Day57pseudoneutid80cat))



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
