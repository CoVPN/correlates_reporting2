# Sys.setenv(TRIAL = "janssen_pooled_realbAb") # just so that _common.R can run
renv::activate(project = here::here(".."))    
    # There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
    if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))
    
library(kyotil)
library(Hmisc)
library(marginalizedRisk)
library(tools) # toTitleCase
library(survey)
library(xtable) # this is a dependency of kyotil
Sys.setenv("TRIAL"="janssen_pooled_realbAb") # value does not matter since we just need to load the common functions in _common.R
source(here::here("..", "_common.R"))
source(here::here("code", "params.R"))
Sys.setenv(VERBOSE=1)
print("collate ...")

# path for figures and tables etc
save.results.to = here::here("output")
if (!dir.exists(save.results.to))  dir.create(save.results.to)
save.results.to = paste0(here::here("output"), "/meta");
if (!dir.exists(save.results.to))  dir.create(save.results.to)


# reading in data for COVE and ENSEMBLE
dat.vac.seroneg.bAb=load.data("janssen_pooled_realbAb", "D29IncludeNotMolecConfirmedstart1")
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

dat.ense.0=load.data("janssen_pooled_realbAb", "D29IncludeNotMolecConfirmedstart1", trt=0)
    
dat.cove.1=load.data("moderna_real", "D57")
dat.cove.0=load.data("moderna_real", "D57", trt=0)

dat.vac.seroneg.id50.na=load.data("janssen_na_realPsV", "D29IncludeNotMolecConfirmedstart1")
dat.pla.seroneg.id50.na=load.data("janssen_na_realPsV", "D29IncludeNotMolecConfirmedstart1", trt=0)
dat.vac.seroneg.id50.la=load.data("janssen_la_realPsV", "D29IncludeNotMolecConfirmedstart1")
dat.pla.seroneg.id50.la=load.data("janssen_la_realPsV", "D29IncludeNotMolecConfirmedstart1", trt=0)
dat.vac.seroneg.id50.sa=load.data("janssen_sa_realPsV", "D29IncludeNotMolecConfirmedstart1")
dat.pla.seroneg.id50.sa=load.data("janssen_sa_realPsV", "D29IncludeNotMolecConfirmedstart1", trt=0)


###################################################################################################
# make tables of HR that contain all markers

for (region in c("pooled","na","la","sa")) {
        
for (ind in 1:ifelse(region=="pooled",2,1)) {
# 1 include ADCP and 2 does not
    
    if(verbose) myprint(region)
    
    trials=c("janssen_"%.%region%.%"_realbAb", "janssen_"%.%region%.%"_realPsV", if(ind==1) "janssen_"%.%region%.%"_realADCP") 
        
    for (COR in c("D29IncludeNotMolecConfirmed", "D29IncludeNotMolecConfirmedstart1")) {
        if(verbose) myprint(COR)
        out=lapply(trials, function (trial) {
            print(trial)
            load (file=paste0(here::here("output", trial, COR), "/coxph_slopes.Rdata")) 
            list(tab.cont, tab.cat, save.s.1, save.s.2, pvals.adj, tab.cont.scaled)
        })
    
        pvals.adj=do.call(rbind, sapply(out, function (res) res[[5]]))
        pvals.adj.fdr=p.adjust(pvals.adj[,1], method="fdr")
        pvals.adj.hol=p.adjust(pvals.adj[,1], method="holm")
        
        tab.cont=do.call(rbind, sapply(out, function (res) res[[1]]))
        p.1=formatDouble(pvals.adj.hol[startsWith(names(pvals.adj.hol), "cont.")], 3, remove.leading0=F); p.1=sub("0.000","<0.001",p.1)
        p.2=formatDouble(pvals.adj.fdr[startsWith(names(pvals.adj.fdr), "cont.")], 3, remove.leading0=F); p.2=sub("0.000","<0.001",p.2)
        tab.cont[,"p.1"]=p.1
        tab.cont[,"p.2"]=p.2
        
        tab.cont.scaled=do.call(rbind, sapply(out, function (res) res[[6]]))
        p.1=formatDouble(pvals.adj.hol[startsWith(names(pvals.adj.hol), "cont.")], 3, remove.leading0=F); p.1=sub("0.000","<0.001",p.1)
        p.2=formatDouble(pvals.adj.fdr[startsWith(names(pvals.adj.fdr), "cont.")], 3, remove.leading0=F); p.2=sub("0.000","<0.001",p.2)
        tab.cont.scaled[,"p.1"]=p.1
        tab.cont.scaled[,"p.2"]=p.2
        
        tab.cat=do.call(rbind, sapply(out, function (res) res[[2]]))
        overall.p.1=formatDouble(pvals.adj.hol[startsWith(names(pvals.adj.hol), "tri.")], 3, remove.leading0=F); p.1=sub("0.000","<0.001",p.1)
        overall.p.2=formatDouble(pvals.adj.fdr[startsWith(names(pvals.adj.fdr), "tri.")], 3, remove.leading0=F); p.2=sub("0.000","<0.001",p.2)
        tab.cat[,"overall.p.1"]=c(rbind(overall.p.1, NA,NA))
        tab.cat[,"overall.p.2"]=c(rbind(overall.p.2, NA,NA))
                
        mytex(tab.cont, file.name="CoR_univariable"%.%ifelse(ind==2,"NoADCP","")%.%"_svycoxph_pretty_ENSEMBLE_"%.%region%.%"_"%.%COR, align="c", include.colnames = F, save2input.only=T, input.foldername=here::here("output/meta"),
            col.headers=paste0("\\hline\n 
                 \\multicolumn{1}{l}{} & \\multicolumn{1}{c}{No. cases /}   & \\multicolumn{2}{c}{HR per 10-fold incr.}                     & \\multicolumn{1}{c}{P-value}   & \\multicolumn{1}{c}{q-value}   & \\multicolumn{1}{c}{FWER} \\\\ 
                 \\multicolumn{1}{l}{Immunologic Marker}            & \\multicolumn{1}{c}{No. at-risk**} & \\multicolumn{1}{c}{Pt. Est.} & \\multicolumn{1}{c}{95\\% CI} & \\multicolumn{1}{c}{(2-sided)} & \\multicolumn{1}{c}{***} & \\multicolumn{1}{c}{} \\\\ 
                 \\hline\n 
            ")
        )    
        
        mytex(tab.cont.scaled, file.name="CoR_univariable"%.%ifelse(ind==2,"NoADCP","")%.%"_svycoxph_pretty_scaled_ENSEMBLE_"%.%region%.%"_"%.%COR, align="c", include.colnames = F, save2input.only=T, input.foldername=here::here("output/meta"),
            col.headers=paste0("\\hline\n 
                 \\multicolumn{1}{l}{} & \\multicolumn{1}{c}{No. cases /}   & \\multicolumn{2}{c}{HR per SD incr.}                     & \\multicolumn{1}{c}{P-value}   & \\multicolumn{1}{c}{q-value}   & \\multicolumn{1}{c}{FWER} \\\\ 
                 \\multicolumn{1}{l}{Immunologic Marker}            & \\multicolumn{1}{c}{No. at-risk**} & \\multicolumn{1}{c}{Pt. Est.} & \\multicolumn{1}{c}{95\\% CI} & \\multicolumn{1}{c}{(2-sided)} & \\multicolumn{1}{c}{***} & \\multicolumn{1}{c}{} \\\\ 
                 \\hline\n 
            ")
        )    
        
        mytex(tab.cat, file.name="CoR_univariable"%.%ifelse(ind==2,"NoADCP","")%.%"_svycoxph_cat_pretty_ENSEMBLE_"%.%region%.%"_"%.%COR, align="c", include.colnames = F, save2input.only=T, input.foldername=here::here("output/meta"),
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
    
}





###################################################################################################
# cross tabulation (cases x non-cases) x (have bnAb data) x (have PsV ID50 data)

with(dat.ense.1, table(!is.na(Day29bindSpike) & !is.na(BbindSpike), !is.na(Day29pseudoneutid50) & !is.na(Bpseudoneutid50), EventIndPrimaryIncludeNotMolecConfirmedD29))

with(dat.ense.1, table(!is.na(Day29bindSpike) & !is.na(BbindSpike), !is.na(Day29pseudoneutid50) & !is.na(Bpseudoneutid50), EventIndPrimaryD29))

tmp=subset(dat.ense.1, !(!is.na(Day29bindSpike) & !is.na(BbindSpike)) & !(!is.na(Day29pseudoneutid50) & !is.na(Bpseudoneutid50)) & EventIndPrimary==1)

###################################################################################################
# correlations between markers

# ENSEMBLE
mypdf(mfrow=c(2,2), file="output/meta/ensemble_corplot")
with(subset(dat.ense.1, SubcohortInd==1), {
    corplot(Day29pseudoneutid50, Day29ADCP, xlab="ID50", ylab="ADCP", method="s")
    corplot(Day29pseudoneutid50, Day29bindRBD, xlab="ID50", ylab="bAb RBD", method="s")
    corplot(Day29bindSpike, Day29bindRBD, xlab="bAb Spike", ylab="bAb RBD", method="s")
    corplot(Day29bindRBD, Day29ADCP, xlab="bAb RBD", ylab="ADCP", method="s")
})
dev.off()

with(dat.ense.1, table(Day29pseudoneutid50cat, Day29ADCPcat))
with(dat.ense.1, table(Day29pseudoneutid50cat, Day29bindRBDcat))
with(dat.ense.1, table(Day29bindSpikecat, Day29bindRBDcat))
with(dat.ense.1, table(Day29bindRBDcat, Day29ADCPcat))

# COVE
mypdf(mfrow=c(2,2), file="output/meta/cove_corplot")
with(subset(dat.cove.1, SubcohortInd==1), {
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
get.marginalized.risk.no.marker(Surv(EventTimePrimary, EventIndPrimary) ~ risk_score + MinorityInd + HighRiskInd, dat.cove.0, day=70)
# 0.039606
get.marginalized.risk.no.marker(Surv(EventTimePrimary, EventIndPrimary) ~ risk_score + as.factor(Region), dat.ense.0, day=66)
# 0.0525176
get.marginalized.risk.no.marker(Surv(EventTimePrimary, EventIndPrimary) ~ risk_score, dat.pla.seroneg.id50.la, day=48)
# 0.0393271

get.marginalized.risk.no.marker(Surv(EventTimePrimary, EventIndPrimary) ~ risk_score + MinorityInd + HighRiskInd, dat.cove.1, day=89)
# 0.00328746
get.marginalized.risk.no.marker(Surv(EventTimePrimary, EventIndPrimary) ~ risk_score + as.factor(Region), dat.ense.1, day=66)
# 0.0196666
get.marginalized.risk.no.marker(Surv(EventTimePrimary, EventIndPrimary) ~ risk_score, dat.vac.seroneg.id50.la, day=48)
# 0.0165755


get.ve=function(form.0, dat.1, dat.0, marker.name, cuts=c(30,100), t) {
    dat.1$discrete.marker=factor(cut(dat.1[[marker.name]], breaks=c(-Inf,log10(cuts),Inf) ))
    print(table(dat.1$discrete.marker, dat.1$EventIndPrimary)             )
    f1=update(form.0, ~ . + discrete.marker)
    tmp.design=twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dat.1)
    fit.risk=try(svycoxph(f1, design=tmp.design))
    if(inherits(fit.risk,"try-error")) stop("syvcoxph fit failed")
    prob.vac=marginalized.risk(fit.risk, "discrete.marker", data=subset(dat.1, ph2=1), ss=NULL, weights=subset(dat.1, ph2=1)$wt, t=t, categorical.s=T, verbose=F)        
    #print(prob.vac)
    prev.plac = get.marginalized.risk.no.marker(form.0, dat.0, t)
    ve = c(1 - prob.vac/prev.plac)    
    ve
}

# COVE
tab.1=get.ve(Surv(EventTimePrimary, EventIndPrimary) ~ risk_score + MinorityInd + HighRiskInd, dat.cove.1, dat.cove.0, "Day57pseudoneutid50", cuts=c(30,100), t=100)              
#                0   1
#  (-Inf,1.48]  21   3
#  (1.48,2]    108   9
#  (2, Inf]    919  28

#(-Inf,1.48]    (1.48,2]    (2, Inf] 
#   0.843034    0.852950    0.955792 
        
# ENSEMBLE, pooled
tab.2=get.ve(Surv(EventTimePrimary, EventIndPrimary) ~ risk_score + as.factor(Region), dat.vac.seroneg.id50, dat.ense.0, "Day29pseudoneutid50", cuts=c(30,100), t=66) 
#                0   1
#  (-Inf,1.48] 590  76
#  (1.48,2]    131  11
#  (2, Inf]     57   2
#(-Inf,1.48]    (1.48,2]    (2, Inf] 
#   0.471382    0.707984    0.871407 
get.ve(Surv(EventTimePrimary, EventIndPrimary) ~ risk_score + as.factor(Region), dat.vac.seroneg.id50, dat.ense.0, "Day29pseudoneutid50", cuts=c(30,100), t=48) 
#(-Inf,1.48]    (1.48,2]    (2, Inf] 
#   0.580821    0.769225    0.898612
   
# ENSEMBLE NA. Cuts are chosen as 10, 30 because there are no data above 100
tab.3=get.ve(Surv(EventTimePrimary, EventIndPrimary) ~ risk_score, dat.vac.seroneg.id50.na, dat.pla.seroneg.id50.na, "Day29pseudoneutid50", cuts=c(30), t=66) 
#               0   1
#  (-Inf,1.48] 299  22
#  (1.48, Inf]  84   2
   
#ENSEMBLE LA
tab.4=get.ve(Surv(EventTimePrimary, EventIndPrimary) ~ risk_score, dat.vac.seroneg.id50.la, dat.pla.seroneg.id50.la, "Day29pseudoneutid50", cuts=c(30,100), t=48) 
#                0   1
#  (-Inf,1.48] 153  39
#  (1.48,2]     37   7
#  (2, Inf]     18   2
   
# ENSEMBLE SA. Cuts are chosen as 10, 30 because there are no data above 100
tab.5=get.ve(Surv(EventTimePrimary, EventIndPrimary) ~ risk_score, dat.vac.seroneg.id50.sa, dat.pla.seroneg.id50.sa, "Day29pseudoneutid50", cuts=c(30), t=26) 
#                0   1
#  (-Inf,1]    100  13
#  (1,1.48]     38   2
#  (1.48, Inf]  49   2

tab=rbind(COVE=tab.1, "ENSEMBLE pooled"=tab.2, "ENSEMBLE NA"=c(tab.3,NA), "ENSEMBLE LA"=tab.4, "ENSEMBLE SA"=c(tab.5,NA))
colnames(tab)=c("<=30",">30,<=100",">100")
tab.1
mytex(tab, file="id50_ve", align="c", save2input.only=T, input.foldername="output/meta")
