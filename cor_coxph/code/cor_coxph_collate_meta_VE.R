# 2022.3.3 collate is not as necessary as before because we combined bAb and PsV into one data file

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
library(WeightedROC)
Sys.setenv("TRIAL"="janssen_pooled_EUA") # value does not matter since we just need to load the common functions in _common.R
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
dat.vac.seroneg.bAbID50=load.data("janssen_pooled_EUA", "D29IncludeNotMolecConfirmedstart1")
# combine markers into one data frame
dat.ense.1=dat.vac.seroneg.bAbID50

dat.ense.0=load.data("janssen_pooled_EUA", "D29IncludeNotMolecConfirmedstart1", trt=0)
    
dat.cove.1=load.data("moderna_real", "D57")
dat.cove.0=load.data("moderna_real", "D57", trt=0)

dat.vac.seroneg.id50.na=load.data("janssen_na_EUA", "D29IncludeNotMolecConfirmedstart1")
dat.pla.seroneg.id50.na=load.data("janssen_na_EUA", "D29IncludeNotMolecConfirmedstart1", trt=0)
dat.vac.seroneg.id50.la=load.data("janssen_la_EUA", "D29IncludeNotMolecConfirmedstart1")
dat.pla.seroneg.id50.la=load.data("janssen_la_EUA", "D29IncludeNotMolecConfirmedstart1", trt=0)
dat.vac.seroneg.id50.sa=load.data("janssen_sa_EUA", "D29IncludeNotMolecConfirmedstart1")
dat.pla.seroneg.id50.sa=load.data("janssen_sa_EUA", "D29IncludeNotMolecConfirmedstart1", trt=0)




###################################################################################################
# make tables of HR that contain all markers

#for (region in c("pooled","na","la","sa")) {
for (region in c("pooled")) {
        
    if(verbose) myprint(region)
    
    trials=c("janssen_"%.%region%.%"_EUA") 
        
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
                
        mytex(tab.cont, file.name="CoR_univariable_svycoxph_pretty_ENSEMBLE_"%.%region%.%"_"%.%COR, align="c", include.colnames = F, save2input.only=T, input.foldername=here::here("output/meta"),
            col.headers=paste0("\\hline\n 
                 \\multicolumn{1}{l}{} & \\multicolumn{1}{c}{No. cases /}   & \\multicolumn{2}{c}{HR per 10-fold incr.}                     & \\multicolumn{1}{c}{P-value}   & \\multicolumn{1}{c}{q-value}   & \\multicolumn{1}{c}{FWER} \\\\ 
                 \\multicolumn{1}{l}{Immunologic Marker}            & \\multicolumn{1}{c}{No. at-risk**} & \\multicolumn{1}{c}{Pt. Est.} & \\multicolumn{1}{c}{95\\% CI} & \\multicolumn{1}{c}{(2-sided)} & \\multicolumn{1}{c}{***} & \\multicolumn{1}{c}{} \\\\ 
                 \\hline\n 
            ")
        )    
        
        mytex(tab.cont.scaled, file.name="CoR_univariable_svycoxph_pretty_scaled_ENSEMBLE_"%.%region%.%"_"%.%COR, align="c", include.colnames = F, save2input.only=T, input.foldername=here::here("output/meta"),
            col.headers=paste0("\\hline\n 
                 \\multicolumn{1}{l}{} & \\multicolumn{1}{c}{No. cases /}   & \\multicolumn{2}{c}{HR per SD incr.}                     & \\multicolumn{1}{c}{P-value}   & \\multicolumn{1}{c}{q-value}   & \\multicolumn{1}{c}{FWER} \\\\ 
                 \\multicolumn{1}{l}{Immunologic Marker}            & \\multicolumn{1}{c}{No. at-risk**} & \\multicolumn{1}{c}{Pt. Est.} & \\multicolumn{1}{c}{95\\% CI} & \\multicolumn{1}{c}{(2-sided)} & \\multicolumn{1}{c}{***} & \\multicolumn{1}{c}{} \\\\ 
                 \\hline\n 
            ")
        )    
        
        mytex(tab.cat, file.name="CoR_univariable_svycoxph_cat_pretty_ENSEMBLE_"%.%region%.%"_"%.%COR, align="c", include.colnames = F, save2input.only=T, input.foldername=here::here("output/meta"),
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
# cross tabulation (cases x non-cases) x (have bnAb data) x (have PsV ID50 data)

with(dat.ense.1, table(!is.na(Day29bindSpike) & !is.na(BbindSpike), !is.na(Day29pseudoneutid50) & !is.na(Bpseudoneutid50), EventIndPrimaryIncludeNotMolecConfirmedD29))

with(dat.ense.1, table(!is.na(Day29bindSpike) & !is.na(BbindSpike), !is.na(Day29pseudoneutid50) & !is.na(Bpseudoneutid50), EventIndPrimaryD29))

tmp=subset(dat.ense.1, !(!is.na(Day29bindSpike) & !is.na(BbindSpike)) & !(!is.na(Day29pseudoneutid50) & !is.na(Bpseudoneutid50)) & EventIndPrimary==1)


###################################################################################################
# correlations between markers

# ENSEMBLE
mypdf(mfrow=c(2,3), file="output/meta/ensemble_corplot")
with(subset(dat.ense.1, SubcohortInd==1), {
    corplot(Day29pseudoneutid50, Day29ADCP, xlab="ID50", ylab="ADCP", method="s", add.diagonal.line=T)
    corplot(Day29pseudoneutid50, Day29bindRBD, xlab="ID50", ylab="bAb RBD", method="s", add.diagonal.line=T)
    corplot(Day29pseudoneutid50, Day29bindSpike, xlab="ID50", ylab="bAb Spike", method="s", add.diagonal.line=T)
    corplot(Day29bindSpike, Day29bindRBD, xlab="bAb Spike", ylab="bAb RBD", method="s", add.diagonal.line=T)
    corplot(Day29bindSpike, Day29ADCP, xlab="bAb Spike", ylab="ADCP", method="s", add.diagonal.line=T)
    corplot(Day29bindRBD, Day29ADCP, xlab="bAb RBD", ylab="ADCP", method="s", add.diagonal.line=T)
})
dev.off()

with(dat.ense.1, table(Day29pseudoneutid50cat, Day29ADCPcat))
with(dat.ense.1, table(Day29pseudoneutid50cat, Day29bindRBDcat))
with(dat.ense.1, table(Day29bindSpikecat, Day29bindRBDcat))
with(dat.ense.1, table(Day29bindRBDcat, Day29ADCPcat))

with(subset(dat.ense.1, ph2.D29start1==1), WeightedAUC(WeightedROC(-Day29bindSpike, EventIndPrimary, weight = wt.D29start1)))
#0.55435
with(subset(dat.ense.1, ph2.D29start1==1), WeightedAUC(WeightedROC(-Day29bindRBD, EventIndPrimary, weight = wt.D29start1)))
#0.560554
with(subset(dat.ense.1, ph2.D29start1ADCP==1), WeightedAUC(WeightedROC(-Day29ADCP, EventIndPrimary, weight = wt.D29start1ADCP)))
#0.590071

# COVE
mypdf(mfrow=c(2,2), file="output/meta/cove_corplot")
with(subset(dat.cove.1, SubcohortInd==1), {
    corplot(Day57pseudoneutid50, Day57bindSpike, xlab="ID50", ylab="bAb Spike", add.diagonal.line=F)
    corplot(Day57pseudoneutid50, Day57bindRBD, xlab="ID50", ylab="bAb RBD", add.diagonal.line=F)
    corplot(Day57bindSpike, Day57bindRBD, xlab="bAb Spike", ylab="bAb RBD", add.diagonal.line=T)
    corplot(Day57pseudoneutid50, Day57pseudoneutid80, xlab="ID50", ylab="ID80", add.diagonal.line=T)
})
dev.off()

with(dat.cove.1, table(Day57pseudoneutid50cat, Day57bindSpikecat))
with(dat.cove.1, table(Day57pseudoneutid50cat, Day57bindRBDcat))
with(dat.cove.1, table(Day57bindSpikecat, Day57bindRBDcat))
with(dat.cove.1, table(Day57pseudoneutid50cat, Day57pseudoneutid80cat))

with(subset(dat.cove.1, ph2.D57==1), WeightedAUC(WeightedROC(-Day57bindSpike, EventIndPrimary, weight = wt.D57)))
#0.658508
with(subset(dat.cove.1, ph2.D57==1), WeightedAUC(WeightedROC(-Day57bindRBD, EventIndPrimary, weight = wt.D57)))
#0.647224
with(subset(dat.cove.1, ph2.D57==1), WeightedAUC(WeightedROC(-Day57pseudoneutid50, EventIndPrimary, weight = wt.D57)))
#0.628095
with(subset(dat.cove.1, ph2.D57==1), WeightedAUC(WeightedROC(-Day57pseudoneutid80, EventIndPrimary, weight = wt.D57)))
#0.618422






####################################################################
# compare controlled VE for ID50 <30, 100> between cove and ensemble

10**wtd.quantile(dat.cove.1[["Day57pseudoneutid50"]], dat.cove.1$wt, probs = c(0.025,0.05,0.1,0.95,0.975))
#    2.5%      5.0%     10.0%     95.0%     97.5% 
#  31.5829   60.4836   88.0196 1090.9402 1308.6377 

10**wtd.quantile(dat.vac.seroneg.id50.na[["Day29pseudoneutid50"]], dat.vac.seroneg.id50.na$wt, c(0.025, 0.05,0.90,0.95,.975))
#   2.5%     5.0%    90.0%    95.0%    97.5% 
#  2.9988   2.9988  51.2652  88.1076 213.3432 



lcut=31.5829; ucut=213.3432

get.ve=function(form.0, dat.1, dat.0, marker.name, cuts=c(lcut,ucut), t) {
    dat.1$discrete.marker=factor(cut(dat.1[[marker.name]], breaks=c(-Inf,log10(cuts),Inf) ))
    print(table(dat.1$discrete.marker, dat.1$EventIndPrimary, useNA="ifany")             )
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
tab.1=get.ve(Surv(EventTimePrimary, EventIndPrimary) ~ risk_score + MinorityInd + HighRiskInd, dat.cove.1, dat.cove.0, "Day57pseudoneutid50", t=100); tab.1
#                0   1
#  (-Inf,1.51]  22   3
#  (1.51,2.01] 111   9
#  (2.01, Inf] 915  28
#(-Inf,1.51] (1.51,2.01] (2.01, Inf] 
#   0.836539    0.827545    0.946230 
   

# ENSEMBLE, pooled
tab.2=get.ve(Surv(EventTimePrimary, EventIndPrimary) ~ risk_score + as.factor(Region), dat.vac.seroneg.bAbID50, dat.ense.0, "Day29pseudoneutid50", t=66); tab.2
#                0   1
#  (-Inf,1.51] 601  76
#  (1.51,2.01] 120   9
#  (2.01, Inf]  54   3
#(-Inf,1.51] (1.51,2.01] (2.01, Inf] 
#   0.473873    0.727231    0.793189 
        
   
# ENSEMBLE NA. Only one cut b/c no data above 102
tab.3=get.ve(Surv(EventTimePrimary, EventIndPrimary) ~ risk_score, dat.vac.seroneg.id50.na, dat.pla.seroneg.id50.na, "Day29pseudoneutid50", cuts=31.5829, t=66); tab.3
#                0   1
#  (-Inf,1.51] 304  22
#  (1.51, Inf]  76   2
#(-Inf,1.51] (1.51, Inf] 
#   0.764128    0.815223 
   
#ENSEMBLE LA
tab.4=get.ve(Surv(EventTimePrimary, EventIndPrimary) ~ risk_score, dat.vac.seroneg.id50.la, dat.pla.seroneg.id50.la, "Day29pseudoneutid50", t=48); tab.4
#                0   1
#  (-Inf,1.51] 158  39
#  (1.51,2.01]  32   5
#  (2.01, Inf]  18   3
#(-Inf,1.51] (1.51,2.01] (2.01, Inf] 
#   0.515871    0.766733    0.708168 
      
# ENSEMBLE SA. Cuts are chosen as 10, 30 because there are no data above 100
tab.5=get.ve(Surv(EventTimePrimary, EventIndPrimary) ~ risk_score, dat.vac.seroneg.id50.sa, dat.pla.seroneg.id50.sa, "Day29pseudoneutid50", cuts=c(32), t=26); tab.5
#                0   1
#  (-Inf,1.51] 139  15
#  (1.51, Inf]  48   2
#(-Inf,1.51] (1.51, Inf] 
#   0.350872    0.768693


tab=rbind(COVE=tab.1, "ENSEMBLE pooled"=tab.2, "ENSEMBLE NA"=c(tab.3,NA), "ENSEMBLE LA"=tab.4, "ENSEMBLE SA"=c(tab.5,NA))
colnames(tab)=c("<=32",">32,<=102",">102")
tab
mytex(tab, file="id50_ve", align="c", save2input.only=T, input.foldername="output/meta")


# missingness

with(dat.ense.1, table(!is.na(Day29bindSpike), !is.na(Day29bindRBD), EventIndPrimaryIncludeNotMolecConfirmedD29)) # same missingness
with(dat.ense.1, table(!is.na(Day29pseudoneutid50), !is.na(Day29bindRBD), EventIndPrimaryIncludeNotMolecConfirmedD29))
with(dat.ense.1, table(!is.na(Day29pseudoneutid50), !is.na(Day29bindRBD), Region))
with(dat.ense.1, table(!is.na(Day29pseudoneutid50), !is.na(Day29bindRBD)))

with(dat.ense.1, table(!is.na(BbindSpike), !is.na(BbindRBD))) # same missingness
with(dat.ense.1, table(!is.na(Bpseudoneutid50), !is.na(BbindRBD)))
with(dat.ense.1, table(!is.na(Bpseudoneutid50), !is.na(Day29pseudoneutid50), EventIndPrimaryIncludeNotMolecConfirmedD29))
with(dat.ense.1, table(!is.na(BbindSpike), !is.na(Day29bindSpike), EventIndPrimaryIncludeNotMolecConfirmedD29)) 

with(dat.cove.1, table(!is.na(Day57pseudoneutid50), !is.na(Day57bindRBD), EventIndPrimaryD57))
