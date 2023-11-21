#Sys.setenv(TRIAL = "id27hpv"); COR="M18"; Sys.setenv(VERBOSE = 1) 

print(date())
renv::activate(project = here::here(".."))     
source(here::here("..", "_common.R")) # dat.mock is made

library(kyotil) # p.adj.perm, getFormattedSummary
library(marginalizedRisk)
library(survey)
library(parallel)
library(Hmisc) # wtd.quantile, cut2
library(xtable) # this is a dependency of kyotil
library (osDesign)

source(here::here("code", "params.R"))
time.start=Sys.time()
myprint(study_name)
myprint(verbose)

# redefine form.0
form.0 = update (EventIndPrimary~1, as.formula(config$covariates_riskscore))
print(form.0)


# path for figures and tables etc
save.results.to = here::here("output");                                if (!dir.exists(save.results.to))  dir.create(save.results.to)
save.results.to = paste0(save.results.to, "/", attr(config,"config")); if (!dir.exists(save.results.to))  dir.create(save.results.to)
save.results.to = paste0(save.results.to, "/", COR,"/");               if (!dir.exists(save.results.to))  dir.create(save.results.to)
print(paste0("save.results.to equals ", save.results.to))

# append to file names for figures and tables
# defined differently in cor_coxph_xx.R
fname.suffix = study_name


# B=1e3 and numPerm=1e4 take 10 min to run with 30 CPUS for one analysis
B <-       config$num_boot_replicates 
numPerm <- config$num_perm_replicates # number permutation replicates 1e4
myprint(B)
myprint(numPerm)

# define an alias for EventIndPrimaryDxx
dat.mock$yy=dat.mock[[config.cor$EventIndPrimary]]

# there is only one analysis population
dat.vac.seroneg=subset(dat.mock, ph1)


# define trichotomized markers
dat.vac.seroneg = add.trichotomized.markers (dat.vac.seroneg, all.markers, wt.col.name="wt")

# for hpv33 and above, too few positive responses, change to cut into to low and high
for (a in c("M18bindL1L2_HPV33", "M18bindL1L2_HPV45", "M18bindL1L2_HPV52", "M18bindL1L2_HPV58")) {        
    cutpoint=min(dat.vac.seroneg[[a]], na.rm = T)
    dat.vac.seroneg[[a%.%'cat']] = cut(dat.vac.seroneg[[a]], breaks = c(-Inf, cutpoint, Inf))
    attr(dat.vac.seroneg, "marker.cutpoints")[[a]] = cutpoint
}


marker.cutpoints=attr(dat.vac.seroneg, "marker.cutpoints")
for (a in all.markers) {        
    q.a=marker.cutpoints[[a]]
    if (startsWith(a, "Day")) {
        # not fold change
        write(paste0(labels.axis[1,marker.name.to.assay(a)], " [", concatList(round(q.a, 2), ", "), ")%"), file=paste0(save.results.to, "cutpoints_", a, "_"%.%study_name))
    } else {
        # fold change
        # gsub("_", "\\\_", a, fixed = TRUE) is a bandaid to escape the marker name for latex, which may have _
        write(paste0(gsub("_", "\\_", a, fixed = TRUE), " [", concatList(round(q.a, 2), ", "), ")%"), file=paste0(save.results.to, "cutpoints_", a, "_"%.%study_name))
    }
}

#create twophase design object
design.vacc.seroneg<-twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dat.vac.seroneg)
with(dat.vac.seroneg, table(Wstratum, ph2))
    
# create verification object to be populated by the following scripts
rv=list() 
rv$marker.cutpoints=marker.cutpoints

# getting some quantiles
#10**wtd.quantile(dat.vac.seroneg$Day57pseudoneutid50, dat.vac.seroneg$wt, c(0.025, 0.05, seq(.2,.9,by=0.01),seq(.9,.99,by=0.005)))

# table of ph1 and ph2 cases
tab=with(dat.vac.seroneg, table(ph2, EventIndPrimary))
names(dimnames(tab))[2]="Event Indicator"
print(tab)
mytex(tab, file.name="tab1", save2input.only=T, input.foldername=save.results.to)

# for use in competing risk estimation
dat.vac.seroneg.ph2=subset(dat.vac.seroneg, ph2)

begin=Sys.time()



###################################################################################################
# get OR on continuous markers
###################################################################################################
    
tab=with(dat.vac.seroneg, table(tps.stratum, EventIndPrimary))
nn0=tab[,1]
nn1=tab[,2]

fits=list()
fits.scaled=list()
for (i in 1:2) { # 1: not scaled, 2: scaled
    for (a in all.markers) {
        fit = tps(update (form.0,  as.formula('~. + '%.%ifelse(i==2,"scale("%.%a%.%")", a))), 
                  dat.vac.seroneg.ph2, nn0 = nn0, nn1 = nn1, 
                  group = dat.vac.seroneg.ph2$tps.stratum) 
            # method = "WL"# to be more robust
        
        ## getFormattedSummary() has this logic, so no necessary
        # res = summary(fit)$coef
        # robust = T # better for small samples
        # idx = ifelse(robust, "Emp ", "Mod ")
        # out = cbind(res[, c("Value", idx %.% "p")],
        #             `lower bound` = res[, "Value"] - 1.96 * res[, idx %.% "SE"],
        #             `upper bound` = res[, "Value"] + 1.96 * res[, idx %.% "SE"])
        # 
        # out[, c(1, 3, 4)] = exp(out[, c(1, 3, 4)])
        # colnames(out) = c("OR", "p.value", "(lower", "upper)")
        # fit$coefficients = out
        
        if (i==1) fits[[a]]=fit else fits.scaled[[a]]=fit
    }
}

assertthat::assert_that(all(abs(fits$M18bindL1L2_mdw$coef-c(-4.3268181,-0.0285102,-0.1003947))<1e-6), msg = "failed cor_logistic unit testing")    
    
natrisk=nrow(dat.vac.seroneg)
nevents=sum(dat.vac.seroneg$yy==1)

# make pretty table
rows=length(coef(fits[[1]]))
est=getFormattedSummary(fits, exp=T, robust=T, rows=rows, type=1)
ci= getFormattedSummary(fits, exp=T, robust=T, rows=rows, type=13)
p=  getFormattedSummary(fits, exp=T, robust=T, rows=rows, type=10)
est.scaled=getFormattedSummary(fits.scaled, exp=T, robust=T, rows=rows, type=1)
ci.scaled= getFormattedSummary(fits.scaled, exp=T, robust=T, rows=rows, type=13)

pvals.cont = sapply(fits, function(x) {
    tmp=getFixedEf(x)
    p.val.col=which(startsWith(tolower(colnames(tmp)),"p"))
    tmp[nrow(tmp),p.val.col]
})




###################################################################################################
if(verbose) print("regression for trichotomized markers")

marker.levels = sapply(all.markers, function(a) length(table(dat.vac.seroneg[[a%.%"cat"]])))

fits.tri=list()
for (a in all.markers) {
    if(verbose) myprint(a)
    f= update(form.0, as.formula(paste0("~.+", a, "cat")))
    
    fit = tps(f, dat.vac.seroneg.ph2, nn0 = nn0, nn1 = nn1, 
              group = dat.vac.seroneg.ph2$tps.stratum) 
    # method = "WL"# to be more robust
    
    # res = summary(fit)$coef
    # robust = FALSE # better for small samples
    # idx = ifelse(robust, "Emp ", "Mod ")
    # out = cbind(res[, c("Value", idx %.% "p")], 
    #             `lower bound` = res[, "Value"] - 1.96 * res[, idx %.% "SE"], 
    #             `upper bound` = res[, "Value"] + 1.96 * res[, idx %.% "SE"])
    # 
    # out[, c(1, 3, 4)] = exp(out[, c(1, 3, 4)])
    # colnames(out) = c("OR", "p.value", "(lower", "upper)")
    # fit$coefficients = out
    
    fits.tri[[a]]=fit
}

fits.tri.coef.ls= lapply(fits.tri, function (fit) getFixedEf(fit, robust=T))


# get generalized Wald p values
overall.p.tri=sapply(all.markers, function(a) {
    fit=fits.tri[[a]]
    rows=length(fit$coef) - (marker.levels[a]-2):0
    
    if (length(fit)==1) NA else {
        stat=coef(fit)[rows] %*% solve(vcov(fit,robust=T)[rows,rows]) %*% coef(fit)[rows]
        pchisq(stat, length(rows), lower.tail = FALSE)
    }
})
#
overall.p.0=formatDouble(c(rbind(overall.p.tri, NA,NA)), digits=3, remove.leading0 = F);   overall.p.0=sub("0.000","<0.001",overall.p.0)



###################################################################################################
if(verbose) print("# multitesting adjustment for continuous and trichotomized markers together")

# If primary_assays is not defined in config, multitesting adjustment is over all assays. 
# If primary_assays defined, multitesting adjustment is only over this subset. If this set is empty, then no multitesting adjustment is done

p.unadj=c(cont=pvals.cont, tri=overall.p.tri)
# save a copy for later use
p.unadj.1 = p.unadj 
# pick out a subset based on config
if (!is.null(config$primary_assays)) {
    if (length(config$primary_assays)>0) {
        p.unadj = p.unadj[c(paste0("cont.",DayPrefix,tpeak,primary_assays), paste0("tri.",DayPrefix,tpeak,primary_assays))]
    } else {
        p.unadj=c()
    }
}

#### Holm and FDR adjustment
pvals.adj.fdr=p.adjust(p.unadj, method="fdr")
pvals.adj.hol=p.adjust(p.unadj, method="holm")

print("not doing Westfall and Young")
pvals.adj=cbind(p.unadj, p.FWER=pvals.adj.hol, p.FDR=pvals.adj.fdr)
write(NA, file=paste0(save.results.to, "permutation_replicates_"%.%study_name))     # so the rmd file can compile


# since we take ID80 out earlier, we may need to add it back for the table and we do it with the help of p.unadj.1
pvals.adj = cbind(p.unadj=p.unadj.1, pvals.adj[match(names(p.unadj.1), rownames(pvals.adj)),2:3, drop=F])



###################################################################################################
# make continuous markers table

p.1=formatDouble(pvals.adj["cont."%.%names(pvals.cont),"p.FWER"], 3, remove.leading0=F); p.1=sub("0.000","<0.001",p.1)
p.2=formatDouble(pvals.adj["cont."%.%names(pvals.cont),"p.FDR" ], 3, remove.leading0=F); p.2=sub("0.000","<0.001",p.2)


tab.1=cbind(paste0(nevents, "/", format(natrisk, big.mark=",")), t(est), t(ci), t(p), p.2, p.1)
rownames(tab.1)=all.markers.names.short
tab.1
mytex(tab.1, file.name="CoR_univariable_logistic_pretty_"%.%study_name, align="c", include.colnames = F, save2input.only=T, input.foldername=save.results.to,
      col.headers=paste0("\\hline\n 
         \\multicolumn{1}{l}{", study_name, "} & \\multicolumn{1}{c}{No. cases /}   & \\multicolumn{2}{c}{OR per 10-fold incr.}                     & \\multicolumn{1}{c}{P-value}   & \\multicolumn{1}{c}{q-value}   & \\multicolumn{1}{c}{FWER} \\\\ 
         \\multicolumn{1}{l}{Immunologic Marker}            & \\multicolumn{1}{c}{No. at-risk**} & \\multicolumn{1}{c}{Pt. Est.} & \\multicolumn{1}{c}{95\\% CI} & \\multicolumn{1}{c}{(2-sided)} & \\multicolumn{1}{c}{***} & \\multicolumn{1}{c}{} \\\\ 
         \\hline\n 
    "),
      longtable=T, 
      label=paste0("tab:CoR_univariable_logistic_pretty"), 
      caption.placement = "top", 
      caption=paste0("Inference for M ", tpeak, "antibody marker covariate-adjusted correlates of risk of ", config.cor$txt.endpoint, " pooled over treatment arms: Odds ratios per 10-fold increment in the marker*")
)
tab.cont=tab.1

tab.1.nop12=cbind(paste0(nevents, "/", format(natrisk, big.mark=",")), t(est), t(ci), t(p))
rownames(tab.1.nop12)=all.markers.names.short
rv$tab.1=tab.1.nop12

# scaled markers
tab.1.scaled=cbind(paste0(nevents, "/", format(natrisk, big.mark=",")), t(est.scaled), t(ci.scaled), t(p), p.2, p.1)
rownames(tab.1.scaled)=all.markers.names.short
tab.1.scaled
mytex(tab.1.scaled, file.name="CoR_univariable_logistic_pretty_scaled_"%.%study_name, align="c", include.colnames = F, save2input.only=T, input.foldername=save.results.to,
      col.headers=paste0("\\hline\n 
         \\multicolumn{1}{l}{", (study_name), "} & \\multicolumn{1}{c}{No. cases /}   & \\multicolumn{2}{c}{OR per SD incr.}                     & \\multicolumn{1}{c}{P-value}   & \\multicolumn{1}{c}{q-value}   & \\multicolumn{1}{c}{FWER} \\\\ 
         \\multicolumn{1}{l}{Immunologic Marker}            & \\multicolumn{1}{c}{No. at-risk**} & \\multicolumn{1}{c}{Pt. Est.} & \\multicolumn{1}{c}{95\\% CI} & \\multicolumn{1}{c}{(2-sided)} & \\multicolumn{1}{c}{***} & \\multicolumn{1}{c}{} \\\\ 
         \\hline\n 
    "),
      longtable=T, 
      label=paste0("tab:CoR_univariable_logistic_pretty_scaled"), 
      caption.placement = "top", 
      caption=paste0("Inference for Day ", tpeak, "antibody marker covariate-adjusted correlates of risk of ", config.cor$txt.endpoint, " pooled over treatment groups: Odds ratios per SD increment in the marker*")
)
tab.cont.scaled=tab.1.scaled



###################################################################################################
# make trichotomized markers table


get.est=function(a) {
  fit=fits.tri[[a]]
  out = getFormattedSummary(list(fit), exp=T, robust=T, rows=rows, type=1)
  if (length(out)==1) c(NA,out) else out
}
get.ci =function(a) {
  fit=fits.tri[[a]]
  out = getFormattedSummary(list(fit), exp=T, robust=T, rows=rows, type=13)
  if (length(out)==1) c(NA,out) else out
}
get.p  =function(a) {
  fit=fits.tri[[a]]
  out = getFormattedSummary(list(fit), exp=T, robust=T, rows=rows, type=10)
  if (length(out)==1) c(NA,out) else out
}
get.est(all.markers[1]); get.ci (all.markers[1]); get.p  (all.markers[1])
get.est(all.markers[6]); get.ci (all.markers[6]); get.p  (all.markers[6])

# regression parameters
est=c(rbind(1.00,  sapply(all.markers, function (a) get.est(a))))
ci= c(rbind("N/A", sapply(all.markers, function (a) get.ci (a))))
p=  c(rbind("N/A", sapply(all.markers, function (a) get.p  (a))))

overall.p.1=formatDouble(pvals.adj["tri."%.%names(pvals.cont),"p.FWER"], 3, remove.leading0=F);   overall.p.1=sub("0.000","<0.001",overall.p.1)
overall.p.2=formatDouble(pvals.adj["tri."%.%names(pvals.cont),"p.FDR" ], 3, remove.leading0=F);   overall.p.2=sub("0.000","<0.001",overall.p.2)

# add space
overall.p.1=c(rbind(overall.p.1, NA,NA))
overall.p.2=c(rbind(overall.p.2, NA,NA))


# n cases and n at risk
natrisk = round(c(sapply (all.markers%.%"cat", function(a) {
  out = aggregate(subset(dat.vac.seroneg,ph2==1)        [["wt"]], subset(dat.vac.seroneg,ph2==1        )[a], sum, na.rm=T, drop=F)[,2]
  if (length(out)==3) out else c(out,NA)
} )))
nevents = round(c(sapply (all.markers%.%"cat", function(a) {
  out = aggregate(subset(dat.vac.seroneg,yy==1 & ph2==1)[["wt"]], subset(dat.vac.seroneg,yy==1 & ph2==1)[a], sum, na.rm=T, drop=F)[,2]
  if (length(out)==3) out else c(out,NA)
} )))

tab=cbind(
    rep(c("Lower","Middle","Upper"), length(p)/3),
    paste0(nevents, "/", format(natrisk, big.mark=",",digit=0, scientific=F)), 
    formatDouble(nevents/natrisk, digit=4, remove.leading0=F),
    est, ci, p, overall.p.0, overall.p.2, overall.p.1
)

tmp=rbind(all.markers.names.short, "", "")
rownames(tab)=c(tmp)
tab

# use longtable because this table could be long, e.g. in hvtn705second
mytex(tab[1:(nrow(tab)),], file.name="CoR_univariable_logistic_cat_pretty_"%.%study_name, align="c", include.colnames = F, save2input.only=T, input.foldername=save.results.to,
      col.headers=paste0("\\hline\n 
         \\multicolumn{1}{l}{", (study_name), "} & \\multicolumn{1}{c}{Tertile}   & \\multicolumn{1}{c}{No. cases /}   & \\multicolumn{1}{c}{Attack}   & \\multicolumn{2}{c}{Haz. Ratio}                     & \\multicolumn{1}{c}{P-value}   & \\multicolumn{1}{c}{Overall P-}      & \\multicolumn{1}{c}{Overall q-}   & \\multicolumn{1}{c}{Overall} \\\\ 
         \\multicolumn{1}{l}{Immunologic Marker}            & \\multicolumn{1}{c}{}          & \\multicolumn{1}{c}{No. at-risk**} & \\multicolumn{1}{c}{rate}   & \\multicolumn{1}{c}{Pt. Est.} & \\multicolumn{1}{c}{95\\% CI} & \\multicolumn{1}{c}{(2-sided)} & \\multicolumn{1}{c}{value***} & \\multicolumn{1}{c}{value $\\dagger$} & \\multicolumn{1}{c}{FWER} \\\\ 
         \\hline\n 
      "),        
      longtable=T, 
      label=paste0("tab:CoR_univariable_logistic_cat_pretty_", study_name), 
      caption.placement = "top", 
      caption=paste0("Inference for Day ", tpeak, "antibody marker covariate-adjusted correlates of risk of ", config.cor$txt.endpoint, " in the vaccine group: Hazard ratios for Middle vs. Upper tertile vs. Lower tertile*")
)


tab.nop12=cbind(
    rep(c("Lower","Middle","Upper"), length(p)/3), 
    paste0(nevents, "/", format(natrisk, big.mark=",",digit=0, scientific=F)), 
    formatDouble(nevents/natrisk, digit=4, remove.leading0=F),
    est, ci, p, overall.p.0
)
rownames(tab.nop12)=c(rbind(all.markers.names.short, "", ""))
rv$tab.2=tab.nop12




print(date())
print("cor_coxph run time: "%.%format(Sys.time()-time.start, digits=1))
