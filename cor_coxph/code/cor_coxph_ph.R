# mandatory input: 
  # fname.suffix, which is used in the file names to save results

# optional inputs
{
if (is.null(show.q)) show.q=T # control whether fwer and q values are shown in tables

# controls whether we are using survey package to handle two-phase samples or coxph for cohort
# for svy, we expect design.dat
# for coxph, we expect dat
if (is.null(use.svy)) use.svy=T 

# for T, we expect dat.pla.seroneg
if (is.null(has.plac)) has.plac=T # control whether there are placebo data
}

# makes one table for continuous markers and one table for discrete markers

myprint(fname.suffix)


###################################################################################################
if(verbose) print("Regression for continuous markers")

fits=list()
for (a in all.markers) {
    f= update(form.0, as.formula(paste0("~.+", a)))
    if (use.svy) {
      fits[[a]]=svycoxph(f, design=design.dat) 
    } else {
      fits[[a]]=coxph(f, dat) 
    }
}

# scaled marker
fits.scaled=list()
for (a in all.markers) {
    f= update(form.0, as.formula(paste0("~.+scale(", a, ")")))
    if (use.svy) {
      fits.scaled[[a]]=svycoxph(f, design=design.dat) 
    } else {
      fits.scaled[[a]]=coxph(f, dat) 
    }
}

# put coxph model coef together to save
fits.cont.coef.ls = lapply(fits, function (fit) getFixedEf(fit, robust=use.svy))

natrisk=nrow(dat)
nevents=sum(dat$yy==1)

# make pretty table
rows=length(coef(fits[[1]]))
est=getFormattedSummary(fits, exp=T, robust=use.svy, rows=rows, type=1)
ci= getFormattedSummary(fits, exp=T, robust=use.svy, rows=rows, type=7)
p=  getFormattedSummary(fits, exp=T, robust=use.svy, rows=rows, type=10)
est.scaled=getFormattedSummary(fits.scaled, exp=T, robust=use.svy, rows=rows, type=1)
ci.scaled= getFormattedSummary(fits.scaled, exp=T, robust=use.svy, rows=rows, type=7)

pvals.cont = sapply(fits, function(x) {
    tmp=getFixedEf(x)
    p.val.col=which(startsWith(tolower(colnames(tmp)),"p"))
    tmp[nrow(tmp),p.val.col]
})


###################################################################################################
if(verbose) print("regression for trichotomized markers")

fits.tri=list()
for (a in all.markers) {
    if(verbose>=2) myprint(a)
    f= update(form.0, as.formula(paste0("~.+", a, "cat")))
    if (use.svy) {
      fits.tri[[a]]=svycoxph(f, design=design.dat) 
    } else {
      fits.tri[[a]]=coxph(f, dat) 
    }
}
fits.tri=fits.tri

fits.tri.coef.ls= lapply(fits.tri, function (fit) getFixedEf(fit, robust=use.svy))


rows=length(coef(fits.tri[[1]]))-1:0
# get generalized Wald p values
overall.p.tri=sapply(fits.tri, function(fit) {
    if (length(fit)==1) NA else {
        stat=coef(fit)[rows] %*% solve(vcov(fit,robust=use.svy)[rows,rows]) %*% coef(fit)[rows]
        pchisq(stat, length(rows), lower.tail = FALSE)
    }
})
#
overall.p.0=formatDouble(c(rbind(overall.p.tri, NA,NA)), digits=3, remove.leading0 = F);   overall.p.0=sub("0.000","<0.001",overall.p.0)




###################################################################################################
if(verbose) print("multitesting adjustment for continuous and trichotomized markers together")

# If primary_assays is not defined in config, multitesting adjustment is over all assays. 
# If primary_assays defined, multitesting adjustment is only over this subset. If this set is empty, then no multitesting adjustment is done

p.unadj=c(cont=pvals.cont, tri=overall.p.tri)
# save a copy for later use
p.unadj.1 = p.unadj 
# pick out a subset based on config
if (!is.null(config$primary_assays)) {
    if (length(config$primary_assays)>0) {
        p.unadj = p.unadj[c("cont.Day"%.%tpeak%.%primary_assays, "tri.Day"%.%tpeak%.%primary_assays)]
    } else {
        p.unadj=c()
    }
}
if (TRIAL=="prevent19" | TRIAL=='prevent19_stage2') {
    # bindSpike tertiary has no cases in the upper tertile, cannot do P value
    p.unadj = p.unadj[startsWith(names(p.unadj), "cont."), drop=F]
}

#### Holm and FDR adjustment
pvals.adj.fdr=p.adjust(p.unadj, method="fdr")
pvals.adj.hol=p.adjust(p.unadj, method="holm")

if (length(p.unadj)>1) {
    print("doing Westfall and Young")
  
    #### Westfall and Young permutation-based adjustment
    if(!file.exists(paste0(save.results.to, "pvals.perm.",fname.suffix,".Rdata"))) {
        
        dat.ph2 = design.dat$phase1$sample$variables
        design.dat.perm=design.dat
        #design.dat.perm$phase1$full$variables
    
    #    # if want to only do multitesting when liveneutmn50 is included
    #    if (!"liveneutmn50" %in% assays) numPerm=5
        
        # TODO: there is no need to permute all.markers
        
        out=mclapply(1:numPerm, mc.cores = numCores, FUN=function(seed) {   
            # store the current rng state 
            save.seed <- try(get(".Random.seed", .GlobalEnv), silent=TRUE) 
            if (class(save.seed)=="try-error") {set.seed(1); save.seed <- get(".Random.seed", .GlobalEnv) }          
            set.seed(seed)        
            
            # permute markers in design.dat.perm
            new.idx=sample(1:nrow(dat.ph2))
            tmp=dat.ph2
            for (a in all.markers) {
                tmp[[a]]=tmp[[a]][new.idx]
                tmp[[a%.%"cat"]]=tmp[[a%.%"cat"]][new.idx]
            }
            design.dat.perm$phase1$sample$variables = tmp
            
            # rename all.markers so that out has the proper names. this is only effective within permutation
            names(all.markers)=all.markers
            out=c(
                cont=sapply (all.markers, function(a) {
                    f= update(form.0, as.formula(paste0("~.+", a)))
                    if (use.svy) {
                      fit=run.svycoxph(f, design=design.dat.perm) 
                    } else {
                      # TODO
                    }
                    if (length(fit)==1) NA else last(c(getFixedEf(fit)))
                })        
                ,    
                tri=sapply (all.markers, function(a) {
                    f= update(form.0, as.formula(paste0("~.+", a, "cat")))
                    if (use.svy) {
                      fit=run.svycoxph(f, design=design.dat.perm) 
                    } else {
                      # TODO
                    }
                    if (length(fit)==1) NA else last(c(getFixedEf(fit)))
                })
            )
            
            # restore rng state 
            assign(".Random.seed", save.seed, .GlobalEnv)    
            
            out
        })
        pvals.perm=do.call(rbind, out)
        save(pvals.perm, file=paste0(save.results.to, "pvals.perm."%.%fname.suffix%.%".Rdata"))
        
    } else {
        load(file=paste0(save.results.to, "pvals.perm."%.%fname.suffix%.%".Rdata"))
    }
    # save number of permutation replicates
    write(nrow(pvals.perm), file=paste0(save.results.to, "permutation_replicates_"%.%fname.suffix))
    
    
    if(any(is.na(p.unadj))) {
        pvals.adj = cbind(p.unadj=p.unadj, p.FWER=NA, p.FDR=NA)
    } else {
        pvals.adj = p.adj.perm (p.unadj, pvals.perm[,names(p.unadj)], alpha=1)  
    }
    if(verbose) print(pvals.adj)
    
} else {
    print("not doing Westfall and Young")
    pvals.adj=cbind(p.unadj, p.FWER=pvals.adj.hol, p.FDR=pvals.adj.fdr)
    write(NA, file=paste0(save.results.to, "permutation_replicates_"%.%fname.suffix))     # so the rmd file can compile
}



# since we take ID80 out earlier, we may need to add it back for the table and we do it with the help of p.unadj.1
pvals.adj = cbind(p.unadj=p.unadj.1, pvals.adj[match(names(p.unadj.1), rownames(pvals.adj)),2:3, drop=F])

if (TRIAL=="prevent19" | TRIAL=='prevent19_stage2') {
  # bindSpike tertiary has no cases in the upper tertile, cannot do P value
    # this code somehow works even though there are also RBD and ID50
    pvals.adj=rbind(pvals.adj, tri.Day35bindSpike=c(NA,NA,NA))
}



###################################################################################################
# make continuous markers table

p.1=formatDouble(pvals.adj["cont."%.%names(pvals.cont),"p.FWER"], 3, remove.leading0=F); p.1=sub("0.000","<0.001",p.1)
p.2=formatDouble(pvals.adj["cont."%.%names(pvals.cont),"p.FDR" ], 3, remove.leading0=F); p.2=sub("0.000","<0.001",p.2)


tab.1=cbind(paste0(nevents, "/", format(natrisk, big.mark=",")), t(est), t(ci), t(p), if(show.q) p.2, if(show.q) p.1)
rownames(tab.1)=all.markers.names.short
tab.1

if (show.q) {
  header=paste0("\\hline\n 
         \\multicolumn{1}{l}{", '', "} & \\multicolumn{1}{c}{No. cases /}   & \\multicolumn{2}{c}{HR per 10-fold incr.}                     & \\multicolumn{1}{c}{P-value}   & \\multicolumn{1}{c}{q-value}   & \\multicolumn{1}{c}{FWER} \\\\ 
         \\multicolumn{1}{l}{Immunologic Marker}            & \\multicolumn{1}{c}{No. at-risk**} & \\multicolumn{1}{c}{Pt. Est.} & \\multicolumn{1}{c}{95\\% CI} & \\multicolumn{1}{c}{} & \\multicolumn{1}{c}{***} & \\multicolumn{1}{c}{} \\\\ 
         \\hline\n 
    ")
} else {
  header=paste0("\\hline\n 
         \\multicolumn{1}{l}{", '', "} & \\multicolumn{1}{c}{No. cases /}   & \\multicolumn{2}{c}{HR per 10-fold incr.}                     & \\multicolumn{1}{c}{P-value}    \\\\ 
         \\multicolumn{1}{l}{Immunologic Marker}            & \\multicolumn{1}{c}{No. at-risk**} & \\multicolumn{1}{c}{Pt. Est.} & \\multicolumn{1}{c}{95\\% CI} & \\multicolumn{1}{c}{}  \\\\ 
         \\hline\n 
    ")
}

mytex(tab.1, file.name="CoR_univariable_svycoxph_pretty_"%.%fname.suffix, align="c", include.colnames = F, save2input.only=T, input.foldername=save.results.to,
    col.headers=header,
    longtable=T, 
    label=paste0("tab:CoR_univariable_svycoxph_pretty"), 
    caption.placement = "top", 
    caption=paste0("Inference for Day ", tpeak, " antibody marker covariate-adjusted correlates of risk of ", config.cor$txt.endpoint, " in the ", escape(fname.suffix), " group: Hazard ratios per 10-fold increment in the marker*")
)
tab.cont=tab.1

tab.1.nop12=cbind(paste0(nevents, "/", format(natrisk, big.mark=",")), t(est), t(ci), t(p))
rownames(tab.1.nop12)=all.markers.names.short

# scaled markers
tab.1.scaled=cbind(paste0(nevents, "/", format(natrisk, big.mark=",")), t(est.scaled), t(ci.scaled), t(p), if(show.q) p.2, if(show.q) p.1)
rownames(tab.1.scaled)=all.markers.names.short
tab.1.scaled

if (show.q) {
  header=paste0("\\hline\n 
         \\multicolumn{1}{l}{", toTitleCase(study_name), "} & \\multicolumn{1}{c}{No. cases /}   & \\multicolumn{2}{c}{HR per SD incr.}                     & \\multicolumn{1}{c}{P-value}   & \\multicolumn{1}{c}{q-value}   & \\multicolumn{1}{c}{FWER} \\\\ 
         \\multicolumn{1}{l}{Immunologic Marker}            & \\multicolumn{1}{c}{No. at-risk**} & \\multicolumn{1}{c}{Pt. Est.} & \\multicolumn{1}{c}{95\\% CI} & \\multicolumn{1}{c}{} & \\multicolumn{1}{c}{***} & \\multicolumn{1}{c}{} \\\\ 
         \\hline\n 
    ")
} else {
  header=paste0("\\hline\n 
         \\multicolumn{1}{l}{", toTitleCase(study_name), "} & \\multicolumn{1}{c}{No. cases /}   & \\multicolumn{2}{c}{HR per SD incr.}                     & \\multicolumn{1}{c}{P-value}   \\\\ 
         \\multicolumn{1}{l}{Immunologic Marker}            & \\multicolumn{1}{c}{No. at-risk**} & \\multicolumn{1}{c}{Pt. Est.} & \\multicolumn{1}{c}{95\\% CI} & \\multicolumn{1}{c}{}  \\\\ 
         \\hline\n 
    ")
}

mytex(tab.1.scaled, file.name="CoR_univariable_svycoxph_pretty_scaled_"%.%fname.suffix, align="c", include.colnames = F, save2input.only=T, input.foldername=save.results.to,
    col.headers=header,
    longtable=T, 
    label=paste0("tab:CoR_univariable_svycoxph_pretty_scaled"), 
    caption.placement = "top", 
    caption=paste0("Inference for Day ", tpeak, " antibody marker covariate-adjusted correlates of risk of ", config.cor$txt.endpoint, " in the ", escape(fname.suffix), " group: Hazard ratios per SD increment in the marker*")
)
tab.cont.scaled=tab.1.scaled



###################################################################################################
# make trichotomized markers table

overall.p.1=formatDouble(pvals.adj["tri."%.%names(pvals.cont),"p.FWER"], 3, remove.leading0=F);   overall.p.1=sub("0.000","<0.001",overall.p.1)
overall.p.2=formatDouble(pvals.adj["tri."%.%names(pvals.cont),"p.FDR" ], 3, remove.leading0=F);   overall.p.2=sub("0.000","<0.001",overall.p.2)

# add space
overall.p.1=c(rbind(overall.p.1, NA,NA))
overall.p.2=c(rbind(overall.p.2, NA,NA))


# if "Delta"%.%tpeak%.%"overB" is included, nevents have a problem because some markers may have only two category in the cases

# n cases and n at risk
natrisk = round(c(sapply (all.markers%.%"cat", function(a) aggregate(subset(dat,ph2==1)        [["wt"]], subset(dat,ph2==1        )[a], sum, na.rm=T, drop=F)[,2] )))
nevents = round(c(sapply (all.markers%.%"cat", function(a) aggregate(subset(dat,yy==1 & ph2==1)[["wt"]], subset(dat,yy==1 & ph2==1)[a], sum, na.rm=T, drop=F)[,2] )))
natrisk[is.na(natrisk)]=0
nevents[is.na(nevents)]=0
colSums(matrix(natrisk, nrow=3))
# regression parameters
est=c(rbind(1.00,  sapply(fits.tri, function (fit) if(length(fit)==1) rep(NA,2) else getFormattedSummary(list(fit), exp=T, robust=use.svy, rows=rows, type=1))  ))
ci= c(rbind("N/A", sapply(fits.tri, function (fit) if(length(fit)==1) rep(NA,2) else getFormattedSummary(list(fit), exp=T, robust=use.svy, rows=rows, type=7)) ))
p=  c(rbind("N/A", sapply(fits.tri, function (fit) if(length(fit)==1) rep(NA,2) else getFormattedSummary(list(fit), exp=T, robust=use.svy, rows=rows, type=10)) ))

tab=cbind(
    rep(c("Lower","Middle","Upper"), length(p)/3), 
    paste0(nevents, "/", format(natrisk, big.mark=",",digit=0, scientific=F)), 
    formatDouble(nevents/natrisk, digit=4, remove.leading0=F),
    est, ci, p, overall.p.0, if(show.q) overall.p.2, if(show.q) overall.p.1
)
tmp=rbind(all.markers.names.short, "", "")
rownames(tab)=c(tmp)
tab
tab.cat=tab[1:(nrow(tab)),]
#cond.plac=dat.pla.seroneg[[config.cor$EventTimePrimary]]<=tfinal.tpeak # not used anymore

if(show.q) {
  header=paste0("\\hline\n 
         \\multicolumn{1}{l}{", '', "} & \\multicolumn{1}{c}{Tertile}   & \\multicolumn{1}{c}{No. cases /}   & \\multicolumn{1}{c}{Attack}   & \\multicolumn{2}{c}{Haz. Ratio}                     & \\multicolumn{1}{c}{P-value}   & \\multicolumn{1}{c}{Overall P-}      & \\multicolumn{1}{c}{Overall q-}   & \\multicolumn{1}{c}{Overall} \\\\ 
         \\multicolumn{1}{l}{Immunologic Marker}            & \\multicolumn{1}{c}{}          & \\multicolumn{1}{c}{No. at-risk**} & \\multicolumn{1}{c}{rate}   & \\multicolumn{1}{c}{Pt. Est.} & \\multicolumn{1}{c}{95\\% CI} & \\multicolumn{1}{c}{} & \\multicolumn{1}{c}{value***} & \\multicolumn{1}{c}{value $\\dagger$} & \\multicolumn{1}{c}{FWER} \\\\ 
         \\hline\n 
    ")
} else {
  header=paste0("\\hline\n 
         \\multicolumn{1}{l}{", '', "} & \\multicolumn{1}{c}{Tertile}   & \\multicolumn{1}{c}{No. cases /}   & \\multicolumn{1}{c}{Attack}   & \\multicolumn{2}{c}{Haz. Ratio}                     & \\multicolumn{1}{c}{P-value}   & \\multicolumn{1}{c}{Overall P-}     \\\\ 
         \\multicolumn{1}{l}{Immunologic Marker}            & \\multicolumn{1}{c}{}          & \\multicolumn{1}{c}{No. at-risk**} & \\multicolumn{1}{c}{rate}   & \\multicolumn{1}{c}{Pt. Est.} & \\multicolumn{1}{c}{95\\% CI} & \\multicolumn{1}{c}{} & \\multicolumn{1}{c}{value***} \\\\ 
         \\hline\n 
    ")
}

if (has.plac) {
  add.to.row=list(list(nrow(tab)), # insert at the beginning of table, and at the end of, say, the first table
                  c(paste0(" \n \\multicolumn{8}{l}{} \\\\ \n", 
                           "\n \\multicolumn{2}{l}{Placebo} & ", 
                           paste0(sum(dat.pla.seroneg$yy), "/", format(nrow(dat.pla.seroneg), big.mark=",")), "&",  
                           formatDouble(sum(dat.pla.seroneg$yy)/nrow(dat.pla.seroneg), digit=4, remove.leading0=F), "&",  
                           "\\multicolumn{4}{l}{}  \\\\ \n")
                    #"\\hline\n \\multicolumn{4}{l}{Standard Deviation 1 mcg/mL}\\\\ \n"
                  )
  )
} else {
  add.to.row = NULL
}

# use longtable because this table could be long, e.g. in hvtn705second
mytex(tab[1:(nrow(tab)),], file.name="CoR_univariable_svycoxph_cat_pretty_"%.%fname.suffix, align="c", include.colnames = F, save2input.only=T, input.foldername=save.results.to,
    col.headers=header,        
    add.to.row=add.to.row,
    longtable=T, 
    label=paste0("tab:CoR_univariable_svycoxph_cat_pretty_", fname.suffix), 
    caption.placement = "top", 
    caption=paste0("Inference for Day ", tpeak, " antibody marker covariate-adjusted correlates of risk of ", config.cor$txt.endpoint, " in the ", escape(fname.suffix), " group: Hazard ratios for Middle vs. Upper tertile vs. Lower tertile*")
)


# save two subjects for collate
if (has.plac) {
  save.s.1=paste0(sum(dat.pla.seroneg$yy), "/", format(nrow(dat.pla.seroneg), big.mark=","))
  save.s.2=formatDouble(sum(dat.pla.seroneg$yy)/nrow(dat.pla.seroneg), digit=4, remove.leading0=F)
}


tab.nop12=cbind(
    rep(c("Lower","Middle","Upper"), length(p)/3), 
    paste0(nevents, "/", format(natrisk, big.mark=",",digit=0, scientific=F)), 
    formatDouble(nevents/natrisk, digit=4, remove.leading0=F),
    est, ci, p, overall.p.0
)
rownames(tab.nop12)=c(rbind(all.markers.names.short, "", ""))




###################################################################################################
# multivariate_assays models

if (!is.null(config$multivariate_assays)) {
    if(verbose) print("Multiple regression")
    
    for (a in config$multivariate_assays) {
      for (i in 1:2) {
      # 1: per SD; 2: per 10-fold
        a.tmp=a
        aa=trim(strsplit(a, "\\+")[[1]])
        for (x in aa[!contain(aa, "\\*")]) {
            # replace every variable with Day210x, with or without scale
            a.tmp=gsub(x, paste0(if(i==1) "scale","(Day",tpeak,x,")"), a.tmp) 
        }
        f= update(form.0, as.formula(paste0("~.+", a.tmp)))
        if (use.svy) {
          fits[[a]]=svycoxph(f, design=design.dat) 
        } else {
          fits[[a]]=coxph(f, dat) 
        }
        var.ind=length(coef(fit)) - length(aa):1 + 1
        
        fits=list(fit)
        est=getFormattedSummary(fits, exp=T, robust=use.svy, rows=var.ind, type=1)
        ci= getFormattedSummary(fits, exp=T, robust=use.svy, rows=var.ind, type=7)
        est = paste0(est, " ", ci)
        p=  getFormattedSummary(fits, exp=T, robust=use.svy, rows=var.ind, type=10)
        
        #generalized Wald test for whether the set of markers has any correlation (rejecting the complete null)
        stat=coef(fit)[var.ind] %*% solve(vcov(fit)[var.ind,var.ind]) %*% coef(fit)[var.ind] 
        p.gwald=pchisq(stat, length(var.ind), lower.tail = FALSE)
        
        tab=cbind(est, p)
        tmp=match(aa, colnames(labels.axis))
        tmp[is.na(tmp)]=1 # otherwise, labels.axis["Day"%.%tpeak, tmp] would throw an error when tmp is NA
        rownames(tab)=ifelse(aa %in% colnames(labels.axis), labels.axis["Day"%.%tpeak, tmp], aa)
        colnames(tab)=c(paste0("HR per ",ifelse(i==1,"sd","10 fold")," incr."), "P value")
        tab
        tab=rbind(tab, "Generalized Wald Test"=c("", formatDouble(p.gwald,3, remove.leading0 = F)))
        
        mytex(tab, file.name=paste0("CoR_multivariable_svycoxph_pretty", match(a, config$multivariate_assays), if(i==2) "_per10fold", fname.suffix), align="c", include.colnames = T, save2input.only=T, 
            input.foldername=save.results.to)
      }
    }
    
}



###################################################################################################
# additional_models

if (!is.null(config$additional_models)) {
    if(verbose) print("Additional models")
    
    for (a in config$additional_models) {
        tmp=gsub("tpeak",tpeak,a)
        f= update(Surv(EventTimePrimary, EventIndPrimary) ~1, as.formula(paste0("~.+", tmp)))
        if (use.svy) {
          fits[[a]]=svycoxph(f, design=design.dat) 
        } else {
          fits[[a]]=coxph(f, dat) 
        }

        fits=list(fit)
        est=getFormattedSummary(fits, exp=T, robust=use.svy, type=1)
        ci= getFormattedSummary(fits, exp=T, robust=use.svy, type=7)
        est = paste0(est, " ", ci)
        p=  getFormattedSummary(fits, exp=T, robust=use.svy, type=10)
        
        tab=cbind(est, p)
        colnames(tab)=c("HR", "P value")
        tab
    
        mytex(tab, file.name=paste0("CoR_add_models", match(a, config$additional_models)), align="c", include.colnames = T, save2input.only=T, 
            input.foldername=save.results.to)
    }
    
}


if (attr(config,"config")=="janssen_pooled_EUA") {
    f=Surv(EventTimePrimary, EventIndPrimary) ~ risk_score + as.factor(Region) * Day29pseudoneutid50    
    if (use.svy) {
      fits[[a]]=svycoxph(f, design=design.dat) 
    } else {
      fits[[a]]=coxph(f, dat) 
    }
    var.ind=5:6
    
    fits=list(fit)
    est=getFormattedSummary(fits, exp=T, robust=use.svy, rows=1:6, type=1)
    ci= getFormattedSummary(fits, exp=T, robust=use.svy, rows=1:6, type=7)
    est = paste0(est, " ", ci)
    p=  getFormattedSummary(fits, exp=T, robust=use.svy, rows=1:6, type=10)
    
    #generalized Wald test for whether the set of markers has any correlation (rejecting the complete null)
    stat=coef(fit)[var.ind] %*% solve(vcov(fit)[var.ind,var.ind]) %*% coef(fit)[var.ind] 
    p.gwald=pchisq(stat, length(var.ind), lower.tail = FALSE)
    
    tab=cbind(est, p)
    colnames(tab)=c("HR", "P value")
    tab=rbind(tab, "Generalized Wald Test for Itxn"=c("", formatDouble(p.gwald,3, remove.leading0 = F)))
    tab
    
    mytex(tab, file.name="CoR_region_itxn", align="c", include.colnames = T, save2input.only=T, 
        input.foldername=save.results.to)
        
}


###################################################################################################
# interaction_models

if (!is.null(config$interaction)) {
    if(verbose) print("Interaction models Cox models")    
    itxn.pvals=c()      
    for (ab in config$interaction) {
        tmp=trim(strsplit(ab, " *\\* *")[[1]])
        aold=tmp[1]
        bold=tmp[2]            
        a=paste0("Day",tpeak,aold)
        b=paste0("Day",tpeak,bold)
                
        # fit the interaction model and save regression results to a table
        f= update(form.0, as.formula(paste0("~.+", a," + ",b," + ",a,":",b)))
        if (use.svy) {
          fit=svycoxph(f, design=twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dat)) 
        } else {
          fits[[a]]=coxph(f, dat) 
        }
        
        fits=list(fit)
        est=getFormattedSummary(fits, exp=T, robust=use.svy, type=1)
        ci= getFormattedSummary(fits, exp=T, robust=use.svy, type=7)
        est = paste0(est, " ", ci)
        p=  getFormattedSummary(fits, exp=T, robust=use.svy, type=10)
        # generalized Wald test for whether the set of markers has any correlation (rejecting the complete null)
        var.ind=length(coef(fit))-2:0
        stat=coef(fit)[var.ind] %*% solve(vcov(fit)[var.ind,var.ind]) %*% coef(fit)[var.ind] 
        p.gwald=pchisq(stat, length(var.ind), lower.tail = FALSE)
        # put together the table
        tab=cbind(est, p)
        colnames(tab)=c("HR", "P value")
        tab=rbind(tab, "Generalized Wald Test for Markers"=c("", formatDouble(p.gwald,3, remove.leading0 = F))); tab
        mytex(tab, file.name=paste0("CoR_itxn_",aold,"_",bold), align="c", include.colnames = T, save2input.only=T, input.foldername=save.results.to)
                
        itxn.pvals=c(itxn.pvals, last(getFixedEf(fit)[,"p.value"]))
    }    
    
    names(itxn.pvals)=config$interaction
    itxn.pvals=itxn.pvals[!contain(config$interaction, "ICS4AnyEnv")] # remove the ones with ICS4AnyEnv
    itx.pvals.adj.fdr=p.adjust(itxn.pvals, method="fdr")
    itx.pvals.adj.hol=p.adjust(itxn.pvals, method="holm")
    tab=cbind(itxn.pvals, itx.pvals.adj.hol, itx.pvals.adj.fdr)
    colnames(tab)=c("interaction P value", "FWER", "FDR")
    mytex(tab, file.name="CoR_itxn_multitesting", align="c", include.colnames = T, save2input.only=T, input.foldername=save.results.to)
    
}



###################################################################################################

save(fits.cont.coef.ls, fits.tri.coef.ls, file=paste0(save.results.to, paste0("coxph_fits", fname.suffix, ".Rdata")))

# save.s.1, save.s.2
save (tab.cont, tab.cat, tab.cont.scaled, pvals.adj, file=paste0(save.results.to, paste0("coxph_slopes", fname.suffix, ".Rdata")))
