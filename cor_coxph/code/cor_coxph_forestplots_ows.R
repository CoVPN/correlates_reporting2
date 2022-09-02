if(verbose) print("forest plots")

###################################################################################################
# forest plots for different phase one baseline strata subgroups
#  "Age >= 65",
#  "Age < 65, At risk",
#  "Age < 65, Not at risk"
fits.all=list()
for (a in assays) {
    fits=list()
    
    fit=svycoxph(update(form.0, as.formula(paste0("~.+Day",tpeak, a))), design=design.vacc.seroneg) 
    fits[[1]]=fit
    
    for (k in 1:max(dat.mock$Bstratum)) {
        if(sum (subset(dat.vac.seroneg, Bstratum==k, yy))<=2) {
            # 0-2 cases
            fits[[k+1]]=NA
        } else {
            design<-try(twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=subset(dat.vac.seroneg, Bstratum==k)),silent=TRUE)
            if (inherits(design,"try-error")) design=NULL
            if (k==1) {
                f = update(form.0, as.formula(paste0("~.+Day",tpeak, a)))
            } else {
                f = update(update(form.0, ~.-HighRiskInd), as.formula(paste0("~.+Day",tpeak, a)))
            }
            fits[[k+1]]=run.svycoxph(f, design=design)
        }
    }       
    
    fits.all[[a]]=fits
}

nevents=c(nrow(subset(dat.vac.seroneg, yy==1)),
          sapply(1:max(dat.mock$Bstratum), function (k) nrow(subset(dat.vac.seroneg, yy==1 & Bstratum==k))) 
)

rv$fr.1=list(nevents=nevents)
for (a in assays) {
    #width and height decide margin
    # to make CI wider, make width bigger and graphwidth larger # onefile has to be F otherwise there will be an empty page inserted
    mypdf(onefile=F, width=10,height=3, file=paste0(save.results.to, "hr_forest_", a, "_", study_name)) 
        fits = fits.all[[a]]
        names(fits)=c("All baseline negative, vaccine", "      "%.%Bstratum.labels)
        est.ci = sapply(fits, function (fit) {
            if (length(fit)==1) return (rep(NA,4))
            tmp=getFixedEf(fit, exp=T, robust=T)
            tmp[nrow(tmp),c("HR", "(lower", "upper)", "p.value")]
        })
        theforestplot(point.estimates=est.ci[1,], lower.bounds=est.ci[2,], upper.bounds=est.ci[3,], p.values=NA, graphwidth=unit(120, "mm"), fontsize=1.2,
            table.labels = c("Group", "HR (95% CI)","No. Events"), group=colnames(est.ci), decimal.places=2, nEvents=nevents, title=paste0(labels.assays.long["Day"%.%tpeak,a]))
    dev.off()
    
    rv$fr.1[[a]]=est.ci
}


###################################################################################################
# forest plots for subpopulations

# different fits may use different formulae
fits.all.2=vector("list", length(assays));names(fits.all.2)=assays
i=1

# whole population
for (a in assays) {
    f= update(form.0, as.formula(paste0("~.+Day",tpeak, a)))
    fits.all.2[[a]]=list()        
    fits.all.2[[a]][[i]]=run.svycoxph(f, design=design.vacc.seroneg) 
}
i=i+1

# senior
design.1<-twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=get.dat.with.no.empty(subset(dat.vac.seroneg, Senior==1)))
design.2<-twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=get.dat.with.no.empty(subset(dat.vac.seroneg, Senior==0)))        
for (a in assays) {
    f= update(form.0, as.formula(paste0("~.+Day",tpeak, a)))
    fits.all.2[[a]][[i]]=run.svycoxph(f, design=design.1) 
    fits.all.2[[a]][[i+1]]=run.svycoxph(f, design=design.2) 
}
i=i+2

# high risk
design.1<-twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=get.dat.with.no.empty(subset(dat.vac.seroneg, HighRiskInd==1)))
design.2<-twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=get.dat.with.no.empty(subset(dat.vac.seroneg, HighRiskInd==0)))
for (a in assays) {
    f= update(update(form.0, ~.-HighRiskInd), as.formula(paste0("~.+Day",tpeak, a)))
    fits.all.2[[a]][[i]]=run.svycoxph(f, design=design.1) 
    fits.all.2[[a]][[i+1]]=run.svycoxph(f, design=design.2) 
}
i=i+2

# MinorityInd makes sense for COVE and US in ensemble
if(study_name %in% c("COVE","MockCOVE") | startsWith(attr(config, "config"), "janssen_pooled") | startsWith(attr(config, "config"), "janssen_na")) {
    if(study_name %in% c("COVE","MockCOVE")) {
        design.1<-twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=get.dat.with.no.empty(subset(dat.vac.seroneg, MinorityInd==1)))
        design.2<-twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=get.dat.with.no.empty(subset(dat.vac.seroneg, MinorityInd==0)))
    } else {
        design.1<-twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=get.dat.with.no.empty(subset(dat.vac.seroneg, MinorityInd==1 & Region==0)))
        design.2<-twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=get.dat.with.no.empty(subset(dat.vac.seroneg, MinorityInd==0 & Region==0)))
    }
    
    for (a in assays) {
        f= update(update(form.0, ~.-as.factor(Region)-MinorityInd), as.formula(paste0("~.+Day",tpeak, a)))
        fits.all.2[[a]][[i]]=run.svycoxph(f, design=design.1) 
        fits.all.2[[a]][[i+1]]=run.svycoxph(f, design=design.2) 
    }
    i=i+2
}

# gender
design.1<-twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=get.dat.with.no.empty(subset(dat.vac.seroneg, Sex==1)))
design.2<-twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=get.dat.with.no.empty(subset(dat.vac.seroneg, Sex==0)))
for (a in assays) {
    f= update(form.0, as.formula(paste0("~.+Day",tpeak, a)))
    fits.all.2[[a]][[i]]=run.svycoxph(f, design=design.1) 
    fits.all.2[[a]][[i+1]]=run.svycoxph(f, design=design.2) 
}
i=i+2

# HIV infection
if (study_name %in% c("ENSEMBLE", "MockENSEMBLE")) {
    design.1<-twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=get.dat.with.no.empty(subset(dat.vac.seroneg, HIVinfection==1)))
    design.2<-twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=get.dat.with.no.empty(subset(dat.vac.seroneg, HIVinfection==0)))        
    for (a in assays) {
        if (startsWith(a, "pseudoneut")) next
        f= update(form.0, as.formula(paste0("~.+Day",tpeak, a)))
        fits.all.2[[a]][[i]]=run.svycoxph(f, design=design.1) 
        fits.all.2[[a]][[i+1]]=run.svycoxph(f, design=design.2) 
    }
}


age.threshold=switch(study_name,COVE=65,MockCOVE=65,AZD1222=65,ENSEMBLE=60,MockENSEMBLE=60,stop("age threshold undefined"))

for (a in assays) {    
    names(fits.all.2[[a]])=c("All Vaccine", 
                             "Age >= "%.%age.threshold, "Age < "%.%age.threshold, 
                             "At risk", "Not at risk", 
                             if (study_name %in% c("COVE","MockCOVE") ) c("Comm. of color", "White Non-Hispanic"),
                             if (startsWith(attr(config, "config"), "janssen_pooled") | startsWith(attr(config, "config"), "janssen_na")) c("Comm. of color (US)", "White Non-Hispanic (US)"),
                             "Men", "Women",
                             if (study_name %in% c("ENSEMBLE","MockENSEMBLE") & !startsWith(a, "pseudoneut")) c("HIV infection Yes", "HIV infection No")
    )
}    


rv$fr.2=list(nevents=nevents)
for (a in assays) {
    nevents=c(nrow(subset(dat.vac.seroneg, yy==1)),
              nrow(subset(dat.vac.seroneg, yy==1 & Senior==1)), 
              nrow(subset(dat.vac.seroneg, yy==1 & Senior==0)), 
              nrow(subset(dat.vac.seroneg, yy==1 & HighRiskInd==1)), 
              nrow(subset(dat.vac.seroneg, yy==1 & HighRiskInd==0)), 
              # MinorityInd also makes sense for US in ensemble
              if(study_name %in% c("COVE","MockCOVE") ) nrow(subset(dat.vac.seroneg, yy==1 & MinorityInd==1)), 
              if(study_name %in% c("COVE","MockCOVE") ) nrow(subset(dat.vac.seroneg, yy==1 & MinorityInd==0)), 
              if(startsWith(attr(config, "config"), "janssen_pooled") | startsWith(attr(config, "config"), "janssen_na")) nrow(subset(dat.vac.seroneg, yy==1 & MinorityInd==1 & Region==0)), 
              if(startsWith(attr(config, "config"), "janssen_pooled") | startsWith(attr(config, "config"), "janssen_na")) nrow(subset(dat.vac.seroneg, yy==1 & MinorityInd==0 & Region==0)), 
              nrow(subset(dat.vac.seroneg, yy==1 & Sex==1)), 
              nrow(subset(dat.vac.seroneg, yy==1 & Sex==0)),
              if (study_name %in% c("ENSEMBLE","MockENSEMBLE") & !startsWith(a, "pseudoneut")) { c(
                  nrow(subset(dat.vac.seroneg, yy==1 & HIVinfection==1)), 
                  nrow(subset(dat.vac.seroneg, yy==1 & HIVinfection==0)))
              }
    )
    
    fits = fits.all.2[[a]]
    est.ci = sapply(fits, function (fit) {
        if (length(fit)==1) return (rep(NA,4))
        if (last(coef(fit))>log(1000)) return (rep(NA,4)) # cannot be true
        tmp=getFixedEf(fit, exp=T, robust=T)
        tmp[nrow(tmp),c("HR", "(lower", "upper)", "p.value")]
    })
    mypdf(onefile=F, file=paste0(save.results.to, "hr_forest_marginal_", a, "_", study_name), width=10,height=.45*ncol(est.ci)) #width and height decide margin
        theforestplot (lineheight=unit(.75,"cm"), point.estimates=est.ci[1,], lower.bounds=est.ci[2,], upper.bounds=est.ci[3,], p.values=NA, graphwidth=unit(120, "mm"), fontsize=1.2,
            table.labels = c("Group (Baseline Negative)", "HR (95% CI)","No. Events"), group=colnames(est.ci), decimal.places=2, nEvents=nevents, title=paste0(labels.assays.long["Day"%.%tpeak,a]))
    dev.off()    
    
    rv$fr.2[[a]]=est.ci
}



###################################################################################################
# forest plots for different countries and regions

if (startsWith(attr(config, "config"), "janssen_pooled") | startsWith(attr(config, "config"), "janssen_la")) {
# la or pooled
# not na or sa
    
    is.pooled=startsWith(attr(config, "config"), "janssen_pooled")
    
    if(verbose) print("forest plots for different countries and regions")
    
    regions=  get("regions.ENSEMBLE")
    countries=get("countries.ENSEMBLE")
    labels.regions=  get("labels.regions.ENSEMBLE")
    labels.countries=get("labels.countries.ENSEMBLE")
#    regions=  get("regions."  %.%study_name_code)
#    countries=get("countries."%.%study_name_code)
#    labels.regions=  get("labels.regions."  %.%study_name_code)
#    labels.countries=get("labels.countries."%.%study_name_code)
    
    countries.1 = countries[countries %in% unique(dat.vac.seroneg$Country)]
    
    designs=list(design.vacc.seroneg); names(designs)="All Vaccine"
    designs=append(designs, lapply(countries.1, function (i) {
        out=try(twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=get.dat.with.no.empty(subset(dat.vac.seroneg, Country==i))), silent=TRUE)
        if (!inherits(out,"try-error")) out else NULL
    }))
    # add Latin America after united states if pooled
    if (is.pooled) {
        designs = append(designs, lapply(regions[2], function (i) twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=get.dat.with.no.empty(subset(dat.vac.seroneg, Region==i)))), after=2)
    }
    fits.all.3=lapply(assays, function(a) {
        f=  update(update(form.0, ~.-as.factor(Region)), as.formula(paste0("~.+Day",tpeak, a)))
        f.1=update(form.0, as.formula(paste0("~.+Day",tpeak, a))) # keep Region for All vaccine
        out=lapply(1:length(designs), function (i) if(i==1) run.svycoxph(f.1, design=designs[[1]]) else run.svycoxph(f, design=designs[[i]]) )
        names(out)=names(designs)
        out
    })
    
    nevents=nrow(subset(dat.vac.seroneg, yy==1))
    nevents=c(nevents, sapply(countries.1, function(i) nrow(subset(dat.vac.seroneg, yy==1 & Country==i))))
    if (is.pooled) {
        nevents=append(nevents, sapply(regions[2],   function(i) nrow(subset(dat.vac.seroneg, yy==1 & Region==i ))), after=2)
    }
    
    select=sapply(designs,function(x) !is.null(x))
    designs=designs[select]
    nevents=nevents[select]
    
    
    rv$fr.3=list(nevents=nevents)
    for (a in assays) {
        fits = fits.all.3[[a]]
        est.ci = sapply(fits[select], function (fit) {
            if (length(fit)==1) return (rep(NA,4)) # fit is NA
            tmp=getFixedEf(fit, exp=T, robust=T)
            if (tmp[nrow(tmp),1]>20)  return (rep(NA,4)) # coefficient is basically infinite
            tmp[nrow(tmp),c("HR", "(lower", "upper)", "p.value")]
        })
        # move latin american country names to the right by two spaces
        if (is.pooled) colnames(est.ci)[4:9] = "    "%.%colnames(est.ci)[4:9]
        mypdf(onefile=F, file=paste0(save.results.to, "hr_forest_countries_", a, "_", study_name), width=10,height=4.5) #width and height decide margin
            theforestplot (lineheight=unit(.75,"cm"), point.estimates=est.ci[1,], lower.bounds=est.ci[2,], upper.bounds=est.ci[3,], p.values=NA, graphwidth=unit(120, "mm"), fontsize=1.2,
                table.labels = c("Group (Baseline Negative)", "HR (95% CI)","No. Events"), group=colnames(est.ci), decimal.places=2, nEvents=nevents, title=paste0(labels.assays.long["Day"%.%tpeak,a]))
        dev.off()        
        rv$fr.3[[a]]=est.ci
    }

} # end if for countries/regions forest plots
