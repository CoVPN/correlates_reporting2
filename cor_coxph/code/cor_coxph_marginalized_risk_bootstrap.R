###################################################################################################
# bootstrap marginalized risks
# type: 
#    1 for S=s
#    2 for S>=s
#    3 for categorical S
# data: ph1 data
# t: a time point near to the time of the last observed outcome will be defined
marginalized.risk.svycoxph.boot=function(formula, marker.name, type, data, t, B, ci.type="quantile", numCores=1) {  
# formula=form.0; marker.name=a; type=1; data=dat.vac.seroneg; t=tfinal.tpeak; B=B; ci.type="quantile"; numCores=1
    
    # store the current rng state 
    save.seed <- try(get(".Random.seed", .GlobalEnv), silent=TRUE) 
    if (class(save.seed)=="try-error") {set.seed(1); save.seed <- get(".Random.seed", .GlobalEnv) } 
    
    data.ph2=subset(data, ph2==1)     
    
    if (type==1) {
    # conditional on S=s (quantitative)
        ss=sort(c(
            # Lars quantiles so that to be consistent with his analyses 
            # every 5% to include s1 and s2 for sensitivity analyses
            report.assay.values(data[[marker.name]][data$EventIndPrimary==1], marker.name.to.assay(marker.name)), 
            # 2.5% and 97.5% as the leftmost and rightmost points 
            wtd.quantile(data[[marker.name]], data$wt, c(0.025,0.05,0.95,0.975)),
            # equally spaced values so that the curves look good  
            seq(min(data[[marker.name]], na.rm=TRUE), max(data[[marker.name]], na.rm=TRUE), length=100)[-c(1,100)]
        ))
        f1=update(formula, as.formula(paste0("~.+",marker.name)))        
        tmp.design=twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=data)
        fit.risk=try(svycoxph(f1, design=tmp.design)) # since we don't need se, we could use coxph, but the weights computed by svycoxph are a little different from the coxph due to fpc
        prob=marginalized.risk(fit.risk, marker.name, data=data.ph2, ss=ss, weights=data.ph2$wt, t=t, categorical.s=F)    
        # Follmann (2018) ratio of sample sizes
        n.dean = last(coef(fit.risk)/sqrt(diag(fit.risk$var))) * sqrt(1/fit.risk$n + 1/fit.risk$nevent)
        
    } else if (type==2) {
    # conditional on S>=s
        ss=quantile(data[[marker.name]], seq(0,.9,by=0.05), na.rm=TRUE); 
        if(verbose) myprint(ss)
        prob=marginalized.risk.threshold (formula, marker.name, data=data.ph2, weights=data.ph2$wt, t=t, ss=ss)
       
    } else if (type==3) {
    # conditional on S=s (categorical)
        f1=update(formula, as.formula(paste0("~.+",marker.name)))        
        tmp.design=twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=data)
        fit.risk=try(svycoxph(f1, design=tmp.design)) # since we don't need se, we could use coxph, but the weights computed by svycoxph are a little different from the coxph due to fpc
        prob=marginalized.risk(fit.risk, marker.name, data=data.ph2, ss=NULL, weights=data.ph2$wt, t=t, categorical.s=T, verbose=F)        
        
    } else if (type==4) {
    # conditional on S=s (quantitative)
        f1=update(formula, as.formula(paste0("~.+",marker.name)))        
        tmp.design=twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=data)
        fit.risk=try(svycoxph(f1, design=tmp.design)) # since we don't need se, we could use coxph, but the weights computed by svycoxph are a little different from the coxph due to fpc
    
    } else stop("wrong type")
    
    # for use in bootstrap
    if(config$case_cohort) ptids.by.stratum=get.ptids.by.stratum.for.bootstrap (data)     
    
    # bootstrap
    seeds=1:B; names(seeds)=seeds
    out=mclapply(seeds, mc.cores = numCores, FUN=function(seed) {   
    seed=seed+560
        if (verbose>=2) myprint(seed)
    
        if(config$case_cohort) {
            dat.b = get.bootstrap.data.cor (data, ptids.by.stratum, seed) 
        } else {
            dat.b = bootstrap.case.control.samples(data, seed, delta.name="EventIndPrimary", strata.name="tps.stratum", ph2.name="ph2") 
        }        
        dat.b.ph2=subset(dat.b, ph2==1)     
        #hist(dat.b$EventTimePrimaryD14)
        #hist(dat.b$EventTimePrimaryD14[dat.b$EventIndPrimaryD14==1])
        #hist(dat.vac.seroneg$EventTimePrimaryD14)
        #hist(dat.vac.seroneg$EventTimePrimaryD14[dat.vac.seroneg$EventIndPrimaryD14==1])
        #get.marginalized.risk.no.marker(dat.b)
        #get.marginalized.risk.no.marker(dat.vac.seroneg)
           
        if(type==1) {
        # conditional on s
            # inline design object b/c it may also throw an error
            fit.risk.1=try(svycoxph(f1, design=twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dat.b)))
#                    summary(survfit(fit.risk.1))
#                    par(mfrow=c(2,2))
#                    hist(dat.b.ph2$EventTimePrimaryD14[dat.b.ph2$Region==1])
#                    sort(dat.b.ph2$EventTimePrimaryD14[dat.b.ph2$EventIndPrimaryD14==1 & dat.b.ph2$Region==1])
#                    hist(data.ph2$EventTimePrimaryD14[data.ph2$Region==1])
#                    sort(data.ph2$EventTimePrimaryD14[data.ph2$EventIndPrimaryD14==1 & data.ph2$Region==1])
    
            #fit.s=svyglm(f2, tmp.design)      
            if ( class (fit.risk.1)[1] != "try-error" ) {
                n.dean = last(coef(fit.risk.1)/sqrt(diag(fit.risk.1$var))) * sqrt(1/fit.risk.1$n + 1/fit.risk.1$nevent)
                c(n.dean, marginalized.risk(fit.risk.1, marker.name, dat.b.ph2, t=t, ss=ss, weights=dat.b.ph2$wt, categorical.s=F))
            } else {
                rep(NA, 1+length(ss))
            }
            
        } else if (type==2) {
        # conditional on S>=s
            tmp=try(marginalized.risk.threshold (formula, marker.name, data=dat.b.ph2, weights=dat.b.ph2$wt, t=t, ss=ss))
            if (class(tmp) != "try-error" ) tmp else rep(NA,length(ss))
            
        } else if (type==3) {
        # conditional on a categorical S
            fit.risk=try(svycoxph(f1, design=twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dat.b)))
            if ( class (fit.risk)[1] != "try-error" ) {
                marginalized.risk(fit.risk, marker.name, dat.b.ph2, t=t, ss=NULL, weights=dat.b.ph2$wt, categorical.s=T)
            } else {
                rep(NA, 3)
            }
            
        } else if (type==4) {
        # conditional on S=s (quantitative)
            fit.risk.b=try(svycoxph(f1, design=twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dat.b)))
            if ( class (fit.risk.b)[1] != "try-error" ) {
            } else {
                NA
            }
            
        } else stop("wrong type")
        
    })
    res=do.call(cbind, out)
    if (type==1) {
        # the first row is n.dean
        boot.n.dean=res[1,]
        res=res[-1,]
    }
    res=res[,!is.na(res[1,])] # remove NA's
    if (verbose) str(res)
    
    # restore rng state 
    assign(".Random.seed", save.seed, .GlobalEnv)    
    
    if (ci.type=="quantile") {
        ci.band=t(apply(res, 1, function(x) quantile(x, c(.025,.975))))
    } else {
        stop("only quantile bootstrap CI supported for now")
    }
    
    ret = list(marker=if(type==3) names(prob) else ss, prob=prob, boot=res, lb=ci.band[,1], ub=ci.band[,2], if(type==1) n.dean=c(n.dean, boot.n.dean))   
    if (type==1) names(ret)[length(ret)]="n.dean" # this is necessary because when using if, that element won't have a name
    ret  
}    



if(!file.exists(paste0(save.results.to, "marginalized.risk.Rdata"))) {    
    cat("make marginalized.risk\n")
    
    # vaccine arm, conditional on continuous S=s
    if (verbose) print("create risks.all.1")
    risks.all.1=lapply(all.markers, function (a) {
        if(verbose) myprint(a)
        marginalized.risk.svycoxph.boot(formula=form.0, marker.name=a, type=1, data=dat.vac.seroneg, tfinal.tpeak, B=B, ci.type="quantile", numCores=numCores)                
    })
    
    # vaccine arm, conditional on S>=s
    if (verbose) print("create risks.all.2")
    risks.all.2=lapply(all.markers, function (a) {
        if(verbose) myprint(a)
        marginalized.risk.svycoxph.boot(formula=form.0, marker.name=a, type=2, data=dat.vac.seroneg, tfinal.tpeak, B=B, ci.type="quantile", numCores=numCores)        
    }) 
    
    # vaccine arm, conditional on categorical S
    if (verbose) print("create risks.all.3")
    risks.all.3=lapply(all.markers, function (a) {
        if(verbose) myprint(a)
        marginalized.risk.svycoxph.boot(formula=form.0, marker.name=a%.%"cat", type=3, data=dat.vac.seroneg, tfinal.tpeak, B=B, ci.type="quantile", numCores=numCores)                
    })    
    
    save(risks.all.1, risks.all.2, risks.all.3, file=paste0(save.results.to, "marginalized.risk.Rdata"))
    
} else {
    load(paste0(save.results.to, "marginalized.risk.Rdata"))
}
write(ncol(risks.all.1[[1]]$boot), file=paste0(save.results.to, "bootstrap_replicates_"%.%study_name))
#rv$marginalized.risk.S.eq.s=list()
#for (a in assays) rv$marginalized.risk.S.eq.s[[a]] = risks.all.1[[a]][c("marker","prob")]
#rv$marginalized.risk.S.geq.s=list()
#for (a in assays) rv$marginalized.risk.S.geq.s[[a]] = risks.all.2[[a]][c("marker","prob")]



###################################################################################################
# trial-specific ad hoc bootstrapping

if (attr(config, "config") == "hvtn705second") {
    
    #### marginalize interaction model
    if(!file.exists(paste0(save.results.to, "itxn.marginalized.risk.Rdata"))) {    
        cat("make itxn.marginalized.risk\n")
        
        #aa=assays[contain(assays,"V2")][-(1:2)]  
        risks.itxn.1=list()      
        risks.itxn.2=list()      
        aa=c("IgG340mdw_V1V2", "IgG3gp70.001428.2.42.V1V240delta")
        bb=c("ICS4AnyEnvIFNg_OR_IL2", "ICS4JMos1gp120IFNg_OR_IL2", "ICS4JMos1gp41IFNg_OR_IL2", "ICS4JMos2GagIFNg_OR_IL2", "ICS4JMos2RNAseIntIFNg_OR_IL2", "ICS4JMos2Sgp120IFNg_OR_IL2", "ICS4JMos2Sgp41IFNg_OR_IL2")
        for(b in bb) { b=paste0("Day210",b)
        for(a in aa) { a=paste0("Day210",a)
                        
            for (idx in 2:2) { # idx 1 uses placebo data. we are not using that anymore
                if (idx==1) {
                    # idx:1 code can only be run for Day210IgG3gp70.001428.2.42.V1V240delta and ICS4AnyEnvIFNg_OR_IL2
                    # for other a and b, min.ls needs to be defined again
                    
                    # treat placebo as having known marker values (undetectable)
                    
                    # okay to start with dat.mock here b/c no need for trichotomized markers
                    dat.ph1=subset(dat.mock, ph1)
                    
    #                par(mfrow=c(1,2))
    #                    myboxplot(Day210IgG3gp70.001428.2.42.V1V240delta ~ Trt, dat.ph1, main="IgG3 V1V2")
    #                    myboxplot(Day210ICS4AnyEnvIFNg_OR_IL2 ~ Trt, dat.ph1, main="CD4")
    #                
    #                # there are two non-Zero'0
    #                head(sort(dat.ph1$Day210IgG3gp70.001428.2.42.V1V240delta[dat.ph1$Trt==0], T))
    #                head(sort(dat.ph1$Day210ICS4AnyEnvIFNg_OR_IL2[dat.ph1$Trt==0], T))
    #                subset(dat.ph1, Trt==0 & Day210IgG3gp70.001428.2.42.V1V240delta>0, c(Day210IgG3gp70.001428.2.42.V1V240delta, Day210ICS4AnyEnvIFNg_OR_IL2))
    #                #     Day210IgG3gp70.001428.2.42.V1V240delta Day210ICS4AnyEnvIFNg_OR_IL2
    #                #905                                0.875061                   -1.707427
    #                #1654                               0.977724                   -0.939479
                            
                    # set the marker value in the placebo arm to the values representing undetectable
                    min.ls=c(0, -2); names(min.ls)=c("Day210IgG3gp70.001428.2.42.V1V240delta","ICS4AnyEnvIFNg_OR_IL2") # it is verified that 0 is the minimum value for all V1V2 bAb markers
                    dat.ph1[[a]][dat.ph1$Trt==0] = min.ls[a]
                    dat.ph1[[b]][dat.ph1$Trt==0] = min.ls[b]
                    # set all placebo to a single stratum
                    dat.ph1$Wstratum[dat.ph1$Trt==0] = min(dat.ph1$Wstratum[dat.ph1$Trt==0])  
                    # set ph2 and wt to 1 for all placebo
                    dat.ph1$ph2[dat.ph1$Trt==0] = 1
                    dat.ph1$wt[dat.ph1$Trt==0] = 1
                    
                } else {
                    # only use vaccine arm
                    dat.ph1=dat.vac.seroneg
                }
                
                data.ph2=subset(dat.ph1, ph2==1)     
                
                # fit the interaction model
                f=as.formula(paste("Surv(EventTimePrimary, EventIndPrimary) ~ RSA + Age + BMI + Riskscore + ",a," + ",b," + ",a,":",b))            
                fit=svycoxph(f, design=twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dat.ph1)) 
                
                # save regression results to a table
                fits=list(fit)
                est=getFormattedSummary(fits, exp=T, robust=T, type=1)
                ci= getFormattedSummary(fits, exp=T, robust=T, type=13)
                est = paste0(est, " ", ci)
                p=  getFormattedSummary(fits, exp=T, robust=T, type=10)
                # generalized Wald test for whether the set of markers has any correlation (rejecting the complete null)
                var.ind=5:7
                stat=coef(fit)[var.ind] %*% solve(vcov(fit)[var.ind,var.ind]) %*% coef(fit)[var.ind] 
                p.gwald=pchisq(stat, length(var.ind), lower.tail = FALSE)
                # put together the table
                tab=cbind(est, p)
                colnames(tab)=c("HR", "P value")
                tab=rbind(tab, "Generalized Wald Test for Itxn"=c("", formatDouble(p.gwald,3, remove.leading0 = F))); tab
                mytex(tab, file.name=paste0("CoR_itxn_",idx,"_",a,"_",b), align="c", include.colnames = T, save2input.only=T, input.foldername=save.results.to)
                
                # once for each marker on the x axis
                for (inner.id in 1:2) {
                    if (inner.id == 1) {
                        # compute risks at three values of ics
                        three.val=c(min=-2, wtd.quantile(dat.vac.seroneg[[b]][dat.vac.seroneg$Trt==1], dat.vac.seroneg$wt[dat.vac.seroneg$Trt==1], c(.5, .9)))
                        # compute risks at a sequence of cd4 values for each of the three ics values
                        ss=sort(c(
                            wtd.quantile(dat.vac.seroneg[[a]], dat.vac.seroneg$wt, c(0.025,0.05,0.95,0.975)), # will be included in the table
                            seq(min(dat.vac.seroneg[[a]], na.rm=TRUE), max(dat.vac.seroneg[[a]], na.rm=TRUE), length=100) # equally spaced between min and max so that the curves look good
                        ))    
                    } else {
                        # compute risks at three values of v1v2
                        three.val=wtd.quantile(dat.vac.seroneg[[a]][dat.vac.seroneg$Trt==1], dat.vac.seroneg$wt[dat.vac.seroneg$Trt==1], c(.165, .5, .825))
                        # compute risks at a sequence of ics cd4 values for each of the three v1v2 values
                        ss=sort(c(
                            wtd.quantile(dat.vac.seroneg[[b]], dat.vac.seroneg$wt, c(0.025,0.05,0.95,0.975)), # will be included in the table
                            seq(min(dat.vac.seroneg[[b]], na.rm=TRUE), max(dat.vac.seroneg[[b]], na.rm=TRUE), length=100) # equally spaced between min and max so that the curves look good
                        ))    
                        
                    }        
                    
                    # estimate marginalized risks, return a matrix
                    prob.ls=sapply (three.val, function(ics) {
                        marginalized.risk.cont.2(fit, marker.name=ifelse(inner.id==1,a,b), data=data.ph2, weights=data.ph2$wt, t=tfinal.tpeak, ss=ss, marker.name.2=b, s.2=ics)
                    })
                    
                    
                    #### bootstrap
                    
                    # store the current rng state
                    save.seed <- try(get(".Random.seed", .GlobalEnv), silent=TRUE) 
                    if (class(save.seed)=="try-error") {set.seed(1); save.seed <- get(".Random.seed", .GlobalEnv) }         
                    
                    seeds=1:B; names(seeds)=seeds
                    out=mclapply(seeds, mc.cores = numCores, FUN=function(seed) {   
                        seed=seed+560
                        if (verbose>=2) myprint(seed)
                        
                        if (idx==1) {
                            # bootstrap vaccine and placebo arm separately
                            dat.b = rbind(bootstrap.case.control.samples(subset(dat.ph1, Trt==1), seed, delta.name="EventIndPrimary", strata.name="tps.stratum", ph2.name="ph2"),
                                          subset(dat.ph1, Trt==0)[sample.int(nrow(subset(dat.ph1, Trt==0)), r=TRUE),])         
                        } else {
                            dat.b = bootstrap.case.control.samples(dat.ph1, seed, delta.name="EventIndPrimary", strata.name="tps.stratum", ph2.name="ph2")
                        }
                        dat.b.ph2=subset(dat.b, ph2==1)
                        with(dat.b, table(Wstratum, ph2))     
                           
                        # inline design object b/c it may also throw an error
                        fit.b=try(svycoxph(f, design=twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dat.b)))
                
                        if ( class (fit.b)[1] != "try-error" ) {
                            probs=sapply (three.val, function(ics) {
                                marginalized.risk.cont.2(fit.b, marker.name=ifelse(inner.id==1,a,b), data=dat.b.ph2, weights=dat.b.ph2$wt, t=tfinal.tpeak, ss=ss, marker.name.2=ifelse(inner.id==1,b,a), s.2=ics)
                            })
                        } else {
                            matrix(NA, length(ss), length(three.val))
                        }
                        
                    })
                    
                    # restore rng state 
                    assign(".Random.seed", save.seed, .GlobalEnv)    
                    
                    # organize bootstrap results into a list of n.ics, each element of which is a matrix of n.s by n.seeds
                    res.ls=lapply (1:length(three.val), function(i) {
                        res=sapply(out, function (x) x[,i])
                        res[,!is.na(res[1,])] # remove NA's
                    })
                    if (verbose) str(res.ls)
                    # put lb and ub into matrices
                    lb.ls=sapply(res.ls, function (res) t(apply(res, 1, function(x) quantile(x, c(.025)))) )
                    ub.ls=sapply(res.ls, function (res) t(apply(res, 1, function(x) quantile(x, c(.975)))) )
                    
                    if(idx==1) risks.itxn.1[[ifelse(inner.id==1,paste0(a,b),paste0(b,a))]]=list(marker=ss, prob=prob.ls, boot=res.ls, lb=lb.ls, ub=ub.ls, marker.2=three.val)
                    if(idx==2) risks.itxn.2[[ifelse(inner.id==1,paste0(a,b),paste0(b,a))]]=list(marker=ss, prob=prob.ls, boot=res.ls, lb=lb.ls, ub=ub.ls, marker.2=three.val)                    
                } # end inner.id

            } # end idx
        }
        }
        save(risks.itxn.1, risks.itxn.2, file=paste0(save.results.to, "itxn.marginalized.risk.Rdata"))
        
    } else {
        load(paste0(save.results.to, "itxn.marginalized.risk.Rdata"))
    }
}
