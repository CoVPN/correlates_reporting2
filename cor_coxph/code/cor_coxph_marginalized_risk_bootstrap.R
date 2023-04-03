###################################################################################################
# bootstrap marginalized risks
# type: 
#    1 for S=s
#    2 for S>=s
#    3 for categorical S
# data: ph1 data
# t: a time point near to the time of the last observed outcome will be defined
marginalized.risk.svycoxph.boot=function(marker.name, type, data, t, B, ci.type="quantile", numCores=1) {  
# marker.name=a; type=2; data=dat.vac.seroneg; t=tfinal.tpeak; B=B; ci.type="quantile"; numCores=1
    
    # store the current rng state 
    save.seed <- try(get(".Random.seed", .GlobalEnv), silent=TRUE) 
    if (class(save.seed)=="try-error") {set.seed(1); save.seed <- get(".Random.seed", .GlobalEnv) } 
    
    data.ph2=subset(data, ph2==1)     
        
    if (comp.risk) {
        f1=lapply(form.0, function(x) update(x, as.formula(paste0("~.+",marker.name))))
    } else {
        f1=update(form.0, as.formula(paste0("~.+",marker.name)))        
    }
    
    # used in both point est and bootstrap
    # many variables are not passed but defined in the scope of marginalized.risk.svycoxph.boot
    fc.1=function(data.ph2, data, categorical.s, n.dean=FALSE){
        if (comp.risk) {
        # competing risk implementation
            newdata=data.ph2
            sapply(ss, function(x) {
                newdata[[marker.name]]=x
                risks = try(pcr2(f1, data.ph2, t, weights=data.ph2$wt, newdata=newdata))
                ifelse (inherits(risks, "try-error"), NA, weighted.mean(risks, data.ph2$wt))
            })
        
        } else {        
        # non-competing risk implementation
            # inline design object b/c it may also throw an error
            fit.risk.1=try(svycoxph(f1, design=twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=data)))        
            if ( !inherits(fit.risk.1, "try-error" )) {
                out=marginalized.risk(fit.risk.1, marker.name, data.ph2, t=t, ss=ss, weights=data.ph2$wt, categorical.s=categorical.s)
                if (n.dean) c(n.dean= last(coef(fit.risk.1)/sqrt(diag(fit.risk.1$var))) * sqrt(1/fit.risk.1$n + 1/fit.risk.1$nevent), out) else out
            } else {
                rep(NA, ifelse(n.dean,1,0)+length(ss))
            }
        }           
    }
    
    fc.2=function(data.ph2){
        if (comp.risk) {
            sapply(ss, function(x) {
                newdata=data.ph2[data.ph2[[marker.name]]>=x, ]
                risks=try(pcr2(f1, newdata, t, weights=newdata$wt))
                ifelse (inherits(risks, "try-error"), NA, weighted.mean(risks, newdata$wt))
            })
        } else {
            out = try(marginalized.risk.threshold (form.0, marker.name, data=data.ph2, weights=data.ph2$wt, t=t, ss=ss))
            if ( !inherits(out, "try-error" )) {
                out
            } else {
                NA # no need to rep, b/c results will be a list when called in bootstrap. for the point est, it is unlikely to be NA
            }
        }
    }    
    
    
    if (type==1) {
        # conditional on S=s (quantitative)
        # don't sort ss or do ss=ss[!duplicated(ss)] because e.g. 15% will be lost and later code depends on that
        ss=sort(c(
            # Lars quantiles so that to be consistent with his analyses, also add every 5% to include s1 and s2 for sensitivity analyses
            report.assay.values(data[[marker.name]][data$EventIndPrimary==1], marker.name.to.assay(marker.name)), 
            # 2.5% and 97.5% as the leftmost and rightmost points 
            wtd.quantile(data[[marker.name]], data$wt, c(0.025,0.05,0.95,0.975)),
            # equally spaced values so that the curves look good  
            seq(min(data[[marker.name]], na.rm=TRUE), max(data[[marker.name]], na.rm=TRUE), length=100)[-c(1,100)],
            # useful for reports
            if (log10(100)>min(data[[marker.name]], na.rm=TRUE) & log10(100)<max(data[[marker.name]], na.rm=TRUE)) log10(100)
        ))
        
        prob = fc.1(data.ph2, data, n.dean=TRUE, categorical.s=F)
        if (!comp.risk) {
            n.dean=prob[1]
            prob=prob[-1]
        } 
        
    } else if (type==2) {
        # conditional on S>=s
        ss=quantile(data[[marker.name]], seq(0,.9,by=0.05), na.rm=TRUE); if(verbose) myprint(ss)
        prob = fc.2(data.ph2)        
       
    } else if (type==3) {
        # conditional on S=s (categorical)
        ss=unique(data[[marker.name]]); ss=sort(ss[!is.na(ss)]); if(verbose) myprint(ss)        
        prob = fc.1(data.ph2, data, n.dean=F, categorical.s=T)
        
    } else if (type==4) {
        # conditional on S=s (quantitative)
        if (comp.risk) {
            stop("need to implement this (like type 1 but coef only)") 
        }else {
            tmp.design=twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=data)
            fit.risk=try(svycoxph(f1, design=tmp.design)) # since we don't need se, we could use coxph, but the weights computed by svycoxph are a little different from the coxph due to fpc
        }
    
    } else stop("wrong type")
    
    # bootstrap
    if(config$case_cohort) ptids.by.stratum=get.ptids.by.stratum.for.bootstrap (data)     
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
           
        if(type==1) {
            # conditional on s
            fc.1(dat.b.ph2, dat.b, categorical.s=F, n.dean=T)
                
        } else if (type==2) {
            # conditional on S>=s
            fc.2(dat.b.ph2)        
            
        } else if (type==3) {
            # conditional on a categorical S
            fc.1(dat.b.ph2, dat.b, n.dean=F, categorical.s=T)
            
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
    if (type==1 & !comp.risk) {
        # the first row is n.dean
        boot.n.dean=res[1,]
        res=res[-1,]
    }
    res=res[,!is.na(res[1,])] # remove NA's
    if (verbose) str(res)
    
    # restore rng state 
    assign(".Random.seed", save.seed, .GlobalEnv)    
        
    if (ci.type=="quantile") {
        ci.band=t(apply(res, 1, function(x) quantile(x, c(.025,.975), na.rm=T)))
    } else {
        stop("only quantile bootstrap CI supported for now")
    }
    
    ret = list(marker=if(type==3) names(prob) else ss, prob=prob, boot=res, lb=ci.band[,1], ub=ci.band[,2], if(type==1 & !comp.risk) n.dean=c(n.dean, boot.n.dean))   
    if (type==1 & !comp.risk) names(ret)[length(ret)]="n.dean" # this is necessary because when using if, that element won't have a name
    ret  
}    



if(!file.exists(paste0(save.results.to, "marginalized.risk.Rdata"))) {    
    cat("make marginalized.risk\n")
    
    # vaccine arm, conditional on continuous S=s
    if (verbose) print("create risks.all.1")
    risks.all.1=lapply(all.markers, function (a) {
        if(verbose) myprint(a)
        marginalized.risk.svycoxph.boot(marker.name=a, type=1, data=dat.vac.seroneg, tfinal.tpeak, B=B, ci.type="quantile", numCores=numCores)                
    })
    
    # vaccine arm, conditional on S>=s
    if (verbose) print("create risks.all.2")
    risks.all.2=lapply(all.markers, function (a) {
        if(verbose) myprint(a)
        marginalized.risk.svycoxph.boot(marker.name=a, type=2, data=dat.vac.seroneg, tfinal.tpeak, B=B, ci.type="quantile", numCores=numCores)        
    }) 
    
    # vaccine arm, conditional on categorical S
    if (verbose) print("create risks.all.3")
    risks.all.3=lapply(all.markers, function (a) {
        if(verbose) myprint(a)
        marginalized.risk.svycoxph.boot(marker.name=a%.%"cat", type=3, data=dat.vac.seroneg, tfinal.tpeak, B=B, ci.type="quantile", numCores=numCores)                
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
# interaction models 

if (!is.null(config$interaction)) {
    if(verbose) print("Interaction models")
    
    if(!file.exists(paste0(save.results.to, "itxn.marginalized.risk.Rdata"))) {    
        cat("make itxn.marginalized.risk\n")
        
        risks.itxn=list()      
        for (ab in config$interaction) {
            tmp=trim(strsplit(ab, " *\\* *")[[1]])
            aold=tmp[1]
            bold=tmp[2]            
            a=paste0("Day",tpeak,aold)
            b=paste0("Day",tpeak,bold)
                    
            # idx=2: only use vaccine arm. idx=1 uses placebo data and structural knowledge; it is commented out and moved to the end of the file
            dat.ph1=dat.vac.seroneg            
            data.ph2=subset(dat.ph1, ph2==1)     
            
            # fit the interaction model and save regression results to a table
            f=as.formula(paste("Surv(EventTimePrimary, EventIndPrimary) ~ RSA + Age + BMI + Riskscore + ",a," + ",b," + ",a,":",b))            
            fit=svycoxph(f, design=twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dat.ph1)) 
            
            # first treat a as the x axis variable, second treat b as the x axis variable
            for (inner.id in 1:2) {
                if (inner.id == 1) {
                    vx=a; vthree=b
                } else {
                    vx=b; vthree=a
                }        
                
                # compute risks at three values of vthree
                three.val=wtd.quantile(dat.ph1[[vthree]][dat.ph1$Trt==1], dat.ph1$wt[dat.ph1$Trt==1], c(.15, .5, .85))                    
                 # compute risks at a sequence of vx values for each of the three vthree values
                ss=sort(c(
                    wtd.quantile(dat.ph1[[vx]], dat.ph1$wt, c(0.025,0.05,0.95,0.975)), # will be included in the table
                    seq(min(dat.ph1[[vx]], na.rm=TRUE), max(dat.ph1[[vx]], na.rm=TRUE), length=100) # equally spaced between min and max so that the curves look good
                ))    
                
#                # special handling code
#                if (attr(config, "config")=="hvtn705second") {
#                    if (inner.id == 1) {
#                        three.val=c(min=-2, wtd.quantile(dat.ph1[[vthree]][dat.ph1$Trt==1], dat.ph1$wt[dat.ph1$Trt==1], c(.5, .9)))
#                    }
#                }
    
    
                # estimate marginalized risks, return a matrix
                prob.ls=sapply (three.val, function(val) {
                        marginalized.risk.cont.2(fit, marker.name  =vx, data=data.ph2, weights=data.ph2$wt, t=tfinal.tpeak, ss=ss, 
                                                      marker.name.2=vthree, s.2=val)
                })
                                
                #### bootstrap
                # store the current rng state
                save.seed <- try(get(".Random.seed", .GlobalEnv), silent=TRUE) 
                if (class(save.seed)=="try-error") {set.seed(1); save.seed <- get(".Random.seed", .GlobalEnv) }         
                
                seeds=1:B; names(seeds)=seeds
                out=mclapply(seeds, mc.cores = numCores, FUN=function(seed) {   
                    seed=seed+560
                    if (verbose>=2) myprint(seed)
                    
#                        if (idx==1) {
#                            # bootstrap vaccine and placebo arm separately
#                            dat.b = rbind(bootstrap.case.control.samples(subset(dat.ph1, Trt==1), seed, delta.name="EventIndPrimary", strata.name="tps.stratum", ph2.name="ph2"),
#                                          subset(dat.ph1, Trt==0)[sample.int(nrow(subset(dat.ph1, Trt==0)), r=TRUE),])         
                    dat.b = bootstrap.case.control.samples(dat.ph1, seed, delta.name="EventIndPrimary", strata.name="tps.stratum", ph2.name="ph2")
                    
                    dat.b.ph2=subset(dat.b, ph2==1)
                    with(dat.b, table(Wstratum, ph2))     
                       
                    # inline design object b/c it may also throw an error
                    fit.b=try(svycoxph(f, design=twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dat.b)))
            
                    if (class (fit.b)[1] != "try-error") {
                        probs=sapply (three.val, function(val) {
                            marginalized.risk.cont.2(fit.b, marker.name  =vx, data=dat.b.ph2, weights=dat.b.ph2$wt, t=tfinal.tpeak, ss=ss, 
                                                            marker.name.2=vthree, s.2=val)
                        })
                    } else {
                        matrix(NA, length(ss), length(three.val))
                    }
                    
                })
                
                # restore rng state 
                assign(".Random.seed", save.seed, .GlobalEnv)    
                
                # organize bootstrap results into a list of 3, each element of which is a matrix of n.s by n.seeds
                res.ls=lapply (1:length(three.val), function(i) {
                    res=sapply(out, function (x) x[,i])
                    res[,!is.na(res[1,])] # remove NA's
                })
                if (verbose) str(res.ls)
                # put lb and ub into matrices
                lb.ls=sapply(res.ls, function (res) t(apply(res, 1, function(x) quantile(x, c(.025)))) )
                ub.ls=sapply(res.ls, function (res) t(apply(res, 1, function(x) quantile(x, c(.975)))) )
                
                risks.itxn[[paste0(vx,vthree)]]=list(marker=ss, prob=prob.ls, boot=res.ls, lb=lb.ls, ub=ub.ls, marker.2=three.val)
            } # end inner.id
    
        }
        save(risks.itxn, file=paste0(save.results.to, "itxn.marginalized.risk.Rdata"))
        
    } else {
        load(paste0(save.results.to, "itxn.marginalized.risk.Rdata"))
    }
}


#                    # idx=1: treat placebo as having known marker values (undetectable)
#                    
#                    # okay to start with dat.mock here b/c no need for trichotomized markers
#                    dat.ph1=subset(dat.mock, ph1)
#                    
#    #                par(mfrow=c(1,2))
#    #                    myboxplot(Day210IgG3gp70.001428.2.42.V1V240delta ~ Trt, dat.ph1, main="IgG3 V1V2")
#    #                    myboxplot(Day210ICS4AnyEnvIFNg_OR_IL2 ~ Trt, dat.ph1, main="CD4")
#    #                
#    #                # there are two non-Zero'0
#    #                head(sort(dat.ph1$Day210IgG3gp70.001428.2.42.V1V240delta[dat.ph1$Trt==0], T))
#    #                head(sort(dat.ph1$Day210ICS4AnyEnvIFNg_OR_IL2[dat.ph1$Trt==0], T))
#    #                subset(dat.ph1, Trt==0 & Day210IgG3gp70.001428.2.42.V1V240delta>0, c(Day210IgG3gp70.001428.2.42.V1V240delta, Day210ICS4AnyEnvIFNg_OR_IL2))
#    #                #     Day210IgG3gp70.001428.2.42.V1V240delta Day210ICS4AnyEnvIFNg_OR_IL2
#    #                #905                                0.875061                   -1.707427
#    #                #1654                               0.977724                   -0.939479
#                            
#                    # set the marker value in the placebo arm to the values representing undetectable
#                    min.ls=c(0, -2); names(min.ls)=c("Day210IgG3gp70.001428.2.42.V1V240delta","ICS4AnyEnvIFNg_OR_IL2") # it is verified that 0 is the minimum value for all V1V2 bAb markers
#                    dat.ph1[[a]][dat.ph1$Trt==0] = min.ls[a]
#                    dat.ph1[[b]][dat.ph1$Trt==0] = min.ls[b]
#                    # set all placebo to a single stratum
#                    dat.ph1$Wstratum[dat.ph1$Trt==0] = min(dat.ph1$Wstratum[dat.ph1$Trt==0])  
#                    # set ph2 and wt to 1 for all placebo
#                    dat.ph1$ph2[dat.ph1$Trt==0] = 1
#                    dat.ph1$wt[dat.ph1$Trt==0] = 1
            
