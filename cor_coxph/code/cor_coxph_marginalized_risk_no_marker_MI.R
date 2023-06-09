
#fit.ve = coxph(Surv(EventTimePrimary, EventIndPrimary) ~ Trt, subset(dat.mock, ph1==1)) 
#summary(fit.ve)



## these results are close to bootstrap results. they are not used later and only for sanity check
## compute overall risk regardless of markers in both arms by integrating over form.0. 
## the point estimate matche the results from bootstrap
## the variance is asymptotic and still needs to be figured out
#prevs=sapply (c(placebo=0, vaccine=1), function(i) {
#    dat.tmp=subset(dat.mock, Trt==i & Bserostatus==0 & ph1)
#    fit.tmp = coxph(form.0, dat.tmp, model=T) # model=T to make predict possible
#    dat.tmp[[config.cor$EventTimePrimary]]=tfinal.tpeak
#    pred.tmp=predict(fit.tmp, newdata=dat.tmp, type="expected", se.fit=T)    
#    sd.tmp=exp(mean(log(pred.tmp$se.fit)))
#    prev=c(est=NA, "2.5%"=NA, "97.5%"=NA)
#    prev[1] = mean (1 - exp(-pred.tmp$fit))    
#    #prev[2:3] = prev[1] + c(-1,1)*1.96*sd.tmp
#    prev        
#})
#prevs

if(!file.exists(paste0(save.results.to, "marginalized.risk.no.marker.Rdata"))) {    
    if (verbose) print("bootstrap marginalized.risk.no.marker Rdata")
    
    for (i in 1:2) {
        
        # point estimate
        prob=mean(sapply(if(i==1) datasets.vac else datasets.pla, function (dataset) {
            get.marginalized.risk.no.marker(form.0, dataset, tfinal.tpeak)
        }))
        
        # store the current rng state 
        save.seed <- try(get(".Random.seed", .GlobalEnv), silent=TRUE) 
        if (class(save.seed)=="try-error") {set.seed(1); save.seed <- get(".Random.seed", .GlobalEnv) }   
        
        if(config$case_cohort) ptids.by.stratum=get.ptids.by.stratum.for.bootstrap (if(i==1) dat.vac.seroneg else dat.pla.seroneg) 
        
        # if mc.cores is >1 here, the process will be stuck in coxph for some unknown reason
        out=mclapply(1:B, mc.cores = 1, FUN=function(seed) {  
            if (verbose>=2) myprint(seed) 
            if(config$case_cohort) {
                dat.b = get.bootstrap.data.cor (if(i==1) dat.vac.seroneg else dat.pla.seroneg, ptids.by.stratum, seed) 
            } else {
                stop("not supported in marginalised_risk_no_marker")
            }
            
            mean(sapply(1:10, function(imp) {
                dat.b$EventIndOfInterest = ifelse(dat.b$EventIndPrimary==1 & dat.b[["seq1.variant.hotdeck"%.%imp]]==variant, 1, 0)
                dat.b$EventIndCompeting  = ifelse(dat.b$EventIndPrimary==1 & dat.b[["seq1.variant.hotdeck"%.%imp]]!=variant, 1, 0)
                get.marginalized.risk.no.marker(form.0, dat.b, tfinal.tpeak)
            }))
        })
        boot=do.call(cbind, out)
        
        # restore rng state 
        assign(".Random.seed", save.seed, .GlobalEnv)    
        
        if(i==1) {
            res.vacc.cont=c(est=prob, boot)
            prev.vacc=c(res.vacc.cont[1], quantile(res.vacc.cont[-1], c(.025,.975)))
        } else {
            res.plac.cont=c(est=prob, boot)
            prev.plac=c(res.plac.cont[1], quantile(res.plac.cont[-1], c(.025,.975)))
        }
        
    }
    
    overall.ve = c(1 - res.vacc.cont["est"]/res.plac.cont["est"], quantile(1 - res.vacc.cont[-1]/res.plac.cont[-1], c(0.025, 0.975)))

    print(cbind(prev.plac, prev.vacc, overall.ve))
    
    save(res.plac.cont, res.vacc.cont, prev.plac, prev.vacc, overall.ve, file=paste0(save.results.to, "marginalized.risk.no.marker.",region,".Rdata"))
    
} else {
    load(paste0(save.results.to, "marginalized.risk.no.marker.",region,".Rdata"))
}
