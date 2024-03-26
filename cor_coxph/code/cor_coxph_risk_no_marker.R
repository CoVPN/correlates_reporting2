
#fit.ve = coxph(Surv(EventTimePrimary, EventIndPrimary) ~ Trt, subset(dat_proc, ph1==1)) 
#summary(fit.ve)



## these results are close to bootstrap results. they are not used later and only for sanity check
## compute overall risk regardless of markers in both arms by integrating over form.0. 
## the point estimate matche the results from bootstrap
## the variance is asymptotic and still needs to be figured out
#prevs=sapply (c(placebo=0, vaccine=1), function(i) {
#    dat.tmp=subset(dat_proc, Trt==i & Bserostatus==0 & ph1)
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

fname=paste0(save.results.to, "marginalized.risk.no.marker.",fname.suffix,".Rdata")

if(!file.exists(fname)) {    
    cat("Bootstrap marginalized risks using models without markers ...\n")
    
    vacc.only=nrow(dat.pla.seroneg)==0
    
    for (.trt in ifelse(vacc.only, 1, 0):1) {
        dat.tmp=if(.trt==1) dat.vac.seroneg else dat.pla.seroneg
                
        prob = if (TRIAL %in% c("janssen_partA_VL")) {
            mean(sapply(1:10, function(imp) {
                dat.tmp$EventIndOfInterest = ifelse(dat.tmp$EventIndPrimary==1 & dat.tmp[["seq1.variant.hotdeck"%.%imp]]==variant, 1, 0)
                dat.tmp$EventIndCompeting  = ifelse(dat.tmp$EventIndPrimary==1 & dat.tmp[["seq1.variant.hotdeck"%.%imp]]!=variant, 1, 0)
                get.marginalized.risk.no.marker(form.0, dat.tmp, tfinal.tpeak)
            }))
        } else {
            get.marginalized.risk.no.marker(form.0, dat.tmp, tfinal.tpeak)
        }
        
        # bootstrapping
        # store the current rng state 
        save.seed <- try(get(".Random.seed", .GlobalEnv), silent=TRUE) 
        if (class(save.seed)=="try-error") {set.seed(1); save.seed <- get(".Random.seed", .GlobalEnv) }   
        
        if(config$sampling_scheme == 'case_cohort') ptids.by.stratum=get.ptids.by.stratum.for.bootstrap (dat.tmp) 
    
        # if mc.cores is >1 here, the process will be stuck in coxph for some unknown reason
        out=mclapply(1:config$num_boot_replicates, mc.cores = 1, FUN=function(seed) {  
            if (verbose>=2) myprint(seed) 
            if(config$sampling_scheme == 'case_cohort') {
                dat.b = get.bootstrap.data.cor (dat.tmp, ptids.by.stratum, seed) 
            } else if(config$sampling_scheme == 'case_control') {
                dat.b = bootstrap.case.control.samples(dat.ph1=dat.tmp, seed, delta.name="EventIndPrimary", strata.name="tps.stratum", ph2.name="ph2", min.cell.size=0) 
            } else stop("sampling_scheme not supported: "%.%config$sampling_scheme)
            
            prob = if (TRIAL %in% c("janssen_partA_VL")) {
                # if there is no missing variant info in a bootstrap dataset, only need to run the MI code once
                nImp=ifelse(any(with(subset(dat.b, EventIndPrimary==1), is.na(seq1.variant))), 10, 1)
                mean(sapply(1:nImp, function(imp) {
                    dat.b$EventIndOfInterest = ifelse(dat.b$EventIndPrimary==1 & dat.b[["seq1.variant.hotdeck"%.%imp]]==variant, 1, 0)
                    dat.b$EventIndCompeting  = ifelse(dat.b$EventIndPrimary==1 & dat.b[["seq1.variant.hotdeck"%.%imp]]!=variant, 1, 0)
                    get.marginalized.risk.no.marker(form.0, dat.b, tfinal.tpeak)
                }))  
            } else {
                    get.marginalized.risk.no.marker(form.0, dat.b, tfinal.tpeak)
            }
        })
        boot=do.call(cbind, out)
        
        # restore rng state 
        assign(".Random.seed", save.seed, .GlobalEnv)    
        
        if (.trt==0) {
            res.plac.cont=c(est=prob, boot)
            prev.plac=c(res.plac.cont[1], quantile(res.plac.cont[-1], c(.025,.975)))
        } else {
            res.vacc.cont=c(est=prob, boot)
            prev.vacc=c(res.vacc.cont[1], quantile(res.vacc.cont[-1], c(.025,.975)))
        }
    }    
    
    if(!vacc.only) {
        overall.ve = c(1 - res.vacc.cont["est"]/res.plac.cont["est"], quantile(1 - res.vacc.cont[-1]/res.plac.cont[-1], c(0.025, 0.975)))
    } else {
        prev.plac=NA
        res.plac.cont = NA
        overall.ve = NA
    }

    print(cbind(prev.plac, prev.vacc, overall.ve))
    
    save(res.plac.cont, res.vacc.cont, prev.plac, prev.vacc, overall.ve, file=fname)
    
} else {
    load(fname)
}
