cat("Bootstrap controlled risks ...\n")


# vaccine arm, conditional on continuous S=s
if (verbose) print("create risks.all.1")

if (TRIAL=="janssen_partA_VL") {
  fname = paste0(save.results.to, "risks.all.1.", region, ".", variant, ".Rdata")
} else fname = paste0(save.results.to, "risks.all.1.Rdata")

if(!file.exists(fname)) {    
  risks.all.1=lapply(all.markers, function (a) {
    if(verbose) myprint(a)
    marginalized.risk.svycoxph.boot(form.0, marker.name=a, type=1, data=dat.vac.seroneg, t=tfinal.tpeak, B=B, ci.type="quantile", numCores=numCores)                
  })
  save(risks.all.1, file=fname)
} else {
  load(fname)
}
  
  
# vaccine arm, conditional on continuous S>=s
if (verbose) print("create risks.all.2")

if (TRIAL=="janssen_partA_VL") {
  fname = paste0(save.results.to, "risks.all.2.", region, ".", variant, ".Rdata")
} else fname = paste0(save.results.to, "risks.all.2.Rdata")

if(!file.exists(fname)) {    
  risks.all.2=lapply(all.markers, function (a) {
    if(verbose) myprint(a)
    marginalized.risk.svycoxph.boot(form.0, marker.name=a, type=2, data=dat.vac.seroneg, tfinal.tpeak, B=B, ci.type="quantile", numCores=numCores)        
  }) 
  save(risks.all.2, file=fname)
} else {
  load(fname)
}


# vaccine arm, conditional on categorical S
if (verbose) print("create risks.all.3")

if (TRIAL=="janssen_partA_VL") {
  fname = paste0(save.results.to, "risks.all.3.", region, ".", variant, ".Rdata")
} else fname = paste0(save.results.to, "risks.all.3.Rdata")

if(!file.exists(fname)) {    
  risks.all.3=lapply(all.markers, function (a) {
    if(verbose) myprint(a)
    marginalized.risk.svycoxph.boot(form.0, marker.name=a%.%"cat", type=3, data=dat.vac.seroneg, tfinal.tpeak, B=B, ci.type="quantile", numCores=numCores)                
  })    
  save(risks.all.3, file=fname)
} else {
  load(fname)
}

write(ncol(risks.all.1[[1]]$boot), file=paste0(save.results.to, "bootstrap_replicates"))
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
      f= update(form.0, as.formula(paste0("~.+", a," + ",b," + ",a,":",b)))
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
