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
  
  # used in both point est and bootstrap
  # many variables are not passed but defined in the scope of marginalized.risk.svycoxph.boot
  fc.1=function(data.ph2, data, f1, categorical.s, n.dean=FALSE){
    # This is no longer necessary b/c pcr2 is updated to handle the situation when there are no competing events
    # if (is.list(f1)) 
    #   if(all(model.frame(f1[[2]], data.ph2)[[1]][,2]==0)) 
    #     # if there are no competing events, drop competing risk formula
    #     f1 = f1[[1]]
    # 
    if (is.list(f1)) {
      # competing risk implementation
      newdata=data.ph2
      out=sapply(ss, function(x) {
        newdata[[marker.name]]=x
        risks = try(pcr2(f1, data.ph2, t, weights=data.ph2$wt, newdata=newdata))
        ifelse (inherits(risks, "try-error"), NA, weighted.mean(risks, data.ph2$wt))
      })
      if (n.dean) c(NA,out) else out

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
  
  fc.2=function(data.ph2, form.0){
    # This is not enough to ensure there are no competing events b/c subsets of datasets are used in pcr2
    # This is no longer necessary b/c pcr2 is updated to handle the situation when there are no competing events
    # if (is.list(form.0)) 
    #   if(all(model.frame(form.0[[2]], data.ph2)[[1]][,2]==0)) 
    #     # if there are no competing events, drop competing risk formula
    #     form.0 = form.0[[1]]
    # 
    if (is.list(form.0)) {
      sapply(ss, function(x) {
        newdata=data.ph2[data.ph2[[marker.name]]>=x, ]
        risks=try(pcr2(form.0, newdata, t, weights=newdata$wt))
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
  
  if (is.list(form.0)) {
    f1=lapply(form.0, function(x) update(x, as.formula(paste0("~.+",marker.name))))
  } else {
    f1=update(form.0, as.formula(paste0("~.+",marker.name)))        
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
    
    prob = if (TRIAL %in% c("janssen_partA_VL")) {
      # do MI
      rowMeans(sapply(1:10, function(imp) {
        data.ph2$EventIndOfInterest = ifelse(data.ph2$EventIndPrimary==1 & data.ph2[["seq1.variant.hotdeck"%.%imp]]==variant, 1, 0)
        data.ph2$EventIndCompeting  = ifelse(data.ph2$EventIndPrimary==1 & data.ph2[["seq1.variant.hotdeck"%.%imp]]!=variant, 1, 0)
        data$EventIndOfInterest = ifelse(data$EventIndPrimary==1 & data[["seq1.variant.hotdeck"%.%imp]]==variant, 1, 0)
        data$EventIndCompeting  = ifelse(data$EventIndPrimary==1 & data[["seq1.variant.hotdeck"%.%imp]]!=variant, 1, 0)
        fc.1(data.ph2, data, f1, n.dean=TRUE, categorical.s=FALSE)
      }))
    } else {
        fc.1(data.ph2, data, f1, n.dean=TRUE, categorical.s=FALSE)
    }
    
    n.dean=prob[1]
    prob=prob[-1]

  } else if (type==2) {
    # conditional on S>=s
    ss=quantile(data[[marker.name]], seq(0,.9,by=0.05), na.rm=TRUE); if(verbose) myprint(ss)
    prob = if (TRIAL %in% c("janssen_partA_VL")) {
      # do MI
      rowMeans(sapply(1:10, function(imp) {
        data.ph2$EventIndOfInterest = ifelse(data.ph2$EventIndPrimary==1 & data.ph2[["seq1.variant.hotdeck"%.%imp]]==variant, 1, 0)
        data.ph2$EventIndCompeting  = ifelse(data.ph2$EventIndPrimary==1 & data.ph2[["seq1.variant.hotdeck"%.%imp]]!=variant, 1, 0)
        fc.2(data.ph2, form.0)
      }))
    } else {
        fc.2(data.ph2, form.0)
    }
    
    
  } else if (type==3) {
    # conditional on S=s (categorical)
    ss=unique(data[[marker.name]]); ss=sort(ss[!is.na(ss)]); if(verbose) myprint(ss)        
    prob = if (TRIAL %in% c("janssen_partA_VL")) {
      # do MI
      rowMeans(sapply(1:10, function(imp) {
        data.ph2$EventIndOfInterest = ifelse(data.ph2$EventIndPrimary==1 & data.ph2[["seq1.variant.hotdeck"%.%imp]]==variant, 1, 0)
        data.ph2$EventIndCompeting  = ifelse(data.ph2$EventIndPrimary==1 & data.ph2[["seq1.variant.hotdeck"%.%imp]]!=variant, 1, 0)
        data$EventIndOfInterest = ifelse(data$EventIndPrimary==1 & data[["seq1.variant.hotdeck"%.%imp]]==variant, 1, 0)
        data$EventIndCompeting  = ifelse(data$EventIndPrimary==1 & data[["seq1.variant.hotdeck"%.%imp]]!=variant, 1, 0)
        fc.1(data.ph2, data, f1, n.dean=FALSE, categorical.s=TRUE)
      }))
    } else {
      fc.1(data.ph2, data, f1, n.dean=FALSE, categorical.s=TRUE)
    }
    
    
  } else if (type==4) {
    # conditional on S=s (quantitative)
    if (is.list(f1)) {
      stop("need to implement this (like type 1 but coef only)") 
    } else {
      tmp.design=twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=data)
      fit.risk=try(svycoxph(f1, design=tmp.design)) # since we don't need se, we could use coxph, but the weights computed by svycoxph are a little different from the coxph due to fpc
    }
    
  } else stop("wrong type")
  
  
  ###############
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
    
    # if there is no missing variant info in a bootstrap dataset, only need to run the MI code once
    if (TRIAL %in% c("janssen_partA_VL")) 
      nImp = ifelse(any(with(subset(dat.b.ph2, EventIndPrimary==1), is.na(seq1.variant))), 10, 1)
    
    if(type==1) {
      # conditional on s
      if (TRIAL %in% c("janssen_partA_VL")) {
        rowMeans(sapply(1:nImp, function(imp) {
          dat.b.ph2$EventIndOfInterest = ifelse(dat.b.ph2$EventIndPrimary==1 & dat.b.ph2[["seq1.variant.hotdeck"%.%imp]]==variant, 1, 0)
          dat.b.ph2$EventIndCompeting  = ifelse(dat.b.ph2$EventIndPrimary==1 & dat.b.ph2[["seq1.variant.hotdeck"%.%imp]]!=variant, 1, 0)
          dat.b$EventIndOfInterest = ifelse(dat.b$EventIndPrimary==1 & dat.b[["seq1.variant.hotdeck"%.%imp]]==variant, 1, 0)
          dat.b$EventIndCompeting  = ifelse(dat.b$EventIndPrimary==1 & dat.b[["seq1.variant.hotdeck"%.%imp]]!=variant, 1, 0)
          fc.1(dat.b.ph2, dat.b, f1, n.dean=TRUE, categorical.s=FALSE)
        }))
      } else {
        fc.1(dat.b.ph2, dat.b, f1, n.dean=TRUE, categorical.s=FALSE)
      }
      
    } else if (type==2) {
      # conditional on S>=s
      if (TRIAL %in% c("janssen_partA_VL")) {
        rowMeans(sapply(1:nImp, function(imp) {
          dat.b.ph2$EventIndOfInterest = ifelse(dat.b.ph2$EventIndPrimary==1 & dat.b.ph2[["seq1.variant.hotdeck"%.%imp]]==variant, 1, 0)
          dat.b.ph2$EventIndCompeting  = ifelse(dat.b.ph2$EventIndPrimary==1 & dat.b.ph2[["seq1.variant.hotdeck"%.%imp]]!=variant, 1, 0)
          fc.2(dat.b.ph2, form.0)
        }))
      } else {
        fc.2(dat.b.ph2, form.0)
      }
      
    } else if (type==3) {
      # conditional on a categorical S
      if (TRIAL %in% c("janssen_partA_VL")) {
        rowMeans(sapply(1:nImp, function(imp) {
          dat.b.ph2$EventIndOfInterest = ifelse(dat.b.ph2$EventIndPrimary==1 & dat.b.ph2[["seq1.variant.hotdeck"%.%imp]]==variant, 1, 0)
          dat.b.ph2$EventIndCompeting  = ifelse(dat.b.ph2$EventIndPrimary==1 & dat.b.ph2[["seq1.variant.hotdeck"%.%imp]]!=variant, 1, 0)
          dat.b$EventIndOfInterest = ifelse(dat.b$EventIndPrimary==1 & dat.b[["seq1.variant.hotdeck"%.%imp]]==variant, 1, 0)
          dat.b$EventIndCompeting  = ifelse(dat.b$EventIndPrimary==1 & dat.b[["seq1.variant.hotdeck"%.%imp]]!=variant, 1, 0)
          fc.1(dat.b.ph2, dat.b, f1, n.dean=FALSE, categorical.s=TRUE)
        }))
      } else {
        fc.1(dat.b.ph2, dat.b, f1, n.dean=FALSE, categorical.s=TRUE)
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
    ci.band=t(apply(res, 1, function(x) quantile(x, c(.025,.975), na.rm=T)))
  } else {
    stop("only quantile bootstrap CI supported for now")
  }
  
  ret = list(marker=if(type==3) names(prob) else ss, prob=prob, boot=res, lb=ci.band[,1], ub=ci.band[,2], if(type==1) n.dean=c(n.dean, boot.n.dean))   
  if (type==1) names(ret)[length(ret)]="n.dean" # this is necessary because when using if, that element won't have a name
  ret  
}    


################################################################################

cat("make marginalized.risk\n")


# vaccine arm, conditional on continuous S=s
if (TRIAL=="janssen_partA_VL") {
  fname = paste0(save.results.to, "risks.all.1.", region, ".", variant, ".Rdata")
} else fname = paste0(save.results.to, "risks.all.1.Rdata")
if(!file.exists(fname)) {    
  if (verbose) print("create risks.all.1")
  risks.all.1=lapply(all.markers, function (a) {
    if(verbose) myprint(a)
    marginalized.risk.svycoxph.boot(marker.name=a, type=1, data=dat.vac.seroneg, tfinal.tpeak, B=B, ci.type="quantile", numCores=numCores)                
  })
  save(risks.all.1, file=fname)
} else {
  load(fname)
}
  
  
# vaccine arm, conditional on S>=s
if (TRIAL=="janssen_partA_VL") {
  fname = paste0(save.results.to, "risks.all.2.", region, ".", variant, ".Rdata")
} else fname = paste0(save.results.to, "risks.all.2.Rdata")
if(!file.exists(fname)) {    
  if (verbose) print("create risks.all.2")
  risks.all.2=lapply(all.markers, function (a) {
    if(verbose) myprint(a)
    marginalized.risk.svycoxph.boot(marker.name=a, type=2, data=dat.vac.seroneg, tfinal.tpeak, B=B, ci.type="quantile", numCores=numCores)        
  }) 
  save(risks.all.2, file=fname)
} else {
  load(fname)
}


# vaccine arm, conditional on categorical S
if (TRIAL=="janssen_partA_VL") {
  fname = paste0(save.results.to, "risks.all.3.", region, ".", variant, ".Rdata")
} else fname = paste0(save.results.to, "risks.all.3.Rdata")
if(!file.exists(fname)) {    
  if (verbose) print("create risks.all.3")
  risks.all.3=lapply(all.markers, function (a) {
    if(verbose) myprint(a)
    marginalized.risk.svycoxph.boot(marker.name=a%.%"cat", type=3, data=dat.vac.seroneg, tfinal.tpeak, B=B, ci.type="quantile", numCores=numCores)                
  })    
  save(risks.all.3, file=fname)
} else {
  load(fname)
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
