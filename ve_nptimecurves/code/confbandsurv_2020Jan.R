## ---------------------------------------------------------------------------------------
##  Code File: confbandsurv_2019Nov.R
##
##  Description: function to compute a pointwise confidence interval and simultaneous 
##               confidence bands (Parzen et al 1997) for the contrast in Nelson-Aalen based 
##               estimates of the survival curves between two treatment groups
##
##  Function Inputs:
## ---------------------------------------------------------------------------------------
##  time    = minimum of censoring time and failure time
##  delta   = 0 if censored, 1 if event
##  Rx      = 1 or 2 (2 = vaccine, 1 = placebo)
##  N       = number of simulated Gaussian variates for computing simultaneous bands
##  [t1,t2] = interval on which the confidence band is computed
##  ngrid   = number of time-points forming a grid from t1 to t2
##  
##  Choice of contrast is specified by the contrast option
##     contrast = "additive": additive difference y - x
##     contrast = "logRR"   : log hazard ratio log([1-y]/[1-x])
##     contrast = "RR"      : cumulative relative risk (1-y)/(1-x)
##     contrast = "VE"      : vaccine efficacy 1 - (1-y)/(1-x) 
##     
##  NOTE: the contrast function is gsurv(x,y) where x is 
##  P(T>t|placebo) and y is P(T>t|vaccine)
##  
##  bootstrap
##     = TRUE : use bootstrap to compute variance of gsurv(x,y)
##     = FALSE: use Delta method and Greenwood variance estimators
##
##  nsims = the number of bootstrap simulations
##
##  Function Outputs:
## ---------------------------------------------------------------------------------------
##  the output object is a matrix with the following columns:
##
##  timesunique = time points at which confidence bands are computed
##  pointests   = point estimates of gsurv(x,y)
##  lowint      = lower limit of 95% pointwise confidence interval
##  upint       = upper limit of 95% pointwise confidence interval
##  lowband95   = lower limit of 95% simultaneous confidence band
##  upband95    = upper limit of 95% simultaneous confidence band
##  lowband90   = lower limit of 90% simultaneous confidence band
##  upband90    = upper limit of 90% simultaneous confidence band
##  lowband80   = lower limit of 80% simultaneous confidence band
##  upband80    = upper limit of 80% simultaneous confidence band
##
##  Code Edits:
## ---------------------------------------------------------------------------------------
##  22Nov2019  B Prigmore
##     - made minor edits/added comments to improve readability
##     - added nsims as a function argument since it is needed when bootstrap=TRUE
##     - require the survival library when bootstrap=TRUE
##
##  30Jan202  Doug Grove
##     - Fixed ordering of data prior to computation of "Nelson-Aalen" estimator
##     - Fixed the linking of the derived time-grid to the point-estimates/variances
##       computed on the observed times
## ---------------------------------------------------------------------------------------

Confbandsurv <- 
  function(time, delta, Rx, N, t1, t2, ngrid, 
           contrast="additive", bootstrap=FALSE, nsims = NULL)
{

  ## Define the contrast function, derivatives of the contrast function 
  ## in both x and y, and variance functions for both pointwise and 
  ## simultaneous confidence bands.
    
  if (contrast=="additive") {
    # contrast function g
    gsurv <- function(x,y) {
      return(y-x)
    }
    
    # function g1 (derivative of gsurv in x)
    g1surv <- function(x,y) {
      return(-1)
    }
    
    # function g2 (derivative of gsurv in y)
    g2surv <- function(x,y) {
      return(1)
    }
    
    # function v
    vsurv <- function(nsamp1,nsamp2,
                      sa1,sa2,
                      varsa1,varsa2) {
      return((nsamp1*varsa1 + nsamp2*varsa2)**(-0.5))
    }
    
    # IMPORTANT NOTE: vsurv is the variance used for Simultaneous confidence bands, following
    # Parzen et al.
    # for pointwise confidence bands, need vsurv2
    
    vsurv2 <- function(nsamp1,nsamp2,
                       sa1,sa2,
                       varsa1,varsa2) {
      return((varsa1 + varsa2)**(-0.5))
    }
  } # end contrast=="additive"
  
  if (contrast!="additive") {
    # contrast function g
    gsurv <- function(x,y) {
      return(log((1-y)/(1-x)))
    }
    
    # function g1 (derivative of gsurv in x)
    g1surv <- function(x,y) {
      return(1/(1-x))
    }
    
    # function g2 (derivative of gsurv in y)
    g2surv <- function(x,y) {
      return(-1/(1-y))
    }
    
    # function v
    vsurv <- function(nsamp1,nsamp2,
                      sa1,sa2,
                      varsa1,varsa2) {
      return(((nsamp1*varsa1/((1-sa1)^2)) + (nsamp2*varsa2/((1-sa2)^2)))**(-0.5))
    }
    
    vsurv2 <- function(nsamp1,nsamp2,
                       sa1,sa2,
                       varsa1,varsa2) {
      return(((varsa1/((1-sa1)^2)) + (varsa2/((1-sa2)^2)))**(-0.5)) # by the delta method
    }
  
  } # end contrast!="additive"
  
  ### B Prigmore
  ### Not sure why code below is commented out
    
  #if (contrast=="RR")
  #gsurv <- function(x,y) {
  #return((1-y)/(1-x))
  #}
  #
  #if (contrast=="VE")
  #gsurv <- function(x,y) {
  #return(1-(1-y)/(1-x))
  #}
  
  ################################################
  
  ## Order the data by time
    
  nsamp <- length(time)
  nsamp1 <- length(time[Rx==1])
  nsamp2 <- length(time[Rx==2]) 
  
  ## 2020Jan30 Doug Grove - corrected ordering
  ## -----------------------------------------
  ## Ordering needs to be done by both time and delta (decreasing delta, actually) so
  ## observations with events are placed first within a set of observations with the
  ## same follow-up times  
  #ordinds <- order(time)     <- Old code, incorrect
  ordinds <- order(time, -delta)
  
  time <- time[ordinds]
  delta <- delta[ordinds]
  Rx <- Rx[ordinds]
  
  indsRx1 <- c(1:nsamp)[Rx==1]
  indsRx2 <- c(1:nsamp)[Rx==2]
  
  delta1 <- delta[indsRx1]
  delta2 <- delta[indsRx2]
  
  time1 <- time[indsRx1]
  time2 <- time[indsRx2]
  
  
  ## Calculate the indices at which there is a jump (event)
  
  jumpinds1 <- c(1:nsamp1)[delta1==1]
  jumpinds2 <- c(1:nsamp2)[delta2==1]
  
  na1 <- rep(0,nsamp1)
  na2 <- rep(0,nsamp2)
  
  oldsa1 <- rep(0,nsamp1)
  oldsa2 <- rep(0,nsamp2)
  
  varoldsa1 <- rep(0,nsamp1)
  varoldsa2 <- rep(0,nsamp2)
  
  ## Compute Nelson-Aalen ests (na1, na2) for each treatment group 
  ## at each failure time
  
  # Rx = 1:
    risk <- nsamp1
          jna1 <- ifelse(delta1[1]==1,1.0/risk,0.0)
    na1[1] <- jna1
  
  for (i in 2:nsamp1) {
          risk<-nsamp1-i+1
          jna1 <- ifelse(delta1[i]==1,1.0/risk,0.0)
    na1[i] <- na1[i-1]+jna1
  
  }
  
  # Rx = 2:
    risk <- nsamp2
          jna2 <- ifelse(delta2[1]==1,1.0/risk,0.0)
    na2[1] <- jna2
  
  for (i in 2:nsamp2) {
          risk<-nsamp2-i+1
          jna2 <- ifelse(delta2[i]==1,1.0/risk,0.0)
    na2[i] <- na2[i-1]+jna2
  
  }
  
    oldsa1 <- exp(-na1)
    oldsa2 <- exp(-na2)
  
  ## Compute Greenwood variance estimates at each failure time
  
  # Rx = 1:
    risk <- nsamp1
       jna1 <- ifelse(delta1[1]==1,((risk*(risk - delta1[1]))**(-1))*delta1[1],0)
    varoldsa1[1] <- jna1
  
  for (i in 2:nsamp1) {
          risk<-nsamp1-i+1
          jna1 <- ifelse(delta1[i]==1,((risk*(risk - delta1[i]))**(-1))*delta1[i],0)
    varoldsa1[i] <- varoldsa1[i-1]+jna1
  
  }
  
  # Rx = 2:
    risk <- nsamp2
       jna2 <- ifelse(delta2[1]==1,((risk*(risk - delta2[1]))**(-1))*delta2[1],0)
    varoldsa2[1] <- jna2
  
  for (i in 2:nsamp2) {
          risk<-nsamp2-i+1
          jna2 <- ifelse(delta2[i]==1,((risk*(risk - delta2[i]))**(-1))*delta2[i],0)
    varoldsa2[i] <- varoldsa2[i-1]+jna2
  
  }
  
    varoldsa1 <- (oldsa1**2)*varoldsa1
    varoldsa2 <- (oldsa2**2)*varoldsa2
  
  
  ## Create a time grid with ngrid time points between [t1, t2], 
  ## the interval on which the confidence bands will be computed.
  
  ## Find the point estimates and variance estimates for each 
  ## treatment group that correspond to each time point in the grid.
  
  timegrid <- rep(0,ngrid)
  
  sa1 <- rep(0,ngrid)
  sa2 <- rep(0,ngrid)
  varsa1 <- rep(0,ngrid)
  varsa2 <- rep(0,ngrid)
  
        for (i in 1:ngrid) {
  
        tt <- t1+(i-1)*(t2-t1)/(ngrid-1)
        timegrid[i] <- tt
  
        ## 2020Jan30 Doug Grove - corrected computation of 'ind1' and 'ind2'
        ## ------------------------------------------------------------------
        ## 'ind1' and 'ind2' are used to identify which record to pull estimates from
        ## within 'oldsa1', 'oldsa2' and 'varoldsa1', 'varoldsa2'.  This was being done
        ## in such a way that for time point with many observations (ties) the obs. 
        ## identified was NOT the correct one.  It was pointing to the first when it 
        ## should point to the last.  The corrected code immediately follows the incorrect
        ## code below
        #ind1 <- c(1:nsamp1)[time1==max(time1[time1 <= tt])][1]
        #if (is.na(ind1)) { ind1 <- c(1:nsamp1)[time1==min(time1)] }
        w1 <- which( time1 <= tt )
        ind1 <- if ( length(w1) > 0 ) max(w1)
                else  c(1:nsamp1)[time1==min(time1)]
  
        sa1[i] <- oldsa1[ind1]
        varsa1[i] <- varoldsa1[ind1]
  
        #ind2 <- c(1:nsamp2)[time2==max(time2[time2 <= tt])][1]
        #if (is.na(ind2)) { ind2 <- c(1:nsamp2)[time2==min(time2)] }
        w2 <- which( time2 <= tt )
        ind2 <- if ( length(w2) > 0 ) max(w2)
                else c(1:nsamp2)[time2==min(time2)]
  
        sa2[i] <- oldsa2[ind2]
        varsa2[i] <- varoldsa2[ind2]
  
  }
  
  ## Modify the time grid to only include time points at which 
  ## variance estimates for both treatment groups (varsa1 and 
  ## varsa2) are positive
  
  timesunique <- timegrid[varsa1 > 0 & varsa2 > 0] 
  indsunique <- c(1:length(timegrid))[varsa1 > 0 & varsa2 > 0]
  
  sa1 <- sa1[indsunique]
  varsa1 <- varsa1[indsunique]
  sa2 <- sa2[indsunique]
  varsa2 <- varsa2[indsunique]
  
  lenunique <- length(indsunique)
  
  ## Compute the reciprocal of the standard error of gsurv(x,y) 
  ## [the vtilde function in Parzen et al.]
  
  vargsurv <- rep(NA,lenunique)
  vargsurv2 <- rep(NA,lenunique)
  
  ## Analytic variance calculation (the default, performed in any case): 
  for (i in 1:lenunique) {
    vargsurv[i] <- vsurv(nsamp1,nsamp2,sa1[i],sa2[i],varsa1[i],varsa2[i]) 
    vargsurv2[i] <- vsurv2(nsamp1,nsamp2,sa1[i],sa2[i],varsa1[i],varsa2[i])
  }
  vargsurvanalytic <- vargsurv
  vargsurvanalytic2 <- vargsurv2

  ## Boostrap step
  
  if(bootstrap) { 
    require(survival)
    
    Fpp1boot <- matrix(rep(NA,nsims*lenunique),nrow=nsims)
    Fpp2boot <- matrix(rep(NA,nsims*lenunique),nrow=nsims)
    
    diff.ests.boot <- matrix(rep(NA,nsims*lenunique),nrow=nsims)
    logRR.ests.boot <- matrix(rep(NA,nsims*lenunique),nrow=nsims)
    
    for (ii in 1:nsims) {
        boot<-sample(1:nsamp,nsamp,replace=TRUE)
        z.boot<-Rx[boot]
        d.boot<-delta[boot]
        y.boot<-time[boot]
      
      survfit1b<-summary(survfit(Surv(y.boot[z.boot==1],d.boot[z.boot==1])~1))
      survfit2b<-summary(survfit(Surv(y.boot[z.boot==2],d.boot[z.boot==2])~1))
      
      for (i in 1:lenunique) {
        Fpp1b.time.point<-1-min(survfit1b$surv[survfit1b$time<=timesunique[i]])
        Fpp2b.time.point<-1-min(survfit2b$surv[survfit2b$time<=timesunique[i]])
      #  if (Fpp1b.time.point < -999) { Fpp1b.time.point <- 0.000001}
      #  if (Fpp2b.time.point < -999) { Fpp2b.time.point <- 0.000001}
        Fpp1boot[ii,i] <- Fpp1b.time.point
        Fpp2boot[ii,i] <- Fpp2b.time.point
        diff.ests.boot[ii,i]<-Fpp2b.time.point-Fpp1b.time.point
        logRR.ests.boot[ii,i]<-log(Fpp2b.time.point/Fpp1b.time.point)
        
      } # end lenunique
    } # end nsims
    
    for (i in 1:lenunique) {
      if(contrast=="additive") {
        kp <- Fpp1boot[,i] >= 0 & Fpp2boot[,i] >= 0
        #vargsurv[i] <- 1/sqrt(nsamp*var(diff.ests.boot[kp,i]))
        vargsurv[i] <- 1/sqrt(nsamp1**var(Fpp1boot[kp,i]) + nsamp2*var(Fpp2boot[kp,i]))
        vargsurv2[i] <- 1/sqrt(var(Fpp1boot[kp,i]) + var(Fpp2boot[kp,i]))
        # If there is no variability or it is huge, use the analytic variance estimate
        if (is.na(vargsurv[i]) | (!is.na(vargsurv[i]) & (vargsurv[i]==0 | vargsurv[i] > 100))) { vargsurv[i] <- vargsurvanalytic[i] }
        if (is.na(vargsurv2[i]) | (!is.na(vargsurv2[i]) & (vargsurv2[i]==0 | vargsurv2[i] > 100))) { vargsurv2[i] <- vargsurvanalytic2[i] }
        
      } # end if(contrast=="additive")
      
      if(contrast!="additive") {
        kp <- Fpp1boot[,i] >= 0 & Fpp2boot[,i] >= 0
        #vargsurv[i] <- 1/sqrt(nsamp*var(logRR.ests.boot[kp,i]))
        vargsurv[i] <- 1/sqrt(nsamp1**var(log(Fpp1boot[kp,i])) + nsamp2*var(log(Fpp2boot[kp,i])))
        vargsurv2[i] <- 1/sqrt(var(log(Fpp1boot[kp,i])) + var(log(Fpp2boot[kp,i])))
        if (is.na(vargsurv[i]) | (!is.na(vargsurv[i]) & (vargsurv[i]==0 | vargsurv[i] > 100))) { vargsurv[i] <- vargsurvanalytic[i] }
        if (is.na(vargsurv2[i]) | (!is.na(vargsurv2[i]) & (vargsurv2[i]==0 | vargsurv2[i] > 100))) { vargsurv2[i] <- vargsurvanalytic2[i] }
        } # end if(contrast!="additive")
    } # end lenunique
  } # end if(bootstrap)
  
  #vargsurv[vargsurv==Inf] <- vargsurv[min(c(1:lenunique)[vargsurv < 9999999])]
  
  pointests <- rep(0,lenunique)
  lowint <- rep(0,lenunique)
  upint <- rep(0,lenunique)
  lowband95 <- rep(0,lenunique)
  upband95 <- rep(0,lenunique)
  lowband90 <- rep(0,lenunique)
  upband90 <- rep(0,lenunique)
  lowband80 <- rep(0,lenunique)
  upband80 <- rep(0,lenunique)
  
  alpha <- c(0.05,0.10,0.20)
  
  critvalint <- -qnorm(alpha[1]/2)
  
  ## 'critvalband' is a vector with 3 components, each being a critical value for
  ## construction of simultaneous CI at (1-alpha)*100% confidence level
  
  critvalband <- Critvaluesurv(alpha,N,
                               jumpinds1,jumpinds2,
                               nsamp,nsamp1,nsamp2,
                               timesunique,time1,time2,
                               delta1,delta2,
                               sa1,sa2,
                               varsa1,varsa2,
                               lenunique,
                               t1,t2,
                               g1surv,g2surv,vargsurv)
  
  critvalband95 <- critvalband[1]
  critvalband90 <- critvalband[2]
  critvalband80 <- critvalband[3]
  
  ## Compute point estimates for gsurv(x,y) at time points in the time grid.
  ## Compute a 95% pointwise confidence band and simultaneous confidence bands 
  ## at (1-alpha)*100% confidence levels for alpha = 0.05, 0.10, and 0.20.
  
  for (i in 1:lenunique) {
  
    tt<-timesunique[i]
  
    pointests[i] <- gsurv(sa1[i],sa2[i])
  
    x <- pointests[i]
  #  if (contrast=="RR" | contrast=="VE") {x <- log((1-sa2[i])/(1-sa1[i])) }
   
   if(!bootstrap) { 
    vval <- vsurv(nsamp1,nsamp2,sa1[i],sa2[i],varsa1[i],varsa2[i]) 
    vval2 <- vsurv2(nsamp1,nsamp2,sa1[i],sa2[i],varsa1[i],varsa2[i]) 
  }
   if(bootstrap) { vval <- vargsurv[i] 
                   vval2 <- vargsurv2[i]  }
  
    lowint[i] <-    x - critvalint/vval2
    upint[i] <-     x + critvalint/vval2
    lowband95[i] <- x - ((nsamp**(-1/2))*critvalband95)/vval
    upband95[i] <-  x + ((nsamp**(-1/2))*critvalband95)/vval
    lowband90[i] <- x - ((nsamp**(-1/2))*critvalband90)/vval
    upband90[i] <-  x + ((nsamp**(-1/2))*critvalband90)/vval
    lowband80[i] <- x - ((nsamp**(-1/2))*critvalband80)/vval
    upband80[i] <-  x + ((nsamp**(-1/2))*critvalband80)/vval
  
    if(contrast=="RR") {
    pointests[i] <- exp(pointests[i]) 
    lowint[i] <- exp(lowint[i])
    upint[i] <- exp(upint[i])
    lowband95[i] <- exp(lowband95[i])
    upband95[i] <- exp(upband95[i])
    lowband90[i] <- exp(lowband90[i])
    upband90[i] <- exp(upband90[i])
    lowband80[i] <- exp(lowband80[i])
    upband80[i] <- exp(upband80[i])
    }
  
    if(contrast=="VE") {
    pointests[i] <- 1-exp(pointests[i])
    lowint.tmp <- lowint[i]
    lowint[i] <- 1-exp(upint[i])
    upint[i] <- 1-exp(lowint.tmp)
    lowband95.tmp <- lowband95[i]
    lowband95[i] <- 1-exp(upband95[i])
    upband95[i] <- 1-exp(lowband95.tmp)
    lowband90.tmp <- lowband90[i]
    lowband90[i] <- 1-exp(upband90[i])
    upband90[i] <- 1-exp(lowband90.tmp)
    lowband80.tmp <- lowband80[i]
    lowband80[i] <- 1-exp(upband80[i])
    upband80[i] <- 1-exp(lowband80.tmp)
    }
  }
  
  mat <- cbind(timesunique,pointests,lowint,upint,lowband95,upband95,
              lowband90,upband90,lowband80,upband80)
  
  return(mat)

} # end Confbandsurv( )


##########################################################################################
# Sub-functions:
# Parzen, Wei and Ying: Simultaneous Confidence Bands for the Difference of Two Survival
# Functions (SJS, 1997)c
##########################################################################################

##############
# function U
##############

Usurv <- function(N,t,inds,nsamp,time,sa,t1,t2) {
  
  sum <- rep(0,N)
  
  for (j in inds) {
  
    x <- length(time[time >= time[j]])
    
    rec.atrisk <- ifelse(x > 0,1/x,0)
    
    sum <- sum + 
    (ifelse(time[j] >= t1 & time[j] <= t2,1,0)*rec.atrisk*
    ifelse(time[j] <= t,1,0)*rnorm(N))

  } # end inds

  return(-(nsamp**(0.5))*sa*sum)

} # end Usurv


###################
# function Vtilde
###################

Vtildesurv <- function(N,i,t,
                       inds1,inds2,
                       nsamp,nsamp1,nsamp2,
                       time1,time2,
                       sa1,sa2,
                       varsa1,varsa2,
                       t1,t2,
                       g1surv,g2surv,vargsurv)  {

  return(vargsurv[i]*
         (g2surv(sa1[i],sa2[i])*Usurv(N,t,inds2,nsamp,time2,sa2[i],t1,t2)
        + g1surv(sa1[i],sa2[i])*Usurv(N,t,inds1,nsamp,time1,sa1[i],t1,t2)))

} # end vtildesurv

#####################
# function Gtilde
#####################

Gtildesurv <- function(N,inds1,inds2,
                       nsamp,nsamp1,nsamp2,
                       timesunique,time1,time2,
                       delta1,delta2,
                       sa1,sa2,
                       varsa1,varsa2,
                       lenunique,
                       t1,t2,
                       g1surv,g2surv,vargsurv) {

  max <- rep(0,N)
  
  for (i in 1:lenunique) {
  
    tt<-timesunique[i]
    
    x <- abs(Vtildesurv(N,i,tt,inds1,inds2,
                        nsamp,nsamp1,nsamp2,
                        time1,time2,sa1,sa2,
                        varsa1,varsa2,t1,t2,
                        g1surv,g2surv,vargsurv))
    
    max <- ifelse(x > max,x,max)
  
  } # end lenunique
  
  return(max)

} # end Gtildesurv

#########################
# function Critvaluesurv
#########################

Critvaluesurv <- function(alpha,N,
                          inds1,inds2,
                          nsamp,nsamp1,nsamp2,
                          timesunique,time1,time2,
                          delta1,delta2,
                          sa1,sa2,
                          varsa1,varsa2,
                          lenunique,
                          t1,t2,
                          g1surv,g2surv,vargsurv) {

  Gtildevect <- Gtildesurv(N,inds1,inds2,
                           nsamp,nsamp1,nsamp2,
                           timesunique,time1,time2,
                           delta1,delta2,
                           sa1,sa2,
                           varsa1,varsa2,
                           lenunique,
                           t1,t2,
                           g1surv,g2surv,vargsurv)
  
  return(sort(Gtildevect)[floor((1-alpha)*N)])

} # end Critvaluesurv

