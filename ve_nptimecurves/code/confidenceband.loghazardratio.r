#####################################################################################
#  Function confidenceband.loghazardratio
#
#  The function estimates the log hazard ratio (vaccine/placebo) over time and
#  computes 95% pointwise and 80%, 90%, and 95% simultaneous confidence bands
#  for the log hazard ratio over time.  The method is described in Gilbert, Wei, 
#  Kosorok, and Clemens, Biometrics 2002, in press.  Kernel estimation as described 
#  in ABGK (Andersen, Gill, Borgan, and Keiding, 1993) is used to calculate
#  the log hazard ratio estimate.
######################################################################################
#
# INPUT variables into the function confidenceband.loghazardratio.r:
#
# time     = n-vector of mimimum of censoring time and failure time (n subjects in the trial)
# delta    = n-vector of failure indi#cators: 0 if censored, 1 if event
# Rx       = n-vector of treatment arm: 1 = treatment 1 (e.g., control), 2 = treatment 2 (e.g., active)
#            This variable must be coded as 1 and 2.
# [u1,u2]  = the time interval over which the confidence band is computed
# [t1,t2]  = the time interval over which the confidence band could be computed (at least
#            as wide as [u1,u2]). t1 should be chosen larger than both of the smallest
#            infection times for each group. t2 should be chosen smaller than both of 
#            the largest infection times for each group.
# N        = number of simulated Gaussian variates used for computing the
#            simultaneous confidence bands (N = 1000 is adequate)
# band11   = second derivative bandwidth for treatment 1, as on page 249 of Andersen, Borgan, 
#            Gill, and Keiding (1993), used to calculate the optimal bandwidth band1
#            for the treatment 1 sample. 
#            If band11 = 0, the optimal band11 is calculated using a bootstrap procedure. 
# band1    = the bandwidth for the first sample. The default value is 
#            the optimal bandwidth, calculated using the formula  in ABGK (1993)
# band12   = second derivative bandwidth for treatment 2, as on page 249 of Andersen, Borgan, 
#            Gill, and Keiding (1993), used to calculate the optimal bandwidth band2
#            for the treatment 2 sample. 
#            If band12 = 0, the optimal band12 is calculated using a bootstrap procedure. 
# band2    = the bandwidth for the second sample. The default value is
#            the optimal bandwidth, calculated using the formula  in ABGK (1993)
# ngrid    = the number of time-points spanning [t1,t2] for computing
#            the log hazard ratio and the confidence bands. The default value is 40.
# biasadjust = T means that the hazard rate estimates are adjusted by the asymptotic
#            bias correction term (4.2.27) in ABGK.
# tailsl   = T means that the method described on page 251 of
#            ABGK is used to calculate the log hazard ratio estimate
#            and the confidence intervals/bands in the lower tail [t1,u1]
# tailsu   = T means that the method described on page 251 of
#	     ABGK is used to calculate the log hazard ratio estimate
#	     and the confidence intervals/bands in the upper tail [u2,t2]
# nboot    = Number of bootstrap samples used for the cross-validation
#            procedure for estimating the bandwidth parameters used to
#            estimate the second derivatives band11 of the hazard functions
# Nv       = 0 implies that an asymptotic variance formula will be used for the 
#            variance function v()
#          > 0 implies that the simulation procedure will be used to calculate
#            the variance function v(). In this case, Nv is the number of iterations
#            used to calculate the simulated estimate of v.
###################################################################################
#
# OUTPUT of confidenceband.loghazardratio.r, which is used for plotting 
#
# The function outputs an ngrid by 14 matrix, with the following columns:
#
# Column 1: Grid of follow-up times
# Column 2: Log hazard ratio estimate (treatment 2/treatment 1) on the grid of times
# Column 3: lower limit of pointwise 95% CI for the log hazard ratio
# Column 4: upper limit of pointwise 95% CI for the log hazard ratio
# Column 5: lower limit of simultaneous 95% CI for the log hazard ratio
# Column 6: upper limit of simultaneous 95% CI for the log hazard ratio
# Column 7: lower limit of simultaneous 90% CI for the log hazard ratio
# Column 8: upper limit of simultaneous 90% CI for the log hazard ratio
# Column 9: lower limit of simultaneous 80% CI for the log hazard ratio
# Column 10: upper limit of simultaneous 80% CI for the log hazard ratio
# Column 11: Point estimate of the Rx=1 hazard function over the grid of times
# Column 12: Point estimate of the Rx=2 hazard function over the grid of times
# Column 13: Estimated variance of the Rx=1 hazard estimate over the grid of times
# Column 14: Estimated variance of the Rx=2 hazard estimate over the grid of times
###################################################################################
#
# The function Plotconfidenceband.loghazratio is used to call confidenceband.loghazratio
# and to plot the confidence bands
###################################################################################

confidenceband.loghazardratio <- function(time,delta,Rx,u1,u2,t1,t2,N=1000,band11=0,
band1=0,band12=0,band2=0,ngrid=50,biasadjust=T,tailsl=T,tailsu=T,nboot=500,Nv=0) {

nsamp <- length(time)
nsamp1 <- length(time[Rx==1])
nsamp2 <- length(time[Rx==2]) 

# Order the data by time

ordinds <- order(time)

time <- time[ordinds]
delta <- delta[ordinds]
Rx <- Rx[ordinds]

indsRx1 <- c(1:nsamp)[Rx==1]
indsRx2 <- c(1:nsamp)[Rx==2]

delta1 <- delta[indsRx1]
delta2 <- delta[indsRx2]

time1 <- time[indsRx1]
time2 <- time[indsRx2]

# Calculate the constant k2K, used for bias adjustment and for
# calculation of the optimal bandwidths

#k2K <- 0.020202*sum((seq(-1,1,length=100)**2)*epan(seq(-1,1,length=100)))
# k2K = 1/5 for the Epanechnikov kernel
k2K <- 0.20

# Calculate the constant k2Klmod and k2Kumod, used for bias 
# adjustment in the tails

k2Klmod <- function(q) {

return(0.020202*sum((seq(-1,1,length=100)**2)*epanltail(seq(-1,1,length=100),q)))

}

k2Kumod <- function(q) {

return(0.020202*sum((seq(-1,1,length=100)**2)*epanutail(seq(-1,1,length=100),q)))

}


# calculate the indices at which there is a jump (event)

jumpinds1 <- c(1:nsamp1)[delta1==1]
jumpinds2 <- c(1:nsamp2)[delta2==1]

na1 <- rep(0,nsamp1)
na2 <- rep(0,nsamp2)

# Compute Nelson-Aalen ests (na1, na2) at each failure time

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

###################
#  Calculate the optimal bandwidth for each sample

# Calculate int K^2(t)

#intK2 <- 0.02020202*sum(epan(seq(-1,1,length=100))**2)
# Epan kernel special case:
intK2 <- 0.6

# Calculate the integral term with alpha and y, for
# each sample

at.risks1 <- rep(0,length(jumpinds1))

for (j in 1:length(jumpinds1)) {

at.risks1[j] <- length(time1[time1 >= time1[jumpinds1[j]]])

}

delta11 <- sum(ifelse(at.risks1 > 0 & time1[jumpinds1] <= u1,1/(at.risks1**2),0))
delta12 <- sum(ifelse(at.risks1 > 0 & time1[jumpinds1] <= u2,1/(at.risks1**2),0))

at.risks2 <- rep(0,length(jumpinds2))
for (j in 1:length(jumpinds2)) {

at.risks2[j] <- length(time2[time2 >= time2[jumpinds2[j]]])

}

delta21 <- sum(ifelse(at.risks2 > 0 & time2[jumpinds2] <= u1,1/(at.risks2**2),0))
delta22 <- sum(ifelse(at.risks2 > 0 & time2[jumpinds2] <= u2,1/(at.risks2**2),0))

alpy1 <- nsamp1*(delta12 - delta11)
alpy2 <- nsamp2*(delta22 - delta21)

# Calculate the integrated second derivative term, for each sample
# Do this using the formulae on the top of page 249 of ABGK

# Iterate multiple times, for different values of band11, to
# select band11 via cross-validation based on the MISE of the 
# second derivative term

# This bootstrap step is optional, alternatively I've found it reliable
# to simply set band111 = band112 = (u2-u1)/3.  

findband11 <- function(inu1,inu2,innsamp, intimevect, injumpinds, intime, indelta, 
inna, inintK2, inalpy) {

f <- function(band11, nsamp, timevect, jumpinds, time, delta, na, intK2, alpy, u1, u2) 
{
ind <- 1
inds <- rep(NA,length(time))
newtime <- rep(NA,length(time))
newtime[1] <- time[1]
inds[1] <- 1
for (i in 2:nsamp) {
if (time[i]!=time[i-1]) {
newtime[ind] <- time[i-1] 
inds[ind] <- 1
ind <- ind + 1 
}
}

time <- newtime[!is.na(newtime)]
na <- na[inds==1]
na <- na[!is.na(na)]
delta <- delta[inds==1]
delta <- delta[!is.na(delta)]

jumpinds <- c(1:length(na))[delta==1]

timevect <- time[time <= u2 & time >= u1]

temporig <- alphadblprime(timevect,
                     band11,jumpinds,time,na)
MISE <- 0

for (iboot in 1:nboot) { 

# Draw bootstrap samples

boottime <- sample(time,length(time),replace=TRUE)
bootdelta <- sample(delta,length(delta),replace=TRUE)

bootordinds <- order(boottime)
boottime <- boottime[bootordinds]
bootdelta <- bootdelta[bootordinds]

bootjumpinds <- c(1:nsamp)[bootdelta==1]

bootna <- rep(0,length(boottime))

  risk <- nsamp
        jna <- ifelse(bootdelta[1]==1,1.0/risk,0.0)
  bootna[1] <- jna

for (i in 2:length(boottime)) {
        risk<-nsamp-i+1
        jna <- ifelse(bootdelta[i]==1,1.0/risk,0.0)
  bootna[i] <- bootna[i-1]+jna

}

ind <- 1

inds <- rep(NA,length(boottime))

newboottime <- rep(NA,length(boottime))

newboottime[1] <- boottime[1]
inds[1] <- 1

for (i in 2:length(boottime)) {

if (boottime[i]!=boottime[i-1]) {

newboottime[ind] <- boottime[i-1] 
inds[ind] <- 1
ind <- ind + 1 

}
}

boottime <- newboottime[!is.na(newboottime)]
bootna <- bootna[inds==1]
bootna <- bootna[!is.na(bootna)]
bootdelta <- bootdelta[inds==1]
bootdelta <- bootdelta[!is.na(bootdelta)]
 
bootjumpinds <- c(1:length(bootna))[bootdelta==1]
boottimevect <- boottime[boottime <= u2 & boottime >= u1]
temp <- alphadblprime(boottimevect,
                     band11,bootjumpinds,boottime,bootna)
bootsum <- 0

times <- sort(unique(c(timevect,boottimevect)))

indstimevect <- rep(0,length(timevect))
indsboottimevect <- rep(0,length(boottimevect))

ind <- 1
indbt <- 1

for (k in 1:length(times)) {

if (!is.na(timevect[ind]) & !is.na(boottimevect[indbt]) & 
length((times[k] - timevect[1:k])[times[k] - timevect[1:k]==0]) > 0 &
length((times[k] - boottimevect[1:k])[times[k] - boottimevect[1:k]==0]) > 0) {
     
indstimevect[ind] <- 1
indsboottimevect[indbt] <- 1
      
ind <- ind + 1
indbt <- indbt + 1
       
}

if (!is.na(timevect[ind]) & !is.na(boottimevect[indbt]) & 
length((times[k] - timevect[1:k])[times[k] - timevect[1:k]==0]) > 0 &
length((times[k] - boottimevect[1:k])[times[k] - boottimevect[1:k]==0])==0) {

indstimevect[ind] <- 0

ind <- ind + 1

}

if (!is.na(timevect[ind]) & !is.na(boottimevect[indbt]) & 
length((times[k] - timevect[1:k])[times[k] - timevect[1:k]==0])==0 &
length((times[k] - boottimevect[1:k])[times[k] - boottimevect[1:k]==0]) > 0) {

indsboottimevect[indbt] <- 0

indbt <- indbt + 1

}


}

newtemporig <- temporig[indstimevect==1]
temp <- temp[indsboottimevect==1]

for (k in 1:length(temp)) {
bootsum <- bootsum + (temp[k]-newtemporig[k])**2
}

MISE <- MISE + bootsum/length(temp)

}       
return(MISE/nboot) 
}

func1 <- optimize(f, c(u1,(u2-u1)/3), nsamp=innsamp, timevect=intimevect, jumpinds=injumpinds, 
time=intime, delta=indelta, na=inna, intK2=inintK2, alpy=inalpyi, u1=inu1, u2=inu2)

##cat(paste("maximum possible band11 = ",(u2-u1)/3),"\n")

#cat(paste("optimal band11 = ",func1$minimum),"\n")
##cat(paste("corresponding minimum MISE = ",func1$objective),"\n")

return(func1$minimum)

}


##
timevect <- time1[time1 <= u2 & time1 >= u1]

band111 <- band11

if (band11==0) {

band111 <- findband11(u1,u2,nsamp1,timevect,jumpinds1,time1,delta1,na1,intK2,alpy1)

##cat(paste("band111 = ",band111),"\n")

}

temp <- alphadblprime(timevect,
                     band111,jumpinds1,time1,na1)**2

intalpdblprime1 <- temp[1]*timevect[1]

for (k in 2:length(temp)) {

intalpdblprime1 <- intalpdblprime1 + temp[k]*
                               (timevect[k] - timevect[k-1])
}

band1opt <- band1

if (band1==0) {

band1opt <- (k2K**(-2/5))*((intK2*alpy1)**(1/5))*
            (intalpdblprime1**(-1/5))*(nsamp1**(-1/5))
}
##

timevect <- time2[time2 <= u2 & time2 >= u1]

band112 <- band12

if (band12==0) {

band112 <- findband11(u1,u2,nsamp2,timevect,jumpinds2,time2,delta2,na2,intK2,alpy2)

##cat(paste("band112 = ",band112),"\n")

}

temp <- alphadblprime(timevect,
                     band112,jumpinds2,time2,na2)**2

intalpdblprime2 <- temp[1]*timevect[1]

for (k in 2:length(temp)) {

intalpdblprime2 <- intalpdblprime2 + temp[k]*
                               (timevect[k] - timevect[k-1])
}

band2opt <- band2

if (band2==0) {

band2opt <- (k2K**(-2/5))*((intK2*alpy2)**(1/5))*
            (intalpdblprime2**(-1/5))*(nsamp2**(-1/5))

}


#cat(paste("Optimal bandwidth placebo group  = ",band1opt),"\n")

#cat(paste("Optimal bandwidth vaccine group = ",band2opt),"\n")

#    Estimate the hazards on the grid of unique failure times in [u1,u2]
#    by Ramlau-Hansen smoothing of the cumulative CSHR's

#    The variances are calculated as displayed in formula (4.2.7) of ABGK

         ee1<-rep(0,ngrid)
         ee2<-rep(0,ngrid) 
         vee1<-rep(0,ngrid)
         vee2<-rep(0,ngrid)
         timegrid <- rep(0,ngrid)

if (jumpinds1[1]==1) {
at.risks1 <- at.risks1[-1] }

if (jumpinds2[1]==1) {
at.risks2 <- at.risks2[-1] }

      for (i in 1:ngrid) {

         tt <- u1+(i-1)*(u2-u1)/(ngrid-1)
         timegrid[i] <- tt

# Rx = 1:
# Use the modified kernel on the interval [u1,t1+band1opt]

         if (tt < t1 + band1opt) {

         if (tailsl==T) {

         ee1[i]<-epanltail((tt-time1[1])/band1opt,(tt-t1)/band1opt)*(na1[1]-0)/band1opt

         at.risk <- length(time1[time1 >= time1[1]])

         vee1[i] <- ifelse(at.risk > 0,
         ((epanltail((tt-time1[1])/band1opt,(tt-t1)/band1opt)/at.risk)**2)*
         (delta1[1]/(band1opt**2)),0)
         
         ee1[i]<-ee1[i]+sum(epanltail((tt-time1[jumpinds1[jumpinds1!=1]])/
                  band1opt,(tt-t1)/band1opt)*
                  (na1[jumpinds1[jumpinds1!=1]]-
                   na1[jumpinds1[jumpinds1!=1]-1]))/band1opt

         vee1[i] <- vee1[i] + sum(ifelse(at.risks1 > 0,1,0)*
         ((epanltail((tt-time1[jumpinds1[jumpinds1!=1]])/band1opt,
         (tt-t1)/band1opt)/at.risks1)**2)*
         (delta1[jumpinds1[jumpinds1!=1]]/(band1opt**2)),0)

# Adjust ee1 for bias

if (biasadjust==T) {

ee1[i] <- ee1[i] - 0.5*(nsamp1**(-2/5))*
         alphadblprime(tt,band111,jumpinds1,time1,na1)*k2Klmod((tt-t1)/band1opt)

}

}

         if (tailsl==F) {

         ee1[i]<-epan((tt-time1[1])/band1opt)*(na1[1]-0)/band1opt

         at.risk <- length(time1[time1 >= time1[1]])

         vee1[i] <- ifelse(at.risk > 0,
         ((epan((tt-time1[1])/band1opt)/at.risk)**2)*(delta1[1]/(band1opt**2)),0)

         ee1[i] <- ee1[i] + sum(epan((tt-time1[jumpinds1[jumpinds1!=1]])
                     /band1opt)*
                     (na1[jumpinds1[jumpinds1!=1]]
                     -na1[jumpinds1[jumpinds1!=1]-1]))/band1opt

         vee1[i] <- vee1[i] + sum(ifelse(at.risks1 > 0,1,0)*
         ((epan((tt-time1[jumpinds1[jumpinds1!=1]])/band1opt)
         /at.risks1)**2)*
         (delta1[jumpinds1[jumpinds1!=1]]/(band1opt**2)),0)

# Adjust ee1 for bias

if (biasadjust==T) {

ee1[i] <- ee1[i] - 0.5*(nsamp1**(-2/5))*
                  alphadblprime(tt,band111,jumpinds1,time1,na1)*k2K

}

}

}

# Use the modified kernel on the interval [u2-band1opt,u2]

         if (tt > t2 - band1opt) {

         if (tailsu==T) {

         ee1[i]<-epanutail((tt-time1[1])/band1opt,-(t2-tt)/band1opt)*(na1[1]-0)/band1opt

         at.risk <- length(time1[time1 >= time1[1]])

         vee1[i] <- ifelse(at.risk > 0,  
         ((epanutail((tt-time1[1])/band1opt,-(t2-tt)/band1opt)/at.risk)**2)*
         (delta1[1]/(band1opt**2)),0)

         ee1[i] <- ee1[i] + sum(epanutail((tt-time1[jumpinds1[jumpinds1!=1]])
                    /band1opt,-(t2-tt)/band1opt)*
                     (na1[jumpinds1[jumpinds1!=1]]
                     -na1[jumpinds1[jumpinds1!=1]-1]))/band1opt

         vee1[i] <- vee1[i] + sum(ifelse(at.risks1 > 0,1,0)*
         ((epanutail((tt-time1[jumpinds1[jumpinds1!=1]])/band1opt,
         -(t2-tt)/band1opt)/at.risks1)**2)*
         (delta1[jumpinds1[jumpinds1!=1]]/(band1opt**2)),0)

# Adjust ee1 for bias

if (biasadjust==T) {

ee1[i] <- ee1[i] - 0.5*(nsamp1**(-2/5))*
         alphadblprime(tt,band111,jumpinds1,time1,na1)*k2Kumod(-(t2-tt)/band1opt)

}

}

         if (tailsu==F) {

         ee1[i]<-epan((tt-time1[1])/band1opt)*(na1[1]-0)/band1opt

         at.risk <- length(time1[time1 >= time1[1]])

         vee1[i] <- ifelse(at.risk > 0,
         ((epan((tt-time1[1])/band1opt)/at.risk)**2)*(delta1[1]/(band1opt**2)),0)

         ee1[i] <- ee1[i] + sum(epan((tt-time1[jumpinds1[jumpinds1!=1]])
                     /band1opt)*
                     (na1[jumpinds1[jumpinds1!=1]]
                     -na1[jumpinds1[jumpinds1!=1]-1]))/band1opt

         vee1[i] <- vee1[i] + sum(ifelse(at.risks1 > 0,1,0)*
         ((epan((tt-time1[jumpinds1[jumpinds1!=1]])/band1opt)
         /at.risks1)**2)*
         (delta1[jumpinds1[jumpinds1!=1]]/(band1opt**2)),0)

# Adjust ee1 for bias

if (biasadjust==T) {

ee1[i] <- ee1[i] - 0.5*(nsamp1**(-2/5))*
                  alphadblprime(tt,band111,jumpinds1,time1,na1)*k2K

}

}

}

if (tt >= t1 + band1opt & tt <= t2 - band1opt) {

         ee1[i]<-epan((tt-time1[1])/band1opt)*(na1[1]-0)/band1opt

         at.risk <- length(time1[time1 >= time1[1]])

         vee1[i] <- ifelse(at.risk > 0,
         ((epan((tt-time1[1])/band1opt)/at.risk)**2)*(delta1[1]/(band1opt**2)),0)

         ee1[i] <- ee1[i] + sum(epan((tt-time1[jumpinds1[jumpinds1!=1]])
                     /band1opt)*
                     (na1[jumpinds1[jumpinds1!=1]]
                     -na1[jumpinds1[jumpinds1!=1]-1]))/band1opt

         vee1[i] <- vee1[i] + sum(ifelse(at.risks1 > 0,1,0)*
         ((epan((tt-time1[jumpinds1[jumpinds1!=1]])/band1opt)
         /at.risks1)**2)*
         (delta1[jumpinds1[jumpinds1!=1]]/(band1opt**2)),0)


# Adjust ee1 for bias

if (biasadjust==T) {

ee1[i] <- ee1[i] - 0.5*(nsamp1**(-2/5))*
                  alphadblprime(tt,band111,jumpinds1,time1,na1)*k2K

}

} 



# Rx = 2:
# Use the modified kernel on the interval [u1,t1+band2opt]

         if (tt < t1 + band2opt) {

         if (tailsl==T) {

         ee2[i]<-epanltail((tt-time2[1])/band2opt,(tt-t1)/band2opt)*(na2[1]-0)/band2opt

         at.risk <- length(time2[time2 >= time2[1]])

         vee2[i] <- ifelse(at.risk > 0,
         ((epanltail((tt-time2[1])/band2opt,(tt-t1)/band2opt)/at.risk)**2)*
         (delta2[1]/(band2opt**2)),0)

         ee2[i] <- ee2[i] + sum(epanltail((tt-time2[jumpinds2[jumpinds2!=1]])
                    /band2opt,(tt-t1)/band2opt)*
                     (na2[jumpinds2[jumpinds2!=1]]
                     -na2[jumpinds2[jumpinds2!=1]-1]))/band2opt

         vee2[i] <- vee2[i] + sum(ifelse(at.risks2 > 0,1,0)*
         ((epanltail((tt-time2[jumpinds2[jumpinds2!=1]])/band2opt,
         (tt-t1)/band2opt)/at.risks2)**2)*
         (delta2[jumpinds2[jumpinds2!=1]]/(band2opt**2)),0)

# Adjust ee2 for bias

if (biasadjust==T) {

ee2[i] <- ee2[i] - 0.5*(nsamp2**(-2/5))*
         alphadblprime(tt,band112,jumpinds2,time2,na2)*k2Klmod((tt-t1)/band2opt)

}

}

         if (tailsl==F) {

         ee2[i]<-epan((tt-time2[1])/band2opt)*(na2[1]-0)/band2opt

         at.risk <- length(time2[time2 >= time2[1]])

         vee2[i] <- ifelse(at.risk > 0,
         ((epan((tt-time2[1])/band2opt)/at.risk)**2)*(delta2[1]/(band2opt**2)),0)

         ee2[i] <- ee2[i] + sum(epan((tt-time2[jumpinds2[jumpinds2!=1]])
                     /band2opt)*
                     (na2[jumpinds2[jumpinds2!=1]]
                     -na2[jumpinds2[jumpinds2!=1]-1]))/band2opt

         vee2[i] <- vee2[i] + sum(ifelse(at.risks2 > 0,1,0)*
         ((epan((tt-time2[jumpinds2[jumpinds2!=1]])/band2opt)
         /at.risks2)**2)*
         (delta2[jumpinds2[jumpinds2!=1]]/(band2opt**2)),0)

# Adjust ee2 for bias

if (biasadjust==T) {

ee2[i] <- ee2[i] - 0.5*(nsamp2**(-2/5))*
                  alphadblprime(tt,band112,jumpinds2,time2,na2)*k2K

}

}

}

# Use the modified kernel on the interval [u2-band2opt,u2]

         if (tt > t2 - band2opt) {

         if (tailsu==T) {

         ee2[i]<-epanutail((tt-time2[1])/band2opt,-(t2-tt)/band2opt)*(na2[1]-0)/band2opt

         at.risk <- length(time2[time2 >= time2[1]])

         vee2[i] <- ifelse(at.risk > 0,  
         ((epanutail((tt-time2[1])/band2opt,-(t2-tt)/band2opt)/at.risk)**2)*
         (delta2[1]/(band2opt**2)),0)

         ee2[i] <- ee2[i] + sum(epanutail((tt-time2[jumpinds2[jumpinds2!=1]])
                    /band2opt,-(t2-tt)/band2opt)*
                     (na2[jumpinds2[jumpinds2!=1]]
                     -na2[jumpinds2[jumpinds2!=1]-1]))/band2opt

         vee2[i] <- vee2[i] + sum(ifelse(at.risks2 > 0,1,0)*
         ((epanutail((tt-time2[jumpinds2[jumpinds2!=1]])/band2opt,
         -(t2-tt)/band2opt)/at.risks2)**2)*
         (delta2[jumpinds2[jumpinds2!=1]]/(band2opt**2)),0)

# Adjust ee2 for bias

if (biasadjust==T) {

ee2[i] <- ee2[i] - 0.5*(nsamp2**(-2/5))*
         alphadblprime(tt,band112,jumpinds2,time2,na2)*k2Kumod(-(t2-tt)/band2opt)

}

}

         if (tailsu==F) {

         ee2[i]<-epan((tt-time2[1])/band2opt)*(na2[1]-0)/band2opt

         at.risk <- length(time2[time2 >= time2[1]])

         vee2[i] <- ifelse(at.risk > 0,
         ((epan((tt-time2[1])/band2opt)/at.risk)**2)*(delta2[1]/(band2opt**2)),0)

         ee2[i] <- ee2[i] + sum(epan((tt-time2[jumpinds2[jumpinds2!=1]])
                     /band2opt)*
                     (na2[jumpinds2[jumpinds2!=1]]
                     -na2[jumpinds2[jumpinds2!=1]-1]))/band2opt

         vee2[i] <- vee2[i] + sum(ifelse(at.risks2 > 0,1,0)*
         ((epan((tt-time2[jumpinds2[jumpinds2!=1]])/band2opt)
         /at.risks2)**2)*
         (delta2[jumpinds2[jumpinds2!=1]]/(band2opt**2)),0)

# Adjust ee2 for bias

if (biasadjust==T) {

ee2[i] <- ee2[i] - 0.5*(nsamp2**(-2/5))*
                  alphadblprime(tt,band112,jumpinds2,time2,na2)*k2K

}

}

}

if (tt >= t1 + band2opt & tt <= t2 - band2opt) {

         ee2[i]<-epan((tt-time2[1])/band2opt)*(na2[1]-0)/band2opt

         at.risk <- length(time2[time2 >= time2[1]])

         vee2[i] <- ifelse(at.risk > 0,
         ((epan((tt-time2[1])/band2opt)/at.risk)**2)*(delta2[1]/(band2opt**2)),0)

         ee2[i] <- ee2[i] + sum(epan((tt-time2[jumpinds2[jumpinds2!=1]])
                     /band2opt)*
                     (na2[jumpinds2[jumpinds2!=1]]
                     -na2[jumpinds2[jumpinds2!=1]-1]))/band2opt

         vee2[i] <- vee2[i] + sum(ifelse(at.risks2 > 0,1,0)*
         ((epan((tt-time2[jumpinds2[jumpinds2!=1]])/band2opt)
         /at.risks2)**2)*
         (delta2[jumpinds2[jumpinds2!=1]]/(band2opt**2)),0)


# Adjust ee2 for bias

if (biasadjust==T) {

ee2[i] <- ee2[i] - 0.5*(nsamp2**(-2/5))*
                  alphadblprime(tt,band112,jumpinds2,time2,na2)*k2K

}

} 

}

# Modify the time grid to only include time points at which vee1 and 
# vee2 are positive

timesunique <- timegrid[vee1 > 0 & vee2 > 0] 
indsunique <- c(1:length(timegrid))[vee1 > 0 & vee2 > 0]
#browser()
ee1 <- ee1[indsunique]
vee1 <- vee1[indsunique]
ee2 <- ee2[indsunique]
vee2 <- vee2[indsunique]

lenunique <- length(indsunique)

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

# Calculate the function v over the grid of time points

vval <- rep(0,lenunique)

if (Nv==0)  {
vval <- v(nsamp1,nsamp2,ee1,ee2,vee1,vee2)
}

if (Nv > 0) {  

U1simulmat <- matrix(rep(NA,Nv*lenunique),nrow=lenunique)
U2simulmat <- matrix(rep(NA,Nv*lenunique),nrow=lenunique)

for (i in 1:lenunique) {

tt <- timesunique[i]

U1simulmat[i,] <- U(Nv,tt,jumpinds1,nsamp,time1,delta1,na1,band1opt,band111,u1,
   u2,t1,t2,tailsl,tailsu,biasadjust)

U2simulmat[i,] <- U(Nv,tt,jumpinds2,nsamp,time2,delta2,na2,band2opt,band112,u1,
   u2,t1,t2,tailsl,tailsu,biasadjust)


  temp <- (var(U1simulmat[i,])*g1(ee1[i],ee2[i])**2) +
         (var(U2simulmat[i,])*g2(ee1[i],ee2[i])**2)

  vval[i] <- temp**(-1/2)

#paste(#cat("g1(ee1[i],ee2[i])**2=",g1(ee1[i],ee2[i])**2),"\n")
#paste(#cat("g2(ee1[i],ee2[i])**2=",g2(ee1[i],ee2[i])**2),"\n")
#paste(#cat("var(U1simulmat[i,])=",var(U1simulmat[i,])),"\n")
#paste(#cat("var(U2simulmat[i,])=",var(U2simulmat[i,])),"\n")
#paste(#cat("vval[i]=",vval[i]),"\n")

}

}


critvalband <- Critvalue(alpha,N,jumpinds1,jumpinds2,nsamp,timesunique,time1,
time2,delta1,delta2,na1,na2,ee1,ee2,vee1,vee2,vval,lenunique,band1opt,band2opt,band111,band112,
u1,u2,t1,t2,tailsl,tailsu,biasadjust)

critvalband95 <- critvalband[1]

critvalband90 <- critvalband[2]

critvalband80 <- critvalband[3]

#cat(paste("Critical value for 95% bands = ",critvalband95),"\n")

#cat(paste("Critical value for 90% bands = ",critvalband90),"\n")

#cat(paste("Critical value for 80% bands = ",critvalband80),"\n")

  lowint <- g(ee1,ee2) - ((nsamp**(-2/5))*critvalint)/vval
  upint <- g(ee1,ee2) + ((nsamp**(-2/5))*critvalint)/vval
  lowband95 <- g(ee1,ee2) - ((nsamp**(-2/5))*critvalband95)/vval
  upband95 <- g(ee1,ee2) + ((nsamp**(-2/5))*critvalband95)/vval
  lowband90 <- g(ee1,ee2) - ((nsamp**(-2/5))*critvalband90)/vval
  upband90 <- g(ee1,ee2) + ((nsamp**(-2/5))*critvalband90)/vval
  lowband80 <- g(ee1,ee2) - ((nsamp**(-2/5))*critvalband80)/vval
  upband80 <- g(ee1,ee2) + ((nsamp**(-2/5))*critvalband80)/vval


###############
###############

answermat <- cbind(timesunique,log(ee2/ee1),lowint,upint,lowband95,upband95,
                   lowband90,upband90,lowband80,upband80,ee1,ee2,vee1,vee2)

return(answermat)

}


##############################################################################
# The notation of the following functions used in the main program is as in 
# Gilbert et al. (2002, Biometrics)

# Sub-functions:
######################  

######################
# the Epanechnikov function

epan <- function(x) {

return(ifelse(abs(x) <= 1,0.75*(1-x**2),0))

}

#######################

######################
# the derivative of the Epanechnikov function

epander <- function(x) {

return(ifelse(abs(x) <= 1,-1.50*x,0))

}

#######################

######################
# the Modified Epanechnikov function, for the lower tail

epanltail <- function(x,q) {

gammaq <- (3/4)*(((2/15)+(q^3/3)-(q^5/5))*((2/3)+q-(q^3/3))
         - (1/16)*((1 - q^2)**4))

gammaq <- 1/gammaq

alphaq <- ((2/15)+(q^3/3)-(q^5/5))*gammaq

betaq <- ((1 - q^2)**2)*(gammaq/4)

return(ifelse(x <= q & x >= -1,epan(x)*(alphaq + betaq*x),0))

}

#######################

######################
# the derivative of the Modified Epanechnikov function, for the lower tail

epanltailder <- function(x,q) {

gammaq <- (3/4)*(((2/15)+(q^3/3)-(q^5/5))*((2/3)+q-(q^3/3))
         - (1/16)*((1 - q^2)**4))

gammaq <- 1/gammaq

alphaq <- ((2/15)+(q^3/3)-(q^5/5))*gammaq

betaq <- ((1 - q^2)**2)*(gammaq/4)

return(ifelse(x <= q & x >= -1,epander(x)*(alphaq + betaq*x) + betaq*epan(x),0))

}

#######################


######################
# the Modified Epanechnikov function, for the upper tail

epanutail <- function(x,q) {

tempq<-2/15-(q^3/3)+(q^5/5)

gammaq <-1/((3/4)*(tempq*((4/5)-q+(q^3/3)) - ((1 - q^2)^4)/16))

alphaq <- tempq*gammaq

betaq <- -((1 - q^2)^2)*gammaq/4

return(as.numeric(x >=q & x <=1)*epan(x)*(alphaq + betaq*x))
}


#######################

######################
# the derivative of the Modified Epanechnikov function, for the upper tail

epanutailder <- function(x,q) {

gammaq <- (3/4)*(((-4/5)+q-(q^3/3))*((2/15)-(q^3/3)+(q^5/5))/
         ((1/4)-(q^2/2)+(q^4/4)) + ((3/4)-(q^2/2)-(q^4/4)))

betaq <- ifelse(((1/4)-(q^2/2)+(q^4/4))==0,0,(4/3)*(((-4/5)+q-(q^3/3))*
((2/15)-(q^3/3)+(q^5/5))/((1/4)-(q^2/2)+(q^4/4)) + ((3/4)-(q^2/2)-(q^4/4)))**(-1))

alphaq <- ((4/3) - betaq*((3/4)-(q^2/2)-(q^4/4)))/((4/5)-q+(q^3/3))

return(ifelse(x <= 1 & x >= q,epander(x)*(alphaq + betaq*x) + betaq*epan(x),0))

}

#######################



#######################
# function g

g <- function(x,y) {

return(log(y) - log(x))

}

####################### 

#######################
# function g1

g1 <- function(x,y) {

return(-1/x)

}

#######################

#######################
# function g2

g2 <- function(x,y) {

return(1/y)

}

#######################


#######################
# function v

v <- function(n1,n2,ee1,ee2,vee1,vee2) {

return((((n1+n2)**(4/5))*(vee1/(ee1**2) +
                         vee2/(ee2**2)))**(-1/2))
}

#######################

#######################
# function U

U <- function(N,t,inds,nsamp,time,delta,na,band,band11,u1,u2,t1,t2,tailsl,tailsu,
biasadjust) {

k2K <- 0.20

# Calculate the constant k2Klmod and k2Kumod, used for bias
# adjustment in the tails

k2Klmod <- function(q) {

return(0.020202*sum((seq(-1,1,length=100)**2)*epanltail(seq(-1,1,length=100),q)))

}

k2Kumod <- function(q) {

return(0.020202*sum((seq(-1,1,length=100)**2)*epanutail(seq(-1,1,length=100),q)))

}

sum <- rep(0,N)

for (j in inds) {

x <- length(time[time >= time[j]])

rec.atrisk <- ifelse(x > 0,1/x,0)

if (t >= u1 & t < t1 + band & tailsl==T) {

sum <- sum + (ifelse(time[j] >= t1 & time[j] <= t2,1,0)*rec.atrisk*
epanltail((t - time[j])/band,(t-t1)/band)*
delta[j]*rnorm(N))

}

if (t <= u2 & t > t2 - band & tailsu==T) {

sum <- sum + (ifelse(time[j] >= t1 & time[j] <= t2,1,0)*rec.atrisk*
epanutail((t - time[j])/band,-(t2-t)/band)*
delta[j]*rnorm(N))

}

else {

sum <- sum + 
(ifelse(time[j] >= t1 & time[j] <= t2,1,0)*rec.atrisk*
epan((t - time[j])/band)*
delta[j]*rnorm(N))

}
}

if (biasadjust==T) {

if (t < t1 + band & tailsl==T) {

sum <- sum + 0.5*nsamp**(-2/5)*alphadblprime(t,band11,inds,time,na)*
            k2Klmod((t-t1)/band) 

}

if (t <= u2 & t > t2 - band & tailsu==T) {

sum <- sum + 0.5*nsamp**(-2/5)*alphadblprime(t,band11,inds,time,na)*
            k2Kumod(-(t2-t)/band)

}

else {

sum <- sum + 0.5*nsamp**(-2/5)*alphadblprime(t,band11,inds,time,na)*
            k2K
}
}

return(((nsamp**(2/5))/band)*sum)

}

#######################


#######################
# function Vtilde

Vtilde <- function(N,i,t,inds1,inds2,nsamp,time1,time2,delta1,delta2,na1,na2,
         ee1,ee2,vee1,vee2,vval,band1,band2,band111,band112,u1,u2,t1,t2,
         tailsl,tailsu,biasadjust)  {

return(vval[i]*
       (g2(ee1[i],ee2[i])*U(N,t,inds2,nsamp,time2,delta2,na2,band2,band111,u1,
u2,t1,t2,tailsl,tailsu,biasadjust)
      + g1(ee1[i],ee2[i])*U(N,t,inds1,nsamp,time1,delta1,na1,band1,band112,u1,
u2,t1,t2,tailsl,tailsu,biasadjust)))

}



#######################

# The simulated test process L*(t)

#######################

########################
# function Gtilde

Gtilde <- function(N,inds1,inds2,nsamp,timesunique,time1,time2,delta1,delta2,
na1,na2,ee1,ee2,vee1,vee2,vval,lenunique,band1,band2,band111,band112,u1,u2,t1,t2,
tailsl,tailsu,biasadjust)  {

max <- rep(0,N)

for (i in 1:lenunique) {

tt<-timesunique[i]

x <- abs(Vtilde(N,i,tt,inds1,inds2,nsamp,time1,time2,delta1,delta2,na1,na2,
   ee1,ee2,vee1,vee2,vval,band1,band2,band111,band112,u1,u2,t1,t2,tailsl,tailsu,biasadjust))

max <- ifelse(x > max,x,max)

}

return(max)

}

#########################

#########################
# function Critvalue

Critvalue <- function(alpha,N,inds1,inds2,nsamp,timesunique,time1,time2,
delta1,delta2,na1,na2,ee1,ee2,vee1,vee2,vval,lenunique,band1,band2,band111,band112,
u1,u2,t1,t2,tailsl,tailsu,biasadjust)  {

Gtildevect <- Gtilde(N,inds1,inds2,nsamp,timesunique,time1,time2,delta1,delta2,
na1,na2,ee1,ee2,vee1,vee2,vval,lenunique,band1,band2,band111,band112,u1,u2,t1,t2,tailsl,tailsu,
biasadjust)

return(sort(Gtildevect)[floor((1-alpha)*N)])

}

#########################


###################
# function K1dblprime

K1dblprime <- function(x) {

y <- (-105/16)*(1-6*x**2+5*x**4)

return(ifelse(abs(x) <= 1,y,0))

}

###################

###################
# function alphadblprime 

alphadblprime <- function(t,band,inds,time,na) {

ans <- 0

ans <- K1dblprime((t - time[1])/band)*(na[1] - 0)

for (j in inds[inds!=1]) {

ans <- ans + K1dblprime((t - time[j])/band)*(na[j] - na[j-1])

}

return(ans/(band**3))

}

