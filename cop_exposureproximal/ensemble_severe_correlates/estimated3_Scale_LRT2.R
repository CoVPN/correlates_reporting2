
ind.event<-as.integer(commandArgs(trailingOnly=T)[1])
ind.marker<-as.integer(commandArgs(trailingOnly=T)[2])
ind.region<-as.integer(commandArgs(trailingOnly=T)[3])

# ind.event=1; ind.marker=1; ind.region=1

library(WeMix)
library(knitr)
library(readr)
library(tidyr)
library(extraDistr)
library(dplyr)
library(lubridate)
library(purrr)
#library(patchwork)
library(survival)
library(survey)
#library(marginalizedRisk)

library(foreach)
library(doRNG)
library(lme4)

oc<-3

excludeOutlierC<-function(x,coef){
 Q3<-quantile(x,.75,na.rm=T)
 Q1<-quantile(x,.25,na.rm=T)
 IQR.x<-IQR(x,na.rm=T)
 x.out<-ifelse(x>Q3+coef*IQR.x,NA,x)
 x.out<-ifelse(x<Q1-coef*IQR.x,NA,x.out)
 return(x.out)
}

# setwd("~/COVPN_P3003")

source(file="../FunctionCall.R")
library(lme4)

#setwd("C:/Users/yhuang/OneDrive - Fred Hutchinson Cancer Research Center/Documents/All_Files/1yingsstuff/COVPN_P3003")
#source(file="../Pepe/RA/Functions/FunctionCall.R")


vv.endi<-"EventIndPrimaryIncludeNotMolecConfirmedD29_NoRegionCens" 
vv.endt<-"EventTimePrimaryIncludeNotMolecConfirmedD29_NoRegionCens"

# library(config)# don't call this because then merge becomes ambiguous
config.reporting <- config::get(config = "janssen_pooled_partA", file="../../config.yml") 
dat<-read.csv(config.reporting$data_cleaned,na.strings=c("n/a","NA","N/A","","."))
# dat<-read.csv("adata/janssen_pooled_partA_data_processed_with_riskscore_20240305.csv",na.strings=c("n/a","NA","N/A","","."))

vv.delta<-c("EventIndPrimaryIncludeNotMolecConfirmedD29","SevereEventIndPrimaryIncludeNotMolecConfirmedD29",
 "ModerateEventIndPrimaryIncludeNotMolecConfirmedD29")

vv.time<-c("EventTimePrimaryIncludeNotMolecConfirmedD29","SevereEventTimePrimaryIncludeNotMolecConfirmedD29",
"ModerateEventTimePrimaryIncludeNotMolecConfirmedD29")
 
vv.D29<-c("Day29bindSpike","Day29bindRBD","Day29pseudoneutid50")
vv.D71<-c("Day71bindSpike","Day71bindRBD","Day71pseudoneutid50")
vv.M6<-c("Mon6bindSpike","Mon6bindRBD","Mon6pseudoneutid50")


oo.wdate<-dat$NumberdaysD1toD29>=dat$NumberdaysD1toD71 | dat$NumberdaysD1toD29>=dat$NumberdaysD1toM6 |
dat$NumberdaysD1toD71>=dat$NumberdaysD1toM6

oo.wdate<-(1:length(oo.wdate))[oo.wdate]
oo.wdate<-na.omit(oo.wdate)
length(oo.wdate)
#[1] 0

dat[oo.wdate,c('NumberdaysD1toD29','NumberdaysD1toD71','NumberdaysD1toM6')]

oo.immune<-dat$ph1.D29=='TRUE' & (dat$Trt==0 | dat$ph2.D29=='TRUE') & dat$Bserostatus==0



#summary(dat$CalendarDateEnrollment)

#dat$D29toD71<-dat$NumberdaysD1toD71-dat$NumberdaysD1toD29
#dat$D29toM6<-dat$NumberdaysD1toM6-dat$NumberdaysD1toD29
#
#summary(dat$D29toD71)
#summary(dat$D29toM6)
#
#dat$D29toD71<-excludeOutlierC(dat$D29toD71,3)
#dat$D29toM6<-excludeOutlierC(dat$D29toM6,3)
#summary(dat$D29toD71)
#summary(dat$D29toM6)
##
#dat$NumberdaysD1toD71<-ifelse(is.na(dat$D29toD71),NA,dat$NumberdaysD1toD71)
#dat$NumberdaysD1toM6<-ifelse(is.na(dat$D29toM6),NA,dat$NumberdaysD1toM6)
##

####

data<-dat[oo.immune,]
X0.min<-min(data[,vv.D29[ind.marker]],na.rm=T)
#> table(!is.na(data$NumberdaysD1toD29),!is.na(data$NumberdaysD1toD71))
#
#       FALSE  TRUE
#  TRUE   409 18615
#
#> table(!is.na(data$NumberdaysD1toD29),!is.na(data$NumberdaysD1toM6))
#
#       FALSE  TRUE
#  TRUE  1373 17651

#############



####

by(data[,vv.D29[1]],data$Trt,summary)
#data$Trt: 0
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.734   0.734   0.734   0.736   0.734   1.226   17620 
#--------------------------------------------------------------------------------------------------------------------- 
#data$Trt: 1
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.7341  1.2088  1.5330  1.5158  1.8448  3.6505 

by(data[,vv.D29[2]],data$Trt,summary)
#data$Trt: 0
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.848   0.848   0.848   0.856   0.848   1.411   17620 
#--------------------------------------------------------------------------------------------------------------------- 
#data$Trt: 1
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.8478  1.2127  1.5061  1.5106  1.8080  3.8632 


by(data[,vv.D29[3]],data$Trt,summary)
#data$Trt: 0
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.389   0.389   0.389   0.393   0.389   1.337   17610 
#--------------------------------------------------------------------------------------------------------------------- 
#data$Trt: 
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.3889  0.3889  0.3889  0.7489  1.0324  3.5170 

###
data$Sex<-ifelse(data$Sex==2,0.5,data$Sex)


data$weight<-ifelse(data$Trt==1,data$wt.D29,1)

data$ID<-1:nrow(data)


###########3

oo.Z1<-(1:nrow(data))[data$Trt==1]


#############
dd<-data[oo.Z1,]
ddd<-dd[!is.na(dd[,vv.M6[3]]),]
check=cbind(ddd[,'NumberdaysD1toM6']-ddd$NumberdaysD1toD29,ddd[,vv.time[1]],ddd[,vv.endt])
#> sum(check[,1]<=check[,3],na.rm=T)
#[1] 78
#> sum(check[,1]<=check[,3],na.rm=T)
#[1] 78
tt<-ww<-matrix(NA,nrow(data),3)

datlong<-NULL
for (i in oo.Z1){
  
   T.visit<-unlist(data[i,c('NumberdaysD1toD29','NumberdaysD1toD71','NumberdaysD1toM6')]-data$NumberdaysD1toD29[i])
   XE.long<-unlist(data[i,c(vv.D29[ind.marker],vv.D71[ind.marker],vv.M6[ind.marker])])
   #oo.visit<-(1:3)[T.visit<TE[i] & !is.na(T.visit) & !is.na(XE.long)]
   oo.visit<-(1:3)[T.visit<data[i,vv.endt] & !is.na(T.visit) & !is.na(XE.long)]
   #if (delta[i]==0) oo.visit<-(1:3)[!is.na(T.visit) & !is.na(XE.long)]
   TE.truc<-T.visit[oo.visit]
   XE.truc<-XE.long[oo.visit]
    #TE.truc<-c(TE.truc,TE[i])
    #XE.truc<-c(XE.truc,XE[i])
    RE.temp<-data$CalendarDateEnrollment[i]+data$NumberdaysD1toD29[i]
    TE.temp<-data[i,vv.time[ind.event]]
    SE.temp<-TE.temp+RE.temp
    SE.truc<-TE.truc+RE.temp
    delta.temp<-data[i,vv.delta[ind.event]]
    ll=length(TE.truc)
  #  W1.truc<-rep(WE1[i],ll)
  #  W2.truc<-rep(WE2[i],ll)
  #  W3.truc<-rep(WE3[i],ll)
    datlong<-rbind(datlong,cbind(rep(data[i,]$ID,ll),rep(RE.temp,ll),TE.truc,XE.truc,SE.truc,rep(data$Trt[i],ll),rep(SE.temp,ll),rep(0,ll),1:ll,rep(delta.temp,ll),rep(data$weight[i],ll),
    rep(data$Age[i],ll),rep(data$BMI[i],ll),rep(data$Sex[i],ll),rep(data$Region[i],ll)))
    
    ########3
    
    tt.temp<-ww.temp<-rep(NA,3)
    #ifelse(T.visit<TE[i],T.visit,NA)
    ifelse(T.visit<data[i,vv.endt],T.visit,NA)
    tt.temp[oo.visit]<-TE.truc
    ww.temp[oo.visit]<-XE.truc
    tt[i,]<-tt.temp
    ww[i,]<-ww.temp
    
    rm(RE.temp,TE.temp,SE.temp,delta.temp)
    
}

tt[data$Trt==0,]<-c(0,NA,NA)
ww[data$Trt==0,]<-c(0,NA,NA)


  
datlong<-as.data.frame(datlong)

### delta is censoring indicator ##
colnames(datlong)<-c('id','R','time','X','S','Z','eventime','delta','visit','event','weight','Age','BMI','Sex','Region')
datlong<-datlong[datlong$Z==1 & !is.na(datlong$X),]

datlong$X<-datlong$X-X0.min
################



##################

X.long<-datlong$X
T.long<-datlong$time
id.long<-datlong$id

###############

qq=tapply(rep(1,nrow(datlong)),datlong$id,sum)
qq.id<-as.numeric(names(qq)[qq>1])

datlong.2<-datlong[datlong$id%in%qq.id,]
####


#datlong$X.star=datlong$X+rnorm(nrow(datlong),0,0.2)
### Here we just use random subcohort to estimate LME
#fit<-lmer(X~time +factor(Region)+Age+Sex+(1|id), data=datlong.2[datlong.2$Z==1 & datlong.2$delta==0,])
#

#fit<-lmer(X~time+Age+Sex+(1|id), data=datlong.2[datlong.2$Z==1 & datlong.2$delta==0,])
#coef.fix<-getME(fit,name=c("fixef"))
#coef.rd<-as.matrix(VarCorr(fit)[[1]])
#sigma.e<-getME(fit,name=c("sigma"))

#coef.fix<-getME(fit,name=c("fixef"))
#coef.rd<-as.matrix(VarCorr(fit)[[1]])
#coef.rd<-c(coef.rd,0,0,0)
#coef.rd<-matrix(coef.rd,byrow=T,nrow=2)
#sigma.e<-getME(fit,name=c("sigma"))

fac0=length(unique(datlong[datlong$event==0,]$id))/length(unique(datlong.2[datlong.2$event==0,]$id))
fac1=length(unique(datlong[datlong$event==1,]$id))/length(unique(datlong.2[datlong.2$event==1,]$id))

datlong.2$weightScale<-ifelse(datlong.2$event==0,datlong.2$weight*fac0,datlong.2$weight*fac1)
#


### here account for sampling weight

dd<-datlong.2[datlong.2$Z==1,]
dd$W1<-1
#fit<-mix(X~time + Age+Sex + (1|id), data=dd,weights=c('W1','weight'),cWeights=TRUE)
fit<-mix(X~time + Age+Sex + (1|id), data=dd,weights=c('W1','weightScale'),cWeights=TRUE)

#coef.fix<-fit$coef
#coef.rd<-unlist(fit$varVC[2])
#coef.rd<-matrix(coef.rd,byrow=2,nrow=2)
#sigma.e<-fit$sigma
#
coef.fix<-fit$coef
coef.rd<-unlist(fit$varVC[2])
coef.rd<-c(coef.rd,0,0,0)
coef.rd<-matrix(coef.rd,byrow=2,nrow=2)
sigma.e<-fit$sigma

#### create a wide-form data with one observation each participant that include multiple time point with immune response measure

T.visit<-c(1,2)

datwide.delta0 <- reshape(datlong[c('id','visit','X')], idvar = "id", 
timevar = "visit",v.names="X",direction = "wide")
datwide.delta0$groupall<-apply(!is.na(datwide.delta0[,-1]),1,sum,na.rm=T)

Vname<-paste('X.',T.visit,sep='')
datwide.delta0$group<-apply(!is.na(datwide.delta0[,Vname]),1,sum,na.rm=T)


datwide.delta0.time <- reshape(datlong[c('id','visit','time')], idvar = "id", 
timevar = "visit",v.names="time",direction = "wide")

datwide.delta0<-merge(datwide.delta0,datwide.delta0.time,by="id",all=T)
Tname<-paste('time.',T.visit,sep='')


#########
#### if Region specific, only take the data from the specific region

if (ind.region<=2) data<-data[data$Region==ind.region,]

RE<-data$CalendarDateEnrollment+data$NumberdaysD1toD29
TE<-data[,vv.time[ind.event]]
SE<-RE+TE
ZE=data$Trt
XE0=ifelse(data$Trt==1, data[,vv.D29[ind.marker]]-X0.min,NA)
XE0=ifelse(data$Trt==0,0,XE0)
XE=XE0
WE=cbind(data$Age,data$Sex)
WEc=cbind(data$risk_score,data$Region==1,data$Region==2)
#WE1=data$risk_score
#WE2=data$Region==1
#WE3=data$Region==2
#WE3=corr_dat$MinorityInd
#WE4=corr_dat$HighRiskInd
#WE5=corr_dat$risk_score
delta=data[,vv.delta[ind.event]]
weight=data$weight

ID<-data$ID
#incoh<-corr_dat$Trt==1 & !is.na(corr_dat$Day57pseudoneutid50)
incoh<-data$Trt==1 & data$ph2.D29 & !is.na(data[,vv.D29[ind.marker]])












 

###########

SEE<-ifelse(delta==1, SE, SE+max(SE))
#oo<-order(SE)
oo<-order(SEE)

SE<-SE[oo]
TE<-TE[oo]
XE<-XE[oo]
XE0<-XE0[oo]
WE=WE[oo,]
WEc<-WEc[oo,]
#WE1<-WE1[oo]
#WE2<-WE2[oo]
#WE3<-WE3[oo]
#WE4<-WE4[oo]
#WE5<-WE5[oo]
ZE<-ZE[oo]
RE<-RE[oo]
delta=delta[oo]
ll=sum(delta)
incoh=incoh[oo]
weight=weight[oo]
ID=ID[oo]


#### 

oo<-ZE==0 | (ZE==1 & incoh==1)
SE<-SE[oo]
TE<-TE[oo]
XE<-XE[oo]
XE0<-XE0[oo]
WE<-WE[oo,]
WEc<-WEc[oo,]
ZE<-ZE[oo]
RE<-RE[oo]
#mu.EXE<-mu.EXE[oo]
#sigma.EXE<-sigma.EXE[oo]
delta=delta[oo]
ll=sum(delta)
incoh=incoh[oo]
ID=ID[oo]
weight=weight[oo]


## unique value for SE[delta==1]

SE.u<-unique(SE[delta==1])

ll.u<-length(SE.u)
#### Get corresponding data from vaccinees only to enter the cox model based on only vaccine recipients ####################

oo.v=ZE==1
SE.v=SE[oo.v]
TE.v=TE[oo.v]
XE.v=XE[oo.v]
XE0.v=XE0[oo.v]
WE.v=WE[oo.v,]
WEc.v=WEc[oo.v,]
#WE1.v=WE1[oo.v]
#WE2.v=WE2[oo.v]
#WE3.v=WE3[oo.v]
#WE4.v=WE4[oo.v]
#WE5.v=WE5[oo.v]
ZE.v=ZE[oo.v]
RE.v=RE[oo.v]
delta.v=delta[oo.v]
ll.v=sum(delta.v)
incoh.v=incoh[oo.v]
ID.v=ID[oo.v]
weight.v=weight[oo.v]


### unique value for SE.v[delta.v==1]

SE.uv<-sort(unique(SE.v[delta.v==1]))

  
### use those with XE(0)
### 

weights=weight
################
  dw<-2 ### dimension of W
  if (ind.region==3){
  dc<-3
  } else {
  
  WEc=as.matrix(WEc[,1])
  WEc.v=as.matrix(WEc.v[,1])
  dc<-1
  }
  
######

load(file=paste("output/ES_event",ind.event,"_marker",ind.marker,"_region",ind.region,".Rdata",sep=''))

  
  low=-10;up=10
  
  if (ind.region==3){
  source(file="MethodFund_forEnsemble3TieFull.R")
#  B1E=-.8;B2E=-.8;B3E=0
#  parm1=B1E;parm2=B2E;parm3=B3E
 
  
  ### fit a cox model with coefficients with x(0),x(T), and T
#  fit123<-optim(c(B1E,B2E,B3E),nloglik.M1.a.Cal)
  
  ### fit a cox model with coefficient with x(0) only 
#   B00=fit1.vp$coef[1]; B1E=fit1.vp$coef[2];parm2=0;parm3=0; B4E=fit1.vp$coef[3];B5E=fit1.vp$coef[4];B6E=fit1.vp$coef[5]
#
# 
# fit1<-optim(c(B00,B1E,B4E,B5E,B6E),nloglik.M1.b.Cal.1)
#  
#  ### etc
  B00=fit2.vp$coef[1];parm1=0;B2E=fit2.vp$coef[2];parm3=0; B4E=fit2.vp$coef[3];B5E=fit2.vp$coef[4];B6E=fit2.vp$coef[5]
#
  fit2<-optim(c(B00,B2E,B4E,B5E,B6E),nloglik.M1.b.Cal.2)
#  
#  B00=fit3.vp$coef[1];parm1=0;parm2=0;B3E=fit3.vp$coef[2]; B4E=fit3.vp$coef[3];B5E=fit3.vp$coef[4];B6E=fit3.vp$coef[5]
##
#  fit3<-optim(c(B00,B3E,B4E,B5E,B6E),nloglik.M1.b.Cal.3)
##  
#  B00=fit12.vp$coef[1]; B1E=fit12.vp$coef[2];B2E=fit12.vp$coef[3];parm3=0;B4E=fit12.vp$coef[4];B5E=fit12.vp$coef[5];B6E=fit12.vp$coef[6]
##  
#  fit12<-optim(c(B00,B1E,B2E,B4E,B5E,B6E),nloglik.M1.b.Cal.12)
##  
#  B00=fit13.vp$coef[1]; B1E=fit13.vp$coef[2];parm2=0; B3E=fit13.vp$coef[3];B4E=fit13.vp$coef[4];B5E=fit13.vp$coef[5];B6E=fit13.vp$coef[6]
##  
#  fit13<-optim(c(B00,B1E,B3E,B4E,B5E,B6E),nloglik.M1.b.Cal.13)
##  
#  B00=fit23.vp$coef[1];parm1=0; B2E=fit23.vp$coef[2];B3E=fit23.vp$coef[3];B4E=fit23.vp$coef[4];B5E=fit23.vp$coef[5];B6E=fit23.vp$coef[6]
##  
#  fit23<-optim(c(B00,B2E,B3E,B4E,B5E,B6E),nloglik.M1.b.Cal.23)
} else {
 source(file="MethodFund_forEnsemble3STieFull.R")
#  B1E=-.8;B2E=-.8;B3E=0
#  parm1=B1E;parm2=B2E;parm3=B3E
 
  
  ### fit a cox model with coefficients with x(0),x(T), and T

#  B00=fit1.vp$coef[1]; B1E=fit1.vp$coef[2];parm2=0;parm3=0; B4E=fit1.vp$coef[3];

 
 #fit1<-optim(c(B00,B1E,B4E),nloglik.M1.b.Cal.1)
#  
#  ### etc
  B00=fit2.vp$coef[1];parm1=0;B2E=fit2.vp$coef[2];parm3=0; B4E=fit2.vp$coef[3];
#  
  fit2<-optim(c(B00,B2E,B4E),nloglik.M1.b.Cal.2)
#  
#  B00=fit3.vp$coef[1];parm1=0;parm2=0;B3E=fit3.vp$coef[2]; B4E=fit3.vp$coef[3];
##
#  fit3<-optim(c(B00,B3E,B4E),nloglik.M1.b.Cal.3)
##  
#  B00=fit12.vp$coef[1]; B1E=fit12.vp$coef[2];B2E=fit12.vp$coef[3];parm3=0;B4E=fit12.vp$coef[4];
##  
#  fit12<-optim(c(B00,B1E,B2E,B4E),nloglik.M1.b.Cal.12)
##  
#  B00=fit13.vp$coef[1]; B1E=fit13.vp$coef[2];parm2=0; B3E=fit13.vp$coef[3];B4E=fit13.vp$coef[4];
##  
#  fit13<-optim(c(B00,B1E,B3E,B4E),nloglik.M1.b.Cal.13)
##  
#  B00=fit23.vp$coef[1];parm1=0; B2E=fit23.vp$coef[2];B3E=fit23.vp$coef[3];B4E=fit23.vp$coef[4];
##  
#  fit23<-optim(c(B00,B2E,B3E,B4E),nloglik.M1.b.Cal.23)
###  

}
 save(X0.min,
 #fit1,fit2,fit3,fit12,fit13,fit23,
 fit2,
file=paste("output/outd3_Scale_LRT2_event",ind.event,"_marker",ind.marker,"_region",ind.region,".Rdata",sep=''))
