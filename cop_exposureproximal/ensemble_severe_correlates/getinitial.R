
ind.event<-as.integer(commandArgs(trailingOnly=T)[1])
ind.marker<-as.integer(commandArgs(trailingOnly=T)[2])
ind.region<-as.integer(commandArgs(trailingOnly=T)[3])

#print(ind.event, ind.marker,ind.region)
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

oc<-3

excludeOutlierC<-function(x,coef){
 Q3<-quantile(x,.75,na.rm=T)
 Q1<-quantile(x,.25,na.rm=T)
 IQR.x<-IQR(x,na.rm=T)
 x.out<-ifelse(x>Q3+coef*IQR.x,NA,x)
 x.out<-ifelse(x<Q1-coef*IQR.x,NA,x.out)
 return(x.out)
}

setwd("~/COVPN_P3003")

source(file="~/Pepe/RA/Functions/FunctionCall.R")
library(lme4)

#setwd("C:/Users/yhuang/OneDrive - Fred Hutchinson Cancer Research Center/Documents/All_Files/1yingsstuff/COVPN_P3003")
#source(file="../Pepe/RA/Functions/FunctionCall.R")


vv.endi<-"EventIndPrimaryIncludeNotMolecConfirmedD29_NoRegionCens" 
vv.endt<-"EventTimePrimaryIncludeNotMolecConfirmedD29_NoRegionCens"

dat<-read.csv("adata/janssen_pooled_partA_data_processed_with_riskscore_20240305.csv",na.strings=c("n/a","NA","N/A","","."))

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

#dat[oo.wdate,c('NumberdaysD1toD29','NumberdaysD1toD71','NumberdaysD1toM6')]

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


###
data$Sex<-ifelse(data$Sex==2,0.5,data$Sex)


data$weight<-ifelse(data$Trt==1,data$wt.D29,1)

data$Eventime<-data[,vv.time[ind.event]]
data$Eventind<-data[,vv.delta[ind.event]]
data$X0<-data[,vv.D29[ind.marker]]


X0.min<-min(data$X0,na.rm=T)
data$X0=data$X0-X0.min

RE<-data$CalendarDateEnrollment+data$NumberdaysD1toD29
TE<-data[,vv.time[ind.event]]
SE<-RE+TE
ZE=data$Trt
XE0=ifelse(data$Trt==1, data[,vv.D29[ind.marker]],NA)
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

ID<-1:nrow(data)
#incoh<-corr_dat$Trt==1 & !is.na(corr_dat$Day57pseudoneutid50)
incoh<-data$Trt==1 & data$ph2.D29 & !is.na(data[,vv.D29[ind.marker]])


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
   XE.truc<-XE.truc-X0.min

    #TE.truc<-c(TE.truc,TE[i])
    #XE.truc<-c(XE.truc,XE[i])
    SE.truc<-TE.truc+RE[i]
    ll=length(TE.truc)
  #  W1.truc<-rep(WE1[i],ll)
  #  W2.truc<-rep(WE2[i],ll)
  #  W3.truc<-rep(WE3[i],ll)
    datlong<-rbind(datlong,cbind(rep(i,ll),rep(RE[i],ll),TE.truc,XE.truc,SE.truc,rep(ZE[i],ll),rep(SE[i],ll),rep(0,ll),1:ll,rep(delta[i],ll),rep(weight[i],ll),
    rep(data$Age[i],ll),rep(data$BMI[i],ll),rep(data$Sex[i],ll),rep(data$Region[i],ll)))
    
    ########3
    
    tt.temp<-ww.temp<-rep(NA,3)
    ifelse(T.visit<TE[i],T.visit,NA)
    tt.temp[oo.visit]<-TE.truc
    ww.temp[oo.visit]<-XE.truc
    tt[i,]<-tt.temp
    ww[i,]<-ww.temp
    
}

tt[ZE==0,]<-c(0,NA,NA)
ww[ZE==0,]<-c(0,NA,NA)



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

  
  
datlong<-as.data.frame(datlong)

### delta is censoring indicator ##
colnames(datlong)<-c('id','R','time','X','S','Z','eventime','delta','visit','event','weight','Age','BMI','Sex','Region')
datlong<-datlong[datlong$Z==1 & !is.na(datlong$X),]





##################

X.long<-datlong$X
T.long<-datlong$time
id.long<-datlong$id

###############

qq=tapply(rep(1,nrow(datlong)),datlong$id,sum)
qq.id<-as.numeric(names(qq)[qq>1])

datlong.2<-datlong[datlong$id%in%qq.id,]
####


#
### Here we just use random subcohort to estimate LME


#fit<-lmer(X~time+Age+Sex+(1|id), data=datlong.2[datlong.2$Z==1 & datlong.2$delta==0,])
#coef.fix<-getME(fit,name=c("fixef"))
#coef.rd<-as.matrix(VarCorr(fit)[[1]])
#sigma.e<-getME(fit,name=c("sigma"))

#coef.fix<-getME(fit,name=c("fixef"))
#slope<-coef.fix[2]

### here account for sampling weight

dd<-datlong.2[datlong.2$Z==1,]
dd$W1<-1
fit<-mix(X~time + Age+Sex + (1|id), data=dd,weights=c('W1','weight'),cWeights=TRUE)
coef.fix<-fit$coef
slope<-coef.fix[2]


##########

data$daystart<-data$CalendarDateEnrollment+data$NumberdaysD1toD29
data$daystop<-data$Eventime+data$daystart

data<-data[,c('Ptid','X0','Eventime','Eventind','daystart','daystop','Region','risk_score','weight','Trt')]

make_dat_tveff = function(dat, slope = -0.0043) {
        dat %>% 
        group_by(Ptid) %>% 
        summarise(tstart = seq(daystart, daystop - 1, by = 1),
                  tstop = tstart + 1,
                  daystart = unique(daystart),
                  daystop=unique(daystop),
                  event = case_when(tstop == daystop ~ Eventind,
                                    TRUE ~ 0), ### event takes value of EventIndPrimaryD57 if tstop==daystop, otherwise event=0
                  Trt = unique(Trt),
                #  tps_stratum = unique(tps_stratum),
                #  Wstratum = unique(Wstratum),
                #  SubcohortInd = unique(SubcohortInd),
                  Region = unique(Region),
                #  HighRiskInd = unique(HighRiskInd),
                  risk_score = unique(risk_score),
                #  wt_D57 = unique(wt_D57),
                #  day57_titre = unique(Day57pseudoneutid50),
                X0=unique(X0),
                weight=unique(weight),
                XE = unique(X0) + slope * (tstart - daystart),
                Eventind=unique(Eventind),
                Eventime=unique(Eventime),
                  .groups = "keep") %>% 
        ungroup()
}

 dat_tveff= make_dat_tveff(data, slope = slope)
 
 dat_tveff$XE<-ifelse(dat_tveff$Trt==0,0,dat_tveff$XE)
 dat_tveff$X0<-ifelse(dat_tveff$Trt==0,0,dat_tveff$X0)
 dat_tveff$weight<-ifelse(dat_tveff$Trt==0,1,dat_tveff$weight)
 
   
    #dat_tveff<-merge(dat_tveff,corr_dat[,c('Ptid','daystart')],by="Ptid",all.x=T)
    
    dat_tveff$te<-dat_tveff$tstop-dat_tveff$daystart
    
    # vaccine model
    if (ind.region==3){
    fit1.v = 
        dat_tveff[dat_tveff$Trt==1,] %>% 
        coxph(Surv(tstart, tstop, event) ~ X0+risk_score+ factor(Region),
              data = .,
              id = Ptid,
              weights = dat_tveff[dat_tveff$Trt==1,]$weight,
              model = TRUE)
    
     fit2.v = 
        dat_tveff[dat_tveff$Trt==1,] %>% 
        coxph(Surv(tstart, tstop, event) ~ XE+risk_score+ factor(Region),
              data = .,
              id = Ptid,
              weights = dat_tveff[dat_tveff$Trt==1,]$weight,
              model = TRUE)
    
   fit3.v = 
        dat_tveff[dat_tveff$Trt==1,] %>% 
        coxph(Surv(tstart, tstop, event) ~ te+risk_score+ factor(Region),
              data = .,
              id = Ptid,
              weights = dat_tveff[dat_tveff$Trt==1,]$weight,
              model = TRUE)

     fit12.v = 
        dat_tveff[dat_tveff$Trt==1,] %>% 
        coxph(Surv(tstart, tstop, event) ~ X0+XE+risk_score+ factor(Region),
              data = .,
              id = Ptid,
              weights = dat_tveff[dat_tveff$Trt==1,]$weight,
              model = TRUE)

    
       fit13.v = 
        dat_tveff[dat_tveff$Trt==1,] %>% 
        coxph(Surv(tstart, tstop, event) ~ X0+te+risk_score+ factor(Region),
              data = .,
              id = Ptid,
              weights = dat_tveff[dat_tveff$Trt==1,]$weight,
              model = TRUE)
    
    
        fit23.v = 
        dat_tveff[dat_tveff$Trt==1,] %>% 
        coxph(Surv(tstart, tstop, event) ~ XE+te+risk_score+ factor(Region),
              data = .,
              id = Ptid,
              weights = dat_tveff[dat_tveff$Trt==1,]$weight,
              model = TRUE)
    
     
   fit123.v = 
        dat_tveff[dat_tveff$Trt==1,] %>% 
        coxph(Surv(tstart, tstop, event) ~ X0+XE+te+risk_score+ factor(Region),
              data = .,
              id = Ptid,
              weights = dat_tveff[dat_tveff$Trt==1,]$weight,
              model = TRUE)

 
      
    fit1.vp = 
        dat_tveff %>% 
        coxph(Surv(tstart, tstop, event) ~ Trt+X0+risk_score+ factor(Region),
              data = .,
              id = Ptid,
              weights = dat_tveff$weight,
              model = TRUE)
    
  
  fit2.vp = 
        dat_tveff %>% 
        coxph(Surv(tstart, tstop, event) ~ Trt+XE+risk_score+ factor(Region),
              data = .,
              id = Ptid,
              weights = dat_tveff$weight,
              model = TRUE)
   
   
  fit3.vp = 
        dat_tveff %>% 
        coxph(Surv(tstart, tstop, event) ~ Trt+I(Trt*te)+risk_score+ factor(Region),
              data = .,
              id = Ptid,
              weights = dat_tveff$weight,
              model = TRUE)

     fit12.vp = 
        dat_tveff %>% 
        coxph(Surv(tstart, tstop, event) ~ Trt+X0+XE+risk_score+ factor(Region),
              data = .,
              id = Ptid,
              weights = dat_tveff$weight,
              model = TRUE)

              
          fit13.vp = 
        dat_tveff %>% 
        coxph(Surv(tstart, tstop, event) ~ Trt+X0+I(Trt*te)+risk_score+ factor(Region),
              data = .,
              id = Ptid,
              weights = dat_tveff$weight,
              model = TRUE)
              
        fit23.vp = 
        dat_tveff %>% 
        coxph(Surv(tstart, tstop, event) ~ Trt+XE+I(Trt*te)+risk_score+ factor(Region),
              data = .,
              id = Ptid,
              weights = dat_tveff$weight,
              model = TRUE)
     fit123.vp = 
        dat_tveff %>% 
        coxph(Surv(tstart, tstop, event) ~ Trt+X0+XE+I(Trt*te)+risk_score+ factor(Region),
              data = .,
              id = Ptid,
              weights = dat_tveff$weight,
              model = TRUE)
     } else {
     
     fit1.v = 
        dat_tveff[dat_tveff$Trt==1 & dat_tveff$Region==ind.region,] %>% 
        coxph(Surv(tstart, tstop, event) ~ X0+risk_score,
              data = .,
              id = Ptid,
              weights = dat_tveff[dat_tveff$Trt==1 & dat_tveff$Region==ind.region, ]$weight,
              model = TRUE)
    
     fit2.v = 
        dat_tveff[dat_tveff$Trt==1 & dat_tveff$Region==ind.region,] %>% 
        coxph(Surv(tstart, tstop, event) ~ XE+risk_score,
              data = .,
              id = Ptid,
              weights = dat_tveff[dat_tveff$Trt==1 & dat_tveff$Region==ind.region,]$weight,
              model = TRUE)
    
   fit3.v = 
        dat_tveff[dat_tveff$Trt==1 & dat_tveff$Region==ind.region,] %>% 
        coxph(Surv(tstart, tstop, event) ~ te+risk_score,
              data = .,
              id = Ptid,
              weights = dat_tveff[dat_tveff$Trt==1& dat_tveff$Region==ind.region,]$weight,
              model = TRUE)

     fit12.v = 
        dat_tveff[dat_tveff$Trt==1& dat_tveff$Region==ind.region,] %>% 
        coxph(Surv(tstart, tstop, event) ~ X0+XE+risk_score,
              data = .,
              id = Ptid,
              weights = dat_tveff[dat_tveff$Trt==1 & dat_tveff$Region==ind.region,]$weight,
              model = TRUE)

    
       fit13.v = 
        dat_tveff[dat_tveff$Trt==1& dat_tveff$Region==ind.region,] %>% 
        coxph(Surv(tstart, tstop, event) ~ X0+te+risk_score,
              data = .,
              id = Ptid,
              weights = dat_tveff[dat_tveff$Trt==1 & dat_tveff$Region==ind.region,]$weight,
              model = TRUE)
    
    
        fit23.v = 
        dat_tveff[dat_tveff$Trt==1 & dat_tveff$Region==ind.region,] %>% 
        coxph(Surv(tstart, tstop, event) ~ XE+te+risk_score,
              data = .,
              id = Ptid,
              weights = dat_tveff[dat_tveff$Trt==1& dat_tveff$Region==ind.region,]$weight,
              model = TRUE)
    
     
   fit123.v = 
        dat_tveff[dat_tveff$Trt==1 & dat_tveff$Region==ind.region,] %>% 
        coxph(Surv(tstart, tstop, event) ~ X0+XE+te+risk_score,
              data = .,
              id = Ptid,
              weights = dat_tveff[dat_tveff$Trt==1& dat_tveff$Region==ind.region,]$weight,
              model = TRUE)

 
      
    fit1.vp = 
        dat_tveff[dat_tveff$Region==ind.region,] %>% 
        coxph(Surv(tstart, tstop, event) ~ Trt+X0+risk_score,
              data = .,
              id = Ptid,
              weights = dat_tveff[dat_tveff$Region==ind.region,]$weight,
              model = TRUE)
    
  
  fit2.vp = 
        dat_tveff[dat_tveff$Region==ind.region,] %>% 
        coxph(Surv(tstart, tstop, event) ~ Trt+XE+risk_score,
              data = .,
              id = Ptid,
              weights = dat_tveff[dat_tveff$Region==ind.region,]$weight,
              model = TRUE)
   
   
  fit3.vp = 
        dat_tveff[dat_tveff$Region==ind.region,] %>% 
        coxph(Surv(tstart, tstop, event) ~ Trt+I(Trt*te)+risk_score,
              data = .,
              id = Ptid,
              weights = dat_tveff[dat_tveff$Region==ind.region,]$weight,
              model = TRUE)

     fit12.vp = 
        dat_tveff[dat_tveff$Region==ind.region,] %>% 
        coxph(Surv(tstart, tstop, event) ~ Trt+X0+XE+risk_score,
              data = .,
              id = Ptid,
              weights = dat_tveff[dat_tveff$Region==ind.region,]$weight,
              model = TRUE)

              
          fit13.vp = 
        dat_tveff[dat_tveff$Region==ind.region,] %>% 
        coxph(Surv(tstart, tstop, event) ~ Trt+X0+I(Trt*te)+risk_score,
              data = .,
              id = Ptid,
              weights = dat_tveff[dat_tveff$Region==ind.region,]$weight,
              model = TRUE)
              
        fit23.vp = 
        dat_tveff[dat_tveff$Region==ind.region,] %>% 
        coxph(Surv(tstart, tstop, event) ~ Trt+XE+I(Trt*te)+risk_score,
              data = .,
              id = Ptid,
              weights = dat_tveff[dat_tveff$Region==ind.region,]$weight,
              model = TRUE)
     fit123.vp = 
        dat_tveff[dat_tveff$Region==ind.region,] %>% 
        coxph(Surv(tstart, tstop, event) ~ Trt+X0+XE+I(Trt*te)+risk_score,
              data = .,
              id = Ptid,
              weights = dat_tveff[dat_tveff$Region==ind.region,]$weight,
              model = TRUE)
     
     
     
     }
    
save(slope, fit1.v,fit2.v,fit3.v,fit12.v,fit13.v,fit23.v,fit123.v,fit1.vp,fit2.vp,fit3.vp,fit12.vp,fit13.vp,fit23.vp,
fit123.vp,file=paste("~/TND/Ensemble/Result/ES_event",ind.event,"_marker",ind.marker,"_region",ind.region,".Rdata",sep=''))
