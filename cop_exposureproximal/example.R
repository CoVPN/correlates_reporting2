library(WeMix)

### load the trial data (data) and the immune correlate study data in long form (datlong)

load(file="~/TND/example.Rdata")

## Using vaccine recipients with >1 observation in the immune correlates study to fit the LME model for immune response trajectory
qq=tapply(rep(1,nrow(datlong)),datlong$id,sum)
qq.id<-as.numeric(names(qq)[qq>1])

datlong.2<-datlong[datlong$id%in%qq.id,]
####
### here account for sampling weight

dd<-datlong.2[datlong.2$Z==1,]

###############

dd$W1<-1
fit<-mix(X~time + W + (1+time|id), data=dd,weights=c('W1','weights'),cWeights=TRUE)
coef.fix<-fit$coef             ### fixed effect 
coef.rd<-unlist(fit$varVC[2])  ### random effect variance matrix             
coef.rd<-matrix(coef.rd,byrow=2,nrow=2)
sigma.e<-fit$sigma    ### measurement error standard deviation


####### Now process the trial data for entering the Cox model 

ID<-data$ID
SE<-data$SE
TE<-data$TE
WE<-data$WE
WEc<-data$WEc
ZE<-data$ZE
RE<-data$RE
XE<-data$XE
delta=data$delta
weight=data$weight
incoh=data$incoh

SEE<-ifelse(delta==1, SE, SE+max(SE))
oo<-order(SEE)

SE<-SE[oo]
TE<-TE[oo]
XE<-XE[oo]
WE<-WE[oo]
WEc<-WEc[oo]
ZE<-ZE[oo]
RE<-RE[oo]
ll=sum(delta)
delta=delta[oo]
incoh=incoh[oo]
ID=ID[oo]
weight=weight[oo]

### 

oo<-ZE==0 | (ZE==1 & incoh==1)
SE<-SE[oo]
TE<-TE[oo]
XE<-XE[oo]
WE<-matrix(WE[oo])
WEc<-matrix(WEc[oo])
ZE<-ZE[oo]
RE<-RE[oo]
delta=delta[oo]
incoh=incoh[oo]
ID=ID[oo]
weight=weight[oo]
## unique value for SE[delta==1]

SE.u<-unique(SE[delta==1])
ll.u<-length(SE.u)

#### Get data from vaccine arm only ####################

oo.v=ZE==1
SE.v=SE[oo.v]
TE.v=TE[oo.v]
XE.v=XE[oo.v]
WE.v=matrix(WE[oo.v])
WEc.v=matrix(WEc[oo.v])
ZE.v=ZE[oo.v]
RE.v=RE[oo.v]
delta.v=delta[oo.v]
ll.v=sum(delta.v)
incoh.v=incoh[oo.v]
ID.v=ID[oo.v]
weight.v=weight[oo.v]

### unique value for SE.v[delta.v==1]

SE.uv<-sort(unique(SE.v[delta.v==1]))


### Data for immune response history from the immune correlate study

X.long<-datlong$X
T.long<-datlong$time
id.long<-datlong$id

#### create a wide-form data with one observation each participant that include multiple time point with immune response measure

T.visit<-c(1,2,3)

datwide.delta0 <- reshape(datlong[c('id','visit','X')], idvar = "id", 
timevar = "visit",v.names="X",direction = "wide")
datwide.delta0$groupall<-apply(!is.na(datwide.delta0[,-1]),1,sum,na.rm=T)

Vname<-paste('X.',T.visit,sep='')
datwide.delta0$group<-apply(!is.na(datwide.delta0[,Vname]),1,sum,na.rm=T)


datwide.delta0.time <- reshape(datlong[c('id','visit','time')], idvar = "id", 
timevar = "visit",v.names="time",direction = "wide")

datwide.delta0<-merge(datwide.delta0,datwide.delta0.time,by="id",all=T)
Tname<-paste('time.',T.visit,sep='')


#################


weights.v=weight.v  ### weight to enter into Cox model for vaccine recipients only 
weights=weight      ### weight to enter into Cox model for both vaccine and placebo recipients


 dw<-1 ### dimension of covariate to adjust for in the LME
 dc<-1 ### dimension of confounders to adjust for in the Cox model
  
  
 B2E<-(-1)  ## initial value for coef for x(t) in the Cox model
 B1E<-0     ## initial value for coef for x(0) in the Cox model
 B3E<-0     ## initial value for coef for t in the Cox model
 B4E<-0     ## initial value for confounder in the Cox model


  
#### Derive vaccine-only estimator for Cox model coefficients  
  
  source(file="~/TND/code/May2021/Revision/Fun_Vacc.R")
  
  parm1=0;parm2=0;parm3=0;parm4=0
 
 ## fit Cox model with x(t) and W
  fit2<-optim(c(B2E,B4E),nloglik.M1.a.Cal.2)
  
        #  > fit2
        #$par
        #[1] -0.82832007 -0.07885602
        #
        #$value
        #[1] 9984.025
        #
        #$counts
        #function gradient
        #      35       NA
        #
        #$convergence
        #[1] 0
        #
        #$message
        #NULL

  ##fit Cox model with x(0), x(t) and W
  
  fit12<-optim(c(B1E,B2E,B4E),nloglik.M1.a.Cal.12)
  
        # > fit12
        #$par
        #[1]  0.2923914 -0.9120382 -0.1009345
        #
        #$value
        #[1] 9982.056
        #
        #$counts
        #function gradient
        #      70       NA
        #
        #$convergence
        #[1] 0
        #
        #$message
        #NULL

  
  ## fit Cox model with x(t), t, and W
  
  fit23<-optim(c(B2E,B3E,B4E),nloglik.M1.a.Cal.23)
  
        #     > fit23
        #$par
        #[1] -0.79119577  0.59681177 -0.08408344
        #
        #$value
        #[1] 9982.567
        #
        #$counts
        #function gradient
        #     108       NA
        #
        #$convergence
        #[1] 0
        #
        #$message
        #NULL
        
        
   #### Derive vaccine+placebo estimator for Cox model coefficients  
     
      source(file="~/TND/code/May2021/Revision/Fun_VaccPlac.R")
      
       B00=-0.1  ## initial value for Z
       B1E<-0    ## initial value for Z*x(0)
       B2E<-(-1) ## initial value for Z*x(t) 
       B3E<-0    ## initial value for Z*t
       B4E=0     ## initial value for confounder
        
       parm0=B00; parm1=parm2=parm3=parm4=0;
 
        
   ## fit Cox model with Z,Z*x(t), and W
   
       fit2<-optim(c(B00,B2E,B4E),nloglik.M1.b.Cal.2)
       
        #      > fit2
        #$par
        #[1]  1.87535464 -0.98543414 -0.03119901
        #
        #$value
        #[1] 27531.17
        #
        #$counts
        #function gradient
        #     100       NA
        #
        #$convergence
        #[1] 0
        #
        #$message
        #NULL


       
   ## fit Cox model with Z, Z*x(0),Z*x(t), and W    
     fit12<-optim(c(B00,B1E,B2E,B4E),nloglik.M1.b.Cal.12)


        #> fit12
        #$par
        #[1]  0.68068529  0.39252299 -1.05778330 -0.04391891
        #
        #$value
        #[1] 27526.58
        #
        #$counts
        #function gradient
        #     251       NA
        #
        #$convergence
        #[1] 0
        #
        #$message
        #NULL


   ## fit Cox model with Z, Z*x(t),Z*t, and W    
    fit23<-optim(c(B00,B2E,B3E,B4E),nloglik.M1.b.Cal.23)
  
        # 
        #  > fit23
        #$par
        #[1]  0.96106315 -0.80391443  0.74776767 -0.03903916
        #
        #$value
        #[1] 27521.4
        #
        #$counts
        #function gradient
        #     277       NA
        #
        #$convergence
        #[1] 0
        #
        #$message
        #NULL
