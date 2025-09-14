
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

dat[oo.wdate,c('NumberdaysD1toD29','NumberdaysD1toD71','NumberdaysD1toM6')]

oo.immune<-dat$ph1.D29=='TRUE' & (dat$Trt==0 | dat$ph2.D29=='TRUE') & dat$Bserostatus==0

data<-dat[oo.immune,]

#ind.event<-1
#ind.marker<-1
#ind.region<-3


for (ind.event in 1:3){
for (ind.marker in 1:3){
for (ind.region in 0:3){

SS<-sort(unique(data[,vv.D29[ind.marker]]))


load(file=paste("~/TND/Ensemble/Result/outd3_Scale_LRT2_event",ind.event,"_marker",ind.marker,"_region",ind.region,".Rdata",sep=''))

calVE.2<-function(beta0,beta2,x) {
 return(1-exp(beta0+beta2*x))
}

VE.SS<-calVE.2(fit2$par[1],fit2$par[2],SS-X0.min)


out.VE.b<-NULL
for (iter in 1:500){

ff=try.error(load(file=paste("~/TND/Ensemble/Result/outd3_Scale_LRT2_event",ind.event,"_marker",ind.marker,"_region",ind.region,"_iter",iter,".Rdata",sep='')))
if (!inherits(ff,'try-error')){
 VE.SS.b<-calVE.2(fit2.b$par[1],fit2.b$par[2],SS-X0.min)
 out.VE.b<-rbind(out.VE.b,VE.SS.b)
 rm(ff)
 } else {
 print(c(ind.event,ind.marker,ind.region,iter))
 }
}
load(file=paste("~/TND/Ensemble/ResultVE/XEQ_Scale_event",ind.event,"_marker",ind.marker,"_region",ind.region,".Rdata",sep=''))
SS.max=q[3]
save(SS,VE.SS,out.VE.b,SS.max,file=paste("~/TND/Ensemble/ResultVE/outVE_Scale_LRT2_event",ind.event,"_marker",ind.marker,"_region",ind.region,".Rdata",sep=''))
rm(fit2,SS,VE.SS,SS.max)
}
}
}
