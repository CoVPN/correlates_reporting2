# setwd("C:/Users/yhuang/OneDrive - Fred Hutchinson Cancer Research Center/Documents/All_Files/1yingsstuff/TND/Ensemble")
source(file="../FunctionCall.R")


transf=function(y) -log(1-y) 

yat=c(seq(0,.9,by=.1),.95,.99)

SS.max.3<-c(log10(266),log10(279),log10(123))
SS.max.0<-c(log10(209),log10(201),log10(157))
SS.max.1<-c(log10(258),log10(279),log10(97.8))
SS.max.2<-c(log10(326),log10(363),log10(89.7))

vv1=c('Anti Spike IgG (BAU/ml)','Anti RBD IgG (BAU/ml)','nAb-ID50 (IU50/ml)')
vv2=c('COVID-19','Severe-Critical COVID-19','Moderate COVID-19')
ind.event<-1
ind.marker<-1
ind.region<-0

pdf(file=paste0("Graph/VE2_Scale_LRT2_event.pdf"),width=7, height=7)


par(mfcol=c(2,2))
par(mar=c(3.5, 4, 1.8, 1) + 0.1)
par(mar=c(3.5, 3, 1.8, 1) + 0.1)

ind.region<-3

for (ind.event in 2:3){
#for (ind.region in 3:3){


#pdf(file=paste0("Graph/VE2_LRT2_event",ind.event,"_marker",ind.marker,"_region",ind.region,".pdf"),width=7, height=7)



#SS.maxx=SS.max.0*(ind.region==0)+SS.max.1*(ind.region==1)+SS.max.2*(ind.region==2)+SS.max.3*(ind.region==3)


for (ind.marker in c(3,1)){

#SS.max<-SS.maxx[ind.marker]


load(file=paste("output/outVE_Scale_LRT2_event",ind.event,"_marker",ind.marker,"_region",ind.region,".Rdata",sep=''))


SS.min<-min(SS)

oo<-SS<=SS.max
   
    sde<-apply(log(1-out.VE.b),2,sd,na.rm=T)
    
    low.ci.per<-1-exp(log(1-VE.SS)+qnorm(0.975)*sde)
    high.ci.per<-1-exp(log(1-VE.SS)-qnorm(0.975)*sde)
    
    low.ci<-apply(out.VE.b,2,quantile,0.025)
    high.ci<-apply(out.VE.b,2,quantile,0.975)


#plot(SS,VE.SS,ylim=c(0,1),type='n')
#lines(SS,VE.SS)
#lines(SS,low.ci.per,col=2)
#lines(SS,high.ci.per,col=2)
#lines(SS,low.ci,col=3)
#lines(SS,high.ci,col=3)


low.est<-low.ci
high.est<-high.ci

x.vv=vv1[ind.marker]
y.vv=vv2[ind.event]
#plot(SS[oo],transf(low.est[oo]),type='n',xlim=c(SS.min,SS.max*1.02),ylim=transf(c(0,min(.99,max(high.est[oo])*1.02))), xaxt="n", yaxt="n", lwd=2.5,xlab='',ylab='',main='',cex.main=1.2)
plot(c(SS[oo],max(SS[oo])*1.01),c(transf(low.est[oo]),transf(.95)),type='n',xlim=c(SS.min,SS.max*1.02),ylim=transf(c(0,.95)), xaxt="n", yaxt="n", lwd=2.5,xlab='',ylab=vv2[ind.event],main='',cex.main=1.2)

#lines(Su,high.un)
#axis(side=1, at=seq(0,8,by=0.5),labels=F,cex.axis=0.8)
  axis(side=1, at=c(0,log10(5),1,log10(30)),labels=c(1,5,10,30),cex.axis=1)
#if (ind.marker==3) axis(side=1, at=c(0,1,log10(30)),labels=F,cex.axis=0.8)
  #axis(1, at=seq(0,8,by=0.5), labels=seq(0,8,by=0.5), tick=F, line=0,cex.axis=1,las=0)
  #axis(side=1, at=seq(0,5,by=1), labels=c(1,10,expression(10^2),10^3,10^4,10^5)), tick=F, line=0,cex.axis=1,las=0)
  #axis(side=1, at=1, labels=expression(10^2), tick=F, line=0,cex.axis=1,las=0)
 # axis(side=1, at=c(0,log10(5),1,log10(30)), labels=c(1,5,10,30), tick=F, line=0,cex.axis=1,las=0)
#if (ind.marker==3)  axis(side=1, at=c(log10(30)), labels=c(30), tick=F, line=0,cex.axis=1,las=0)
  axis(side=1, at=2, labels=expression(10^2), line=0,cex.axis=1,las=0)
  axis(side=1, at=3, labels=expression(10^3), line=0,cex.axis=1,las=0)
  axis(side=1, at=4, labels=expression(10^4), line=0,cex.axis=1,las=0)
  axis(side=1, at=5, labels=expression(10^5), line=0,cex.axis=1,las=0)
   
  #mtext(expression(bold(s[1])),side=1, las=0, line=2.5, cex=1)
  #mtext(bquote(s[1]==.(x.vv)),side=1, las=0, line=2.5, cex=1)
  mtext(paste(x.vv,' at exposure',sep=''),side=1, las=0, line=2.5, cex=.8)
  #axis(side=2, at=seq(0,1,by=0.25), labels=paste0(seq(0,100,by=25)),cex.axis=0.8)
  axis(side=2,at=transf(yat),labels=(yat*100)%.%"%",cex=.8) 

  #mtext(expression(bold(paste('VE(',s[1],') (%)'))),side=2, las=0, line=2, cex=1,at=transf(0.9))
 # mtext(expression(paste('VE againt',bquote(.(y.vv)),'by 170 Days Post Day 29 Visit')),side=2, las=0, line=2, cex=1,at=transf(0.9))
 #mtext(paste('VE againt',y.vv,'by 170 Days Post Day 29 Visit'),side=2, las=0, line=2, cex=.8)
 mtext(paste('VE againt',y.vv),side=2, las=0, line=2, cex=.8)
#polygon(c(Su,rev(SS)), c(transf(high.un) ,rev(transf(low.un))), col = 'lightskyblue1',border=NA )
#polygon(c(SS,rev(SS)), c(transf(high.est) ,rev(transf(low.est))), col = 'lightskyblue1',border=NA )


#polygon(c(Su,rev(Su)), c(high.ci.per.A3,rev(low.ci.per.A3)), col = 'green',border=NA )
#polygon(c(Su,rev(Su)), c(high.est ,rev(low.est)), col = 'black',border=NA )

lines(SS[oo],transf(low.est[oo]),col=4,lty=2,lwd=2)
lines(SS[oo],transf(high.est[oo]),col=4,lty=2,lwd=2)
polygon(c(SS[oo],rev(SS[oo])), c(transf(high.est[oo]) ,rev(transf(low.est[oo]))), col = 'lightskyblue1',border=NA )
lines(SS[oo],transf(VE.SS[oo]),col=1,lty=1,lwd=2)

}
}
dev.off()
