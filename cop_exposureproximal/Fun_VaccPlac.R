
##### This file include functions for estimating exposure-proximal correlate model based on data from vaccine+placebo recipients.
##### parm0 parm1 parm2 parm3 are logHR for Z, Z*x(0), Z*x(t), and Z*t respectively, parm4 is logHR for confounder WEc that is 
##### adjusted for in the Cox model. Here parm4 is a scalar assuming univariate WEc. Multivariate WEc can be handled by specifying parm4 to 
##### be a vector.


######## Fit a Cox hazard model on calendar time with coefficient for Z, Z*x(0), Z*x(t), Z*t, and W
### It requires the following input ###
    ### SE.u: ordered unique failure time
    ### ll: number of events
    ### datwide.delta0: a data in wide format showing immunogenicity measurements for vaccinees in the immune correlate study 
    ### Followings are vectors from the immune correlate study that is arranged in long format
        ### X.long: immune response measurement
        ### T.long: t (time post-peak) for immune response measurement
        ### id: participant id
    ### Follwing vectors correspond to information from vaccine recipient included in the analysis 
        ### ID: Unique ID that can be matched to "id" in datwide.delta0
        ### delta: case indicator
        ### SE: Minimum of Censoring and end of followup time, measured one calendar scale
        ### RE: Peak time point after vaccination, measured on calendar time scale
        ### TE: SE-RE: Time post-peak to SE
        ### WE: Covariates used to model immune response trajectory
        ### WEc: confounders to adjust for in Cox model

nloglik.M1.b.Cal<- function(parm){

 out=0
 
 outi<-NULL
   ### coef for Z, Z*x(0), Z*x(t), Z*t, and confounders respectively
    parm0=parm[1]
    parm1=parm[2]
    parm2=parm[3]
    parm3=parm[4]
    parm4=parm[5]
   
   jj.4.full<-(1:length(RE))[datwide.delta0[match(ID,datwide.delta0$id),]$groupall > datwide.delta0[match(ID,datwide.delta0$id),]$group]
   jj.4.full<-as.numeric(na.omit(jj.4.full))
   out.jj.4.full<-NULL
   for (j in jj.4.full){
      
    keep=1:ll
 
    VVj.11<-coef.rd[1,1]
    VVj.12<-VVj.21<-coef.rd[1,1]+(SE[keep]-RE[j])*coef.rd[1,2]
    VVj.22<-coef.rd[1,1]+2*(SE[keep]-RE[j])*coef.rd[1,2]+(SE[keep]-RE[j])^2*coef.rd[2,2]
    
    muj.1=cbind(1,0,matrix(WE[j,],length(j),dw))%*%coef.fix
    muj.2=cbind(1,SE[keep]-RE[j],matrix(WE[j,],length(j),dw))%*%coef.fix
    
   
    oo.j<-id.long==ID[j]
    Xj<-X.long[oo.j]
    Tj<-T.long[oo.j]
    
    
    VVj.c<-cbind(1,Tj)%*%coef.rd%*%rbind(1,Tj)+diag(rep(sigma.e^2,length(Tj)),nrow=length(Tj),ncol=length(Tj))               
   
    COVVj.1<-coef.rd[1,1]+coef.rd[1,2]*Tj
    COVVj.2<-coef.rd[1,1]+coef.rd[1,2]*outer(SE[keep]-RE[j],Tj,'+')+coef.rd[2,2]*outer(SE[keep]-RE[j],Tj,'*')
    COVVj.1=matrix(COVVj.1,1,length(COVVj.1))
    if (length(Tj)==1) COVVj.2=matrix(COVVj.2,length(COVVj.2),1)
    
    muj.c<-cbind(1,Tj,matrix(WE[j,],length(j),dw))%*%coef.fix
        
    comp1j<-cbind(as.numeric(muj.1),muj.2)+cbind(as.numeric(COVVj.1%*%solve(VVj.c)%*%(Xj-muj.c)),COVVj.2%*%solve(VVj.c)%*%(Xj-muj.c))
    
     VVjminus.11<-COVVj.1%*%solve(VVj.c)%*%(t(COVVj.1))
     VVjminus.12<-COVVj.2%*%solve(VVj.c)%*%(t(COVVj.1))
     VVjminus.22<-diag(COVVj.2%*%solve(VVj.c)%*%(t(COVVj.2)))
     
     comp2j.11<-VVj.11-VVjminus.11
     comp2j.12<-VVj.12-VVjminus.12
     comp2j.22<-VVj.22-VVjminus.22
   
     out.jj.4.full<-cbind(out.jj.4.full,as.numeric(c(parm1,parm2)%*%t(comp1j))+parm1^2*as.numeric(comp2j.11)/2+parm1*parm2*comp2j.12+parm2^2*comp2j.22/2+parm3*(SE[keep]-RE[j])
     +parm0*ZE[j]+matrix(WEc[j,],length(j),dc)%*%parm4)
     
    }
        
    S.tie=SE.u[1]  ### unique SE value
    ind.tie=0       ### position within the unique value
    
    for (i in 1:ll){
    
    ## track #of ties
    
        di<-sum(SE==SE[i] & delta==1)
        
        if (SE[i]==S.tie){  ## in the same tie as previous value
          ind.tie=ind.tie+1
        } else { ## go to a new SE value
           S.tie=SE[i]
           ind.tie=1
        }
     ### compute average weight among ties
        weights.bar=mean(weights[SE==SE[i] & delta==1])
     
    ss=0
    comp12i=0 ### placebo contribute nothing
    if (ZE[i]==1) { 
    VVi<-cbind(1,c(0,TE[i]))%*%coef.rd%*%rbind(1,c(0,TE[i]))
    mui<-rbind(c(1,0,matrix(WE[i,],length(i),dw)),c(1,TE[i],matrix(WE[i,],length(i),dw)))%*%coef.fix
   
    oo.i<-id.long==ID[i]
        
       #### (i) if case has at least one X(Ti)* measured, i.e. case ID[i] belong to id.long
        if (sum(oo.i)>0){ 
        Xi<-X.long[oo.i]
        Ti<-T.long[oo.i]
        
        VVi.c<-cbind(1,Ti)%*%coef.rd%*%rbind(1,Ti)+diag(rep(sigma.e^2,length(Ti)),nrow=length(Ti),ncol=length(Ti))               
        COVVi<-cbind(1,c(0,TE[i]))%*%coef.rd%*%rbind(1,Ti)
        mui.c<-cbind(1,Ti,matrix(WE[i,],length(Ti),dw,byrow=TRUE))%*%coef.fix
        
        comp1i<-mui+COVVi%*%solve(VVi.c)%*%(Xi-mui.c)      ## E(X(0),X(t)|X(Ti)*,W)
        comp2i<-VVi-COVVi%*%solve(VVi.c)%*%t(COVVi)        ## Var(X(0),X(t)|X(Ti)*,W)
        
        } else {
        #### (ii) if case does not have X(Ti)* measure, i.e. case ID[i] not in id.long
        
       comp1i<-mui            
       comp2i<-VVi            
       } 
 
     comp12i<-as.numeric(t(c(parm1,parm2))%*%comp1i+t(c(parm1,parm2))%*%comp2i%*%c(parm1,parm2)/2)
    
    } 
    
    jj.r<-((i+1):(length(RE)))[RE[(i+1):length(RE)]<SE[i] & SE[(i+1):length(RE)]>=SE[i] & ID[(i+1):length(RE)]%in%id.long]
    jj.l<-NULL
    if (i>1) jj.l<-(1:(i-1))[RE[1:(i-1)]<SE[i] & SE[1:(i-1)]>=SE[i] & ID[1:(i-1)]%in%id.long]
    jj<-c(jj.l,jj.r) ## among risk set  
    
    jj.tie<-(1:(length(RE)))[RE[1:length(RE)]<SE[i] & ID[1:length(RE)]%in%id.long & SE[1:length(SE)]==SE[i] & delta[1:length(delta)]==1]
   
    jj.W.r<-((i+1):(length(RE)))[RE[(i+1):length(RE)]<SE[i] & ZE[(i+1):length(ZE)]==1 & SE[(i+1):length(RE)]>=SE[i] & !ID[(i+1):length(RE)]%in%id.long]
    jj.W.l<-NULL
    if (i>1) jj.W.l<-(1:(i-1))[RE[1:(i-1)]<SE[i] & ZE[1:(i-1)]==1 & !ID[1:(i-1)]%in%id.long & SE[1:(i-1)]>=SE[i]]
    jj.W<-c(jj.W.l,jj.W.r)    

    
    jj.W.tie<-(1:(length(RE)))[RE[1:length(RE)]<SE[i] & ZE[1:length(ZE)]==1 & SE[1:length(RE)]==SE[i] & !ID[1:length(RE)]%in%id.long  & delta[1:length(delta)]==1]
 
  
    jj.Z0.r<-((i+1):length(RE))[RE[(i+1):length(RE)]<SE[i] & ZE[(i+1):length(ZE)]==0 & SE[(i+1):length(SE)]>=SE[i]]  
    jj.Z0.l<-NULL
    if (i>1) jj.Z0.l<-(1:(i-1))[RE[1:(i-1)]<SE[i] & ZE[1:(i-1)]==0 & SE[1:(i-1)]>=SE[i]]
    jj.Z0=c(jj.Z0.l,jj.Z0.r)

    jj.Z0.tie<-(1:length(RE))[RE[1:length(RE)]<SE[i] & ZE[1:length(ZE)]==0 & SE[1:length(SE)]==SE[i] & delta[1:length(delta)]==1]  
  
    jj.1.r<-((i+1):(length(RE)))[RE[(i+1):length(RE)]<SE[i] & datwide.delta0[match(ID[(i+1):length(RE)],datwide.delta0$id),]$group==1 &
       datwide.delta0[match(ID[(i+1):length(RE)],datwide.delta0$id),]$groupall==1 & SE[(i+1):length(RE)]>=SE[i]]
    jj.1.l<-NULL
    if (i>1) jj.1.l<-(1:(i-1))[RE[1:(i-1)]<SE[i] & datwide.delta0[match(ID[1:(i-1)],datwide.delta0$id),]$group==1 &
       datwide.delta0[match(ID[1:(i-1)],datwide.delta0$id),]$groupall==1 & SE[1:(i-1)]>=SE[i]]
    jj.1<-c(jj.1.l,jj.1.r) ## among risk set

    
    jj.1.tie<-(1:(length(RE)))[RE[1:length(RE)]<SE[i] & datwide.delta0[match(ID[1:length(RE)],datwide.delta0$id),]$group==1 &
       datwide.delta0[match(ID[1:length(RE)],datwide.delta0$id),]$groupall==1 & SE[1:length(SE)]==SE[i] & delta[1:length(delta)]==1]
   
   
    #######
 
    jj.2.r<-((i+1):(length(RE)))[RE[(i+1):length(RE)]<SE[i] & datwide.delta0[match(ID[(i+1):length(RE)],datwide.delta0$id),]$group==2 & 
       datwide.delta0[match(ID[(i+1):length(RE)],datwide.delta0$id),]$groupall==2 & SE[(i+1):length(SE)]>=SE[i]]
    jj.2.l<-NULL
    if (i>1)   jj.2.l<-(1:(i-1))[RE[1:(i-1)]<SE[i] & datwide.delta0[match(ID[1:(i-1)],datwide.delta0$id),]$group==2 & 
       datwide.delta0[match(ID[1:(i-1)],datwide.delta0$id),]$groupall==2 & SE[1:(i-1)]>=SE[i]]
    jj.2<-c(jj.2.l,jj.2.r)

    jj.2.tie<-(1:(length(RE)))[RE[1:length(RE)]<SE[i] & datwide.delta0[match(ID[1:length(RE)],datwide.delta0$id),]$group==2 & 
       datwide.delta0[match(ID[1:length(RE)],datwide.delta0$id),]$groupall==2 & SE[1:length(SE)]==SE[i] & delta[1:length(delta)]==1]
   

    ###
    jj.3.r<-((i+1):(length(RE)))[RE[(i+1):length(RE)]<SE[i] & datwide.delta0[match(ID[(i+1):length(RE)],datwide.delta0$id),]$group==3 & 
       datwide.delta0[match(ID[(i+1):length(RE)],datwide.delta0$id),]$groupall==3 & SE[(i+1):length(SE)]>=SE[i]]
    jj.3.l<-NULL
    if (i>1) jj.3.l<-(1:(i-1))[RE[1:(i-1)]<SE[i] & datwide.delta0[match(ID[1:(i-1)],datwide.delta0$id),]$group==3 & 
       datwide.delta0[match(ID[1:(i-1)],datwide.delta0$id),]$groupall==3 & SE[1:(i-1)]>=SE[i]]
    jj.3<-c(jj.3.l,jj.3.r)
   

    jj.3.tie<-(1:(length(RE)))[RE[1:length(RE)]<SE[i] & datwide.delta0[match(ID[1:length(RE)],datwide.delta0$id),]$group==3 & 
       datwide.delta0[match(ID[1:length(RE)],datwide.delta0$id),]$groupall==3 & SE[1:length(SE)]==SE[i] & delta[1:length(delta)]==1]
  
     jj.1<-as.numeric(na.omit(jj.1))
     jj.2<-as.numeric(na.omit(jj.2))
     jj.3<-as.numeric(na.omit(jj.3))
   
       jj.1.tie<-as.numeric(na.omit(jj.1.tie))
     jj.2.tie<-as.numeric(na.omit(jj.2.tie))
     jj.3.tie<-as.numeric(na.omit(jj.3.tie))
     
    ss=ss+sum((RE[jj.4.full]<SE[i])*(SE[jj.4.full]>=SE[i])*(!jj.4.full==i)*weights[jj.4.full]*exp(out.jj.4.full[i,]-comp12i-parm3*(ZE[i]*TE[i])
    -parm0*ZE[i]-matrix(WEc[i,],length(i),dc)%*%parm4))  ## risk set from jj.4.full
     
    ss=ss-(ind.tie-1)/di*sum((RE[jj.4.full]<SE[i])*(SE[jj.4.full]==SE[i])*(!jj.4.full==i)*weights[jj.4.full]*exp(out.jj.4.full[i,]-comp12i-parm3*(ZE[i]*TE[i])
    -parm0*ZE[i]-matrix(WEc[i,],length(i),dc)%*%parm4))  ## tie from jj.4.full
     

   if (length(jj.W)>0){
      ### No X measure 
    mu0.EXE<-cbind(1,0,matrix(as.numeric(WE[jj.W,]),length(jj.W),dw))%*%coef.fix
    mu.EXE<-cbind(1,SE[i]-RE[jj.W],matrix(as.numeric(WE[jj.W,]),length(jj.W),dw))%*%coef.fix
  
   comp1.b<-parm1*mu0.EXE+parm2*mu.EXE
       
   comp2.b<-parm1^2*coef.rd[1,1]+2*parm1*parm2*(coef.rd[1,1]+(SE[i]-RE[jj.W])*coef.rd[1,2])+
   parm2^2*(coef.rd[1,1]+2*(SE[i]-RE[jj.W])*coef.rd[1,2]+(SE[i]-RE[jj.W])^2*coef.rd[2,2])
   
   ss=ss+sum(weights[jj.W]*exp(as.numeric(comp1.b)+0.5*comp2.b-comp12i+parm3*(SE[i]-RE[jj.W]-ZE[i]*TE[i])
   +parm0*(ZE[jj.W]-ZE[i])-t(as.numeric(WEc[i,])-matrix(as.numeric(t(WEc[jj.W,])),dc,length(jj.W)))%*%parm4) )  
    
   } 
   ### minus ties
   if (length(jj.W.tie)>0){
   
        mu0.EXE<-cbind(1,0,matrix(as.numeric(WE[jj.W.tie,]),length(jj.W.tie),dw))%*%coef.fix
        mu.EXE<-cbind(1,SE[i]-RE[jj.W.tie],matrix(as.numeric(WE[jj.W.tie,]),length(jj.W.tie),dw))%*%coef.fix
      
       comp1.b<-parm1*mu0.EXE+parm2*mu.EXE
           
       comp2.b<-parm1^2*coef.rd[1,1]+2*parm1*parm2*(coef.rd[1,1]+(SE[i]-RE[jj.W.tie])*coef.rd[1,2])+
       parm2^2*(coef.rd[1,1]+2*(SE[i]-RE[jj.W.tie])*coef.rd[1,2]+(SE[i]-RE[jj.W.tie])^2*coef.rd[2,2])
       
       ss=ss-(ind.tie-1)/di*sum(weights[jj.W.tie]*exp(as.numeric(comp1.b)+0.5*comp2.b-comp12i+parm3*(SE[i]-RE[jj.W.tie]-ZE[i]*TE[i])
       +parm0*(ZE[jj.W.tie]-ZE[i])-t(as.numeric(WEc[i,])-matrix(as.numeric(t(WEc[jj.W.tie,])),dc,length(jj.W.tie)))%*%parm4) ) 
   }      

 
   ### Measure of X(0)* only
   
    #account for different t levels
   
   if (length(jj.1)>0){
   Xt1=datwide.delta0[match(ID[jj.1],datwide.delta0$id),Vname[1]]
    t1=datwide.delta0[match(ID[jj.1],datwide.delta0$id),Tname[1]]
 
    ## make t1 vector
   COVV.a<-matrix(c(coef.rd[1,1]+t1*coef.rd[1,2],coef.rd[1,1]+(t1+SE[i]-RE[jj.1])*coef.rd[1,2]+t1*(SE[i]-RE[jj.1])*coef.rd[2,2]),length(jj.1),2,byrow=F)
   SS.a<-sapply(coef.rd[1,1]+2*t1*coef.rd[1,2]+t1^2*coef.rd[2,2]+sigma.e^2,solve)
   
   slope.a<-COVV.a*as.numeric(SS.a)
   Sigma1.2.a=outer(t(slope.a),(COVV.a),'*'); 
   kappa1=diag(as.matrix(Sigma1.2.a[1,,,1])); kappa2=diag(as.matrix(Sigma1.2.a[1,,,2])); 
   kappa3=diag(as.matrix(Sigma1.2.a[2,,,1])); kappa4=diag(as.matrix(Sigma1.2.a[2,,,2]))
   
     mu0.EXE<-cbind(1,0,matrix(as.numeric(WE[jj.1,]),length(jj.1),dw))%*%coef.fix
      mu.EXE<-cbind(1,SE[i]-RE[jj.1],matrix(as.numeric(WE[jj.1,]),length(jj.1),dw))%*%coef.fix
  
       comp1<-parm1*(mu0.EXE+slope.a[,1]*(Xt1-mu0.EXE))+
       parm2*(mu.EXE+slope.a[,2]*(Xt1-mu0.EXE))
       
       comp2<-parm1^2*(coef.rd[1,1]-kappa1)+
       parm2*parm1*(2*coef.rd[1,1]+2*(SE[i]-RE[jj.1])*coef.rd[1,2]-kappa2-kappa3)+
       parm2^2*(coef.rd[1,1]+(SE[i]-RE[jj.1])^2*coef.rd[2,2]+2*(SE[i]-RE[jj.1])*coef.rd[1,2]-kappa4)
    
    ss=ss+sum(weights[jj.1]*exp(as.numeric(comp1)+0.5*comp2-comp12i+parm3*(SE[i]-RE[jj.1]-ZE[i]*TE[i])
    +parm0*(ZE[jj.1]-ZE[i])-t(as.numeric(WEc[i,])-matrix(as.numeric(t(WEc[jj.1,])),dc,length(jj.1)))%*%parm4) )  
   
   rm(Xt1)
   }

   ## minus ties

   if (length(jj.1.tie)>0){
   Xt1=datwide.delta0[match(ID[jj.1.tie],datwide.delta0$id),Vname[1]]
   
   t1=datwide.delta0[match(ID[jj.1.tie],datwide.delta0$id),Tname[1]]
   
    # make t1 vector
   COVV.a<-matrix(c(coef.rd[1,1]+t1*coef.rd[1,2],coef.rd[1,1]+(t1+SE[i]-RE[jj.1.tie])*coef.rd[1,2]+t1*(SE[i]-RE[jj.1.tie])*coef.rd[2,2]),length(jj.1.tie),2,byrow=F)
   SS.a<-sapply(coef.rd[1,1]+2*t1*coef.rd[1,2]+t1^2*coef.rd[2,2]+sigma.e^2,solve)
 
   
   slope.a<-COVV.a*as.numeric(SS.a)
   Sigma1.2.a=outer(t(slope.a),(COVV.a),'*'); 
   kappa1=diag(as.matrix(Sigma1.2.a[1,,,1])); kappa2=diag(as.matrix(Sigma1.2.a[1,,,2])); 
   kappa3=diag(as.matrix(Sigma1.2.a[2,,,1])); kappa4=diag(as.matrix(Sigma1.2.a[2,,,2]))
   
     mu0.EXE<-cbind(1,0,matrix(as.numeric(WE[jj.1.tie,]),length(jj.1.tie),dw))%*%coef.fix
      mu.EXE<-cbind(1,SE[i]-RE[jj.1.tie],matrix(as.numeric(WE[jj.1.tie,]),length(jj.1.tie),dw))%*%coef.fix
  
       comp1<-parm1*(mu0.EXE+slope.a[,1]*(Xt1-mu0.EXE))+
       parm2*(mu.EXE+slope.a[,2]*(Xt1-mu0.EXE))
       
       comp2<-parm1^2*(coef.rd[1,1]-kappa1)+
       parm2*parm1*(2*coef.rd[1,1]+2*(SE[i]-RE[jj.1.tie])*coef.rd[1,2]-kappa2-kappa3)+
       parm2^2*(coef.rd[1,1]+(SE[i]-RE[jj.1.tie])^2*coef.rd[2,2]+2*(SE[i]-RE[jj.1.tie])*coef.rd[1,2]-kappa4)
    
    ss=ss-(ind.tie-1)/di*sum(weights[jj.1.tie]*exp(as.numeric(comp1)+0.5*comp2-comp12i+parm3*(SE[i]-RE[jj.1.tie]-ZE[i]*TE[i])
    +parm0*(ZE[jj.1.tie]-ZE[i])-t(as.numeric(WEc[i,])-matrix(as.numeric(t(WEc[jj.1.tie,])),dc,length(jj.1.tie)))%*%parm4 ) )  
   
   rm(Xt1)
   }



   ###### Measure of X(0)* and X(t2)* 
     if (length(jj.2)>0){
     Xt1=datwide.delta0[match(ID[jj.2],datwide.delta0$id),Vname[1]]
     Xt2=datwide.delta0[match(ID[jj.2],datwide.delta0$id),Vname[2]]
     
       ## make t1, t2 vectors
      t1=datwide.delta0[match(ID[jj.2],datwide.delta0$id),Tname[1]]
      t2=datwide.delta0[match(ID[jj.2],datwide.delta0$id),Tname[2]]
         
     COVV.b1<-coef.rd[1,1]+t1*coef.rd[1,2]
     COVV.b2<-coef.rd[1,1]+t2*coef.rd[1,2]
     COVV.b3<-coef.rd[1,1]+(t1+SE[i]-RE[jj.2])*coef.rd[1,2]+t1*(SE[i]-RE[jj.2])*coef.rd[2,2]
     COVV.b4<-coef.rd[1,1]+t2*(SE[i]-RE[jj.2])*coef.rd[2,2]+(t2+SE[i]-RE[jj.2])*coef.rd[1,2]
     

     SS.b.1<-coef.rd[1,1]+2*t1*coef.rd[1,2]+t1*2*coef.rd[2,2]+sigma.e^2
     SS.b.2<-coef.rd[1,1]+2*t2*coef.rd[1,2]+t2^2*coef.rd[2,2]+sigma.e^2
     SS.b.12<-coef.rd[1,1]+(t1+t2)*coef.rd[1,2]+t1*t2*coef.rd[2,2]
     
     SS.b<-array(NA,dim=c(2,2,length(jj.2)))
     SS.b[1,1,]<-SS.b.1
     SS.b[1,2,]<-SS.b[2,1,]<-SS.b.12
     SS.b[2,2,]<-SS.b.2
     
     SS.b<-apply(SS.b,3,solve)  ### SS.b is a 4 by length(jj.2) matrix with each column the inverse of a matrix
 
     ### slope
    
     ### When t1 t2 are vectors, SS.b[1,], SS.b[2,], SS.b[3,], SS.b[4,] 
     ### corresponds [1,1],[2,1],[1,2],[2,2] of a 2x2 matrix
     
     gamma1=COVV.b1*SS.b[1,]+COVV.b2*SS.b[2,]
     gamma2=COVV.b1*SS.b[3,]+COVV.b2*SS.b[4,]
     gamma3=COVV.b3*SS.b[1,]+COVV.b4*SS.b[2,]
     gamma4=COVV.b3*SS.b[3,]+COVV.b4*SS.b[4,]
     

     ### Sigma1.2.b
     eta1=gamma1*COVV.b1+gamma2*COVV.b2
     eta2=gamma1*COVV.b3+gamma2*COVV.b4
     eta3=gamma3*COVV.b1+gamma4*COVV.b2
     eta4=gamma3*COVV.b3+gamma4*COVV.b4
     
      mu0.EXE<-cbind(1,0,matrix(as.numeric(WE[jj.2,]),length(jj.2),dw))%*%coef.fix
       mu.EXE<-cbind(1,SE[i]-RE[jj.2],matrix(as.numeric(WE[jj.2,]),length(jj.2),dw))%*%coef.fix
       mut1.EXE<-cbind(1,t1,matrix(as.numeric(WE[jj.2,]),length(jj.2),dw))%*%coef.fix
       mut2.EXE<-cbind(1,t2,matrix(as.numeric(WE[jj.2,]),length(jj.2),dw))%*%coef.fix
       
      
      comp1<-parm1*(mu0.EXE+gamma1*Xt1+gamma2*Xt2-gamma1*mut1.EXE-gamma2*mut2.EXE)+
        parm2*(mu.EXE+gamma3*Xt1+gamma4*Xt2-gamma3*mut1.EXE-gamma4*mut2.EXE)
      
      comp2<-parm1^2*(coef.rd[1,1]-eta1)+parm2*parm1*(2*coef.rd[1,1]+2*(SE[i]-RE[jj.2])*coef.rd[1,2]-eta2-eta3)+
      parm2^2*(coef.rd[1,1]+(SE[i]-RE[jj.2])^2*coef.rd[2,2]+2*(SE[i]-RE[jj.2])*coef.rd[1,2]-eta4)
      
    
    ss=ss+sum(weights[jj.2]*exp(as.numeric(comp1)+0.5*comp2-comp12i+parm3*(SE[i]-RE[jj.2]-ZE[i]*TE[i])
    +parm0*(ZE[jj.2]-ZE[i])-t(as.numeric(WEc[i,])-matrix(as.numeric(t(WEc[jj.2,])),dc,length(jj.2)))%*%parm4) )  
    rm(Xt1,Xt2)
    }

    ## minus ties

     if (length(jj.2.tie)>0){
     Xt1=datwide.delta0[match(ID[jj.2.tie],datwide.delta0$id),Vname[1]]
     Xt2=datwide.delta0[match(ID[jj.2.tie],datwide.delta0$id),Vname[2]]
     
         ## make t1, t2 vectors
     t1=datwide.delta0[match(ID[jj.2.tie],datwide.delta0$id),Tname[1]]
     t2=datwide.delta0[match(ID[jj.2.tie],datwide.delta0$id),Tname[2]]
    

   ### when t1 t2 are vectors 
     COVV.b1<-coef.rd[1,1]+t1*coef.rd[1,2]
     COVV.b2<-coef.rd[1,1]+t2*coef.rd[1,2]
     COVV.b3<-coef.rd[1,1]+(t1+SE[i]-RE[jj.2.tie])*coef.rd[1,2]+t1*(SE[i]-RE[jj.2.tie])*coef.rd[2,2]
     COVV.b4<-coef.rd[1,1]+t2*(SE[i]-RE[jj.2.tie])*coef.rd[2,2]+(t2+SE[i]-RE[jj.2.tie])*coef.rd[1,2]
     

     SS.b.1<-coef.rd[1,1]+2*t1*coef.rd[1,2]+t1*2*coef.rd[2,2]+sigma.e^2
     SS.b.2<-coef.rd[1,1]+2*t2*coef.rd[1,2]+t2^2*coef.rd[2,2]+sigma.e^2
     SS.b.12<-coef.rd[1,1]+(t1+t2)*coef.rd[1,2]+t1*t2*coef.rd[2,2]
     
     SS.b<-array(NA,dim=c(2,2,length(jj.2.tie)))
     SS.b[1,1,]<-SS.b.1
     SS.b[1,2,]<-SS.b[2,1,]<-SS.b.12
     SS.b[2,2,]<-SS.b.2
     
     SS.b<-apply(SS.b,3,solve)  ### SS.b is a 4 by length(jj.2) matrix with each column the inverse of a matrix
 
     ### slope
    
     ### When t1 t2 are vectors , SS.b[1,], SS.b[2,], SS.b[3,], SS.b[4,] 
     ### corresponds [1,1],[2,1],[1,2],[2,2] of a 2x2 matrix
     
     gamma1=COVV.b1*SS.b[1,]+COVV.b2*SS.b[2,]
     gamma2=COVV.b1*SS.b[3,]+COVV.b2*SS.b[4,]
     gamma3=COVV.b3*SS.b[1,]+COVV.b4*SS.b[2,]
     gamma4=COVV.b3*SS.b[3,]+COVV.b4*SS.b[4,]
     
     ### Sigma1.2.b
     eta1=gamma1*COVV.b1+gamma2*COVV.b2
     eta2=gamma1*COVV.b3+gamma2*COVV.b4
     eta3=gamma3*COVV.b1+gamma4*COVV.b2
     eta4=gamma3*COVV.b3+gamma4*COVV.b4
     
      mu0.EXE<-cbind(1,0,matrix(as.numeric(WE[jj.2.tie,]),length(jj.2.tie),dw))%*%coef.fix
       mu.EXE<-cbind(1,SE[i]-RE[jj.2.tie],matrix(as.numeric(WE[jj.2.tie,]),length(jj.2.tie),dw))%*%coef.fix
       mut1.EXE<-cbind(1,t1,matrix(as.numeric(WE[jj.2.tie,]),length(jj.2.tie),dw))%*%coef.fix
       mut2.EXE<-cbind(1,t2,matrix(as.numeric(WE[jj.2.tie,]),length(jj.2.tie),dw))%*%coef.fix
       
      
      comp1<-parm1*(mu0.EXE+gamma1*Xt1+gamma2*Xt2-gamma1*mut1.EXE-gamma2*mut2.EXE)+
        parm2*(mu.EXE+gamma3*Xt1+gamma4*Xt2-gamma3*mut1.EXE-gamma4*mut2.EXE)
      
      comp2<-parm1^2*(coef.rd[1,1]-eta1)+parm2*parm1*(2*coef.rd[1,1]+2*(SE[i]-RE[jj.2.tie])*coef.rd[1,2]-eta2-eta3)+
      parm2^2*(coef.rd[1,1]+(SE[i]-RE[jj.2.tie])^2*coef.rd[2,2]+2*(SE[i]-RE[jj.2.tie])*coef.rd[1,2]-eta4)
      
    
    ss=ss-(ind.tie-1)/di*sum(weights[jj.2.tie]*exp(as.numeric(comp1)+0.5*comp2-comp12i+parm3*(SE[i]-RE[jj.2.tie]-ZE[i]*TE[i])
    +parm0*(ZE[jj.2.tie]-ZE[i])-t(as.numeric(WEc[i,])-matrix(as.numeric(t(WEc[jj.2.tie,])),dc,length(jj.2.tie)))%*%parm4) )  
    rm(Xt1,Xt2)
    }
    #####  ###### Measure of X(t1)*, X(t2)*, and X(t3)*
    
     if (length(jj.3)>0){
     Xt1=datwide.delta0[match(ID[jj.3],datwide.delta0$id),Vname[1]]
     Xt2=datwide.delta0[match(ID[jj.3],datwide.delta0$id),Vname[2]]
     Xt3=datwide.delta0[match(ID[jj.3],datwide.delta0$id),Vname[3]]
     
       ### t1, t2, t3 vector
     t1=datwide.delta0[match(ID[jj.3],datwide.delta0$id),Tname[1]]
     t2=datwide.delta0[match(ID[jj.3],datwide.delta0$id),Tname[2]]
     t3=datwide.delta0[match(ID[jj.3],datwide.delta0$id),Tname[3]]
     
   ####  t1,t2,t3 vector

     COVV.c1=coef.rd[1,1]+t1*coef.rd[1,2]
     COVV.c2=coef.rd[1,1]+t2*coef.rd[1,2]
     COVV.c3=coef.rd[1,1]+t3*coef.rd[1,2]
     COVV.c4=coef.rd[1,1]+(t1+SE[i]-RE[jj.3])*coef.rd[1,2]+t1*(SE[i]-RE[jj.3])*coef.rd[2,2]
     COVV.c5=coef.rd[1,1]+(t2+SE[i]-RE[jj.3])*coef.rd[1,2]+t2*(SE[i]-RE[jj.3])*coef.rd[2,2]
     COVV.c6=coef.rd[1,1]+(t3+SE[i]-RE[jj.3])*coef.rd[1,2]+t3*(SE[i]-RE[jj.3])*coef.rd[2,2]
     
     ##########
    
     SS.c.11<-coef.rd[1,1]+2*t1*coef.rd[1,2]+t1^2*coef.rd[2,2]+sigma.e^2
     SS.c.22<-coef.rd[1,1]+2*t2*coef.rd[1,2]+t2^2*coef.rd[2,2]+sigma.e^2
     SS.c.33<-coef.rd[1,1]+2*t3*coef.rd[1,2]+t3^2*coef.rd[2,2]+sigma.e^2
     SS.c.12=coef.rd[1,1]+(t1+t2)*coef.rd[1,2]+t1*t2*coef.rd[2,2]
     SS.c.13=coef.rd[1,1]+(t1+t3)*coef.rd[1,2]+t1*t3*coef.rd[2,2]
     SS.c.23=coef.rd[1,1]+(t2+t3)*coef.rd[1,2]+t2*t3*coef.rd[2,2]
     
     SS.c<-array(NA,dim=c(3,3,length(jj.3)))
     SS.c[1,1,]<-SS.c.11
     SS.c[2,2,]<-SS.c.22
     SS.c[3,3,]<-SS.c.33
     SS.c[1,2,]<-SS.c[2,1,]<-SS.c.12
     SS.c[1,3,]<-SS.c[3,1,]<-SS.c.13
     SS.c[2,3,]<-SS.c[3,2,]<-SS.c.23
     
    SS.c<-apply(SS.c,3,solve)
    

     #### SS.c is 9xlength(jj.3), elements in each column corresponds to [1,1],[2,1],[3,1],[1,2],[2,2],[3,2],[1,3],[2,3],[3,3] of a matrix

     delta1=COVV.c1*SS.c[1,]+COVV.c2*SS.c[2,]+COVV.c3*SS.c[3,]
     delta2=COVV.c1*SS.c[4,]+COVV.c2*SS.c[5,]+COVV.c3*SS.c[6,]
     delta3=COVV.c1*SS.c[7,]+COVV.c2*SS.c[8,]+COVV.c3*SS.c[9,]
     delta4=COVV.c4*SS.c[1,]+COVV.c5*SS.c[2,]+COVV.c6*SS.c[3,]
     delta5=COVV.c4*SS.c[4,]+COVV.c5*SS.c[5,]+COVV.c6*SS.c[6,]
     delta6=COVV.c4*SS.c[7,]+COVV.c5*SS.c[8,]+COVV.c6*SS.c[9,]

     xi1=delta1*COVV.c1+delta2*COVV.c2+delta3*COVV.c3
     xi2=delta1*COVV.c4+delta2*COVV.c5+delta3*COVV.c6
     xi3=delta4*COVV.c1+delta5*COVV.c2+delta6*COVV.c3
     xi4=delta4*COVV.c4+delta5*COVV.c5+delta6*COVV.c6
      
      mu0.EXE<-cbind(1,0,matrix(as.numeric(WE[jj.3,]),length(jj.3),dw))%*%coef.fix
       mu.EXE<-cbind(1,SE[i]-RE[jj.3],matrix(as.numeric(WE[jj.3,]),length(jj.3),dw))%*%coef.fix
       mut1.EXE<-cbind(1,t1,matrix(as.numeric(WE[jj.3,]),length(jj.3),dw))%*%coef.fix
       mut2.EXE<-cbind(1,t2,matrix(as.numeric(WE[jj.3,]),length(jj.3),dw))%*%coef.fix
       mut3.EXE<-cbind(1,t3,matrix(as.numeric(WE[jj.3,]),length(jj.3),dw))%*%coef.fix
       
      
      comp1<-parm1*(mu0.EXE+delta1*Xt1+delta2*Xt2+delta3*Xt3-delta1*mut1.EXE-delta2*mut2.EXE-delta3*mut3.EXE)+
        parm2*(mu.EXE+delta4*Xt1+delta5*Xt2+delta6*Xt3-delta4*mut1.EXE-delta5*mut2.EXE-delta6*mut3.EXE)
      
      comp2<-parm1^2*(coef.rd[1,1]-xi1)+parm2*parm1*(2*coef.rd[1,1]+2*(SE[i]-RE[jj.3])*coef.rd[1,2]-xi2-xi3)+
      parm2^2*(coef.rd[1,1]+(SE[i]-RE[jj.3])^2*coef.rd[2,2]+2*(SE[i]-RE[jj.3])*coef.rd[1,2]-xi4)
      
    
    ss=ss+sum(weights[jj.3]*exp(as.numeric(comp1)+0.5*comp2-comp12i+parm3*(SE[i]-RE[jj.3]-ZE[i]*TE[i])
    +parm0*(ZE[jj.3]-ZE[i])-t(as.numeric(WEc[i,])-matrix(as.numeric(t(WEc[jj.3,])),dc,length(jj.3)))%*%parm4) )  
   
    rm(Xt1,Xt2,Xt3)
    }

    #### Minus ties
 
     if (length(jj.3.tie)>0){
     Xt1=datwide.delta0[match(ID[jj.3.tie],datwide.delta0$id),Vname[1]]
     Xt2=datwide.delta0[match(ID[jj.3.tie],datwide.delta0$id),Vname[2]]
     Xt3=datwide.delta0[match(ID[jj.3.tie],datwide.delta0$id),Vname[3]]
     
      ### t1, t2, t3 vector
     t1=datwide.delta0[match(ID[jj.3.tie],datwide.delta0$id),Tname[1]]
     t2=datwide.delta0[match(ID[jj.3.tie],datwide.delta0$id),Tname[2]]
     t3=datwide.delta0[match(ID[jj.3.tie],datwide.delta0$id),Tname[3]]
 
     COVV.c1=coef.rd[1,1]+t1*coef.rd[1,2]
     COVV.c2=coef.rd[1,1]+t2*coef.rd[1,2]
     COVV.c3=coef.rd[1,1]+t3*coef.rd[1,2]
     COVV.c4=coef.rd[1,1]+(t1+SE[i]-RE[jj.3.tie])*coef.rd[1,2]+t1*(SE[i]-RE[jj.3.tie])*coef.rd[2,2]
     COVV.c5=coef.rd[1,1]+(t2+SE[i]-RE[jj.3.tie])*coef.rd[1,2]+t2*(SE[i]-RE[jj.3.tie])*coef.rd[2,2]
     COVV.c6=coef.rd[1,1]+(t3+SE[i]-RE[jj.3.tie])*coef.rd[1,2]+t3*(SE[i]-RE[jj.3.tie])*coef.rd[2,2]
     
     ##########
    
     SS.c.11<-coef.rd[1,1]+2*t1*coef.rd[1,2]+t1^2*coef.rd[2,2]+sigma.e^2
     SS.c.22<-coef.rd[1,1]+2*t2*coef.rd[1,2]+t2^2*coef.rd[2,2]+sigma.e^2
     SS.c.33<-coef.rd[1,1]+2*t3*coef.rd[1,2]+t3^2*coef.rd[2,2]+sigma.e^2
     SS.c.12=coef.rd[1,1]+(t1+t2)*coef.rd[1,2]+t1*t2*coef.rd[2,2]
     SS.c.13=coef.rd[1,1]+(t1+t3)*coef.rd[1,2]+t1*t3*coef.rd[2,2]
     SS.c.23=coef.rd[1,1]+(t2+t3)*coef.rd[1,2]+t2*t3*coef.rd[2,2]
     
     SS.c<-array(NA,dim=c(3,3,length(jj.3.tie)))
     SS.c[1,1,]<-SS.c.11
     SS.c[2,2,]<-SS.c.22
     SS.c[3,3,]<-SS.c.33
     SS.c[1,2,]<-SS.c[2,1,]<-SS.c.12
     SS.c[1,3,]<-SS.c[3,1,]<-SS.c.13
     SS.c[2,3,]<-SS.c[3,2,]<-SS.c.23
         
    SS.c<-apply(SS.c,3,solve)
    

     #### SS.c is 9xlength(jj.3), elements in each column corresponds to [1,1],[2,1],[3,1],[1,2],[2,2],[3,2],[1,3],[2,3],[3,3] of a matrix

     delta1=COVV.c1*SS.c[1,]+COVV.c2*SS.c[2,]+COVV.c3*SS.c[3,]
     delta2=COVV.c1*SS.c[4,]+COVV.c2*SS.c[5,]+COVV.c3*SS.c[6,]
     delta3=COVV.c1*SS.c[7,]+COVV.c2*SS.c[8,]+COVV.c3*SS.c[9,]
     delta4=COVV.c4*SS.c[1,]+COVV.c5*SS.c[2,]+COVV.c6*SS.c[3,]
     delta5=COVV.c4*SS.c[4,]+COVV.c5*SS.c[5,]+COVV.c6*SS.c[6,]
     delta6=COVV.c4*SS.c[7,]+COVV.c5*SS.c[8,]+COVV.c6*SS.c[9,]

     
     xi1=delta1*COVV.c1+delta2*COVV.c2+delta3*COVV.c3
     xi2=delta1*COVV.c4+delta2*COVV.c5+delta3*COVV.c6
     xi3=delta4*COVV.c1+delta5*COVV.c2+delta6*COVV.c3
     xi4=delta4*COVV.c4+delta5*COVV.c5+delta6*COVV.c6
      
      mu0.EXE<-cbind(1,0,matrix(as.numeric(WE[jj.3.tie,]),length(jj.3.tie),dw))%*%coef.fix
       mu.EXE<-cbind(1,SE[i]-RE[jj.3.tie],matrix(as.numeric(WE[jj.3.tie,]),length(jj.3.tie),dw))%*%coef.fix
       mut1.EXE<-cbind(1,t1,matrix(as.numeric(WE[jj.3.tie,]),length(jj.3.tie),dw))%*%coef.fix
       mut2.EXE<-cbind(1,t2,matrix(as.numeric(WE[jj.3.tie,]),length(jj.3.tie),dw))%*%coef.fix
       mut3.EXE<-cbind(1,t3,matrix(as.numeric(WE[jj.3.tie,]),length(jj.3.tie),dw))%*%coef.fix
            
      comp1<-parm1*(mu0.EXE+delta1*Xt1+delta2*Xt2+delta3*Xt3-delta1*mut1.EXE-delta2*mut2.EXE-delta3*mut3.EXE)+
        parm2*(mu.EXE+delta4*Xt1+delta5*Xt2+delta6*Xt3-delta4*mut1.EXE-delta5*mut2.EXE-delta6*mut3.EXE)
      
      comp2<-parm1^2*(coef.rd[1,1]-xi1)+parm2*parm1*(2*coef.rd[1,1]+2*(SE[i]-RE[jj.3.tie])*coef.rd[1,2]-xi2-xi3)+
      parm2^2*(coef.rd[1,1]+(SE[i]-RE[jj.3.tie])^2*coef.rd[2,2]+2*(SE[i]-RE[jj.3.tie])*coef.rd[1,2]-xi4)
      
    ss=ss-(ind.tie-1)/di*sum(weights[jj.3.tie]*exp(as.numeric(comp1)+0.5*comp2-comp12i+parm3*(SE[i]-RE[jj.3.tie]-ZE[i]*TE[i])
    +parm0*(ZE[jj.3.tie]-ZE[i])-t(as.numeric(WEc[i,])-matrix(as.numeric(t(WEc[jj.3.tie,])),dc,length(jj.3.tie)))%*%parm4) )  
   
    rm(Xt1,Xt2,Xt3)
    }

    ss.new=ss
    ### placebo
   
   ### need to use jj.Z0 and jj.Z0.tie when there are confounders in the risk model
   if (length(jj.Z0)>0)
   ss=ss+sum(weights[jj.Z0]*exp(-comp12i-parm3*(ZE[i]*TE[i])-parm0*(ZE[i])
   -t(as.numeric(WEc[i,])-matrix(as.numeric(t(WEc[jj.Z0,])),dc,length(jj.Z0)))%*%parm4))
   
   if (length(jj.Z0.tie)>0)
   ss=ss-(ind.tie-1)/di*sum(weights[jj.Z0.tie]*exp(-comp12i-parm3*(ZE[i]*TE[i])-parm0*(ZE[i])-t(as.numeric(WEc[i,])-matrix(as.numeric(t(WEc[jj.Z0.tie,])),dc,length(jj.Z0.tie)))%*%parm4))
   
    out=out+weights.bar*log(weights[i]+ss)+(weights.bar-weights[i])*(comp12i+parm3*(SE[i]-RE[i])*ZE[i]+parm0*(ZE[i])
    +matrix(WEc[i,],length(i),dc)%*%parm4)
    
    outi=c(outi,ss)
 
    }
 return(out)
 }

########

######## Fit a Cox hazard model on calendar time with coefficient for Z, Z*x(0), and W

nloglik.M1.b.Cal.1<- function(parm){

 out=0
 
 outi<-NULL
   ### coef for Z, Z*x(0), and confounders respectively
    parm0=parm[1]
    parm1=parm[2]
    #parm2=parm[3]
    #parm3=parm[4]
    parm4=parm[3]
   
   jj.4.full<-(1:length(RE))[datwide.delta0[match(ID,datwide.delta0$id),]$groupall > datwide.delta0[match(ID,datwide.delta0$id),]$group]
   jj.4.full<-as.numeric(na.omit(jj.4.full))
   out.jj.4.full<-NULL
   for (j in jj.4.full){
      
    keep=1:ll
 
    VVj.11<-coef.rd[1,1]
    VVj.12<-VVj.21<-coef.rd[1,1]+(SE[keep]-RE[j])*coef.rd[1,2]
    VVj.22<-coef.rd[1,1]+2*(SE[keep]-RE[j])*coef.rd[1,2]+(SE[keep]-RE[j])^2*coef.rd[2,2]
    
    muj.1=cbind(1,0,matrix(WE[j,],length(j),dw))%*%coef.fix
    muj.2=cbind(1,SE[keep]-RE[j],matrix(WE[j,],length(j),dw))%*%coef.fix
    
   
    oo.j<-id.long==ID[j]
    Xj<-X.long[oo.j]
    Tj<-T.long[oo.j]
    
    
    VVj.c<-cbind(1,Tj)%*%coef.rd%*%rbind(1,Tj)+diag(rep(sigma.e^2,length(Tj)),nrow=length(Tj),ncol=length(Tj))               
   
    COVVj.1<-coef.rd[1,1]+coef.rd[1,2]*Tj
    COVVj.2<-coef.rd[1,1]+coef.rd[1,2]*outer(SE[keep]-RE[j],Tj,'+')+coef.rd[2,2]*outer(SE[keep]-RE[j],Tj,'*')
    COVVj.1=matrix(COVVj.1,1,length(COVVj.1))
    if (length(Tj)==1) COVVj.2=matrix(COVVj.2,length(COVVj.2),1)
    
    muj.c<-cbind(1,Tj,matrix(WE[j,],length(j),dw))%*%coef.fix
        
    comp1j<-cbind(as.numeric(muj.1),muj.2)+cbind(as.numeric(COVVj.1%*%solve(VVj.c)%*%(Xj-muj.c)),COVVj.2%*%solve(VVj.c)%*%(Xj-muj.c))
    
     VVjminus.11<-COVVj.1%*%solve(VVj.c)%*%(t(COVVj.1))
     VVjminus.12<-COVVj.2%*%solve(VVj.c)%*%(t(COVVj.1))
     VVjminus.22<-diag(COVVj.2%*%solve(VVj.c)%*%(t(COVVj.2)))
     
     comp2j.11<-VVj.11-VVjminus.11
     comp2j.12<-VVj.12-VVjminus.12
     comp2j.22<-VVj.22-VVjminus.22
   
     out.jj.4.full<-cbind(out.jj.4.full,as.numeric(c(parm1,parm2)%*%t(comp1j))+parm1^2*as.numeric(comp2j.11)/2+parm1*parm2*comp2j.12+parm2^2*comp2j.22/2+parm3*(SE[keep]-RE[j])
     +parm0*ZE[j]+matrix(WEc[j,],length(j),dc)%*%parm4)
     
    }
        
    S.tie=SE.u[1]  ### unique SE value
    ind.tie=0       ### position within the unique value
    
    for (i in 1:ll){
    
    ## track #of ties
    
        di<-sum(SE==SE[i] & delta==1)
        
        if (SE[i]==S.tie){  ## in the same tie as previous value
          ind.tie=ind.tie+1
        } else { ## go to a new SE value
           S.tie=SE[i]
           ind.tie=1
        }
     ### compute average weight among ties
        weights.bar=mean(weights[SE==SE[i] & delta==1])
     
    ss=0
    comp12i=0 ### placebo contribute nothing
    if (ZE[i]==1) { 
    VVi<-cbind(1,c(0,TE[i]))%*%coef.rd%*%rbind(1,c(0,TE[i]))
    mui<-rbind(c(1,0,matrix(WE[i,],length(i),dw)),c(1,TE[i],matrix(WE[i,],length(i),dw)))%*%coef.fix
   
    oo.i<-id.long==ID[i]
        
       #### (i) if case has at least one X(Ti)* measured, i.e. case ID[i] belong to id.long
        if (sum(oo.i)>0){ 
        Xi<-X.long[oo.i]
        Ti<-T.long[oo.i]
        
        VVi.c<-cbind(1,Ti)%*%coef.rd%*%rbind(1,Ti)+diag(rep(sigma.e^2,length(Ti)),nrow=length(Ti),ncol=length(Ti))               
        COVVi<-cbind(1,c(0,TE[i]))%*%coef.rd%*%rbind(1,Ti)
        mui.c<-cbind(1,Ti,matrix(WE[i,],length(Ti),dw,byrow=TRUE))%*%coef.fix
        
        comp1i<-mui+COVVi%*%solve(VVi.c)%*%(Xi-mui.c)      ## E(X(0),X(t)|X(Ti)*,W)
        comp2i<-VVi-COVVi%*%solve(VVi.c)%*%t(COVVi)        ## Var(X(0),X(t)|X(Ti)*,W)
        
        } else {
        #### (ii) if case does not have X(Ti)* measure, i.e. case ID[i] not in id.long
        
       comp1i<-mui            
       comp2i<-VVi            
       } 
 
     comp12i<-as.numeric(t(c(parm1,parm2))%*%comp1i+t(c(parm1,parm2))%*%comp2i%*%c(parm1,parm2)/2)
    
    } 
    
    jj.r<-((i+1):(length(RE)))[RE[(i+1):length(RE)]<SE[i] & SE[(i+1):length(RE)]>=SE[i] & ID[(i+1):length(RE)]%in%id.long]
    jj.l<-NULL
    if (i>1) jj.l<-(1:(i-1))[RE[1:(i-1)]<SE[i] & SE[1:(i-1)]>=SE[i] & ID[1:(i-1)]%in%id.long]
    jj<-c(jj.l,jj.r) ## among risk set  
    
    jj.tie<-(1:(length(RE)))[RE[1:length(RE)]<SE[i] & ID[1:length(RE)]%in%id.long & SE[1:length(SE)]==SE[i] & delta[1:length(delta)]==1]
   
    jj.W.r<-((i+1):(length(RE)))[RE[(i+1):length(RE)]<SE[i] & ZE[(i+1):length(ZE)]==1 & SE[(i+1):length(RE)]>=SE[i] & !ID[(i+1):length(RE)]%in%id.long]
    jj.W.l<-NULL
    if (i>1) jj.W.l<-(1:(i-1))[RE[1:(i-1)]<SE[i] & ZE[1:(i-1)]==1 & !ID[1:(i-1)]%in%id.long & SE[1:(i-1)]>=SE[i]]
    jj.W<-c(jj.W.l,jj.W.r)    

    
    jj.W.tie<-(1:(length(RE)))[RE[1:length(RE)]<SE[i] & ZE[1:length(ZE)]==1 & SE[1:length(RE)]==SE[i] & !ID[1:length(RE)]%in%id.long  & delta[1:length(delta)]==1]
 
  
    jj.Z0.r<-((i+1):length(RE))[RE[(i+1):length(RE)]<SE[i] & ZE[(i+1):length(ZE)]==0 & SE[(i+1):length(SE)]>=SE[i]]  
    jj.Z0.l<-NULL
    if (i>1) jj.Z0.l<-(1:(i-1))[RE[1:(i-1)]<SE[i] & ZE[1:(i-1)]==0 & SE[1:(i-1)]>=SE[i]]
    jj.Z0=c(jj.Z0.l,jj.Z0.r)

    jj.Z0.tie<-(1:length(RE))[RE[1:length(RE)]<SE[i] & ZE[1:length(ZE)]==0 & SE[1:length(SE)]==SE[i] & delta[1:length(delta)]==1]  
  
    jj.1.r<-((i+1):(length(RE)))[RE[(i+1):length(RE)]<SE[i] & datwide.delta0[match(ID[(i+1):length(RE)],datwide.delta0$id),]$group==1 &
       datwide.delta0[match(ID[(i+1):length(RE)],datwide.delta0$id),]$groupall==1 & SE[(i+1):length(RE)]>=SE[i]]
    jj.1.l<-NULL
    if (i>1) jj.1.l<-(1:(i-1))[RE[1:(i-1)]<SE[i] & datwide.delta0[match(ID[1:(i-1)],datwide.delta0$id),]$group==1 &
       datwide.delta0[match(ID[1:(i-1)],datwide.delta0$id),]$groupall==1 & SE[1:(i-1)]>=SE[i]]
    jj.1<-c(jj.1.l,jj.1.r) ## among risk set

    
    jj.1.tie<-(1:(length(RE)))[RE[1:length(RE)]<SE[i] & datwide.delta0[match(ID[1:length(RE)],datwide.delta0$id),]$group==1 &
       datwide.delta0[match(ID[1:length(RE)],datwide.delta0$id),]$groupall==1 & SE[1:length(SE)]==SE[i] & delta[1:length(delta)]==1]
   
   
    #######
 
    jj.2.r<-((i+1):(length(RE)))[RE[(i+1):length(RE)]<SE[i] & datwide.delta0[match(ID[(i+1):length(RE)],datwide.delta0$id),]$group==2 & 
       datwide.delta0[match(ID[(i+1):length(RE)],datwide.delta0$id),]$groupall==2 & SE[(i+1):length(SE)]>=SE[i]]
    jj.2.l<-NULL
    if (i>1)   jj.2.l<-(1:(i-1))[RE[1:(i-1)]<SE[i] & datwide.delta0[match(ID[1:(i-1)],datwide.delta0$id),]$group==2 & 
       datwide.delta0[match(ID[1:(i-1)],datwide.delta0$id),]$groupall==2 & SE[1:(i-1)]>=SE[i]]
    jj.2<-c(jj.2.l,jj.2.r)

    jj.2.tie<-(1:(length(RE)))[RE[1:length(RE)]<SE[i] & datwide.delta0[match(ID[1:length(RE)],datwide.delta0$id),]$group==2 & 
       datwide.delta0[match(ID[1:length(RE)],datwide.delta0$id),]$groupall==2 & SE[1:length(SE)]==SE[i] & delta[1:length(delta)]==1]
   

    ###
    jj.3.r<-((i+1):(length(RE)))[RE[(i+1):length(RE)]<SE[i] & datwide.delta0[match(ID[(i+1):length(RE)],datwide.delta0$id),]$group==3 & 
       datwide.delta0[match(ID[(i+1):length(RE)],datwide.delta0$id),]$groupall==3 & SE[(i+1):length(SE)]>=SE[i]]
    jj.3.l<-NULL
    if (i>1) jj.3.l<-(1:(i-1))[RE[1:(i-1)]<SE[i] & datwide.delta0[match(ID[1:(i-1)],datwide.delta0$id),]$group==3 & 
       datwide.delta0[match(ID[1:(i-1)],datwide.delta0$id),]$groupall==3 & SE[1:(i-1)]>=SE[i]]
    jj.3<-c(jj.3.l,jj.3.r)
   

    jj.3.tie<-(1:(length(RE)))[RE[1:length(RE)]<SE[i] & datwide.delta0[match(ID[1:length(RE)],datwide.delta0$id),]$group==3 & 
       datwide.delta0[match(ID[1:length(RE)],datwide.delta0$id),]$groupall==3 & SE[1:length(SE)]==SE[i] & delta[1:length(delta)]==1]
  
     jj.1<-as.numeric(na.omit(jj.1))
     jj.2<-as.numeric(na.omit(jj.2))
     jj.3<-as.numeric(na.omit(jj.3))
   
       jj.1.tie<-as.numeric(na.omit(jj.1.tie))
     jj.2.tie<-as.numeric(na.omit(jj.2.tie))
     jj.3.tie<-as.numeric(na.omit(jj.3.tie))
     
    ss=ss+sum((RE[jj.4.full]<SE[i])*(SE[jj.4.full]>=SE[i])*(!jj.4.full==i)*weights[jj.4.full]*exp(out.jj.4.full[i,]-comp12i-parm3*(ZE[i]*TE[i])
    -parm0*ZE[i]-matrix(WEc[i,],length(i),dc)%*%parm4))  ## risk set from jj.4.full
     
    ss=ss-(ind.tie-1)/di*sum((RE[jj.4.full]<SE[i])*(SE[jj.4.full]==SE[i])*(!jj.4.full==i)*weights[jj.4.full]*exp(out.jj.4.full[i,]-comp12i-parm3*(ZE[i]*TE[i])
    -parm0*ZE[i]-matrix(WEc[i,],length(i),dc)%*%parm4))  ## tie from jj.4.full
     

   if (length(jj.W)>0){
      ### No X measure 
    mu0.EXE<-cbind(1,0,matrix(as.numeric(WE[jj.W,]),length(jj.W),dw))%*%coef.fix
    mu.EXE<-cbind(1,SE[i]-RE[jj.W],matrix(as.numeric(WE[jj.W,]),length(jj.W),dw))%*%coef.fix
  
   comp1.b<-parm1*mu0.EXE+parm2*mu.EXE
       
   comp2.b<-parm1^2*coef.rd[1,1]+2*parm1*parm2*(coef.rd[1,1]+(SE[i]-RE[jj.W])*coef.rd[1,2])+
   parm2^2*(coef.rd[1,1]+2*(SE[i]-RE[jj.W])*coef.rd[1,2]+(SE[i]-RE[jj.W])^2*coef.rd[2,2])
   
   ss=ss+sum(weights[jj.W]*exp(as.numeric(comp1.b)+0.5*comp2.b-comp12i+parm3*(SE[i]-RE[jj.W]-ZE[i]*TE[i])
   +parm0*(ZE[jj.W]-ZE[i])-t(as.numeric(WEc[i,])-matrix(as.numeric(t(WEc[jj.W,])),dc,length(jj.W)))%*%parm4) )  
    
   } 
   ### minus ties
   if (length(jj.W.tie)>0){
   
        mu0.EXE<-cbind(1,0,matrix(as.numeric(WE[jj.W.tie,]),length(jj.W.tie),dw))%*%coef.fix
        mu.EXE<-cbind(1,SE[i]-RE[jj.W.tie],matrix(as.numeric(WE[jj.W.tie,]),length(jj.W.tie),dw))%*%coef.fix
      
       comp1.b<-parm1*mu0.EXE+parm2*mu.EXE
           
       comp2.b<-parm1^2*coef.rd[1,1]+2*parm1*parm2*(coef.rd[1,1]+(SE[i]-RE[jj.W.tie])*coef.rd[1,2])+
       parm2^2*(coef.rd[1,1]+2*(SE[i]-RE[jj.W.tie])*coef.rd[1,2]+(SE[i]-RE[jj.W.tie])^2*coef.rd[2,2])
       
       ss=ss-(ind.tie-1)/di*sum(weights[jj.W.tie]*exp(as.numeric(comp1.b)+0.5*comp2.b-comp12i+parm3*(SE[i]-RE[jj.W.tie]-ZE[i]*TE[i])
       +parm0*(ZE[jj.W.tie]-ZE[i])-t(as.numeric(WEc[i,])-matrix(as.numeric(t(WEc[jj.W.tie,])),dc,length(jj.W.tie)))%*%parm4) ) 
   }      

 
   ### Measure of X(0)* only
   
    #account for different t levels
   
   if (length(jj.1)>0){
   Xt1=datwide.delta0[match(ID[jj.1],datwide.delta0$id),Vname[1]]
    t1=datwide.delta0[match(ID[jj.1],datwide.delta0$id),Tname[1]]
 
    ## make t1 vector
   COVV.a<-matrix(c(coef.rd[1,1]+t1*coef.rd[1,2],coef.rd[1,1]+(t1+SE[i]-RE[jj.1])*coef.rd[1,2]+t1*(SE[i]-RE[jj.1])*coef.rd[2,2]),length(jj.1),2,byrow=F)
   SS.a<-sapply(coef.rd[1,1]+2*t1*coef.rd[1,2]+t1^2*coef.rd[2,2]+sigma.e^2,solve)
   
   slope.a<-COVV.a*as.numeric(SS.a)
   Sigma1.2.a=outer(t(slope.a),(COVV.a),'*'); 
   kappa1=diag(as.matrix(Sigma1.2.a[1,,,1])); kappa2=diag(as.matrix(Sigma1.2.a[1,,,2])); 
   kappa3=diag(as.matrix(Sigma1.2.a[2,,,1])); kappa4=diag(as.matrix(Sigma1.2.a[2,,,2]))
   
     mu0.EXE<-cbind(1,0,matrix(as.numeric(WE[jj.1,]),length(jj.1),dw))%*%coef.fix
      mu.EXE<-cbind(1,SE[i]-RE[jj.1],matrix(as.numeric(WE[jj.1,]),length(jj.1),dw))%*%coef.fix
  
       comp1<-parm1*(mu0.EXE+slope.a[,1]*(Xt1-mu0.EXE))+
       parm2*(mu.EXE+slope.a[,2]*(Xt1-mu0.EXE))
       
       comp2<-parm1^2*(coef.rd[1,1]-kappa1)+
       parm2*parm1*(2*coef.rd[1,1]+2*(SE[i]-RE[jj.1])*coef.rd[1,2]-kappa2-kappa3)+
       parm2^2*(coef.rd[1,1]+(SE[i]-RE[jj.1])^2*coef.rd[2,2]+2*(SE[i]-RE[jj.1])*coef.rd[1,2]-kappa4)
    
    ss=ss+sum(weights[jj.1]*exp(as.numeric(comp1)+0.5*comp2-comp12i+parm3*(SE[i]-RE[jj.1]-ZE[i]*TE[i])
    +parm0*(ZE[jj.1]-ZE[i])-t(as.numeric(WEc[i,])-matrix(as.numeric(t(WEc[jj.1,])),dc,length(jj.1)))%*%parm4) )  
   
   rm(Xt1)
   }

   ## minus ties

   if (length(jj.1.tie)>0){
   Xt1=datwide.delta0[match(ID[jj.1.tie],datwide.delta0$id),Vname[1]]
   
   t1=datwide.delta0[match(ID[jj.1.tie],datwide.delta0$id),Tname[1]]
   
    # make t1 vector
   COVV.a<-matrix(c(coef.rd[1,1]+t1*coef.rd[1,2],coef.rd[1,1]+(t1+SE[i]-RE[jj.1.tie])*coef.rd[1,2]+t1*(SE[i]-RE[jj.1.tie])*coef.rd[2,2]),length(jj.1.tie),2,byrow=F)
   SS.a<-sapply(coef.rd[1,1]+2*t1*coef.rd[1,2]+t1^2*coef.rd[2,2]+sigma.e^2,solve)
 
   
   slope.a<-COVV.a*as.numeric(SS.a)
   Sigma1.2.a=outer(t(slope.a),(COVV.a),'*'); 
   kappa1=diag(as.matrix(Sigma1.2.a[1,,,1])); kappa2=diag(as.matrix(Sigma1.2.a[1,,,2])); 
   kappa3=diag(as.matrix(Sigma1.2.a[2,,,1])); kappa4=diag(as.matrix(Sigma1.2.a[2,,,2]))
   
     mu0.EXE<-cbind(1,0,matrix(as.numeric(WE[jj.1.tie,]),length(jj.1.tie),dw))%*%coef.fix
      mu.EXE<-cbind(1,SE[i]-RE[jj.1.tie],matrix(as.numeric(WE[jj.1.tie,]),length(jj.1.tie),dw))%*%coef.fix
  
       comp1<-parm1*(mu0.EXE+slope.a[,1]*(Xt1-mu0.EXE))+
       parm2*(mu.EXE+slope.a[,2]*(Xt1-mu0.EXE))
       
       comp2<-parm1^2*(coef.rd[1,1]-kappa1)+
       parm2*parm1*(2*coef.rd[1,1]+2*(SE[i]-RE[jj.1.tie])*coef.rd[1,2]-kappa2-kappa3)+
       parm2^2*(coef.rd[1,1]+(SE[i]-RE[jj.1.tie])^2*coef.rd[2,2]+2*(SE[i]-RE[jj.1.tie])*coef.rd[1,2]-kappa4)
    
    ss=ss-(ind.tie-1)/di*sum(weights[jj.1.tie]*exp(as.numeric(comp1)+0.5*comp2-comp12i+parm3*(SE[i]-RE[jj.1.tie]-ZE[i]*TE[i])
    +parm0*(ZE[jj.1.tie]-ZE[i])-t(as.numeric(WEc[i,])-matrix(as.numeric(t(WEc[jj.1.tie,])),dc,length(jj.1.tie)))%*%parm4 ) )  
   
   rm(Xt1)
   }



   ###### Measure of X(0)* and X(t2)* 
     if (length(jj.2)>0){
     Xt1=datwide.delta0[match(ID[jj.2],datwide.delta0$id),Vname[1]]
     Xt2=datwide.delta0[match(ID[jj.2],datwide.delta0$id),Vname[2]]
     
       ## make t1, t2 vectors
      t1=datwide.delta0[match(ID[jj.2],datwide.delta0$id),Tname[1]]
      t2=datwide.delta0[match(ID[jj.2],datwide.delta0$id),Tname[2]]
         
     COVV.b1<-coef.rd[1,1]+t1*coef.rd[1,2]
     COVV.b2<-coef.rd[1,1]+t2*coef.rd[1,2]
     COVV.b3<-coef.rd[1,1]+(t1+SE[i]-RE[jj.2])*coef.rd[1,2]+t1*(SE[i]-RE[jj.2])*coef.rd[2,2]
     COVV.b4<-coef.rd[1,1]+t2*(SE[i]-RE[jj.2])*coef.rd[2,2]+(t2+SE[i]-RE[jj.2])*coef.rd[1,2]
     

     SS.b.1<-coef.rd[1,1]+2*t1*coef.rd[1,2]+t1*2*coef.rd[2,2]+sigma.e^2
     SS.b.2<-coef.rd[1,1]+2*t2*coef.rd[1,2]+t2^2*coef.rd[2,2]+sigma.e^2
     SS.b.12<-coef.rd[1,1]+(t1+t2)*coef.rd[1,2]+t1*t2*coef.rd[2,2]
     
     SS.b<-array(NA,dim=c(2,2,length(jj.2)))
     SS.b[1,1,]<-SS.b.1
     SS.b[1,2,]<-SS.b[2,1,]<-SS.b.12
     SS.b[2,2,]<-SS.b.2
     
     SS.b<-apply(SS.b,3,solve)  ### SS.b is a 4 by length(jj.2) matrix with each column the inverse of a matrix
 
     ### slope
    
     ### When t1 t2 are vectors, SS.b[1,], SS.b[2,], SS.b[3,], SS.b[4,] 
     ### corresponds [1,1],[2,1],[1,2],[2,2] of a 2x2 matrix
     
     gamma1=COVV.b1*SS.b[1,]+COVV.b2*SS.b[2,]
     gamma2=COVV.b1*SS.b[3,]+COVV.b2*SS.b[4,]
     gamma3=COVV.b3*SS.b[1,]+COVV.b4*SS.b[2,]
     gamma4=COVV.b3*SS.b[3,]+COVV.b4*SS.b[4,]
     

     ### Sigma1.2.b
     eta1=gamma1*COVV.b1+gamma2*COVV.b2
     eta2=gamma1*COVV.b3+gamma2*COVV.b4
     eta3=gamma3*COVV.b1+gamma4*COVV.b2
     eta4=gamma3*COVV.b3+gamma4*COVV.b4
     
      mu0.EXE<-cbind(1,0,matrix(as.numeric(WE[jj.2,]),length(jj.2),dw))%*%coef.fix
       mu.EXE<-cbind(1,SE[i]-RE[jj.2],matrix(as.numeric(WE[jj.2,]),length(jj.2),dw))%*%coef.fix
       mut1.EXE<-cbind(1,t1,matrix(as.numeric(WE[jj.2,]),length(jj.2),dw))%*%coef.fix
       mut2.EXE<-cbind(1,t2,matrix(as.numeric(WE[jj.2,]),length(jj.2),dw))%*%coef.fix
       
      
      comp1<-parm1*(mu0.EXE+gamma1*Xt1+gamma2*Xt2-gamma1*mut1.EXE-gamma2*mut2.EXE)+
        parm2*(mu.EXE+gamma3*Xt1+gamma4*Xt2-gamma3*mut1.EXE-gamma4*mut2.EXE)
      
      comp2<-parm1^2*(coef.rd[1,1]-eta1)+parm2*parm1*(2*coef.rd[1,1]+2*(SE[i]-RE[jj.2])*coef.rd[1,2]-eta2-eta3)+
      parm2^2*(coef.rd[1,1]+(SE[i]-RE[jj.2])^2*coef.rd[2,2]+2*(SE[i]-RE[jj.2])*coef.rd[1,2]-eta4)
      
    
    ss=ss+sum(weights[jj.2]*exp(as.numeric(comp1)+0.5*comp2-comp12i+parm3*(SE[i]-RE[jj.2]-ZE[i]*TE[i])
    +parm0*(ZE[jj.2]-ZE[i])-t(as.numeric(WEc[i,])-matrix(as.numeric(t(WEc[jj.2,])),dc,length(jj.2)))%*%parm4) )  
    rm(Xt1,Xt2)
    }

    ## minus ties

     if (length(jj.2.tie)>0){
     Xt1=datwide.delta0[match(ID[jj.2.tie],datwide.delta0$id),Vname[1]]
     Xt2=datwide.delta0[match(ID[jj.2.tie],datwide.delta0$id),Vname[2]]
     
         ## make t1, t2 vectors
     t1=datwide.delta0[match(ID[jj.2.tie],datwide.delta0$id),Tname[1]]
     t2=datwide.delta0[match(ID[jj.2.tie],datwide.delta0$id),Tname[2]]
    

   ### when t1 t2 are vectors 
     COVV.b1<-coef.rd[1,1]+t1*coef.rd[1,2]
     COVV.b2<-coef.rd[1,1]+t2*coef.rd[1,2]
     COVV.b3<-coef.rd[1,1]+(t1+SE[i]-RE[jj.2.tie])*coef.rd[1,2]+t1*(SE[i]-RE[jj.2.tie])*coef.rd[2,2]
     COVV.b4<-coef.rd[1,1]+t2*(SE[i]-RE[jj.2.tie])*coef.rd[2,2]+(t2+SE[i]-RE[jj.2.tie])*coef.rd[1,2]
     

     SS.b.1<-coef.rd[1,1]+2*t1*coef.rd[1,2]+t1*2*coef.rd[2,2]+sigma.e^2
     SS.b.2<-coef.rd[1,1]+2*t2*coef.rd[1,2]+t2^2*coef.rd[2,2]+sigma.e^2
     SS.b.12<-coef.rd[1,1]+(t1+t2)*coef.rd[1,2]+t1*t2*coef.rd[2,2]
     
     SS.b<-array(NA,dim=c(2,2,length(jj.2.tie)))
     SS.b[1,1,]<-SS.b.1
     SS.b[1,2,]<-SS.b[2,1,]<-SS.b.12
     SS.b[2,2,]<-SS.b.2
     
     SS.b<-apply(SS.b,3,solve)  ### SS.b is a 4 by length(jj.2) matrix with each column the inverse of a matrix
 
     ### slope
    
     ### When t1 t2 are vectors , SS.b[1,], SS.b[2,], SS.b[3,], SS.b[4,] 
     ### corresponds [1,1],[2,1],[1,2],[2,2] of a 2x2 matrix
     
     gamma1=COVV.b1*SS.b[1,]+COVV.b2*SS.b[2,]
     gamma2=COVV.b1*SS.b[3,]+COVV.b2*SS.b[4,]
     gamma3=COVV.b3*SS.b[1,]+COVV.b4*SS.b[2,]
     gamma4=COVV.b3*SS.b[3,]+COVV.b4*SS.b[4,]
     
     ### Sigma1.2.b
     eta1=gamma1*COVV.b1+gamma2*COVV.b2
     eta2=gamma1*COVV.b3+gamma2*COVV.b4
     eta3=gamma3*COVV.b1+gamma4*COVV.b2
     eta4=gamma3*COVV.b3+gamma4*COVV.b4
     
      mu0.EXE<-cbind(1,0,matrix(as.numeric(WE[jj.2.tie,]),length(jj.2.tie),dw))%*%coef.fix
       mu.EXE<-cbind(1,SE[i]-RE[jj.2.tie],matrix(as.numeric(WE[jj.2.tie,]),length(jj.2.tie),dw))%*%coef.fix
       mut1.EXE<-cbind(1,t1,matrix(as.numeric(WE[jj.2.tie,]),length(jj.2.tie),dw))%*%coef.fix
       mut2.EXE<-cbind(1,t2,matrix(as.numeric(WE[jj.2.tie,]),length(jj.2.tie),dw))%*%coef.fix
       
      
      comp1<-parm1*(mu0.EXE+gamma1*Xt1+gamma2*Xt2-gamma1*mut1.EXE-gamma2*mut2.EXE)+
        parm2*(mu.EXE+gamma3*Xt1+gamma4*Xt2-gamma3*mut1.EXE-gamma4*mut2.EXE)
      
      comp2<-parm1^2*(coef.rd[1,1]-eta1)+parm2*parm1*(2*coef.rd[1,1]+2*(SE[i]-RE[jj.2.tie])*coef.rd[1,2]-eta2-eta3)+
      parm2^2*(coef.rd[1,1]+(SE[i]-RE[jj.2.tie])^2*coef.rd[2,2]+2*(SE[i]-RE[jj.2.tie])*coef.rd[1,2]-eta4)
      
    
    ss=ss-(ind.tie-1)/di*sum(weights[jj.2.tie]*exp(as.numeric(comp1)+0.5*comp2-comp12i+parm3*(SE[i]-RE[jj.2.tie]-ZE[i]*TE[i])
    +parm0*(ZE[jj.2.tie]-ZE[i])-t(as.numeric(WEc[i,])-matrix(as.numeric(t(WEc[jj.2.tie,])),dc,length(jj.2.tie)))%*%parm4) )  
    rm(Xt1,Xt2)
    }
    #####  ###### Measure of X(t1)*, X(t2)*, and X(t3)*
    
     if (length(jj.3)>0){
     Xt1=datwide.delta0[match(ID[jj.3],datwide.delta0$id),Vname[1]]
     Xt2=datwide.delta0[match(ID[jj.3],datwide.delta0$id),Vname[2]]
     Xt3=datwide.delta0[match(ID[jj.3],datwide.delta0$id),Vname[3]]
     
       ### t1, t2, t3 vector
     t1=datwide.delta0[match(ID[jj.3],datwide.delta0$id),Tname[1]]
     t2=datwide.delta0[match(ID[jj.3],datwide.delta0$id),Tname[2]]
     t3=datwide.delta0[match(ID[jj.3],datwide.delta0$id),Tname[3]]
     
   ####  t1,t2,t3 vector

     COVV.c1=coef.rd[1,1]+t1*coef.rd[1,2]
     COVV.c2=coef.rd[1,1]+t2*coef.rd[1,2]
     COVV.c3=coef.rd[1,1]+t3*coef.rd[1,2]
     COVV.c4=coef.rd[1,1]+(t1+SE[i]-RE[jj.3])*coef.rd[1,2]+t1*(SE[i]-RE[jj.3])*coef.rd[2,2]
     COVV.c5=coef.rd[1,1]+(t2+SE[i]-RE[jj.3])*coef.rd[1,2]+t2*(SE[i]-RE[jj.3])*coef.rd[2,2]
     COVV.c6=coef.rd[1,1]+(t3+SE[i]-RE[jj.3])*coef.rd[1,2]+t3*(SE[i]-RE[jj.3])*coef.rd[2,2]
     
     ##########
    
     SS.c.11<-coef.rd[1,1]+2*t1*coef.rd[1,2]+t1^2*coef.rd[2,2]+sigma.e^2
     SS.c.22<-coef.rd[1,1]+2*t2*coef.rd[1,2]+t2^2*coef.rd[2,2]+sigma.e^2
     SS.c.33<-coef.rd[1,1]+2*t3*coef.rd[1,2]+t3^2*coef.rd[2,2]+sigma.e^2
     SS.c.12=coef.rd[1,1]+(t1+t2)*coef.rd[1,2]+t1*t2*coef.rd[2,2]
     SS.c.13=coef.rd[1,1]+(t1+t3)*coef.rd[1,2]+t1*t3*coef.rd[2,2]
     SS.c.23=coef.rd[1,1]+(t2+t3)*coef.rd[1,2]+t2*t3*coef.rd[2,2]
     
     SS.c<-array(NA,dim=c(3,3,length(jj.3)))
     SS.c[1,1,]<-SS.c.11
     SS.c[2,2,]<-SS.c.22
     SS.c[3,3,]<-SS.c.33
     SS.c[1,2,]<-SS.c[2,1,]<-SS.c.12
     SS.c[1,3,]<-SS.c[3,1,]<-SS.c.13
     SS.c[2,3,]<-SS.c[3,2,]<-SS.c.23
     
    SS.c<-apply(SS.c,3,solve)
    

     #### SS.c is 9xlength(jj.3), elements in each column corresponds to [1,1],[2,1],[3,1],[1,2],[2,2],[3,2],[1,3],[2,3],[3,3] of a matrix

     delta1=COVV.c1*SS.c[1,]+COVV.c2*SS.c[2,]+COVV.c3*SS.c[3,]
     delta2=COVV.c1*SS.c[4,]+COVV.c2*SS.c[5,]+COVV.c3*SS.c[6,]
     delta3=COVV.c1*SS.c[7,]+COVV.c2*SS.c[8,]+COVV.c3*SS.c[9,]
     delta4=COVV.c4*SS.c[1,]+COVV.c5*SS.c[2,]+COVV.c6*SS.c[3,]
     delta5=COVV.c4*SS.c[4,]+COVV.c5*SS.c[5,]+COVV.c6*SS.c[6,]
     delta6=COVV.c4*SS.c[7,]+COVV.c5*SS.c[8,]+COVV.c6*SS.c[9,]

     xi1=delta1*COVV.c1+delta2*COVV.c2+delta3*COVV.c3
     xi2=delta1*COVV.c4+delta2*COVV.c5+delta3*COVV.c6
     xi3=delta4*COVV.c1+delta5*COVV.c2+delta6*COVV.c3
     xi4=delta4*COVV.c4+delta5*COVV.c5+delta6*COVV.c6
      
      mu0.EXE<-cbind(1,0,matrix(as.numeric(WE[jj.3,]),length(jj.3),dw))%*%coef.fix
       mu.EXE<-cbind(1,SE[i]-RE[jj.3],matrix(as.numeric(WE[jj.3,]),length(jj.3),dw))%*%coef.fix
       mut1.EXE<-cbind(1,t1,matrix(as.numeric(WE[jj.3,]),length(jj.3),dw))%*%coef.fix
       mut2.EXE<-cbind(1,t2,matrix(as.numeric(WE[jj.3,]),length(jj.3),dw))%*%coef.fix
       mut3.EXE<-cbind(1,t3,matrix(as.numeric(WE[jj.3,]),length(jj.3),dw))%*%coef.fix
       
      
      comp1<-parm1*(mu0.EXE+delta1*Xt1+delta2*Xt2+delta3*Xt3-delta1*mut1.EXE-delta2*mut2.EXE-delta3*mut3.EXE)+
        parm2*(mu.EXE+delta4*Xt1+delta5*Xt2+delta6*Xt3-delta4*mut1.EXE-delta5*mut2.EXE-delta6*mut3.EXE)
      
      comp2<-parm1^2*(coef.rd[1,1]-xi1)+parm2*parm1*(2*coef.rd[1,1]+2*(SE[i]-RE[jj.3])*coef.rd[1,2]-xi2-xi3)+
      parm2^2*(coef.rd[1,1]+(SE[i]-RE[jj.3])^2*coef.rd[2,2]+2*(SE[i]-RE[jj.3])*coef.rd[1,2]-xi4)
      
    
    ss=ss+sum(weights[jj.3]*exp(as.numeric(comp1)+0.5*comp2-comp12i+parm3*(SE[i]-RE[jj.3]-ZE[i]*TE[i])
    +parm0*(ZE[jj.3]-ZE[i])-t(as.numeric(WEc[i,])-matrix(as.numeric(t(WEc[jj.3,])),dc,length(jj.3)))%*%parm4) )  
   
    rm(Xt1,Xt2,Xt3)
    }

    #### Minus ties
 
     if (length(jj.3.tie)>0){
     Xt1=datwide.delta0[match(ID[jj.3.tie],datwide.delta0$id),Vname[1]]
     Xt2=datwide.delta0[match(ID[jj.3.tie],datwide.delta0$id),Vname[2]]
     Xt3=datwide.delta0[match(ID[jj.3.tie],datwide.delta0$id),Vname[3]]
     
      ### t1, t2, t3 vector
     t1=datwide.delta0[match(ID[jj.3.tie],datwide.delta0$id),Tname[1]]
     t2=datwide.delta0[match(ID[jj.3.tie],datwide.delta0$id),Tname[2]]
     t3=datwide.delta0[match(ID[jj.3.tie],datwide.delta0$id),Tname[3]]
 
     COVV.c1=coef.rd[1,1]+t1*coef.rd[1,2]
     COVV.c2=coef.rd[1,1]+t2*coef.rd[1,2]
     COVV.c3=coef.rd[1,1]+t3*coef.rd[1,2]
     COVV.c4=coef.rd[1,1]+(t1+SE[i]-RE[jj.3.tie])*coef.rd[1,2]+t1*(SE[i]-RE[jj.3.tie])*coef.rd[2,2]
     COVV.c5=coef.rd[1,1]+(t2+SE[i]-RE[jj.3.tie])*coef.rd[1,2]+t2*(SE[i]-RE[jj.3.tie])*coef.rd[2,2]
     COVV.c6=coef.rd[1,1]+(t3+SE[i]-RE[jj.3.tie])*coef.rd[1,2]+t3*(SE[i]-RE[jj.3.tie])*coef.rd[2,2]
     
     ##########
    
     SS.c.11<-coef.rd[1,1]+2*t1*coef.rd[1,2]+t1^2*coef.rd[2,2]+sigma.e^2
     SS.c.22<-coef.rd[1,1]+2*t2*coef.rd[1,2]+t2^2*coef.rd[2,2]+sigma.e^2
     SS.c.33<-coef.rd[1,1]+2*t3*coef.rd[1,2]+t3^2*coef.rd[2,2]+sigma.e^2
     SS.c.12=coef.rd[1,1]+(t1+t2)*coef.rd[1,2]+t1*t2*coef.rd[2,2]
     SS.c.13=coef.rd[1,1]+(t1+t3)*coef.rd[1,2]+t1*t3*coef.rd[2,2]
     SS.c.23=coef.rd[1,1]+(t2+t3)*coef.rd[1,2]+t2*t3*coef.rd[2,2]
     
     SS.c<-array(NA,dim=c(3,3,length(jj.3.tie)))
     SS.c[1,1,]<-SS.c.11
     SS.c[2,2,]<-SS.c.22
     SS.c[3,3,]<-SS.c.33
     SS.c[1,2,]<-SS.c[2,1,]<-SS.c.12
     SS.c[1,3,]<-SS.c[3,1,]<-SS.c.13
     SS.c[2,3,]<-SS.c[3,2,]<-SS.c.23
         
    SS.c<-apply(SS.c,3,solve)
    

     #### SS.c is 9xlength(jj.3), elements in each column corresponds to [1,1],[2,1],[3,1],[1,2],[2,2],[3,2],[1,3],[2,3],[3,3] of a matrix

     delta1=COVV.c1*SS.c[1,]+COVV.c2*SS.c[2,]+COVV.c3*SS.c[3,]
     delta2=COVV.c1*SS.c[4,]+COVV.c2*SS.c[5,]+COVV.c3*SS.c[6,]
     delta3=COVV.c1*SS.c[7,]+COVV.c2*SS.c[8,]+COVV.c3*SS.c[9,]
     delta4=COVV.c4*SS.c[1,]+COVV.c5*SS.c[2,]+COVV.c6*SS.c[3,]
     delta5=COVV.c4*SS.c[4,]+COVV.c5*SS.c[5,]+COVV.c6*SS.c[6,]
     delta6=COVV.c4*SS.c[7,]+COVV.c5*SS.c[8,]+COVV.c6*SS.c[9,]

     
     xi1=delta1*COVV.c1+delta2*COVV.c2+delta3*COVV.c3
     xi2=delta1*COVV.c4+delta2*COVV.c5+delta3*COVV.c6
     xi3=delta4*COVV.c1+delta5*COVV.c2+delta6*COVV.c3
     xi4=delta4*COVV.c4+delta5*COVV.c5+delta6*COVV.c6
      
      mu0.EXE<-cbind(1,0,matrix(as.numeric(WE[jj.3.tie,]),length(jj.3.tie),dw))%*%coef.fix
       mu.EXE<-cbind(1,SE[i]-RE[jj.3.tie],matrix(as.numeric(WE[jj.3.tie,]),length(jj.3.tie),dw))%*%coef.fix
       mut1.EXE<-cbind(1,t1,matrix(as.numeric(WE[jj.3.tie,]),length(jj.3.tie),dw))%*%coef.fix
       mut2.EXE<-cbind(1,t2,matrix(as.numeric(WE[jj.3.tie,]),length(jj.3.tie),dw))%*%coef.fix
       mut3.EXE<-cbind(1,t3,matrix(as.numeric(WE[jj.3.tie,]),length(jj.3.tie),dw))%*%coef.fix
            
      comp1<-parm1*(mu0.EXE+delta1*Xt1+delta2*Xt2+delta3*Xt3-delta1*mut1.EXE-delta2*mut2.EXE-delta3*mut3.EXE)+
        parm2*(mu.EXE+delta4*Xt1+delta5*Xt2+delta6*Xt3-delta4*mut1.EXE-delta5*mut2.EXE-delta6*mut3.EXE)
      
      comp2<-parm1^2*(coef.rd[1,1]-xi1)+parm2*parm1*(2*coef.rd[1,1]+2*(SE[i]-RE[jj.3.tie])*coef.rd[1,2]-xi2-xi3)+
      parm2^2*(coef.rd[1,1]+(SE[i]-RE[jj.3.tie])^2*coef.rd[2,2]+2*(SE[i]-RE[jj.3.tie])*coef.rd[1,2]-xi4)
      
    ss=ss-(ind.tie-1)/di*sum(weights[jj.3.tie]*exp(as.numeric(comp1)+0.5*comp2-comp12i+parm3*(SE[i]-RE[jj.3.tie]-ZE[i]*TE[i])
    +parm0*(ZE[jj.3.tie]-ZE[i])-t(as.numeric(WEc[i,])-matrix(as.numeric(t(WEc[jj.3.tie,])),dc,length(jj.3.tie)))%*%parm4) )  
   
    rm(Xt1,Xt2,Xt3)
    }

    ss.new=ss
    ### placebo
   
   ### need to use jj.Z0 and jj.Z0.tie when there are confounders in the risk model
   if (length(jj.Z0)>0)
   ss=ss+sum(weights[jj.Z0]*exp(-comp12i-parm3*(ZE[i]*TE[i])-parm0*(ZE[i])
   -t(as.numeric(WEc[i,])-matrix(as.numeric(t(WEc[jj.Z0,])),dc,length(jj.Z0)))%*%parm4))
   
   if (length(jj.Z0.tie)>0)
   ss=ss-(ind.tie-1)/di*sum(weights[jj.Z0.tie]*exp(-comp12i-parm3*(ZE[i]*TE[i])-parm0*(ZE[i])-t(as.numeric(WEc[i,])-matrix(as.numeric(t(WEc[jj.Z0.tie,])),dc,length(jj.Z0.tie)))%*%parm4))
   
    out=out+weights.bar*log(weights[i]+ss)+(weights.bar-weights[i])*(comp12i+parm3*(SE[i]-RE[i])*ZE[i]+parm0*(ZE[i])
    +matrix(WEc[i,],length(i),dc)%*%parm4)
    
    outi=c(outi,ss)
 
    }
 return(out)
 }

########

######## Fit a Cox hazard model on calendar time with coefficient for Z, Z*x(t), and W

nloglik.M1.b.Cal.2<- function(parm){

 out=0
 
 outi<-NULL
   ### coef for Z, Z*x(t), and confounders respectively
    parm0=parm[1]
    #parm1=parm[2]
    parm2=parm[2]
    #parm3=parm[4]
    parm4=parm[3]
   
   jj.4.full<-(1:length(RE))[datwide.delta0[match(ID,datwide.delta0$id),]$groupall > datwide.delta0[match(ID,datwide.delta0$id),]$group]
   jj.4.full<-as.numeric(na.omit(jj.4.full))
   out.jj.4.full<-NULL
   for (j in jj.4.full){
      
    keep=1:ll
 
    VVj.11<-coef.rd[1,1]
    VVj.12<-VVj.21<-coef.rd[1,1]+(SE[keep]-RE[j])*coef.rd[1,2]
    VVj.22<-coef.rd[1,1]+2*(SE[keep]-RE[j])*coef.rd[1,2]+(SE[keep]-RE[j])^2*coef.rd[2,2]
    
    muj.1=cbind(1,0,matrix(WE[j,],length(j),dw))%*%coef.fix
    muj.2=cbind(1,SE[keep]-RE[j],matrix(WE[j,],length(j),dw))%*%coef.fix
    
   
    oo.j<-id.long==ID[j]
    Xj<-X.long[oo.j]
    Tj<-T.long[oo.j]
    
    
    VVj.c<-cbind(1,Tj)%*%coef.rd%*%rbind(1,Tj)+diag(rep(sigma.e^2,length(Tj)),nrow=length(Tj),ncol=length(Tj))               
   
    COVVj.1<-coef.rd[1,1]+coef.rd[1,2]*Tj
    COVVj.2<-coef.rd[1,1]+coef.rd[1,2]*outer(SE[keep]-RE[j],Tj,'+')+coef.rd[2,2]*outer(SE[keep]-RE[j],Tj,'*')
    COVVj.1=matrix(COVVj.1,1,length(COVVj.1))
    if (length(Tj)==1) COVVj.2=matrix(COVVj.2,length(COVVj.2),1)
    
    muj.c<-cbind(1,Tj,matrix(WE[j,],length(j),dw))%*%coef.fix
        
    comp1j<-cbind(as.numeric(muj.1),muj.2)+cbind(as.numeric(COVVj.1%*%solve(VVj.c)%*%(Xj-muj.c)),COVVj.2%*%solve(VVj.c)%*%(Xj-muj.c))
    
     VVjminus.11<-COVVj.1%*%solve(VVj.c)%*%(t(COVVj.1))
     VVjminus.12<-COVVj.2%*%solve(VVj.c)%*%(t(COVVj.1))
     VVjminus.22<-diag(COVVj.2%*%solve(VVj.c)%*%(t(COVVj.2)))
     
     comp2j.11<-VVj.11-VVjminus.11
     comp2j.12<-VVj.12-VVjminus.12
     comp2j.22<-VVj.22-VVjminus.22
   
     out.jj.4.full<-cbind(out.jj.4.full,as.numeric(c(parm1,parm2)%*%t(comp1j))+parm1^2*as.numeric(comp2j.11)/2+parm1*parm2*comp2j.12+parm2^2*comp2j.22/2+parm3*(SE[keep]-RE[j])
     +parm0*ZE[j]+matrix(WEc[j,],length(j),dc)%*%parm4)
     
    }
        
    S.tie=SE.u[1]  ### unique SE value
    ind.tie=0       ### position within the unique value
    
    for (i in 1:ll){
    
    ## track #of ties
    
        di<-sum(SE==SE[i] & delta==1)
        
        if (SE[i]==S.tie){  ## in the same tie as previous value
          ind.tie=ind.tie+1
        } else { ## go to a new SE value
           S.tie=SE[i]
           ind.tie=1
        }
     ### compute average weight among ties
        weights.bar=mean(weights[SE==SE[i] & delta==1])
     
    ss=0
    comp12i=0 ### placebo contribute nothing
    if (ZE[i]==1) { 
    VVi<-cbind(1,c(0,TE[i]))%*%coef.rd%*%rbind(1,c(0,TE[i]))
    mui<-rbind(c(1,0,matrix(WE[i,],length(i),dw)),c(1,TE[i],matrix(WE[i,],length(i),dw)))%*%coef.fix
   
    oo.i<-id.long==ID[i]
        
       #### (i) if case has at least one X(Ti)* measured, i.e. case ID[i] belong to id.long
        if (sum(oo.i)>0){ 
        Xi<-X.long[oo.i]
        Ti<-T.long[oo.i]
        
        VVi.c<-cbind(1,Ti)%*%coef.rd%*%rbind(1,Ti)+diag(rep(sigma.e^2,length(Ti)),nrow=length(Ti),ncol=length(Ti))               
        COVVi<-cbind(1,c(0,TE[i]))%*%coef.rd%*%rbind(1,Ti)
        mui.c<-cbind(1,Ti,matrix(WE[i,],length(Ti),dw,byrow=TRUE))%*%coef.fix
        
        comp1i<-mui+COVVi%*%solve(VVi.c)%*%(Xi-mui.c)      ## E(X(0),X(t)|X(Ti)*,W)
        comp2i<-VVi-COVVi%*%solve(VVi.c)%*%t(COVVi)        ## Var(X(0),X(t)|X(Ti)*,W)
        
        } else {
        #### (ii) if case does not have X(Ti)* measure, i.e. case ID[i] not in id.long
        
       comp1i<-mui            
       comp2i<-VVi            
       } 
 
     comp12i<-as.numeric(t(c(parm1,parm2))%*%comp1i+t(c(parm1,parm2))%*%comp2i%*%c(parm1,parm2)/2)
    
    } 
    
    jj.r<-((i+1):(length(RE)))[RE[(i+1):length(RE)]<SE[i] & SE[(i+1):length(RE)]>=SE[i] & ID[(i+1):length(RE)]%in%id.long]
    jj.l<-NULL
    if (i>1) jj.l<-(1:(i-1))[RE[1:(i-1)]<SE[i] & SE[1:(i-1)]>=SE[i] & ID[1:(i-1)]%in%id.long]
    jj<-c(jj.l,jj.r) ## among risk set  
    
    jj.tie<-(1:(length(RE)))[RE[1:length(RE)]<SE[i] & ID[1:length(RE)]%in%id.long & SE[1:length(SE)]==SE[i] & delta[1:length(delta)]==1]
   
    jj.W.r<-((i+1):(length(RE)))[RE[(i+1):length(RE)]<SE[i] & ZE[(i+1):length(ZE)]==1 & SE[(i+1):length(RE)]>=SE[i] & !ID[(i+1):length(RE)]%in%id.long]
    jj.W.l<-NULL
    if (i>1) jj.W.l<-(1:(i-1))[RE[1:(i-1)]<SE[i] & ZE[1:(i-1)]==1 & !ID[1:(i-1)]%in%id.long & SE[1:(i-1)]>=SE[i]]
    jj.W<-c(jj.W.l,jj.W.r)    

    
    jj.W.tie<-(1:(length(RE)))[RE[1:length(RE)]<SE[i] & ZE[1:length(ZE)]==1 & SE[1:length(RE)]==SE[i] & !ID[1:length(RE)]%in%id.long  & delta[1:length(delta)]==1]
 
  
    jj.Z0.r<-((i+1):length(RE))[RE[(i+1):length(RE)]<SE[i] & ZE[(i+1):length(ZE)]==0 & SE[(i+1):length(SE)]>=SE[i]]  
    jj.Z0.l<-NULL
    if (i>1) jj.Z0.l<-(1:(i-1))[RE[1:(i-1)]<SE[i] & ZE[1:(i-1)]==0 & SE[1:(i-1)]>=SE[i]]
    jj.Z0=c(jj.Z0.l,jj.Z0.r)

    jj.Z0.tie<-(1:length(RE))[RE[1:length(RE)]<SE[i] & ZE[1:length(ZE)]==0 & SE[1:length(SE)]==SE[i] & delta[1:length(delta)]==1]  
  
    jj.1.r<-((i+1):(length(RE)))[RE[(i+1):length(RE)]<SE[i] & datwide.delta0[match(ID[(i+1):length(RE)],datwide.delta0$id),]$group==1 &
       datwide.delta0[match(ID[(i+1):length(RE)],datwide.delta0$id),]$groupall==1 & SE[(i+1):length(RE)]>=SE[i]]
    jj.1.l<-NULL
    if (i>1) jj.1.l<-(1:(i-1))[RE[1:(i-1)]<SE[i] & datwide.delta0[match(ID[1:(i-1)],datwide.delta0$id),]$group==1 &
       datwide.delta0[match(ID[1:(i-1)],datwide.delta0$id),]$groupall==1 & SE[1:(i-1)]>=SE[i]]
    jj.1<-c(jj.1.l,jj.1.r) ## among risk set

    
    jj.1.tie<-(1:(length(RE)))[RE[1:length(RE)]<SE[i] & datwide.delta0[match(ID[1:length(RE)],datwide.delta0$id),]$group==1 &
       datwide.delta0[match(ID[1:length(RE)],datwide.delta0$id),]$groupall==1 & SE[1:length(SE)]==SE[i] & delta[1:length(delta)]==1]
   
   
    #######
 
    jj.2.r<-((i+1):(length(RE)))[RE[(i+1):length(RE)]<SE[i] & datwide.delta0[match(ID[(i+1):length(RE)],datwide.delta0$id),]$group==2 & 
       datwide.delta0[match(ID[(i+1):length(RE)],datwide.delta0$id),]$groupall==2 & SE[(i+1):length(SE)]>=SE[i]]
    jj.2.l<-NULL
    if (i>1)   jj.2.l<-(1:(i-1))[RE[1:(i-1)]<SE[i] & datwide.delta0[match(ID[1:(i-1)],datwide.delta0$id),]$group==2 & 
       datwide.delta0[match(ID[1:(i-1)],datwide.delta0$id),]$groupall==2 & SE[1:(i-1)]>=SE[i]]
    jj.2<-c(jj.2.l,jj.2.r)

    jj.2.tie<-(1:(length(RE)))[RE[1:length(RE)]<SE[i] & datwide.delta0[match(ID[1:length(RE)],datwide.delta0$id),]$group==2 & 
       datwide.delta0[match(ID[1:length(RE)],datwide.delta0$id),]$groupall==2 & SE[1:length(SE)]==SE[i] & delta[1:length(delta)]==1]
   

    ###
    jj.3.r<-((i+1):(length(RE)))[RE[(i+1):length(RE)]<SE[i] & datwide.delta0[match(ID[(i+1):length(RE)],datwide.delta0$id),]$group==3 & 
       datwide.delta0[match(ID[(i+1):length(RE)],datwide.delta0$id),]$groupall==3 & SE[(i+1):length(SE)]>=SE[i]]
    jj.3.l<-NULL
    if (i>1) jj.3.l<-(1:(i-1))[RE[1:(i-1)]<SE[i] & datwide.delta0[match(ID[1:(i-1)],datwide.delta0$id),]$group==3 & 
       datwide.delta0[match(ID[1:(i-1)],datwide.delta0$id),]$groupall==3 & SE[1:(i-1)]>=SE[i]]
    jj.3<-c(jj.3.l,jj.3.r)
   

    jj.3.tie<-(1:(length(RE)))[RE[1:length(RE)]<SE[i] & datwide.delta0[match(ID[1:length(RE)],datwide.delta0$id),]$group==3 & 
       datwide.delta0[match(ID[1:length(RE)],datwide.delta0$id),]$groupall==3 & SE[1:length(SE)]==SE[i] & delta[1:length(delta)]==1]
  
     jj.1<-as.numeric(na.omit(jj.1))
     jj.2<-as.numeric(na.omit(jj.2))
     jj.3<-as.numeric(na.omit(jj.3))
   
       jj.1.tie<-as.numeric(na.omit(jj.1.tie))
     jj.2.tie<-as.numeric(na.omit(jj.2.tie))
     jj.3.tie<-as.numeric(na.omit(jj.3.tie))
     
    ss=ss+sum((RE[jj.4.full]<SE[i])*(SE[jj.4.full]>=SE[i])*(!jj.4.full==i)*weights[jj.4.full]*exp(out.jj.4.full[i,]-comp12i-parm3*(ZE[i]*TE[i])
    -parm0*ZE[i]-matrix(WEc[i,],length(i),dc)%*%parm4))  ## risk set from jj.4.full
     
    ss=ss-(ind.tie-1)/di*sum((RE[jj.4.full]<SE[i])*(SE[jj.4.full]==SE[i])*(!jj.4.full==i)*weights[jj.4.full]*exp(out.jj.4.full[i,]-comp12i-parm3*(ZE[i]*TE[i])
    -parm0*ZE[i]-matrix(WEc[i,],length(i),dc)%*%parm4))  ## tie from jj.4.full
     

   if (length(jj.W)>0){
      ### No X measure 
    mu0.EXE<-cbind(1,0,matrix(as.numeric(WE[jj.W,]),length(jj.W),dw))%*%coef.fix
    mu.EXE<-cbind(1,SE[i]-RE[jj.W],matrix(as.numeric(WE[jj.W,]),length(jj.W),dw))%*%coef.fix
  
   comp1.b<-parm1*mu0.EXE+parm2*mu.EXE
       
   comp2.b<-parm1^2*coef.rd[1,1]+2*parm1*parm2*(coef.rd[1,1]+(SE[i]-RE[jj.W])*coef.rd[1,2])+
   parm2^2*(coef.rd[1,1]+2*(SE[i]-RE[jj.W])*coef.rd[1,2]+(SE[i]-RE[jj.W])^2*coef.rd[2,2])
   
   ss=ss+sum(weights[jj.W]*exp(as.numeric(comp1.b)+0.5*comp2.b-comp12i+parm3*(SE[i]-RE[jj.W]-ZE[i]*TE[i])
   +parm0*(ZE[jj.W]-ZE[i])-t(as.numeric(WEc[i,])-matrix(as.numeric(t(WEc[jj.W,])),dc,length(jj.W)))%*%parm4) )  
    
   } 
   ### minus ties
   if (length(jj.W.tie)>0){
   
        mu0.EXE<-cbind(1,0,matrix(as.numeric(WE[jj.W.tie,]),length(jj.W.tie),dw))%*%coef.fix
        mu.EXE<-cbind(1,SE[i]-RE[jj.W.tie],matrix(as.numeric(WE[jj.W.tie,]),length(jj.W.tie),dw))%*%coef.fix
      
       comp1.b<-parm1*mu0.EXE+parm2*mu.EXE
           
       comp2.b<-parm1^2*coef.rd[1,1]+2*parm1*parm2*(coef.rd[1,1]+(SE[i]-RE[jj.W.tie])*coef.rd[1,2])+
       parm2^2*(coef.rd[1,1]+2*(SE[i]-RE[jj.W.tie])*coef.rd[1,2]+(SE[i]-RE[jj.W.tie])^2*coef.rd[2,2])
       
       ss=ss-(ind.tie-1)/di*sum(weights[jj.W.tie]*exp(as.numeric(comp1.b)+0.5*comp2.b-comp12i+parm3*(SE[i]-RE[jj.W.tie]-ZE[i]*TE[i])
       +parm0*(ZE[jj.W.tie]-ZE[i])-t(as.numeric(WEc[i,])-matrix(as.numeric(t(WEc[jj.W.tie,])),dc,length(jj.W.tie)))%*%parm4) ) 
   }      

 
   ### Measure of X(0)* only
   
    #account for different t levels
   
   if (length(jj.1)>0){
   Xt1=datwide.delta0[match(ID[jj.1],datwide.delta0$id),Vname[1]]
    t1=datwide.delta0[match(ID[jj.1],datwide.delta0$id),Tname[1]]
 
    ## make t1 vector
   COVV.a<-matrix(c(coef.rd[1,1]+t1*coef.rd[1,2],coef.rd[1,1]+(t1+SE[i]-RE[jj.1])*coef.rd[1,2]+t1*(SE[i]-RE[jj.1])*coef.rd[2,2]),length(jj.1),2,byrow=F)
   SS.a<-sapply(coef.rd[1,1]+2*t1*coef.rd[1,2]+t1^2*coef.rd[2,2]+sigma.e^2,solve)
   
   slope.a<-COVV.a*as.numeric(SS.a)
   Sigma1.2.a=outer(t(slope.a),(COVV.a),'*'); 
   kappa1=diag(as.matrix(Sigma1.2.a[1,,,1])); kappa2=diag(as.matrix(Sigma1.2.a[1,,,2])); 
   kappa3=diag(as.matrix(Sigma1.2.a[2,,,1])); kappa4=diag(as.matrix(Sigma1.2.a[2,,,2]))
   
     mu0.EXE<-cbind(1,0,matrix(as.numeric(WE[jj.1,]),length(jj.1),dw))%*%coef.fix
      mu.EXE<-cbind(1,SE[i]-RE[jj.1],matrix(as.numeric(WE[jj.1,]),length(jj.1),dw))%*%coef.fix
  
       comp1<-parm1*(mu0.EXE+slope.a[,1]*(Xt1-mu0.EXE))+
       parm2*(mu.EXE+slope.a[,2]*(Xt1-mu0.EXE))
       
       comp2<-parm1^2*(coef.rd[1,1]-kappa1)+
       parm2*parm1*(2*coef.rd[1,1]+2*(SE[i]-RE[jj.1])*coef.rd[1,2]-kappa2-kappa3)+
       parm2^2*(coef.rd[1,1]+(SE[i]-RE[jj.1])^2*coef.rd[2,2]+2*(SE[i]-RE[jj.1])*coef.rd[1,2]-kappa4)
    
    ss=ss+sum(weights[jj.1]*exp(as.numeric(comp1)+0.5*comp2-comp12i+parm3*(SE[i]-RE[jj.1]-ZE[i]*TE[i])
    +parm0*(ZE[jj.1]-ZE[i])-t(as.numeric(WEc[i,])-matrix(as.numeric(t(WEc[jj.1,])),dc,length(jj.1)))%*%parm4) )  
   
   rm(Xt1)
   }

   ## minus ties

   if (length(jj.1.tie)>0){
   Xt1=datwide.delta0[match(ID[jj.1.tie],datwide.delta0$id),Vname[1]]
   
   t1=datwide.delta0[match(ID[jj.1.tie],datwide.delta0$id),Tname[1]]
   
    # make t1 vector
   COVV.a<-matrix(c(coef.rd[1,1]+t1*coef.rd[1,2],coef.rd[1,1]+(t1+SE[i]-RE[jj.1.tie])*coef.rd[1,2]+t1*(SE[i]-RE[jj.1.tie])*coef.rd[2,2]),length(jj.1.tie),2,byrow=F)
   SS.a<-sapply(coef.rd[1,1]+2*t1*coef.rd[1,2]+t1^2*coef.rd[2,2]+sigma.e^2,solve)
 
   
   slope.a<-COVV.a*as.numeric(SS.a)
   Sigma1.2.a=outer(t(slope.a),(COVV.a),'*'); 
   kappa1=diag(as.matrix(Sigma1.2.a[1,,,1])); kappa2=diag(as.matrix(Sigma1.2.a[1,,,2])); 
   kappa3=diag(as.matrix(Sigma1.2.a[2,,,1])); kappa4=diag(as.matrix(Sigma1.2.a[2,,,2]))
   
     mu0.EXE<-cbind(1,0,matrix(as.numeric(WE[jj.1.tie,]),length(jj.1.tie),dw))%*%coef.fix
      mu.EXE<-cbind(1,SE[i]-RE[jj.1.tie],matrix(as.numeric(WE[jj.1.tie,]),length(jj.1.tie),dw))%*%coef.fix
  
       comp1<-parm1*(mu0.EXE+slope.a[,1]*(Xt1-mu0.EXE))+
       parm2*(mu.EXE+slope.a[,2]*(Xt1-mu0.EXE))
       
       comp2<-parm1^2*(coef.rd[1,1]-kappa1)+
       parm2*parm1*(2*coef.rd[1,1]+2*(SE[i]-RE[jj.1.tie])*coef.rd[1,2]-kappa2-kappa3)+
       parm2^2*(coef.rd[1,1]+(SE[i]-RE[jj.1.tie])^2*coef.rd[2,2]+2*(SE[i]-RE[jj.1.tie])*coef.rd[1,2]-kappa4)
    
    ss=ss-(ind.tie-1)/di*sum(weights[jj.1.tie]*exp(as.numeric(comp1)+0.5*comp2-comp12i+parm3*(SE[i]-RE[jj.1.tie]-ZE[i]*TE[i])
    +parm0*(ZE[jj.1.tie]-ZE[i])-t(as.numeric(WEc[i,])-matrix(as.numeric(t(WEc[jj.1.tie,])),dc,length(jj.1.tie)))%*%parm4 ) )  
   
   rm(Xt1)
   }



   ###### Measure of X(0)* and X(t2)* 
     if (length(jj.2)>0){
     Xt1=datwide.delta0[match(ID[jj.2],datwide.delta0$id),Vname[1]]
     Xt2=datwide.delta0[match(ID[jj.2],datwide.delta0$id),Vname[2]]
     
       ## make t1, t2 vectors
      t1=datwide.delta0[match(ID[jj.2],datwide.delta0$id),Tname[1]]
      t2=datwide.delta0[match(ID[jj.2],datwide.delta0$id),Tname[2]]
         
     COVV.b1<-coef.rd[1,1]+t1*coef.rd[1,2]
     COVV.b2<-coef.rd[1,1]+t2*coef.rd[1,2]
     COVV.b3<-coef.rd[1,1]+(t1+SE[i]-RE[jj.2])*coef.rd[1,2]+t1*(SE[i]-RE[jj.2])*coef.rd[2,2]
     COVV.b4<-coef.rd[1,1]+t2*(SE[i]-RE[jj.2])*coef.rd[2,2]+(t2+SE[i]-RE[jj.2])*coef.rd[1,2]
     

     SS.b.1<-coef.rd[1,1]+2*t1*coef.rd[1,2]+t1*2*coef.rd[2,2]+sigma.e^2
     SS.b.2<-coef.rd[1,1]+2*t2*coef.rd[1,2]+t2^2*coef.rd[2,2]+sigma.e^2
     SS.b.12<-coef.rd[1,1]+(t1+t2)*coef.rd[1,2]+t1*t2*coef.rd[2,2]
     
     SS.b<-array(NA,dim=c(2,2,length(jj.2)))
     SS.b[1,1,]<-SS.b.1
     SS.b[1,2,]<-SS.b[2,1,]<-SS.b.12
     SS.b[2,2,]<-SS.b.2
     
     SS.b<-apply(SS.b,3,solve)  ### SS.b is a 4 by length(jj.2) matrix with each column the inverse of a matrix
 
     ### slope
    
     ### When t1 t2 are vectors, SS.b[1,], SS.b[2,], SS.b[3,], SS.b[4,] 
     ### corresponds [1,1],[2,1],[1,2],[2,2] of a 2x2 matrix
     
     gamma1=COVV.b1*SS.b[1,]+COVV.b2*SS.b[2,]
     gamma2=COVV.b1*SS.b[3,]+COVV.b2*SS.b[4,]
     gamma3=COVV.b3*SS.b[1,]+COVV.b4*SS.b[2,]
     gamma4=COVV.b3*SS.b[3,]+COVV.b4*SS.b[4,]
     

     ### Sigma1.2.b
     eta1=gamma1*COVV.b1+gamma2*COVV.b2
     eta2=gamma1*COVV.b3+gamma2*COVV.b4
     eta3=gamma3*COVV.b1+gamma4*COVV.b2
     eta4=gamma3*COVV.b3+gamma4*COVV.b4
     
      mu0.EXE<-cbind(1,0,matrix(as.numeric(WE[jj.2,]),length(jj.2),dw))%*%coef.fix
       mu.EXE<-cbind(1,SE[i]-RE[jj.2],matrix(as.numeric(WE[jj.2,]),length(jj.2),dw))%*%coef.fix
       mut1.EXE<-cbind(1,t1,matrix(as.numeric(WE[jj.2,]),length(jj.2),dw))%*%coef.fix
       mut2.EXE<-cbind(1,t2,matrix(as.numeric(WE[jj.2,]),length(jj.2),dw))%*%coef.fix
       
      
      comp1<-parm1*(mu0.EXE+gamma1*Xt1+gamma2*Xt2-gamma1*mut1.EXE-gamma2*mut2.EXE)+
        parm2*(mu.EXE+gamma3*Xt1+gamma4*Xt2-gamma3*mut1.EXE-gamma4*mut2.EXE)
      
      comp2<-parm1^2*(coef.rd[1,1]-eta1)+parm2*parm1*(2*coef.rd[1,1]+2*(SE[i]-RE[jj.2])*coef.rd[1,2]-eta2-eta3)+
      parm2^2*(coef.rd[1,1]+(SE[i]-RE[jj.2])^2*coef.rd[2,2]+2*(SE[i]-RE[jj.2])*coef.rd[1,2]-eta4)
      
    
    ss=ss+sum(weights[jj.2]*exp(as.numeric(comp1)+0.5*comp2-comp12i+parm3*(SE[i]-RE[jj.2]-ZE[i]*TE[i])
    +parm0*(ZE[jj.2]-ZE[i])-t(as.numeric(WEc[i,])-matrix(as.numeric(t(WEc[jj.2,])),dc,length(jj.2)))%*%parm4) )  
    rm(Xt1,Xt2)
    }

    ## minus ties

     if (length(jj.2.tie)>0){
     Xt1=datwide.delta0[match(ID[jj.2.tie],datwide.delta0$id),Vname[1]]
     Xt2=datwide.delta0[match(ID[jj.2.tie],datwide.delta0$id),Vname[2]]
     
         ## make t1, t2 vectors
     t1=datwide.delta0[match(ID[jj.2.tie],datwide.delta0$id),Tname[1]]
     t2=datwide.delta0[match(ID[jj.2.tie],datwide.delta0$id),Tname[2]]
    

   ### when t1 t2 are vectors 
     COVV.b1<-coef.rd[1,1]+t1*coef.rd[1,2]
     COVV.b2<-coef.rd[1,1]+t2*coef.rd[1,2]
     COVV.b3<-coef.rd[1,1]+(t1+SE[i]-RE[jj.2.tie])*coef.rd[1,2]+t1*(SE[i]-RE[jj.2.tie])*coef.rd[2,2]
     COVV.b4<-coef.rd[1,1]+t2*(SE[i]-RE[jj.2.tie])*coef.rd[2,2]+(t2+SE[i]-RE[jj.2.tie])*coef.rd[1,2]
     

     SS.b.1<-coef.rd[1,1]+2*t1*coef.rd[1,2]+t1*2*coef.rd[2,2]+sigma.e^2
     SS.b.2<-coef.rd[1,1]+2*t2*coef.rd[1,2]+t2^2*coef.rd[2,2]+sigma.e^2
     SS.b.12<-coef.rd[1,1]+(t1+t2)*coef.rd[1,2]+t1*t2*coef.rd[2,2]
     
     SS.b<-array(NA,dim=c(2,2,length(jj.2.tie)))
     SS.b[1,1,]<-SS.b.1
     SS.b[1,2,]<-SS.b[2,1,]<-SS.b.12
     SS.b[2,2,]<-SS.b.2
     
     SS.b<-apply(SS.b,3,solve)  ### SS.b is a 4 by length(jj.2) matrix with each column the inverse of a matrix
 
     ### slope
    
     ### When t1 t2 are vectors , SS.b[1,], SS.b[2,], SS.b[3,], SS.b[4,] 
     ### corresponds [1,1],[2,1],[1,2],[2,2] of a 2x2 matrix
     
     gamma1=COVV.b1*SS.b[1,]+COVV.b2*SS.b[2,]
     gamma2=COVV.b1*SS.b[3,]+COVV.b2*SS.b[4,]
     gamma3=COVV.b3*SS.b[1,]+COVV.b4*SS.b[2,]
     gamma4=COVV.b3*SS.b[3,]+COVV.b4*SS.b[4,]
     
     ### Sigma1.2.b
     eta1=gamma1*COVV.b1+gamma2*COVV.b2
     eta2=gamma1*COVV.b3+gamma2*COVV.b4
     eta3=gamma3*COVV.b1+gamma4*COVV.b2
     eta4=gamma3*COVV.b3+gamma4*COVV.b4
     
      mu0.EXE<-cbind(1,0,matrix(as.numeric(WE[jj.2.tie,]),length(jj.2.tie),dw))%*%coef.fix
       mu.EXE<-cbind(1,SE[i]-RE[jj.2.tie],matrix(as.numeric(WE[jj.2.tie,]),length(jj.2.tie),dw))%*%coef.fix
       mut1.EXE<-cbind(1,t1,matrix(as.numeric(WE[jj.2.tie,]),length(jj.2.tie),dw))%*%coef.fix
       mut2.EXE<-cbind(1,t2,matrix(as.numeric(WE[jj.2.tie,]),length(jj.2.tie),dw))%*%coef.fix
       
      
      comp1<-parm1*(mu0.EXE+gamma1*Xt1+gamma2*Xt2-gamma1*mut1.EXE-gamma2*mut2.EXE)+
        parm2*(mu.EXE+gamma3*Xt1+gamma4*Xt2-gamma3*mut1.EXE-gamma4*mut2.EXE)
      
      comp2<-parm1^2*(coef.rd[1,1]-eta1)+parm2*parm1*(2*coef.rd[1,1]+2*(SE[i]-RE[jj.2.tie])*coef.rd[1,2]-eta2-eta3)+
      parm2^2*(coef.rd[1,1]+(SE[i]-RE[jj.2.tie])^2*coef.rd[2,2]+2*(SE[i]-RE[jj.2.tie])*coef.rd[1,2]-eta4)
      
    
    ss=ss-(ind.tie-1)/di*sum(weights[jj.2.tie]*exp(as.numeric(comp1)+0.5*comp2-comp12i+parm3*(SE[i]-RE[jj.2.tie]-ZE[i]*TE[i])
    +parm0*(ZE[jj.2.tie]-ZE[i])-t(as.numeric(WEc[i,])-matrix(as.numeric(t(WEc[jj.2.tie,])),dc,length(jj.2.tie)))%*%parm4) )  
    rm(Xt1,Xt2)
    }
    #####  ###### Measure of X(t1)*, X(t2)*, and X(t3)*
    
     if (length(jj.3)>0){
     Xt1=datwide.delta0[match(ID[jj.3],datwide.delta0$id),Vname[1]]
     Xt2=datwide.delta0[match(ID[jj.3],datwide.delta0$id),Vname[2]]
     Xt3=datwide.delta0[match(ID[jj.3],datwide.delta0$id),Vname[3]]
     
       ### t1, t2, t3 vector
     t1=datwide.delta0[match(ID[jj.3],datwide.delta0$id),Tname[1]]
     t2=datwide.delta0[match(ID[jj.3],datwide.delta0$id),Tname[2]]
     t3=datwide.delta0[match(ID[jj.3],datwide.delta0$id),Tname[3]]
     
   ####  t1,t2,t3 vector

     COVV.c1=coef.rd[1,1]+t1*coef.rd[1,2]
     COVV.c2=coef.rd[1,1]+t2*coef.rd[1,2]
     COVV.c3=coef.rd[1,1]+t3*coef.rd[1,2]
     COVV.c4=coef.rd[1,1]+(t1+SE[i]-RE[jj.3])*coef.rd[1,2]+t1*(SE[i]-RE[jj.3])*coef.rd[2,2]
     COVV.c5=coef.rd[1,1]+(t2+SE[i]-RE[jj.3])*coef.rd[1,2]+t2*(SE[i]-RE[jj.3])*coef.rd[2,2]
     COVV.c6=coef.rd[1,1]+(t3+SE[i]-RE[jj.3])*coef.rd[1,2]+t3*(SE[i]-RE[jj.3])*coef.rd[2,2]
     
     ##########
    
     SS.c.11<-coef.rd[1,1]+2*t1*coef.rd[1,2]+t1^2*coef.rd[2,2]+sigma.e^2
     SS.c.22<-coef.rd[1,1]+2*t2*coef.rd[1,2]+t2^2*coef.rd[2,2]+sigma.e^2
     SS.c.33<-coef.rd[1,1]+2*t3*coef.rd[1,2]+t3^2*coef.rd[2,2]+sigma.e^2
     SS.c.12=coef.rd[1,1]+(t1+t2)*coef.rd[1,2]+t1*t2*coef.rd[2,2]
     SS.c.13=coef.rd[1,1]+(t1+t3)*coef.rd[1,2]+t1*t3*coef.rd[2,2]
     SS.c.23=coef.rd[1,1]+(t2+t3)*coef.rd[1,2]+t2*t3*coef.rd[2,2]
     
     SS.c<-array(NA,dim=c(3,3,length(jj.3)))
     SS.c[1,1,]<-SS.c.11
     SS.c[2,2,]<-SS.c.22
     SS.c[3,3,]<-SS.c.33
     SS.c[1,2,]<-SS.c[2,1,]<-SS.c.12
     SS.c[1,3,]<-SS.c[3,1,]<-SS.c.13
     SS.c[2,3,]<-SS.c[3,2,]<-SS.c.23
     
    SS.c<-apply(SS.c,3,solve)
    

     #### SS.c is 9xlength(jj.3), elements in each column corresponds to [1,1],[2,1],[3,1],[1,2],[2,2],[3,2],[1,3],[2,3],[3,3] of a matrix

     delta1=COVV.c1*SS.c[1,]+COVV.c2*SS.c[2,]+COVV.c3*SS.c[3,]
     delta2=COVV.c1*SS.c[4,]+COVV.c2*SS.c[5,]+COVV.c3*SS.c[6,]
     delta3=COVV.c1*SS.c[7,]+COVV.c2*SS.c[8,]+COVV.c3*SS.c[9,]
     delta4=COVV.c4*SS.c[1,]+COVV.c5*SS.c[2,]+COVV.c6*SS.c[3,]
     delta5=COVV.c4*SS.c[4,]+COVV.c5*SS.c[5,]+COVV.c6*SS.c[6,]
     delta6=COVV.c4*SS.c[7,]+COVV.c5*SS.c[8,]+COVV.c6*SS.c[9,]

     xi1=delta1*COVV.c1+delta2*COVV.c2+delta3*COVV.c3
     xi2=delta1*COVV.c4+delta2*COVV.c5+delta3*COVV.c6
     xi3=delta4*COVV.c1+delta5*COVV.c2+delta6*COVV.c3
     xi4=delta4*COVV.c4+delta5*COVV.c5+delta6*COVV.c6
      
      mu0.EXE<-cbind(1,0,matrix(as.numeric(WE[jj.3,]),length(jj.3),dw))%*%coef.fix
       mu.EXE<-cbind(1,SE[i]-RE[jj.3],matrix(as.numeric(WE[jj.3,]),length(jj.3),dw))%*%coef.fix
       mut1.EXE<-cbind(1,t1,matrix(as.numeric(WE[jj.3,]),length(jj.3),dw))%*%coef.fix
       mut2.EXE<-cbind(1,t2,matrix(as.numeric(WE[jj.3,]),length(jj.3),dw))%*%coef.fix
       mut3.EXE<-cbind(1,t3,matrix(as.numeric(WE[jj.3,]),length(jj.3),dw))%*%coef.fix
       
      
      comp1<-parm1*(mu0.EXE+delta1*Xt1+delta2*Xt2+delta3*Xt3-delta1*mut1.EXE-delta2*mut2.EXE-delta3*mut3.EXE)+
        parm2*(mu.EXE+delta4*Xt1+delta5*Xt2+delta6*Xt3-delta4*mut1.EXE-delta5*mut2.EXE-delta6*mut3.EXE)
      
      comp2<-parm1^2*(coef.rd[1,1]-xi1)+parm2*parm1*(2*coef.rd[1,1]+2*(SE[i]-RE[jj.3])*coef.rd[1,2]-xi2-xi3)+
      parm2^2*(coef.rd[1,1]+(SE[i]-RE[jj.3])^2*coef.rd[2,2]+2*(SE[i]-RE[jj.3])*coef.rd[1,2]-xi4)
      
    
    ss=ss+sum(weights[jj.3]*exp(as.numeric(comp1)+0.5*comp2-comp12i+parm3*(SE[i]-RE[jj.3]-ZE[i]*TE[i])
    +parm0*(ZE[jj.3]-ZE[i])-t(as.numeric(WEc[i,])-matrix(as.numeric(t(WEc[jj.3,])),dc,length(jj.3)))%*%parm4) )  
   
    rm(Xt1,Xt2,Xt3)
    }

    #### Minus ties
 
     if (length(jj.3.tie)>0){
     Xt1=datwide.delta0[match(ID[jj.3.tie],datwide.delta0$id),Vname[1]]
     Xt2=datwide.delta0[match(ID[jj.3.tie],datwide.delta0$id),Vname[2]]
     Xt3=datwide.delta0[match(ID[jj.3.tie],datwide.delta0$id),Vname[3]]
     
      ### t1, t2, t3 vector
     t1=datwide.delta0[match(ID[jj.3.tie],datwide.delta0$id),Tname[1]]
     t2=datwide.delta0[match(ID[jj.3.tie],datwide.delta0$id),Tname[2]]
     t3=datwide.delta0[match(ID[jj.3.tie],datwide.delta0$id),Tname[3]]
 
     COVV.c1=coef.rd[1,1]+t1*coef.rd[1,2]
     COVV.c2=coef.rd[1,1]+t2*coef.rd[1,2]
     COVV.c3=coef.rd[1,1]+t3*coef.rd[1,2]
     COVV.c4=coef.rd[1,1]+(t1+SE[i]-RE[jj.3.tie])*coef.rd[1,2]+t1*(SE[i]-RE[jj.3.tie])*coef.rd[2,2]
     COVV.c5=coef.rd[1,1]+(t2+SE[i]-RE[jj.3.tie])*coef.rd[1,2]+t2*(SE[i]-RE[jj.3.tie])*coef.rd[2,2]
     COVV.c6=coef.rd[1,1]+(t3+SE[i]-RE[jj.3.tie])*coef.rd[1,2]+t3*(SE[i]-RE[jj.3.tie])*coef.rd[2,2]
     
     ##########
    
     SS.c.11<-coef.rd[1,1]+2*t1*coef.rd[1,2]+t1^2*coef.rd[2,2]+sigma.e^2
     SS.c.22<-coef.rd[1,1]+2*t2*coef.rd[1,2]+t2^2*coef.rd[2,2]+sigma.e^2
     SS.c.33<-coef.rd[1,1]+2*t3*coef.rd[1,2]+t3^2*coef.rd[2,2]+sigma.e^2
     SS.c.12=coef.rd[1,1]+(t1+t2)*coef.rd[1,2]+t1*t2*coef.rd[2,2]
     SS.c.13=coef.rd[1,1]+(t1+t3)*coef.rd[1,2]+t1*t3*coef.rd[2,2]
     SS.c.23=coef.rd[1,1]+(t2+t3)*coef.rd[1,2]+t2*t3*coef.rd[2,2]
     
     SS.c<-array(NA,dim=c(3,3,length(jj.3.tie)))
     SS.c[1,1,]<-SS.c.11
     SS.c[2,2,]<-SS.c.22
     SS.c[3,3,]<-SS.c.33
     SS.c[1,2,]<-SS.c[2,1,]<-SS.c.12
     SS.c[1,3,]<-SS.c[3,1,]<-SS.c.13
     SS.c[2,3,]<-SS.c[3,2,]<-SS.c.23
         
    SS.c<-apply(SS.c,3,solve)
    

     #### SS.c is 9xlength(jj.3), elements in each column corresponds to [1,1],[2,1],[3,1],[1,2],[2,2],[3,2],[1,3],[2,3],[3,3] of a matrix

     delta1=COVV.c1*SS.c[1,]+COVV.c2*SS.c[2,]+COVV.c3*SS.c[3,]
     delta2=COVV.c1*SS.c[4,]+COVV.c2*SS.c[5,]+COVV.c3*SS.c[6,]
     delta3=COVV.c1*SS.c[7,]+COVV.c2*SS.c[8,]+COVV.c3*SS.c[9,]
     delta4=COVV.c4*SS.c[1,]+COVV.c5*SS.c[2,]+COVV.c6*SS.c[3,]
     delta5=COVV.c4*SS.c[4,]+COVV.c5*SS.c[5,]+COVV.c6*SS.c[6,]
     delta6=COVV.c4*SS.c[7,]+COVV.c5*SS.c[8,]+COVV.c6*SS.c[9,]

     
     xi1=delta1*COVV.c1+delta2*COVV.c2+delta3*COVV.c3
     xi2=delta1*COVV.c4+delta2*COVV.c5+delta3*COVV.c6
     xi3=delta4*COVV.c1+delta5*COVV.c2+delta6*COVV.c3
     xi4=delta4*COVV.c4+delta5*COVV.c5+delta6*COVV.c6
      
      mu0.EXE<-cbind(1,0,matrix(as.numeric(WE[jj.3.tie,]),length(jj.3.tie),dw))%*%coef.fix
       mu.EXE<-cbind(1,SE[i]-RE[jj.3.tie],matrix(as.numeric(WE[jj.3.tie,]),length(jj.3.tie),dw))%*%coef.fix
       mut1.EXE<-cbind(1,t1,matrix(as.numeric(WE[jj.3.tie,]),length(jj.3.tie),dw))%*%coef.fix
       mut2.EXE<-cbind(1,t2,matrix(as.numeric(WE[jj.3.tie,]),length(jj.3.tie),dw))%*%coef.fix
       mut3.EXE<-cbind(1,t3,matrix(as.numeric(WE[jj.3.tie,]),length(jj.3.tie),dw))%*%coef.fix
            
      comp1<-parm1*(mu0.EXE+delta1*Xt1+delta2*Xt2+delta3*Xt3-delta1*mut1.EXE-delta2*mut2.EXE-delta3*mut3.EXE)+
        parm2*(mu.EXE+delta4*Xt1+delta5*Xt2+delta6*Xt3-delta4*mut1.EXE-delta5*mut2.EXE-delta6*mut3.EXE)
      
      comp2<-parm1^2*(coef.rd[1,1]-xi1)+parm2*parm1*(2*coef.rd[1,1]+2*(SE[i]-RE[jj.3.tie])*coef.rd[1,2]-xi2-xi3)+
      parm2^2*(coef.rd[1,1]+(SE[i]-RE[jj.3.tie])^2*coef.rd[2,2]+2*(SE[i]-RE[jj.3.tie])*coef.rd[1,2]-xi4)
      
    ss=ss-(ind.tie-1)/di*sum(weights[jj.3.tie]*exp(as.numeric(comp1)+0.5*comp2-comp12i+parm3*(SE[i]-RE[jj.3.tie]-ZE[i]*TE[i])
    +parm0*(ZE[jj.3.tie]-ZE[i])-t(as.numeric(WEc[i,])-matrix(as.numeric(t(WEc[jj.3.tie,])),dc,length(jj.3.tie)))%*%parm4) )  
   
    rm(Xt1,Xt2,Xt3)
    }

    ss.new=ss
    ### placebo
   
   ### need to use jj.Z0 and jj.Z0.tie when there are confounders in the risk model
   if (length(jj.Z0)>0)
   ss=ss+sum(weights[jj.Z0]*exp(-comp12i-parm3*(ZE[i]*TE[i])-parm0*(ZE[i])
   -t(as.numeric(WEc[i,])-matrix(as.numeric(t(WEc[jj.Z0,])),dc,length(jj.Z0)))%*%parm4))
   
   if (length(jj.Z0.tie)>0)
   ss=ss-(ind.tie-1)/di*sum(weights[jj.Z0.tie]*exp(-comp12i-parm3*(ZE[i]*TE[i])-parm0*(ZE[i])-t(as.numeric(WEc[i,])-matrix(as.numeric(t(WEc[jj.Z0.tie,])),dc,length(jj.Z0.tie)))%*%parm4))
   
    out=out+weights.bar*log(weights[i]+ss)+(weights.bar-weights[i])*(comp12i+parm3*(SE[i]-RE[i])*ZE[i]+parm0*(ZE[i])
    +matrix(WEc[i,],length(i),dc)%*%parm4)
    
    outi=c(outi,ss)
 
    }
 return(out)
 }

########

######## Fit a Cox hazard model on calendar time with coefficient for Z, Z*t, and W

nloglik.M1.b.Cal.3<- function(parm){

 out=0
 
 outi<-NULL
   ### coef for Z, Z*t, and confounders respectively
    parm0=parm[1]
    #parm1=parm[2]
    #parm2=parm[2]
    parm3=parm[2]
    parm4=parm[3]
   
   jj.4.full<-(1:length(RE))[datwide.delta0[match(ID,datwide.delta0$id),]$groupall > datwide.delta0[match(ID,datwide.delta0$id),]$group]
   jj.4.full<-as.numeric(na.omit(jj.4.full))
   out.jj.4.full<-NULL
   for (j in jj.4.full){
      
    keep=1:ll
 
    VVj.11<-coef.rd[1,1]
    VVj.12<-VVj.21<-coef.rd[1,1]+(SE[keep]-RE[j])*coef.rd[1,2]
    VVj.22<-coef.rd[1,1]+2*(SE[keep]-RE[j])*coef.rd[1,2]+(SE[keep]-RE[j])^2*coef.rd[2,2]
    
    muj.1=cbind(1,0,matrix(WE[j,],length(j),dw))%*%coef.fix
    muj.2=cbind(1,SE[keep]-RE[j],matrix(WE[j,],length(j),dw))%*%coef.fix
    
   
    oo.j<-id.long==ID[j]
    Xj<-X.long[oo.j]
    Tj<-T.long[oo.j]
    
    
    VVj.c<-cbind(1,Tj)%*%coef.rd%*%rbind(1,Tj)+diag(rep(sigma.e^2,length(Tj)),nrow=length(Tj),ncol=length(Tj))               
   
    COVVj.1<-coef.rd[1,1]+coef.rd[1,2]*Tj
    COVVj.2<-coef.rd[1,1]+coef.rd[1,2]*outer(SE[keep]-RE[j],Tj,'+')+coef.rd[2,2]*outer(SE[keep]-RE[j],Tj,'*')
    COVVj.1=matrix(COVVj.1,1,length(COVVj.1))
    if (length(Tj)==1) COVVj.2=matrix(COVVj.2,length(COVVj.2),1)
    
    muj.c<-cbind(1,Tj,matrix(WE[j,],length(j),dw))%*%coef.fix
        
    comp1j<-cbind(as.numeric(muj.1),muj.2)+cbind(as.numeric(COVVj.1%*%solve(VVj.c)%*%(Xj-muj.c)),COVVj.2%*%solve(VVj.c)%*%(Xj-muj.c))
    
     VVjminus.11<-COVVj.1%*%solve(VVj.c)%*%(t(COVVj.1))
     VVjminus.12<-COVVj.2%*%solve(VVj.c)%*%(t(COVVj.1))
     VVjminus.22<-diag(COVVj.2%*%solve(VVj.c)%*%(t(COVVj.2)))
     
     comp2j.11<-VVj.11-VVjminus.11
     comp2j.12<-VVj.12-VVjminus.12
     comp2j.22<-VVj.22-VVjminus.22
   
     out.jj.4.full<-cbind(out.jj.4.full,as.numeric(c(parm1,parm2)%*%t(comp1j))+parm1^2*as.numeric(comp2j.11)/2+parm1*parm2*comp2j.12+parm2^2*comp2j.22/2+parm3*(SE[keep]-RE[j])
     +parm0*ZE[j]+matrix(WEc[j,],length(j),dc)%*%parm4)
     
    }
        
    S.tie=SE.u[1]  ### unique SE value
    ind.tie=0       ### position within the unique value
    
    for (i in 1:ll){
    
    ## track #of ties
    
        di<-sum(SE==SE[i] & delta==1)
        
        if (SE[i]==S.tie){  ## in the same tie as previous value
          ind.tie=ind.tie+1
        } else { ## go to a new SE value
           S.tie=SE[i]
           ind.tie=1
        }
     ### compute average weight among ties
        weights.bar=mean(weights[SE==SE[i] & delta==1])
     
    ss=0
    comp12i=0 ### placebo contribute nothing
    if (ZE[i]==1) { 
    VVi<-cbind(1,c(0,TE[i]))%*%coef.rd%*%rbind(1,c(0,TE[i]))
    mui<-rbind(c(1,0,matrix(WE[i,],length(i),dw)),c(1,TE[i],matrix(WE[i,],length(i),dw)))%*%coef.fix
   
    oo.i<-id.long==ID[i]
        
       #### (i) if case has at least one X(Ti)* measured, i.e. case ID[i] belong to id.long
        if (sum(oo.i)>0){ 
        Xi<-X.long[oo.i]
        Ti<-T.long[oo.i]
        
        VVi.c<-cbind(1,Ti)%*%coef.rd%*%rbind(1,Ti)+diag(rep(sigma.e^2,length(Ti)),nrow=length(Ti),ncol=length(Ti))               
        COVVi<-cbind(1,c(0,TE[i]))%*%coef.rd%*%rbind(1,Ti)
        mui.c<-cbind(1,Ti,matrix(WE[i,],length(Ti),dw,byrow=TRUE))%*%coef.fix
        
        comp1i<-mui+COVVi%*%solve(VVi.c)%*%(Xi-mui.c)      ## E(X(0),X(t)|X(Ti)*,W)
        comp2i<-VVi-COVVi%*%solve(VVi.c)%*%t(COVVi)        ## Var(X(0),X(t)|X(Ti)*,W)
        
        } else {
        #### (ii) if case does not have X(Ti)* measure, i.e. case ID[i] not in id.long
        
       comp1i<-mui            
       comp2i<-VVi            
       } 
 
     comp12i<-as.numeric(t(c(parm1,parm2))%*%comp1i+t(c(parm1,parm2))%*%comp2i%*%c(parm1,parm2)/2)
    
    } 
    
    jj.r<-((i+1):(length(RE)))[RE[(i+1):length(RE)]<SE[i] & SE[(i+1):length(RE)]>=SE[i] & ID[(i+1):length(RE)]%in%id.long]
    jj.l<-NULL
    if (i>1) jj.l<-(1:(i-1))[RE[1:(i-1)]<SE[i] & SE[1:(i-1)]>=SE[i] & ID[1:(i-1)]%in%id.long]
    jj<-c(jj.l,jj.r) ## among risk set  
    
    jj.tie<-(1:(length(RE)))[RE[1:length(RE)]<SE[i] & ID[1:length(RE)]%in%id.long & SE[1:length(SE)]==SE[i] & delta[1:length(delta)]==1]
   
    jj.W.r<-((i+1):(length(RE)))[RE[(i+1):length(RE)]<SE[i] & ZE[(i+1):length(ZE)]==1 & SE[(i+1):length(RE)]>=SE[i] & !ID[(i+1):length(RE)]%in%id.long]
    jj.W.l<-NULL
    if (i>1) jj.W.l<-(1:(i-1))[RE[1:(i-1)]<SE[i] & ZE[1:(i-1)]==1 & !ID[1:(i-1)]%in%id.long & SE[1:(i-1)]>=SE[i]]
    jj.W<-c(jj.W.l,jj.W.r)    

    
    jj.W.tie<-(1:(length(RE)))[RE[1:length(RE)]<SE[i] & ZE[1:length(ZE)]==1 & SE[1:length(RE)]==SE[i] & !ID[1:length(RE)]%in%id.long  & delta[1:length(delta)]==1]
 
  
    jj.Z0.r<-((i+1):length(RE))[RE[(i+1):length(RE)]<SE[i] & ZE[(i+1):length(ZE)]==0 & SE[(i+1):length(SE)]>=SE[i]]  
    jj.Z0.l<-NULL
    if (i>1) jj.Z0.l<-(1:(i-1))[RE[1:(i-1)]<SE[i] & ZE[1:(i-1)]==0 & SE[1:(i-1)]>=SE[i]]
    jj.Z0=c(jj.Z0.l,jj.Z0.r)

    jj.Z0.tie<-(1:length(RE))[RE[1:length(RE)]<SE[i] & ZE[1:length(ZE)]==0 & SE[1:length(SE)]==SE[i] & delta[1:length(delta)]==1]  
  
    jj.1.r<-((i+1):(length(RE)))[RE[(i+1):length(RE)]<SE[i] & datwide.delta0[match(ID[(i+1):length(RE)],datwide.delta0$id),]$group==1 &
       datwide.delta0[match(ID[(i+1):length(RE)],datwide.delta0$id),]$groupall==1 & SE[(i+1):length(RE)]>=SE[i]]
    jj.1.l<-NULL
    if (i>1) jj.1.l<-(1:(i-1))[RE[1:(i-1)]<SE[i] & datwide.delta0[match(ID[1:(i-1)],datwide.delta0$id),]$group==1 &
       datwide.delta0[match(ID[1:(i-1)],datwide.delta0$id),]$groupall==1 & SE[1:(i-1)]>=SE[i]]
    jj.1<-c(jj.1.l,jj.1.r) ## among risk set

    
    jj.1.tie<-(1:(length(RE)))[RE[1:length(RE)]<SE[i] & datwide.delta0[match(ID[1:length(RE)],datwide.delta0$id),]$group==1 &
       datwide.delta0[match(ID[1:length(RE)],datwide.delta0$id),]$groupall==1 & SE[1:length(SE)]==SE[i] & delta[1:length(delta)]==1]
   
   
    #######
 
    jj.2.r<-((i+1):(length(RE)))[RE[(i+1):length(RE)]<SE[i] & datwide.delta0[match(ID[(i+1):length(RE)],datwide.delta0$id),]$group==2 & 
       datwide.delta0[match(ID[(i+1):length(RE)],datwide.delta0$id),]$groupall==2 & SE[(i+1):length(SE)]>=SE[i]]
    jj.2.l<-NULL
    if (i>1)   jj.2.l<-(1:(i-1))[RE[1:(i-1)]<SE[i] & datwide.delta0[match(ID[1:(i-1)],datwide.delta0$id),]$group==2 & 
       datwide.delta0[match(ID[1:(i-1)],datwide.delta0$id),]$groupall==2 & SE[1:(i-1)]>=SE[i]]
    jj.2<-c(jj.2.l,jj.2.r)

    jj.2.tie<-(1:(length(RE)))[RE[1:length(RE)]<SE[i] & datwide.delta0[match(ID[1:length(RE)],datwide.delta0$id),]$group==2 & 
       datwide.delta0[match(ID[1:length(RE)],datwide.delta0$id),]$groupall==2 & SE[1:length(SE)]==SE[i] & delta[1:length(delta)]==1]
   

    ###
    jj.3.r<-((i+1):(length(RE)))[RE[(i+1):length(RE)]<SE[i] & datwide.delta0[match(ID[(i+1):length(RE)],datwide.delta0$id),]$group==3 & 
       datwide.delta0[match(ID[(i+1):length(RE)],datwide.delta0$id),]$groupall==3 & SE[(i+1):length(SE)]>=SE[i]]
    jj.3.l<-NULL
    if (i>1) jj.3.l<-(1:(i-1))[RE[1:(i-1)]<SE[i] & datwide.delta0[match(ID[1:(i-1)],datwide.delta0$id),]$group==3 & 
       datwide.delta0[match(ID[1:(i-1)],datwide.delta0$id),]$groupall==3 & SE[1:(i-1)]>=SE[i]]
    jj.3<-c(jj.3.l,jj.3.r)
   

    jj.3.tie<-(1:(length(RE)))[RE[1:length(RE)]<SE[i] & datwide.delta0[match(ID[1:length(RE)],datwide.delta0$id),]$group==3 & 
       datwide.delta0[match(ID[1:length(RE)],datwide.delta0$id),]$groupall==3 & SE[1:length(SE)]==SE[i] & delta[1:length(delta)]==1]
  
     jj.1<-as.numeric(na.omit(jj.1))
     jj.2<-as.numeric(na.omit(jj.2))
     jj.3<-as.numeric(na.omit(jj.3))
   
       jj.1.tie<-as.numeric(na.omit(jj.1.tie))
     jj.2.tie<-as.numeric(na.omit(jj.2.tie))
     jj.3.tie<-as.numeric(na.omit(jj.3.tie))
     
    ss=ss+sum((RE[jj.4.full]<SE[i])*(SE[jj.4.full]>=SE[i])*(!jj.4.full==i)*weights[jj.4.full]*exp(out.jj.4.full[i,]-comp12i-parm3*(ZE[i]*TE[i])
    -parm0*ZE[i]-matrix(WEc[i,],length(i),dc)%*%parm4))  ## risk set from jj.4.full
     
    ss=ss-(ind.tie-1)/di*sum((RE[jj.4.full]<SE[i])*(SE[jj.4.full]==SE[i])*(!jj.4.full==i)*weights[jj.4.full]*exp(out.jj.4.full[i,]-comp12i-parm3*(ZE[i]*TE[i])
    -parm0*ZE[i]-matrix(WEc[i,],length(i),dc)%*%parm4))  ## tie from jj.4.full
     

   if (length(jj.W)>0){
      ### No X measure 
    mu0.EXE<-cbind(1,0,matrix(as.numeric(WE[jj.W,]),length(jj.W),dw))%*%coef.fix
    mu.EXE<-cbind(1,SE[i]-RE[jj.W],matrix(as.numeric(WE[jj.W,]),length(jj.W),dw))%*%coef.fix
  
   comp1.b<-parm1*mu0.EXE+parm2*mu.EXE
       
   comp2.b<-parm1^2*coef.rd[1,1]+2*parm1*parm2*(coef.rd[1,1]+(SE[i]-RE[jj.W])*coef.rd[1,2])+
   parm2^2*(coef.rd[1,1]+2*(SE[i]-RE[jj.W])*coef.rd[1,2]+(SE[i]-RE[jj.W])^2*coef.rd[2,2])
   
   ss=ss+sum(weights[jj.W]*exp(as.numeric(comp1.b)+0.5*comp2.b-comp12i+parm3*(SE[i]-RE[jj.W]-ZE[i]*TE[i])
   +parm0*(ZE[jj.W]-ZE[i])-t(as.numeric(WEc[i,])-matrix(as.numeric(t(WEc[jj.W,])),dc,length(jj.W)))%*%parm4) )  
    
   } 
   ### minus ties
   if (length(jj.W.tie)>0){
   
        mu0.EXE<-cbind(1,0,matrix(as.numeric(WE[jj.W.tie,]),length(jj.W.tie),dw))%*%coef.fix
        mu.EXE<-cbind(1,SE[i]-RE[jj.W.tie],matrix(as.numeric(WE[jj.W.tie,]),length(jj.W.tie),dw))%*%coef.fix
      
       comp1.b<-parm1*mu0.EXE+parm2*mu.EXE
           
       comp2.b<-parm1^2*coef.rd[1,1]+2*parm1*parm2*(coef.rd[1,1]+(SE[i]-RE[jj.W.tie])*coef.rd[1,2])+
       parm2^2*(coef.rd[1,1]+2*(SE[i]-RE[jj.W.tie])*coef.rd[1,2]+(SE[i]-RE[jj.W.tie])^2*coef.rd[2,2])
       
       ss=ss-(ind.tie-1)/di*sum(weights[jj.W.tie]*exp(as.numeric(comp1.b)+0.5*comp2.b-comp12i+parm3*(SE[i]-RE[jj.W.tie]-ZE[i]*TE[i])
       +parm0*(ZE[jj.W.tie]-ZE[i])-t(as.numeric(WEc[i,])-matrix(as.numeric(t(WEc[jj.W.tie,])),dc,length(jj.W.tie)))%*%parm4) ) 
   }      

 
   ### Measure of X(0)* only
   
    #account for different t levels
   
   if (length(jj.1)>0){
   Xt1=datwide.delta0[match(ID[jj.1],datwide.delta0$id),Vname[1]]
    t1=datwide.delta0[match(ID[jj.1],datwide.delta0$id),Tname[1]]
 
    ## make t1 vector
   COVV.a<-matrix(c(coef.rd[1,1]+t1*coef.rd[1,2],coef.rd[1,1]+(t1+SE[i]-RE[jj.1])*coef.rd[1,2]+t1*(SE[i]-RE[jj.1])*coef.rd[2,2]),length(jj.1),2,byrow=F)
   SS.a<-sapply(coef.rd[1,1]+2*t1*coef.rd[1,2]+t1^2*coef.rd[2,2]+sigma.e^2,solve)
   
   slope.a<-COVV.a*as.numeric(SS.a)
   Sigma1.2.a=outer(t(slope.a),(COVV.a),'*'); 
   kappa1=diag(as.matrix(Sigma1.2.a[1,,,1])); kappa2=diag(as.matrix(Sigma1.2.a[1,,,2])); 
   kappa3=diag(as.matrix(Sigma1.2.a[2,,,1])); kappa4=diag(as.matrix(Sigma1.2.a[2,,,2]))
   
     mu0.EXE<-cbind(1,0,matrix(as.numeric(WE[jj.1,]),length(jj.1),dw))%*%coef.fix
      mu.EXE<-cbind(1,SE[i]-RE[jj.1],matrix(as.numeric(WE[jj.1,]),length(jj.1),dw))%*%coef.fix
  
       comp1<-parm1*(mu0.EXE+slope.a[,1]*(Xt1-mu0.EXE))+
       parm2*(mu.EXE+slope.a[,2]*(Xt1-mu0.EXE))
       
       comp2<-parm1^2*(coef.rd[1,1]-kappa1)+
       parm2*parm1*(2*coef.rd[1,1]+2*(SE[i]-RE[jj.1])*coef.rd[1,2]-kappa2-kappa3)+
       parm2^2*(coef.rd[1,1]+(SE[i]-RE[jj.1])^2*coef.rd[2,2]+2*(SE[i]-RE[jj.1])*coef.rd[1,2]-kappa4)
    
    ss=ss+sum(weights[jj.1]*exp(as.numeric(comp1)+0.5*comp2-comp12i+parm3*(SE[i]-RE[jj.1]-ZE[i]*TE[i])
    +parm0*(ZE[jj.1]-ZE[i])-t(as.numeric(WEc[i,])-matrix(as.numeric(t(WEc[jj.1,])),dc,length(jj.1)))%*%parm4) )  
   
   rm(Xt1)
   }

   ## minus ties

   if (length(jj.1.tie)>0){
   Xt1=datwide.delta0[match(ID[jj.1.tie],datwide.delta0$id),Vname[1]]
   
   t1=datwide.delta0[match(ID[jj.1.tie],datwide.delta0$id),Tname[1]]
   
    # make t1 vector
   COVV.a<-matrix(c(coef.rd[1,1]+t1*coef.rd[1,2],coef.rd[1,1]+(t1+SE[i]-RE[jj.1.tie])*coef.rd[1,2]+t1*(SE[i]-RE[jj.1.tie])*coef.rd[2,2]),length(jj.1.tie),2,byrow=F)
   SS.a<-sapply(coef.rd[1,1]+2*t1*coef.rd[1,2]+t1^2*coef.rd[2,2]+sigma.e^2,solve)
 
   
   slope.a<-COVV.a*as.numeric(SS.a)
   Sigma1.2.a=outer(t(slope.a),(COVV.a),'*'); 
   kappa1=diag(as.matrix(Sigma1.2.a[1,,,1])); kappa2=diag(as.matrix(Sigma1.2.a[1,,,2])); 
   kappa3=diag(as.matrix(Sigma1.2.a[2,,,1])); kappa4=diag(as.matrix(Sigma1.2.a[2,,,2]))
   
     mu0.EXE<-cbind(1,0,matrix(as.numeric(WE[jj.1.tie,]),length(jj.1.tie),dw))%*%coef.fix
      mu.EXE<-cbind(1,SE[i]-RE[jj.1.tie],matrix(as.numeric(WE[jj.1.tie,]),length(jj.1.tie),dw))%*%coef.fix
  
       comp1<-parm1*(mu0.EXE+slope.a[,1]*(Xt1-mu0.EXE))+
       parm2*(mu.EXE+slope.a[,2]*(Xt1-mu0.EXE))
       
       comp2<-parm1^2*(coef.rd[1,1]-kappa1)+
       parm2*parm1*(2*coef.rd[1,1]+2*(SE[i]-RE[jj.1.tie])*coef.rd[1,2]-kappa2-kappa3)+
       parm2^2*(coef.rd[1,1]+(SE[i]-RE[jj.1.tie])^2*coef.rd[2,2]+2*(SE[i]-RE[jj.1.tie])*coef.rd[1,2]-kappa4)
    
    ss=ss-(ind.tie-1)/di*sum(weights[jj.1.tie]*exp(as.numeric(comp1)+0.5*comp2-comp12i+parm3*(SE[i]-RE[jj.1.tie]-ZE[i]*TE[i])
    +parm0*(ZE[jj.1.tie]-ZE[i])-t(as.numeric(WEc[i,])-matrix(as.numeric(t(WEc[jj.1.tie,])),dc,length(jj.1.tie)))%*%parm4 ) )  
   
   rm(Xt1)
   }



   ###### Measure of X(0)* and X(t2)* 
     if (length(jj.2)>0){
     Xt1=datwide.delta0[match(ID[jj.2],datwide.delta0$id),Vname[1]]
     Xt2=datwide.delta0[match(ID[jj.2],datwide.delta0$id),Vname[2]]
     
       ## make t1, t2 vectors
      t1=datwide.delta0[match(ID[jj.2],datwide.delta0$id),Tname[1]]
      t2=datwide.delta0[match(ID[jj.2],datwide.delta0$id),Tname[2]]
         
     COVV.b1<-coef.rd[1,1]+t1*coef.rd[1,2]
     COVV.b2<-coef.rd[1,1]+t2*coef.rd[1,2]
     COVV.b3<-coef.rd[1,1]+(t1+SE[i]-RE[jj.2])*coef.rd[1,2]+t1*(SE[i]-RE[jj.2])*coef.rd[2,2]
     COVV.b4<-coef.rd[1,1]+t2*(SE[i]-RE[jj.2])*coef.rd[2,2]+(t2+SE[i]-RE[jj.2])*coef.rd[1,2]
     

     SS.b.1<-coef.rd[1,1]+2*t1*coef.rd[1,2]+t1*2*coef.rd[2,2]+sigma.e^2
     SS.b.2<-coef.rd[1,1]+2*t2*coef.rd[1,2]+t2^2*coef.rd[2,2]+sigma.e^2
     SS.b.12<-coef.rd[1,1]+(t1+t2)*coef.rd[1,2]+t1*t2*coef.rd[2,2]
     
     SS.b<-array(NA,dim=c(2,2,length(jj.2)))
     SS.b[1,1,]<-SS.b.1
     SS.b[1,2,]<-SS.b[2,1,]<-SS.b.12
     SS.b[2,2,]<-SS.b.2
     
     SS.b<-apply(SS.b,3,solve)  ### SS.b is a 4 by length(jj.2) matrix with each column the inverse of a matrix
 
     ### slope
    
     ### When t1 t2 are vectors, SS.b[1,], SS.b[2,], SS.b[3,], SS.b[4,] 
     ### corresponds [1,1],[2,1],[1,2],[2,2] of a 2x2 matrix
     
     gamma1=COVV.b1*SS.b[1,]+COVV.b2*SS.b[2,]
     gamma2=COVV.b1*SS.b[3,]+COVV.b2*SS.b[4,]
     gamma3=COVV.b3*SS.b[1,]+COVV.b4*SS.b[2,]
     gamma4=COVV.b3*SS.b[3,]+COVV.b4*SS.b[4,]
     

     ### Sigma1.2.b
     eta1=gamma1*COVV.b1+gamma2*COVV.b2
     eta2=gamma1*COVV.b3+gamma2*COVV.b4
     eta3=gamma3*COVV.b1+gamma4*COVV.b2
     eta4=gamma3*COVV.b3+gamma4*COVV.b4
     
      mu0.EXE<-cbind(1,0,matrix(as.numeric(WE[jj.2,]),length(jj.2),dw))%*%coef.fix
       mu.EXE<-cbind(1,SE[i]-RE[jj.2],matrix(as.numeric(WE[jj.2,]),length(jj.2),dw))%*%coef.fix
       mut1.EXE<-cbind(1,t1,matrix(as.numeric(WE[jj.2,]),length(jj.2),dw))%*%coef.fix
       mut2.EXE<-cbind(1,t2,matrix(as.numeric(WE[jj.2,]),length(jj.2),dw))%*%coef.fix
       
      
      comp1<-parm1*(mu0.EXE+gamma1*Xt1+gamma2*Xt2-gamma1*mut1.EXE-gamma2*mut2.EXE)+
        parm2*(mu.EXE+gamma3*Xt1+gamma4*Xt2-gamma3*mut1.EXE-gamma4*mut2.EXE)
      
      comp2<-parm1^2*(coef.rd[1,1]-eta1)+parm2*parm1*(2*coef.rd[1,1]+2*(SE[i]-RE[jj.2])*coef.rd[1,2]-eta2-eta3)+
      parm2^2*(coef.rd[1,1]+(SE[i]-RE[jj.2])^2*coef.rd[2,2]+2*(SE[i]-RE[jj.2])*coef.rd[1,2]-eta4)
      
    
    ss=ss+sum(weights[jj.2]*exp(as.numeric(comp1)+0.5*comp2-comp12i+parm3*(SE[i]-RE[jj.2]-ZE[i]*TE[i])
    +parm0*(ZE[jj.2]-ZE[i])-t(as.numeric(WEc[i,])-matrix(as.numeric(t(WEc[jj.2,])),dc,length(jj.2)))%*%parm4) )  
    rm(Xt1,Xt2)
    }

    ## minus ties

     if (length(jj.2.tie)>0){
     Xt1=datwide.delta0[match(ID[jj.2.tie],datwide.delta0$id),Vname[1]]
     Xt2=datwide.delta0[match(ID[jj.2.tie],datwide.delta0$id),Vname[2]]
     
         ## make t1, t2 vectors
     t1=datwide.delta0[match(ID[jj.2.tie],datwide.delta0$id),Tname[1]]
     t2=datwide.delta0[match(ID[jj.2.tie],datwide.delta0$id),Tname[2]]
    

   ### when t1 t2 are vectors 
     COVV.b1<-coef.rd[1,1]+t1*coef.rd[1,2]
     COVV.b2<-coef.rd[1,1]+t2*coef.rd[1,2]
     COVV.b3<-coef.rd[1,1]+(t1+SE[i]-RE[jj.2.tie])*coef.rd[1,2]+t1*(SE[i]-RE[jj.2.tie])*coef.rd[2,2]
     COVV.b4<-coef.rd[1,1]+t2*(SE[i]-RE[jj.2.tie])*coef.rd[2,2]+(t2+SE[i]-RE[jj.2.tie])*coef.rd[1,2]
     

     SS.b.1<-coef.rd[1,1]+2*t1*coef.rd[1,2]+t1*2*coef.rd[2,2]+sigma.e^2
     SS.b.2<-coef.rd[1,1]+2*t2*coef.rd[1,2]+t2^2*coef.rd[2,2]+sigma.e^2
     SS.b.12<-coef.rd[1,1]+(t1+t2)*coef.rd[1,2]+t1*t2*coef.rd[2,2]
     
     SS.b<-array(NA,dim=c(2,2,length(jj.2.tie)))
     SS.b[1,1,]<-SS.b.1
     SS.b[1,2,]<-SS.b[2,1,]<-SS.b.12
     SS.b[2,2,]<-SS.b.2
     
     SS.b<-apply(SS.b,3,solve)  ### SS.b is a 4 by length(jj.2) matrix with each column the inverse of a matrix
 
     ### slope
    
     ### When t1 t2 are vectors , SS.b[1,], SS.b[2,], SS.b[3,], SS.b[4,] 
     ### corresponds [1,1],[2,1],[1,2],[2,2] of a 2x2 matrix
     
     gamma1=COVV.b1*SS.b[1,]+COVV.b2*SS.b[2,]
     gamma2=COVV.b1*SS.b[3,]+COVV.b2*SS.b[4,]
     gamma3=COVV.b3*SS.b[1,]+COVV.b4*SS.b[2,]
     gamma4=COVV.b3*SS.b[3,]+COVV.b4*SS.b[4,]
     
     ### Sigma1.2.b
     eta1=gamma1*COVV.b1+gamma2*COVV.b2
     eta2=gamma1*COVV.b3+gamma2*COVV.b4
     eta3=gamma3*COVV.b1+gamma4*COVV.b2
     eta4=gamma3*COVV.b3+gamma4*COVV.b4
     
      mu0.EXE<-cbind(1,0,matrix(as.numeric(WE[jj.2.tie,]),length(jj.2.tie),dw))%*%coef.fix
       mu.EXE<-cbind(1,SE[i]-RE[jj.2.tie],matrix(as.numeric(WE[jj.2.tie,]),length(jj.2.tie),dw))%*%coef.fix
       mut1.EXE<-cbind(1,t1,matrix(as.numeric(WE[jj.2.tie,]),length(jj.2.tie),dw))%*%coef.fix
       mut2.EXE<-cbind(1,t2,matrix(as.numeric(WE[jj.2.tie,]),length(jj.2.tie),dw))%*%coef.fix
       
      
      comp1<-parm1*(mu0.EXE+gamma1*Xt1+gamma2*Xt2-gamma1*mut1.EXE-gamma2*mut2.EXE)+
        parm2*(mu.EXE+gamma3*Xt1+gamma4*Xt2-gamma3*mut1.EXE-gamma4*mut2.EXE)
      
      comp2<-parm1^2*(coef.rd[1,1]-eta1)+parm2*parm1*(2*coef.rd[1,1]+2*(SE[i]-RE[jj.2.tie])*coef.rd[1,2]-eta2-eta3)+
      parm2^2*(coef.rd[1,1]+(SE[i]-RE[jj.2.tie])^2*coef.rd[2,2]+2*(SE[i]-RE[jj.2.tie])*coef.rd[1,2]-eta4)
      
    
    ss=ss-(ind.tie-1)/di*sum(weights[jj.2.tie]*exp(as.numeric(comp1)+0.5*comp2-comp12i+parm3*(SE[i]-RE[jj.2.tie]-ZE[i]*TE[i])
    +parm0*(ZE[jj.2.tie]-ZE[i])-t(as.numeric(WEc[i,])-matrix(as.numeric(t(WEc[jj.2.tie,])),dc,length(jj.2.tie)))%*%parm4) )  
    rm(Xt1,Xt2)
    }
    #####  ###### Measure of X(t1)*, X(t2)*, and X(t3)*
    
     if (length(jj.3)>0){
     Xt1=datwide.delta0[match(ID[jj.3],datwide.delta0$id),Vname[1]]
     Xt2=datwide.delta0[match(ID[jj.3],datwide.delta0$id),Vname[2]]
     Xt3=datwide.delta0[match(ID[jj.3],datwide.delta0$id),Vname[3]]
     
       ### t1, t2, t3 vector
     t1=datwide.delta0[match(ID[jj.3],datwide.delta0$id),Tname[1]]
     t2=datwide.delta0[match(ID[jj.3],datwide.delta0$id),Tname[2]]
     t3=datwide.delta0[match(ID[jj.3],datwide.delta0$id),Tname[3]]
     
   ####  t1,t2,t3 vector

     COVV.c1=coef.rd[1,1]+t1*coef.rd[1,2]
     COVV.c2=coef.rd[1,1]+t2*coef.rd[1,2]
     COVV.c3=coef.rd[1,1]+t3*coef.rd[1,2]
     COVV.c4=coef.rd[1,1]+(t1+SE[i]-RE[jj.3])*coef.rd[1,2]+t1*(SE[i]-RE[jj.3])*coef.rd[2,2]
     COVV.c5=coef.rd[1,1]+(t2+SE[i]-RE[jj.3])*coef.rd[1,2]+t2*(SE[i]-RE[jj.3])*coef.rd[2,2]
     COVV.c6=coef.rd[1,1]+(t3+SE[i]-RE[jj.3])*coef.rd[1,2]+t3*(SE[i]-RE[jj.3])*coef.rd[2,2]
     
     ##########
    
     SS.c.11<-coef.rd[1,1]+2*t1*coef.rd[1,2]+t1^2*coef.rd[2,2]+sigma.e^2
     SS.c.22<-coef.rd[1,1]+2*t2*coef.rd[1,2]+t2^2*coef.rd[2,2]+sigma.e^2
     SS.c.33<-coef.rd[1,1]+2*t3*coef.rd[1,2]+t3^2*coef.rd[2,2]+sigma.e^2
     SS.c.12=coef.rd[1,1]+(t1+t2)*coef.rd[1,2]+t1*t2*coef.rd[2,2]
     SS.c.13=coef.rd[1,1]+(t1+t3)*coef.rd[1,2]+t1*t3*coef.rd[2,2]
     SS.c.23=coef.rd[1,1]+(t2+t3)*coef.rd[1,2]+t2*t3*coef.rd[2,2]
     
     SS.c<-array(NA,dim=c(3,3,length(jj.3)))
     SS.c[1,1,]<-SS.c.11
     SS.c[2,2,]<-SS.c.22
     SS.c[3,3,]<-SS.c.33
     SS.c[1,2,]<-SS.c[2,1,]<-SS.c.12
     SS.c[1,3,]<-SS.c[3,1,]<-SS.c.13
     SS.c[2,3,]<-SS.c[3,2,]<-SS.c.23
     
    SS.c<-apply(SS.c,3,solve)
    

     #### SS.c is 9xlength(jj.3), elements in each column corresponds to [1,1],[2,1],[3,1],[1,2],[2,2],[3,2],[1,3],[2,3],[3,3] of a matrix

     delta1=COVV.c1*SS.c[1,]+COVV.c2*SS.c[2,]+COVV.c3*SS.c[3,]
     delta2=COVV.c1*SS.c[4,]+COVV.c2*SS.c[5,]+COVV.c3*SS.c[6,]
     delta3=COVV.c1*SS.c[7,]+COVV.c2*SS.c[8,]+COVV.c3*SS.c[9,]
     delta4=COVV.c4*SS.c[1,]+COVV.c5*SS.c[2,]+COVV.c6*SS.c[3,]
     delta5=COVV.c4*SS.c[4,]+COVV.c5*SS.c[5,]+COVV.c6*SS.c[6,]
     delta6=COVV.c4*SS.c[7,]+COVV.c5*SS.c[8,]+COVV.c6*SS.c[9,]

     xi1=delta1*COVV.c1+delta2*COVV.c2+delta3*COVV.c3
     xi2=delta1*COVV.c4+delta2*COVV.c5+delta3*COVV.c6
     xi3=delta4*COVV.c1+delta5*COVV.c2+delta6*COVV.c3
     xi4=delta4*COVV.c4+delta5*COVV.c5+delta6*COVV.c6
      
      mu0.EXE<-cbind(1,0,matrix(as.numeric(WE[jj.3,]),length(jj.3),dw))%*%coef.fix
       mu.EXE<-cbind(1,SE[i]-RE[jj.3],matrix(as.numeric(WE[jj.3,]),length(jj.3),dw))%*%coef.fix
       mut1.EXE<-cbind(1,t1,matrix(as.numeric(WE[jj.3,]),length(jj.3),dw))%*%coef.fix
       mut2.EXE<-cbind(1,t2,matrix(as.numeric(WE[jj.3,]),length(jj.3),dw))%*%coef.fix
       mut3.EXE<-cbind(1,t3,matrix(as.numeric(WE[jj.3,]),length(jj.3),dw))%*%coef.fix
       
      
      comp1<-parm1*(mu0.EXE+delta1*Xt1+delta2*Xt2+delta3*Xt3-delta1*mut1.EXE-delta2*mut2.EXE-delta3*mut3.EXE)+
        parm2*(mu.EXE+delta4*Xt1+delta5*Xt2+delta6*Xt3-delta4*mut1.EXE-delta5*mut2.EXE-delta6*mut3.EXE)
      
      comp2<-parm1^2*(coef.rd[1,1]-xi1)+parm2*parm1*(2*coef.rd[1,1]+2*(SE[i]-RE[jj.3])*coef.rd[1,2]-xi2-xi3)+
      parm2^2*(coef.rd[1,1]+(SE[i]-RE[jj.3])^2*coef.rd[2,2]+2*(SE[i]-RE[jj.3])*coef.rd[1,2]-xi4)
      
    
    ss=ss+sum(weights[jj.3]*exp(as.numeric(comp1)+0.5*comp2-comp12i+parm3*(SE[i]-RE[jj.3]-ZE[i]*TE[i])
    +parm0*(ZE[jj.3]-ZE[i])-t(as.numeric(WEc[i,])-matrix(as.numeric(t(WEc[jj.3,])),dc,length(jj.3)))%*%parm4) )  
   
    rm(Xt1,Xt2,Xt3)
    }

    #### Minus ties
 
     if (length(jj.3.tie)>0){
     Xt1=datwide.delta0[match(ID[jj.3.tie],datwide.delta0$id),Vname[1]]
     Xt2=datwide.delta0[match(ID[jj.3.tie],datwide.delta0$id),Vname[2]]
     Xt3=datwide.delta0[match(ID[jj.3.tie],datwide.delta0$id),Vname[3]]
     
      ### t1, t2, t3 vector
     t1=datwide.delta0[match(ID[jj.3.tie],datwide.delta0$id),Tname[1]]
     t2=datwide.delta0[match(ID[jj.3.tie],datwide.delta0$id),Tname[2]]
     t3=datwide.delta0[match(ID[jj.3.tie],datwide.delta0$id),Tname[3]]
 
     COVV.c1=coef.rd[1,1]+t1*coef.rd[1,2]
     COVV.c2=coef.rd[1,1]+t2*coef.rd[1,2]
     COVV.c3=coef.rd[1,1]+t3*coef.rd[1,2]
     COVV.c4=coef.rd[1,1]+(t1+SE[i]-RE[jj.3.tie])*coef.rd[1,2]+t1*(SE[i]-RE[jj.3.tie])*coef.rd[2,2]
     COVV.c5=coef.rd[1,1]+(t2+SE[i]-RE[jj.3.tie])*coef.rd[1,2]+t2*(SE[i]-RE[jj.3.tie])*coef.rd[2,2]
     COVV.c6=coef.rd[1,1]+(t3+SE[i]-RE[jj.3.tie])*coef.rd[1,2]+t3*(SE[i]-RE[jj.3.tie])*coef.rd[2,2]
     
     ##########
    
     SS.c.11<-coef.rd[1,1]+2*t1*coef.rd[1,2]+t1^2*coef.rd[2,2]+sigma.e^2
     SS.c.22<-coef.rd[1,1]+2*t2*coef.rd[1,2]+t2^2*coef.rd[2,2]+sigma.e^2
     SS.c.33<-coef.rd[1,1]+2*t3*coef.rd[1,2]+t3^2*coef.rd[2,2]+sigma.e^2
     SS.c.12=coef.rd[1,1]+(t1+t2)*coef.rd[1,2]+t1*t2*coef.rd[2,2]
     SS.c.13=coef.rd[1,1]+(t1+t3)*coef.rd[1,2]+t1*t3*coef.rd[2,2]
     SS.c.23=coef.rd[1,1]+(t2+t3)*coef.rd[1,2]+t2*t3*coef.rd[2,2]
     
     SS.c<-array(NA,dim=c(3,3,length(jj.3.tie)))
     SS.c[1,1,]<-SS.c.11
     SS.c[2,2,]<-SS.c.22
     SS.c[3,3,]<-SS.c.33
     SS.c[1,2,]<-SS.c[2,1,]<-SS.c.12
     SS.c[1,3,]<-SS.c[3,1,]<-SS.c.13
     SS.c[2,3,]<-SS.c[3,2,]<-SS.c.23
         
    SS.c<-apply(SS.c,3,solve)
    

     #### SS.c is 9xlength(jj.3), elements in each column corresponds to [1,1],[2,1],[3,1],[1,2],[2,2],[3,2],[1,3],[2,3],[3,3] of a matrix

     delta1=COVV.c1*SS.c[1,]+COVV.c2*SS.c[2,]+COVV.c3*SS.c[3,]
     delta2=COVV.c1*SS.c[4,]+COVV.c2*SS.c[5,]+COVV.c3*SS.c[6,]
     delta3=COVV.c1*SS.c[7,]+COVV.c2*SS.c[8,]+COVV.c3*SS.c[9,]
     delta4=COVV.c4*SS.c[1,]+COVV.c5*SS.c[2,]+COVV.c6*SS.c[3,]
     delta5=COVV.c4*SS.c[4,]+COVV.c5*SS.c[5,]+COVV.c6*SS.c[6,]
     delta6=COVV.c4*SS.c[7,]+COVV.c5*SS.c[8,]+COVV.c6*SS.c[9,]

     
     xi1=delta1*COVV.c1+delta2*COVV.c2+delta3*COVV.c3
     xi2=delta1*COVV.c4+delta2*COVV.c5+delta3*COVV.c6
     xi3=delta4*COVV.c1+delta5*COVV.c2+delta6*COVV.c3
     xi4=delta4*COVV.c4+delta5*COVV.c5+delta6*COVV.c6
      
      mu0.EXE<-cbind(1,0,matrix(as.numeric(WE[jj.3.tie,]),length(jj.3.tie),dw))%*%coef.fix
       mu.EXE<-cbind(1,SE[i]-RE[jj.3.tie],matrix(as.numeric(WE[jj.3.tie,]),length(jj.3.tie),dw))%*%coef.fix
       mut1.EXE<-cbind(1,t1,matrix(as.numeric(WE[jj.3.tie,]),length(jj.3.tie),dw))%*%coef.fix
       mut2.EXE<-cbind(1,t2,matrix(as.numeric(WE[jj.3.tie,]),length(jj.3.tie),dw))%*%coef.fix
       mut3.EXE<-cbind(1,t3,matrix(as.numeric(WE[jj.3.tie,]),length(jj.3.tie),dw))%*%coef.fix
            
      comp1<-parm1*(mu0.EXE+delta1*Xt1+delta2*Xt2+delta3*Xt3-delta1*mut1.EXE-delta2*mut2.EXE-delta3*mut3.EXE)+
        parm2*(mu.EXE+delta4*Xt1+delta5*Xt2+delta6*Xt3-delta4*mut1.EXE-delta5*mut2.EXE-delta6*mut3.EXE)
      
      comp2<-parm1^2*(coef.rd[1,1]-xi1)+parm2*parm1*(2*coef.rd[1,1]+2*(SE[i]-RE[jj.3.tie])*coef.rd[1,2]-xi2-xi3)+
      parm2^2*(coef.rd[1,1]+(SE[i]-RE[jj.3.tie])^2*coef.rd[2,2]+2*(SE[i]-RE[jj.3.tie])*coef.rd[1,2]-xi4)
      
    ss=ss-(ind.tie-1)/di*sum(weights[jj.3.tie]*exp(as.numeric(comp1)+0.5*comp2-comp12i+parm3*(SE[i]-RE[jj.3.tie]-ZE[i]*TE[i])
    +parm0*(ZE[jj.3.tie]-ZE[i])-t(as.numeric(WEc[i,])-matrix(as.numeric(t(WEc[jj.3.tie,])),dc,length(jj.3.tie)))%*%parm4) )  
   
    rm(Xt1,Xt2,Xt3)
    }

    ss.new=ss
    ### placebo
   
   ### need to use jj.Z0 and jj.Z0.tie when there are confounders in the risk model
   if (length(jj.Z0)>0)
   ss=ss+sum(weights[jj.Z0]*exp(-comp12i-parm3*(ZE[i]*TE[i])-parm0*(ZE[i])
   -t(as.numeric(WEc[i,])-matrix(as.numeric(t(WEc[jj.Z0,])),dc,length(jj.Z0)))%*%parm4))
   
   if (length(jj.Z0.tie)>0)
   ss=ss-(ind.tie-1)/di*sum(weights[jj.Z0.tie]*exp(-comp12i-parm3*(ZE[i]*TE[i])-parm0*(ZE[i])-t(as.numeric(WEc[i,])-matrix(as.numeric(t(WEc[jj.Z0.tie,])),dc,length(jj.Z0.tie)))%*%parm4))
   
    out=out+weights.bar*log(weights[i]+ss)+(weights.bar-weights[i])*(comp12i+parm3*(SE[i]-RE[i])*ZE[i]+parm0*(ZE[i])
    +matrix(WEc[i,],length(i),dc)%*%parm4)
    
    outi=c(outi,ss)
 
    }
 return(out)
 }

########

######## Fit a Cox hazard model on calendar time with coefficient for Z, Z*x(0), Z*x(t), and W

nloglik.M1.b.Cal.12<- function(parm){

 out=0
 
 outi<-NULL
   ### coef for Z, Z*x(0), Z*x(t), and confounders respectively
    parm0=parm[1]
    parm1=parm[2]
    parm2=parm[3]
    #parm3=parm[2]
    parm4=parm[4]
   
   jj.4.full<-(1:length(RE))[datwide.delta0[match(ID,datwide.delta0$id),]$groupall > datwide.delta0[match(ID,datwide.delta0$id),]$group]
   jj.4.full<-as.numeric(na.omit(jj.4.full))
   out.jj.4.full<-NULL
   for (j in jj.4.full){
      
    keep=1:ll
 
    VVj.11<-coef.rd[1,1]
    VVj.12<-VVj.21<-coef.rd[1,1]+(SE[keep]-RE[j])*coef.rd[1,2]
    VVj.22<-coef.rd[1,1]+2*(SE[keep]-RE[j])*coef.rd[1,2]+(SE[keep]-RE[j])^2*coef.rd[2,2]
    
    muj.1=cbind(1,0,matrix(WE[j,],length(j),dw))%*%coef.fix
    muj.2=cbind(1,SE[keep]-RE[j],matrix(WE[j,],length(j),dw))%*%coef.fix
    
   
    oo.j<-id.long==ID[j]
    Xj<-X.long[oo.j]
    Tj<-T.long[oo.j]
    
    
    VVj.c<-cbind(1,Tj)%*%coef.rd%*%rbind(1,Tj)+diag(rep(sigma.e^2,length(Tj)),nrow=length(Tj),ncol=length(Tj))               
   
    COVVj.1<-coef.rd[1,1]+coef.rd[1,2]*Tj
    COVVj.2<-coef.rd[1,1]+coef.rd[1,2]*outer(SE[keep]-RE[j],Tj,'+')+coef.rd[2,2]*outer(SE[keep]-RE[j],Tj,'*')
    COVVj.1=matrix(COVVj.1,1,length(COVVj.1))
    if (length(Tj)==1) COVVj.2=matrix(COVVj.2,length(COVVj.2),1)
    
    muj.c<-cbind(1,Tj,matrix(WE[j,],length(j),dw))%*%coef.fix
        
    comp1j<-cbind(as.numeric(muj.1),muj.2)+cbind(as.numeric(COVVj.1%*%solve(VVj.c)%*%(Xj-muj.c)),COVVj.2%*%solve(VVj.c)%*%(Xj-muj.c))
    
     VVjminus.11<-COVVj.1%*%solve(VVj.c)%*%(t(COVVj.1))
     VVjminus.12<-COVVj.2%*%solve(VVj.c)%*%(t(COVVj.1))
     VVjminus.22<-diag(COVVj.2%*%solve(VVj.c)%*%(t(COVVj.2)))
     
     comp2j.11<-VVj.11-VVjminus.11
     comp2j.12<-VVj.12-VVjminus.12
     comp2j.22<-VVj.22-VVjminus.22
   
     out.jj.4.full<-cbind(out.jj.4.full,as.numeric(c(parm1,parm2)%*%t(comp1j))+parm1^2*as.numeric(comp2j.11)/2+parm1*parm2*comp2j.12+parm2^2*comp2j.22/2+parm3*(SE[keep]-RE[j])
     +parm0*ZE[j]+matrix(WEc[j,],length(j),dc)%*%parm4)
     
    }
        
    S.tie=SE.u[1]  ### unique SE value
    ind.tie=0       ### position within the unique value
    
    for (i in 1:ll){
    
    ## track #of ties
    
        di<-sum(SE==SE[i] & delta==1)
        
        if (SE[i]==S.tie){  ## in the same tie as previous value
          ind.tie=ind.tie+1
        } else { ## go to a new SE value
           S.tie=SE[i]
           ind.tie=1
        }
     ### compute average weight among ties
        weights.bar=mean(weights[SE==SE[i] & delta==1])
     
    ss=0
    comp12i=0 ### placebo contribute nothing
    if (ZE[i]==1) { 
    VVi<-cbind(1,c(0,TE[i]))%*%coef.rd%*%rbind(1,c(0,TE[i]))
    mui<-rbind(c(1,0,matrix(WE[i,],length(i),dw)),c(1,TE[i],matrix(WE[i,],length(i),dw)))%*%coef.fix
   
    oo.i<-id.long==ID[i]
        
       #### (i) if case has at least one X(Ti)* measured, i.e. case ID[i] belong to id.long
        if (sum(oo.i)>0){ 
        Xi<-X.long[oo.i]
        Ti<-T.long[oo.i]
        
        VVi.c<-cbind(1,Ti)%*%coef.rd%*%rbind(1,Ti)+diag(rep(sigma.e^2,length(Ti)),nrow=length(Ti),ncol=length(Ti))               
        COVVi<-cbind(1,c(0,TE[i]))%*%coef.rd%*%rbind(1,Ti)
        mui.c<-cbind(1,Ti,matrix(WE[i,],length(Ti),dw,byrow=TRUE))%*%coef.fix
        
        comp1i<-mui+COVVi%*%solve(VVi.c)%*%(Xi-mui.c)      ## E(X(0),X(t)|X(Ti)*,W)
        comp2i<-VVi-COVVi%*%solve(VVi.c)%*%t(COVVi)        ## Var(X(0),X(t)|X(Ti)*,W)
        
        } else {
        #### (ii) if case does not have X(Ti)* measure, i.e. case ID[i] not in id.long
        
       comp1i<-mui            
       comp2i<-VVi            
       } 
 
     comp12i<-as.numeric(t(c(parm1,parm2))%*%comp1i+t(c(parm1,parm2))%*%comp2i%*%c(parm1,parm2)/2)
    
    } 
    
    jj.r<-((i+1):(length(RE)))[RE[(i+1):length(RE)]<SE[i] & SE[(i+1):length(RE)]>=SE[i] & ID[(i+1):length(RE)]%in%id.long]
    jj.l<-NULL
    if (i>1) jj.l<-(1:(i-1))[RE[1:(i-1)]<SE[i] & SE[1:(i-1)]>=SE[i] & ID[1:(i-1)]%in%id.long]
    jj<-c(jj.l,jj.r) ## among risk set  
    
    jj.tie<-(1:(length(RE)))[RE[1:length(RE)]<SE[i] & ID[1:length(RE)]%in%id.long & SE[1:length(SE)]==SE[i] & delta[1:length(delta)]==1]
   
    jj.W.r<-((i+1):(length(RE)))[RE[(i+1):length(RE)]<SE[i] & ZE[(i+1):length(ZE)]==1 & SE[(i+1):length(RE)]>=SE[i] & !ID[(i+1):length(RE)]%in%id.long]
    jj.W.l<-NULL
    if (i>1) jj.W.l<-(1:(i-1))[RE[1:(i-1)]<SE[i] & ZE[1:(i-1)]==1 & !ID[1:(i-1)]%in%id.long & SE[1:(i-1)]>=SE[i]]
    jj.W<-c(jj.W.l,jj.W.r)    

    
    jj.W.tie<-(1:(length(RE)))[RE[1:length(RE)]<SE[i] & ZE[1:length(ZE)]==1 & SE[1:length(RE)]==SE[i] & !ID[1:length(RE)]%in%id.long  & delta[1:length(delta)]==1]
 
  
    jj.Z0.r<-((i+1):length(RE))[RE[(i+1):length(RE)]<SE[i] & ZE[(i+1):length(ZE)]==0 & SE[(i+1):length(SE)]>=SE[i]]  
    jj.Z0.l<-NULL
    if (i>1) jj.Z0.l<-(1:(i-1))[RE[1:(i-1)]<SE[i] & ZE[1:(i-1)]==0 & SE[1:(i-1)]>=SE[i]]
    jj.Z0=c(jj.Z0.l,jj.Z0.r)

    jj.Z0.tie<-(1:length(RE))[RE[1:length(RE)]<SE[i] & ZE[1:length(ZE)]==0 & SE[1:length(SE)]==SE[i] & delta[1:length(delta)]==1]  
  
    jj.1.r<-((i+1):(length(RE)))[RE[(i+1):length(RE)]<SE[i] & datwide.delta0[match(ID[(i+1):length(RE)],datwide.delta0$id),]$group==1 &
       datwide.delta0[match(ID[(i+1):length(RE)],datwide.delta0$id),]$groupall==1 & SE[(i+1):length(RE)]>=SE[i]]
    jj.1.l<-NULL
    if (i>1) jj.1.l<-(1:(i-1))[RE[1:(i-1)]<SE[i] & datwide.delta0[match(ID[1:(i-1)],datwide.delta0$id),]$group==1 &
       datwide.delta0[match(ID[1:(i-1)],datwide.delta0$id),]$groupall==1 & SE[1:(i-1)]>=SE[i]]
    jj.1<-c(jj.1.l,jj.1.r) ## among risk set

    
    jj.1.tie<-(1:(length(RE)))[RE[1:length(RE)]<SE[i] & datwide.delta0[match(ID[1:length(RE)],datwide.delta0$id),]$group==1 &
       datwide.delta0[match(ID[1:length(RE)],datwide.delta0$id),]$groupall==1 & SE[1:length(SE)]==SE[i] & delta[1:length(delta)]==1]
   
   
    #######
 
    jj.2.r<-((i+1):(length(RE)))[RE[(i+1):length(RE)]<SE[i] & datwide.delta0[match(ID[(i+1):length(RE)],datwide.delta0$id),]$group==2 & 
       datwide.delta0[match(ID[(i+1):length(RE)],datwide.delta0$id),]$groupall==2 & SE[(i+1):length(SE)]>=SE[i]]
    jj.2.l<-NULL
    if (i>1)   jj.2.l<-(1:(i-1))[RE[1:(i-1)]<SE[i] & datwide.delta0[match(ID[1:(i-1)],datwide.delta0$id),]$group==2 & 
       datwide.delta0[match(ID[1:(i-1)],datwide.delta0$id),]$groupall==2 & SE[1:(i-1)]>=SE[i]]
    jj.2<-c(jj.2.l,jj.2.r)

    jj.2.tie<-(1:(length(RE)))[RE[1:length(RE)]<SE[i] & datwide.delta0[match(ID[1:length(RE)],datwide.delta0$id),]$group==2 & 
       datwide.delta0[match(ID[1:length(RE)],datwide.delta0$id),]$groupall==2 & SE[1:length(SE)]==SE[i] & delta[1:length(delta)]==1]
   

    ###
    jj.3.r<-((i+1):(length(RE)))[RE[(i+1):length(RE)]<SE[i] & datwide.delta0[match(ID[(i+1):length(RE)],datwide.delta0$id),]$group==3 & 
       datwide.delta0[match(ID[(i+1):length(RE)],datwide.delta0$id),]$groupall==3 & SE[(i+1):length(SE)]>=SE[i]]
    jj.3.l<-NULL
    if (i>1) jj.3.l<-(1:(i-1))[RE[1:(i-1)]<SE[i] & datwide.delta0[match(ID[1:(i-1)],datwide.delta0$id),]$group==3 & 
       datwide.delta0[match(ID[1:(i-1)],datwide.delta0$id),]$groupall==3 & SE[1:(i-1)]>=SE[i]]
    jj.3<-c(jj.3.l,jj.3.r)
   

    jj.3.tie<-(1:(length(RE)))[RE[1:length(RE)]<SE[i] & datwide.delta0[match(ID[1:length(RE)],datwide.delta0$id),]$group==3 & 
       datwide.delta0[match(ID[1:length(RE)],datwide.delta0$id),]$groupall==3 & SE[1:length(SE)]==SE[i] & delta[1:length(delta)]==1]
  
     jj.1<-as.numeric(na.omit(jj.1))
     jj.2<-as.numeric(na.omit(jj.2))
     jj.3<-as.numeric(na.omit(jj.3))
   
       jj.1.tie<-as.numeric(na.omit(jj.1.tie))
     jj.2.tie<-as.numeric(na.omit(jj.2.tie))
     jj.3.tie<-as.numeric(na.omit(jj.3.tie))
     
    ss=ss+sum((RE[jj.4.full]<SE[i])*(SE[jj.4.full]>=SE[i])*(!jj.4.full==i)*weights[jj.4.full]*exp(out.jj.4.full[i,]-comp12i-parm3*(ZE[i]*TE[i])
    -parm0*ZE[i]-matrix(WEc[i,],length(i),dc)%*%parm4))  ## risk set from jj.4.full
     
    ss=ss-(ind.tie-1)/di*sum((RE[jj.4.full]<SE[i])*(SE[jj.4.full]==SE[i])*(!jj.4.full==i)*weights[jj.4.full]*exp(out.jj.4.full[i,]-comp12i-parm3*(ZE[i]*TE[i])
    -parm0*ZE[i]-matrix(WEc[i,],length(i),dc)%*%parm4))  ## tie from jj.4.full
     

   if (length(jj.W)>0){
      ### No X measure 
    mu0.EXE<-cbind(1,0,matrix(as.numeric(WE[jj.W,]),length(jj.W),dw))%*%coef.fix
    mu.EXE<-cbind(1,SE[i]-RE[jj.W],matrix(as.numeric(WE[jj.W,]),length(jj.W),dw))%*%coef.fix
  
   comp1.b<-parm1*mu0.EXE+parm2*mu.EXE
       
   comp2.b<-parm1^2*coef.rd[1,1]+2*parm1*parm2*(coef.rd[1,1]+(SE[i]-RE[jj.W])*coef.rd[1,2])+
   parm2^2*(coef.rd[1,1]+2*(SE[i]-RE[jj.W])*coef.rd[1,2]+(SE[i]-RE[jj.W])^2*coef.rd[2,2])
   
   ss=ss+sum(weights[jj.W]*exp(as.numeric(comp1.b)+0.5*comp2.b-comp12i+parm3*(SE[i]-RE[jj.W]-ZE[i]*TE[i])
   +parm0*(ZE[jj.W]-ZE[i])-t(as.numeric(WEc[i,])-matrix(as.numeric(t(WEc[jj.W,])),dc,length(jj.W)))%*%parm4) )  
    
   } 
   ### minus ties
   if (length(jj.W.tie)>0){
   
        mu0.EXE<-cbind(1,0,matrix(as.numeric(WE[jj.W.tie,]),length(jj.W.tie),dw))%*%coef.fix
        mu.EXE<-cbind(1,SE[i]-RE[jj.W.tie],matrix(as.numeric(WE[jj.W.tie,]),length(jj.W.tie),dw))%*%coef.fix
      
       comp1.b<-parm1*mu0.EXE+parm2*mu.EXE
           
       comp2.b<-parm1^2*coef.rd[1,1]+2*parm1*parm2*(coef.rd[1,1]+(SE[i]-RE[jj.W.tie])*coef.rd[1,2])+
       parm2^2*(coef.rd[1,1]+2*(SE[i]-RE[jj.W.tie])*coef.rd[1,2]+(SE[i]-RE[jj.W.tie])^2*coef.rd[2,2])
       
       ss=ss-(ind.tie-1)/di*sum(weights[jj.W.tie]*exp(as.numeric(comp1.b)+0.5*comp2.b-comp12i+parm3*(SE[i]-RE[jj.W.tie]-ZE[i]*TE[i])
       +parm0*(ZE[jj.W.tie]-ZE[i])-t(as.numeric(WEc[i,])-matrix(as.numeric(t(WEc[jj.W.tie,])),dc,length(jj.W.tie)))%*%parm4) ) 
   }      

 
   ### Measure of X(0)* only
   
    #account for different t levels
   
   if (length(jj.1)>0){
   Xt1=datwide.delta0[match(ID[jj.1],datwide.delta0$id),Vname[1]]
    t1=datwide.delta0[match(ID[jj.1],datwide.delta0$id),Tname[1]]
 
    ## make t1 vector
   COVV.a<-matrix(c(coef.rd[1,1]+t1*coef.rd[1,2],coef.rd[1,1]+(t1+SE[i]-RE[jj.1])*coef.rd[1,2]+t1*(SE[i]-RE[jj.1])*coef.rd[2,2]),length(jj.1),2,byrow=F)
   SS.a<-sapply(coef.rd[1,1]+2*t1*coef.rd[1,2]+t1^2*coef.rd[2,2]+sigma.e^2,solve)
   
   slope.a<-COVV.a*as.numeric(SS.a)
   Sigma1.2.a=outer(t(slope.a),(COVV.a),'*'); 
   kappa1=diag(as.matrix(Sigma1.2.a[1,,,1])); kappa2=diag(as.matrix(Sigma1.2.a[1,,,2])); 
   kappa3=diag(as.matrix(Sigma1.2.a[2,,,1])); kappa4=diag(as.matrix(Sigma1.2.a[2,,,2]))
   
     mu0.EXE<-cbind(1,0,matrix(as.numeric(WE[jj.1,]),length(jj.1),dw))%*%coef.fix
      mu.EXE<-cbind(1,SE[i]-RE[jj.1],matrix(as.numeric(WE[jj.1,]),length(jj.1),dw))%*%coef.fix
  
       comp1<-parm1*(mu0.EXE+slope.a[,1]*(Xt1-mu0.EXE))+
       parm2*(mu.EXE+slope.a[,2]*(Xt1-mu0.EXE))
       
       comp2<-parm1^2*(coef.rd[1,1]-kappa1)+
       parm2*parm1*(2*coef.rd[1,1]+2*(SE[i]-RE[jj.1])*coef.rd[1,2]-kappa2-kappa3)+
       parm2^2*(coef.rd[1,1]+(SE[i]-RE[jj.1])^2*coef.rd[2,2]+2*(SE[i]-RE[jj.1])*coef.rd[1,2]-kappa4)
    
    ss=ss+sum(weights[jj.1]*exp(as.numeric(comp1)+0.5*comp2-comp12i+parm3*(SE[i]-RE[jj.1]-ZE[i]*TE[i])
    +parm0*(ZE[jj.1]-ZE[i])-t(as.numeric(WEc[i,])-matrix(as.numeric(t(WEc[jj.1,])),dc,length(jj.1)))%*%parm4) )  
   
   rm(Xt1)
   }

   ## minus ties

   if (length(jj.1.tie)>0){
   Xt1=datwide.delta0[match(ID[jj.1.tie],datwide.delta0$id),Vname[1]]
   
   t1=datwide.delta0[match(ID[jj.1.tie],datwide.delta0$id),Tname[1]]
   
    # make t1 vector
   COVV.a<-matrix(c(coef.rd[1,1]+t1*coef.rd[1,2],coef.rd[1,1]+(t1+SE[i]-RE[jj.1.tie])*coef.rd[1,2]+t1*(SE[i]-RE[jj.1.tie])*coef.rd[2,2]),length(jj.1.tie),2,byrow=F)
   SS.a<-sapply(coef.rd[1,1]+2*t1*coef.rd[1,2]+t1^2*coef.rd[2,2]+sigma.e^2,solve)
 
   
   slope.a<-COVV.a*as.numeric(SS.a)
   Sigma1.2.a=outer(t(slope.a),(COVV.a),'*'); 
   kappa1=diag(as.matrix(Sigma1.2.a[1,,,1])); kappa2=diag(as.matrix(Sigma1.2.a[1,,,2])); 
   kappa3=diag(as.matrix(Sigma1.2.a[2,,,1])); kappa4=diag(as.matrix(Sigma1.2.a[2,,,2]))
   
     mu0.EXE<-cbind(1,0,matrix(as.numeric(WE[jj.1.tie,]),length(jj.1.tie),dw))%*%coef.fix
      mu.EXE<-cbind(1,SE[i]-RE[jj.1.tie],matrix(as.numeric(WE[jj.1.tie,]),length(jj.1.tie),dw))%*%coef.fix
  
       comp1<-parm1*(mu0.EXE+slope.a[,1]*(Xt1-mu0.EXE))+
       parm2*(mu.EXE+slope.a[,2]*(Xt1-mu0.EXE))
       
       comp2<-parm1^2*(coef.rd[1,1]-kappa1)+
       parm2*parm1*(2*coef.rd[1,1]+2*(SE[i]-RE[jj.1.tie])*coef.rd[1,2]-kappa2-kappa3)+
       parm2^2*(coef.rd[1,1]+(SE[i]-RE[jj.1.tie])^2*coef.rd[2,2]+2*(SE[i]-RE[jj.1.tie])*coef.rd[1,2]-kappa4)
    
    ss=ss-(ind.tie-1)/di*sum(weights[jj.1.tie]*exp(as.numeric(comp1)+0.5*comp2-comp12i+parm3*(SE[i]-RE[jj.1.tie]-ZE[i]*TE[i])
    +parm0*(ZE[jj.1.tie]-ZE[i])-t(as.numeric(WEc[i,])-matrix(as.numeric(t(WEc[jj.1.tie,])),dc,length(jj.1.tie)))%*%parm4 ) )  
   
   rm(Xt1)
   }



   ###### Measure of X(0)* and X(t2)* 
     if (length(jj.2)>0){
     Xt1=datwide.delta0[match(ID[jj.2],datwide.delta0$id),Vname[1]]
     Xt2=datwide.delta0[match(ID[jj.2],datwide.delta0$id),Vname[2]]
     
       ## make t1, t2 vectors
      t1=datwide.delta0[match(ID[jj.2],datwide.delta0$id),Tname[1]]
      t2=datwide.delta0[match(ID[jj.2],datwide.delta0$id),Tname[2]]
         
     COVV.b1<-coef.rd[1,1]+t1*coef.rd[1,2]
     COVV.b2<-coef.rd[1,1]+t2*coef.rd[1,2]
     COVV.b3<-coef.rd[1,1]+(t1+SE[i]-RE[jj.2])*coef.rd[1,2]+t1*(SE[i]-RE[jj.2])*coef.rd[2,2]
     COVV.b4<-coef.rd[1,1]+t2*(SE[i]-RE[jj.2])*coef.rd[2,2]+(t2+SE[i]-RE[jj.2])*coef.rd[1,2]
     

     SS.b.1<-coef.rd[1,1]+2*t1*coef.rd[1,2]+t1*2*coef.rd[2,2]+sigma.e^2
     SS.b.2<-coef.rd[1,1]+2*t2*coef.rd[1,2]+t2^2*coef.rd[2,2]+sigma.e^2
     SS.b.12<-coef.rd[1,1]+(t1+t2)*coef.rd[1,2]+t1*t2*coef.rd[2,2]
     
     SS.b<-array(NA,dim=c(2,2,length(jj.2)))
     SS.b[1,1,]<-SS.b.1
     SS.b[1,2,]<-SS.b[2,1,]<-SS.b.12
     SS.b[2,2,]<-SS.b.2
     
     SS.b<-apply(SS.b,3,solve)  ### SS.b is a 4 by length(jj.2) matrix with each column the inverse of a matrix
 
     ### slope
    
     ### When t1 t2 are vectors, SS.b[1,], SS.b[2,], SS.b[3,], SS.b[4,] 
     ### corresponds [1,1],[2,1],[1,2],[2,2] of a 2x2 matrix
     
     gamma1=COVV.b1*SS.b[1,]+COVV.b2*SS.b[2,]
     gamma2=COVV.b1*SS.b[3,]+COVV.b2*SS.b[4,]
     gamma3=COVV.b3*SS.b[1,]+COVV.b4*SS.b[2,]
     gamma4=COVV.b3*SS.b[3,]+COVV.b4*SS.b[4,]
     

     ### Sigma1.2.b
     eta1=gamma1*COVV.b1+gamma2*COVV.b2
     eta2=gamma1*COVV.b3+gamma2*COVV.b4
     eta3=gamma3*COVV.b1+gamma4*COVV.b2
     eta4=gamma3*COVV.b3+gamma4*COVV.b4
     
      mu0.EXE<-cbind(1,0,matrix(as.numeric(WE[jj.2,]),length(jj.2),dw))%*%coef.fix
       mu.EXE<-cbind(1,SE[i]-RE[jj.2],matrix(as.numeric(WE[jj.2,]),length(jj.2),dw))%*%coef.fix
       mut1.EXE<-cbind(1,t1,matrix(as.numeric(WE[jj.2,]),length(jj.2),dw))%*%coef.fix
       mut2.EXE<-cbind(1,t2,matrix(as.numeric(WE[jj.2,]),length(jj.2),dw))%*%coef.fix
       
      
      comp1<-parm1*(mu0.EXE+gamma1*Xt1+gamma2*Xt2-gamma1*mut1.EXE-gamma2*mut2.EXE)+
        parm2*(mu.EXE+gamma3*Xt1+gamma4*Xt2-gamma3*mut1.EXE-gamma4*mut2.EXE)
      
      comp2<-parm1^2*(coef.rd[1,1]-eta1)+parm2*parm1*(2*coef.rd[1,1]+2*(SE[i]-RE[jj.2])*coef.rd[1,2]-eta2-eta3)+
      parm2^2*(coef.rd[1,1]+(SE[i]-RE[jj.2])^2*coef.rd[2,2]+2*(SE[i]-RE[jj.2])*coef.rd[1,2]-eta4)
      
    
    ss=ss+sum(weights[jj.2]*exp(as.numeric(comp1)+0.5*comp2-comp12i+parm3*(SE[i]-RE[jj.2]-ZE[i]*TE[i])
    +parm0*(ZE[jj.2]-ZE[i])-t(as.numeric(WEc[i,])-matrix(as.numeric(t(WEc[jj.2,])),dc,length(jj.2)))%*%parm4) )  
    rm(Xt1,Xt2)
    }

    ## minus ties

     if (length(jj.2.tie)>0){
     Xt1=datwide.delta0[match(ID[jj.2.tie],datwide.delta0$id),Vname[1]]
     Xt2=datwide.delta0[match(ID[jj.2.tie],datwide.delta0$id),Vname[2]]
     
         ## make t1, t2 vectors
     t1=datwide.delta0[match(ID[jj.2.tie],datwide.delta0$id),Tname[1]]
     t2=datwide.delta0[match(ID[jj.2.tie],datwide.delta0$id),Tname[2]]
    

   ### when t1 t2 are vectors 
     COVV.b1<-coef.rd[1,1]+t1*coef.rd[1,2]
     COVV.b2<-coef.rd[1,1]+t2*coef.rd[1,2]
     COVV.b3<-coef.rd[1,1]+(t1+SE[i]-RE[jj.2.tie])*coef.rd[1,2]+t1*(SE[i]-RE[jj.2.tie])*coef.rd[2,2]
     COVV.b4<-coef.rd[1,1]+t2*(SE[i]-RE[jj.2.tie])*coef.rd[2,2]+(t2+SE[i]-RE[jj.2.tie])*coef.rd[1,2]
     

     SS.b.1<-coef.rd[1,1]+2*t1*coef.rd[1,2]+t1*2*coef.rd[2,2]+sigma.e^2
     SS.b.2<-coef.rd[1,1]+2*t2*coef.rd[1,2]+t2^2*coef.rd[2,2]+sigma.e^2
     SS.b.12<-coef.rd[1,1]+(t1+t2)*coef.rd[1,2]+t1*t2*coef.rd[2,2]
     
     SS.b<-array(NA,dim=c(2,2,length(jj.2.tie)))
     SS.b[1,1,]<-SS.b.1
     SS.b[1,2,]<-SS.b[2,1,]<-SS.b.12
     SS.b[2,2,]<-SS.b.2
     
     SS.b<-apply(SS.b,3,solve)  ### SS.b is a 4 by length(jj.2) matrix with each column the inverse of a matrix
 
     ### slope
    
     ### When t1 t2 are vectors , SS.b[1,], SS.b[2,], SS.b[3,], SS.b[4,] 
     ### corresponds [1,1],[2,1],[1,2],[2,2] of a 2x2 matrix
     
     gamma1=COVV.b1*SS.b[1,]+COVV.b2*SS.b[2,]
     gamma2=COVV.b1*SS.b[3,]+COVV.b2*SS.b[4,]
     gamma3=COVV.b3*SS.b[1,]+COVV.b4*SS.b[2,]
     gamma4=COVV.b3*SS.b[3,]+COVV.b4*SS.b[4,]
     
     ### Sigma1.2.b
     eta1=gamma1*COVV.b1+gamma2*COVV.b2
     eta2=gamma1*COVV.b3+gamma2*COVV.b4
     eta3=gamma3*COVV.b1+gamma4*COVV.b2
     eta4=gamma3*COVV.b3+gamma4*COVV.b4
     
      mu0.EXE<-cbind(1,0,matrix(as.numeric(WE[jj.2.tie,]),length(jj.2.tie),dw))%*%coef.fix
       mu.EXE<-cbind(1,SE[i]-RE[jj.2.tie],matrix(as.numeric(WE[jj.2.tie,]),length(jj.2.tie),dw))%*%coef.fix
       mut1.EXE<-cbind(1,t1,matrix(as.numeric(WE[jj.2.tie,]),length(jj.2.tie),dw))%*%coef.fix
       mut2.EXE<-cbind(1,t2,matrix(as.numeric(WE[jj.2.tie,]),length(jj.2.tie),dw))%*%coef.fix
       
      
      comp1<-parm1*(mu0.EXE+gamma1*Xt1+gamma2*Xt2-gamma1*mut1.EXE-gamma2*mut2.EXE)+
        parm2*(mu.EXE+gamma3*Xt1+gamma4*Xt2-gamma3*mut1.EXE-gamma4*mut2.EXE)
      
      comp2<-parm1^2*(coef.rd[1,1]-eta1)+parm2*parm1*(2*coef.rd[1,1]+2*(SE[i]-RE[jj.2.tie])*coef.rd[1,2]-eta2-eta3)+
      parm2^2*(coef.rd[1,1]+(SE[i]-RE[jj.2.tie])^2*coef.rd[2,2]+2*(SE[i]-RE[jj.2.tie])*coef.rd[1,2]-eta4)
      
    
    ss=ss-(ind.tie-1)/di*sum(weights[jj.2.tie]*exp(as.numeric(comp1)+0.5*comp2-comp12i+parm3*(SE[i]-RE[jj.2.tie]-ZE[i]*TE[i])
    +parm0*(ZE[jj.2.tie]-ZE[i])-t(as.numeric(WEc[i,])-matrix(as.numeric(t(WEc[jj.2.tie,])),dc,length(jj.2.tie)))%*%parm4) )  
    rm(Xt1,Xt2)
    }
    #####  ###### Measure of X(t1)*, X(t2)*, and X(t3)*
    
     if (length(jj.3)>0){
     Xt1=datwide.delta0[match(ID[jj.3],datwide.delta0$id),Vname[1]]
     Xt2=datwide.delta0[match(ID[jj.3],datwide.delta0$id),Vname[2]]
     Xt3=datwide.delta0[match(ID[jj.3],datwide.delta0$id),Vname[3]]
     
       ### t1, t2, t3 vector
     t1=datwide.delta0[match(ID[jj.3],datwide.delta0$id),Tname[1]]
     t2=datwide.delta0[match(ID[jj.3],datwide.delta0$id),Tname[2]]
     t3=datwide.delta0[match(ID[jj.3],datwide.delta0$id),Tname[3]]
     
   ####  t1,t2,t3 vector

     COVV.c1=coef.rd[1,1]+t1*coef.rd[1,2]
     COVV.c2=coef.rd[1,1]+t2*coef.rd[1,2]
     COVV.c3=coef.rd[1,1]+t3*coef.rd[1,2]
     COVV.c4=coef.rd[1,1]+(t1+SE[i]-RE[jj.3])*coef.rd[1,2]+t1*(SE[i]-RE[jj.3])*coef.rd[2,2]
     COVV.c5=coef.rd[1,1]+(t2+SE[i]-RE[jj.3])*coef.rd[1,2]+t2*(SE[i]-RE[jj.3])*coef.rd[2,2]
     COVV.c6=coef.rd[1,1]+(t3+SE[i]-RE[jj.3])*coef.rd[1,2]+t3*(SE[i]-RE[jj.3])*coef.rd[2,2]
     
     ##########
    
     SS.c.11<-coef.rd[1,1]+2*t1*coef.rd[1,2]+t1^2*coef.rd[2,2]+sigma.e^2
     SS.c.22<-coef.rd[1,1]+2*t2*coef.rd[1,2]+t2^2*coef.rd[2,2]+sigma.e^2
     SS.c.33<-coef.rd[1,1]+2*t3*coef.rd[1,2]+t3^2*coef.rd[2,2]+sigma.e^2
     SS.c.12=coef.rd[1,1]+(t1+t2)*coef.rd[1,2]+t1*t2*coef.rd[2,2]
     SS.c.13=coef.rd[1,1]+(t1+t3)*coef.rd[1,2]+t1*t3*coef.rd[2,2]
     SS.c.23=coef.rd[1,1]+(t2+t3)*coef.rd[1,2]+t2*t3*coef.rd[2,2]
     
     SS.c<-array(NA,dim=c(3,3,length(jj.3)))
     SS.c[1,1,]<-SS.c.11
     SS.c[2,2,]<-SS.c.22
     SS.c[3,3,]<-SS.c.33
     SS.c[1,2,]<-SS.c[2,1,]<-SS.c.12
     SS.c[1,3,]<-SS.c[3,1,]<-SS.c.13
     SS.c[2,3,]<-SS.c[3,2,]<-SS.c.23
     
    SS.c<-apply(SS.c,3,solve)
    

     #### SS.c is 9xlength(jj.3), elements in each column corresponds to [1,1],[2,1],[3,1],[1,2],[2,2],[3,2],[1,3],[2,3],[3,3] of a matrix

     delta1=COVV.c1*SS.c[1,]+COVV.c2*SS.c[2,]+COVV.c3*SS.c[3,]
     delta2=COVV.c1*SS.c[4,]+COVV.c2*SS.c[5,]+COVV.c3*SS.c[6,]
     delta3=COVV.c1*SS.c[7,]+COVV.c2*SS.c[8,]+COVV.c3*SS.c[9,]
     delta4=COVV.c4*SS.c[1,]+COVV.c5*SS.c[2,]+COVV.c6*SS.c[3,]
     delta5=COVV.c4*SS.c[4,]+COVV.c5*SS.c[5,]+COVV.c6*SS.c[6,]
     delta6=COVV.c4*SS.c[7,]+COVV.c5*SS.c[8,]+COVV.c6*SS.c[9,]

     xi1=delta1*COVV.c1+delta2*COVV.c2+delta3*COVV.c3
     xi2=delta1*COVV.c4+delta2*COVV.c5+delta3*COVV.c6
     xi3=delta4*COVV.c1+delta5*COVV.c2+delta6*COVV.c3
     xi4=delta4*COVV.c4+delta5*COVV.c5+delta6*COVV.c6
      
      mu0.EXE<-cbind(1,0,matrix(as.numeric(WE[jj.3,]),length(jj.3),dw))%*%coef.fix
       mu.EXE<-cbind(1,SE[i]-RE[jj.3],matrix(as.numeric(WE[jj.3,]),length(jj.3),dw))%*%coef.fix
       mut1.EXE<-cbind(1,t1,matrix(as.numeric(WE[jj.3,]),length(jj.3),dw))%*%coef.fix
       mut2.EXE<-cbind(1,t2,matrix(as.numeric(WE[jj.3,]),length(jj.3),dw))%*%coef.fix
       mut3.EXE<-cbind(1,t3,matrix(as.numeric(WE[jj.3,]),length(jj.3),dw))%*%coef.fix
       
      
      comp1<-parm1*(mu0.EXE+delta1*Xt1+delta2*Xt2+delta3*Xt3-delta1*mut1.EXE-delta2*mut2.EXE-delta3*mut3.EXE)+
        parm2*(mu.EXE+delta4*Xt1+delta5*Xt2+delta6*Xt3-delta4*mut1.EXE-delta5*mut2.EXE-delta6*mut3.EXE)
      
      comp2<-parm1^2*(coef.rd[1,1]-xi1)+parm2*parm1*(2*coef.rd[1,1]+2*(SE[i]-RE[jj.3])*coef.rd[1,2]-xi2-xi3)+
      parm2^2*(coef.rd[1,1]+(SE[i]-RE[jj.3])^2*coef.rd[2,2]+2*(SE[i]-RE[jj.3])*coef.rd[1,2]-xi4)
      
    
    ss=ss+sum(weights[jj.3]*exp(as.numeric(comp1)+0.5*comp2-comp12i+parm3*(SE[i]-RE[jj.3]-ZE[i]*TE[i])
    +parm0*(ZE[jj.3]-ZE[i])-t(as.numeric(WEc[i,])-matrix(as.numeric(t(WEc[jj.3,])),dc,length(jj.3)))%*%parm4) )  
   
    rm(Xt1,Xt2,Xt3)
    }

    #### Minus ties
 
     if (length(jj.3.tie)>0){
     Xt1=datwide.delta0[match(ID[jj.3.tie],datwide.delta0$id),Vname[1]]
     Xt2=datwide.delta0[match(ID[jj.3.tie],datwide.delta0$id),Vname[2]]
     Xt3=datwide.delta0[match(ID[jj.3.tie],datwide.delta0$id),Vname[3]]
     
      ### t1, t2, t3 vector
     t1=datwide.delta0[match(ID[jj.3.tie],datwide.delta0$id),Tname[1]]
     t2=datwide.delta0[match(ID[jj.3.tie],datwide.delta0$id),Tname[2]]
     t3=datwide.delta0[match(ID[jj.3.tie],datwide.delta0$id),Tname[3]]
 
     COVV.c1=coef.rd[1,1]+t1*coef.rd[1,2]
     COVV.c2=coef.rd[1,1]+t2*coef.rd[1,2]
     COVV.c3=coef.rd[1,1]+t3*coef.rd[1,2]
     COVV.c4=coef.rd[1,1]+(t1+SE[i]-RE[jj.3.tie])*coef.rd[1,2]+t1*(SE[i]-RE[jj.3.tie])*coef.rd[2,2]
     COVV.c5=coef.rd[1,1]+(t2+SE[i]-RE[jj.3.tie])*coef.rd[1,2]+t2*(SE[i]-RE[jj.3.tie])*coef.rd[2,2]
     COVV.c6=coef.rd[1,1]+(t3+SE[i]-RE[jj.3.tie])*coef.rd[1,2]+t3*(SE[i]-RE[jj.3.tie])*coef.rd[2,2]
     
     ##########
    
     SS.c.11<-coef.rd[1,1]+2*t1*coef.rd[1,2]+t1^2*coef.rd[2,2]+sigma.e^2
     SS.c.22<-coef.rd[1,1]+2*t2*coef.rd[1,2]+t2^2*coef.rd[2,2]+sigma.e^2
     SS.c.33<-coef.rd[1,1]+2*t3*coef.rd[1,2]+t3^2*coef.rd[2,2]+sigma.e^2
     SS.c.12=coef.rd[1,1]+(t1+t2)*coef.rd[1,2]+t1*t2*coef.rd[2,2]
     SS.c.13=coef.rd[1,1]+(t1+t3)*coef.rd[1,2]+t1*t3*coef.rd[2,2]
     SS.c.23=coef.rd[1,1]+(t2+t3)*coef.rd[1,2]+t2*t3*coef.rd[2,2]
     
     SS.c<-array(NA,dim=c(3,3,length(jj.3.tie)))
     SS.c[1,1,]<-SS.c.11
     SS.c[2,2,]<-SS.c.22
     SS.c[3,3,]<-SS.c.33
     SS.c[1,2,]<-SS.c[2,1,]<-SS.c.12
     SS.c[1,3,]<-SS.c[3,1,]<-SS.c.13
     SS.c[2,3,]<-SS.c[3,2,]<-SS.c.23
         
    SS.c<-apply(SS.c,3,solve)
    

     #### SS.c is 9xlength(jj.3), elements in each column corresponds to [1,1],[2,1],[3,1],[1,2],[2,2],[3,2],[1,3],[2,3],[3,3] of a matrix

     delta1=COVV.c1*SS.c[1,]+COVV.c2*SS.c[2,]+COVV.c3*SS.c[3,]
     delta2=COVV.c1*SS.c[4,]+COVV.c2*SS.c[5,]+COVV.c3*SS.c[6,]
     delta3=COVV.c1*SS.c[7,]+COVV.c2*SS.c[8,]+COVV.c3*SS.c[9,]
     delta4=COVV.c4*SS.c[1,]+COVV.c5*SS.c[2,]+COVV.c6*SS.c[3,]
     delta5=COVV.c4*SS.c[4,]+COVV.c5*SS.c[5,]+COVV.c6*SS.c[6,]
     delta6=COVV.c4*SS.c[7,]+COVV.c5*SS.c[8,]+COVV.c6*SS.c[9,]

     
     xi1=delta1*COVV.c1+delta2*COVV.c2+delta3*COVV.c3
     xi2=delta1*COVV.c4+delta2*COVV.c5+delta3*COVV.c6
     xi3=delta4*COVV.c1+delta5*COVV.c2+delta6*COVV.c3
     xi4=delta4*COVV.c4+delta5*COVV.c5+delta6*COVV.c6
      
      mu0.EXE<-cbind(1,0,matrix(as.numeric(WE[jj.3.tie,]),length(jj.3.tie),dw))%*%coef.fix
       mu.EXE<-cbind(1,SE[i]-RE[jj.3.tie],matrix(as.numeric(WE[jj.3.tie,]),length(jj.3.tie),dw))%*%coef.fix
       mut1.EXE<-cbind(1,t1,matrix(as.numeric(WE[jj.3.tie,]),length(jj.3.tie),dw))%*%coef.fix
       mut2.EXE<-cbind(1,t2,matrix(as.numeric(WE[jj.3.tie,]),length(jj.3.tie),dw))%*%coef.fix
       mut3.EXE<-cbind(1,t3,matrix(as.numeric(WE[jj.3.tie,]),length(jj.3.tie),dw))%*%coef.fix
            
      comp1<-parm1*(mu0.EXE+delta1*Xt1+delta2*Xt2+delta3*Xt3-delta1*mut1.EXE-delta2*mut2.EXE-delta3*mut3.EXE)+
        parm2*(mu.EXE+delta4*Xt1+delta5*Xt2+delta6*Xt3-delta4*mut1.EXE-delta5*mut2.EXE-delta6*mut3.EXE)
      
      comp2<-parm1^2*(coef.rd[1,1]-xi1)+parm2*parm1*(2*coef.rd[1,1]+2*(SE[i]-RE[jj.3.tie])*coef.rd[1,2]-xi2-xi3)+
      parm2^2*(coef.rd[1,1]+(SE[i]-RE[jj.3.tie])^2*coef.rd[2,2]+2*(SE[i]-RE[jj.3.tie])*coef.rd[1,2]-xi4)
      
    ss=ss-(ind.tie-1)/di*sum(weights[jj.3.tie]*exp(as.numeric(comp1)+0.5*comp2-comp12i+parm3*(SE[i]-RE[jj.3.tie]-ZE[i]*TE[i])
    +parm0*(ZE[jj.3.tie]-ZE[i])-t(as.numeric(WEc[i,])-matrix(as.numeric(t(WEc[jj.3.tie,])),dc,length(jj.3.tie)))%*%parm4) )  
   
    rm(Xt1,Xt2,Xt3)
    }

    ss.new=ss
    ### placebo
   
   ### need to use jj.Z0 and jj.Z0.tie when there are confounders in the risk model
   if (length(jj.Z0)>0)
   ss=ss+sum(weights[jj.Z0]*exp(-comp12i-parm3*(ZE[i]*TE[i])-parm0*(ZE[i])
   -t(as.numeric(WEc[i,])-matrix(as.numeric(t(WEc[jj.Z0,])),dc,length(jj.Z0)))%*%parm4))
   
   if (length(jj.Z0.tie)>0)
   ss=ss-(ind.tie-1)/di*sum(weights[jj.Z0.tie]*exp(-comp12i-parm3*(ZE[i]*TE[i])-parm0*(ZE[i])-t(as.numeric(WEc[i,])-matrix(as.numeric(t(WEc[jj.Z0.tie,])),dc,length(jj.Z0.tie)))%*%parm4))
   
    out=out+weights.bar*log(weights[i]+ss)+(weights.bar-weights[i])*(comp12i+parm3*(SE[i]-RE[i])*ZE[i]+parm0*(ZE[i])
    +matrix(WEc[i,],length(i),dc)%*%parm4)
    
    outi=c(outi,ss)
 
    }
 return(out)
 }

########

######## Fit a Cox hazard model on calendar time with coefficient for Z, Z*x(0), Z*t, and W

nloglik.M1.b.Cal.13<- function(parm){

 out=0
 
 outi<-NULL
   ### coef for Z, Z*x(0), Z*t, and confounders respectively
    parm0=parm[1]
    parm1=parm[2]
    #parm2=parm[3]
    parm3=parm[3]
    parm4=parm[4]
   
   jj.4.full<-(1:length(RE))[datwide.delta0[match(ID,datwide.delta0$id),]$groupall > datwide.delta0[match(ID,datwide.delta0$id),]$group]
   jj.4.full<-as.numeric(na.omit(jj.4.full))
   out.jj.4.full<-NULL
   for (j in jj.4.full){
      
    keep=1:ll
 
    VVj.11<-coef.rd[1,1]
    VVj.12<-VVj.21<-coef.rd[1,1]+(SE[keep]-RE[j])*coef.rd[1,2]
    VVj.22<-coef.rd[1,1]+2*(SE[keep]-RE[j])*coef.rd[1,2]+(SE[keep]-RE[j])^2*coef.rd[2,2]
    
    muj.1=cbind(1,0,matrix(WE[j,],length(j),dw))%*%coef.fix
    muj.2=cbind(1,SE[keep]-RE[j],matrix(WE[j,],length(j),dw))%*%coef.fix
    
   
    oo.j<-id.long==ID[j]
    Xj<-X.long[oo.j]
    Tj<-T.long[oo.j]
    
    
    VVj.c<-cbind(1,Tj)%*%coef.rd%*%rbind(1,Tj)+diag(rep(sigma.e^2,length(Tj)),nrow=length(Tj),ncol=length(Tj))               
   
    COVVj.1<-coef.rd[1,1]+coef.rd[1,2]*Tj
    COVVj.2<-coef.rd[1,1]+coef.rd[1,2]*outer(SE[keep]-RE[j],Tj,'+')+coef.rd[2,2]*outer(SE[keep]-RE[j],Tj,'*')
    COVVj.1=matrix(COVVj.1,1,length(COVVj.1))
    if (length(Tj)==1) COVVj.2=matrix(COVVj.2,length(COVVj.2),1)
    
    muj.c<-cbind(1,Tj,matrix(WE[j,],length(j),dw))%*%coef.fix
        
    comp1j<-cbind(as.numeric(muj.1),muj.2)+cbind(as.numeric(COVVj.1%*%solve(VVj.c)%*%(Xj-muj.c)),COVVj.2%*%solve(VVj.c)%*%(Xj-muj.c))
    
     VVjminus.11<-COVVj.1%*%solve(VVj.c)%*%(t(COVVj.1))
     VVjminus.12<-COVVj.2%*%solve(VVj.c)%*%(t(COVVj.1))
     VVjminus.22<-diag(COVVj.2%*%solve(VVj.c)%*%(t(COVVj.2)))
     
     comp2j.11<-VVj.11-VVjminus.11
     comp2j.12<-VVj.12-VVjminus.12
     comp2j.22<-VVj.22-VVjminus.22
   
     out.jj.4.full<-cbind(out.jj.4.full,as.numeric(c(parm1,parm2)%*%t(comp1j))+parm1^2*as.numeric(comp2j.11)/2+parm1*parm2*comp2j.12+parm2^2*comp2j.22/2+parm3*(SE[keep]-RE[j])
     +parm0*ZE[j]+matrix(WEc[j,],length(j),dc)%*%parm4)
     
    }
        
    S.tie=SE.u[1]  ### unique SE value
    ind.tie=0       ### position within the unique value
    
    for (i in 1:ll){
    
    ## track #of ties
    
        di<-sum(SE==SE[i] & delta==1)
        
        if (SE[i]==S.tie){  ## in the same tie as previous value
          ind.tie=ind.tie+1
        } else { ## go to a new SE value
           S.tie=SE[i]
           ind.tie=1
        }
     ### compute average weight among ties
        weights.bar=mean(weights[SE==SE[i] & delta==1])
     
    ss=0
    comp12i=0 ### placebo contribute nothing
    if (ZE[i]==1) { 
    VVi<-cbind(1,c(0,TE[i]))%*%coef.rd%*%rbind(1,c(0,TE[i]))
    mui<-rbind(c(1,0,matrix(WE[i,],length(i),dw)),c(1,TE[i],matrix(WE[i,],length(i),dw)))%*%coef.fix
   
    oo.i<-id.long==ID[i]
        
       #### (i) if case has at least one X(Ti)* measured, i.e. case ID[i] belong to id.long
        if (sum(oo.i)>0){ 
        Xi<-X.long[oo.i]
        Ti<-T.long[oo.i]
        
        VVi.c<-cbind(1,Ti)%*%coef.rd%*%rbind(1,Ti)+diag(rep(sigma.e^2,length(Ti)),nrow=length(Ti),ncol=length(Ti))               
        COVVi<-cbind(1,c(0,TE[i]))%*%coef.rd%*%rbind(1,Ti)
        mui.c<-cbind(1,Ti,matrix(WE[i,],length(Ti),dw,byrow=TRUE))%*%coef.fix
        
        comp1i<-mui+COVVi%*%solve(VVi.c)%*%(Xi-mui.c)      ## E(X(0),X(t)|X(Ti)*,W)
        comp2i<-VVi-COVVi%*%solve(VVi.c)%*%t(COVVi)        ## Var(X(0),X(t)|X(Ti)*,W)
        
        } else {
        #### (ii) if case does not have X(Ti)* measure, i.e. case ID[i] not in id.long
        
       comp1i<-mui            
       comp2i<-VVi            
       } 
 
     comp12i<-as.numeric(t(c(parm1,parm2))%*%comp1i+t(c(parm1,parm2))%*%comp2i%*%c(parm1,parm2)/2)
    
    } 
    
    jj.r<-((i+1):(length(RE)))[RE[(i+1):length(RE)]<SE[i] & SE[(i+1):length(RE)]>=SE[i] & ID[(i+1):length(RE)]%in%id.long]
    jj.l<-NULL
    if (i>1) jj.l<-(1:(i-1))[RE[1:(i-1)]<SE[i] & SE[1:(i-1)]>=SE[i] & ID[1:(i-1)]%in%id.long]
    jj<-c(jj.l,jj.r) ## among risk set  
    
    jj.tie<-(1:(length(RE)))[RE[1:length(RE)]<SE[i] & ID[1:length(RE)]%in%id.long & SE[1:length(SE)]==SE[i] & delta[1:length(delta)]==1]
   
    jj.W.r<-((i+1):(length(RE)))[RE[(i+1):length(RE)]<SE[i] & ZE[(i+1):length(ZE)]==1 & SE[(i+1):length(RE)]>=SE[i] & !ID[(i+1):length(RE)]%in%id.long]
    jj.W.l<-NULL
    if (i>1) jj.W.l<-(1:(i-1))[RE[1:(i-1)]<SE[i] & ZE[1:(i-1)]==1 & !ID[1:(i-1)]%in%id.long & SE[1:(i-1)]>=SE[i]]
    jj.W<-c(jj.W.l,jj.W.r)    

    
    jj.W.tie<-(1:(length(RE)))[RE[1:length(RE)]<SE[i] & ZE[1:length(ZE)]==1 & SE[1:length(RE)]==SE[i] & !ID[1:length(RE)]%in%id.long  & delta[1:length(delta)]==1]
 
  
    jj.Z0.r<-((i+1):length(RE))[RE[(i+1):length(RE)]<SE[i] & ZE[(i+1):length(ZE)]==0 & SE[(i+1):length(SE)]>=SE[i]]  
    jj.Z0.l<-NULL
    if (i>1) jj.Z0.l<-(1:(i-1))[RE[1:(i-1)]<SE[i] & ZE[1:(i-1)]==0 & SE[1:(i-1)]>=SE[i]]
    jj.Z0=c(jj.Z0.l,jj.Z0.r)

    jj.Z0.tie<-(1:length(RE))[RE[1:length(RE)]<SE[i] & ZE[1:length(ZE)]==0 & SE[1:length(SE)]==SE[i] & delta[1:length(delta)]==1]  
  
    jj.1.r<-((i+1):(length(RE)))[RE[(i+1):length(RE)]<SE[i] & datwide.delta0[match(ID[(i+1):length(RE)],datwide.delta0$id),]$group==1 &
       datwide.delta0[match(ID[(i+1):length(RE)],datwide.delta0$id),]$groupall==1 & SE[(i+1):length(RE)]>=SE[i]]
    jj.1.l<-NULL
    if (i>1) jj.1.l<-(1:(i-1))[RE[1:(i-1)]<SE[i] & datwide.delta0[match(ID[1:(i-1)],datwide.delta0$id),]$group==1 &
       datwide.delta0[match(ID[1:(i-1)],datwide.delta0$id),]$groupall==1 & SE[1:(i-1)]>=SE[i]]
    jj.1<-c(jj.1.l,jj.1.r) ## among risk set

    
    jj.1.tie<-(1:(length(RE)))[RE[1:length(RE)]<SE[i] & datwide.delta0[match(ID[1:length(RE)],datwide.delta0$id),]$group==1 &
       datwide.delta0[match(ID[1:length(RE)],datwide.delta0$id),]$groupall==1 & SE[1:length(SE)]==SE[i] & delta[1:length(delta)]==1]
   
   
    #######
 
    jj.2.r<-((i+1):(length(RE)))[RE[(i+1):length(RE)]<SE[i] & datwide.delta0[match(ID[(i+1):length(RE)],datwide.delta0$id),]$group==2 & 
       datwide.delta0[match(ID[(i+1):length(RE)],datwide.delta0$id),]$groupall==2 & SE[(i+1):length(SE)]>=SE[i]]
    jj.2.l<-NULL
    if (i>1)   jj.2.l<-(1:(i-1))[RE[1:(i-1)]<SE[i] & datwide.delta0[match(ID[1:(i-1)],datwide.delta0$id),]$group==2 & 
       datwide.delta0[match(ID[1:(i-1)],datwide.delta0$id),]$groupall==2 & SE[1:(i-1)]>=SE[i]]
    jj.2<-c(jj.2.l,jj.2.r)

    jj.2.tie<-(1:(length(RE)))[RE[1:length(RE)]<SE[i] & datwide.delta0[match(ID[1:length(RE)],datwide.delta0$id),]$group==2 & 
       datwide.delta0[match(ID[1:length(RE)],datwide.delta0$id),]$groupall==2 & SE[1:length(SE)]==SE[i] & delta[1:length(delta)]==1]
   

    ###
    jj.3.r<-((i+1):(length(RE)))[RE[(i+1):length(RE)]<SE[i] & datwide.delta0[match(ID[(i+1):length(RE)],datwide.delta0$id),]$group==3 & 
       datwide.delta0[match(ID[(i+1):length(RE)],datwide.delta0$id),]$groupall==3 & SE[(i+1):length(SE)]>=SE[i]]
    jj.3.l<-NULL
    if (i>1) jj.3.l<-(1:(i-1))[RE[1:(i-1)]<SE[i] & datwide.delta0[match(ID[1:(i-1)],datwide.delta0$id),]$group==3 & 
       datwide.delta0[match(ID[1:(i-1)],datwide.delta0$id),]$groupall==3 & SE[1:(i-1)]>=SE[i]]
    jj.3<-c(jj.3.l,jj.3.r)
   

    jj.3.tie<-(1:(length(RE)))[RE[1:length(RE)]<SE[i] & datwide.delta0[match(ID[1:length(RE)],datwide.delta0$id),]$group==3 & 
       datwide.delta0[match(ID[1:length(RE)],datwide.delta0$id),]$groupall==3 & SE[1:length(SE)]==SE[i] & delta[1:length(delta)]==1]
  
     jj.1<-as.numeric(na.omit(jj.1))
     jj.2<-as.numeric(na.omit(jj.2))
     jj.3<-as.numeric(na.omit(jj.3))
   
       jj.1.tie<-as.numeric(na.omit(jj.1.tie))
     jj.2.tie<-as.numeric(na.omit(jj.2.tie))
     jj.3.tie<-as.numeric(na.omit(jj.3.tie))
     
    ss=ss+sum((RE[jj.4.full]<SE[i])*(SE[jj.4.full]>=SE[i])*(!jj.4.full==i)*weights[jj.4.full]*exp(out.jj.4.full[i,]-comp12i-parm3*(ZE[i]*TE[i])
    -parm0*ZE[i]-matrix(WEc[i,],length(i),dc)%*%parm4))  ## risk set from jj.4.full
     
    ss=ss-(ind.tie-1)/di*sum((RE[jj.4.full]<SE[i])*(SE[jj.4.full]==SE[i])*(!jj.4.full==i)*weights[jj.4.full]*exp(out.jj.4.full[i,]-comp12i-parm3*(ZE[i]*TE[i])
    -parm0*ZE[i]-matrix(WEc[i,],length(i),dc)%*%parm4))  ## tie from jj.4.full
     

   if (length(jj.W)>0){
      ### No X measure 
    mu0.EXE<-cbind(1,0,matrix(as.numeric(WE[jj.W,]),length(jj.W),dw))%*%coef.fix
    mu.EXE<-cbind(1,SE[i]-RE[jj.W],matrix(as.numeric(WE[jj.W,]),length(jj.W),dw))%*%coef.fix
  
   comp1.b<-parm1*mu0.EXE+parm2*mu.EXE
       
   comp2.b<-parm1^2*coef.rd[1,1]+2*parm1*parm2*(coef.rd[1,1]+(SE[i]-RE[jj.W])*coef.rd[1,2])+
   parm2^2*(coef.rd[1,1]+2*(SE[i]-RE[jj.W])*coef.rd[1,2]+(SE[i]-RE[jj.W])^2*coef.rd[2,2])
   
   ss=ss+sum(weights[jj.W]*exp(as.numeric(comp1.b)+0.5*comp2.b-comp12i+parm3*(SE[i]-RE[jj.W]-ZE[i]*TE[i])
   +parm0*(ZE[jj.W]-ZE[i])-t(as.numeric(WEc[i,])-matrix(as.numeric(t(WEc[jj.W,])),dc,length(jj.W)))%*%parm4) )  
    
   } 
   ### minus ties
   if (length(jj.W.tie)>0){
   
        mu0.EXE<-cbind(1,0,matrix(as.numeric(WE[jj.W.tie,]),length(jj.W.tie),dw))%*%coef.fix
        mu.EXE<-cbind(1,SE[i]-RE[jj.W.tie],matrix(as.numeric(WE[jj.W.tie,]),length(jj.W.tie),dw))%*%coef.fix
      
       comp1.b<-parm1*mu0.EXE+parm2*mu.EXE
           
       comp2.b<-parm1^2*coef.rd[1,1]+2*parm1*parm2*(coef.rd[1,1]+(SE[i]-RE[jj.W.tie])*coef.rd[1,2])+
       parm2^2*(coef.rd[1,1]+2*(SE[i]-RE[jj.W.tie])*coef.rd[1,2]+(SE[i]-RE[jj.W.tie])^2*coef.rd[2,2])
       
       ss=ss-(ind.tie-1)/di*sum(weights[jj.W.tie]*exp(as.numeric(comp1.b)+0.5*comp2.b-comp12i+parm3*(SE[i]-RE[jj.W.tie]-ZE[i]*TE[i])
       +parm0*(ZE[jj.W.tie]-ZE[i])-t(as.numeric(WEc[i,])-matrix(as.numeric(t(WEc[jj.W.tie,])),dc,length(jj.W.tie)))%*%parm4) ) 
   }      

 
   ### Measure of X(0)* only
   
    #account for different t levels
   
   if (length(jj.1)>0){
   Xt1=datwide.delta0[match(ID[jj.1],datwide.delta0$id),Vname[1]]
    t1=datwide.delta0[match(ID[jj.1],datwide.delta0$id),Tname[1]]
 
    ## make t1 vector
   COVV.a<-matrix(c(coef.rd[1,1]+t1*coef.rd[1,2],coef.rd[1,1]+(t1+SE[i]-RE[jj.1])*coef.rd[1,2]+t1*(SE[i]-RE[jj.1])*coef.rd[2,2]),length(jj.1),2,byrow=F)
   SS.a<-sapply(coef.rd[1,1]+2*t1*coef.rd[1,2]+t1^2*coef.rd[2,2]+sigma.e^2,solve)
   
   slope.a<-COVV.a*as.numeric(SS.a)
   Sigma1.2.a=outer(t(slope.a),(COVV.a),'*'); 
   kappa1=diag(as.matrix(Sigma1.2.a[1,,,1])); kappa2=diag(as.matrix(Sigma1.2.a[1,,,2])); 
   kappa3=diag(as.matrix(Sigma1.2.a[2,,,1])); kappa4=diag(as.matrix(Sigma1.2.a[2,,,2]))
   
     mu0.EXE<-cbind(1,0,matrix(as.numeric(WE[jj.1,]),length(jj.1),dw))%*%coef.fix
      mu.EXE<-cbind(1,SE[i]-RE[jj.1],matrix(as.numeric(WE[jj.1,]),length(jj.1),dw))%*%coef.fix
  
       comp1<-parm1*(mu0.EXE+slope.a[,1]*(Xt1-mu0.EXE))+
       parm2*(mu.EXE+slope.a[,2]*(Xt1-mu0.EXE))
       
       comp2<-parm1^2*(coef.rd[1,1]-kappa1)+
       parm2*parm1*(2*coef.rd[1,1]+2*(SE[i]-RE[jj.1])*coef.rd[1,2]-kappa2-kappa3)+
       parm2^2*(coef.rd[1,1]+(SE[i]-RE[jj.1])^2*coef.rd[2,2]+2*(SE[i]-RE[jj.1])*coef.rd[1,2]-kappa4)
    
    ss=ss+sum(weights[jj.1]*exp(as.numeric(comp1)+0.5*comp2-comp12i+parm3*(SE[i]-RE[jj.1]-ZE[i]*TE[i])
    +parm0*(ZE[jj.1]-ZE[i])-t(as.numeric(WEc[i,])-matrix(as.numeric(t(WEc[jj.1,])),dc,length(jj.1)))%*%parm4) )  
   
   rm(Xt1)
   }

   ## minus ties

   if (length(jj.1.tie)>0){
   Xt1=datwide.delta0[match(ID[jj.1.tie],datwide.delta0$id),Vname[1]]
   
   t1=datwide.delta0[match(ID[jj.1.tie],datwide.delta0$id),Tname[1]]
   
    # make t1 vector
   COVV.a<-matrix(c(coef.rd[1,1]+t1*coef.rd[1,2],coef.rd[1,1]+(t1+SE[i]-RE[jj.1.tie])*coef.rd[1,2]+t1*(SE[i]-RE[jj.1.tie])*coef.rd[2,2]),length(jj.1.tie),2,byrow=F)
   SS.a<-sapply(coef.rd[1,1]+2*t1*coef.rd[1,2]+t1^2*coef.rd[2,2]+sigma.e^2,solve)
 
   
   slope.a<-COVV.a*as.numeric(SS.a)
   Sigma1.2.a=outer(t(slope.a),(COVV.a),'*'); 
   kappa1=diag(as.matrix(Sigma1.2.a[1,,,1])); kappa2=diag(as.matrix(Sigma1.2.a[1,,,2])); 
   kappa3=diag(as.matrix(Sigma1.2.a[2,,,1])); kappa4=diag(as.matrix(Sigma1.2.a[2,,,2]))
   
     mu0.EXE<-cbind(1,0,matrix(as.numeric(WE[jj.1.tie,]),length(jj.1.tie),dw))%*%coef.fix
      mu.EXE<-cbind(1,SE[i]-RE[jj.1.tie],matrix(as.numeric(WE[jj.1.tie,]),length(jj.1.tie),dw))%*%coef.fix
  
       comp1<-parm1*(mu0.EXE+slope.a[,1]*(Xt1-mu0.EXE))+
       parm2*(mu.EXE+slope.a[,2]*(Xt1-mu0.EXE))
       
       comp2<-parm1^2*(coef.rd[1,1]-kappa1)+
       parm2*parm1*(2*coef.rd[1,1]+2*(SE[i]-RE[jj.1.tie])*coef.rd[1,2]-kappa2-kappa3)+
       parm2^2*(coef.rd[1,1]+(SE[i]-RE[jj.1.tie])^2*coef.rd[2,2]+2*(SE[i]-RE[jj.1.tie])*coef.rd[1,2]-kappa4)
    
    ss=ss-(ind.tie-1)/di*sum(weights[jj.1.tie]*exp(as.numeric(comp1)+0.5*comp2-comp12i+parm3*(SE[i]-RE[jj.1.tie]-ZE[i]*TE[i])
    +parm0*(ZE[jj.1.tie]-ZE[i])-t(as.numeric(WEc[i,])-matrix(as.numeric(t(WEc[jj.1.tie,])),dc,length(jj.1.tie)))%*%parm4 ) )  
   
   rm(Xt1)
   }



   ###### Measure of X(0)* and X(t2)* 
     if (length(jj.2)>0){
     Xt1=datwide.delta0[match(ID[jj.2],datwide.delta0$id),Vname[1]]
     Xt2=datwide.delta0[match(ID[jj.2],datwide.delta0$id),Vname[2]]
     
       ## make t1, t2 vectors
      t1=datwide.delta0[match(ID[jj.2],datwide.delta0$id),Tname[1]]
      t2=datwide.delta0[match(ID[jj.2],datwide.delta0$id),Tname[2]]
         
     COVV.b1<-coef.rd[1,1]+t1*coef.rd[1,2]
     COVV.b2<-coef.rd[1,1]+t2*coef.rd[1,2]
     COVV.b3<-coef.rd[1,1]+(t1+SE[i]-RE[jj.2])*coef.rd[1,2]+t1*(SE[i]-RE[jj.2])*coef.rd[2,2]
     COVV.b4<-coef.rd[1,1]+t2*(SE[i]-RE[jj.2])*coef.rd[2,2]+(t2+SE[i]-RE[jj.2])*coef.rd[1,2]
     

     SS.b.1<-coef.rd[1,1]+2*t1*coef.rd[1,2]+t1*2*coef.rd[2,2]+sigma.e^2
     SS.b.2<-coef.rd[1,1]+2*t2*coef.rd[1,2]+t2^2*coef.rd[2,2]+sigma.e^2
     SS.b.12<-coef.rd[1,1]+(t1+t2)*coef.rd[1,2]+t1*t2*coef.rd[2,2]
     
     SS.b<-array(NA,dim=c(2,2,length(jj.2)))
     SS.b[1,1,]<-SS.b.1
     SS.b[1,2,]<-SS.b[2,1,]<-SS.b.12
     SS.b[2,2,]<-SS.b.2
     
     SS.b<-apply(SS.b,3,solve)  ### SS.b is a 4 by length(jj.2) matrix with each column the inverse of a matrix
 
     ### slope
    
     ### When t1 t2 are vectors, SS.b[1,], SS.b[2,], SS.b[3,], SS.b[4,] 
     ### corresponds [1,1],[2,1],[1,2],[2,2] of a 2x2 matrix
     
     gamma1=COVV.b1*SS.b[1,]+COVV.b2*SS.b[2,]
     gamma2=COVV.b1*SS.b[3,]+COVV.b2*SS.b[4,]
     gamma3=COVV.b3*SS.b[1,]+COVV.b4*SS.b[2,]
     gamma4=COVV.b3*SS.b[3,]+COVV.b4*SS.b[4,]
     

     ### Sigma1.2.b
     eta1=gamma1*COVV.b1+gamma2*COVV.b2
     eta2=gamma1*COVV.b3+gamma2*COVV.b4
     eta3=gamma3*COVV.b1+gamma4*COVV.b2
     eta4=gamma3*COVV.b3+gamma4*COVV.b4
     
      mu0.EXE<-cbind(1,0,matrix(as.numeric(WE[jj.2,]),length(jj.2),dw))%*%coef.fix
       mu.EXE<-cbind(1,SE[i]-RE[jj.2],matrix(as.numeric(WE[jj.2,]),length(jj.2),dw))%*%coef.fix
       mut1.EXE<-cbind(1,t1,matrix(as.numeric(WE[jj.2,]),length(jj.2),dw))%*%coef.fix
       mut2.EXE<-cbind(1,t2,matrix(as.numeric(WE[jj.2,]),length(jj.2),dw))%*%coef.fix
       
      
      comp1<-parm1*(mu0.EXE+gamma1*Xt1+gamma2*Xt2-gamma1*mut1.EXE-gamma2*mut2.EXE)+
        parm2*(mu.EXE+gamma3*Xt1+gamma4*Xt2-gamma3*mut1.EXE-gamma4*mut2.EXE)
      
      comp2<-parm1^2*(coef.rd[1,1]-eta1)+parm2*parm1*(2*coef.rd[1,1]+2*(SE[i]-RE[jj.2])*coef.rd[1,2]-eta2-eta3)+
      parm2^2*(coef.rd[1,1]+(SE[i]-RE[jj.2])^2*coef.rd[2,2]+2*(SE[i]-RE[jj.2])*coef.rd[1,2]-eta4)
      
    
    ss=ss+sum(weights[jj.2]*exp(as.numeric(comp1)+0.5*comp2-comp12i+parm3*(SE[i]-RE[jj.2]-ZE[i]*TE[i])
    +parm0*(ZE[jj.2]-ZE[i])-t(as.numeric(WEc[i,])-matrix(as.numeric(t(WEc[jj.2,])),dc,length(jj.2)))%*%parm4) )  
    rm(Xt1,Xt2)
    }

    ## minus ties

     if (length(jj.2.tie)>0){
     Xt1=datwide.delta0[match(ID[jj.2.tie],datwide.delta0$id),Vname[1]]
     Xt2=datwide.delta0[match(ID[jj.2.tie],datwide.delta0$id),Vname[2]]
     
         ## make t1, t2 vectors
     t1=datwide.delta0[match(ID[jj.2.tie],datwide.delta0$id),Tname[1]]
     t2=datwide.delta0[match(ID[jj.2.tie],datwide.delta0$id),Tname[2]]
    

   ### when t1 t2 are vectors 
     COVV.b1<-coef.rd[1,1]+t1*coef.rd[1,2]
     COVV.b2<-coef.rd[1,1]+t2*coef.rd[1,2]
     COVV.b3<-coef.rd[1,1]+(t1+SE[i]-RE[jj.2.tie])*coef.rd[1,2]+t1*(SE[i]-RE[jj.2.tie])*coef.rd[2,2]
     COVV.b4<-coef.rd[1,1]+t2*(SE[i]-RE[jj.2.tie])*coef.rd[2,2]+(t2+SE[i]-RE[jj.2.tie])*coef.rd[1,2]
     

     SS.b.1<-coef.rd[1,1]+2*t1*coef.rd[1,2]+t1*2*coef.rd[2,2]+sigma.e^2
     SS.b.2<-coef.rd[1,1]+2*t2*coef.rd[1,2]+t2^2*coef.rd[2,2]+sigma.e^2
     SS.b.12<-coef.rd[1,1]+(t1+t2)*coef.rd[1,2]+t1*t2*coef.rd[2,2]
     
     SS.b<-array(NA,dim=c(2,2,length(jj.2.tie)))
     SS.b[1,1,]<-SS.b.1
     SS.b[1,2,]<-SS.b[2,1,]<-SS.b.12
     SS.b[2,2,]<-SS.b.2
     
     SS.b<-apply(SS.b,3,solve)  ### SS.b is a 4 by length(jj.2) matrix with each column the inverse of a matrix
 
     ### slope
    
     ### When t1 t2 are vectors , SS.b[1,], SS.b[2,], SS.b[3,], SS.b[4,] 
     ### corresponds [1,1],[2,1],[1,2],[2,2] of a 2x2 matrix
     
     gamma1=COVV.b1*SS.b[1,]+COVV.b2*SS.b[2,]
     gamma2=COVV.b1*SS.b[3,]+COVV.b2*SS.b[4,]
     gamma3=COVV.b3*SS.b[1,]+COVV.b4*SS.b[2,]
     gamma4=COVV.b3*SS.b[3,]+COVV.b4*SS.b[4,]
     
     ### Sigma1.2.b
     eta1=gamma1*COVV.b1+gamma2*COVV.b2
     eta2=gamma1*COVV.b3+gamma2*COVV.b4
     eta3=gamma3*COVV.b1+gamma4*COVV.b2
     eta4=gamma3*COVV.b3+gamma4*COVV.b4
     
      mu0.EXE<-cbind(1,0,matrix(as.numeric(WE[jj.2.tie,]),length(jj.2.tie),dw))%*%coef.fix
       mu.EXE<-cbind(1,SE[i]-RE[jj.2.tie],matrix(as.numeric(WE[jj.2.tie,]),length(jj.2.tie),dw))%*%coef.fix
       mut1.EXE<-cbind(1,t1,matrix(as.numeric(WE[jj.2.tie,]),length(jj.2.tie),dw))%*%coef.fix
       mut2.EXE<-cbind(1,t2,matrix(as.numeric(WE[jj.2.tie,]),length(jj.2.tie),dw))%*%coef.fix
       
      
      comp1<-parm1*(mu0.EXE+gamma1*Xt1+gamma2*Xt2-gamma1*mut1.EXE-gamma2*mut2.EXE)+
        parm2*(mu.EXE+gamma3*Xt1+gamma4*Xt2-gamma3*mut1.EXE-gamma4*mut2.EXE)
      
      comp2<-parm1^2*(coef.rd[1,1]-eta1)+parm2*parm1*(2*coef.rd[1,1]+2*(SE[i]-RE[jj.2.tie])*coef.rd[1,2]-eta2-eta3)+
      parm2^2*(coef.rd[1,1]+(SE[i]-RE[jj.2.tie])^2*coef.rd[2,2]+2*(SE[i]-RE[jj.2.tie])*coef.rd[1,2]-eta4)
      
    
    ss=ss-(ind.tie-1)/di*sum(weights[jj.2.tie]*exp(as.numeric(comp1)+0.5*comp2-comp12i+parm3*(SE[i]-RE[jj.2.tie]-ZE[i]*TE[i])
    +parm0*(ZE[jj.2.tie]-ZE[i])-t(as.numeric(WEc[i,])-matrix(as.numeric(t(WEc[jj.2.tie,])),dc,length(jj.2.tie)))%*%parm4) )  
    rm(Xt1,Xt2)
    }
    #####  ###### Measure of X(t1)*, X(t2)*, and X(t3)*
    
     if (length(jj.3)>0){
     Xt1=datwide.delta0[match(ID[jj.3],datwide.delta0$id),Vname[1]]
     Xt2=datwide.delta0[match(ID[jj.3],datwide.delta0$id),Vname[2]]
     Xt3=datwide.delta0[match(ID[jj.3],datwide.delta0$id),Vname[3]]
     
       ### t1, t2, t3 vector
     t1=datwide.delta0[match(ID[jj.3],datwide.delta0$id),Tname[1]]
     t2=datwide.delta0[match(ID[jj.3],datwide.delta0$id),Tname[2]]
     t3=datwide.delta0[match(ID[jj.3],datwide.delta0$id),Tname[3]]
     
   ####  t1,t2,t3 vector

     COVV.c1=coef.rd[1,1]+t1*coef.rd[1,2]
     COVV.c2=coef.rd[1,1]+t2*coef.rd[1,2]
     COVV.c3=coef.rd[1,1]+t3*coef.rd[1,2]
     COVV.c4=coef.rd[1,1]+(t1+SE[i]-RE[jj.3])*coef.rd[1,2]+t1*(SE[i]-RE[jj.3])*coef.rd[2,2]
     COVV.c5=coef.rd[1,1]+(t2+SE[i]-RE[jj.3])*coef.rd[1,2]+t2*(SE[i]-RE[jj.3])*coef.rd[2,2]
     COVV.c6=coef.rd[1,1]+(t3+SE[i]-RE[jj.3])*coef.rd[1,2]+t3*(SE[i]-RE[jj.3])*coef.rd[2,2]
     
     ##########
    
     SS.c.11<-coef.rd[1,1]+2*t1*coef.rd[1,2]+t1^2*coef.rd[2,2]+sigma.e^2
     SS.c.22<-coef.rd[1,1]+2*t2*coef.rd[1,2]+t2^2*coef.rd[2,2]+sigma.e^2
     SS.c.33<-coef.rd[1,1]+2*t3*coef.rd[1,2]+t3^2*coef.rd[2,2]+sigma.e^2
     SS.c.12=coef.rd[1,1]+(t1+t2)*coef.rd[1,2]+t1*t2*coef.rd[2,2]
     SS.c.13=coef.rd[1,1]+(t1+t3)*coef.rd[1,2]+t1*t3*coef.rd[2,2]
     SS.c.23=coef.rd[1,1]+(t2+t3)*coef.rd[1,2]+t2*t3*coef.rd[2,2]
     
     SS.c<-array(NA,dim=c(3,3,length(jj.3)))
     SS.c[1,1,]<-SS.c.11
     SS.c[2,2,]<-SS.c.22
     SS.c[3,3,]<-SS.c.33
     SS.c[1,2,]<-SS.c[2,1,]<-SS.c.12
     SS.c[1,3,]<-SS.c[3,1,]<-SS.c.13
     SS.c[2,3,]<-SS.c[3,2,]<-SS.c.23
     
    SS.c<-apply(SS.c,3,solve)
    

     #### SS.c is 9xlength(jj.3), elements in each column corresponds to [1,1],[2,1],[3,1],[1,2],[2,2],[3,2],[1,3],[2,3],[3,3] of a matrix

     delta1=COVV.c1*SS.c[1,]+COVV.c2*SS.c[2,]+COVV.c3*SS.c[3,]
     delta2=COVV.c1*SS.c[4,]+COVV.c2*SS.c[5,]+COVV.c3*SS.c[6,]
     delta3=COVV.c1*SS.c[7,]+COVV.c2*SS.c[8,]+COVV.c3*SS.c[9,]
     delta4=COVV.c4*SS.c[1,]+COVV.c5*SS.c[2,]+COVV.c6*SS.c[3,]
     delta5=COVV.c4*SS.c[4,]+COVV.c5*SS.c[5,]+COVV.c6*SS.c[6,]
     delta6=COVV.c4*SS.c[7,]+COVV.c5*SS.c[8,]+COVV.c6*SS.c[9,]

     xi1=delta1*COVV.c1+delta2*COVV.c2+delta3*COVV.c3
     xi2=delta1*COVV.c4+delta2*COVV.c5+delta3*COVV.c6
     xi3=delta4*COVV.c1+delta5*COVV.c2+delta6*COVV.c3
     xi4=delta4*COVV.c4+delta5*COVV.c5+delta6*COVV.c6
      
      mu0.EXE<-cbind(1,0,matrix(as.numeric(WE[jj.3,]),length(jj.3),dw))%*%coef.fix
       mu.EXE<-cbind(1,SE[i]-RE[jj.3],matrix(as.numeric(WE[jj.3,]),length(jj.3),dw))%*%coef.fix
       mut1.EXE<-cbind(1,t1,matrix(as.numeric(WE[jj.3,]),length(jj.3),dw))%*%coef.fix
       mut2.EXE<-cbind(1,t2,matrix(as.numeric(WE[jj.3,]),length(jj.3),dw))%*%coef.fix
       mut3.EXE<-cbind(1,t3,matrix(as.numeric(WE[jj.3,]),length(jj.3),dw))%*%coef.fix
       
      
      comp1<-parm1*(mu0.EXE+delta1*Xt1+delta2*Xt2+delta3*Xt3-delta1*mut1.EXE-delta2*mut2.EXE-delta3*mut3.EXE)+
        parm2*(mu.EXE+delta4*Xt1+delta5*Xt2+delta6*Xt3-delta4*mut1.EXE-delta5*mut2.EXE-delta6*mut3.EXE)
      
      comp2<-parm1^2*(coef.rd[1,1]-xi1)+parm2*parm1*(2*coef.rd[1,1]+2*(SE[i]-RE[jj.3])*coef.rd[1,2]-xi2-xi3)+
      parm2^2*(coef.rd[1,1]+(SE[i]-RE[jj.3])^2*coef.rd[2,2]+2*(SE[i]-RE[jj.3])*coef.rd[1,2]-xi4)
      
    
    ss=ss+sum(weights[jj.3]*exp(as.numeric(comp1)+0.5*comp2-comp12i+parm3*(SE[i]-RE[jj.3]-ZE[i]*TE[i])
    +parm0*(ZE[jj.3]-ZE[i])-t(as.numeric(WEc[i,])-matrix(as.numeric(t(WEc[jj.3,])),dc,length(jj.3)))%*%parm4) )  
   
    rm(Xt1,Xt2,Xt3)
    }

    #### Minus ties
 
     if (length(jj.3.tie)>0){
     Xt1=datwide.delta0[match(ID[jj.3.tie],datwide.delta0$id),Vname[1]]
     Xt2=datwide.delta0[match(ID[jj.3.tie],datwide.delta0$id),Vname[2]]
     Xt3=datwide.delta0[match(ID[jj.3.tie],datwide.delta0$id),Vname[3]]
     
      ### t1, t2, t3 vector
     t1=datwide.delta0[match(ID[jj.3.tie],datwide.delta0$id),Tname[1]]
     t2=datwide.delta0[match(ID[jj.3.tie],datwide.delta0$id),Tname[2]]
     t3=datwide.delta0[match(ID[jj.3.tie],datwide.delta0$id),Tname[3]]
 
     COVV.c1=coef.rd[1,1]+t1*coef.rd[1,2]
     COVV.c2=coef.rd[1,1]+t2*coef.rd[1,2]
     COVV.c3=coef.rd[1,1]+t3*coef.rd[1,2]
     COVV.c4=coef.rd[1,1]+(t1+SE[i]-RE[jj.3.tie])*coef.rd[1,2]+t1*(SE[i]-RE[jj.3.tie])*coef.rd[2,2]
     COVV.c5=coef.rd[1,1]+(t2+SE[i]-RE[jj.3.tie])*coef.rd[1,2]+t2*(SE[i]-RE[jj.3.tie])*coef.rd[2,2]
     COVV.c6=coef.rd[1,1]+(t3+SE[i]-RE[jj.3.tie])*coef.rd[1,2]+t3*(SE[i]-RE[jj.3.tie])*coef.rd[2,2]
     
     ##########
    
     SS.c.11<-coef.rd[1,1]+2*t1*coef.rd[1,2]+t1^2*coef.rd[2,2]+sigma.e^2
     SS.c.22<-coef.rd[1,1]+2*t2*coef.rd[1,2]+t2^2*coef.rd[2,2]+sigma.e^2
     SS.c.33<-coef.rd[1,1]+2*t3*coef.rd[1,2]+t3^2*coef.rd[2,2]+sigma.e^2
     SS.c.12=coef.rd[1,1]+(t1+t2)*coef.rd[1,2]+t1*t2*coef.rd[2,2]
     SS.c.13=coef.rd[1,1]+(t1+t3)*coef.rd[1,2]+t1*t3*coef.rd[2,2]
     SS.c.23=coef.rd[1,1]+(t2+t3)*coef.rd[1,2]+t2*t3*coef.rd[2,2]
     
     SS.c<-array(NA,dim=c(3,3,length(jj.3.tie)))
     SS.c[1,1,]<-SS.c.11
     SS.c[2,2,]<-SS.c.22
     SS.c[3,3,]<-SS.c.33
     SS.c[1,2,]<-SS.c[2,1,]<-SS.c.12
     SS.c[1,3,]<-SS.c[3,1,]<-SS.c.13
     SS.c[2,3,]<-SS.c[3,2,]<-SS.c.23
         
    SS.c<-apply(SS.c,3,solve)
    

     #### SS.c is 9xlength(jj.3), elements in each column corresponds to [1,1],[2,1],[3,1],[1,2],[2,2],[3,2],[1,3],[2,3],[3,3] of a matrix

     delta1=COVV.c1*SS.c[1,]+COVV.c2*SS.c[2,]+COVV.c3*SS.c[3,]
     delta2=COVV.c1*SS.c[4,]+COVV.c2*SS.c[5,]+COVV.c3*SS.c[6,]
     delta3=COVV.c1*SS.c[7,]+COVV.c2*SS.c[8,]+COVV.c3*SS.c[9,]
     delta4=COVV.c4*SS.c[1,]+COVV.c5*SS.c[2,]+COVV.c6*SS.c[3,]
     delta5=COVV.c4*SS.c[4,]+COVV.c5*SS.c[5,]+COVV.c6*SS.c[6,]
     delta6=COVV.c4*SS.c[7,]+COVV.c5*SS.c[8,]+COVV.c6*SS.c[9,]

     
     xi1=delta1*COVV.c1+delta2*COVV.c2+delta3*COVV.c3
     xi2=delta1*COVV.c4+delta2*COVV.c5+delta3*COVV.c6
     xi3=delta4*COVV.c1+delta5*COVV.c2+delta6*COVV.c3
     xi4=delta4*COVV.c4+delta5*COVV.c5+delta6*COVV.c6
      
      mu0.EXE<-cbind(1,0,matrix(as.numeric(WE[jj.3.tie,]),length(jj.3.tie),dw))%*%coef.fix
       mu.EXE<-cbind(1,SE[i]-RE[jj.3.tie],matrix(as.numeric(WE[jj.3.tie,]),length(jj.3.tie),dw))%*%coef.fix
       mut1.EXE<-cbind(1,t1,matrix(as.numeric(WE[jj.3.tie,]),length(jj.3.tie),dw))%*%coef.fix
       mut2.EXE<-cbind(1,t2,matrix(as.numeric(WE[jj.3.tie,]),length(jj.3.tie),dw))%*%coef.fix
       mut3.EXE<-cbind(1,t3,matrix(as.numeric(WE[jj.3.tie,]),length(jj.3.tie),dw))%*%coef.fix
            
      comp1<-parm1*(mu0.EXE+delta1*Xt1+delta2*Xt2+delta3*Xt3-delta1*mut1.EXE-delta2*mut2.EXE-delta3*mut3.EXE)+
        parm2*(mu.EXE+delta4*Xt1+delta5*Xt2+delta6*Xt3-delta4*mut1.EXE-delta5*mut2.EXE-delta6*mut3.EXE)
      
      comp2<-parm1^2*(coef.rd[1,1]-xi1)+parm2*parm1*(2*coef.rd[1,1]+2*(SE[i]-RE[jj.3.tie])*coef.rd[1,2]-xi2-xi3)+
      parm2^2*(coef.rd[1,1]+(SE[i]-RE[jj.3.tie])^2*coef.rd[2,2]+2*(SE[i]-RE[jj.3.tie])*coef.rd[1,2]-xi4)
      
    ss=ss-(ind.tie-1)/di*sum(weights[jj.3.tie]*exp(as.numeric(comp1)+0.5*comp2-comp12i+parm3*(SE[i]-RE[jj.3.tie]-ZE[i]*TE[i])
    +parm0*(ZE[jj.3.tie]-ZE[i])-t(as.numeric(WEc[i,])-matrix(as.numeric(t(WEc[jj.3.tie,])),dc,length(jj.3.tie)))%*%parm4) )  
   
    rm(Xt1,Xt2,Xt3)
    }

    ss.new=ss
    ### placebo
   
   ### need to use jj.Z0 and jj.Z0.tie when there are confounders in the risk model
   if (length(jj.Z0)>0)
   ss=ss+sum(weights[jj.Z0]*exp(-comp12i-parm3*(ZE[i]*TE[i])-parm0*(ZE[i])
   -t(as.numeric(WEc[i,])-matrix(as.numeric(t(WEc[jj.Z0,])),dc,length(jj.Z0)))%*%parm4))
   
   if (length(jj.Z0.tie)>0)
   ss=ss-(ind.tie-1)/di*sum(weights[jj.Z0.tie]*exp(-comp12i-parm3*(ZE[i]*TE[i])-parm0*(ZE[i])-t(as.numeric(WEc[i,])-matrix(as.numeric(t(WEc[jj.Z0.tie,])),dc,length(jj.Z0.tie)))%*%parm4))
   
    out=out+weights.bar*log(weights[i]+ss)+(weights.bar-weights[i])*(comp12i+parm3*(SE[i]-RE[i])*ZE[i]+parm0*(ZE[i])
    +matrix(WEc[i,],length(i),dc)%*%parm4)
    
    outi=c(outi,ss)
 
    }
 return(out)
 }

########

######## Fit a Cox hazard model on calendar time with coefficient for Z, Z*x(t), Z*t, and W

nloglik.M1.b.Cal.23<- function(parm){

 out=0
 
 outi<-NULL
   ### coef for Z, Z*x(t), Z*t, and confounders respectively
    parm0=parm[1]
    #parm1=parm[2]
    parm2=parm[2]
    parm3=parm[3]
    parm4=parm[4]
   
   jj.4.full<-(1:length(RE))[datwide.delta0[match(ID,datwide.delta0$id),]$groupall > datwide.delta0[match(ID,datwide.delta0$id),]$group]
   jj.4.full<-as.numeric(na.omit(jj.4.full))
   out.jj.4.full<-NULL
   for (j in jj.4.full){
      
    keep=1:ll
 
    VVj.11<-coef.rd[1,1]
    VVj.12<-VVj.21<-coef.rd[1,1]+(SE[keep]-RE[j])*coef.rd[1,2]
    VVj.22<-coef.rd[1,1]+2*(SE[keep]-RE[j])*coef.rd[1,2]+(SE[keep]-RE[j])^2*coef.rd[2,2]
    
    muj.1=cbind(1,0,matrix(WE[j,],length(j),dw))%*%coef.fix
    muj.2=cbind(1,SE[keep]-RE[j],matrix(WE[j,],length(j),dw))%*%coef.fix
    
   
    oo.j<-id.long==ID[j]
    Xj<-X.long[oo.j]
    Tj<-T.long[oo.j]
    
    
    VVj.c<-cbind(1,Tj)%*%coef.rd%*%rbind(1,Tj)+diag(rep(sigma.e^2,length(Tj)),nrow=length(Tj),ncol=length(Tj))               
   
    COVVj.1<-coef.rd[1,1]+coef.rd[1,2]*Tj
    COVVj.2<-coef.rd[1,1]+coef.rd[1,2]*outer(SE[keep]-RE[j],Tj,'+')+coef.rd[2,2]*outer(SE[keep]-RE[j],Tj,'*')
    COVVj.1=matrix(COVVj.1,1,length(COVVj.1))
    if (length(Tj)==1) COVVj.2=matrix(COVVj.2,length(COVVj.2),1)
    
    muj.c<-cbind(1,Tj,matrix(WE[j,],length(j),dw))%*%coef.fix
        
    comp1j<-cbind(as.numeric(muj.1),muj.2)+cbind(as.numeric(COVVj.1%*%solve(VVj.c)%*%(Xj-muj.c)),COVVj.2%*%solve(VVj.c)%*%(Xj-muj.c))
    
     VVjminus.11<-COVVj.1%*%solve(VVj.c)%*%(t(COVVj.1))
     VVjminus.12<-COVVj.2%*%solve(VVj.c)%*%(t(COVVj.1))
     VVjminus.22<-diag(COVVj.2%*%solve(VVj.c)%*%(t(COVVj.2)))
     
     comp2j.11<-VVj.11-VVjminus.11
     comp2j.12<-VVj.12-VVjminus.12
     comp2j.22<-VVj.22-VVjminus.22
   
     out.jj.4.full<-cbind(out.jj.4.full,as.numeric(c(parm1,parm2)%*%t(comp1j))+parm1^2*as.numeric(comp2j.11)/2+parm1*parm2*comp2j.12+parm2^2*comp2j.22/2+parm3*(SE[keep]-RE[j])
     +parm0*ZE[j]+matrix(WEc[j,],length(j),dc)%*%parm4)
     
    }
        
    S.tie=SE.u[1]  ### unique SE value
    ind.tie=0       ### position within the unique value
    
    for (i in 1:ll){
    
    ## track #of ties
    
        di<-sum(SE==SE[i] & delta==1)
        
        if (SE[i]==S.tie){  ## in the same tie as previous value
          ind.tie=ind.tie+1
        } else { ## go to a new SE value
           S.tie=SE[i]
           ind.tie=1
        }
     ### compute average weight among ties
        weights.bar=mean(weights[SE==SE[i] & delta==1])
     
    ss=0
    comp12i=0 ### placebo contribute nothing
    if (ZE[i]==1) { 
    VVi<-cbind(1,c(0,TE[i]))%*%coef.rd%*%rbind(1,c(0,TE[i]))
    mui<-rbind(c(1,0,matrix(WE[i,],length(i),dw)),c(1,TE[i],matrix(WE[i,],length(i),dw)))%*%coef.fix
   
    oo.i<-id.long==ID[i]
        
       #### (i) if case has at least one X(Ti)* measured, i.e. case ID[i] belong to id.long
        if (sum(oo.i)>0){ 
        Xi<-X.long[oo.i]
        Ti<-T.long[oo.i]
        
        VVi.c<-cbind(1,Ti)%*%coef.rd%*%rbind(1,Ti)+diag(rep(sigma.e^2,length(Ti)),nrow=length(Ti),ncol=length(Ti))               
        COVVi<-cbind(1,c(0,TE[i]))%*%coef.rd%*%rbind(1,Ti)
        mui.c<-cbind(1,Ti,matrix(WE[i,],length(Ti),dw,byrow=TRUE))%*%coef.fix
        
        comp1i<-mui+COVVi%*%solve(VVi.c)%*%(Xi-mui.c)      ## E(X(0),X(t)|X(Ti)*,W)
        comp2i<-VVi-COVVi%*%solve(VVi.c)%*%t(COVVi)        ## Var(X(0),X(t)|X(Ti)*,W)
        
        } else {
        #### (ii) if case does not have X(Ti)* measure, i.e. case ID[i] not in id.long
        
       comp1i<-mui            
       comp2i<-VVi            
       } 
 
     comp12i<-as.numeric(t(c(parm1,parm2))%*%comp1i+t(c(parm1,parm2))%*%comp2i%*%c(parm1,parm2)/2)
    
    } 
    
    jj.r<-((i+1):(length(RE)))[RE[(i+1):length(RE)]<SE[i] & SE[(i+1):length(RE)]>=SE[i] & ID[(i+1):length(RE)]%in%id.long]
    jj.l<-NULL
    if (i>1) jj.l<-(1:(i-1))[RE[1:(i-1)]<SE[i] & SE[1:(i-1)]>=SE[i] & ID[1:(i-1)]%in%id.long]
    jj<-c(jj.l,jj.r) ## among risk set  
    
    jj.tie<-(1:(length(RE)))[RE[1:length(RE)]<SE[i] & ID[1:length(RE)]%in%id.long & SE[1:length(SE)]==SE[i] & delta[1:length(delta)]==1]
   
    jj.W.r<-((i+1):(length(RE)))[RE[(i+1):length(RE)]<SE[i] & ZE[(i+1):length(ZE)]==1 & SE[(i+1):length(RE)]>=SE[i] & !ID[(i+1):length(RE)]%in%id.long]
    jj.W.l<-NULL
    if (i>1) jj.W.l<-(1:(i-1))[RE[1:(i-1)]<SE[i] & ZE[1:(i-1)]==1 & !ID[1:(i-1)]%in%id.long & SE[1:(i-1)]>=SE[i]]
    jj.W<-c(jj.W.l,jj.W.r)    

    
    jj.W.tie<-(1:(length(RE)))[RE[1:length(RE)]<SE[i] & ZE[1:length(ZE)]==1 & SE[1:length(RE)]==SE[i] & !ID[1:length(RE)]%in%id.long  & delta[1:length(delta)]==1]
 
  
    jj.Z0.r<-((i+1):length(RE))[RE[(i+1):length(RE)]<SE[i] & ZE[(i+1):length(ZE)]==0 & SE[(i+1):length(SE)]>=SE[i]]  
    jj.Z0.l<-NULL
    if (i>1) jj.Z0.l<-(1:(i-1))[RE[1:(i-1)]<SE[i] & ZE[1:(i-1)]==0 & SE[1:(i-1)]>=SE[i]]
    jj.Z0=c(jj.Z0.l,jj.Z0.r)

    jj.Z0.tie<-(1:length(RE))[RE[1:length(RE)]<SE[i] & ZE[1:length(ZE)]==0 & SE[1:length(SE)]==SE[i] & delta[1:length(delta)]==1]  
  
    jj.1.r<-((i+1):(length(RE)))[RE[(i+1):length(RE)]<SE[i] & datwide.delta0[match(ID[(i+1):length(RE)],datwide.delta0$id),]$group==1 &
       datwide.delta0[match(ID[(i+1):length(RE)],datwide.delta0$id),]$groupall==1 & SE[(i+1):length(RE)]>=SE[i]]
    jj.1.l<-NULL
    if (i>1) jj.1.l<-(1:(i-1))[RE[1:(i-1)]<SE[i] & datwide.delta0[match(ID[1:(i-1)],datwide.delta0$id),]$group==1 &
       datwide.delta0[match(ID[1:(i-1)],datwide.delta0$id),]$groupall==1 & SE[1:(i-1)]>=SE[i]]
    jj.1<-c(jj.1.l,jj.1.r) ## among risk set

    
    jj.1.tie<-(1:(length(RE)))[RE[1:length(RE)]<SE[i] & datwide.delta0[match(ID[1:length(RE)],datwide.delta0$id),]$group==1 &
       datwide.delta0[match(ID[1:length(RE)],datwide.delta0$id),]$groupall==1 & SE[1:length(SE)]==SE[i] & delta[1:length(delta)]==1]
   
   
    #######
 
    jj.2.r<-((i+1):(length(RE)))[RE[(i+1):length(RE)]<SE[i] & datwide.delta0[match(ID[(i+1):length(RE)],datwide.delta0$id),]$group==2 & 
       datwide.delta0[match(ID[(i+1):length(RE)],datwide.delta0$id),]$groupall==2 & SE[(i+1):length(SE)]>=SE[i]]
    jj.2.l<-NULL
    if (i>1)   jj.2.l<-(1:(i-1))[RE[1:(i-1)]<SE[i] & datwide.delta0[match(ID[1:(i-1)],datwide.delta0$id),]$group==2 & 
       datwide.delta0[match(ID[1:(i-1)],datwide.delta0$id),]$groupall==2 & SE[1:(i-1)]>=SE[i]]
    jj.2<-c(jj.2.l,jj.2.r)

    jj.2.tie<-(1:(length(RE)))[RE[1:length(RE)]<SE[i] & datwide.delta0[match(ID[1:length(RE)],datwide.delta0$id),]$group==2 & 
       datwide.delta0[match(ID[1:length(RE)],datwide.delta0$id),]$groupall==2 & SE[1:length(SE)]==SE[i] & delta[1:length(delta)]==1]
   

    ###
    jj.3.r<-((i+1):(length(RE)))[RE[(i+1):length(RE)]<SE[i] & datwide.delta0[match(ID[(i+1):length(RE)],datwide.delta0$id),]$group==3 & 
       datwide.delta0[match(ID[(i+1):length(RE)],datwide.delta0$id),]$groupall==3 & SE[(i+1):length(SE)]>=SE[i]]
    jj.3.l<-NULL
    if (i>1) jj.3.l<-(1:(i-1))[RE[1:(i-1)]<SE[i] & datwide.delta0[match(ID[1:(i-1)],datwide.delta0$id),]$group==3 & 
       datwide.delta0[match(ID[1:(i-1)],datwide.delta0$id),]$groupall==3 & SE[1:(i-1)]>=SE[i]]
    jj.3<-c(jj.3.l,jj.3.r)
   

    jj.3.tie<-(1:(length(RE)))[RE[1:length(RE)]<SE[i] & datwide.delta0[match(ID[1:length(RE)],datwide.delta0$id),]$group==3 & 
       datwide.delta0[match(ID[1:length(RE)],datwide.delta0$id),]$groupall==3 & SE[1:length(SE)]==SE[i] & delta[1:length(delta)]==1]
  
     jj.1<-as.numeric(na.omit(jj.1))
     jj.2<-as.numeric(na.omit(jj.2))
     jj.3<-as.numeric(na.omit(jj.3))
   
       jj.1.tie<-as.numeric(na.omit(jj.1.tie))
     jj.2.tie<-as.numeric(na.omit(jj.2.tie))
     jj.3.tie<-as.numeric(na.omit(jj.3.tie))
     
    ss=ss+sum((RE[jj.4.full]<SE[i])*(SE[jj.4.full]>=SE[i])*(!jj.4.full==i)*weights[jj.4.full]*exp(out.jj.4.full[i,]-comp12i-parm3*(ZE[i]*TE[i])
    -parm0*ZE[i]-matrix(WEc[i,],length(i),dc)%*%parm4))  ## risk set from jj.4.full
     
    ss=ss-(ind.tie-1)/di*sum((RE[jj.4.full]<SE[i])*(SE[jj.4.full]==SE[i])*(!jj.4.full==i)*weights[jj.4.full]*exp(out.jj.4.full[i,]-comp12i-parm3*(ZE[i]*TE[i])
    -parm0*ZE[i]-matrix(WEc[i,],length(i),dc)%*%parm4))  ## tie from jj.4.full
     

   if (length(jj.W)>0){
      ### No X measure 
    mu0.EXE<-cbind(1,0,matrix(as.numeric(WE[jj.W,]),length(jj.W),dw))%*%coef.fix
    mu.EXE<-cbind(1,SE[i]-RE[jj.W],matrix(as.numeric(WE[jj.W,]),length(jj.W),dw))%*%coef.fix
  
   comp1.b<-parm1*mu0.EXE+parm2*mu.EXE
       
   comp2.b<-parm1^2*coef.rd[1,1]+2*parm1*parm2*(coef.rd[1,1]+(SE[i]-RE[jj.W])*coef.rd[1,2])+
   parm2^2*(coef.rd[1,1]+2*(SE[i]-RE[jj.W])*coef.rd[1,2]+(SE[i]-RE[jj.W])^2*coef.rd[2,2])
   
   ss=ss+sum(weights[jj.W]*exp(as.numeric(comp1.b)+0.5*comp2.b-comp12i+parm3*(SE[i]-RE[jj.W]-ZE[i]*TE[i])
   +parm0*(ZE[jj.W]-ZE[i])-t(as.numeric(WEc[i,])-matrix(as.numeric(t(WEc[jj.W,])),dc,length(jj.W)))%*%parm4) )  
    
   } 
   ### minus ties
   if (length(jj.W.tie)>0){
   
        mu0.EXE<-cbind(1,0,matrix(as.numeric(WE[jj.W.tie,]),length(jj.W.tie),dw))%*%coef.fix
        mu.EXE<-cbind(1,SE[i]-RE[jj.W.tie],matrix(as.numeric(WE[jj.W.tie,]),length(jj.W.tie),dw))%*%coef.fix
      
       comp1.b<-parm1*mu0.EXE+parm2*mu.EXE
           
       comp2.b<-parm1^2*coef.rd[1,1]+2*parm1*parm2*(coef.rd[1,1]+(SE[i]-RE[jj.W.tie])*coef.rd[1,2])+
       parm2^2*(coef.rd[1,1]+2*(SE[i]-RE[jj.W.tie])*coef.rd[1,2]+(SE[i]-RE[jj.W.tie])^2*coef.rd[2,2])
       
       ss=ss-(ind.tie-1)/di*sum(weights[jj.W.tie]*exp(as.numeric(comp1.b)+0.5*comp2.b-comp12i+parm3*(SE[i]-RE[jj.W.tie]-ZE[i]*TE[i])
       +parm0*(ZE[jj.W.tie]-ZE[i])-t(as.numeric(WEc[i,])-matrix(as.numeric(t(WEc[jj.W.tie,])),dc,length(jj.W.tie)))%*%parm4) ) 
   }      

 
   ### Measure of X(0)* only
   
    #account for different t levels
   
   if (length(jj.1)>0){
   Xt1=datwide.delta0[match(ID[jj.1],datwide.delta0$id),Vname[1]]
    t1=datwide.delta0[match(ID[jj.1],datwide.delta0$id),Tname[1]]
 
    ## make t1 vector
   COVV.a<-matrix(c(coef.rd[1,1]+t1*coef.rd[1,2],coef.rd[1,1]+(t1+SE[i]-RE[jj.1])*coef.rd[1,2]+t1*(SE[i]-RE[jj.1])*coef.rd[2,2]),length(jj.1),2,byrow=F)
   SS.a<-sapply(coef.rd[1,1]+2*t1*coef.rd[1,2]+t1^2*coef.rd[2,2]+sigma.e^2,solve)
   
   slope.a<-COVV.a*as.numeric(SS.a)
   Sigma1.2.a=outer(t(slope.a),(COVV.a),'*'); 
   kappa1=diag(as.matrix(Sigma1.2.a[1,,,1])); kappa2=diag(as.matrix(Sigma1.2.a[1,,,2])); 
   kappa3=diag(as.matrix(Sigma1.2.a[2,,,1])); kappa4=diag(as.matrix(Sigma1.2.a[2,,,2]))
   
     mu0.EXE<-cbind(1,0,matrix(as.numeric(WE[jj.1,]),length(jj.1),dw))%*%coef.fix
      mu.EXE<-cbind(1,SE[i]-RE[jj.1],matrix(as.numeric(WE[jj.1,]),length(jj.1),dw))%*%coef.fix
  
       comp1<-parm1*(mu0.EXE+slope.a[,1]*(Xt1-mu0.EXE))+
       parm2*(mu.EXE+slope.a[,2]*(Xt1-mu0.EXE))
       
       comp2<-parm1^2*(coef.rd[1,1]-kappa1)+
       parm2*parm1*(2*coef.rd[1,1]+2*(SE[i]-RE[jj.1])*coef.rd[1,2]-kappa2-kappa3)+
       parm2^2*(coef.rd[1,1]+(SE[i]-RE[jj.1])^2*coef.rd[2,2]+2*(SE[i]-RE[jj.1])*coef.rd[1,2]-kappa4)
    
    ss=ss+sum(weights[jj.1]*exp(as.numeric(comp1)+0.5*comp2-comp12i+parm3*(SE[i]-RE[jj.1]-ZE[i]*TE[i])
    +parm0*(ZE[jj.1]-ZE[i])-t(as.numeric(WEc[i,])-matrix(as.numeric(t(WEc[jj.1,])),dc,length(jj.1)))%*%parm4) )  
   
   rm(Xt1)
   }

   ## minus ties

   if (length(jj.1.tie)>0){
   Xt1=datwide.delta0[match(ID[jj.1.tie],datwide.delta0$id),Vname[1]]
   
   t1=datwide.delta0[match(ID[jj.1.tie],datwide.delta0$id),Tname[1]]
   
    # make t1 vector
   COVV.a<-matrix(c(coef.rd[1,1]+t1*coef.rd[1,2],coef.rd[1,1]+(t1+SE[i]-RE[jj.1.tie])*coef.rd[1,2]+t1*(SE[i]-RE[jj.1.tie])*coef.rd[2,2]),length(jj.1.tie),2,byrow=F)
   SS.a<-sapply(coef.rd[1,1]+2*t1*coef.rd[1,2]+t1^2*coef.rd[2,2]+sigma.e^2,solve)
 
   
   slope.a<-COVV.a*as.numeric(SS.a)
   Sigma1.2.a=outer(t(slope.a),(COVV.a),'*'); 
   kappa1=diag(as.matrix(Sigma1.2.a[1,,,1])); kappa2=diag(as.matrix(Sigma1.2.a[1,,,2])); 
   kappa3=diag(as.matrix(Sigma1.2.a[2,,,1])); kappa4=diag(as.matrix(Sigma1.2.a[2,,,2]))
   
     mu0.EXE<-cbind(1,0,matrix(as.numeric(WE[jj.1.tie,]),length(jj.1.tie),dw))%*%coef.fix
      mu.EXE<-cbind(1,SE[i]-RE[jj.1.tie],matrix(as.numeric(WE[jj.1.tie,]),length(jj.1.tie),dw))%*%coef.fix
  
       comp1<-parm1*(mu0.EXE+slope.a[,1]*(Xt1-mu0.EXE))+
       parm2*(mu.EXE+slope.a[,2]*(Xt1-mu0.EXE))
       
       comp2<-parm1^2*(coef.rd[1,1]-kappa1)+
       parm2*parm1*(2*coef.rd[1,1]+2*(SE[i]-RE[jj.1.tie])*coef.rd[1,2]-kappa2-kappa3)+
       parm2^2*(coef.rd[1,1]+(SE[i]-RE[jj.1.tie])^2*coef.rd[2,2]+2*(SE[i]-RE[jj.1.tie])*coef.rd[1,2]-kappa4)
    
    ss=ss-(ind.tie-1)/di*sum(weights[jj.1.tie]*exp(as.numeric(comp1)+0.5*comp2-comp12i+parm3*(SE[i]-RE[jj.1.tie]-ZE[i]*TE[i])
    +parm0*(ZE[jj.1.tie]-ZE[i])-t(as.numeric(WEc[i,])-matrix(as.numeric(t(WEc[jj.1.tie,])),dc,length(jj.1.tie)))%*%parm4 ) )  
   
   rm(Xt1)
   }



   ###### Measure of X(0)* and X(t2)* 
     if (length(jj.2)>0){
     Xt1=datwide.delta0[match(ID[jj.2],datwide.delta0$id),Vname[1]]
     Xt2=datwide.delta0[match(ID[jj.2],datwide.delta0$id),Vname[2]]
     
       ## make t1, t2 vectors
      t1=datwide.delta0[match(ID[jj.2],datwide.delta0$id),Tname[1]]
      t2=datwide.delta0[match(ID[jj.2],datwide.delta0$id),Tname[2]]
         
     COVV.b1<-coef.rd[1,1]+t1*coef.rd[1,2]
     COVV.b2<-coef.rd[1,1]+t2*coef.rd[1,2]
     COVV.b3<-coef.rd[1,1]+(t1+SE[i]-RE[jj.2])*coef.rd[1,2]+t1*(SE[i]-RE[jj.2])*coef.rd[2,2]
     COVV.b4<-coef.rd[1,1]+t2*(SE[i]-RE[jj.2])*coef.rd[2,2]+(t2+SE[i]-RE[jj.2])*coef.rd[1,2]
     

     SS.b.1<-coef.rd[1,1]+2*t1*coef.rd[1,2]+t1*2*coef.rd[2,2]+sigma.e^2
     SS.b.2<-coef.rd[1,1]+2*t2*coef.rd[1,2]+t2^2*coef.rd[2,2]+sigma.e^2
     SS.b.12<-coef.rd[1,1]+(t1+t2)*coef.rd[1,2]+t1*t2*coef.rd[2,2]
     
     SS.b<-array(NA,dim=c(2,2,length(jj.2)))
     SS.b[1,1,]<-SS.b.1
     SS.b[1,2,]<-SS.b[2,1,]<-SS.b.12
     SS.b[2,2,]<-SS.b.2
     
     SS.b<-apply(SS.b,3,solve)  ### SS.b is a 4 by length(jj.2) matrix with each column the inverse of a matrix
 
     ### slope
    
     ### When t1 t2 are vectors, SS.b[1,], SS.b[2,], SS.b[3,], SS.b[4,] 
     ### corresponds [1,1],[2,1],[1,2],[2,2] of a 2x2 matrix
     
     gamma1=COVV.b1*SS.b[1,]+COVV.b2*SS.b[2,]
     gamma2=COVV.b1*SS.b[3,]+COVV.b2*SS.b[4,]
     gamma3=COVV.b3*SS.b[1,]+COVV.b4*SS.b[2,]
     gamma4=COVV.b3*SS.b[3,]+COVV.b4*SS.b[4,]
     

     ### Sigma1.2.b
     eta1=gamma1*COVV.b1+gamma2*COVV.b2
     eta2=gamma1*COVV.b3+gamma2*COVV.b4
     eta3=gamma3*COVV.b1+gamma4*COVV.b2
     eta4=gamma3*COVV.b3+gamma4*COVV.b4
     
      mu0.EXE<-cbind(1,0,matrix(as.numeric(WE[jj.2,]),length(jj.2),dw))%*%coef.fix
       mu.EXE<-cbind(1,SE[i]-RE[jj.2],matrix(as.numeric(WE[jj.2,]),length(jj.2),dw))%*%coef.fix
       mut1.EXE<-cbind(1,t1,matrix(as.numeric(WE[jj.2,]),length(jj.2),dw))%*%coef.fix
       mut2.EXE<-cbind(1,t2,matrix(as.numeric(WE[jj.2,]),length(jj.2),dw))%*%coef.fix
       
      
      comp1<-parm1*(mu0.EXE+gamma1*Xt1+gamma2*Xt2-gamma1*mut1.EXE-gamma2*mut2.EXE)+
        parm2*(mu.EXE+gamma3*Xt1+gamma4*Xt2-gamma3*mut1.EXE-gamma4*mut2.EXE)
      
      comp2<-parm1^2*(coef.rd[1,1]-eta1)+parm2*parm1*(2*coef.rd[1,1]+2*(SE[i]-RE[jj.2])*coef.rd[1,2]-eta2-eta3)+
      parm2^2*(coef.rd[1,1]+(SE[i]-RE[jj.2])^2*coef.rd[2,2]+2*(SE[i]-RE[jj.2])*coef.rd[1,2]-eta4)
      
    
    ss=ss+sum(weights[jj.2]*exp(as.numeric(comp1)+0.5*comp2-comp12i+parm3*(SE[i]-RE[jj.2]-ZE[i]*TE[i])
    +parm0*(ZE[jj.2]-ZE[i])-t(as.numeric(WEc[i,])-matrix(as.numeric(t(WEc[jj.2,])),dc,length(jj.2)))%*%parm4) )  
    rm(Xt1,Xt2)
    }

    ## minus ties

     if (length(jj.2.tie)>0){
     Xt1=datwide.delta0[match(ID[jj.2.tie],datwide.delta0$id),Vname[1]]
     Xt2=datwide.delta0[match(ID[jj.2.tie],datwide.delta0$id),Vname[2]]
     
         ## make t1, t2 vectors
     t1=datwide.delta0[match(ID[jj.2.tie],datwide.delta0$id),Tname[1]]
     t2=datwide.delta0[match(ID[jj.2.tie],datwide.delta0$id),Tname[2]]
    

   ### when t1 t2 are vectors 
     COVV.b1<-coef.rd[1,1]+t1*coef.rd[1,2]
     COVV.b2<-coef.rd[1,1]+t2*coef.rd[1,2]
     COVV.b3<-coef.rd[1,1]+(t1+SE[i]-RE[jj.2.tie])*coef.rd[1,2]+t1*(SE[i]-RE[jj.2.tie])*coef.rd[2,2]
     COVV.b4<-coef.rd[1,1]+t2*(SE[i]-RE[jj.2.tie])*coef.rd[2,2]+(t2+SE[i]-RE[jj.2.tie])*coef.rd[1,2]
     

     SS.b.1<-coef.rd[1,1]+2*t1*coef.rd[1,2]+t1*2*coef.rd[2,2]+sigma.e^2
     SS.b.2<-coef.rd[1,1]+2*t2*coef.rd[1,2]+t2^2*coef.rd[2,2]+sigma.e^2
     SS.b.12<-coef.rd[1,1]+(t1+t2)*coef.rd[1,2]+t1*t2*coef.rd[2,2]
     
     SS.b<-array(NA,dim=c(2,2,length(jj.2.tie)))
     SS.b[1,1,]<-SS.b.1
     SS.b[1,2,]<-SS.b[2,1,]<-SS.b.12
     SS.b[2,2,]<-SS.b.2
     
     SS.b<-apply(SS.b,3,solve)  ### SS.b is a 4 by length(jj.2) matrix with each column the inverse of a matrix
 
     ### slope
    
     ### When t1 t2 are vectors , SS.b[1,], SS.b[2,], SS.b[3,], SS.b[4,] 
     ### corresponds [1,1],[2,1],[1,2],[2,2] of a 2x2 matrix
     
     gamma1=COVV.b1*SS.b[1,]+COVV.b2*SS.b[2,]
     gamma2=COVV.b1*SS.b[3,]+COVV.b2*SS.b[4,]
     gamma3=COVV.b3*SS.b[1,]+COVV.b4*SS.b[2,]
     gamma4=COVV.b3*SS.b[3,]+COVV.b4*SS.b[4,]
     
     ### Sigma1.2.b
     eta1=gamma1*COVV.b1+gamma2*COVV.b2
     eta2=gamma1*COVV.b3+gamma2*COVV.b4
     eta3=gamma3*COVV.b1+gamma4*COVV.b2
     eta4=gamma3*COVV.b3+gamma4*COVV.b4
     
      mu0.EXE<-cbind(1,0,matrix(as.numeric(WE[jj.2.tie,]),length(jj.2.tie),dw))%*%coef.fix
       mu.EXE<-cbind(1,SE[i]-RE[jj.2.tie],matrix(as.numeric(WE[jj.2.tie,]),length(jj.2.tie),dw))%*%coef.fix
       mut1.EXE<-cbind(1,t1,matrix(as.numeric(WE[jj.2.tie,]),length(jj.2.tie),dw))%*%coef.fix
       mut2.EXE<-cbind(1,t2,matrix(as.numeric(WE[jj.2.tie,]),length(jj.2.tie),dw))%*%coef.fix
       
      
      comp1<-parm1*(mu0.EXE+gamma1*Xt1+gamma2*Xt2-gamma1*mut1.EXE-gamma2*mut2.EXE)+
        parm2*(mu.EXE+gamma3*Xt1+gamma4*Xt2-gamma3*mut1.EXE-gamma4*mut2.EXE)
      
      comp2<-parm1^2*(coef.rd[1,1]-eta1)+parm2*parm1*(2*coef.rd[1,1]+2*(SE[i]-RE[jj.2.tie])*coef.rd[1,2]-eta2-eta3)+
      parm2^2*(coef.rd[1,1]+(SE[i]-RE[jj.2.tie])^2*coef.rd[2,2]+2*(SE[i]-RE[jj.2.tie])*coef.rd[1,2]-eta4)
      
    
    ss=ss-(ind.tie-1)/di*sum(weights[jj.2.tie]*exp(as.numeric(comp1)+0.5*comp2-comp12i+parm3*(SE[i]-RE[jj.2.tie]-ZE[i]*TE[i])
    +parm0*(ZE[jj.2.tie]-ZE[i])-t(as.numeric(WEc[i,])-matrix(as.numeric(t(WEc[jj.2.tie,])),dc,length(jj.2.tie)))%*%parm4) )  
    rm(Xt1,Xt2)
    }
    #####  ###### Measure of X(t1)*, X(t2)*, and X(t3)*
    
     if (length(jj.3)>0){
     Xt1=datwide.delta0[match(ID[jj.3],datwide.delta0$id),Vname[1]]
     Xt2=datwide.delta0[match(ID[jj.3],datwide.delta0$id),Vname[2]]
     Xt3=datwide.delta0[match(ID[jj.3],datwide.delta0$id),Vname[3]]
     
       ### t1, t2, t3 vector
     t1=datwide.delta0[match(ID[jj.3],datwide.delta0$id),Tname[1]]
     t2=datwide.delta0[match(ID[jj.3],datwide.delta0$id),Tname[2]]
     t3=datwide.delta0[match(ID[jj.3],datwide.delta0$id),Tname[3]]
     
   ####  t1,t2,t3 vector

     COVV.c1=coef.rd[1,1]+t1*coef.rd[1,2]
     COVV.c2=coef.rd[1,1]+t2*coef.rd[1,2]
     COVV.c3=coef.rd[1,1]+t3*coef.rd[1,2]
     COVV.c4=coef.rd[1,1]+(t1+SE[i]-RE[jj.3])*coef.rd[1,2]+t1*(SE[i]-RE[jj.3])*coef.rd[2,2]
     COVV.c5=coef.rd[1,1]+(t2+SE[i]-RE[jj.3])*coef.rd[1,2]+t2*(SE[i]-RE[jj.3])*coef.rd[2,2]
     COVV.c6=coef.rd[1,1]+(t3+SE[i]-RE[jj.3])*coef.rd[1,2]+t3*(SE[i]-RE[jj.3])*coef.rd[2,2]
     
     ##########
    
     SS.c.11<-coef.rd[1,1]+2*t1*coef.rd[1,2]+t1^2*coef.rd[2,2]+sigma.e^2
     SS.c.22<-coef.rd[1,1]+2*t2*coef.rd[1,2]+t2^2*coef.rd[2,2]+sigma.e^2
     SS.c.33<-coef.rd[1,1]+2*t3*coef.rd[1,2]+t3^2*coef.rd[2,2]+sigma.e^2
     SS.c.12=coef.rd[1,1]+(t1+t2)*coef.rd[1,2]+t1*t2*coef.rd[2,2]
     SS.c.13=coef.rd[1,1]+(t1+t3)*coef.rd[1,2]+t1*t3*coef.rd[2,2]
     SS.c.23=coef.rd[1,1]+(t2+t3)*coef.rd[1,2]+t2*t3*coef.rd[2,2]
     
     SS.c<-array(NA,dim=c(3,3,length(jj.3)))
     SS.c[1,1,]<-SS.c.11
     SS.c[2,2,]<-SS.c.22
     SS.c[3,3,]<-SS.c.33
     SS.c[1,2,]<-SS.c[2,1,]<-SS.c.12
     SS.c[1,3,]<-SS.c[3,1,]<-SS.c.13
     SS.c[2,3,]<-SS.c[3,2,]<-SS.c.23
     
    SS.c<-apply(SS.c,3,solve)
    

     #### SS.c is 9xlength(jj.3), elements in each column corresponds to [1,1],[2,1],[3,1],[1,2],[2,2],[3,2],[1,3],[2,3],[3,3] of a matrix

     delta1=COVV.c1*SS.c[1,]+COVV.c2*SS.c[2,]+COVV.c3*SS.c[3,]
     delta2=COVV.c1*SS.c[4,]+COVV.c2*SS.c[5,]+COVV.c3*SS.c[6,]
     delta3=COVV.c1*SS.c[7,]+COVV.c2*SS.c[8,]+COVV.c3*SS.c[9,]
     delta4=COVV.c4*SS.c[1,]+COVV.c5*SS.c[2,]+COVV.c6*SS.c[3,]
     delta5=COVV.c4*SS.c[4,]+COVV.c5*SS.c[5,]+COVV.c6*SS.c[6,]
     delta6=COVV.c4*SS.c[7,]+COVV.c5*SS.c[8,]+COVV.c6*SS.c[9,]

     xi1=delta1*COVV.c1+delta2*COVV.c2+delta3*COVV.c3
     xi2=delta1*COVV.c4+delta2*COVV.c5+delta3*COVV.c6
     xi3=delta4*COVV.c1+delta5*COVV.c2+delta6*COVV.c3
     xi4=delta4*COVV.c4+delta5*COVV.c5+delta6*COVV.c6
      
      mu0.EXE<-cbind(1,0,matrix(as.numeric(WE[jj.3,]),length(jj.3),dw))%*%coef.fix
       mu.EXE<-cbind(1,SE[i]-RE[jj.3],matrix(as.numeric(WE[jj.3,]),length(jj.3),dw))%*%coef.fix
       mut1.EXE<-cbind(1,t1,matrix(as.numeric(WE[jj.3,]),length(jj.3),dw))%*%coef.fix
       mut2.EXE<-cbind(1,t2,matrix(as.numeric(WE[jj.3,]),length(jj.3),dw))%*%coef.fix
       mut3.EXE<-cbind(1,t3,matrix(as.numeric(WE[jj.3,]),length(jj.3),dw))%*%coef.fix
       
      
      comp1<-parm1*(mu0.EXE+delta1*Xt1+delta2*Xt2+delta3*Xt3-delta1*mut1.EXE-delta2*mut2.EXE-delta3*mut3.EXE)+
        parm2*(mu.EXE+delta4*Xt1+delta5*Xt2+delta6*Xt3-delta4*mut1.EXE-delta5*mut2.EXE-delta6*mut3.EXE)
      
      comp2<-parm1^2*(coef.rd[1,1]-xi1)+parm2*parm1*(2*coef.rd[1,1]+2*(SE[i]-RE[jj.3])*coef.rd[1,2]-xi2-xi3)+
      parm2^2*(coef.rd[1,1]+(SE[i]-RE[jj.3])^2*coef.rd[2,2]+2*(SE[i]-RE[jj.3])*coef.rd[1,2]-xi4)
      
    
    ss=ss+sum(weights[jj.3]*exp(as.numeric(comp1)+0.5*comp2-comp12i+parm3*(SE[i]-RE[jj.3]-ZE[i]*TE[i])
    +parm0*(ZE[jj.3]-ZE[i])-t(as.numeric(WEc[i,])-matrix(as.numeric(t(WEc[jj.3,])),dc,length(jj.3)))%*%parm4) )  
   
    rm(Xt1,Xt2,Xt3)
    }

    #### Minus ties
 
     if (length(jj.3.tie)>0){
     Xt1=datwide.delta0[match(ID[jj.3.tie],datwide.delta0$id),Vname[1]]
     Xt2=datwide.delta0[match(ID[jj.3.tie],datwide.delta0$id),Vname[2]]
     Xt3=datwide.delta0[match(ID[jj.3.tie],datwide.delta0$id),Vname[3]]
     
      ### t1, t2, t3 vector
     t1=datwide.delta0[match(ID[jj.3.tie],datwide.delta0$id),Tname[1]]
     t2=datwide.delta0[match(ID[jj.3.tie],datwide.delta0$id),Tname[2]]
     t3=datwide.delta0[match(ID[jj.3.tie],datwide.delta0$id),Tname[3]]
 
     COVV.c1=coef.rd[1,1]+t1*coef.rd[1,2]
     COVV.c2=coef.rd[1,1]+t2*coef.rd[1,2]
     COVV.c3=coef.rd[1,1]+t3*coef.rd[1,2]
     COVV.c4=coef.rd[1,1]+(t1+SE[i]-RE[jj.3.tie])*coef.rd[1,2]+t1*(SE[i]-RE[jj.3.tie])*coef.rd[2,2]
     COVV.c5=coef.rd[1,1]+(t2+SE[i]-RE[jj.3.tie])*coef.rd[1,2]+t2*(SE[i]-RE[jj.3.tie])*coef.rd[2,2]
     COVV.c6=coef.rd[1,1]+(t3+SE[i]-RE[jj.3.tie])*coef.rd[1,2]+t3*(SE[i]-RE[jj.3.tie])*coef.rd[2,2]
     
     ##########
    
     SS.c.11<-coef.rd[1,1]+2*t1*coef.rd[1,2]+t1^2*coef.rd[2,2]+sigma.e^2
     SS.c.22<-coef.rd[1,1]+2*t2*coef.rd[1,2]+t2^2*coef.rd[2,2]+sigma.e^2
     SS.c.33<-coef.rd[1,1]+2*t3*coef.rd[1,2]+t3^2*coef.rd[2,2]+sigma.e^2
     SS.c.12=coef.rd[1,1]+(t1+t2)*coef.rd[1,2]+t1*t2*coef.rd[2,2]
     SS.c.13=coef.rd[1,1]+(t1+t3)*coef.rd[1,2]+t1*t3*coef.rd[2,2]
     SS.c.23=coef.rd[1,1]+(t2+t3)*coef.rd[1,2]+t2*t3*coef.rd[2,2]
     
     SS.c<-array(NA,dim=c(3,3,length(jj.3.tie)))
     SS.c[1,1,]<-SS.c.11
     SS.c[2,2,]<-SS.c.22
     SS.c[3,3,]<-SS.c.33
     SS.c[1,2,]<-SS.c[2,1,]<-SS.c.12
     SS.c[1,3,]<-SS.c[3,1,]<-SS.c.13
     SS.c[2,3,]<-SS.c[3,2,]<-SS.c.23
         
    SS.c<-apply(SS.c,3,solve)
    

     #### SS.c is 9xlength(jj.3), elements in each column corresponds to [1,1],[2,1],[3,1],[1,2],[2,2],[3,2],[1,3],[2,3],[3,3] of a matrix

     delta1=COVV.c1*SS.c[1,]+COVV.c2*SS.c[2,]+COVV.c3*SS.c[3,]
     delta2=COVV.c1*SS.c[4,]+COVV.c2*SS.c[5,]+COVV.c3*SS.c[6,]
     delta3=COVV.c1*SS.c[7,]+COVV.c2*SS.c[8,]+COVV.c3*SS.c[9,]
     delta4=COVV.c4*SS.c[1,]+COVV.c5*SS.c[2,]+COVV.c6*SS.c[3,]
     delta5=COVV.c4*SS.c[4,]+COVV.c5*SS.c[5,]+COVV.c6*SS.c[6,]
     delta6=COVV.c4*SS.c[7,]+COVV.c5*SS.c[8,]+COVV.c6*SS.c[9,]

     
     xi1=delta1*COVV.c1+delta2*COVV.c2+delta3*COVV.c3
     xi2=delta1*COVV.c4+delta2*COVV.c5+delta3*COVV.c6
     xi3=delta4*COVV.c1+delta5*COVV.c2+delta6*COVV.c3
     xi4=delta4*COVV.c4+delta5*COVV.c5+delta6*COVV.c6
      
      mu0.EXE<-cbind(1,0,matrix(as.numeric(WE[jj.3.tie,]),length(jj.3.tie),dw))%*%coef.fix
       mu.EXE<-cbind(1,SE[i]-RE[jj.3.tie],matrix(as.numeric(WE[jj.3.tie,]),length(jj.3.tie),dw))%*%coef.fix
       mut1.EXE<-cbind(1,t1,matrix(as.numeric(WE[jj.3.tie,]),length(jj.3.tie),dw))%*%coef.fix
       mut2.EXE<-cbind(1,t2,matrix(as.numeric(WE[jj.3.tie,]),length(jj.3.tie),dw))%*%coef.fix
       mut3.EXE<-cbind(1,t3,matrix(as.numeric(WE[jj.3.tie,]),length(jj.3.tie),dw))%*%coef.fix
            
      comp1<-parm1*(mu0.EXE+delta1*Xt1+delta2*Xt2+delta3*Xt3-delta1*mut1.EXE-delta2*mut2.EXE-delta3*mut3.EXE)+
        parm2*(mu.EXE+delta4*Xt1+delta5*Xt2+delta6*Xt3-delta4*mut1.EXE-delta5*mut2.EXE-delta6*mut3.EXE)
      
      comp2<-parm1^2*(coef.rd[1,1]-xi1)+parm2*parm1*(2*coef.rd[1,1]+2*(SE[i]-RE[jj.3.tie])*coef.rd[1,2]-xi2-xi3)+
      parm2^2*(coef.rd[1,1]+(SE[i]-RE[jj.3.tie])^2*coef.rd[2,2]+2*(SE[i]-RE[jj.3.tie])*coef.rd[1,2]-xi4)
      
    ss=ss-(ind.tie-1)/di*sum(weights[jj.3.tie]*exp(as.numeric(comp1)+0.5*comp2-comp12i+parm3*(SE[i]-RE[jj.3.tie]-ZE[i]*TE[i])
    +parm0*(ZE[jj.3.tie]-ZE[i])-t(as.numeric(WEc[i,])-matrix(as.numeric(t(WEc[jj.3.tie,])),dc,length(jj.3.tie)))%*%parm4) )  
   
    rm(Xt1,Xt2,Xt3)
    }

    ss.new=ss
    ### placebo
   
   ### need to use jj.Z0 and jj.Z0.tie when there are confounders in the risk model
   if (length(jj.Z0)>0)
   ss=ss+sum(weights[jj.Z0]*exp(-comp12i-parm3*(ZE[i]*TE[i])-parm0*(ZE[i])
   -t(as.numeric(WEc[i,])-matrix(as.numeric(t(WEc[jj.Z0,])),dc,length(jj.Z0)))%*%parm4))
   
   if (length(jj.Z0.tie)>0)
   ss=ss-(ind.tie-1)/di*sum(weights[jj.Z0.tie]*exp(-comp12i-parm3*(ZE[i]*TE[i])-parm0*(ZE[i])-t(as.numeric(WEc[i,])-matrix(as.numeric(t(WEc[jj.Z0.tie,])),dc,length(jj.Z0.tie)))%*%parm4))
   
    out=out+weights.bar*log(weights[i]+ss)+(weights.bar-weights[i])*(comp12i+parm3*(SE[i]-RE[i])*ZE[i]+parm0*(ZE[i])
    +matrix(WEc[i,],length(i),dc)%*%parm4)
    
    outi=c(outi,ss)
 
    }
 return(out)
 }

########
