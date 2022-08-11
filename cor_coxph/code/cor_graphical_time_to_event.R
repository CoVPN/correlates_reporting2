#Sys.setenv(TRIAL = "vat08m_nonnaive"); COR="D22D43"; Sys.setenv(VERBOSE = 1)
renv::activate(project = here::here(".."))    
    # There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
    if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))
source(here::here("..", "_common.R"))
#-----------------------------------------------


library(kyotil) # p.adj.perm, getFormattedSummary


# path for figures and tables etc
save.results.to = here::here("output")
if (!dir.exists(save.results.to))  dir.create(save.results.to)
save.results.to = paste0(here::here("output"), "/", attr(config,"config"));
if (!dir.exists(save.results.to))  dir.create(save.results.to)
save.results.to = paste0(save.results.to, "/", COR,"/");
if (!dir.exists(save.results.to))  dir.create(save.results.to)
print(paste0("save.results.to equals ", save.results.to))
    


mypdf (mfrow=c(1,3), file=paste0(save.results.to, "barplot_mixed"))     
    tmp.1=table(subset(dat.mock, ph1.intercurrent.cases       & Trt==1, "EventTimePrimaryD"%.%tinterm, drop=T))
    tmp.2=table(subset(dat.mock, dat.mock[["ph1.D"%.%tpeak]] & dat.mock[["EventIndPrimaryD"%.%tpeak]] & Trt==1, "EventTimePrimaryD"%.%tinterm, drop=T))
    tmp.3=table(subset(dat.mock, dat.mock[["ph1.D"%.%tpeak]] & dat.mock[["EventIndPrimaryD"%.%tpeak]] & Trt==1, "EventTimePrimaryD"%.%tpeak  , drop=T))
    
    minmax=range(as.numeric(names(tmp.1)), as.numeric(names(tmp.2)), as.numeric(names(tmp.3)))
    tmp.4=rep(0, minmax[2]-minmax[1]+1)    
    names(tmp.4)=minmax[1]:minmax[2]
    tmp.5=as.table(tmp.4)
        
    tmp=cbinduneven(list(tmp.1, tmp.2, tmp.3, tmp.5))
    tmp=tmp[order(as.numeric(rownames(tmp))),]
    tmp.1=tmp[,1]; names(tmp.1)=rownames(tmp)
    tmp.2=tmp[,2]; names(tmp.2)=rownames(tmp)
    tmp.3=tmp[,3]; names(tmp.3)=rownames(tmp)
    
    if(all(is.na(tmp.1))) {
        empty.plot()
    } else {
        barplot(tmp.1, main="D"%.%tinterm%.%" to COVID", xlab="Days", yaxt="n", xaxt="n"); title(line=3, main="Intercurrent Cases"); axis(2, at=0:10); axis(1, at=seq(0,minmax[2],by=10)); 
    }
    barplot(tmp.2, main="D"%.%tinterm%.%" to COVID", xlab="Days", yaxt="n", xaxt="n"); title(line=3, main="Post Day "%.%tpeak%.%" Cases");  axis(2, at=0:10); axis(1, at=seq(0,minmax[2],by=10)); 
    barplot(tmp.3, main="D"%.%tpeak%.%"   to COVID", xlab="Days", yaxt="n", xaxt="n"); title(line=3, main="Post Day "%.%tpeak%.%" Cases");  axis(2, at=0:10); axis(1, at=seq(0,minmax[2],by=10)); 
dev.off()  
