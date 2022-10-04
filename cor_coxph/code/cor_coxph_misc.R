#with(subset(dat.mock, Trt==1 & ph1), table(EventIndPrimary, useNA="ifany"))
#with(subset(dat.mock, Trt==1 & ph1), summary(EventTimePrimary, useNA="ifany"))

#par(mfrow=c(2,2))
#with(subset(dat.mock, Trt==1 & ph2 & Region==0), hist(EventTimePrimary[EventIndPrimary==0], xlim=c(0,100), main="NA"))
#with(subset(dat.mock, Trt==1 & ph2 & Region==1), hist(EventTimePrimary[EventIndPrimary==0], xlim=c(0,100), main="LA"))
#with(subset(dat.mock, Trt==1 & ph2 & Region==0), hist(EventTimePrimary[EventIndPrimary==1], xlim=c(0,100), main=""))
#with(subset(dat.mock, Trt==1 & ph2 & Region==1), hist(EventTimePrimary[EventIndPrimary==1], xlim=c(0,100), main=""))

#with(subset(dat.mock, Trt==1 & ph1), summary(EventTimePrimary[EventIndPrimary==0], useNA="ifany"))
#with(subset(dat.mock, Trt==1 & ph1), summary(EventTimePrimary[EventIndPrimary==1], useNA="ifany"))

#with(subset(dat.mock, ph1==1), table(EventIndPrimary, Trt))
#with(subset(dat.mock, ph2==1), table(EventIndPrimary, Trt))
#with(subset(dat.mock, Trt==1 & ph1), table(Wstratum, ph2))
#with(subset(dat.mock, Trt==1 & ph1), table(tps.stratum, ph2))
#with(subset(dat.mock, Trt==1), table(Wstratum, ph2))
#with(subset(dat.mock, Trt==1), table(Wstratum, wt.D29))
#with(subset(dat.mock, ph1==1), table(Delta.D210, Trt))
#with(subset(dat.mock, ph2==1), corplot(Day210ELCZ,Day210ADCPgp140C97ZAfib))
#with(subset(dat.mock, ph2==1), corr(cbind(Day210ELCZ,Day210ADCPgp140C97ZAfib), w = wt))
    
#with(subset(dat.mock, Trt==1 & ph2 & EventIndPrimary), mean(Day29pseudoneutid50>log10(3)) )
#with(subset(dat.mock, Trt==1 & ph2 & !EventIndPrimary), mean(Day29pseudoneutid50>log10(3)) )



# Average follow-up of vaccine recipients starting at tpeaklag days post visit
write(round(mean(subset(dat.mock, Trt==1 & ph1, EventTimePrimary, drop=T), na.rm=T)-tpeaklag), file=paste0(save.results.to, "avg_followup_"%.%study_name))


## experimenting with multitesting
## based on moderna_mock
## make their correlation 0.98
#dat.vac.seroneg$Day57bindRBD = dat.vac.seroneg$Day57bindSpike + rnorm(nrow(dat.vac.seroneg), sd=.1)
#dat.vac.seroneg$Day57pseudoneutid50 = dat.vac.seroneg$Day57pseudoneutid80 + rnorm(nrow(dat.vac.seroneg), sd=.1)
## switch 50 and 80
#tmp=dat.vac.seroneg$Day57pseudoneutid50
#dat.vac.seroneg$Day57pseudoneutid50=dat.vac.seroneg$Day57pseudoneutid80
#dat.vac.seroneg$Day57pseudoneutid80=tmp



## commented out now because it won't run with the real processed moderna dataset, which did not include baseline pos
## Decile immunogenicity table
## do this before uloq censoring
#if (config$is_ows_trial) {
#    
#    if ("Day"%.%tpeak%.%"pseudoneutid50" %in% names(dat.mock)) {    
#        res=lapply (0:1, function(ii) {
#            dat.immuno.seroneg=subset(dat.mock, Trt==1 & Bserostatus==ii & ph2.immuno)    
#            ww=sort(unique(dat.immuno.seroneg$demo.stratum))
#            myprint(ww)
#            stopifnot(min(ww)==1)
#            stopifnot(max(ww)==length(ww))
#            names(ww)=demo.stratum.labels
#            mysapply (c(All=0,ww), function(w) { 
#                if(verbose) myprint(w)
#                dat.tmp= if (w==0) dat.immuno.seroneg else subset(dat.immuno.seroneg, demo.stratum==w)
#                10**wtd.quantile(dat.tmp[["Day"%.%tpeak%.%"pseudoneutid50"]], weights = dat.tmp$wt.subcohort, probs = 0:10/10)
#            })
#        })
#        tab=rbind(res[[1]], res[[2]])
#        colnames(tab)[1]="min"
#        colnames(tab)[ncol(tab)]="max"
#        tab
#        mytex(tab, file.name="cID50_deciles_"%.%study_name, align="r", include.colnames = T, save2input.only=T, input.foldername=save.results.to, digits=0,
#            add.to.row=list(list(0,nrow(res[[1]])), # insert at the beginning of table, and at the end of, say, the first table
#                c("       \n \\multicolumn{12}{l}{Baseline negative} \\\\ \n",
#                  "\\hline\n \\multicolumn{12}{l}{Baseline positive} \\\\ \n"
#                 )
#            )    
#        )
#    }
#        
#}


#Median and IQR and range of days from dose 1 to Day 29 visit, and from dose 1 to Day 57 visit (requested by Lindsey B).  
#subsetting by (a) the whole immunogenicity subcohort, (2) non-cases in the immunogenicity subcohort, (3) intercurrent cases, (4) primary cases.
if (study_name=="COVE" | study_name=="MockCOVE") {
        
    # D1 to tpeak
    tab=sapply(1:4, function (i) {
        idx=with(dat.mock, {
            tmp = (if(i==1) ph1.immuno else if(i==2) (ph1.immuno & EventIndPrimary) else if(i==3) ph1.intercurrent.cases else if(i==4) (ph1.D57 & EventIndPrimaryD57)) & Trt==1 & Bserostatus==0
            tmp [!is.na(tmp)]
        })
        res=c(quantile(dat.mock[idx, "NumberdaysD1toD"%.%tpeak], c(0, 1:3/4, 1), na.rm=T))
        res
    })
    tab=t(tab)
    rownames(tab)=c("(a)","(b)","(c)","(d)")
    colnames(tab)=c("min", "1st quartile", "median", "3d quartile", "max")
    # the file name is misnomer for D57
    mytex(tab, file.name="dose1_interval_summary_"%.%study_name, align="c", include.colnames = T, save2input.only=T, input.foldername=save.results.to, digits=0)    
}

    

if ((study_name=="COVE" | study_name=="MockCOVE") & tpeak=="57") {

    # Number of breakthrough vaccine cases with Day 57 ID80 > 660 IU
    res=nrow(subset(dat.mock, Trt==1 & Bserostatus==0 & ph1 & EventIndPrimary & Day57pseudoneutid80>log10(660)))
    write(res, file=paste0(save.results.to, "num_vacc_cases_highid80_"%.%study_name))
    
    dat.mock$NumberdaysD29toD57=dat.mock$NumberdaysD1toD57-dat.mock$NumberdaysD1toD29
    
    # D29 to D57
    tab=sapply(1:4, function (i) {
        idx=with(dat.mock, {
            tmp = (if(i==1) ph1.immuno else if(i==2) (ph1.immuno & EventIndPrimary) else if(i==3) ph1.intercurrent.cases else if(i==4) (ph1.D57 & EventIndPrimaryD57)) & Trt==1 & Bserostatus==0
            tmp [!is.na(tmp)]
        })
        res=c(quantile(dat.mock[idx, "NumberdaysD29toD57"], c(0, 1:3/4, 1), na.rm=T))
        res
    })
    tab=t(tab)
    rownames(tab)=c("(a)","(b)","(c)","(d)")
    colnames(tab)=c("min", "1st quartile", "median", "3d quartile", "max")
    mytex(tab, file.name="dose1_dose2_interval_summary_"%.%study_name, align="c", include.colnames = T, save2input.only=T, input.foldername=save.results.to, digits=0)

    
    # For Day 29 vaccine breakthrough cases: 
    # (1) for Intercurrent cases number of days from Day 1 to Day 29, 
    # (2) for Intercurrent cases number of days from Day 29 to COVID endpoint, 
    # (3) for Post Day 57 cases number of days from Day 1 to Day 29, 
    # (4) for Post Day 57 cases number of days from Day 29 to Day 57, 
    # (5) for Post Day 57 cases number of days from Day 57 to COVID endpoint (e.g. with 2 panels on a first row and 3 panels on a bottom row).
    mypdf (mfrow=c(1,2), file=paste0(save.results.to, "barplot_visit_intervals_intercurrentcases"))    
        tmp.1=table(subset(dat.mock, ph1.intercurrent.cases & Trt==1 & Bserostatus==0, NumberdaysD1toD29, drop=T))
        tmp.2=table(subset(dat.mock, ph1.intercurrent.cases & Trt==1 & Bserostatus==0, EventTimePrimaryD29, drop=T))
        
        tmp=cbinduneven(list(tmp.1, tmp.2))
        tmp=tmp[order(as.numeric(rownames(tmp))),]
        tmp.1=tmp[,1]; names(tmp.1)=rownames(tmp)
        tmp.2=tmp[,2]; names(tmp.2)=rownames(tmp)
        
        barplot(tmp.1, main="D1 to D29",  xlab="Days")
        barplot(tmp.2, main="D29 to COVID", xlab="Days")
    dev.off()  
    
    
    mypdf (mfrow=c(1,3), file=paste0(save.results.to, "barplot_visit_intervals_D57cases"))    
        tmp.1=table(subset(dat.mock, ph1.D57 & EventIndPrimaryD57 & Trt==1 & Bserostatus==0, NumberdaysD1toD29, drop=T))
        tmp.2=table(subset(dat.mock, ph1.D57 & EventIndPrimaryD57 & Trt==1 & Bserostatus==0, NumberdaysD29toD57, drop=T))
        tmp.3=table(subset(dat.mock, ph1.D57 & EventIndPrimaryD57 & Trt==1 & Bserostatus==0, EventTimePrimaryD57, drop=T))

        tmp=cbinduneven(list(tmp.1, tmp.2, tmp.3))
        tmp=tmp[order(as.numeric(rownames(tmp))),]
        tmp.1=tmp[,1]; names(tmp.1)=rownames(tmp)
        tmp.2=tmp[,2]; names(tmp.2)=rownames(tmp)
        tmp.3=tmp[,3]; names(tmp.3)=rownames(tmp)
        
        barplot(tmp.1, main="D1 to D29",  xlab="Days")
        barplot(tmp.2, main="D29 to D57", xlab="Days")
        barplot(tmp.3, main="D57 to COVID", xlab="Days")
    dev.off()  




} 


if (startsWith(attr(config, "config"), "janssen_pooled_")) {
    # the distribution of the number of days from D29 until COVID-19 endpoint occurrence or until right-censoring
    # stratified by COVID-19 endpoint or right-censoring and by geographic region.  Based on this figure,        

    tmp.1=table(subset(dat.mock, ph1.D29start1 & EventIndPrimary==1 & Region==0, EventTimePrimary, drop=T))
    tmp.2=table(subset(dat.mock, ph1.D29start1 & EventIndPrimary==1 & Region==1, EventTimePrimary, drop=T))
    tmp.3=table(subset(dat.mock, ph1.D29start1 & EventIndPrimary==1 & Region==2, EventTimePrimary, drop=T))
    tmp.4=table(subset(dat.mock, ph1.D29start1 & EventIndPrimary==0 & Region==0, EventTimePrimary, drop=T))
    tmp.5=table(subset(dat.mock, ph1.D29start1 & EventIndPrimary==0 & Region==1, EventTimePrimary, drop=T))
    tmp.6=table(subset(dat.mock, ph1.D29start1 & EventIndPrimary==0 & Region==2, EventTimePrimary, drop=T))
    
    minmax=range(c(as.numeric(names(tmp.1)),as.numeric(names(tmp.2)),as.numeric(names(tmp.3)),as.numeric(names(tmp.4)),as.numeric(names(tmp.5)),as.numeric(names(tmp.6))))
    tmp.7=rep(0, minmax[2]-minmax[1]+1)    
    names(tmp.7)=minmax[1]:minmax[2]
    tmp.7=as.table(tmp.7)
        
    tmp=cbinduneven(list(tmp.1, tmp.2, tmp.3, tmp.4, tmp.5, tmp.6, tmp.7))
    tmp=tmp[order(as.numeric(rownames(tmp))),]
    tmp.1=tmp[,1]; names(tmp.1)=rownames(tmp)
    tmp.2=tmp[,2]; names(tmp.2)=rownames(tmp)
    tmp.3=tmp[,3]; names(tmp.3)=rownames(tmp)
    tmp.4=tmp[,4]; names(tmp.4)=rownames(tmp)
    tmp.5=tmp[,5]; names(tmp.5)=rownames(tmp)
    tmp.6=tmp[,6]; names(tmp.6)=rownames(tmp)
    
    mypdf (mfrow=c(2,3), file=paste0(save.results.to, "revision_fig_1"))     
        barplot(tmp.1, xlab="Days", main="Cases, United States")
        barplot(tmp.2, xlab="Days", main="Cases, Latin America")
        barplot(tmp.3, xlab="Days", main="Cases, South Africa") 
        barplot(tmp.4, xlab="Days", main="Non-cases, United States")
        barplot(tmp.5, xlab="Days", main="Non-cases, Latin America")
        barplot(tmp.6, xlab="Days", main="Non-cases, South Africa") 
    dev.off()  


    # last case in vaccine arm
    tmp=subset(dat.vac.seroneg, ph2==1 & EventIndPrimary==1, select=c(EventIndPrimary, EventTimePrimary))
    tmp[order(tmp[,2]),] # last 144
    
    tmp=subset(dat.vac.seroneg, ph1==1 & EventIndPrimary==1 & Region==0, select=c(EventIndPrimary, EventTimePrimary))
    tmp[order(tmp[,2]),] # last 47
    
    tmp=subset(dat.vac.seroneg, ph2==1 & EventIndPrimary==1 & Region==1, select=c(EventIndPrimary, EventTimePrimary))
    tmp[order(tmp[,2]),] # last 101
    
    tmp=subset(dat.vac.seroneg, ph2==1 & EventIndPrimary==1 & Region==2, select=c(EventIndPrimary, EventTimePrimary))
    tmp[order(tmp[,2]),] # last 144
    
    
    # last case in vaccine arm
    tmp=subset(dat.vac.seroneg, EventIndPrimary==1, select=c(EventIndPrimary, EventTimePrimary))
    tmp[order(tmp[,2]),] # last 195
    
    tmp=subset(dat.vac.seroneg, EventIndPrimary==1 & Region==0, select=c(EventIndPrimary, EventTimePrimary))
    tmp[order(tmp[,2]),] # last 111
    
    tmp=subset(dat.vac.seroneg, EventIndPrimary==1 & Region==1, select=c(EventIndPrimary, EventTimePrimary))
    tmp[order(tmp[,2]),] # last 191
    
    tmp=subset(dat.vac.seroneg, EventIndPrimary==1 & Region==2, select=c(EventIndPrimary, EventTimePrimary))
    tmp[order(tmp[,2]),] # last 195
    
    # with 15 subjects at risk
    tmp=subset(dat.vac.seroneg, ph2==1 & SubcohortInd==1, select=c(EventIndPrimary, EventTimePrimary))
    tail(tmp[order(tmp[,2]),], 15) # last 220, first 190
    
    tmp=subset(dat.vac.seroneg, ph2==1 & SubcohortInd==1 & Region==0, select=c(EventIndPrimary, EventTimePrimary))
    tail(tmp[order(tmp[,2]),], 15) # last 220, first 158
    
    tmp=subset(dat.vac.seroneg, ph2==1 & SubcohortInd==1 & Region==1, select=c(EventIndPrimary, EventTimePrimary))
    tail(tmp[order(tmp[,2]),], 15) # last 211, first 182
    
    tmp=subset(dat.vac.seroneg, ph2==1 & SubcohortInd==1 & Region==2, select=c(EventIndPrimary, EventTimePrimary))
    tail(tmp[order(tmp[,2]),], 15) # last 207, first 183
    
    
    # placebo, note that the following days are all shorter than vaccine arm
    tmp=subset(dat.pla.seroneg,          SubcohortInd==1,             select=c(EventIndPrimary, EventTimePrimary))
    tail(tmp[order(tmp[,2]),], 15) # last 198, first 157
    
    tmp=subset(dat.pla.seroneg,          SubcohortInd==1 & Region==0, select=c(EventIndPrimary, EventTimePrimary))
    tail(tmp[order(tmp[,2]),], 15) # last 155, first 90
    
    tmp=subset(dat.pla.seroneg,          SubcohortInd==1 & Region==1, select=c(EventIndPrimary, EventTimePrimary))
    tail(tmp[order(tmp[,2]),], 15) # last 198, first 84
    
    tmp=subset(dat.pla.seroneg,          SubcohortInd==1 & Region==2, select=c(EventIndPrimary, EventTimePrimary))
    tail(tmp[order(tmp[,2]),], 15) # last 182, first 99
}

if (tpeak=="57") {
        
    # barplots for number of days from Day 1 to Day 29, and number of days from Day 29 to Day 57 in the immunogenicity subcohort
    mypdf (mfrow=c(1,2), file=paste0(save.results.to, "barplot_visit_intervals_immuno"))    
        tmp.1=table(subset(dat.mock, ph1.immuno & Trt==1 & Bserostatus==0, NumberdaysD1toD29,  drop=T))
        tmp.2=table(subset(dat.mock, ph1.immuno & Trt==1 & Bserostatus==0, NumberdaysD29toD57, drop=T))
        tmp=cbinduneven(list(tmp.1, tmp.2))
        tmp=tmp[order(rownames(tmp)),]
        tmp.1=tmp[,1]; names(tmp.1)=rownames(tmp)
        tmp.2=tmp[,2]; names(tmp.2)=rownames(tmp)
        barplot(tmp.1, main="D1 to D29",  xlab="Days"); title(line=3, main="Immunogenicity Subcohort")# , cex.names=.7
        barplot(tmp.2, main="D29 to D57", xlab="Days"); title(line=3, main="Immunogenicity Subcohort")
    dev.off()  
    
    #
    res=round(quantile(subset(dat.mock, ph1.immuno & Trt==1 & Bserostatus==0, NumberdaysD1toD29, drop=T), c(1/4,1/2,3/4), na.rm=T)); res
    write(paste0(res[2], " (", res[1], "-", res[3], ")"), file=paste0(save.results.to, "quartiles_visit_intervals_immuno_1to29"))
    #
    res=round(quantile(subset(dat.mock, ph1.immuno & Trt==1 & Bserostatus==0, NumberdaysD29toD57, drop=T), c(1/4,1/2,3/4), na.rm=T)); res
    write(paste0(res[2], " (", res[1], "-", res[3], ")"), file=paste0(save.results.to, "quartiles_visit_intervals_immuno_29to57"))
        
} else if (tpeak=="29") {
    
    # barplots for number of days from Day 1 to Day 29 in the immunogenicity subcohort
    mypdf (mfrow=c(1,1), file=paste0(save.results.to, "barplot_visit_intervals_immuno"))    
        tmp.1=table(subset(dat.mock, ph1.immuno & Trt==1 & Bserostatus==0, NumberdaysD1toD29,  drop=T))
        barplot(tmp.1, main="D1 to D29",  xlab="Days")# , cex.names=.7
    dev.off()  
    
    res=round(quantile(subset(dat.mock, ph1.immuno & Trt==1 & Bserostatus==0, NumberdaysD1toD29, drop=T), c(1/4,1/2,3/4), na.rm=T)); res
    write(paste0(res[2], " (", res[1], "-", res[3], ")"), file=paste0(save.results.to, "quartiles_visit_intervals_immuno_1to29"))
    
}


# understand different follow-up periods in different regions
if (attr(config,"config")=="janssen_pooled_EUA") {   
    form = Surv(EventTimePrimary, EventIndPrimary) ~ Region
    
    par(mfrow=c(1,2))
    
    km <- survfit(form, subset(dat.mock, Trt==1 & ph1))
    plot(km, ylim=c(.96,1), col=1:3, main="ph1", xlim=c(0,100))
    
    km <- survfit(form, subset(dat.mock, Trt==1 & ph2))
    plot(km, ylim=c(0.6,1), col=1:3, main="ph2", xlim=c(0,100))
    
}
