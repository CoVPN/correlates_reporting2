################################################################################
# prevent19_stage2 Delta COVID correlates

library(survey)
library(kyotil)

dat_mapped=read.csv('/trials/covpn/p3004/analysis/mapping_immune_correlates/stage2/adata/COVID_Novavax_stage2_mapped_20240116.csv')
dat_proc=read.csv('/trials/covpn/p3004/analysis/correlates/Part_A_Blinded_Phase_Data/adata/prevent19_stage2_data_processed_20230928.csv')

assay_metadata=read.csv('~/correlates_reporting2/assay_metadata/prevent19_stage2_assay_metadata.csv')
assays=assay_metadata$assay
assays

times=c("Day35","C1","BD1","DD1")


##############################################
# data availability

# dataset contains vaccine arm only and seroneg only
mytable(dat_mapped$Trt)
mytable(dat_mapped$Bserostatus)

# who are these 126 ptids not measured a second time?
with(subset(dat_mapped, D35ancestralnAbBatch!=23), mytable(EventIndPrimaryD1))

# missingness across time points for ancestral id50
for (j in 1:2){
  # change the comments to focus on different subsets
  dat=subset(dat_mapped)
  # dat=subset(dat_mapped, D35ancestralnAbBatch==23) # remove batch 1
  # dat=subset(dat_mapped, VARTYPE=="Delta") # delta cases
  
  v=ifelse(j==1,"D614G","Delta")
  print(v)
  
  for (i in 1:2) {
    dat[['tmp'%.%i]]=ifelse(is.na(dat[[times[i]%.%'pseudoneutid50_'%.%v]]), 0, 1)
  }
  dat$tmp=with(dat,paste0(tmp1,tmp2))
  print(table(dat$tmp))
}


# correlation between ancestral, delta ID50s
dat=dat_mapped[,c(times%.%"pseudoneutid50_D614G",times%.%"pseudoneutid50_Delta")]
names(dat)=sub("pseudoneutid50","ID50",names(dat))
mypairs(dat)

dat=dat_mapped[,c(rbind(times%.%"pseudoneutid50_D614G",times%.%"pseudoneutid50_Delta"))]
names(dat)=sub("pseudoneutid50","ID50",names(dat))
mypairs(dat)

# cor highers at BD1 and DD1
# delta ID50 less sensitive, i.e. less pos rate, than ancestral ID50

dat=dat_mapped[,c(rbind(times%.%"pseudoneutid50_D614G",times%.%"pseudoneutid50_Delta"))]
names(dat)=sub("pseudoneutid50","ID50",names(dat))
myboxplot(dat)



##############################################
# batch effect
# There are 679 samples at D35 measured twice, once in Dec 21 (Test.Results.1) and once in July 23 (Test.Results.2)
# The Pearson correlation is 0.88. Two-sample test by wilcoxon paired test p value is 6.046e-16. 
# Applications that are sensitives to a shift of ~20% fold change in titers will be impacted

dat.2=read.csv('/trials/covpn/p3004/download_data/2019nCoV2-301_ nAB_D614G_27Oct2023 - 2.csv')

mytable(dat.2$PVC, dat.2$Test.Method)

# find duplicated measurement of the same samples
dat.2.anc=subset(dat.2, Test.Method=='SARS-COV-2 D614G NEUT' & PVC=='D35')
dat.2.anc=dat.2.anc[order(dat.2.anc$Patient.Number ),]

tmp = mytable(dat.2.anc$Patient.Number)
mytable(tmp)
tmp[tmp==3] # US097-0021
subset(dat.2.anc, Patient.Number=='US097-0021')

dat.2.anc=subset(dat.2.anc, Patient.Number %in% names(tmp[tmp==2]))
dat.2.anc$batch = ifelse(startsWith(dat.2.anc$Monogram.Accession,'23'),2,1)

subset(dat.2.anc, Patient.Number=='US152-0064')

dat.n=myreshapewide(Titer.Result~batch, dat.2.anc, 'Accession.Container')

clean=function(x) {
  x[x=='< 40']=20
  x[x=='ZNG40']=20
  as.integer((x))
}

dat.n$Titer.Result.1=clean(dat.n$Titer.Result.1)
dat.n$Titer.Result.2=clean(dat.n$Titer.Result.2)

dat.n$dif=log10(dat.n$Titer.Result.1)-log10(dat.n$Titer.Result.2)
summary(dat.n$dif)
summary(10**dat.n$dif)
myboxplot(list(dat.n$dif))

corplot(Titer.Result.1~Titer.Result.2,dat.n, log='xy', method='p')
# abline(h=(8282), lty=2)
# abline(v=(26859), lty=2)
abline(a=-log10(2),b=1, lty=2) # 2 fold difference
abline(a=log10(2),b=1, lty=2)

myboxplot(dat.n[,c('Titer.Result.1','Titer.Result.2')], log='y', test='w')

wilcox.test(dat.n$Titer.Result.1, dat.n$Titer.Result.2,paired=T)

nrow(dat.n)




################################################################################
# stage 1 - ancestral COVID correlates

library(kyotil)


#binding Ab
dat=read.csv("/trials/covpn/p3004/download_data/Binding Antibody/BARDA_SARS_MSD_2019nCoV2-301_Nexelis_6001_04Jan2022.csv")

with(dat, mytable(Test.Name)) # SARS-CoV-2 Spike Protein (Wild Type) 


x=subset(dat, Visit=="35.00 Day", Test.Result, drop=T)
x.1=ifelse(x=="<150.4", 75.2, as.numeric(x))

hist(log10(x.1), xlab="log10 (bAb Spike)")



#ID50
dat.1=read.csv("/trials/covpn/p3004/download_data/Neutralizing_Antibody/2019nCoV2-301_ nAB_D614G_09Mar2022.csv")
#dat.1=read.csv("/trials/covpn/p3004/download_data/Neutralizing_Antibody/2019nCoV2-301_ nAB_D614G_30Mar2022.csv")
with(dat.1, mytable(Visit.Description, PVC)) # just use PVC D35 and D0
with(dat.1, mytable(Test.Method, Specimen.Type, Unit.of.Results)) # SARS-COV-2 D614G NEUT, SERUM, 1/Dilution

with(dat.1, mytable(Reason.Not.Done.or.Comment.for.Result)) # some are "ZPABG"

# two that have records duplicated on D0 due to input errors
subset(dat.1, Patient.Number=="MX006-0187")
subset(dat.1, Patient.Number=="US179-0077")


tab=with(dat.1, mytable(Patient.Number, PVC))





###################################################################################################
# earlier stuff

dat=as.data.frame(read_excel("../Novavax/COVPN Sample Strata Freqs.xlsx"))

dat$racestrata[dat$country=="MX"]=""
dat$comorbid[dat$country=="MX"]=""

dat$country[dat$country=="MX"]="XM"
dat$N=dat$Frequency

dat=aggregate(N~agestrata+racestrata+comorbid+country+blserostat, dat, sum)

cbind(dat, n=c(rep("100+12",8), rep("50+6",2), rep("35+35",8), rep("13+13",2)))


# mapped data
dat.04=read.csv("/trials/covpn/p3004/analysis/mapping_immune_correlates/adata/COVID_Novavax_realdata_20220404.csv")
dat.18=read.csv("/trials/covpn/p3004/analysis/mapping_immune_correlates/adata/COVID_Novavax_realdata_20220418.csv")

dat.tmp.1=subset(dat.04, Subjectid %in% c("2019nCoV301-US049-0039","2019nCoV301-US159-0063","2019nCoV301-US159-0135","2019nCoV301-US162-0231"
                                          ,"2019nCoV301-US163-0193","2019nCoV301-US172-0511","2019nCoV301-US200-0046","2019nCoV301-US228-0067"
                                          ,"2019nCoV301-US228-0210","2019nCoV301-US239-0100","2019nCoV301-US241-0121","2019nCoV301-US242-0118"
))$Day35pseudoneutid50
dat.tmp.2=subset(dat.18, Subjectid %in% c("2019nCoV301-US049-0039","2019nCoV301-US159-0063","2019nCoV301-US159-0135","2019nCoV301-US162-0231"
                                          ,"2019nCoV301-US163-0193","2019nCoV301-US172-0511","2019nCoV301-US200-0046","2019nCoV301-US228-0067"
                                          ,"2019nCoV301-US228-0210","2019nCoV301-US239-0100","2019nCoV301-US241-0121","2019nCoV301-US242-0118"
))$Day35pseudoneutid50

with(subset(dat.1, Trt==1), mytable(!is.na(Day35bindSpike), !is.na(Day35pseudoneutid50), EventIndPrimaryD35))
#, , EventIndPrimaryD35 = 0
#
#       
#        FALSE  TRUE
#  FALSE 16594     4
#  TRUE     23   987
#
#, , EventIndPrimaryD35 = 1
#
#       
#        FALSE  TRUE
#  FALSE     0     0
#  TRUE      1    11
## This supports imputing ID50 by Spike binding


subset(dat.1, Trt==1 & EventIndPrimaryD35)


with(subset(dat.1, Trt==1), mytable(!is.na(BbindSpike), !is.na(Bpseudoneutid50), EventIndPrimaryD35))
#, , EventIndPrimaryD35 = 0
#
#       
#        FALSE  TRUE
#  FALSE 16591     2
#  TRUE     45   970
#
#, , EventIndPrimaryD35 = 1
#
#       
#        FALSE  TRUE
#  FALSE     0     0
#  TRUE      1    11





# analysis ready data
dat.2.pr.new=read.csv("/trials/covpn/p3004/analysis/correlates/Part_A_Blinded_Phase_Data/adata/prevent19_data_processed_with_riskscore.csv")
dat.2.pr.old=read.csv("/trials/covpn/p3004/analysis/correlates/Part_A_Blinded_Phase_Data/adata/archive/prevent19_data_processed_with_riskscore_2022-04-04 time 18.05.03.csv")

dat.tmp.1=subset(dat.2.pr.new, Ptid %in% c("2019nCoV301-US049-0039","2019nCoV301-US159-0063","2019nCoV301-US159-0135","2019nCoV301-US162-0231"
                                           ,"2019nCoV301-US163-0193","2019nCoV301-US172-0511","2019nCoV301-US200-0046","2019nCoV301-US228-0067"
                                           ,"2019nCoV301-US228-0210","2019nCoV301-US239-0100","2019nCoV301-US241-0121","2019nCoV301-US242-0118"
))$Day35pseudoneutid50
dat.tmp.2=subset(dat.2.pr.old, Ptid %in% c("2019nCoV301-US049-0039","2019nCoV301-US159-0063","2019nCoV301-US159-0135","2019nCoV301-US162-0231"
                                           ,"2019nCoV301-US163-0193","2019nCoV301-US172-0511","2019nCoV301-US200-0046","2019nCoV301-US228-0067"
                                           ,"2019nCoV301-US228-0210","2019nCoV301-US239-0100","2019nCoV301-US241-0121","2019nCoV301-US242-0118"
))$Day35pseudoneutid50


with(subset(dat.2.pr, Trt==1), mytable(EventIndPrimaryD57, Country, ph1.D57))
with(subset(dat.2.pr, Trt==1 & ph1.D57==1), mytable(EventIndPrimaryD57, Country, ph2.D57))

with(subset(dat.2.pr, Trt==1 & Bserostatus==0), mytable(EventIndPrimaryD35, ph1.D35))

with(subset(dat.2.pr.new, Trt==1 & SubcohortInd==1 & Country==0), mytable(!is.na(BbindSpike)&!is.na(Day35bindSpike), ph1.D35))

with(subset(dat.mock, Country==0 & SubcohortInd), mytable(Bserostatus, Trt))
with(subset(dat.mock, Country==0 & SubcohortInd & Bserostatus==0 & Trt==1), mytable(!is.na(BbindSpike), !is.na(Day35bindSpike)))
with(subset(dat.mock, Country==0 & SubcohortInd & Bserostatus==0 & Trt==1), mytable(!is.na(Bpseudoneutid50), !is.na(Day35pseudoneutid50)))
with(subset(dat.2.pr.new, Country==0 & SubcohortInd & Bserostatus==0 & Trt==1 & ph1.D35 & !is.na(BbindSpike) & !is.na(Day35bindSpike)), mytable(is.na(Day35pseudoneutid50), EventIndPrimaryD35))

with(subset(dat.2.pr.new, Country==0 & SubcohortInd==1 & Bserostatus==0 & Trt==1 & Perprotocol ), mytable(!is.na(BbindSpike)&!is.na(Day35bindSpike)) )
with(subset(dat.2.pr.new, Country==0 & SubcohortInd==1 & Bserostatus==0 & Trt==1 & Perprotocol & ph2.immuno ), mytable(ph2.D35) )
subset(dat.2.pr.new, Country==0 & SubcohortInd==1 & Bserostatus==0 & Trt==1 & Perprotocol & ph2.immuno & ph2.D35 & AnyinfectionD1==1)


# downloaded data
dat.new=read.csv("/trials/covpn/p3004/download_data/Neutralizing_Antibody/2019nCoV2-301_ nAB_D614G_14Apr2022.csv")
dat.old=read.csv("/trials/covpn/p3004/download_data/Neutralizing_Antibody/2019nCoV2-301_ nAB_D614G_30Mar2022.csv")

all(dat.old$Monogram.Accession %in% dat.new$Monogram.Accession)

dat.tmp = dat.new[!dat.new$Monogram.Accession %in% dat.old$Monogram.Accession, ]

dat.new$Subjectid = "2019nCoV301-" %.% dat.new$Patient.Number
dat.old$Subjectid = "2019nCoV301-" %.% dat.old$Patient.Number

dat.tmp.1=subset(dat.new, Subjectid %in% c("2019nCoV301-US049-0039","2019nCoV301-US159-0063","2019nCoV301-US159-0135","2019nCoV301-US162-0231"
                                           ,"2019nCoV301-US163-0193","2019nCoV301-US172-0511","2019nCoV301-US200-0046","2019nCoV301-US228-0067"
                                           ,"2019nCoV301-US228-0210","2019nCoV301-US239-0100","2019nCoV301-US241-0121","2019nCoV301-US242-0118"
))[,c(4:10)]
dat.tmp.2=subset(dat.old, Subjectid %in% c("2019nCoV301-US049-0039","2019nCoV301-US159-0063","2019nCoV301-US159-0135","2019nCoV301-US162-0231"
                                           ,"2019nCoV301-US163-0193","2019nCoV301-US172-0511","2019nCoV301-US200-0046","2019nCoV301-US228-0067"
                                           ,"2019nCoV301-US228-0210","2019nCoV301-US239-0100","2019nCoV301-US241-0121","2019nCoV301-US242-0118"
))[,c(4:10)]

dat.tmp.1==dat.tmp.2




