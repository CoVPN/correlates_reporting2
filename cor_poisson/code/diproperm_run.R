# remotes::install_github("elbamos/clusteringdatasets")
# remotes::install_github("jiadongm/diproperm")

library(glue)
library(diproperm)
library(config)

file = config::get(config = "iliad_ib202p", file="~/correlates_reporting2/config.yml")$data_cleaned
dat_proc_202=read.csv(file)
dat_proc_202$Ptid = as.factor(dat_proc_202$Ptid) # need to make it a factor to use as id
dat_proc_202$Trt = factor(dat_proc_202$Trt, levels = c("PBO","BPZE1")) 



dat = subset(dat_proc_202, Trt == "BPZE1" & PPAI==1)

panels=c("Nasal_IgA",      "Norm_Nasal_IgA", "Serum_IgA",      "Serum_IgG")
antigens=c("WCE", "FHA", "FIM", "PRN", "PT")

pvals = 
sapply(c("Day28", "B", "Delta28overB"), function (t) {
  sapply(panels, function(p) {
    m <- matrix(unlist(dat[,glue("{t}{antigens}_{p}")]), nc=5)
    # classifier choice: md is extrememly fast; dwd is very slow
    out=DiProPerm (X=m, y=ifelse(dat$EventIndC9_11_14==1,1,-1), B=10000, classifier="md", univ.stat="md", cores=10)
    out$pvalue
  })
})
pvals
