# one at a time, only used for debugging
get.n.ratio.1=function(beta.1, beta.2, version=c("","L","U")) {   
    version=match.arg(version)
    if        (beta.1<0 & beta.2<0) {
        (beta.1/beta.2)**2
    } else if (beta.1>0 & beta.2<0) {
        0
    } else if (beta.1>0 & beta.2>0) {
        if(version=="") NA else if (version=="L") 0 else if (version=="U") Inf else stop("wrong version in get.n.ratio")
    } else if (beta.1<0 & beta.2>0) {
        Inf
    } else stop("something wrong in get.n.ratio")        
}

# vectorized
get.n.ratio=function(beta.1, beta.2, version=c("","L","U")) {
    version=match.arg(version)
    ifelse (beta.1<0 & beta.2<0, 
        (beta.1/beta.2)**2,
    ifelse (beta.1>0 & beta.2<0, 
        0, 
    ifelse (beta.1>0 & beta.2>0,
        if(version=="") NA else if (version=="L") 0 else if (version=="U") Inf,
    ifelse (beta.1<0 & beta.2>0, 
        Inf, -Inf) # -Inf would only appear when either beta is 0
        )
        )
        )           
}
# test
#x=rnorm(10); y=rnorm(10); get.n.ratio(x,y); sapply(1:length(x), function(i) get.n.ratio.1(x[i],y[i]) )

n.dean=sapply(risks.all.1, function(x) x$n.dean)
assays=colnames(n.dean)
k=ncol(n.dean)  

samplesizeratio <- samplesizeratio.lb <- samplesizeratio.ub <- matrix(NA, k, k, dimnames=list(assays, assays))
for(i in 1:k) {
    for(j in setdiff(1:k, i)) {
        samplesizeratio[i, j]=get.n.ratio(n.dean[1,i],  n.dean[1,j])
                           lb=get.n.ratio(n.dean[-1,i], n.dean[-1,j], "L")
                           ub=get.n.ratio(n.dean[-1,i], n.dean[-1,j], "U")
        samplesizeratio.lb[i, j]=quantile(lb, 2.5/100)
        samplesizeratio.ub[i, j]=quantile(ub, 1-2.5/100)
    }
}

mytex(samplesizeratio, file.name="follmann2018_samplesizeratio", align="c", include.colnames = T, save2input.only=T, input.foldername=save.results.to, digits=2)

tab=matrix(paste0("(", formatDouble(samplesizeratio.lb, 2), ", ", formatDouble(samplesizeratio.ub, 2), ")"), k, k, dimnames=list(assays, assays))
mytex(tab, file.name="follmann2018_samplesizeratio_lbub", align="c", include.colnames = T, save2input.only=T, input.foldername=save.results.to)
