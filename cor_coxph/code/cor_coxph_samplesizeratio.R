renv::activate(project = here::here(".."))     
library(kyotil)
library(Hmisc)
library(xtable) # this is a dependency of kyotil

Args <- commandArgs(trailingOnly=TRUE) 
if (length(Args)!=2) {
    stop("Two command line arguments are required: study_name, COR. For example, COVE D57 or ENSEMBLE D29IncludeNotMolecConfirmedstart1")
}
myprint(Args)
study.name <- Args[1]
COR=Args[2]

if (study.name=="COVE") {
    load(paste0("output/moderna_real/", COR, "/marginalized.risk.Rdata"))
    tmp=lapply(risks.all.1, function(x) x$n.dean)
    
    
} else if (study.name=="ENSEMBLE") {
    
} 
    
assays=  

# one at a time
get.n.ratio.1=function(beta.1, beta.2, version=NULL) {
    
    if        (beta.1>0 & beta.2>0) {
        (beta.1/beta.2)**2
        
    } else if (beta.1<0 & beta.2>0) {
        0
        
    } else if (beta.1<0 & beta.2<0) {
        if(is.null(version)) NA else if (version=="L") 0 else if (version=="U") Inf else stop("wrong version in get.n.ratio")
        
    } else if (beta.1>0 & beta.2<0) {
        Inf
        
    } else stop("something wrong in get.n.ratio")        
}

# vectorized
get.n.ratio=function(beta.1, beta.2, version=c("","L","U")) {
    match(version)
    ifelse (beta.1>0 & beta.2>0, 
        (beta.1/beta.2)**2,
        
    ifelse (beta.1<0 & beta.2>0, 
        0, 
        
    ifelse (beta.1<0 & beta.2<0,
        if(version=="") NA else if (version=="L") 0 else if (version=="U") Inf,
        
    ifelse (beta.1>0 & beta.2<0, 
        Inf, -Inf) # -Inf would only appear when either beta is 0
        
        )
        )
        )
                
}


x=rnorm(10); y=rnorm(10)
get.n.ratio(x,y)
sapply(1:length(x), function(i) get.n.ratio.1(x[i],y[i]) )


  
for(i in 1:(length(assays)-1)) {
    for(j in 2:length(assays)) {
        get.n.ratio(beta.1, beta.2)
        sapply(, simplify="array", function())
        
        
        
    }
}
