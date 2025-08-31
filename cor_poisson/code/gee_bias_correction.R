# Bias corrections to the sandwich variance estimator for small samples (originally in ~adecamp/thesis/general/code/gagStudySize_small_samp.R)

# Small Sample Adjustments for Sandwich Estimator of variance in GEE
# 
# Argument fit.gee is the fit object from a call to geeglm from package geepack.
# The function is currently only implemented for family="gaussian" with link="identity" and corstr="independence","exchangeable", or "ar1"
# Returns a list with elements MBN, MD, KC, and mFG each of which is a corrected sandwich estimate of the variance
# of the beta coefficients from fit.gee.  Additionally, returns the standard sandwich estimate.
#   0) Standard sandwich estimate with no correction (NO), should correspond to fit.gee$geese$vbeta
#   1) Morel, Bokossa, and Neerchal (MBN)
#   2) Mancl and DeRouen (MD)
#   3) Kauermann and Carroll (KC)
#   4) modified Fay and Graubard (mFG)

# change log

# changed code to allow "exchangeable" and "ar1" working correlations, note waves argument is not
# implimented which would impact the results for "ar1" if specified in the call to geeglm 
# Allan deCamp (July 24, 2013)

# NOTES
#
# to grab the waves argument from a call to geeglm consider the code below
#   mf <- fit$call
#   mf[[1]] <- as.name("model.frame")
#   mftmp <- mf
#   mftmp$family <- mftmp$corstr <- mftmp$control <- mftmp$zcor <- mftmp$std.err <- NULL
#   mftmp$scale.fix <- NULL
#   mf <- eval(mftmp, parent.frame())
#   waves <- model.extract(mf, waves)
# 
# Future iterations might consider an implementation that incorporates ideas from 
# "Fast Pure R Implementation of GEE: Application of the Matrix Package"
# http://journal.r-project.org/archive/2013-1/mcdaniel-henderson-rathouz.pdf


ssase = function(fit.gee) {

  stopifnot(is(fit.gee)=="geeglm")
  stopifnot(fit.gee$corstr %in% c("independence","exchangeable","ar1")) 
  # stopifnot(fit.gee$family$family=="gaussian")
  # stopifnot(fit.gee$family$link=="identity")
  stopifnot(is.null(fit.gee$call$waves))
  
  # the standard sandwich estimate uses alpha, gamma, A and B
  alpha = fit.gee$geese$alpha
  gamma = fit.gee$geese$gamma
  p = ncol(fit.gee$geese$X)
  corstr = fit.gee$corstr
  clusz = fit.gee$geese$clusz
  ids   = rep(1:length(clusz), clusz)
  uid = 1:length(clusz)
  
  N = max(clusz)
  if( corstr=="independence" ) {
    SIG = diag(nrow=N)
  } else if( corstr=="exchangeable" ) {
    SIG = matrix(alpha, nrow=N, ncol=N)
    diag(SIG) = 1
  } else if( corstr=="ar1" ) {
    SIG = matrix(NA, nrow=N, ncol=N)
    for( i in 1:N ) {
      for( j in 1:N ) {
        SIG[i,j] = alpha^abs(i-j)
      }
    }
  } else {
    stop(sprintf("unexpected corstr %s",corstr))
  }
  SIG = gamma * SIG
  
  sig = lapply(uid, function(id) {
    n = clusz[id]
    return(SIG[1:n, 1:n])
  })
  siginv = lapply(sig, function(m) solve(m))
    
  B = sapply(uid, function(id) {
    w = which(ids==id)
    r = fit.gee$residuals[w]
    r = matrix(r, nc=1)
    ohat = r %*% t(r)
    D = fit.gee$geese$X[w,,drop=F]
    t(D) %*% siginv[[id]] %*% ohat %*% siginv[[id]] %*% D
  })
  B = matrix(rowSums(B), nr=p)

  A = sapply(uid, function(id) {
    w = which(ids==id)
    D = fit.gee$geese$X[w,,drop=F]
    t(D) %*% siginv[[id]] %*% D
  })
  A = matrix(rowSums(A), nr=p)

  vbeta.NO = solve(A) %*% B %*% solve(A) # standard sandwich estimate, no small sample correction

  # Morel (MBN)
  dhat = min(0.5,p/(length(uid)-p))
  zhat = max(1,sum(diag(solve(A) %*% B))/p)

  vbeta.MBN = solve(A) %*% B %*% solve(A) + dhat*zhat*solve(A)

  # Mancl
  BMD = sapply(uid, function(id) {
    w = which(ids==id)
    r = fit.gee$residuals[w]
    r = matrix(r, nc=1)
    ohat = r %*% t(r)
    D = fit.gee$geese$X[w,,drop=F]
    H = D %*% solve(A) %*% t(D) %*% siginv[[id]]
    I = diag(rep(1,nrow(H)))
    t(D) %*% siginv[[id]] %*% solve(I-H) %*% ohat %*% solve(I-H) %*% siginv[[id]] %*% D
  })
  BMD = matrix(rowSums(BMD), nr=p)

  vbeta.MD = solve(A) %*% BMD %*% solve(A)

  # Kauermann and Carroll (I-H)^-.5 as described in Lu et al. 
  # could replace (I-H)^-1 in Mancl and DeRouen with (I-H)^-.5 but running into 
  # problems taking the inverse of the square root...
  BKC = sapply(uid, function(id) {
    w = which(ids==id)
    r = fit.gee$residuals[w]
    r = matrix(r, nc=1)
    ohat = r %*% t(r)
    D = fit.gee$geese$X[w,,drop=F]
    H = D %*% solve(A) %*% t(D) %*% siginv[[id]]
    I = diag(rep(1,nrow(H)))
    M = solve(I-H)
    e = eigen(M)
    M = e$vectors %*% diag(sqrt(e$values)) %*% solve(e$vectors)
    
    t(D) %*% siginv[[id]] %*% M %*% ohat %*% M %*% siginv[[id]] %*% D
  })
  BKC = matrix(rowSums(BKC), nr=p)

  vbeta.KC = solve(A) %*% BKC %*% solve(A)

  # modified FG
  BmFG = sapply(uid, function(id) {
    w = which(ids==id)
    r = fit.gee$residuals[w]
    r = matrix(r, nc=1)
    ohat = r %*% t(r)
    D = fit.gee$geese$X[w,,drop=F]

    # for inverse of square root of a matrix...
    I = diag(nrow=p)
    M = I-t(D) %*% siginv[[id]] %*% D %*% solve(A)
    e = eigen(M)
    H = solve(e$vectors %*% diag(sqrt(e$values)) %*% solve(e$vectors))
          
    H %*% t(D) %*% siginv[[id]] %*% ohat %*% siginv[[id]] %*% D %*% t(H)
  })
  BmFG = matrix(rowSums(BmFG), nr=p)

  vbeta.mFG = solve(A) %*% BmFG %*% solve(A)

  list(NO=vbeta.NO, MBN=vbeta.MBN, MD=vbeta.MD, KC=vbeta.KC, mFG=vbeta.mFG)
}

####################################################################################################################
### quick sanity checks, either run function or copy and paste code...
sanity_check = function() {
  require(geepack)
  require(mvtnorm)

  # generate clustered data
  np=11
  nv=17
  n = np+nv

  rho = 0.9
  mcor = matrix(rho, nr=12, nc=12)
  diag(mcor) = 1

  m.samp = rpois(np+nv, lambda=6.5)
  m.samp = pmax(1,m.samp)
  m.samp = pmin(12,m.samp)
    
  df.samp = data.frame( seqid=rep(1:(np+nv), times=m.samp), 
                        trt=c(rep("placebo", sum(m.samp[1:np])), rep("vaccine", sum(m.samp[(np+1):(np+nv)]))) )
  dist = rmvnorm(n, sigma=mcor)
  df.samp$dist = unlist(sapply(1:n, function(row) dist[row,1:m.samp[row]]))


  # run GEE
  fit.gee = geeglm(df.samp$dist~df.samp$trt, id=df.samp$seqid, corstr="independence")
  vbeta = ssase(fit.gee)

  # show variance estimates
  fit.gee$geese$vbeta
  vbeta$NO
  vbeta$MBN
  vbeta$MD
  vbeta$KC
  vbeta$mFG

  # No adjustment and variance computed by geeglm should be the same
  max(abs(vbeta$NO - fit.gee$geese$vbeta))


  # run GEE 
  for( corstr in c("ind", "exch", "ar1") ) {
    fit.gee = geeglm(df.samp$dist~df.samp$trt, id=df.samp$seqid, corstr=corstr)
    vbeta = ssase(fit.gee)

    # No adjustment and variance computed by geeglm should be the same
    print(fit.gee$corstr)
    print(fit.gee$geese$alpha)
    print(max(abs(vbeta$NO - fit.gee$geese$vbeta)))
    print(max(abs((vbeta$NO - fit.gee$geese$vbeta)/fit.gee$geese$vbeta)))

    # KC and mFG typically give nearly the same result (mathematically identical, differences due to instability in 
    # matrix manipulations...)
    print(max(abs(vbeta$KC-vbeta$mFG)))

  }



  ### using Michal's data
  data = read.table("~/test_data.txt", header=T)
  fit1 = geeglm(Y.parallel ~ X.parallel + J + X.parallel:S.parallel, id=cluster, family=gaussian("identity"), data=data)
  vbeta = ssase(fit1)

  # run GEE 
  for( corstr in c("ind", "exch", "ar1") ) {
    fit1 = geeglm(Y.parallel ~ X.parallel + J + X.parallel:S.parallel, id=cluster, family=gaussian("identity"), data=data, corstr=corstr)
    vbeta = ssase(fit1)

    # No adjustment and variance computed by geeglm should be the same
    cat(sprintf("corstr %s\n", fit1$corstr))
    cat(sprintf("maximal absolute difference in geese vbeta and ssase vbeta$NO %0.2g\n", max(abs(vbeta$NO - fit1$geese$vbeta))))
    cat(sprintf("maximal relative difference in geese vbeta and ssase vbeta$NO %0.2g\n", max(abs((vbeta$NO - fit1$geese$vbeta)/fit1$geese$vbeta))))
  }

}





