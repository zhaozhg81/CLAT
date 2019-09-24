source("CDF.R")


CLAT <- function(X, q, mu=0, s=1)
  {
    clat.pos <- CLAT.pos(X,q,mu,s)
    clat.neg <- CLAT.neg(X,q,mu,s)

    R <- clat.pos$R + clat.neg$R
    SigInd <- clat.pos$SigInd + clat.neg$SigInd

    y <- list(R=R,SigInd=SigInd)

    return(y)
  }

CLAT.pos <- function(X,q,mu=0,s=1){
  ## X: is the test statistic
  ## q: the FDR level
  ## mu, s: the parameter for the null distribution. By default, under the null, X~ N(0,1). But mu, s can be estimated empirically.
  
  cdfEsti <- CDF(X,mu,s)
  n <- length(X)
  i=1
  j=1

  R=0
  SigInd=rep(0,n)
  Rrange=c(Inf,Inf)
  phiCDF <-pnorm(X,mu,s,lower.tail=TRUE)
    CC <- (q*cdfEsti$empCDF-phiCDF)
    ORDER <- order(CC)
    CSort <- CC[ORDER]
    
    A <- cdfEsti$empCDF
    OriInd <- ORDER

    MAX=0
    iPrime=ORDER[1];
    lower=iPrime;
    upper=iPrime;
  
    for( j in 1:n)
      {
        jPrime <- OriInd[j]
        if(X[jPrime]<0)
        	next
        if((A[jPrime]-A[iPrime]>MAX) & (abs(X[jPrime]-X[iPrime])>0.01)){
          MAX=A[jPrime]-A[iPrime]
          lower=iPrime
          upper=jPrime
        }
        if(A[jPrime]<A[iPrime]){
          iPrime=jPrime
        }
      }
    Rrange=c(X[lower],X[upper])
    
    SigInd <- (X>=Rrange[1])*(X<=Rrange[2])
    
  R <- sum(SigInd)
  R <- R*(R>1)
  y <- list(R=R,SigInd=SigInd,Rrange=Rrange)

  return(y)
}


CLAT.neg <- function(X,q,mu=0,s=1){
  return( CLAT.pos (-X,q,-mu,s) )
}
