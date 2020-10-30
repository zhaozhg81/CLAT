                                        #' The main routine for CLAT
                                        #'
                                        #' This function loads a vector of p-values and q, the designated mFDR level. 
                                        #'
                                        #' @param pv A vector of pvalues.
                                        #' @param q The designated mFDR level
                                        #' @return R: The total number of rejection.
                                        #' @return SigInd: A vector of 0-1, indicating whether the corresponding hypothesis is rejected (1) or is not rejected (0).
                                        #' @return 
                                        #' @examples
                                        #' 
                                        #' 
                                        #' 
                                        #' p <- 1000
                                        #' 
                                        #' X <- c( rnorm(100, -2, 1), rnorm(100, 2, 1), rnorm(800, 0, ) )
                                        #' 
                                        #' q <- 0.05
                                        #' 
                                        #' ## Right-sided test
                                        #' p.right <- pnorm( X, lower.tail=FALSE )
                                        #' clat.right <- CLAT( p.right, q)
                                        #' rejind <- which( clat.right$SigInd!=0 )
                                        #' 
                                        #' ## Right-sided test
                                        #' p.left <- pnorm( X, lower.tail=TRUE )
                                        #' clat.left <- CLAT( p.left, q)
                                        #' rejind <- which( clat.left$SigInd!=0 )
                                        #' 
                                        #' 
                                        #' ## Two-sided test
                                        #' p.right <- pnorm( X, lower.tail=FALSE )
                                        #' clat.right <- CLAT( p.right, q)
                                        #' 
                                        #' p.left <- pnorm( X, lower.tail=TRUE )
                                        #' clat.left <- CLAT( p.left, q)
                                        #' rejind <- c( which( clat.left$SigInd!=0 ), which(clat.right$SigInd!=0) ) 
                                        #' @export


CLAT <- function(pv, q){## Positive Side
  
  n <- length(pv)
  i=1
  j=1

  R=0
  SigInd=rep(0,n)
  Rrange=c(Inf,Inf)
  pv.sort <- sort(pv, decreasing=FALSE)
  X.sort <- qnorm( pv.sort )##*1.53
  Ti <- q*c(1:n)/n - pv.sort
  li <- order(Ti, decreasing=FALSE)
  
  
  MAX=0
  iPrime=li[1]
  lower=0; ## lower is I in the paper
  upper=0; ## upper is J in the paper
  
  for( j in 1:n)
    {
      jPrime <- li[j]

      if( pv.sort[jPrime]>0.5)
        next
      if(   (jPrime - iPrime> MAX ) & ( ( (abs( X.sort[jPrime]-X.sort[iPrime])>(2*log(n)/sqrt(n)) ))|( abs(X.sort[jPrime])==Inf)|( abs(X.sort[iPrime])==Inf) )  ){
      ## if(   (jPrime - iPrime> MAX ) & ( (abs( X.sort[jPrime]-X.sort[iPrime])>0.5 ))|( abs(X.sort[jPrime])==Inf)|( abs(X.sort[iPrime])==Inf)   ){
          
        MAX=jPrime - iPrime
        lower=iPrime
        upper=jPrime
      }
      if( jPrime < iPrime){
        iPrime=jPrime
      }
      if( ( Ti[jPrime] >=0 ) & (jPrime>MAX) ){
        lower=0
        upper=jPrime
        MAX=jPrime
      }
    }
  if( upper - lower > 1 ){
    if(lower>0){
      Rrange=c(pv.sort[lower],pv.sort[upper])
    }else{
      Rrange=c(0, pv.sort[upper] )
    }
    SigInd <- (pv>=Rrange[1])*(pv<=Rrange[2])
    R <- sum(SigInd)
  }else{
    Rrange=c(0,0)
    SigInd <- 0*pv
    R <- 0
  }
  
  y <- list(R=R,SigInd=SigInd)
  
  return(y)
}
