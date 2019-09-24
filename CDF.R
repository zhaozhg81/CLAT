CDF <- function(X,mu=0,s=1){
  p <- length(X)
  empCDF <- X
  phiCDF <- X
  YY=ecdf(X)
  empCDF<- YY(X)
  #phiCDF <- pnorm(X,mu,s,lower.tail=TRUE)
  y <- list(empCDF=empCDF)
  return(y)
}
