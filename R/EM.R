norm.mix.den <- function(data, L, probL, muL, sigmasqL )
  {
    temp <- 0
    for(l in 1:L)
      temp <- temp + probL[l] * dnorm( data, muL[l], sqrt(sigmasqL[l]) )

    temp
  }

EM.Esti <- function(x, L, DELTA=0.001,max.iter=500)
  {
    n <- length(x)
    probL <- runif(L)
    probL <- probL/sum(probL)
    
    muL <- rnorm(L)
    sigmaSqL <- rchisq(L, 2, 2)
    pi1 <- 0.5

    muL.new <- muL
    sigmaSqL.new <- sigmaSqL
    pi1.new <- pi1
    probL.new <- probL
    
    delta <- 1
    omega.mat <- array(0, c( n, L) )
    lfdr.tmp <- array(0, n)
    iter <- 1
    
    while( delta > DELTA)
      {
        if( iter  > max.iter )
          break
        
        for(i in 1:n)
          {
            mixden <- norm.mix.den(x[i], L, probL, muL, sigmaSqL )
            lfdr.tmp[i] <- (1-pi1)*dnorm(x[i],0,1)/( (1-pi1)*dnorm(x[i],0,1) + pi1 * mixden )
            for(l in 1:L )
              omega.mat[i,l] <- probL[l]*dnorm(x[i], muL[l],sqrt(sigmaSqL[l]))/mixden
          }
        pi1.new <- 1 - mean( lfdr.tmp )
        for(l in 1:L)
          {
            probL.new[l] <- sum( (1-lfdr.tmp) * omega.mat[,l] )/sum(1- lfdr.tmp)
            muL.new[l] <- sum( (1-lfdr.tmp)*omega.mat[,l]*x )/sum( (1-lfdr.tmp)*omega.mat[,l] )
            sigmaSqL.new[l] <- sum( (1-lfdr.tmp)*omega.mat[,l]*(x-muL.new[l])^2)/sum( (1-lfdr.tmp)*omega.mat[,l] ) 
          }
        delta <- 0
        delta <- delta + (pi1.new-pi1)^2
        for(l in 1:L)
          delta <- delta + (probL.new[l]-probL[l])^2 + (muL.new[l]-muL[l])^2 + (sigmaSqL.new[l]-sigmaSqL[l])^2

        if( is.na(delta) )
          {            
            probL.new <- runif(L)
            probL.new <- probL/sum(probL)
            
            muL.new <- rnorm(L)
            sigmaSqL.new <- rchisq(L, 2, 2)
            pi1.new <- 0.5
            delta <- 1
            iter <- 0
          }

        pi1 <- pi1.new
        probL <- probL.new
        muL <- muL.new
        sigmaSqL <- sigmaSqL.new
        iter <- iter + 1

      }

    est <- list( pi1=pi1, probL=probL, muL=muL, sigmaSqL=sigmaSqL )
  }


LFDR.EM <- function(X, L=2, DELTA=0.001)
  {
    n <- length(X)
    em.esti <- EM.Esti(X, L=2, DELTA=0.001)
    lfdr <- array(0,n)
    for(i in 1:n)
      {
        lfdr[i] <- ( 1-em.esti$pi1 ) * dnorm(X[i])/( (1-em.esti$pi1) * dnorm(X[i])  + em.esti$pi1 * norm.mix.den(X[i], L, em.esti$probL, em.esti$muL, em.esti$sigmaSqL) )
      }
    lfdr
  }

LFDR.EM.output <- function(X, L, DELTA=0.001)
  {
    n <- length(X)
    em.esti <- EM.Esti(X, L, DELTA=0.001)
    lfdr <- array(0,n)
    for(i in 1:n)
      {
        lfdr[i] <- ( 1-em.esti$pi1 ) * dnorm(X[i])/( (1-em.esti$pi1) * dnorm(X[i])  + em.esti$pi1 * norm.mix.den(X[i], L, em.esti$probL, em.esti$muL, em.esti$sigmaSqL) )
      }
    list( lfdr=lfdr, pi1=em.esti$pi1, probL = em.esti$probL, muL=em.esti$muL, sigmaSqL=em.esti$sigmaSqL)
  }
