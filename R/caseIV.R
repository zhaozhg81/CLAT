                                        #' This function runs the simulation and plots the figure for case IV.
                                        #'
                                        #' @param beta The value of beta, which is between 0 and 1. In the paper, beta is chosen as 0.3 and 0.4 respectively.
                                        #' @export

caseIV <- function(beta)
  {
    
    numSim <- 2
    
    n <- 5000
    
    sigma <- 0.7    
    p1 <- 0.9
    q <- 0.1
    
    
    NONzero=floor(n^{1-beta})
    epsilon <- NONzero/n;
    
    qPrime <- epsilon*q/(1-(1-epsilon)*q)
    if(sigma<1){
      muLB <- sqrt(2*(1-sigma^2)*log(sigma/(p1*qPrime)))
      mu <- c( muLB+0.2+c(4:20)/40, muLB + c(14:20)/20, muLB + c(6:8)/5)
    }
    if(sigma>1){
      mu=c(1.5,2,2.5,3,3.5)+1.5
      mu=mu-1.5
    }
    
    rNum <- length(mu)
    
    set.seed(3)
    
    FDP <- array(0,c(rNum,numSim)); FNP <- FDP; FDPBH <- FDP; R <- FDP
    FDPSunCai <- FDP; FNPSunCai <- FDP; RSunCai <- FDP
    FDPBH <- FDP; FNPBH <- FDP; RBH <- FDP
    FDPIdeal <- FDP; FNPIdeal <- FDP; RIdeal <- FDP
    FDPlocfdr <- FDP; FNPlocfdr <- FDP; Rlocfdr <- FDP
    FDPlocfdr.em <- FDP; FNPlocfdr.em <- FDP; Rlocfdr.em <- FDP
    FDPoracle <- FDP; FNPoracle <- FDP; Roracle <- FDP
    
    
    NoFalseSunCai <- array(0,c(rNum,numSim,2))
    
    FDR <- array(0,c(1,rNum)); FNR <- FDR; AveRej <- FDR
    FDRSunCai <- FDR; FNRSunCai <- FDR; AveRejSunCai <- FDR
    FDRBH <- FDR; FNRBH <- FDR; AveRejBH <- FDR
    FDRIdeal <- FDR; FNRIdeal <- FDR; AveRejIdeal <- FDR
    FDRlocfdr <- FDR; FNRlocfdr <- FDR; AveRejlocfdr <- FDR
    FDRlocfdr.em <- FDR; FNRlocfdr.em <- FDR; AveRejlocfdr.em <- FDR
    FDRoracle <- FDR; FNRoracle <- FNR; AveRejoracle <- FDR
    
    MaxX <- array(0,c(rNum,numSim))
    
    sigma2 <- 0.5
    
    
    EM.esti <- array(0, c(rNum, numSim, 7 ) )
    
    for (l in 1:rNum){
      for (k in 1:numSim){
        
        ## ########################################################
        ## #####################################################
        ## ## Generate the random numbers
        X <- array(n,c(1,n))
        trueInd <- array(n,c(1,n))
        trueInd[1:(n-NONzero)] <- 0
        trueInd[(n-NONzero+1):n] <- 1
        X[1:(n-NONzero)] <- rnorm(n-NONzero,0,1)
        X[(n-NONzero+1):n] <- rnorm(NONzero,c(rep(mu[l],floor(NONzero*p1)),rep(-mu[l],NONzero-floor(NONzero*p1))),sigma)
        ZZ <- rnorm(n, 0, sigma2)
        X <- (X+ZZ)/sqrt( 1 + sigma2^2 )
        
        f0x <- dnorm(X)
        f1x <- p1 * dnorm( X, mu[l], sigma ) + (1-p1) * dnorm(X, -mu[l], sigma)
        
        ## ###########################################################
        ## ###############################################################
        ## ## Sun and Cai 's approach
        adaptiveZ<-adaptZ.func(X, q)
        ## the threshold
        threshold<-adaptiveZ$th
        ## number of rejected hypotheses
        RSunCai[l,k]<-adaptiveZ$nr
        ## the rejected hypotheses
        rh<-adaptiveZ$re
        FDPSunCai[l,k] <- sum(rh<(n*(1-epsilon)))/(RSunCai[l,k]+(RSunCai[l,k]==0))
        FNPSunCai[l,k] <- (n*epsilon-sum(rh>(n*(1-epsilon))))/(n-RSunCai[l,k])
        ## ##########################################################
        ## ##########################################################
        ## ## CLAT 
        
        ## # Positve side
        if(p1!=1){
          pv <- 1 - pnorm( X )
          CLATPos <- CLAT(pv, q)
          
          pv.neg <- pnorm(X)
          CLATNeg <- CLAT(pv.neg, q)
          
          R[l,k] <- CLATPos$R + CLATNeg$R   
          SigInd = CLATPos$SigInd + CLATNeg$SigInd
          FDP[l,k] <- sum(SigInd*(1-trueInd))/(R[l,k]+(R[l,k]==0))
          FNP[l,k] <- sum((1-SigInd)*trueInd)/(n-R[l,k])
        }
        if(p1==1){
          CLATres <- CLAT(X,q,'Normal')
          R[l,k] <- CLATres$R
          FDP[l,k] <- sum( CLATres$SigInd*(1-trueInd))/(R[l,k]+(R[l,k]==0))
          FNP[l,k] <- sum((1-CLATres$SigInd)*trueInd)/(n-R[l,k])
        }
        

        ## ################################################################
        ## #########################################################
        ##    oracle
        LFDR.oracle <- (1-epsilon)*f0x/( (1-epsilon)*f0x + epsilon * f1x )
        
        LFDR.sort <- sort(LFDR.oracle, decreasing=FALSE)
        Rtmp <- max( c(1:n) * ( ( 1/c(1:n) * cumsum( LFDR.sort) )< q ) )
        
        Roracle[l,k] <- Rtmp
        
        if( Rtmp> 0){
          Indoracle <- (LFDR.oracle <= LFDR.sort[Rtmp]);
          FDPoracle[l,k] <- sum(Indoracle*(1-trueInd))/(Roracle[l,k]+(Roracle[l,k]==0));
          FNPoracle[l,k] <- sum((1-Indoracle)*trueInd)/(n-Roracle[l,k]);
        }else{
          FDPoracle[l,k] <- 0
          FNPoracle[l,k] <- sum(trueInd)/n;
        }
        
        
        ## ################################################################
        ## #########################################################
        ##    localfdr (Efron)
        LFDR <- locfdr(X, plot=0);
        LFDR.sort <- sort(LFDR$fdr, decreasing=FALSE)
        R.lfdr <- max( c(1:n) * ( ( 1/c(1:n) * cumsum( LFDR.sort) )< q ) )
        
        Rlocfdr[l,k] <- R.lfdr
        
        if( R.lfdr> 0){
          Indlocfdr <- (LFDR$fdr <= LFDR.sort[R.lfdr]);
          FDPlocfdr[l,k] <- sum(Indlocfdr*(1-trueInd))/(Rlocfdr[l,k]+(Rlocfdr[l,k]==0));
          FNPlocfdr[l,k] <- sum((1-Indlocfdr)*trueInd)/(n-Rlocfdr[l,k]);
        }else{
          FDPlocfdr[l,k] <- 0
          FNPlocfdr[l,k] <- sum(trueInd)/n;
        }
    
        ## ################################################################
        ## #########################################################
        ## ##   EM Algorithm for estimating the true parameter
        lfdr.res <- LFDR.EM.output(X[1,], L=2, DELTA=0.001)
        
        lfdr.em <- lfdr.res$lfdr
        
        EM.esti[l,k, 1] <- lfdr.res$pi1
        EM.esti[l,k, c(2:3 )] <- lfdr.res$probL
        EM.esti[l,k, c(4:5)] <- lfdr.res$muL
        EM.esti[l,k, c(6:7)] <- lfdr.res$sigmaSqL
        
        lfdr.em.sort <- sort(lfdr.em, decreasing=FALSE)
        R.lfdr.em <- max( c(1:n) * ( ( 1/c(1:n) * cumsum( lfdr.em.sort) )< q ) )
        
        Rlocfdr.em[l,k] <- R.lfdr.em
        
        if( R.lfdr.em> 0){
          Indlocfdr.em <- (lfdr.em <= lfdr.em.sort[R.lfdr.em]);
          FDPlocfdr.em[l,k] <- sum( Indlocfdr.em*(1-trueInd[1,]))/(Rlocfdr.em[l,k]+(Rlocfdr.em[l,k]==0));
          FNPlocfdr.em[l,k] <- sum((1-Indlocfdr.em)*trueInd[1,])/(n-Rlocfdr.em[l,k]);
        }else{
          FDPlocfdr.em[l,k] <- 0
          FNPlocfdr.em[l,k] <- sum(trueInd)/n;
        }
        
        
        
        ## ################################################################
        ## #########################################################
        ## ##   Benjamini Hochberg procedure
        if(p1==1){
          pv <- 1-pnorm(X)
        }
        if(p1!=1){
          pv <- 2*pnorm(-abs(X))
        }
        BH <- BH.cut(pv,q)
        
        RBH[l,k] <- BH$nr
        if( RBH[l,k]==0)
          {
            FDPBH[l,k] <- 0
            FNPBH[l,k] <- epsilon
          }else{
            FDPBH[l,k] <- sum(BH$re<(n*(1-epsilon)))/(RBH[l,k]+(RBH[l,k]==0))
            FNPBH[l,k] <- (n*epsilon-sum(BH$re>(n*(1-epsilon))))/(n-RBH[l,k])
          }
        
        
    
        
        if(k%%50==0)
          print(paste("k=",k,sep="")) 
      }
      
      print(paste("l=",l,sep=""))
    }
    
    FDR=rowMeans(FDP*R)/ ( rowMeans(R) +(rowMeans(R)==0 ) ) ;  FNR=rowMeans( (n-R)*FNP)/( rowMeans(n-R) + (rowMeans(n-R)==0) );  AveRej=rowMeans(R)
    FDRSunCai=rowMeans(RSunCai*FDPSunCai)/( rowMeans(RSunCai) + (rowMeans(RSunCai)==0) );  FNRSunCai=rowMeans((n-RSunCai)*FNPSunCai)/( rowMeans( n-RSunCai) + (rowMeans( n-RSunCai)==0));  AveRejSunCai=rowMeans(RSunCai)
    FDRBH=rowMeans(RBH*FDPBH)/( rowMeans(RBH) + (rowMeans(RBH)==0) );  FNRBH=rowMeans((n-RBH)*FNPBH)/( rowMeans(n-RBH) + (rowMeans(n-RBH)==0 ));  AveRejBH=rowMeans(RBH)
    FDRIdeal=rowMeans(RIdeal*FDPIdeal)/( rowMeans(RIdeal) + (rowMeans(RIdeal)==0) ); FNRIdeal=rowMeans((n-RIdeal)*FNPIdeal)/( rowMeans(n-RIdeal) + (rowMeans(n-RIdeal)==0) ); AveRejIdeal=rowMeans(RIdeal)
    FDRlocfdr=rowMeans(Rlocfdr*FDPlocfdr)/( rowMeans(Rlocfdr) + (rowMeans(Rlocfdr)==0) ); FNRlocfdr=rowMeans((n-Rlocfdr)*FNPlocfdr)/( rowMeans(n-Rlocfdr) + (rowMeans(n-Rlocfdr)==0) ); AveRejlocfdr=rowMeans(Rlocfdr)
    FDRlocfdr.em=rowMeans(Rlocfdr.em*FDPlocfdr.em)/( rowMeans(Rlocfdr.em) + (rowMeans(Rlocfdr.em)==0) ); FNRlocfdr.em=rowMeans((n-Rlocfdr.em)*FNPlocfdr.em)/( rowMeans(n-Rlocfdr.em) + (rowMeans(n-Rlocfdr.em)==0 ) ); AveRejlocfdr.em=rowMeans(Rlocfdr.em)
    FDRoracle=rowMeans(Roracle*FDPoracle)/( rowMeans(Roracle) + (rowMeans(Roracle)==0 ) ); FNRoracle=rowMeans((n-Roracle)*FNPoracle)/( rowMeans(n-Roracle) + (rowMeans(n-Roracle)==0 )); AveRejoracle=rowMeans(Roracle)
    
    
    y=list(q=q,n=n,epsilon=epsilon,NONzero=NONzero,mu=mu,beta=beta, sigma=sigma, FDR=FDR,FNR=FNR,FDRSunCai=FDRSunCai,FNRSunCai=FNRSunCai,AveRej=AveRej,AveRejSunCai=AveRejSunCai,FDRBH=FDRBH,FNRBH=FNRBH,AveRejBH=AveRejBH,FDP=FDP,FNP=FNP,FDPSunCai=FDPSunCai,FNPSunCai=FNPSunCai,FDPBH=FDPBH,FNPBH=FNPBH,R=R,RSunCai=RSunCai,RBH=RBH,FDRIdeal=FDRIdeal, FNRIdeal=FNRIdeal, AveRejIdeal=AveRejIdeal,FDPIdeal=FDPIdeal, FNPIdeal=FNPIdeal, FDPlocfdr=FDPlocfdr, FNPlocfdr=FNPlocfdr, Rlocfdr=Rlocfdr, FDRlocfdr=FDRlocfdr, FNRlocfdr=FNRlocfdr, AveRejlocfdr=AveRejlocfdr, MaxX=MaxX, NoFalseSunCai=NoFalseSunCai, FDPlocfdr.em=FDPlocfdr.em, FNPlocfdr.em=FNPlocfdr.em, Rlocfdr.em=Rlocfdr.em, FDRlocfdr.em=FDRlocfdr.em, FNRlocfdr.em=FNRlocfdr.em, AveRejlocfdr.em = AveRejlocfdr.em, FDPoracle=FDPoracle, FNPoracle=FNPoracle, Rorcale=Roracle, FDRoracle=FDRoracle, FNRoracle=FNRoracle, AveRejoracle=AveRejoracle, EM.esti=EM.esti)
    
    
    
    ## case <- "Dep_Normal"
    ## file=paste("./result/Dep_Normal_n_",n,"_sigma_",sigma,"_sigma_2_",sigma2,"_q_",q,"_beta_",beta,"_p1_",p1,".Rdata",sep="")
    ## save(y,file=file)
    
    ## fig.fdr <- paste("../../../tex/TEST/Revision2/fig_FDR_",case,"_n_",n,"_sigma_",10*sigma,"_sigma2_",10*sigma2, "_q_",10*q,"_beta_",10*beta,"_p1_",10*p1,".pdf", sep="")
    ## fig.fnr <- paste("../../../tex/TEST/Revision2/fig_FNR_",case,"_n_",n,"_sigma_",10*sigma,"_sigma2_",10*sigma2, "_q_",10*q,"_beta_",10*beta,"_p1_",10*p1,".pdf", sep="")
    ## fig.power <- paste("../../../tex/TEST/Revision2/fig_power_",case,"_n_",n,"_sigma_",10*sigma,"_sigma2_",10*sigma2, "_q_",10*q,"_beta_",10*beta,"_p1_",10*p1,".pdf", sep="")
    

    LWD <- 1.5
    CEX <- 1.5
    
    x <- y$mu
    
    ## mFDR plot
    ## pdf(file=fig.fdr)
    
    plot(x, y$FDRBH, type="b", col="grey", pch=25, ylim=c(0,max( y$FDR, y$FDRBH, y$FDRlocfdr, y$FDRlocfdr.em, y$FDRSunCai, y$FDRoracle)*1.5), xlab=TeX("$\\mu$"), ylab="", main="mFDR", cex=CEX, lwd=LWD )
    lines(x, y$FDRlocfdr, type="b", col="purple", pch=13, cex=CEX, lwd=LWD)
    lines(x, y$FDRSunCai, type="b", col='green', pch=5, cex=CEX, lwd=LWD )
    lines(x, y$FDRlocfdr.em, type="b", col="cyan", pch=0, cex=CEX, lwd=LWD)
    lines(x, y$FDR, type="b", col='red', pch=16, cex=CEX, lwd=LWD)
    lines(x, 0*x + q, type="l", col="black", pch=2, lwd=LWD )
    legend("topleft",c("BH","Lfdr-locfdr","Lfdr-SC","Lfdr-EM","CLAT","q-level"),
           col=c("grey","purple","green","cyan","red","black"),
           pch=c(25,13,5,0, 16,NA),cex=CEX,pt.cex=CEX,lty=1,lwd=LWD) 
    ## dev.off()
    
    ## mFNR plot
    ## pdf(file=fig.fnr)
    
    plot(x, y$FNRBH, type="b", col="grey", pch=25, ylim=c(0,max( y$FNR, y$FNRBH, y$FNRlocfdr, y$FNRlocfdr.em, y$FNRSunCai, y$FNRoracle)), xlab=TeX("$\\mu$"), ylab="", main="FNR", cex=CEX, lwd=LWD )
    lines(x, y$FNRlocfdr, type="b", col="purple", pch=13, cex=CEX, lwd=LWD)
    lines(x, y$FNRSunCai, type="b", col='green', pch=5, cex=CEX, lwd=LWD )
    lines(x, y$FNRlocfdr.em, type="b", col="cyan", pch=0, cex=CEX, lwd=LWD)
    lines(x, y$FNR, type="b", col='red', pch=16, cex=CEX, lwd=LWD)
    legend(x=0.75,y=0.02,c("BH","Lfdr-locfdr","Lfdr-SC","Lfdr-EM","CLAT"),         col=c("grey","purple","green","cyan","red"),         pch=c(25,13,5,0,16),cex=CEX,pt.cex=CEX,lty=1,lwd=LWD) 
    
    ## dev.off()
    
    
    ## Proportion of true rejection
    ## pdf(file=fig.power)
    
    plot(x, apply( y$RBH*(1-y$FDPBH), 1, mean)/y$NONzero, type="b", col="grey", pch=25, ylim=c(0,1), xlab=TeX("$\\mu$"), ylab="", main="Proportion of True Rejection", cex=CEX, lwd=LWD )
    lines(x, apply( y$Rlocfdr*(1-y$FDPlocfdr), 1, mean)/y$NONzero, type="b", col="purple", pch=13, cex=CEX, lwd=LWD)
    lines(x, apply( y$RSunCai*(1-y$FDPSunCai), 1, mean)/y$NONzero, type="b", col='green', pch=5, cex=CEX, lwd=LWD )
lines(x, apply( y$Rlocfdr.em*(1-y$FDPlocfdr.em), 1, mean)/y$NONzero, type="b", col="cyan", pch=0, cex=CEX, lwd=LWD)
    lines(x, apply( y$R*(1-y$FDP), 1, mean)/y$NONzero, type="b", col='red', pch=16, cex=CEX, lwd=LWD)
    legend("topleft",c("BH","Lfdr-locfdr","Lfdr-SC","Lfdr-EM","CLAT"),
           col=c("grey","purple","green","cyan","red"),
       pch=c(25,13,5,0,16),cex=CEX,pt.cex=CEX,lty=1, lwd=LWD) 

    ## dev.off()

  }
