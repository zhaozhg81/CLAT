                                        #' This function runs the simulation and plots the figure for case III.
                                        #'
                                        #' @param beta The value of beta, which is between 0 and 1. In the paper, beta is chosen as 0.3 and 0.4 respectively.
                                        #' @export


caseIII <- function(beta)
  {
    numSim <- 2
    
    n <- 5000
    l <- 1.2
    q <- 0.1
    alpha.LB <- beta + log( (1-q)/q * l)/log(n)
    
    NONzero=floor(n^{1-beta})
    epsilon <- NONzero/n;
    
    
    alphas <- alpha.LB + c(0:40)/400
    
    rNum <- length(alphas)
    
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
    
    
    
    for (lll in 1:length(alphas)){
      
      for (k in 1:numSim){
        
        ## ########################################################
        ## #####################################################
        ## ## Generate the random numbers
        alpha <- alphas[lll]
        
        U <- runif( n, 0, 1)
        X <- U
        trueInd <- array(n,c(1,n))
        trueInd[1:(n-NONzero)] <- 0
        trueInd[(n-NONzero+1):n] <- 1
        X[1:(n-NONzero)] <- U[1:(n-NONzero)]
        X[(n-NONzero+1):n] <- ( l*n^{-alpha}*sqrt(2*U[(n-NONzero+1):n]) ) * ( U[(n-NONzero+1):n]<0.5) +  ( 2*l*n^{-alpha} - l*n^{-alpha}*sqrt( 1- 2*(U[(n-NONzero+1):n]-0.5) ) )*(U[(n-NONzero+1):n]>0.5)
        
        
        f0x <- array(1, n)
        f1x <- n^{2*alpha}*X/l/l*(X <= (l*n^{-alpha})) + ( -n^{2*alpha}*(X-2*l*n^{-alpha}) )/l/l * ( X>( l*n^{-alpha}) ) *( X <= (2*l*n^{-alpha}) ) 
    
        ## ################################################################
        ## #########################################################
        ## ##   Benjamini Hochberg procedure
        pv <- X
        BH <- BH.cut(pv,q)
        
        RBH[lll,k] <- BH$nr
        if( RBH[lll, k] !=0 ){
          FDPBH[lll,k] <- sum(BH$re<(n*(1-epsilon)))/(RBH[lll,k]+(RBH[lll,k]==0))
          FNPBH[lll,k] <- (n*epsilon-sum(BH$re>(n*(1-epsilon))))/(n-RBH[lll,k])
        }else
        {
          FDPBH[lll,k] <- 0
          FNPBH[lll,k] <- epsilon
        }
        
        
        
        
        ## ###########################################################
        ## ###############################################################
        ## ## Sun and Cai 's approach
        adaptiveZ<-adaptZ.func(X, q, nullDist='Unif')
        ## the threshold
        threshold<-adaptiveZ$th
        ## number of rejected hypotheses
        RSunCai[lll,k]<-adaptiveZ$nr
        ## the rejected hypotheses
        if( RSunCai[lll,k]!=0 ){
          rh<-adaptiveZ$re
          FDPSunCai[lll,k] <- sum(rh<(n*(1-epsilon)))/(RSunCai[lll,k]+(RSunCai[lll,k]==0))
          FNPSunCai[lll,k] <- (n*epsilon-sum(rh>(n*(1-epsilon))))/(n-RSunCai[lll,k])
        }else{
          FDPSunCai[lll,k] <- 0
          FNPSunCai[lll,k] <- epsilon
        }   
        
        
        ## ##########################################################
        ## ##########################################################
        ## ## CLAT
        
        ## # Positve side
        pv=X
        CLATres <- CLAT(pv,q)
        R[lll,k] <- CLATres$R
        FDP[lll,k] <- sum(CLATres$SigInd*(1-trueInd))/(R[lll,k]+(R[lll,k]==0))
        FNP[lll,k] <- sum((1-CLATres$SigInd)*trueInd)/(n-R[lll,k])
        
        
        
        ## ################################################################
        ## #########################################################
        ##    oracle
        LFDR.oracle <- (1-epsilon)*f0x/( (1-epsilon)*f0x + epsilon * f1x )
        LFDR.ord <- order(LFDR.oracle,decreasing=FALSE)
        
        LFDR.sort <- sort(LFDR.oracle, decreasing=FALSE)
        Rtmp <- max( c(1:n) * ( ( 1/c(1:n) * cumsum( LFDR.sort) )< q ) )
        
        Roracle[lll,k] <- Rtmp
        
        if( Rtmp> 0){
          Indoracle <- array(0, n)
          Indoracle[ LFDR.ord[1:Rtmp] ] <- 1
          FDPoracle[lll,k] <- sum(Indoracle*(1-trueInd[1,]))/(Roracle[lll,k]+(Roracle[lll,k]==0));
          FNPoracle[lll,k] <- sum((1-Indoracle)*trueInd[1,])/(n-Roracle[lll,k]);
        }else{
          FDPoracle[lll,k] <- 0
          FNPoracle[lll,k] <- sum(trueInd)/n;
        }
        
        ## ################################################################
        ## #########################################################
        ##    localfdr (Efron)
        LFDR <- locfdr( qnorm(X, lower.tail=FALSE), plot=0 );
        LFDR.sort <- sort(LFDR$fdr, decreasing=FALSE)
        R.lfdr <- max( c(1:n) * ( ( 1/c(1:n) * cumsum( LFDR.sort) )< q ) )
        
        Rlocfdr[lll,k] <- R.lfdr
        
        if( R.lfdr> 0){
          Indlocfdr <- (LFDR$fdr <= LFDR.sort[R.lfdr]);
          FDPlocfdr[lll,k] <- sum(Indlocfdr*(1-trueInd))/(Rlocfdr[lll,k]+(Rlocfdr[lll,k]==0));
          FNPlocfdr[lll,k] <- sum((1-Indlocfdr)*trueInd)/(n-Rlocfdr[lll,k]);
        }else{
          FDPlocfdr[lll,k] <- 0
          FNPlocfdr[lll,k] <- sum(trueInd)/n;
        }
        
        ## ################################################################
        ## #########################################################
        ## ##   EM Algorithm for estimating the true parameter
        lfdr.em <- LFDR.EM( qnorm(X, lower.tail=FALSE), L=1, DELTA=0.001)
        lfdr.em.sort <- sort(lfdr.em, decreasing=FALSE)
        R.lfdr.em <- max( c(1:n) * ( ( 1/c(1:n) * cumsum( lfdr.em.sort) )< q ) )
        
        Rlocfdr.em[lll,k] <- R.lfdr.em
        
        if( R.lfdr.em> 0){
          Indlocfdr.em <- (lfdr.em <= lfdr.em.sort[R.lfdr.em]);
          FDPlocfdr.em[lll,k] <- sum( Indlocfdr.em*(1-trueInd[1,]))/(Rlocfdr.em[lll,k]+(Rlocfdr.em[lll,k]==0));
          FNPlocfdr.em[lll,k] <- sum((1-Indlocfdr.em)*trueInd[1,])/(n-Rlocfdr.em[lll,k]);
        }else{
          FDPlocfdr.em[lll,k] <- 0
          FNPlocfdr.em[lll,k] <- sum(trueInd)/n;
        }          
        
        
        if(k%%50==0)
          print(paste("k=",k,sep="")) 
      }
      
      print(paste("l=",lll,sep=""))
    }
    
    
    FDR=rowMeans(FDP*R)/ ( rowMeans(R) +(rowMeans(R)==0 ) ) ;  FNR=rowMeans( (n-R)*FNP)/( rowMeans(n-R) + (rowMeans(n-R)==0) );  AveRej=rowMeans(R)
    FDRSunCai=rowMeans(RSunCai*FDPSunCai)/( rowMeans(RSunCai) + (rowMeans(RSunCai)==0) );  FNRSunCai=rowMeans((n-RSunCai)*FNPSunCai)/( rowMeans( n-RSunCai) + (rowMeans( n-RSunCai)==0));  AveRejSunCai=rowMeans(RSunCai)
    FDRBH=rowMeans(RBH*FDPBH)/( rowMeans(RBH) + (rowMeans(RBH)==0) );  FNRBH=rowMeans((n-RBH)*FNPBH)/( rowMeans(n-RBH) + (rowMeans(n-RBH)==0 ));  AveRejBH=rowMeans(RBH)
    FDRIdeal=rowMeans(RIdeal*FDPIdeal)/( rowMeans(RIdeal) + (rowMeans(RIdeal)==0) ); FNRIdeal=rowMeans((n-RIdeal)*FNPIdeal)/( rowMeans(n-RIdeal) + (rowMeans(n-RIdeal)==0) ); AveRejIdeal=rowMeans(RIdeal)
    FDRlocfdr=rowMeans(Rlocfdr*FDPlocfdr)/( rowMeans(Rlocfdr) + (rowMeans(Rlocfdr)==0) ); FNRlocfdr=rowMeans((n-Rlocfdr)*FNPlocfdr)/( rowMeans(n-Rlocfdr) + (rowMeans(n-Rlocfdr)==0) ); AveRejlocfdr=rowMeans(Rlocfdr)
    FDRlocfdr.em=rowMeans(Rlocfdr.em*FDPlocfdr.em)/( rowMeans(Rlocfdr.em) + (rowMeans(Rlocfdr.em)==0) ); FNRlocfdr.em=rowMeans((n-Rlocfdr.em)*FNPlocfdr.em)/( rowMeans(n-Rlocfdr.em) + (rowMeans(n-Rlocfdr.em)==0 ) ); AveRejlocfdr.em=rowMeans(Rlocfdr.em)
    FDRoracle=rowMeans(Roracle*FDPoracle)/( rowMeans(Roracle) + (rowMeans(Roracle)==0 ) ); FNRoracle=rowMeans((n-Roracle)*FNPoracle)/( rowMeans(n-Roracle) + (rowMeans(n-Roracle)==0 )); AveRejoracle=rowMeans(Roracle)
    
    y=list(q=q,n=n,epsilon=epsilon,NONzero=NONzero,alphas=alphas,l=l,beta=beta, FDR=FDR,FNR=FNR,FDRSunCai=FDRSunCai,FNRSunCai=FNRSunCai,AveRej=AveRej,AveRejSunCai=AveRejSunCai,FDRBH=FDRBH,FNRBH=FNRBH,AveRejBH=AveRejBH,FDP=FDP,FNP=FNP,FDPSunCai=FDPSunCai,FNPSunCai=FNPSunCai,FDPBH=FDPBH,FNPBH=FNPBH,R=R,RSunCai=RSunCai,RBH=RBH,FDRIdeal=FDRIdeal, FNRIdeal=FNRIdeal, AveRejIdeal=AveRejIdeal,FDPIdeal=FDPIdeal, FNPIdeal=FNPIdeal, FDPlocfdr=FDPlocfdr, FNPlocfdr=FNPlocfdr, Rlocfdr=Rlocfdr, FDRlocfdr=FDRlocfdr, FNRlocfdr=FNRlocfdr, AveRejlocfdr=AveRejlocfdr, MaxX=MaxX, NoFalseSunCai=NoFalseSunCai, FDPlocfdr.em=FDPlocfdr.em, FNPlocfdr.em=FNPlocfdr.em, Rlocfdr.em=Rlocfdr.em, FDRlocfdr.em=FDRlocfdr.em, FNRlocfdr.em=FNRlocfdr.em, AveRejlocfdr.em = AveRejlocfdr.em, FDPoracle=FDPoracle, FNPoracle=FNPoracle, Rorcale=Roracle, FDRoracle=FDRoracle, FNRoracle=FNRoracle, AveRejoracle=AveRejoracle)
    
    
    ## file=paste("./result/unif_vary_alpha_n_",n,"_q_",q,"_beta_",beta,"_l_",l,".Rdata",sep="")
    ## save(y,file=file)
    
    ## case <- "Unif"
    ## fig.fdr <- paste("../../../tex/TEST/Revision/fig_FDR_unif_n_",n,"_q_",10*q,"_beta_",10*beta,"_l_",10*l,".pdf", sep="")
    ## fig.fnr <- paste("../../../tex/TEST/Revision/fig_FNR_unif_n_",n,"_q_",10*q,"_beta_",10*beta,"_l_",10*l,".pdf", sep="")
    ## fig.power <- paste("../../../tex/TEST/Revision/fig_power_unif_n_",n,"_q_",10*q,"_beta_",10*beta,"_l_",10*l,".pdf", sep="")
        
    x <- y$alphas
    LWD <- 1.5
    CEX <- 1.5
    
    
    ## mFDR plot
    ## pdf(file=fig.fdr)    
    
    plot(x, y$FDRBH, type="b", col="grey", pch=25, ylim=c(0,max( y$FDR, y$FDRBH, y$FDRlocfdr, y$FDRlocfdr.em, y$FDRSunCai, y$FDRoracle)*1.5), xlab=TeX("$\\alpha$"), ylab="", main="mFDR", cex=CEX, lwd=LWD )
    lines(x, y$FDRlocfdr, type="b", col="purple", pch=13, cex=CEX, lwd=LWD)
    lines(x, y$FDRSunCai, type="b", col='green', pch=5, cex=CEX, lwd=LWD )
    lines(x, y$FDRlocfdr.em, type="b", col="cyan", pch=0, cex=CEX, lwd=LWD)
    lines(x, y$FDRoracle, type="b", col='blue', pch=8, cex=CEX, lwd=LWD)
    lines(x, y$FDR, type="b", col='red', pch=16, cex=CEX, lwd=LWD)
    lines(x, 0*x + q, type="l", col="black", pch=2, lwd=LWD )
    
    legend("topleft",c("BH","Lfdr-locfdr","Lfdr-SC","Lfdr-EM","Lfdr-oracle","CLAT","q-level"),
           col=c("grey","purple","green","cyan","blue","red","black"),
       pch=c(25,13,5,0,8, 16,NA),cex=CEX,pt.cex=CEX,lty=1,lwd=LWD) 
    ## dev.off()
    
    ## mFNR plot
    ## pdf(file=fig.fnr)

    plot(x, y$FNRBH, type="b", col="grey", pch=25, ylim=c(0,max( y$FNR, y$FNRBH, y$FNRlocfdr, y$FNRlocfdr.em, y$FNRSunCai, y$FNRoracle)), xlab=TeX("$\\alpha$"), ylab="", main="mFNR", cex=CEX, lwd=LWD )
    lines(x, y$FNRlocfdr, type="b", col="purple", pch=13, cex=CEX, lwd=LWD)
lines(x, y$FNRSunCai, type="b", col='green', pch=5, cex=CEX, lwd=LWD )
    lines(x, y$FNRlocfdr.em, type="b", col="cyan", pch=0, cex=CEX, lwd=LWD)
    lines(x, y$FNRoracle, type="b", col='blue', pch=8, cex=CEX, lwd=LWD)
    lines(x, y$FNR, type="b", col='red', pch=16, cex=CEX, lwd=LWD)

    
    legend(x=0.75,y=0.02,c("BH","Lfdr-locfdr","Lfdr-SC","Lfdr-EM", "Lfdr-oracle","CLAT"),         col=c("grey","purple","green","cyan", "blue","red"),         pch=c(25,13,5,0,8,16),cex=CEX,pt.cex=CEX,lty=1,lwd=LWD) 
    
    ## dev.off()
    
    
    ## Proportion of true rejection
    ## pdf(file=fig.power)
    
    plot(x, apply( y$RBH*(1-y$FDPBH), 1, mean)/y$NONzero, type="b", col="grey", pch=25, ylim=c(0,1), xlab=TeX("$\\alpha$"), ylab="", main="Proportion of True Rejection", cex=CEX, lwd=LWD )
    lines(x, apply( y$Rlocfdr*(1-y$FDPlocfdr), 1, mean)/y$NONzero, type="b", col="purple", pch=13, cex=CEX, lwd=LWD)
    lines(x, apply( y$RSunCai*(1-y$FDPSunCai), 1, mean)/y$NONzero, type="b", col='green', pch=5, cex=CEX, lwd=LWD )
    lines(x, apply( y$Rlocfdr.em*(1-y$FDPlocfdr.em), 1, mean)/y$NONzero, type="b", col="cyan", pch=0, cex=CEX, lwd=LWD)
    lines(x, apply( y$Rorcale*(1-y$FDPoracle), 1, mean)/y$NONzero, type="b", col='blue', pch=8, cex=CEX, lwd=LWD)
    lines(x, apply( y$R*(1-y$FDP), 1, mean)/y$NONzero, type="b", col='red', pch=16, cex=CEX, lwd=LWD)
    
    legend("topleft",c("BH","Lfdr-locfdr","Lfdr-SC","Lfdr-EM","Lfdr-oracle","CLAT"),
           col=c("grey","purple","green","cyan", "blue","red"),
           pch=c(25,13,5,0,8,16),cex=CEX,pt.cex=CEX,lty=1, lwd=LWD) 
    ## dev.off()

  }
