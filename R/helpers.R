rm(list=ls())
#library(mvtnorm)
############
# simulate data
DoRun <- FALSE
if ( DoRun ) {
  input <- list(nsims = 25,
                n = 100,
                mu_EF = 2.5,
                sd_EF = 0.15,
                mu_B = 500,
                sd_B = 25,
                rho = 0.96,
                rho_B = 0)
  
  ExtraPars <- list(K = 24*365*4,
                    mu_aux = 2,
                    sd_aux = 0.2)
  #
  TSimDat <- SimulateData(Sinput=input,SExtraPars=ExtraPars)
  FitAuxModel(FSimDat=TSimDat,Finput=input,doGraph=TRUE,returnAuxLm=FALSE)
  
  RAuxLm <- FitAuxModel(FSimDat=TSimDat,Finput=input,doGraph=FALSE,returnAuxLm=TRUE)
  Results <- ApplyEstimatorsApp(QAuxLm=RAuxLm,QSimDat=TSimDat,Qinput=input,QExtraPars=ExtraPars,doBoot=TRUE,nboots=250)
  Results
  #
  #to <- Sys.time()
  SimResults <- SimulateResults(Zinput=input,ZExtraPars=ExtraPars,doBoot=TRUE,nboots=500)
  #Sys.time()-to
  #SimResults
  100*SimResults$QCoverage/input$nsims
  TableResults(TSimRes=SimResults, doBoot=TRUE, Tinputs=input)#
  ##
  #par(mfrow=c(1,1))
  #plotResults(VSimResults=SimResults, plotCoverage=TRUE, doBoot = TRUE,Vinputs=input)
  #
}

SimulateData <- function(Sinput=input,SExtraPars=ExtraPars) {
  Crho <- Sinput$rho*Sinput$rho_B 
  covmat <- matrix(data=c(Sinput$sd_EF^2,Sinput$rho*Sinput$sd_EF*SExtraPars$sd_aux,Sinput$rho_B*Sinput$sd_EF*Sinput$sd_B,
                          Sinput$rho*Sinput$sd_EF*SExtraPars$sd_aux,SExtraPars$sd_aux^2,Crho*SExtraPars$sd_aux*Sinput$sd_B,
                          Sinput$rho_B*Sinput$sd_EF*Sinput$sd_B,Crho*SExtraPars$sd_aux*Sinput$sd_B,Sinput$sd_B^2),nrow=3,ncol=3)
  #browser()
  RandDat <- t(chol(covmat)) %*% matrix(rnorm(n=SExtraPars$K*3),ncol=SExtraPars$K,nrow=3)
  #RandDat <- rmvnorm(n=K, mean=c(mu_EF, mu_aux, mu_B), sigma=covmat)
  y <- Sinput$mu_EF + RandDat[1,]
  x <- SExtraPars$mu_aux + RandDat[2,]
  B <- Sinput$mu_B + RandDat[3,]
 
  y[y<=0] <- 0.01
  x[x<=0] <- 0.01
  B[B<=0] <- 0.01
  
  #plot(x,y, xlab="auxiliary variable", ylab="EF")
  #browser()
  #lm_check <- lm(y~x)
  #sd_lm_truth <- sd(residuals(lm_check))
  #coef_truth <- coefficients(lm_check)
  
  #abline(lm_check, col="red")
  #hist(residuals(lm_check), main="residuals")
  
  # 4. Truth
  truth <- sum(B*y/sum(B))
  x_bar <- sum(B*x/sum(B)) # true mean of auxiliary variable over period of interest
  #x_bar <- mean(x)
  ###########################
  # sample 
  II <- sample(x = 1:SExtraPars$K, size = Sinput$n, replace = FALSE)
  y_samp <- y[II]
  x_samp <- x[II]
  B_samp <- B[II]

  return(list(truth=truth,x_bar=x_bar, B=B, x=x, y=y, x_samp=x_samp, y_samp=y_samp, B_samp=B_samp))  
}

#SimDat <- SimulateData(mu_EF=mu_EF,sd_EF=sd_EF,mu_B=mu_B,sd_B=sd_B,n=n,K=K,M=M,mu_aux=mu_aux,sd_aux=sd_aux, Srho=0.94)

FitAuxModel <- function(FSimDat=SimDat,Finput=input,doGraph=TRUE,returnAuxLm=FALSE) {
  lm_samp<- lm(FSimDat$y_samp ~ FSimDat$x_samp)
  sd_lm_samp <- sqrt((1/(Finput$n - 2))*sum(residuals(lm_samp)^2))
  coef_samp <- coefficients(lm_samp)
  
  if ( doGraph ) {
    par(mfrow=c(2,3))
    par(mar=c(5,5,4,3))
    hist(FSimDat$y_samp, main="", 
         xlab=expression(paste("EF (tC",O[2],"/t)")), ylab="numbers of samples", cex.lab=1.5, cex.axis=1.5)
    mtext(side=3, text=bquote("EF process variability"))
    
    hist(FSimDat$x_samp, main="", 
         xlab=paste("auxiliary variable"), ylab="numbers of samples", cex.lab=1.5, cex.axis=1.5)
    mtext(side=3, text=bquote("auxiliary variable process variability (-)"))
    
    mtext(side=3, text="Illustration (single realisation) of scenario", line=2, cex=1.5)
    
    hist(FSimDat$B_samp, main="", 
         xlab=paste("Flow (t/day)"), ylab="numbers of samples", cex.lab=1.5, cex.axis=1.5)
    mtext(side=3, text=bquote("Fuel gas flow process variability"))

    plot(FSimDat$B_samp, FSimDat$x_samp, xlab="Flow (t/day)", ylab="auxiliary variable (-)",
         cex.axis=1.5, cex.lab=1.5)

    plot(FSimDat$B_samp, FSimDat$y_samp, xlab="Flow (t/day)", ylab=expression(paste("EF (tC",O[2],"/t)")),
         cex.axis=1.5, cex.lab=1.5)
    
    plot(FSimDat$y_samp ~ FSimDat$x_samp, xlab="auxiliary variable (-)", ylab=expression(paste("EF (tC",O[2],"/t)")),
         cex.axis=1.5, cex.lab=1.5)
    abline(lm_samp, col="blue", lwd=2)
    qrho <- round(cor(x = FSimDat$y_samp, y = FSimDat$x_samp, method = "pearson"),3)
    
    mtext(side=3, text=bquote("correlation coefficient" ~ rho==.(qrho) ~ " sample size"~n==.(Finput$n)))
    box(which = "outer")
  }
  if (returnAuxLm ) { return(list(a=coef_samp[1],b=coef_samp[2],sd_lm_samp=sd_lm_samp)) }
}

#AuxLm <- FitAuxModel(SimDat=SimDat,mu_EF=mu_EF,sd_EF=sd_EF,mu_B=mu_B,sd_B=sd_B,n=n,doGraph=FALSE,returnAuxLm=TRUE)
#AuxLm <- FitAuxModel(SimDat=SimDat,mu_EF=mu_EF,sd_EF=sd_EF,mu_B=mu_B,sd_B=sd_B,n=n,doGraph=TRUE,returnAuxLm=FALSE)

#Res <- ApplyEstimatorsApp(AuxLm=AuxLm,SimDat=SimDat,mu_EF=mu_EF,sd_EF=sd_EF,mu_B=mu_B,sd_B=sd_B,n=n,K=K)


#############################
# Estimators

# estimators for in App only

ApplyEstimatorsApp <- function(QAuxLm=AuxLm,QSimDat=SimDat,Qinput=input,QExtraPars=ExtraPars,doBoot=TRUE,nboots=250) {
  
  y_samp_mean <- mean(QSimDat$y_samp)
  x_samp_mean <- mean(QSimDat$x_samp)
  
  y_samp_bar <- sum(QSimDat$y_samp*QSimDat$B_samp)/sum(QSimDat$B_samp)
  x_samp_bar <- sum(QSimDat$x_samp*QSimDat$B_samp)/sum(QSimDat$B_samp)
  
  w_samp <- QSimDat$B_samp/sum(QSimDat$B_samp)
  # 1a. Simple Random Sample (SRS)
  # as in NEA (Netherlands Emissions Authority) Excel sheet
  
  srs <- list(estimate = mean(QSimDat$y_samp),
               #se = sd(SimDat$y_samp)/sqrt(n),
               se = sd(QSimDat$y_samp)/sqrt(Qinput$n),
               t_val = qt(p = 0.975, df = Qinput$n-1))
  srs$Low <- srs$estimate-srs$t_val*srs$se
  srs$Up <- srs$estimate+srs$t_val*srs$se
  #srs$U <- abs(srs$Up-srs$Low) # old definition
  srs$U <- abs(srs$Up-srs$estimate) # new definition
  srs$relU <- 100*(srs$U/srs$estimate)
  srs$Contain <- srs$Low < QSimDat$truth & srs$Up > QSimDat$truth
  
  # 1b. Simple Random Sample (SRS)
  # weighted by throughput
  srsW <- list(estimate = sum(w_samp*QSimDat$y_samp),
               #se = sd(SimDat$y_samp)/sqrt(n),
               se = sqrt( (Qinput$n/(Qinput$n-1)) *sum(w_samp^2*(QSimDat$y_samp-y_samp_bar)^2)),
               t_val = qt(p = 0.975, df = Qinput$n-1))
  srsW$Low <- srsW$estimate-srsW$t_val*srsW$se
  srsW$Up <- srsW$estimate+srsW$t_val*srsW$se
  #srsW$U <- abs(srsW$Up-srsW$Low) # old definition
  srsW$U <- abs(srsW$Up-srsW$estimate) # new definition
  srsW$relU <- 100*(srsW$U/srsW$estimate)
  srsW$Contain <- srsW$Low < QSimDat$truth & srsW$Up > QSimDat$truth

  # 2b. regression estimator assuming x_bar is known
  # weighted by throughput
  lrA <- list(estimate = y_samp_bar + QAuxLm$b*(QSimDat$x_bar - x_samp_bar),
              se_bit1 =  sum((QSimDat$y_samp - y_samp_bar)^2),
              se_bit2 = sum((QSimDat$y_samp-y_samp_bar)*(QSimDat$x_samp-x_samp_bar))^2,
              se_bit3 = sum((QSimDat$x_samp-x_samp_bar)^2),
              t_val = qt(p = 0.975, df = Qinput$n-2))
  lrA$se <- sqrt( (1/(Qinput$n*(Qinput$n-2))) * (lrA$se_bit1 - (lrA$se_bit2/lrA$se_bit3)) )
  lrA$Low <- lrA$estimate-lrA$t_val*lrA$se
  lrA$Up <- lrA$estimate+lrA$t_val*lrA$se
  # lrA$U <- abs(lrA$Up-lrA$Low) # old definition
  lrA$U <- abs(lrA$Up-lrA$estimate) # new definition
  lrA$relU <- 100*(lrA$U/lrA$estimate)
  lrA$Contain <- lrA$Low < QSimDat$truth & lrA$Up > QSimDat$truth
  
  # 3a. Van Zanten
  # based on all the x measurements, so assuming auxiliary variable continuously monitored
  vz <- list(estimate = QAuxLm$a + QAuxLm$b*QSimDat$x_bar,
             se_bit1 = sum((QSimDat$x_samp - QSimDat$x_bar)^2)/sum((QSimDat$x_samp - x_samp_mean)^2),
             se_bit2 = 1 + (1/QExtraPars$K) * sum((QSimDat$B - mean(QSimDat$B))^2) / mean(QSimDat$B)^2,
             t_val = qt(p = 0.975, df = Qinput$n-2))
  #browser()
  vz$se <- QAuxLm$sd_lm_samp * sqrt((1/Qinput$n)*vz$se_bit1 + (1/QExtraPars$K)*vz$se_bit2)
  vz$Low <- vz$estimate-vz$t_val*vz$se
  vz$Up <- vz$estimate+vz$t_val*vz$se
  # vz$U <- abs(vz$Up-vz$Low) # old definition
  vz$U <- abs(vz$Up-vz$estimate) # new definition
  vz$relU <- 100*(vz$U/vz$estimate)
  vz$Contain <- vz$Low < QSimDat$truth & vz$Up > QSimDat$truth

  if (doBoot ) {
    # 4. Bootstrap SRS weighted mean
    # 5. Bootstrap regression estimator
    
    m_b <- mlr_b <- rep(NA,nboots)
    for (qq in 1:nboots) {
      BI <- sample(x = 1:Qinput$n, size = Qinput$n, replace = TRUE)
      b_y <- QSimDat$y_samp[BI]
      b_B <- QSimDat$B_samp[BI]
      b_x <- QSimDat$x_samp[BI]
      
      m_b[qq] <- sum(b_y*b_B)/sum(b_B)
      
      if ( abs(Qinput$rho_B)<=0.01 ) { b_lm <- lm(b_y~b_x) }
      if ( abs(Qinput$rho_B)>0.01 ) { b_lm <- lm(b_y~b_x+b_B) }
      b_sd <- sd(residuals(b_lm))
      b_co <- coefficients(b_lm)
      if ( abs(Qinput$rho_B)<=0.01 ) { 
        t_p <- b_co[1] + b_co[2]*QSimDat$x + rnorm(n = QExtraPars$K, mean = 0, sd = b_sd)
        } # if ( abs(Qinput$rho_B)<=0.01 )
      if ( abs(Qinput$rho_B)>0.01 ) { 
        t_p <- b_co[1] + b_co[2]*QSimDat$x  + b_co[3]*QSimDat$B + rnorm(n = QExtraPars$K, mean = 0, sd = b_sd)
      } # if ( abs(Qinput$rho_B)>0.01 )
      
      mlr_b[qq] <- sum(QSimDat$B*t_p)/sum(QSimDat$B)
    } # for (qq in 1:nboots)
    bo <- list(estimate=quantile(x = m_b,probs = 0.5),
               Low=quantile(x = m_b,probs = 0.025),
               Up=quantile(x = m_b,probs = 0.975))
    bo$se <- sd(m_b)
    #bo$U <- abs(bo$Up-bo$Low) # old definition
    bo$U <- abs(bo$Up-bo$estimate) # old definition
    bo$relU <- 100*(bo$U/bo$estimate)
    bo$Contain <- bo$Low < QSimDat$truth & bo$Up > QSimDat$truth
    
    bolr <- list(#estimate=quantile(x = mlr_b,probs = 0.5),
               estimate=mean(mlr_b),
               Low=quantile(x = mlr_b,probs = 0.025),
               Up=quantile(x = mlr_b,probs = 0.975))
    bolr$se <- sd(mlr_b)
    #bolr$U <- abs(bolr$Up-bolr$Low)
    bolr$U <- abs(bolr$Up-bolr$estimate)
    bolr$relU <- 100*(bolr$U/bolr$estimate)
    bolr$Contain <- bolr$Low < QSimDat$truth & bolr$Up > QSimDat$truth
  } # if (doBoot )
  #hist(mlr_b,50)
  #abline(v=QSimDat$truth, col="blue")
  #abline(v=bolr$Low, col="green")
  #abline(v=bolr$Up, col="green")
  
  #browser()
  ##
  if (!doBoot ) {
    Results <- data.frame(NEAE=srs$estimate,NEAse=srs$se,NEArelU=srs$relU,NEAContain=srs$Contain,
                          srsE=srsW$estimate,srsse=srsW$se,srsrelU=srsW$relU,srsContain=srsW$Contain,
                          lrE=lrA$estimate,lrse=lrA$se,lrrelU=lrA$relU,lrContain=lrA$Contain,
                          vzE=vz$estimate,vzse=vz$se,vzrelU=vz$relU,vzContain=vz$Contain)
  } # if (!doBoot )
  if (doBoot ) {
    Results <- data.frame(NEAE=srs$estimate,NEAse=srs$se,NEArelU=srs$relU,NEAContain=srs$Contain,
                          srsE=srsW$estimate,srsse=srsW$se,srsrelU=srsW$relU,srsContain=srsW$Contain,
                          lrE=lrA$estimate,lrse=lrA$se,lrrelU=lrA$relU,lrContain=lrA$Contain,
                          vzE=vz$estimate,vzse=vz$se,vzrelU=vz$relU,vzContain=vz$Contain,
                          boE=bo$estimate,bose=bo$se,borelU=bo$relU,boContain=bo$Contain,
                          bolrE=bolr$estimate,bolrse=bolr$se,bolrrelU=bolr$relU,bolrContain=bolr$Contain)
  } # if (doBoot )
  
  return(Results)
} # ApplyEstimatorsApp
#Results <- ApplyEstimatorsApp(AuxLm=AuxLm,SimDat=SimDat,mu_EF=mu_EF,sd_EF=sd_EF,mu_B=mu_B,sd_B=sd_B,n=n)

SimulateResults <- function(Zinput=input,ZExtraPars=ExtraPars,doBoot=TRUE,nboots=250) {
  
  Uncertainties <- StdErrMean <- matrix(data=NA,nrow=Zinput$nsims,ncol=6)
  Coverage <- rep(0,6)
  for ( i in 1:Zinput$nsims ) {
    
    PSimDat <- SimulateData(Sinput=Zinput,SExtraPars=ZExtraPars)
    PAuxLm <- FitAuxModel(FSimDat=PSimDat,Finput=Zinput,doGraph=FALSE,returnAuxLm=TRUE)
    PResults <- ApplyEstimatorsApp(QAuxLm=PAuxLm,QSimDat=PSimDat,Qinput=Zinput,QExtraPars=ZExtraPars,doBoot=doBoot,nboots=nboots)
    
    #browser()
    Uncertainties[i,1:4] <- c(PResults$NEArelU,PResults$srsrelU,PResults$vzrelU,PResults$lrrelU)
    StdErrMean[i,1:4] <- c(PResults$NEAse,PResults$srsse,PResults$vzse,PResults$lrse)
    Coverage[1] <- Coverage[1] + PResults$NEAContain
    Coverage[2] <- Coverage[2] + PResults$srsContain
    Coverage[3] <- Coverage[3] + PResults$vzContain
    Coverage[4] <- Coverage[4] + PResults$lrContain
    if ( doBoot ){ 
      Uncertainties[i,5] <- PResults$borelU 
      Uncertainties[i,6] <- PResults$bolrrelU
      StdErrMean[i,5] <- PResults$bose 
      StdErrMean[i,6] <- PResults$bolrse
      Coverage[5] <- Coverage[5] + PResults$boContain 
      Coverage[6] <- Coverage[6] + PResults$bolrContain 
    } # if ( doBoot )
  } # for ( i in 1:Zinput$nsims )
  PCoverage <- Coverage#/nsims
  PNEA <- quantile(x = Uncertainties[,1],probs = c(0.025,0.5,0.975))
  Psrs <- quantile(x = Uncertainties[,2],probs = c(0.025,0.5,0.975))
  Pvz <- quantile(x = Uncertainties[,3],probs = c(0.025,0.5,0.975))
  Plr <- quantile(x = Uncertainties[,4],probs = c(0.025,0.5,0.975))

  seNEA <- quantile(x = StdErrMean[,1],probs = c(0.025,0.5,0.975))
  sesrs <- quantile(x = StdErrMean[,2],probs = c(0.025,0.5,0.975))
  sevz <- quantile(x = StdErrMean[,3],probs = c(0.025,0.5,0.975))
  selr <- quantile(x = StdErrMean[,4],probs = c(0.025,0.5,0.975))
  
  if ( !doBoot ){ 
    return(list(QNEA=PNEA,Qsrs=Psrs,Qvz=Pvz,Qlr=Plr,
                seNEA=seNEA,sesrs=sesrs,sevz=sevz,selr=selr,
                QCoverage=PCoverage))
  } # if ( !doBoot )
  if ( doBoot ){ 
    Pbo <- quantile(x = Uncertainties[,5],probs = c(0.025,0.5,0.975))
    Pbolr <- quantile(x = Uncertainties[,6],probs = c(0.025,0.5,0.975))
    sebo <- quantile(x = StdErrMean[,5],probs = c(0.025,0.5,0.975))
    sebolr <- quantile(x = StdErrMean[,6],probs = c(0.025,0.5,0.975))
    return(list(QNEA=PNEA,Qsrs=Psrs,Qvz=Pvz,Qlr=Plr,Qbo=Pbo,Qbolr=Pbolr,
                seNEA=seNEA,sesrs=sesrs,sevz=sevz,selr=selr,sebo=sebo,sebolr=sebolr,
                QCoverage=PCoverage))
  } # if ( doBoot )
  
} # SimulateResults

#SimResults <- SimulateResults(mu_EF=mu_EF,sd_EF=sd_EF,mu_B=mu_B,sd_B=sd_B,n=n,K=K,M=M,mu_aux=mu_aux,sd_aux=sd_aux, Srho=0.94, nsims=50)

TableResults <- function(TSimRes=SimResults, doBoot=FALSE, Tinputs=input) {

  if ( doBoot ) {
    TabDat <- data.frame(Estimator=c("SRS","SRS flow weighted","Bootstrap flow weighted","Cochran","van Zanten","Bootstrap auxiliary"),
                         SEM=c(TSimRes$seNEA[2],TSimRes$sesrs[2],TSimRes$sebo[2],TSimRes$selr[2],TSimRes$sevz[2],TSimRes$sebolr[2]),
                         SEMlower=c(TSimRes$seNEA[1],TSimRes$sesrs[1],TSimRes$sebo[1],TSimRes$selr[1],TSimRes$sevz[1],TSimRes$sebolr[1]),
                         SEMupper=c(TSimRes$seNEA[3],TSimRes$sesrs[3],TSimRes$sebo[3],TSimRes$selr[3],TSimRes$sevz[3],TSimRes$sebolr[3]),
                         RelU=c(TSimRes$QNEA[2],TSimRes$Qsrs[2],TSimRes$Qbo[2],TSimRes$Qlr[2],TSimRes$Qvz[2],TSimRes$Qbolr[2]),
                         RelUlower=c(TSimRes$QNEA[1],TSimRes$Qsrs[1],TSimRes$Qbo[1],TSimRes$Qlr[1],TSimRes$Qvz[1],TSimRes$Qbolr[1]),
                         RelUupper=c(TSimRes$QNEA[3],TSimRes$Qsrs[3],TSimRes$Qbo[3],TSimRes$Qlr[3],TSimRes$Qvz[3],TSimRes$Qbolr[3]),
                         Coverage=round(100*TSimRes$QCoverage[c(1,2,5,4,3,6)]/Tinputs$nsims,2),
                         Coveragelower=round(100*qbeta(p=0.025,TSimRes$QCoverage[c(1,2,5,4,3,6)],Tinputs$nsims-TSimRes$QCoverage[c(1,2,5,4,3,6)]),2),
                         Coverageupper=round(100*qbeta(p=0.975,TSimRes$QCoverage[c(1,2,5,4,3,6)],Tinputs$nsims-TSimRes$QCoverage[c(1,2,5,4,3,6)]),2)  )
  } # if ( doBoot )
  if ( !doBoot ) {
    TabDat <- data.frame(Estimator=c("SRS","SRS flow weighted","Cochran","van Zanten"),
                         SEM=c(TSimRes$seNEA[2],TSimRes$sesrs[2],TSimRes$selr[2],TSimRes$sevz[2]),
                         SEMlower=c(TSimRes$seNEA[1],TSimRes$sesrs[1],TSimRes$selr[1],TSimRes$sevz[1]),
                         SEMupper=c(TSimRes$seNEA[3],TSimRes$sesrs[3],TSimRes$selr[3],TSimRes$sevz[3]),
                         RelU=c(TSimRes$QNEA[2],TSimRes$Qsrs[2],TSimRes$Qlr[2],TSimRes$Qvz[2]),
                         RelUlower=c(TSimRes$QNEA[1],TSimRes$Qsrs[1],TSimRes$Qlr[1],TSimRes$Qvz[1]),
                         RelUupper=c(TSimRes$QNEA[3],TSimRes$Qsrs[3],TSimRes$Qlr[3],TSimRes$Qvz[3]),
                         Coverage=round(100*TSimRes$QCoverage[c(1,2,4,3)]/Tinputs$nsims,2),
                         Coveragelower=round(100*qbeta(p=0.025,TSimRes$QCoverage[c(1,2,4,3)],Tinputs$nsims-TSimRes$QCoverage[c(1,2,4,3)]),2),
                         Coverageupper=round(100*qbeta(p=0.975,TSimRes$QCoverage[c(1,2,4,3)],Tinputs$nsims-TSimRes$QCoverage[c(1,2,4,3)]),2)  )
  } # if ( !doBoot )
  return(TabDat)
  
} # TableResults

plotResults <- function(VSimResults, plotCoverage=FALSE, doBoot=FALSE, Vinputs=inputs) {
  par(mar=c(0,0,2,0))
  plot(1,1,type="n",xlim=c(0,10),ylim=c(0,4),axes=FALSE,xlab="",ylab="")
  
  ####
  # Estimator
  text(x = 1, y = 3.5, labels = "estimator", pos=4, cex=1.25)
  segments(x0=1, x1=9, y0=3.1, y1=3.1, col="grey")
  text(x = 1, y = 3.0, labels ="based on lab measurements only:", pos=4)
  text(x = 1, y = 2.7, labels = "SRS no flow weighting", pos =4)
  text(x = 1, y = 2.5, labels = "SRS flow weighting", pos =4)
  if(doBoot) { text(x = 1, y = 2.3, labels = "Bootstrap SRS", pos =4) }
    
  segments(x0=1, x1=9, y0=2.1, y1=2.1, col="grey")
  text(x = 1, y = 2.0, labels ="based on lab measurements and auxiliary variable:", pos=4)
  text(x = 1, y = 1.7, labels = "van Zanten", pos =4)
  text(x = 1, y = 1.5, labels = "Cochran 1977", pos =4)
  if(doBoot) { text(x = 1, y = 1.3, labels = "Bootstrap regression", pos =4) }

  #####
  # Standard error of the Mean
  text(x = 3, y = 3.5, labels = "Std err of the mean (ton/ton)", pos =4, cex=1.25)
  text(x = 3, y = 2.7, labels = paste(round(VSimResults$Qsrs[1],2),"-",round(VSimResults$Qsrs[3],2)), pos =4)
  if(doBoot) {
    text(x = 3, y = 2.4, labels = paste(round(VSimResults$Qbo[1],2),"-",round(VSimResults$Qbo[3],2)), pos =4)
  } # if doBoot
  
  text(x = 3, y = 1.7, labels = paste(round(VSimResults$Qvz[1],2),"-",round(VSimResults$Qvz[3],2)), pos =4)
  text(x = 3, y = 1.4, labels = paste(round(VSimResults$Qlr[1],2),"-",round(VSimResults$Qlr[3],2)), pos =4)
  if(doBoot) {
    text(x = 3, y = 1.1, labels = paste(round(VSimResults$Qbolr[1],2),"-",round(VSimResults$Qbolr[3],2)), pos =4)
  } # if doBoot
  
  #####
  # Relative Uncertainty
  text(x = 3, y = 3.5, labels = "Relative Uncertainty (%)", pos =4, cex=1.5)
  text(x = 3, y = 3.3, labels = "target = 0.5%", pos =4, cex=0.75, col="blue")
  text(x = 3, y = 2.7, labels = paste(round(VSimResults$Qsrs[1],2),"-",round(VSimResults$Qsrs[3],2)), pos =4)
  if(doBoot) {
    text(x = 3, y = 2.4, labels = paste(round(VSimResults$Qbo[1],2),"-",round(VSimResults$Qbo[3],2)), pos =4)
  } # if doBoot
  
  text(x = 3, y = 1.7, labels = paste(round(VSimResults$Qvz[1],2),"-",round(VSimResults$Qvz[3],2)), pos =4)
  text(x = 3, y = 1.4, labels = paste(round(VSimResults$Qlr[1],2),"-",round(VSimResults$Qlr[3],2)), pos =4)
  if(doBoot) {
    text(x = 3, y = 1.1, labels = paste(round(VSimResults$Qbolr[1],2),"-",round(VSimResults$Qbolr[3],2)), pos =4)
  } # if doBoot
  
  #####
  # Coverage
  if ( plotCoverage ) {
    text(x = 7, y = 3.5, labels = "Coverage probability (%)", pos =4, cex=1.5)
    text(x = 7, y = 3.3, labels = "target = 95%", pos =4, cex=0.75, col="blue")
    #text(x = 7, y = 2.7, labels = round(SimResults$QCoverage[1],3)*100, pos =4)
    text(x = 7, y = 2.7, labels = paste(round(100*qbeta(p=0.025,VSimResults$QCoverage[1],Vinputs$nsims-VSimResults$QCoverage[1]),2),
                                        "-",
                                        round(100*qbeta(p=0.975,VSimResults$QCoverage[1],Vinputs$nsims-VSimResults$QCoverage[1]),2)), pos =4)
    
    if(doBoot) {
      #text(x = 7, y = 2.4, labels = round(SimResults$QCoverage[4],3)*100, pos =4)
      text(x = 7, y = 2.4, labels = paste(round(100*qbeta(p=0.025,VSimResults$QCoverage[4],Vinputs$nsims-VSimResults$QCoverage[4]),2),
                                          "-",
                                          round(100*qbeta(p=0.975,VSimResults$QCoverage[4],Vinputs$nsims-VSimResults$QCoverage[4]),2)), pos =4)
      
    } # if doBoot
    #text(x = 7, y = 1.7, labels = round(SimResults$QCoverage[2],3)*100, pos =4)
    text(x = 7, y = 1.7, labels = paste(round(100*qbeta(p=0.025,VSimResults$QCoverage[2],Vinputs$nsims-VSimResults$QCoverage[2]),2),
                                        "-",
                                        round(100*qbeta(p=0.975,VSimResults$QCoverage[2],Vinputs$nsims-VSimResults$QCoverage[2]),2)), pos =4)
    
    #text(x = 7, y = 1.4, labels = round(SimResults$QCoverage[3],3)*100, pos =4)
    text(x = 7, y = 1.4, labels = paste(round(100*qbeta(p=0.025,VSimResults$QCoverage[3],Vinputs$nsims-VSimResults$QCoverage[3]),2),
                                        "-",
                                        round(100*qbeta(p=0.975,VSimResults$QCoverage[3],Vinputs$nsims-VSimResults$QCoverage[3]),2)), pos =4)
    if(doBoot) {
      #text(x = 7, y = 1.1, labels = round(SimResults$QCoverage[5],3)*100, pos =4)
      text(x = 7, y = 1.1, labels = paste(round(100*qbeta(p=0.025,VSimResults$QCoverage[5],Vinputs$nsims-VSimResults$QCoverage[5]),2),
                                          "-",
                                          round(100*qbeta(p=0.975,VSimResults$QCoverage[5],Vinputs$nsims-VSimResults$QCoverage[5]),2)), pos =4)
    } # if doBoot
    
  } # if plotCoverage
  mtext(side = 3, text=paste("Results based on",Vinputs$nsims,"simulations: relative uncertainty and coverage of the 95% confidence of the estimated annual average EF"),
        line=0, cex=1.5)
  box(which = "outer")
  
}

#plotResults(SimResults, plotCoverage=TRUE, doBoot=TRUE)

############################
# full set of estimators
ApplyEstimators <- function(x_samp, y_samp, K, n, B_samp, truth, x_bar, M, Cheap=FALSE) {
  
  y_samp_mean <- mean(y_samp)
  x_samp_mean <- mean(x_samp)
  
  w_samp <- B_samp/sum(B_samp)
  
  # 1a. Simple Random Sample (SRS)
  # not weighted by throughput
  srs <- list(estimate = y_samp_mean,
              se = sd(y_samp)/sqrt(n),
              t_val = qt(p = 0.975, df = n-1) 
  )
  srs$Low <- srs$estimate-srs$t_val*srs$se
  srs$Up <- srs$estimate+srs$t_val*srs$se
  # srs$U <- abs(srs$Up-srs$Low) # old definition
  srs$U <- abs(srs$Up-srs$estimate) # new definition
  srs$relU <- 100*(srs$U/srs$estimate)
  srs$Contain <- srs$Low < truth & srs$Up > truth
  
  # 1b. Simple Random Sample (SRS)
  # weighted by throughput
  srsW <- list(estimate = sum(w_samp*y_samp),
               se = sd(y_samp)/sqrt(n),
               t_val = qt(p = 0.975, df = n-1) 
  )
  srsW$Low <- srsW$estimate-srsW$t_val*srsW$se
  srsW$Up <- srsW$estimate+srsW$t_val*srsW$se
  # srsW$U <- abs(srsW$Up-srsW$Low) # old definition
  srsW$U <- abs(srsW$Up-srsW$estimate) # old definition
  srsW$relU <- 100*(srsW$U/srsW$estimate)
  srsW$Contain <- srsW$Low < truth & srsW$Up > truth
  
  # 2a. regression estimator assuming x_bar is known
  # not weighted by throughput
  lr <- list(estimate = y_samp_mean + coef_samp[2]*(mean(x) - x_samp_mean),
             se = sqrt((1/(n*(n-2))) * sum( ((y_samp - y_samp_mean)-coef_samp[2]*(x_samp - x_samp_mean))^2 )),
             t_val = qt(p = 0.975, df = n-2) 
  )
  lr$Low <- lr$estimate-lr$t_val*lr$se
  lr$Up <- lr$estimate+lr$t_val*lr$se
  # lr$U <- abs(lr$Up-lr$Low) # old definition
  lr$U <- abs(lr$Up-lr$estimate) # new definition
  lr$relU <- 100*(lr$U/lr$estimate)
  lr$Contain <- lr$Low < truth & lr$Up > truth
  
  # 2b. regression estimator assuming x_bar is known
  # weighted by throughput
  lrA <- list(estimate = sum(w_samp*y_samp) + coef_samp[2]*(x_bar - sum(w_samp*x_samp)),
              se = sqrt((1/(n*(n-2))) * sum( ((y_samp - y_samp_mean)-coef_samp[2]*(x_samp - sum(w_samp*x_samp)))^2 )),
              t_val = qt(p = 0.975, df = n-2) 
  )
  lrA$Low <- lrA$estimate-lrA$t_val*lrA$se
  lrA$Up <- lrA$estimate+lrA$t_val*lrA$se
  # lrA$U <- abs(lrA$Up-lrA$Low) # old definition
  lrA$U <- abs(lrA$Up-lrA$estimate) # new definition
  lrA$relU <- 100*(lrA$U/lrA$estimate)
  lrA$Contain <- lrA$Low < truth & lrA$Up > truth
  
  # 3a. Van Zanten
  # based on all the x measurements, so assuming auxiliary variable continuously monitored
  vz <- list(estimate = coef_samp[1] + coef_samp[2]*x_bar,
             se_bit1 = sum((x_samp - x_bar)^2)/sum((x_samp - x_samp_mean)^2),
             se_bit2 = 1 + (1/K) * sum((B - mean(B))^2) / mean(B)^2,
             t_val = qt(p = 0.975, df = n-2)
  )
  
  vz$se <- sd_lm_samp * sqrt((1/n)*vz$se_bit1 + (1/K)*vz$se_bit2)
  vz$Low <- vz$estimate-vz$t_val*vz$se
  vz$Up <- vz$estimate+vz$t_val*vz$se
  # vz$U <- abs(vz$Up-vz$Low) # old definition
  vz$U <- abs(vz$Up-vz$estimate) # new definition
  vz$relU <- 100*(vz$U/vz$estimate)
  vz$Contain <- vz$Low < truth & vz$Up > truth
  
  # 3b. Van Zanten
  # based on a sample of auxiliary variables
  aux_I <- sample(x = 1:K, size = M, replace = FALSE)
  x_samp_aux <- x[aux_I]
  B_samp_aux <- B[aux_I]
  vz2 <- list(estimate = coef_samp[1] + coef_samp[2]*mean(x_samp_aux),
              se_bit1 = sum((x_samp - mean(x_samp_aux))^2)/sum((x_samp - x_samp_mean)^2),
              se_bit2 = 1 + (1/M) * sum((B_samp_aux - mean(B_samp_aux))^2) / mean(B_samp_aux)^2,
              t_val = qt(p = 0.975, df = n-2)
  )
  
  vz2$se <- sd_lm_samp * sqrt((1/n)*vz2$se_bit1 + (1/M)*vz2$se_bit2)
  vz2$Low <- vz2$estimate-vz2$t_val*vz2$se
  vz2$Up <- vz2$estimate+vz2$t_val*vz2$se
  # vz2$U <- abs(vz2$Up-vz2$Low) # old definition
  vz2$U <- abs(vz2$Up-vz2$estimate) # new definition
  vz2$relU <- 100*(vz2$U/vz2$estimate)
  vz2$Contain <- vz2$Low < truth & vz2$Up > truth
  
  if ( !Cheap ) {
    # 4a. Bootstrap and random draws from error model
    # based on all the x measurements, so assuming auxiliary variable continuously monitored
    bo <- BootEst(x_samp=x_samp, y_samp=y_samp, M=K, x=x, K=K, truth=truth, nboots=2000)
    
    # 4b. Boostrap with a random sample of auxiliary variables
    bob <- BootEst(x_samp=x_samp, y_samp=y_samp, M=M, x=x, K=K, truth=truth, nboots=2000)
    
    Results <- c(truth,
                 srs$estimate,srs$U,srs$relU,srs$Contain,
                 lr$estimate,lr$U,lr$relU,lr$Contain,
                 lrA$estimate,lrA$U,lrA$relU,lrA$Contain,
                 vz$estimate,vz$U,vz$relU,vz$Contain,
                 vz2$estimate,vz2$U,vz2$relU,vz2$Contain,
                 bo$estimate,bo$U,bo$relU,bo$Contain,
                 bob$estimate,bob$U,bob$relU,bob$Contain)
  }  
  
  if ( Cheap ) {
    Results <- c(truth,
                 srs$estimate,srs$U,srs$relU,srs$Contain,
                 lr$estimate,lr$U,lr$relU,lr$Contain,
                 lrA$estimate,lrA$U,lrA$relU,lrA$Contain,
                 vz$estimate,vz$U,vz$relU,vz$Contain,
                 vz2$estimate,vz2$U,vz2$relU,vz2$Contain)
  }
  
  return(Results)
}




###########
# Bootstrap estimate of uncertainty
BootEst <- function(x_samp, y_samp, M, x, K, truth, nboots) {
  
  y_simu_est <- rep(NA,nboots)
  for ( i in 1:nboots ) {
    bootI <- sample(x = 1:length(x_samp), size = length(x_samp), replace = TRUE)
    x_samp_boot <- x_samp[bootI]
    y_samp_boot <- y_samp[bootI]
    
    if ( M < K ) { x_aux <- sample(x = x, size = M, replace = TRUE) }
    if ( M == K ) { x_aux <- x }
    
    lm_boot <- lm(y_samp_boot ~ x_samp_boot)
    sd_lm_boot <- sd(residuals(lm_boot))
    coef_boot <- coefficients(lm_boot)
    
    y_simu <- coef_boot[1] + coef_boot[2]*x_aux + rnorm(n = M, mean = 0, sd = sd_lm_boot)
    y_simu_est[i] <- mean(y_simu)
  }
  
  result <- quantile(x = y_simu_est, probs = c(0.025,0.5,0.975))
  Res <- list(estimate = result[2],
              Lo = result[1],
              Up = result[3],
              #U = abs(result[3]-result[1]), # old definition
              U = abs(result[3]-result[2]), # new definition
              relU = 100*U/result[2],
              Contain = truth > result[1] & truth < result[3])
  return(Res)
}
