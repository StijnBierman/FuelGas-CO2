
# test file
rm(list=ls())
source("C:/Apps/FuelGas-CO2/R/helpers.R")

ExtraPars <- list(K = 24*365*4,
                  mu_aux = 2,
                  sd_aux = 0.2)
AllRes <- NULL
for ( N in c(20,40,60,80,100) ) {
  for ( rho_flow in c(0,0.9) ) {
    input <- list(nsims = 10000,
                  n = N,
                  mu_EF = 2.5,
                  sd_EF = 0.15,
                  mu_B = 10000,
                  sd_B = 2000,
                  rho = 0.9,
                  rho_B = rho_flow)
    
    TSimDat <- SimulateData(Sinput=input,SExtraPars=ExtraPars) # create synthetic data sets
    SimResults <- SimulateResults(Zinput=input,ZExtraPars=ExtraPars,doBoot=FALSE,nboots=500) # sample from synthetic data and estimate parameters
    Res <- TableResults(TSimRes=SimResults, doBoot=FALSE, Tinputs=input) # summarise results
    
    Res$n <- rep(input$n,nrow(Res))
    Res$rho_B <- rep(input$rho_B,nrow(Res))
    Res$rho <- rep(input$rho,nrow(Res))
    Res$mu_EF <- rep(input$mu_EF,nrow(Res))
    Res$mu_B <- rep(input$mu_B,nrow(Res))
    Res$sd_EF <- rep(input$sd_EF,nrow(Res))
    Res$sd_B <- rep(input$sd_B,nrow(Res))
    Res$nsims <- rep(input$nsims,nrow(Res))
    
    AllRes <- rbind(AllRes,Res)
  } # for ( rho_flow in c(0,0.9) )
  print(N)
} # for ( N in c(20,40,60,80,100) )

rm(Res,TSimDat,SimResults,input,N,rho_flow)

###################
# visualise
graphics.off()

for ( rhoB in c(0,0.9) ) {
  png(filename = paste0("C:/Apps/FuelGas-CO2/testOutput/rhoB_",rhoB,"test.png"),
      width = 600, height = 700, pointsize = 14)
  mypch <- c(0,1,2,15)
  mycol <- c("black","red","blue",grey(0.3))
  estimators <- unique(AllRes$Estimator)
  
  layout(mat = matrix(c(1,2,3), ncol=1), widths = c(1,1,1), heights = c(1,4,4))
  layout.show(3)
  par(mar=c(0,0,0,0))
  par(oma=c(0,0,0,0))
  plot(x = 1, y = 1, type = "n", xlim = c(0,1), ylim = c(0,1), xlab = "", ylab = "", axes = FALSE)
  text(x = 0.1, y = 0.75, labels = bquote(mu[EF]==.(unique(AllRes$mu_EF))~"; "~
                                            sigma[EF]==.(unique(AllRes$sd_EF))~"; "~
                                            mu[B]==.(unique(AllRes$mu_B))~"; "~
                                            sigma[B]==.(unique(AllRes$sd_B))~"; "~
                                            rho==.(unique(AllRes$rho))~"; "~
                                            rho[B]==.(unique(AllRes$rho_B))),
       pos = 4, cex = 1.5
  )
  xoff <- 0.01
  points(x = 0.05, y = 0.25, pch = mypch[1], col = mycol[1], cex = 2.5)
  points(x = 0.20, y = 0.25, pch = mypch[2], col = mycol[2], cex = 2.5)
  points(x = 0.55, y = 0.25, pch = mypch[3], col = mycol[3], cex = 2.5)
  points(x = 0.80, y = 0.25, pch = mypch[4], col = mycol[4], cex = 2.5)
  text(x = 0.05+xoff, y = 0.25, pos = 4, labels = estimators[1], col = mycol[1], cex = 1.5)
  text(x = 0.20+xoff, y = 0.25, pos = 4, labels = estimators[2], col = mycol[2], cex = 1.5)
  text(x = 0.55+xoff, y = 0.25, pos = 4, labels = estimators[3], col = mycol[3], cex = 1.5)
  text(x = 0.80+xoff, y = 0.25, pos = 4, labels = estimators[4], col = mycol[4], cex = 1.5)
  box()
  #############
  par(mar=c(1,5,1,1))
  
  offs <- c(-4,-2,2,4)
  plot(x = 1, y = 1, type = "n", 
       xlim = range(AllRes$n+offs,AllRes$n+offs), 
       ylim = range(AllRes$Coveragelower,AllRes$Coverageupper,100),
       xlab = "", ylab = "", axes = FALSE)
  #mtext(text = paste("n=",unique(AllRes$n)), at = unique(AllRes$n),
  #      side = 1, line = 0, cex = 1.25)
  axis(2, cex.axis = 1.5); mtext(text = "coverage (%)", side = 2, line = 3, cex = 1.25)
  abline(h = 95, lty = 2)
  abline(v=unique(AllRes$n)+10)
  cc <- 1
  for ( estim in unique(AllRes$Estimator) ) {
    SubS <- which(AllRes$Estimator==estim & AllRes$rho_B==rhoB)
    points(x = AllRes$n[SubS] + offs[cc], y = AllRes$Coverage[SubS], 
           pch = mypch[cc], col = mycol[cc], cex = 2)
    segments(x0 = AllRes$n[SubS] + offs[cc], y0 = AllRes$Coveragelower[SubS],
             x1 = AllRes$n[SubS] + offs[cc], y1 = AllRes$Coverageupper[SubS], pch = mypch[cc], col = mycol[cc])
    cc <- cc + 1
  } # for ( estim in unique(AllRes$Estimator) )

  plot(x = 1, y = 1, type = "n", 
       xlim = range(AllRes$n+offs,AllRes$n+offs), 
       ylim = range(0,AllRes$RelUlower,AllRes$RelUupper),
       xlab = "", ylab = "", axes = FALSE)
  mtext(text = paste("n=",unique(AllRes$n)), at = unique(AllRes$n),
        side = 1, line = 0, cex = 1.25)
  axis(2, cex.axis = 1.5); mtext(text = "uncertainty (%)", side = 2, line = 3, cex = 1.25)
  abline(h = 0.5, lty = 2)
  abline(v=unique(AllRes$n)+10)
  cc <- 1
  for ( estim in unique(AllRes$Estimator) ) {
    SubS <- which(AllRes$Estimator==estim & AllRes$rho_B==rhoB)
    points(x = AllRes$n[SubS] + offs[cc], y = AllRes$RelU[SubS], 
           pch = mypch[cc], col = mycol[cc], cex = 2)
    segments(x0 = AllRes$n[SubS] + offs[cc], y0 = AllRes$RelUlower[SubS],
             x1 = AllRes$n[SubS] + offs[cc], y1 = AllRes$RelUupper[SubS], pch = mypch[cc], col = mycol[cc])
    cc <- cc + 1
  } # for ( estim in unique(AllRes$Estimator) )
  dev.off()
} # for ( rhoB in c(0,0.9) )
