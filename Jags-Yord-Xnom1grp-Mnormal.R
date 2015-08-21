# Jags-Yord-Xnom1grp-Mnormal.R 
# Accompanies the book:
#   Kruschke, J. K. (2015). Doing Bayesian Data Analysis, Second Edition: 
#   A Tutorial with R, JAGS, and Stan. Academic Press / Elsevier.

source("DBDA2E-utilities.R")

#===============================================================================

genMCMC = function( datFrm, yName , qName,
                    numSavedSteps=50000 , thinSteps = 1,
                    saveName=NULL ,
                    runjagsMethod=runjagsMethodDefault , 
                    nChains=nChainsDefault ) { 
  #-----------------------------------------------------------------------------
  # THE DATA.
  q = as.numeric(as.factor(datFrm[[qName]]))
  nQlevels = max(q)
  
  y = as.numeric(datFrm[[yName]])
  # Do some checking that data make sense:
  if ( any( y!=round(y) ) ) { stop("All y values must be integers (whole numbers).") }
  if ( any( y < 1 ) ) { stop("All y values must be 1 or larger.") }
  # COMPRESS OUT ANY EMPTY VALUES OF Y:
  yOrig=y
  y=as.numeric(factor(y,levels=names(table(y))))
  if ( any(y != yOrig) ) { 
    warning("*** WARNING: Y RE-CODED TO REMOVE EMPTY LEVELS ***")
  }
  Ntotal = length(y)
  # Threshold 1 and nYlevels-1 are fixed; other thresholds are estimated.
  # This allows all parameters to be interpretable on the response scale.
  nYlevels = max(y)  
  thresh = matrix(data = NA, nrow = nQlevels, ncol = nYlevels-1)
  thresh[, 1] = 1 + 0.5
  thresh[, nYlevels-1] = nYlevels-1 + 0.5
  # Specify the data in a list, for later shipment to JAGS:
  dataList = list(
    q = q,
    nQlevels = nQlevels,
    y = y,
    nYlevels = nYlevels,
    thresh = thresh,
    Ntotal = Ntotal 
  )
  #-----------------------------------------------------------------------------
  # THE MODEL.
  modelString = "
  model {
    for (i in 1:Ntotal) {
      y[i] ~ dcat(pr[i,1:nYlevels])
      pr[i,1] <- pnorm(thresh[q[i], 1] , mu[q[i]] , 1/sigma[q[i]]^2)
      for (k in 2:(nYlevels-1)) {
        pr[i,k] <- max( 0 ,  pnorm(thresh[q[i], k] ,   mu[q[i]] , 1/sigma[q[i]]^2 )
                           - pnorm(thresh[q[i], k-1] , mu[q[i]] , 1/sigma[q[i]]^2 ))
      }
      pr[i,nYlevels] <- 1 - pnorm(thresh[q[i], nYlevels-1] , mu[q[i]] , 1/sigma[q[i]]^2)
    }

    # Unique mu, sigma, and thresh vector for each question
    for (j in 1:nQlevels) {
      mu[j] ~ dnorm( (1+nYlevels)/2 , 1/(nYlevels)^2 )
      sigma[j] ~ dunif( nYlevels/1000 , nYlevels*10 )

      # Threshold distributions. 1 and nYlevels-1 are fixed, not stochastic
      for (k in 2:(nYlevels-2)) {  
        thresh[j, k] ~ dnorm( k+0.5 , 1/2^2 )
      }
    }
  }
  " # close quote for modelString
  # Write out modelString to a text file
  writeLines( modelString , con="TEMPmodel.txt" )
  #-----------------------------------------------------------------------------
  # INTIALIZE THE CHAINS.
  #   # Initial values of MCMC chains based on data:
  #   muInit = c( mean(y[x==1]) , mean(y[x==2]) )
  #   sigmaInit = c( sd(y[x==1]) , sd(y[x==2]) )
  #   threshInit = 1:(nYlevels-1)+0.5
  #   threshInit[1] = NA
  #   threshInit[nYlevels-1] = NA
  #   # Regarding initial values in next line: (1) sigma will tend to be too big if 
  #   # the data have outliers, and (2) nu starts at 5 as a moderate value. These
  #   # initial values keep the burn-in period moderate.
  #   initsList = list( mu=muInit, sigma=sigmaInit, nuMinusOne=4, thresh=threshInit )
  initsList = NULL
  #-----------------------------------------------------------------------------
  # RUN THE CHAINS
  parameters = c( "mu" , "sigma" , "thresh" )
  adaptSteps = 500               # Number of steps to "tune" the samplers
  burnInSteps = 1000
  runJagsOut <- run.jags( method=runjagsMethod ,
                          model="TEMPmodel.txt" , 
                          monitor=parameters , 
                          data=dataList ,  
                          #inits=initsList , 
                          n.chains=nChains ,
                          adapt=adaptSteps ,
                          burnin=burnInSteps , 
                          sample=ceiling(numSavedSteps/nChains) ,
                          thin=thinSteps ,
                          summarise=FALSE ,
                          plots=FALSE )
  codaSamples = as.mcmc.list( runJagsOut )
  # resulting codaSamples object has these indices: 
  #   codaSamples[[ chainIdx ]][ stepIdx , paramIdx ]
  if ( !is.null(saveName) ) {
    save( codaSamples , file=paste(saveName,"Mcmc.Rdata",sep="") )
  }
  return( codaSamples )
} # end function

#===============================================================================

plotMCMC = function( codaSamples , datFrm , yName , qName, compVal , #RopeEff=NULL ,
                     showCurve=FALSE , 
                     saveName=NULL , saveType="jpg" ) {
  #-----------------------------------------------------------------------------
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  chainLength = NROW( mcmcMat )
  q = as.numeric(as.factor(datFrm[[qName]]))
  nQlevels = max(q)
  
  #-----------------------------------------------------------------------------
  # Plots for each question
  for (i in 1:nQlevels) {
    mu = mcmcMat[, paste0("mu[", i, "]")]
    sigma = mcmcMat[, paste0("sigma[", i, "]")]
    
    # Set up window and layout:
    openGraph(width=6.0,height=5.0)
    layout( matrix( c(1,2,3,4,5,6) , nrow=3, byrow=TRUE ) )
    par( mar=c(3.5,3.5,2.5,0.5) , mgp=c(2.25,0.7,0) )
    
    # Compute limits for plots of data with posterior pred. distributions
    y = datFrm[[yName]][q == i]
    xLim = c( min(y)-0.5 , max(y)+0.5 )
    xBreaks = seq( xLim[1] , xLim[2] , 1 )  
    histInfo = hist(y,breaks=xBreaks,plot=FALSE)
    yMax = 1.2 * max( c( histInfo$density ) )
    xVec = seq( xLim[1] , xLim[2] , length=501 )
    
    #-----------------------------------------------------------------------------
    # 1. mu
    xlim = range( c( mu ) )
    histInfo = plotPost( mu , cex.lab = 1.75 ,
                         showCurve=showCurve , compVal=compVal ,
                         xlab=bquote(mu) , main=paste("Mean") , 
                         col="skyblue" )
  
    #-----------------------------------------------------------------------------
    # 2. Plot data y and smattering of posterior predictive curves:
    # Data histogram:
    histInfo = hist( y , prob=TRUE , xlim=xLim , ylim=c(0,yMax) , breaks=xBreaks,
                     col="pink" , border="white" , xlab="y" , ylab="" , 
                     yaxt="n" , cex.lab=1.5 , 
                     main=paste("Data w. Post. Pred.") )
    # Posterior predictive probabilities of outcomes:
    outProb=matrix(0,nrow=chainLength,ncol=max(y))
    for ( stepIdx in 1:chainLength ) {
      threshCumProb = ( pnorm( ( mcmcMat[ stepIdx , 
                                          paste("thresh[",i,",",1:(max(y)-1),"]",sep="") ]
                                 - mu[stepIdx] )
                               / sigma[stepIdx] ) )
      outProb[stepIdx,] = c(threshCumProb,1) - c(0,threshCumProb)
    }
    outHdi = apply( outProb , 2 , HDIofMCMC )
    outMean = apply( outProb , 2 , median , na.rm=TRUE )
    show(outMean)
    points( x=1:max(y) , y=outMean  , pch=19 , cex=2 , col="skyblue" )
    segments( x0=1:max(y) , y0=outHdi[1,] , 
              x1=1:max(y) , y1=outHdi[2,] , lwd=4 , col="skyblue" )
    # annotate N:
    text( max(xVec) , yMax , bquote(N==.(length(y))) , adj=c(1.1,1.1) )
    
    #-----------------------------------------------------------------------------
    # 3. sigma
    xlim=range( c( sigma ) )
    histInfo = plotPost( sigma ,  xlim=xlim , cex.lab = 1.75 ,
                         showCurve=showCurve ,
                         xlab=bquote(sigma) , 
                         main=paste("Std. Dev.") , 
                         col="skyblue"  )
    #-----------------------------------------------------------------------------
    # 4. effect size. 
    effectSize = ( mu - compVal ) / sigma
    histInfo = plotPost( effectSize , compVal=0.0 , # ROPE=RopeEff ,
                         showCurve=showCurve ,
                         xlab=bquote( ( mu - .(compVal) ) / sigma ) ,
                         cex.lab=1.75 , main="Effect Size" ,
                         col="skyblue" )
   #-----------------------------------------------------------------------------  
   # 5. Thresholds:
    threshCols = grep(paste0("^thresh\\[", i), colnames(mcmcMat), value=TRUE)
    threshMean = rowMeans( mcmcMat[,threshCols] )
    xLim = range(mcmcMat[,threshCols])
    nPtToPlot = 500
    plotIdx = floor(seq(1,nrow(mcmcMat),length=nPtToPlot))
    plot( mcmcMat[plotIdx,threshCols[1]] , threshMean[plotIdx] , col="skyblue" ,
          xlim=xLim , xlab="Threshold" , ylab="Mean Threshold" )
    abline(v=mean(mcmcMat[plotIdx,threshCols[1]]),lty="dashed",col="skyblue")
    for ( jj in 2:length(threshCols) ) {
      points( mcmcMat[plotIdx,threshCols[jj]] , threshMean[plotIdx] , col="skyblue" )
      abline(v=mean(mcmcMat[plotIdx,threshCols[jj]]),lty="dashed",col="skyblue")
    }
    #-----------------------------------------------------------------------------
    if ( !is.null(saveName) ) {
      saveGraph( file=paste(saveName,"Q",i,"Post",sep=""), type=saveType)
    }
  }
}

#===============================================================================
