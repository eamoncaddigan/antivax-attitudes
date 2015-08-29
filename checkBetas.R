# Just a hack script to check that I got the sum-to-zero stuff right.

betaLevels <- c(5, 3, 2)

# Hey, how do you initialize an empty n-dimentional matrix in R?
mFromBeta <- rep(NA, prod(betaLevels))
dim(mFromBeta) <- betaLevels
mFromMCMC <- rep(NA, prod(betaLevels))
dim(mFromMCMC) <- betaLevels

for (j1 in 1:betaLevels[1]) {
  for (j2 in 1:betaLevels[2]) {
    for (j3 in 1:betaLevels[3]) {
      # Oh God.
      mFromBeta[j1,j2,j3] <- mean(mcmcMat[,"b0"] + 
                                    mcmcMat[,paste0("b1[",j1,"]")] + 
                                    mcmcMat[,paste0("b2[",j2,"]")] + 
                                    mcmcMat[,paste0("b3[",j3,"]")] + 
                                    mcmcMat[,paste0("b1b2[",j1,",",j2,"]")] + 
                                    mcmcMat[,paste0("b1b3[",j1,",",j3,"]")] + 
                                    mcmcMat[,paste0("b2b3[",j2,",",j3,"]")] + 
                                    mcmcMat[,paste0("b1b2b3[",j1,",",j2,",",j3,"]")])
      mFromMCMC[j1,j2,j3] <- mean(mcmcMat[,paste0("m[",j1,",",j2,",",j3,"]")])
    }
  }
}
