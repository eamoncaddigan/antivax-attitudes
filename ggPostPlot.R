# ggplot2-based plots for some of the things

library(ggplot2)
library(dplyr)

# Find the responses and the posterior predictions for the given factor levels
calcPosteriorPredictive <- function(modelData, codaObject,
                                    x1Level, x2Level = NA, x3Level = NA) {
  mcmcMat <- as.matrix(codaObject)
  
  # Separate the observed counts for the given levels
  # XXX - hard-coded factor names because laziness
  cellData <- filter(modelData, as.numeric(as.factor(question)) == x1Level)
  if (!is.na(x2Level)) {
    cellData <- filter(cellData, as.numeric(as.factor(intervention)) == x2Level)
  }
  if (!is.na(x3Level)) {
    cellData <- filter(cellData, as.numeric(as.factor(interval)) == x3Level)
  }
  y <- cellData[["response"]]
  
  # Find the mean, standard deviation, and thresholds
  mu <- mcmcMat[, "b0"] + mcmcMat[, paste0("b1[", x1Level, "]")]
  if (!is.na(x2Level)) {
    mu <- mu + mcmcMat[, paste0("b2[", x2Level, "]")] + 
      mcmcMat[, paste0("b1b2[", x1Level, ",", x2Level, "]")]
  }
  if (!is.na(x3Level)) {
    mu <- mu + mcmcMat[, paste0("b3[", x3Level, "]")] + 
      mcmcMat[, paste0("b1b3[", x1Level, ",", x3Level, "]")]
    if (!is.na(x2Level)) {
      mu <- mu + mcmcMat[, paste0("b1b2b3[", x1Level, ",", x2Level, ",", x3Level, "]")]
    }
  }
  sigma <- mcmcMat[, paste0("sigma[", x1Level, "]")]
  
  # Find the HDI and mean of the posterior probabilities of each of the y levels
  # (Stolen from Kruschke)
  chainLength <- nrow(mcmcMat)
  outProb <- matrix(0, nrow=chainLength, ncol=max(y))
  for (stepIdx in 1:chainLength) {
    threshCumProb <- pnorm((mcmcMat[stepIdx, 
                                    paste0("thresh[", x1Level, ",",1:(max(y)-1),"]")]
                            - mu[stepIdx]) / sigma[stepIdx])
    outProb[stepIdx,] <- c(threshCumProb, 1) - c(0, threshCumProb)
  }
  outHdi <- apply(outProb, 2, HDIofMCMC)
  outMean <- apply(outProb, 2, mean, na.rm=TRUE)
  
  # Create a data.frame with the posterior data
  plotData <- cellData %>%
    count(response) %>% 
    left_join(data_frame(response = 1:max(y)), ., by = "response") %>%
    mutate(n = ifelse(is.na(n), 0, n),
           response_proportion = n / sum(n),
           posterior_mean = outMean,
           posterior_hdi_low = outHdi[1,],
           posterior_hdi_high = outHdi[2,])
  return(plotData)
}


# Plot the responses and posterior predictions for the given factor levels
ggPosteriorPredictive <- function(modelData, codaObject, 
                                  x1Level, x2Level = NA, x3Level = NA) {
  plotData <- calcPosteriorPredictive(modelData, codaObject,
                                      x1Level, x2Level, x3Level)
  numLevels = length(unique(plotData[["response"]]))

  p <- ggplot(plotData, aes(x = response, y = response_proportion)) + 
    geom_bar(stat = "identity", binwidth=1, color="white", fill="skyblue") + 
    geom_point(aes(y = posterior_mean), 
               size=5, color="red") + 
    geom_errorbar(aes(ymin = posterior_hdi_low, ymax = posterior_hdi_high), 
                  size=3, color="red", width=0) + 
    scale_x_continuous(breaks = 1:numLevels) +
    scale_y_continuous(limits = c(0, 1)) +
    xlab("Response Level") + 
    ylab("Proportion")
  return(p)
}
