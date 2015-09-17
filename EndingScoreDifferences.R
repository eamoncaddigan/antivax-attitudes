# Worried about failures of random assignment. Let's check differences in FINAL 
# scores, instead of differences in CHANGES. 
# 
# This assumes that all the chunks from antivax-attitudes.Rmd had been run in
# the current environment.

png(file="bayesian_ending_scores.png", width=900, height=300)
par(mfrow = c(1, 3), mar=c(2, 1, 1, 1), oma=c(0, 0, 4, 0))

# We'll look at all the pairs of levels for x2. Easy for three levels, but this
# code is a bit more generic.
x2Combn <- combn(length(levels(questionnaireData$intervention)), 2)

for (combnIdx in seq_len(dim(x2Combn)[2])) {
  x2Level1 <- x2Combn[2, combnIdx]
  x2Level2 <- x2Combn[1, combnIdx]
  plotPost((mcmcMat[, paste0("b2[", x2Level1, "]")] + mcmcMat[, paste0("b2b3[", x2Level1, ",2]")]) - 
               (mcmcMat[, paste0("b2[", x2Level2, "]")] + mcmcMat[, paste0("b2b3[", x2Level2, ",2]")]),
             main = "",
             compVal = 0.0, ROPE = c(-0.05, 0.05),
             xlab = "")
  mtext(paste0(levels(questionnaireData$intervention)[x2Level1], " vs.\n",
               levels(questionnaireData$intervention)[x2Level2]),
        side = 3, line = -2)
}
title("Ending score differences", outer=TRUE)
dev.off()
