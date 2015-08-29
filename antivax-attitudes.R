# Yeah, the standard Hadley stack. /sigh
library(readxl)
library(tidyr)
library(dplyr)
library(ggplot2)


# Gather and clean the data -----------------------------------------------

# Generates warnings, but only for the Ps who didn't do day 2
expData <- read_excel("Vacc_HPHH_publicDataset.xlsx", sheet = 2)

# Just a bit of data cleaning
expData.clean <- expData %>%
  # Add a subject number
  mutate(subject_number = 1:nrow(.)) %>%
  # Just exclude Ps who didn't do day 2 and failed the attention checks
  filter(Returned == 1,
         `AttentionCheck_PostTest (if = 4 then include)` == 4,
         `AttentionChecks_Sum(include if = 4)` == 4,
         Paid_Attention == 1)

# Get all the dependent measures into a DF
questionnaireData <- expData.clean %>%
  # pull out the columns and use BETTER NAMES (jeez Zach)
  select(subject_number,
         intervention = Condition,
         pretest.healthy = Healthy_VaxscalePretest,
         posttest.healthy = Healthy_VaxscalePosttest,
         pretest.diseases = Diseases_VaxScalePretest,
         posttest.diseases = Diseases_VaxScalePosttest,
         pretest.doctors = Doctors_VaxScalePreTest,
         posttest.doctors = Doctors_VaxScalePostTest,
         pretest.side_effects = Sideeffects_VaxScalePreTest,
         posttest.side_effects = Sideeffects_VaxScalePostTest,
         pretest.plan_to = Planto_VaxScalePreTest,
         posttest.plan_to = Planto_VaxScalePostTest) %>%
  # reverse-code the approrpiate columns
  mutate(pretest.diseases = 7 - pretest.diseases,
         posttest.diseases = 7 - posttest.diseases,
         pretest.side_effects = 7 - pretest.side_effects,
         posttest.side_effects = 7 - posttest.side_effects) %>%
  # "tidy" the data
  gather("question", "response", -subject_number, -intervention) %>%
  separate(question, c("interval", "question"), sep = "\\.") %>% 
  mutate(interval = factor(interval, c("pretest", "posttest"), ordered = TRUE),
         question = factor(question))


# Some plots --------------------------------------------------------------

# Check out the distribution of responses before and after the intervention
p1 <- ggplot(questionnaireData, aes(x = question, y = response, fill = interval)) +
  geom_violin() + 
  facet_grid(intervention ~ .)
print(p1)

# Look at each subject's change for each question
p2 <- ggplot(questionnaireData, aes(x = interval, y = response, group = subject_number)) + 
  geom_line(alpha = 0.2, position = position_jitter(w = 0.15, h = 0.15)) + 
  facet_grid(intervention ~ question)
print(p2)


# Bayesian analysis of survey data ----------------------------------------

# Fit a model to each question using pre-intervention data. 
modelData <- questionnaireData

source("Jags-Yord-Xnom1grp-Mnormal.R")
fileNameRoot = "antivax-mcmc"

mcmcCoda <- genMCMC(datFrm = modelData,
                    yName = "response",
                    x1Name = "question",
                    x2Name = "intervention",
                    x3Name = "interval",
                    numSavedSteps = 15000,
                    thinSteps = 10,
                    saveName = fileNameRoot)

# Display diagnostics of chain, for specified parameters:
# betas, sigmas
parameterNames <- varnames(mcmcCoda)
parameterNames <- parameterNames[grepl("^b", parameterNames) | grepl("^sigma", parameterNames)]
for (parName in parameterNames) {
  diagMCMC(codaObject = mcmcCoda,
           parName = parName, 
           saveName = fileNameRoot,
           saveType = "png")
}

# Display posterior information:
plotMCMC(mcmcCoda,
         datFrm = modelData,
         yName = "response",
         qName = "question",
         compVal = 3.5, 
         saveName = fileNameRoot,
         saveType = "png")

# This right here is the good stuff:
mcmcMat <- as.matrix(mcmcCoda)
plotPost((mcmcMat[, "b2b3[3,2]"] - mcmcMat[, "b2b3[3,1]"]) - 
           (mcmcMat[, "b2b3[1,2]"] - mcmcMat[, "b2b3[1,1]"]))
# It shows that post-pre for the "disease risk" is greater than post-pre for
# "autism correction", supporting the authors' findings.
