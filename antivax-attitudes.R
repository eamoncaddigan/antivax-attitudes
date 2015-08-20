# Yeah, the standard Hadley stack. /sigh
library(readxl)
library(tidyr)
library(dplyr)
library(ggplot2)

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
         pretest_healthy = Healthy_VaxscalePretest,
         posttest_healthy = Healthy_VaxscalePosttest,
         pretest_diseases = Diseases_VaxScalePretest,
         posttest_diseases = Diseases_VaxScalePosttest,
         pretest_doctors = Doctors_VaxScalePreTest,
         posttest_doctors = Doctors_VaxScalePostTest,
         pretest_sideeffects = Sideeffects_VaxScalePreTest,
         posttest_sideeffects = Sideeffects_VaxScalePostTest,
         pretest_planto = Planto_VaxScalePreTest,
         posttest_planto = Planto_VaxScalePostTest,
         pretest_autism = Autism_PreTest,
         posttest_autism = AutismAttitude_PostTest) %>%
  # reverse-code the approrpiate columns
  mutate(pretest_diseases = 7 - pretest_diseases,
         posttest_diseases = 7 - posttest_diseases,
         pretest_sideeffects = 7 - pretest_sideeffects,
         posttest_sideeffects = 7 - posttest_sideeffects) %>%
  # "tidy" the data
  gather("question", "response", -subject_number, -intervention) %>%
  separate(question, c("interval", "question"))

# Let's take a look at the data
ggplot(questionnaireData, aes(x = question, y = response, fill = interval)) +
  geom_violin() + 
  facet_grid(intervention ~ .)
