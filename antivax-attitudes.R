# Yeah, the standard Hadley stack. /sigh
library(readxl)
library(tidyr)
library(dplyr)
library(ggplot2)

# Generates warnings, but only for the Ps who didn't do day 2
expData <- read_excel("Vacc_HPHH_publicDataset.xlsx", sheet = 2)

# Just exclude Ps who didn't do day 2 and failed the attention checks
expData.clean <- expData %>%
  filter(!is.na(StartDate_day2),
         `AttentionCheck_PostTest (if = 4 then include)` == 4,
         `AttentionChecks_Sum(include if = 4)` == 4,
         Paid_Attention == 1)

