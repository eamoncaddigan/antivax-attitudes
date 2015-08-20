# Yeah, the standard Hadley stack. /sigh
library(readxl)
library(tidyr)
library(dplyr)
library(ggplot2)

# Generates warnings, but only for the Ps who didn't do day 2
expData <- read_excel("Vacc_HPHH_publicDataset.xlsx", sheet = 2)
