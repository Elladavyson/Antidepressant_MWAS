## Seeing if there are statistical differences in the demographic variables ##
## Using R version 4.4
library(dplyr)
library(tidyverse)
library(readr)
############ READING IN THE FILES ###############
# Read in the demographics file (copied from datastore in Antidep_methylation) for all of Generation Scotland 
setwd("/exports/eddie/scratch/s2112198/")
demographics <- read.table("GS_demograph.txt", header = T)
demographics <- demographics %>% select(-c(starts_with("FID"), "Batch"))
# Antidepressant exposure phenotypes 
sr <- read.csv("selfrep_pheno3_methyl_03_05.csv", header = T)
pd <- read.csv("antidep_pheno1_clean_appt.csv", header =T)
# Residual phenotypes (ones used in the MWAS)
sr_resid <- read.table("residualised_selfrep_pheno3_nocolnames.pheno")
colnames(sr_resid) <- c("FID", "IID", "antidep_resid")
pd_resid <- read.table("residualised_antidep_pheno1_clean_appt_nocolnames.pheno")
colnames(pd_resid) <- c("FID", "IID", "antidep_resid")
# Subset the phenotypes to those still present after residualising (i.e the individuals in the MWAS)
# Some individuals were removed at this step in the GRM as there were mismatching FIDs
sr %>% filter(IID %in% sr_resid$IID) %>% pull(antidep) %>% table()
pd %>% filter(IID %in% pd_resid$IID) %>% pull(antidep_pheno1) %>% table()
sr <- sr %>% filter(IID %in% sr_resid$IID)
pd <- pd %>% filter(IID %in% pd_resid$IID) 

########## T-tests #########
### Continuous variables 
### Age, Sex, Smoking Pack Years 
## Self-report ##
sr_demographics <- left_join(sr, demographics, by = "IID")
sr_demographics %>% group_by(antidep) %>%
summarise(
mean_age = mean(age, na.rm = T),
age_sd = sd(age, na.rm = T),
mean_bmi = mean(bmi, na.rm = T),
bmi_sd = sd(bmi, na.rm = T),
mean_packyears = mean(pack_years, na.rm = T),
packyears_sd = sd(pack_years, na.rm = T))

sr_age_ttest <- t.test(age ~ antidep, data = sr_demographics, subset = antidep %in% c(0, 1))
sr_bmi_ttest <- t.test(bmi ~ antidep, data = sr_demographics, subset = antidep %in% c(0, 1))
sr_packyears_ttest <- t.test(pack_years ~ antidep, data = sr_demographics, subset = antidep %in% c(0, 1))

## Prescription-derived ##
pd_demographics <- left_join(pd, demographics, by = "IID")
pd_demographics %>% group_by(antidep_pheno1) %>%
summarise(
mean_age = mean(age, na.rm = T),
age_sd = sd(age, na.rm = T),
mean_bmi = mean(bmi, na.rm = T),
bmi_sd = sd(bmi, na.rm = T),
mean_packyears = mean(pack_years, na.rm = T),
packyears_sd = sd(pack_years, na.rm = T))

pd_age_ttest <- t.test(age ~ antidep_pheno1, data = pd_demographics, subset = antidep_pheno1 %in% c(0, 1))
pd_bmi_ttest <- t.test(bmi ~ antidep_pheno1, data = pd_demographics, subset = antidep_pheno1 %in% c(0, 1))
pd_packyears_ttest <- t.test(pack_years ~ antidep_pheno1, data = pd_demographics, subset = antidep_pheno1 %in% c(0, 1))

########## Chi squared tests  #########
## Discrete variables 
## MDD status, sex, smoking behaviours (ordinal)

### Self-Report ####
sr_demographics <- sr_demographics %>%
mutate(antidep_factor = as.factor(antidep),
mdd_factor = as.factor(mdd),
sex_factor = as.factor(sex_coded))
sr_mdd_contig <- table(sr_demographics$antidep_factor, sr_demographics$mdd_factor)
sr_sex_contig <- table(sr_demographics$antidep_factor, sr_demographics$sex_factor)
sr_mdd_chisq <- chisq.test(sr_mdd_contig)
sr_sex_chisq <- chisq.test(sr_sex_contig)