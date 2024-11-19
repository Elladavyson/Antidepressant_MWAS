################################################################################

## LBC methylation pre-processing 

################################################################################

.libPaths('/exports/igmm/eddie/GenScotDepression/users/edavyson/R/x86_64-pc-linux-gnu-library/4.1')
library(lumi) #beta2m function
library(readr)
library(data.table)
library(dplyr)

################################################################################

## Read in data 

################################################################################

targets <- readRDS('targets_3489_bloodonly.rds')
betas <- readRDS('LBC_betas_3489_bloodonly.rds')

################################################################################

## Transform the methylation data 

################################################################################
# transpose the data 
# rows will be IDs, and columns CpGs

betas_t <- t(betas)

# convert the betas to M values

mvals <- beta2m(betas_t) 
rm(betas)
rm(betas_t)

mvals <- as.data.frame(mvals)


################################################################################

## Subsetting to 1936 and 1921

################################################################################

mvals$basename <- row.names(mvals)
mvals <- merge(mvals, targets, by.x = 'basename', by.y = 'Basename')

mvals_36 <- mvals %>%
  filter(cohort == 'LBC36' & WAVE == 1) %>%
  select(ID_raw, starts_with("cg"))

mvals_21 <- mvals %>%
  filter(cohort == 'LBC21' & WAVE == 1) %>%
  select(ID_raw, starts_with("cg"))


################################################################################

## Saving files

################################################################################

write_rds(mvals_36, '/exports/eddie/scratch/s2112198/LBC36_mvals_wave1.rds')
write_rds(mvals_21, '/exports/eddie/scratch/s2112198/LBC21_mvals_wave1.rds')



         