
## Adding missingness for MS script testing in GS 
## Specifically the MRS_calc.R script which now has lots of graphs/measurements for the levels
## of missingness within CpGs

.libPaths('/exports/igmm/eddie/GenScotDepression/users/edavyson/R/x86_64-pc-linux-gnu-library/4.1')
library(tidyverse)
library(data.table)
library(dplyr)
setwd('/exports/eddie/scratch/s2112198')

sink('/exports/eddie/scratch/s2112198/GS_simulate_missingness_09_01_2024.log')
DNAm <- read_table('GS_DNAm_preproc.txt') %>% as.data.frame()

set.seed(123)  # Set seed for reproducibility

# Generate random proportions from uniform dist (runif) of missing values for 1000 columns (CpGs)
print('Randomly sampling the proportion of missingness to each CpG (100)')
missing_proportions <- runif(100, 0, 1)  # Adjust the range as needed

# Introduce random missingness to each column based on the generated proportions
print('Replacing the random proportions of missing values for 100 CpGs')
for (i in seq_along(missing_proportions)) {
  col <- names(DNAm %>% select(-IID))[i]
  print(paste0(col, ':', missing_proportions[i]))
  DNAm[sample(nrow(DNAm), round(missing_proportions[i] * nrow(DNAm))), col] <- NA
}

# Table recording the missingness in each CpG
print('The proportion of missingness in each CpG')
missing_sim <- data.frame(CpG = colnames(DNAm)[2:101], missing = missing_proportions)
print(missing_sim)
# Print a summary to verify the missingness

summary(DNAm)

print('Saving the missingness proportions to missing_GS_proportions_09_01_2024.txt')
print('Saving the new DNAm file with missingness in to missing_GS_DNAm_preproc.txt for testing')
# save the new DNAm object with simulated missingness 
write.table(missing_sim, 'missing_GS_proportions_09_01_2024.txt')
write.table(DNAm, 'missing_GS_DNAm_preproc.txt', row.names = F, quote = F)

sink()

# Testing the DNAm object with missingvalues under the cohort name missingGS in the MRS_calc.R script

# Rscript MRS_calc.R --cohort missingGS --DNAm missing_GS_DNAm_preproc.txt --weights GS_AD_MRS_weights.txt --id_column IID --pheno selfrep_pheno3_methyl_03_05.csv --outdir /exports/eddie/scratch/s2112198/
