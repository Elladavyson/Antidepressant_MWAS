## making OSCA covariate files ## 

library(tidyverse)
library(readr)
library(dplyr)

# awaiting the updated methylation data #

# baseline model 
# MOA
# MOA model : phenotype ~ GRM-unadjusted M vals + age + sex + batch + AHRR + lymphcytes + monocytes + ORM
# quantitative covariates : age, AHRR probe, monocyte cell counts, lymphocyte cell counts
# qualitative covariates: sex, batch

# read in the target file (all covariates - AHRR), AHRR probe, and sample info (for linking SSI to IID)
target_file_meth <- readRDS('edavyson/GS20K_20_06/GS20k_Targets_18869.rds')
AHRR <- read.table('edavyson/Antidep_methylation/OSCA_data/GRM_uncorrected_30_06/mvals_GRM_uncorrected_29_06_OSCA_AHRR', header = T)
sample_info <- read_table('/Volumes/igmm/GenScotDepression/data/genscot/methylation/20k_GRM_corrected/sample_info.txt') %>% as.data.frame()

target_file_meth <- target_file_meth %>% 
  rowwise() %>% 
  mutate(
  lymphocytes = sum(Bcell, NK, Bcell, CD8T, CD4T)
) %>% as.data.frame()

# recode sex variable, F - 0 and M - 1

target_file_meth <- target_file_meth %>% 
  rowwise() %>%
  mutate(
    sex_coded = ifelse(sex == 'F', 0, 1)) %>% 
  as.data.frame()

target_covs <- target_file_meth %>% select(Sample_Sentrix_ID, sex, age, Batch, Mono, lymphocytes, sex_coded, Sample_Name)
target_covs_info <- merge(target_covs, sample_info %>% select(IID, FID), by.x = 'Sample_Name', by.y = 'IID')
target_covs_info <- merge(target_covs_info, AHRR, by.x = c('FID', 'Sample_Name'), by.y = c('FID', 'IID'))
target_covs_info <- target_covs_info %>% rename(IID = Sample_Name)
quantitative_covs <- target_covs_info %>% select(FID, IID, age, Mono, lymphocytes, cg05575921)
qual_covs <- target_covs_info %>% select(FID, IID, sex_coded, Batch)

# update to script (11/07)
# make sure the FIDs match up with the GRM file
# as the phenotype files have been aligned also 

GRM_ids <- read.table('edavyson/Antidep_methylation/antidep_phenotypes/QCdGS20K.grm.id', header = F)
colnames(GRM_ids) <- c('FID_GRM', 'IID')
# merge with the GRM IIDs and FIDs to assess whether there is mismatch
quantitative_covs_GRM <- merge(quantitative_covs, GRM_ids, by = 'IID')
quantitative_covs_GRM <- quantitative_covs_GRM %>% select(FID, everything())
qual_covs_GRM <- merge(qual_covs, GRM_ids, by = 'IID')
qual_covs_GRM <- qual_covs_GRM %>% select(FID, everything())
# table(quantitative_covs_GRM$FID_GRM == quantitative_covs_GRM$FID)
# replace the FIDs with the FIDs from the GRM file if they do not match 

quantitative_covs_GRM[quantitative_covs_GRM$FID != quantitative_covs_GRM$FID_GRM, 'FID'] <- quantitative_covs_GRM[quantitative_covs_GRM$FID != quantitative_covs_GRM$FID_GRM, 'FID_GRM']  
qual_covs_GRM[qual_covs_GRM$FID_GRM != qual_covs_GRM$FID, 'FID'] <- qual_covs_GRM[qual_covs_GRM$FID_GRM != qual_covs_GRM$FID, 'FID_GRM']

# save the covariate files (currently 106 people missing due to a mismatch in sentrix IDs between the target file and the IID/FID)

write.table(quantitative_covs, 'edavyson/Antidep_methylation/OSCA_data/cov_files_30_06/qcovs_30_06.txt', row.names = F, quote = F)
write.table(quantitative_covs, 'edavyson/Antidep_methylation/OSCA_data/cov_files_30_06/qcovs_30_06_nocolnames',col.names = F, row.names = F, quote = F)
write.table(qual_covs, 'edavyson/Antidep_methylation/OSCA_data/cov_files_30_06/covs_30_06.txt', row.names = F, quote = F)
write.table(qual_covs, 'edavyson/Antidep_methylation/OSCA_data/cov_files_30_06/covs_30_06_nocolnames',col.names = F, row.names = F, quote = F)

# with the FIDs aligned with those in the GRM file 

write.table(quantitative_covs_GRM %>% select(-FID_GRM), 'edavyson/Antidep_methylation/OSCA_data/cov_files_30_06/qcovs_11_07_GRMFID.txt', row.names=F, quote = F)
write.table(quantitative_covs_GRM %>% select(-FID_GRM), 'edavyson/Antidep_methylation/OSCA_data/cov_files_30_06/qcovs_11_07_GRMFID_nocolnames', row.names=F, quote = F, col.names = F)
write.table(qual_covs_GRM %>% select(-FID_GRM), 'edavyson/Antidep_methylation/OSCA_data/cov_files_30_06/covs_11_07_GRMFID.txt', row.names=F, quote = F)
write.table(qual_covs_GRM %>% select(-FID_GRM), 'edavyson/Antidep_methylation/OSCA_data/cov_files_30_06/covs_11_07_GRMFID_nocolnames', row.names=F, quote = F, col.names = F)

write.table(quantitative_covs_GRM %>% select(-FID_GRM), '/antidep_methylation/qcovs_11_07_GRMFID_nocolnames', row.names=F, quote = F, col.names = F)
write.table(qual_covs_GRM %>% select(-FID_GRM), '/antidep_methylation/covs_11_07_GRMFID_nocolnames', row.names=F, quote = F, col.names = F)

# move these files to the scratch space when running the MOA files 


#### if running an adjusted model (I am not, but may need to after getting reviews)
## additional covariates for the adjusted model 
# quantitative covariates : age, AHRR probe, monocyte cell counts, lymphocyte cell counts,
#+ BMI
#+ SIMD rank
# qualitative covariates: sex, batch,
#+ MDD status

#SIMD <- read_csv('data/genscot/phenotypes/SIMD.csv') %>% as.data.frame()
#bmi <- read_csv('Data/GenScot/2023_release/bmi_general_phenos/body.gwasp.202302.csv') %>% as.data.frame()
#scid <- read_tsv('edavyson/Antidep_methylation/data/genscot_depression.tsv') %>% as.data.frame()
# the scid column is SCID_Diagnosis column from the SCID QC file
# 0 - No Major Disorder, 1 - Single MDD, 2 - Recurrent MDD and 3 - Bipolar Disorder 
# recoded to 0-0, 1- 1/2 and NA - 3

