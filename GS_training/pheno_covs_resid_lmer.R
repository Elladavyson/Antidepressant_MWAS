######## Regressing the GRM-residualised phenotype against the covariates in the MWAS model ########

#--------------------------------------------------------------------------------

# Reading in files 

#--------------------------------------------------------------------------------
  
phenotype <- read.table('/Volumes/igmm/GenScotDepression/users/edavyson/Antidep_methylation/antidep_phenotypes/resid_phenos/residualised_selfrep_pheno3_nocolnames.pheno', header = F)
colnames(phenotype) <- c('FID', 'IID', 'resid_antidep')

# /exports/igmm/datastore/GenScotDepression/users/edavyson/Antidep_methylation/OSCA_data/cov_files_30_06

print('Reading in covariate files')
qcovs <- read_table('/Volumes/igmm/GenScotDepression/users/edavyson/Antidep_methylation/OSCA_data/cov_files_30_06/qcovs_11_07_GRMFID.txt') %>% as.data.frame()
covs <- read_table('/Volumes/igmm/GenScotDepression/users/edavyson/Antidep_methylation/OSCA_data/cov_files_30_06/covs_11_07_GRMFID.txt') %>% as.data.frame()


#--------------------------------------------------------------------------------

# Merge files 

#--------------------------------------------------------------------------------

all_covs <- merge(qcovs, covs, by = c('FID', 'IID'))
all_covs_pheno <- merge(all_covs, phenotype, by = c('FID', 'IID'))

#--------------------------------------------------------------------------------

# Run regression model of phenotype ~ covariates 

#--------------------------------------------------------------------------------

pheno_cov_model <- lmer(resid_antidep ~ scale(age) + scale(Mono) + scale(lymphocytes) + scale(cg05575921) + as.factor(sex_coded) + (1|Batch), data=all_covs_pheno)

#--------------------------------------------------------------------------------

# Extract and save residuals for MRS analysis 

#--------------------------------------------------------------------------------

pheno_residuals <- residuals(pheno_cov_model)
pheno_residuals <- cbind(all_covs_pheno, pheno_residuals) %>% select(c(FID, IID, pheno_residuals))
write.table(pheno_residuals, '/Volumes/igmm/GenScotDepression/users/edavyson/Antidep_methylation/MRS/selfrep_cov_residuals.txt', row.names = F, quote = F)


