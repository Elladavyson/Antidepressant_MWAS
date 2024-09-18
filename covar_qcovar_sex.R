library(tidyr)
library(dplyr)
library(tibble)
# Read in the phenotype files (residualised)
pd_all <- read.table('residualised_antidep_pheno1_clean_appt_nocolnames.pheno')
colnames(pd_all) <- c("FID", "IID", "antidep")
sr_all <- read.table("residualised_selfrep_pheno3_nocolnames.pheno")
colnames(sr_all) <- c("FID", "IID", "antidep")
# Read in the actual phenotypes (0/1) for the counts (check against log file)
pd_raw <- read.csv("antidep_pheno1_clean_appt.csv", header = T)
sr_raw <- read.csv("selfrep_pheno3_methyl_03_05.csv", header = T)
# Read in the covariate and quantitative covariate file 
covs <- read.table("covs_11_07_GRMFID_nocolnames")
colnames(covs) <- c("FID", "IID", "sex_coded", "Batch")
qcovs <- read.table("qcovs_11_07_GRMFID_nocolnames")
colnames(qcovs) <- c("FID", "IID", "age", "Mono", "lymphocytes", "cg05575921")

# Left join with the phenotype file
pd_covs <- left_join(pd_all %>% rename(FID_pheno=FID), covs, by = "IID")
sr_covs <- left_join(sr_all %>% rename(FID_pheno=FID), covs, by = "IID")
pd_covs_raw <- left_join(pd_raw %>% rename(FID_pheno=FID), covs, by = "IID")
pd_covs_raw <- pd_covs_raw %>% mutate(sex_name = ifelse(sex_coded ==0, "Female", "Male"))
sr_covs_raw <- left_join(sr_raw %>% rename(FID_pheno=FID), covs, by = "IID")
sr_covs_raw <- sr_covs_raw %>% mutate(sex_name = ifelse(sex_coded==0, "Female", "Male"))
# Filter out the individuals with mismatched FIDs to the methylation file
sr_covs_raw <- sr_covs_raw %>% filter(FID_pheno == FID) %>% filter(!is.na(sex_coded))
write.table(sr_covs_raw %>% 
select(FID, IID, antidep), "/exports/igmm/eddie/GenScotDepression/users/edavyson/antidep_project/antidep_phenotypes/selfreport_pheno3_MOA_match_16531.tsv", sep = "\t", row.names = F, quote = F)
# Merge with the phenotype file 

pd_sex <- addmargins(table(pd_covs_raw$antidep_pheno1, pd_covs_raw$sex_name, dnn =c( "Antidepressant Exposure", "Sex")),FUN = list(Total = sum))%>% as.data.frame.matrix()
colnames(pd_sex) <- make.names(colnames(pd_sex), unique = TRUE)
pd_sex <- pd_sex %>%
rownames_to_column(var = "antidep_exposure") %>% 
mutate(phenotype = "Prescription_derived")
pd_sex <- pd_sex %>% select(phenotype, everything())

sr_sex <- addmargins(table(sr_covs_raw$antidep, sr_covs_raw$sex_name, dnn =c( "Antidepressant Exposure", "Sex")),FUN = list(Total = sum)) %>% as.data.frame.matrix()
sr_sex <- sr_sex %>%
rownames_to_column(var = "antidep_exposure") %>% 
mutate(phenotype = "Self_report")
sr_sex <- sr_sex %>% select(phenotype, everything())


# Both  
all_sex_summary <- rbind(sr_sex, pd_sex)
# Reformatting to match paper table 
all_sex_long <- all_sex_summary %>%
  pivot_longer(cols = c(Male, Female, Total),
               names_to = "Category",
               values_to = "Value")
sex_long_percentages <- all_sex_long %>%
group_by(phenotype, antidep_exposure) %>% 
mutate(Total_for_group = Value[Category == "Total"],
Percentage = (Value/Total_for_group)*100) %>%
 ungroup() %>%
 mutate(value_perc= paste0(Value, " (", round(Percentage,1), ")")) %>%
 select(-c(Total_for_group, Value, Percentage))
all_sex_wide <- sex_long_percentages %>%
  unite("group", phenotype, antidep_exposure, sep = "_") %>%  # Create a unique identifier for columns
  pivot_wider(names_from = group, values_from = value_perc) 

outdir <- "/exports/igmm/eddie/GenScotDepression/users/edavyson/antidep_project/revisions/output/MWAS_sex_segregated/"
write.table(all_sex_wide,  paste0(outdir, "sex_phenotype_numbers.tsv"), sep = "\t", row.names = F, quote = F)
# Separate the covariate file into the males and females - OSCA_covariates.R file Females (0) and Male (1)
table(covs$sex_coded, useNA = "always")
female_covs <- covs %>% filter(sex_coded == 0) %>% select(-sex_coded)
male_covs <- covs %>% filter(sex_coded == 1) %>% select(-sex_coded)
pd_female <- pd_covs %>% filter(sex_coded == 0) %>% select(FID, IID, antidep)
pd_male <- pd_covs %>% filter(sex_coded == 1) %>% select(FID, IID, antidep)
sr_female <- sr_covs %>% filter(sex_coded == 0) %>% select(FID, IID, antidep)
sr_male <- sr_covs %>% filter(sex_coded == 1) %>% select(FID, IID, antidep)

write.table(female_covs, paste0(outdir, "covs_11_07_GRMFID_female_nocolnames"), row.names =F, col.names = F, quote = F)
write.table(male_covs, paste0(outdir, "covs_11_07_GRMFID_male_nocolnames"), row.names =F, col.names = F, quote = F)
write.table(pd_female, paste0(outdir, "residualised_antidep_pheno1_clean_appt_female_nocolnames.pheno"), row.names =F, col.names = F, quote = F)
write.table(pd_male, paste0(outdir, "residualised_antidep_pheno1_clean_appt_male_nocolnames.pheno"), row.names =F, col.names = F, quote = F)
write.table(sr_female, paste0(outdir, "residualised_selfrep_pheno3_female_nocolnames.pheno"), row.names =F, col.names = F, quote = F)
write.table(sr_male, paste0(outdir, "residualised_selfrep_pheno3_male_nocolnames.pheno"), row.names =F, col.names = F, quote = F)
