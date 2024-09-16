########### Assess the prevalence of any other psychiatric disorders in the MDD-only group #############
.libPaths(c(.libPaths(), '/exports/eddie3_homes_local/s2112198/R/x86_64-pc-linux-gnu-library/4.4'))
library(dplyr)
library(readr)
library(data.table)
library(UpSetR)
library(tibble)
library(ggplot2)

########### Read in antidepressant phenotypes
### Copied from /exports/igmm/datastore/GenScotDepression/users/edavyson/Antidep_methylation/antidep_phenotypes/

## Self-report phenotype 
selfrep <- read.csv("selfrep_pheno3_methyl_03_05.csv", header = T)
selfrep_resid <- read.table("residualised_selfrep_pheno3_nocolnames.pheno", header = F)
colnames(selfrep_resid) <- c("FID", "IID", "resid_antidep")
selfrep_mrg <- merge(selfrep, selfrep_resid %>% select(-FID), by = c( "IID"))
table(selfrep_mrg$antidep) # 15028 antidepressant unexposed and 1508 antidepressant exposed
## Prescription-derived phenotype 
pd <- read.csv("antidep_pheno1_clean_appt.csv", header =T)
pd_resid <- read.table("residualised_antidep_pheno1_clean_appt_nocolnames.pheno", header =F)
colnames(pd_resid) <- c("FID", "IID", "resid_antidep")
pd_mrg <- merge(pd, pd_resid %>% select(-FID), by = "IID")
table(pd_mrg$antidep) # 7090 unexposed and 861 exposed 

########### Assess the prevalence of other psychiatric disorders in this full sample 
### the MDD phenotype 
fam <- read_table('data/genetics/genotypes/GS20K_PLINK_files/QCd_data/QCdGS20K.fam', col_names=c('FID', 'IID', 'father', 'mother', 'sex', 'pheno'))
# age/sex at baseline
agesex <- read_csv('data/phenotypes/agesex_all.csv')
# baseline SCID
scid_qc <- read_csv('data/phenotypes/SCID_QC_201113_GWASids.csv')
# STRADL
# loads a data.frame called 'x'
stradl <- {load('data/phenotypes/GS_Recontact/STRADL.Rdata'); x}
# Participants to exclude based on health records (SMR)
smr_exclude <- read_tsv('data/phenotypes/SMRexclude.txt', col_names=c('FID', 'IID'))
# recode and merge
gs_depression <- 
fam %>%
select(FID, IID) %>%
# select SCID diagnosis
inner_join(scid_qc %>% select(gwas, SCID_Diagnosis), by=c('IID'='gwas')) %>%
# recode scid diagnosis as 0/1
# controls  0 -> 0, single 1 -> 1, recurrent 2 -> 1, bipolar 3 -> NA
mutate(scid=case_when(SCID_Diagnosis == 0 ~ 0L,
                      SCID_Diagnosis %in% 1:2 ~ 1L,
                      TRUE ~ NA_integer_)) %>%
# STRADL follow up
left_join(stradl %>% select(id, CIDI_MDD), by=c('IID'='id')) %>%
# recode cidi to just 0/1
# Control = 0     MDD = 1     Bipolar = 2     Hypomanic = 3
mutate(cidi=case_when(CIDI_MDD == 0 ~ 0L,
                      CIDI_MDD == 1 ~ 1L,
                      TRUE ~ NA_integer_))


# remove intermediate columns
#select(-SCID_Diagnosis, -CIDI_MDD)

# Merge the antidepressant_exposure phenotypes with the depression phenotypes 
selfrep_mrg <- left_join(selfrep_mrg, gs_depression, by = "IID")
selfrep_mrg <- selfrep_mrg %>% mutate(
    SCID_Diagnosis_text = case_when(
        SCID_Diagnosis == 0 ~ "No Major Disorder",
        SCID_Diagnosis == 1 ~ "Single MDD Episode",
        SCID_Diagnosis == 2 ~ "Recurrent MDD",
        SCID_Diagnosis == 3 ~ "Bipolar"
    ),
    antidep_text = ifelse(antidep == 0, "Unexposed", "Exposed")
)
table(selfrep_mrg$SCID_Diagnosis, useNA = "always") # 
selfrep_scid <- addmargins(table(selfrep_mrg$antidep_text,selfrep_mrg$SCID_Diagnosis_text, useNA = "always", dnn = c("Antidepressant exposure", "SCID Diagnosis")), FUN = list(Total = sum))%>%
as.data.frame.matrix() # Dsitribution of SCID diagnoses in the full sample  
selfrep_scid <- selfrep_scid[-3,]
colnames(selfrep_scid) <- make.names(colnames(selfrep_scid), unique = TRUE)
selfrep_scid <- selfrep_scid %>% 
rownames_to_column(var = "antidep_exposure") %>% 
mutate(phenotype = "Self_report")
selfrep_scid <- selfrep_scid %>% select(phenotype, everything())


selfrep_scid$SCID_Diagnosis <- rownames(selfrep_scid)
selfrep_scid <- selfrep_scid[, c(ncol(selfrep_scid), 1:(ncol(selfrep_scid)-1))] 
outdir <- "/exports/igmm/eddie/GenScotDepression/users/edavyson/antidep_project/revisions/output/SCID_SPQ/"
write.table(selfrep_scid, paste0(outdir, "selfrep_scid_summary.tsv"), sep = "\t", row.names = F, quote =F )      

pd_mrg <- pd_mrg %>% rename(antidep=antidep_pheno1)
pd_mrg <- left_join(pd_mrg, gs_depression, by = "IID")
pd_mrg <- pd_mrg %>% mutate(
    SCID_Diagnosis_text = case_when(
        SCID_Diagnosis == 0 ~ "No Major Disorder",
        SCID_Diagnosis == 1 ~ "Single MDD Episode",
        SCID_Diagnosis == 2 ~ "Recurrent MDD",
        SCID_Diagnosis == 3 ~ "Bipolar"
    ),
    antidep_text = ifelse(antidep == 0, "Unexposed", "Exposed")
)
table(pd_mrg$SCID_Diagnosis, useNA = "always") # 
pd_scid <- addmargins(table(pd_mrg$antidep_text, pd_mrg$SCID_Diagnosis_text, useNA = "always", dnn = c("Antidepressant exposure","SCID Diagnosis")), FUN = list(Total = sum))%>%
as.data.frame.matrix() # Dsitribution of SCID diagnoses in the full sample  
pd_scid <- pd_scid[-3,]
colnames(pd_scid) <- make.names(colnames(pd_scid), unique = TRUE)
pd_scid <- pd_scid %>% 
rownames_to_column(var = "antidep_exposure") %>% 
mutate(phenotype = "Prescription_derived")
pd_scid <- pd_scid %>% select(phenotype, everything())

write.table(pd_scid, paste0(outdir, "pd_scid_summary.tsv"), sep = "\t", row.names = F, quote =F )      
all_scid_summary <- rbind(selfrep_scid, pd_scid)
write.table(all_scid_summary, paste0(outdir, "pd_sr_scid_summary.tsv"), sep = "\t", row.names = F, quote =F )      

# Look at how many SMR exclude participants there are 
smr_diag <- read.csv("SMRDiag.csv", header = T)
smr_exclude <- left_join(smr_exclude, smr_diag, by= c("IID"="gwas"))
smr_exclude$smr_criteria <- apply(smr_exclude, 1, function(x) {
    paste(names(smr_exclude)[which(x == TRUE)], collapse = ":")
})
smr_exclude$smr_cond <- apply(smr_exclude %>% select(-starts_with('SMR')), 1, function(x) {
    paste(names(smr_exclude)[which(x == TRUE)], collapse = ":")
})

smr_lists <- list(depression = smr_exclude %>% filter(depression==TRUE) %>% .$IID, bipolar = smr_exclude %>% filter(bipolar==TRUE) %>% .$IID,
schizophrenia_broad = smr_exclude %>% filter(schiz_broad==TRUE) %>% .$IID, szhizophrenia_narrow = smr_exclude
 %>% filter(schiz_narrow==TRUE) %>% .$IID,
all_outpatient = smr_exclude %>% filter(SMR00_outpatient==TRUE) %>% .$IID, inpatient = smr_exclude %>% filter(SMR04_inpatient==TRUE) %>% .$IID)
names(smr_lists) <- c("Depression ICD10", "Bipolar ICD10", "Broad Schizophrenia ICD10","Narrow Schizophrenia (F20 only)",
"SMR00 Outpatient", "SMR04 Inpatient")
upset_plot <- upset(fromList(smr_lists), nsets = 6, order.by = 'freq', set_size.show = TRUE, mb.ratio = c(0.6, 0.4))
png(paste0(outdir, "smr_exclude_upset.png"),width = 800, height = 600)
upset_plot
dev.off()

# Merge this SMR information with the self-report 
selfrep_mrg_smr <- left_join(selfrep_mrg, smr_exclude, by = c("IID"))
selfrep_mrg_smr <- selfrep_mrg_smr %>%
mutate(bip_scz= ifelse(bipolar==TRUE|schiz_broad == TRUE | schiz_narrow == TRUE, 1,0))

table(selfrep_mrg_smr$bip_scz, selfrep_mrg_smr$antidep_text, useNA = "always")

pd_mrg_smr <- left_join(pd_mrg, smr_exclude, by = c("IID"))
pd_mrg_smr <- pd_mrg_smr %>%
mutate(bip_scz= ifelse(bipolar==TRUE|schiz_broad == TRUE | schiz_narrow == TRUE, 1,0))

table(pd_mrg_smr$bip_scz, pd_mrg_smr$antidep_text, useNA = "always")

# Make a data-frame detailing the number of individuals with Bipolar/Schizophrenia related records using the scottish morbidity records
## BIP only - individuals with a bipolar code and no schizophrenia codes 
## SCZ - individuals with schizophrenia code and no bipolar code 
## BIP_SCZ - individuals with a schizophrenia code and a bipolar code 
SMR_information <- data.frame(Phenotype=c(rep("Self-report",2), rep("Prescription_derived",2)),
Exposure_status = rep(c("Exposed", "Unexposed"),2),
N_BIP_only = c(selfrep_mrg_smr %>% filter(antidep == 1 & bipolar == TRUE & schiz_broad == FALSE & schiz_narrow == FALSE) %>% pull(IID) %>% unique() %>% length(),
selfrep_mrg_smr %>% filter(antidep==0 & bipolar == TRUE & schiz_broad == FALSE & schiz_narrow == FALSE) %>% pull(IID) %>% unique() %>% length(),
pd_mrg_smr %>% filter(antidep == 1 & bipolar == TRUE & schiz_broad == FALSE & schiz_narrow == FALSE) %>% pull(IID) %>% unique() %>% length(),
pd_mrg_smr %>% filter(antidep==0 & bipolar == TRUE & schiz_broad == FALSE & schiz_narrow == FALSE) %>% pull(IID) %>% unique() %>% length()),
N_SCZ_only = c(selfrep_mrg_smr %>% filter(antidep == 1 & (schiz_broad == TRUE | schiz_narrow == TRUE) & bipolar == FALSE) %>% pull(IID) %>% unique() %>% length(),
selfrep_mrg_smr %>% filter(antidep==0 & (schiz_broad == TRUE | schiz_narrow == TRUE) & bipolar == FALSE) %>% pull(IID) %>% unique() %>% length(),
pd_mrg_smr %>% filter(antidep == 1 & (schiz_broad == TRUE | schiz_narrow == TRUE) & bipolar == FALSE) %>% pull(IID) %>% unique() %>% length(),
pd_mrg_smr %>% filter(antidep==0 & (schiz_broad == TRUE | schiz_narrow == TRUE) & bipolar == FALSE) %>% pull(IID) %>% unique()%>% length()),
N_BIP_SCZ = c(selfrep_mrg_smr %>% filter(antidep == 1 & bipolar == TRUE & (schiz_broad == TRUE | schiz_narrow == TRUE)) %>% pull(IID) %>% unique() %>% length(),
selfrep_mrg_smr %>% filter(antidep==0 & bipolar == TRUE & (schiz_broad == TRUE | schiz_narrow == TRUE)) %>% pull(IID) %>% unique() %>% length(),
pd_mrg_smr %>% filter(antidep == 1 & bipolar == TRUE & (schiz_broad == TRUE | schiz_narrow == TRUE)) %>% pull(IID) %>% unique() %>% length(),
pd_mrg_smr %>% filter(antidep==0 & bipolar == TRUE & (schiz_broad == TRUE | schiz_narrow == TRUE)) %>% pull(IID) %>% unique() %>% length()

)
)
# Save this dataframe
write.table(SMR_information, paste0(outdir, "SMR_info_antidep_phenotypes_all.tsv"), sep = "\t", row.names = F, quote = F)

### SPQ information 
# Schizotypal personality questionnaire 
spq <- read.csv("SPQ.csv", header = T)

# Getting the numbers for both phenotypes
pd_mrg_spq <- left_join(pd_mrg, spq, by = c("IID"="ID"))
pd_mrg_spq <- pd_mrg_spq %>% mutate(spq_thres = ifelse(total >= 17, "SPQ", "No SPQ"))
sr_mrg_spq <- left_join(selfrep_mrg, spq, by = c("IID"="ID"))
sr_mrg_spq <- sr_mrg_spq %>% mutate(spq_thres = ifelse(total >= 17, "SPQ", "No SPQ"))
# Table of SPQ diagnoses by antidepressant exposure 
# PD
pd_spq <- addmargins(table(pd_mrg_spq$antidep_text, pd_mrg_spq$spq_thres, useNA = "always", dnn = c("Antidepressant exposure","SPQ")), FUN = list(Total = sum))%>%
as.data.frame.matrix() # Number of SPQ
pd_spq  <- pd_spq[-3,]
colnames(pd_spq) <- make.names(colnames(pd_spq), unique = TRUE)
pd_spq <- pd_spq %>% 
rownames_to_column(var = "antidep_exposure") %>% 
mutate(phenotype = "Prescription_derived")
pd_spq <- pd_spq %>% select(phenotype, everything())
# SR
sr_spq <- addmargins(table(sr_mrg_spq$antidep_text, sr_mrg_spq$spq_thres, useNA = "always", dnn = c("Antidepressant exposure","SPQ")), FUN = list(Total = sum))%>%
as.data.frame.matrix() # Number of SPQ
sr_spq  <- sr_spq[-3,]
colnames(sr_spq) <- make.names(colnames(sr_spq), unique = TRUE)
sr_spq <- sr_spq %>% 
rownames_to_column(var = "antidep_exposure") %>% 
mutate(phenotype = "Self_report")
sr_spq <- sr_spq %>% select(phenotype, everything())
# Combine the two 
spq_both <- rbind(sr_spq, pd_spq)
# Save the output 
write.table(spq_both, paste0(outdir, "SPQ_both_phenos.tsv"), sep = "\t", row.names = F, quote = F)

# Graphs of the distribution of the SPQ scores in the phenotypes 
spq_dist_df <- rbind(pd_mrg_spq %>% 
mutate(phenotype = "Prescription_derived") %>% 
select(IID, antidep, total, phenotype),
sr_mrg_spq %>% 
mutate(phenotype = "Self_report") %>% 
select(IID, antidep, total, phenotype)
)
# Look at SPQ scores across the different antidepressants
spq_dist_plt <- ggplot(spq_dist_df, aes(y = total, fill = as.factor(antidep))) + 
geom_boxplot() +
theme_bw() + 
labs(x = "SPQ score", y = "Count")+
facet_wrap(~phenotype)
ggsave(paste0(outdir, "spq__antidep_dists.png"), plot = spq_dist_plt, width = 8, height = 6, device = "png", dpi = 300)
