## Seeing if there are statistical differences in the demographic variables ##
## Using R version 4.4
library(dplyr)
library(tidyverse)
library(readr)
library(ggplot2)
library(effsize)
# Personal library for 4.4 : /exports/eddie3_homes_local/s2112198/R/x86_64-pc-linux-gnu-library/4.4
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

# Distributions 
outdir <- "/exports/igmm/eddie/GenScotDepression/users/edavyson/antidep_project/revisions/output/"

sr_demographics <- left_join(sr, demographics, by = "IID")
hist_save <- function(label, colname) {
if (label =="self_report") {
    dataset <- sr_demographics
} else if (label == "prescription_derived") {
    dataset <- pd_demographics
} else(
    stop("Please give valid label")
)
ggsave(paste0(outdir,label, "_", colname, "_hist.png"),
ggplot(dataset, aes(x = !!sym(colname))) +
geom_histogram() + 
facet_wrap(~antidep),
width = 8,
height = 6,
device = "png",
dpi= 300)
}
hist_save("self_report", "age")
hist_save("self_report", "bmi")
hist_save("self_report", "pack_years")

sr_demographics %>% group_by(antidep) %>%
summarise(
mean_age = mean(age, na.rm = T),
age_sd = sd(age, na.rm = T),
mean_bmi = mean(bmi, na.rm = T),
bmi_sd = sd(bmi, na.rm = T),
mean_packyears = mean(pack_years, na.rm = T),
packyears_sd = sd(pack_years, na.rm = T))

sr_age_ttest <- t.test(age ~ antidep, data = sr_demographics, subset = antidep %in% c(0, 1))
sr_age_cohen <- cohen.d(age~as.factor(antidep), data = sr_demographics)
sr_bmi_ttest <- t.test(bmi ~ antidep, data = sr_demographics, subset = antidep %in% c(0, 1))
sr_bmi_cohen <- cohen.d(bmi~as.factor(antidep), data = sr_demographics)
sr_packyears_ttest <- t.test(pack_years ~ antidep, data = sr_demographics, subset = antidep %in% c(0, 1))
sr_packyears_cohen <- cohen.d(pack_years~as.factor(antidep), data = sr_demographics)

## Prescription-derived ##
pd_demographics <- left_join(pd, demographics, by = "IID")
pd_demographics <- pd_demographics %>% rename(antidep=antidep_pheno1)
hist_save("prescription_derived", "age")
hist_save("prescription_derived", "bmi")
hist_save("prescription_derived", "pack_years")

pd_demographics %>% group_by(antidep) %>%
summarise(
mean_age = mean(age, na.rm = T),
age_sd = sd(age, na.rm = T),
mean_bmi = mean(bmi, na.rm = T),
bmi_sd = sd(bmi, na.rm = T),
mean_packyears = mean(pack_years, na.rm = T),
packyears_sd = sd(pack_years, na.rm = T))

pd_age_ttest <- t.test(age ~ antidep, data = pd_demographics, subset = antidep %in% c(0, 1))
pd_age_cohen <- cohen.d(age~as.factor(antidep), data = pd_demographics)
pd_bmi_ttest <- t.test(bmi ~ antidep, data = pd_demographics, subset = antidep %in% c(0, 1))
pd_bmi_cohen <- cohen.d(bmi~as.factor(antidep), data = pd_demographics)
pd_packyears_ttest <- t.test(pack_years ~ antidep, data = pd_demographics, subset = antidep %in% c(0, 1))
pd_packyears_cohen <- cohen.d(pack_years~as.factor(antidep), data = pd_demographics)

# Save the t-test results
ttest_res <- data.frame(
    antidep_var = c(rep("Self report", 3), rep("Prescription-derived", 3)),
    variable = rep(c("Age", "BMI", "pack years"),2),
    tstat = c(sr_age_ttest$statistic, sr_bmi_ttest$statistic, sr_packyears_ttest$statistic, pd_age_ttest$statistic, pd_bmi_ttest$statistic, pd_packyears_ttest$statistic),
    df = c(sr_age_ttest$parameter, sr_bmi_ttest$parameter, sr_packyears_ttest$parameter,pd_age_ttest$parameter, pd_bmi_ttest$parameter, pd_packyears_ttest$parameter ),
    p_val = c(sr_age_ttest$p.value, sr_bmi_ttest$p.value, sr_packyears_ttest$p.value, pd_age_ttest$p.value, pd_bmi_ttest$p.value, pd_packyears_ttest$p.value),

    conf_int_low = c(sr_age_ttest$conf.int[[1]], sr_bmi_ttest$conf.int[[1]], sr_packyears_ttest$conf[[1]], pd_age_ttest$conf.int[[1]], pd_bmi_ttest$conf.int[[1]], pd_packyears_ttest$conf[[1]]),
    conf_int_high = c(sr_age_ttest$conf.int[[2]], sr_bmi_ttest$conf.int[[2]], sr_packyears_ttest$conf[[2]], pd_age_ttest$conf.int[[2]], pd_bmi_ttest$conf.int[[2]], pd_packyears_ttest$conf[[2]]),
    cohens_D = c(sr_age_cohen$estimate, sr_bmi_cohen$estimate, sr_packyears_cohen$estimate, pd_age_cohen$estimate, pd_bmi_cohen$estimate, pd_packyears_cohen$estimate)) 
round_numeric <- function(df) {
    df[] <- lapply(df, function(col) {
        if (is.numeric(col)) {
            signif(col, digits = 3)
        } else {
            col
        }
    }) 
    return(df)
}
ttest_res <- round_numeric(ttest_res)
write.table(ttest_res, paste0(outdir, "age_bmi_packyears_ttest.tsv"), sep= "\t", row.names = F, quote = F)

########## Chi squared tests  #########
## Discrete variables 
## MDD status, sex, smoking behaviours (ordinal)

### Self-Report ####
sr_demographics <- sr_demographics %>%
mutate(antidep_factor = as.factor(antidep),
mdd_factor = as.factor(mdd),
sex_factor = as.factor(sex_coded),
smoking_var = case_when(
   ever_smoke == 1 ~ "Current",
   ever_smoke == 2 | ever_smoke == 3 ~ "Former",
   ever_smoke == 4 ~ "Never",
   TRUE  ~ NA_character_ 
),
smoking_var_factor = factor(smoking_var, levels = c("Current", "Former", "Never"), ordered = TRUE))
bar_save <- function(label, colname) {
if (label =="self_report") {
    dataset <- sr_demographics
} else if (label == "prescription_derived") {
    dataset <- pd_demographics
} else {
    stop("Please give valid label 'self-report' or 'prescription-derived'")
}
ggsave(paste0(outdir,label, "_", colname, "_bar.png"),
ggplot(dataset, aes(x = as.factor(antidep), fill = !!sym(colname))) +
geom_bar(stat= "count", position = "dodge")+theme_bw()+
labs(x = paste0(label, "antidepressant exposure"), y = "Count"),
width = 8,
height = 6,
device = "png",
dpi= 300)
}
bar_save("self_report", "mdd_factor")
bar_save("self_report", "sex_factor")
bar_save("self_report", "smoking_var_factor")
sr_mdd_contig <- table(sr_demographics$antidep_factor, sr_demographics$mdd_factor)
sr_sex_contig <- table(sr_demographics$antidep_factor, sr_demographics$sex_factor)
sr_smoke_contig <- table(sr_demographics$antidep_factor, sr_demographics$smoking_var_factor)
sr_mdd_chisq <- chisq.test(sr_mdd_contig)
sr_sex_chisq <- chisq.test(sr_sex_contig)
sr_smoke_chisq <- chisq.test(sr_smoke_contig)

### prescription derived ###
pd_demographics <- pd_demographics %>%
mutate(antidep_factor = as.factor(antidep),
mdd_factor = as.factor(mdd),
sex_factor = as.factor(sex_coded),
smoking_var = case_when(
   ever_smoke == 1 ~ "Current",
   ever_smoke == 2 | ever_smoke == 3 ~ "Former",
   ever_smoke == 4 ~ "Never",
   TRUE  ~ NA_character_ 
),
smoking_var_factor = factor(smoking_var, levels = c("Current", "Former", "Never"), ordered = TRUE))
bar_save("prescription_derived", "mdd_factor")
bar_save("prescription_derived", "sex_factor")
bar_save("prescription_derived", "smoking_var_factor")
pd_mdd_contig <- table(pd_demographics$antidep_factor, pd_demographics$mdd_factor)
pd_sex_contig <- table(pd_demographics$antidep_factor, pd_demographics$sex_factor)
pd_smoke_contig <- table(pd_demographics$antidep_factor, pd_demographics$smoking_var_factor)
pd_mdd_chisq <- chisq.test(pd_mdd_contig)
pd_sex_chisq <- chisq.test(pd_sex_contig)
pd_smoke_chisq <- chisq.test(pd_smoke_contig)

### Saving the chi-squared results 
chisquare_res <- data.frame(
    antidep_var = c(rep("Self report", 3), rep("Prescription-derived", 3)),
    variable = rep(c("Sex", "Smoking status", "Lifetime MDD"),2),
    chi_sqaure = c(sr_sex_chisq$statistic, sr_smoke_chisq$statistic, sr_mdd_chisq$statistic, pd_sex_chisq$statistic, pd_smoke_chisq$statistic, pd_mdd_chisq$statistic),
    df = c(sr_sex_chisq$parameter, sr_smoke_chisq$parameter, sr_mdd_chisq$parameter, pd_sex_chisq$parameter, pd_smoke_chisq$parameter, pd_mdd_chisq$parameter),
    pval = c(sr_sex_chisq$p.value, sr_smoke_chisq$p.value, sr_mdd_chisq$p.value, pd_sex_chisq$p.value, pd_smoke_chisq$p.value, pd_mdd_chisq$p.value))
chisquare_res <- round_numeric(chisquare_res)
write.table(chisquare_res, paste0(outdir, "sex_smoke_mdd_chisquare_res.tsv"), sep = "\t", row.names = F, quote = F)