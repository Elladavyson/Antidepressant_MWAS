
###############################################################################

# Demographic information about sample from the cohort 
# age, sex, bmi, smoking behaviors, and lifetime MDD status (where applicable)

###############################################################################
.libPaths('/exports/igmm/eddie/GenScotDepression/users/edavyson/R/x86_64-pc-linux-gnu-library/4.1')

###############################################################################

# Set up libraries and options/files

###############################################################################

library(data.table)
library(dplyr)
library(optparse)
library(readr)
library(tidyr)

parse <- OptionParser()

# setting up options for the filepaths to the correct files
option_list <- list(
  make_option('--cohort', type='character', help="Cohort, ideally no spaces (for graphs and documentation)", action='store'),
  make_option('--id_column', type = 'character', default="IID", help = "Column names of identifier column in phenotype and covariate files", action = 'store'),
  make_option('--mrs', type = 'character', help = 'File path to antidepressant exposure MRS file (made using MRS_calc.R), column names (IID, AD_MRS and (optional) FID'),
  make_option('--pheno', type = 'character', help = 'File path to antidepressant exposure phenotype file, column names (IID, antidep and (optional) FID'),
  make_option('--demo', type = 'character', help = 'File path to demographics file, column names (age, sex_coded, bmi, ever_smoke and (optional) mdd'),
  make_option('--outdir', type = 'character', help = 'The filepath for output directory', action = 'store')
)

args = commandArgs(trailingOnly=TRUE)
opt <- parse_args(OptionParser(option_list=option_list), args=args)
cohort <- opt$cohort
id_col <- opt$id_column # Vector of identifier columns 
pheno_fp=opt$pheno # AD exposure (phenotype of cohort)
MRS_fp=opt$mrs # AD MRS (predictor) 
demo_fp=opt$demo # Demographics 
outdir <- opt$outdir # File path of output directory

# sinking all output to a log file 

sink(paste0(outdir, cohort, "_MRS_demographics.log"))


###############################################################################

# Read in files

###############################################################################

# read in file with additional demographic information for manuscript 

demographics <- read.table(demo_fp, header = T)
MRS <- read.table(MRS_fp, header = T)

# support phenotype .csv files or .txt files 

if (endsWith(pheno_fp, '.csv')){
  ad_pheno <- read.csv(pheno_fp, header = T)
} else if (endsWith(pheno_fp, '.txt')) {
  ad_pheno <- read.table(pheno_fp, header = T)
} else {
  stop('Unsupported phenotype file, please provide the phenotype as a .csv or .txt file')
}


###############################################################################

# File checks 

###############################################################################

if('AD_MRS' %in% colnames(MRS) == FALSE){
  stop('No AD_MRS column in the MRS file')
} else {
  print('AD_MRS column in the MRS file')
}


# check that there is an antidep column in the file 

if('antidep' %in% colnames(ad_pheno) == FALSE){
  stop('No antidep column in the phenotype file, please change name')
} else {
  print('antidep column in the phenotype file')
}

# filter out any missing phenotype values (if any?)
print('Missing values in the antidep phenotype? ')
table(is.na(ad_pheno$antidep))

ad_pheno <- ad_pheno %>% filter(!is.na(antidep))

# Demographics file 
# age, sex, bmi, smoking behaviours and MDD (optional/if applicable)
# load in vector of required demographic variables 

req_demo_vars <- c('age', 'sex_coded', 'bmi', 'ever_smoke')
if(all(req_demo_vars %in% colnames(demographics))){
  print('All demographic variables loaded in and named correctly')
} else {
    print('Not all demographic variables present / not loaded correctly')
    cols_missing = req_demo_vars[!(req_demo_vars %in% colnames(demographics))]
    print(paste0('Incorrect/missing columns: ', cols_missing))
}

# smoking 
# recode the smoking as current, former and never smoker 

# Question: Have you ever smoked tobacco? 
# 1- Yes Current
# 2 - Yes but stopped in last 12 months 
# 3 - Yes but stopped more than 12 months ago
# 4 No Never smoked 
# 5 (NA values? )


print('Smoking variable')
table(demographics$ever_smoke)

if (!("character" %in% class(demographics$ever_smoke))) {
print('Formatting the smoking variable to Current/Former/Never Smokers')
demographics <- demographics %>% 
  mutate(smoking_var = case_when(
    ever_smoke == 1 ~ 'Current', 
    ever_smoke == 2 ~ 'Former', 
    ever_smoke == 3 ~ 'Former', 
    ever_smoke == 4 ~ 'Never', 
    TRUE ~ NA_character_))

print('Smoking variable named')
table(demographics$smoking_var)

} else {
  print('Smoking variable is already formatted')
  table(demographics$smoking_var)
}

# sex 

print('Sex variable')
table(demographics$sex_coded)
demographics <- demographics %>%
  mutate(sex_name = ifelse(
    sex_coded == 0, 'Female', 'Male'))
print('Following coding of sex variable, 0 = Female and 1 = Male')
table(demographics$sex_name)

# merge with the MRS file 

demographics <- merge(demographics, MRS, by= id_col, all = TRUE)

###############################################################################

# Summarising demographics function 

###############################################################################

# create a summary of this phenotype 
# summarise function with the continuous covs 
# ratios for the discrete covs (sex, smoking)

# summary function for cases and controls 
# phenotype is the phenotype column name in the phenotype file 
# file is the dataframe with the phenotype and covariate information 


phenotype_summary <- function(phenotype, file) {
  file <- file[!is.na(file[,phenotype]),]
  print(colnames(file))
  summary <- file %>% 
    group_by(across(all_of(phenotype)))  %>% 
    summarise(age = paste0(signif(mean(age, na.rm = T),3), 
                           ' (', signif(sd(age, na.rm = T),3), ')'), 
              bmi = paste0(signif(mean(bmi, na.rm = T),3), ' (', 
                           signif(sd(bmi, na.rm = T),3),')'), 
              AD_MRS= paste0(signif(mean(AD_MRS, na.rm = T),3),' (', 
                                signif(sd(AD_MRS, na.rm = T),3), ')'), 
              
              n = n()) %>%
    
    mutate(current_smoking= c(file %>% group_by(across(all_of(c(phenotype, 'smoking_var')))) %>%
                                summarise(count = n()) %>%
                                mutate(percentage = (count/sum(count))*100) %>%
                                mutate(val = paste0(count, ' (', signif(percentage,2), '%)')) %>%
                                filter(!!sym(phenotype) == 0 & smoking_var == 'Current') %>% 
                                pull(val),
                              file %>% group_by(across(all_of(c(phenotype, 'smoking_var')))) %>% 
                                summarise(count = n()) %>%
                                mutate(percentage = (count/sum(count))*100) %>%
                                mutate(val = paste0(count, ' (', signif(percentage,2), '%)')) %>%
                                filter(!!sym(phenotype) == 1 & smoking_var == 'Current') %>% 
                                pull(val)), 
           
           former_smoking= c(file %>% group_by(across(all_of(c(phenotype, 'smoking_var')))) %>% 
                               summarise(count = n()) %>%
                               mutate(percentage = (count/sum(count))*100) %>%
                               mutate(val = paste0(count, ' (', signif(percentage,2), '%)')) %>%
                               filter(!!sym(phenotype) == 0 & smoking_var == 'Former') %>% 
                               pull(val),
                             file %>% group_by(across(all_of(c(phenotype, 'smoking_var')))) %>% 
                               summarise(count = n()) %>%
                               mutate(percentage = (count/sum(count))*100) %>%
                               mutate(val = paste0(count, ' (', signif(percentage,2), '%)')) %>%
                               filter(!!sym(phenotype) == 1 & smoking_var == 'Former') %>% 
                               pull(val)) , 
           
           never_smoking= c(file %>% group_by(across(all_of(c(phenotype, 'smoking_var')))) %>% 
                              summarise(count = n()) %>%
                              mutate(percentage = (count/sum(count))*100) %>%
                              mutate(val = paste0(count, ' (', signif(percentage,2), '%)')) %>%
                              filter(!!sym(phenotype) == 0 & smoking_var == 'Never') %>% 
                              pull(val),
                            file %>% group_by(across(all_of(c(phenotype, 'smoking_var')))) %>% 
                              summarise(count = n()) %>%
                              mutate(percentage = (count/sum(count))*100) %>%
                              mutate(val = paste0(count, ' (', signif(percentage,2), '%)')) %>%
                              filter(!!sym(phenotype) == 1 & smoking_var == 'Never') %>% 
                              pull(val)), 
           
           Female = c(file %>% 
                        group_by(across(all_of(c(phenotype, 'sex_name')))) %>% 
                        summarise(count = n()) %>% 
                        mutate(percentage = (count/sum(count))*100) %>%
                        mutate(val = paste0(count, ' (', signif(percentage,2), '%)')) %>%
                        filter(!!sym(phenotype)==0 & sex_name == 'Female') %>% 
                        pull(val), 
                      file %>% group_by(across(all_of(c(phenotype, 'sex_name')))) %>% 
                        summarise(count = n()) %>% 
                        mutate(percentage = (count/sum(count))*100) %>%
                        mutate(val = paste0(count, ' (', signif(percentage,2), '%)')) %>%
                        filter(!!sym(phenotype)==1 & sex_name == 'Female') %>% 
                        pull(val)), 
           Male = c(file %>% 
                      group_by(across(all_of(c(phenotype, 'sex_name')))) %>% 
                      summarise(count = n()) %>% 
                      mutate(percentage = (count/sum(count))*100) %>%
                      mutate(val = paste0(count, ' (', signif(percentage,2), '%)')) %>%
                      filter(!!sym(phenotype)==0 & sex_name == 'Male') %>% 
                      pull(val), 
                    file %>% group_by(across(all_of(c(phenotype, 'sex_name')))) %>% 
                      summarise(count = n()) %>% 
                      mutate(percentage = (count/sum(count))*100) %>%
                      mutate(val = paste0(count, ' (', signif(percentage,2), '%)')) %>%
                      filter(!!sym(phenotype)==1 & sex_name == 'Male') %>% 
                      pull(val))) %>% 
    as.data.frame()
  
  summary <- summary %>% select(phenotype, n, age ,bmi,  Female, Male, current_smoking, former_smoking, never_smoking, packyears, AD_MRS)
  colnames(summary) <- c(phenotype, 'N', 'Age (%)','BMI (%)', 'Sex Female', 'Sex Male', 'Current smoker', 'Former smoker', 'Never smoker','Pack years', 'AD MRS (%)')
  return(summary)
}

###############################################################################

# Merging with phenotype and summarising

###############################################################################

demographics_pheno <- merge(ad_pheno, demographics, by = id_col, all = TRUE)
demo_summary <- phenotype_summary(phenotype = 'antidep', file = demographics_pheno)

# add MDD summary if that info is available 

if('mdd' %in% colnames(demographics_pheno)){
  print('MDD data available for the phenotype')
  # summarise the MDD data 
  # print the table of mdd variable to log (to sanity check the conversion)
  print('The mdd variable:')
  table(demographics_pheno$mdd)
  demographics_pheno <- demographics_pheno %>%
    mutate(mdd_pheno = case_when(mdd == 0 ~ 'Control', mdd==1 ~ 'Case', TRUE ~ NA))
  
  print('The mdd named variable:')
  table(demographics_pheno$mdd_pheno)
  MDD_summary <- data.table(antidep = c(0, 1),
    Cases = c(demographics_pheno %>%
                             group_by(across(all_of(c('antidep', 'mdd_pheno')))) %>% 
                             summarise(count = n()) %>% 
                             mutate(percentage = (count/sum(count))*100) %>%
                             mutate(val = paste0(count, ' (', signif(percentage,2), '%)')) %>%
                             filter(antidep==0 & mdd_pheno == 'Case') %>% 
                             pull(val), 
              demographics_pheno  %>%
                             group_by(across(all_of(c('antidep', 'mdd_pheno')))) %>% 
                             summarise(count = n()) %>% 
                             mutate(percentage = (count/sum(count))*100) %>%
                             mutate(val = paste0(count, ' (', signif(percentage,2), '%)')) %>%
                             filter(antidep==1 & mdd_pheno == 'Case') %>% 
                             pull(val)),
  Controls = c(demographics_pheno  %>%
                 group_by(across(all_of(c('antidep', 'mdd_pheno')))) %>% 
                 summarise(count = n()) %>% 
                 mutate(percentage = (count/sum(count))*100) %>%
                 mutate(val = paste0(count, ' (', signif(percentage,2), '%)')) %>%
                 filter(antidep==0 & mdd_pheno == 'Control') %>% 
                 pull(val), 
               demographics_pheno  %>%
                 group_by(across(all_of(c('antidep', 'mdd_pheno')))) %>% 
                 summarise(count = n()) %>% 
                 mutate(percentage = (count/sum(count))*100) %>%
                 mutate(val = paste0(count, ' (', signif(percentage,2), '%)')) %>%
                 filter(antidep==1 & mdd_pheno == 'Control') %>% 
                 pull(val)
  )) %>% as.data.frame()
  colnames(MDD_summary) <- c('antidep','MDD cases', 'MDD controls')

  demo_summary <- merge(demo_summary, MDD_summary, by = 'antidep')

} else{
  print('No MDD phenotype data available')
}

###############################################################################

# Save the demographic summary 

###############################################################################

write.table(demo_summary, paste0(outdir, cohort, '_demo_summary.txt'), quote = F, row.names = F)
print(paste0('Saved the demographic summary to ', cohort, '_demo_summary.txt'))

sink()





