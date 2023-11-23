###############################################################################

# Set up libraries and options/files

###############################################################################

.libPaths('/exports/igmm/eddie/GenScotDepression/users/edavyson/R/x86_64-pc-linux-gnu-library/4.1')
library(data.table)
library(dplyr)
library(optparse)
library(readr)
library(tidyr)
library(ggplot2)
library(tools)
library(lme4)

parse <- OptionParser()

# setting up options for the filepaths to the correct files
option_list <- list(
  make_option('--cohort', type='character', help="Cohort, ideally no spaces (for graphs and documentation)", action='store'),
  make_option('--id_column', type = 'character', default="IID", help = "Column names of identifier column in phenotype and covariate files", action = 'store'),
  make_option('--mrs', type = 'character', help = 'File path to antidepressant exposure MRS file (made using MRS_calc.R)'),
  make_option('--pheno', type = 'character', help = 'File path to antidepressant exposure phenotype file'),
  make_option('--covs', type = 'character', help = 'File path to covariate file'),
  make_option('--outdir', type = 'character', help = 'The filepath for output directory', action = 'store')
)

args = commandArgs(trailingOnly=TRUE)
opt <- parse_args(OptionParser(option_list=option_list), args=args)
cohort <- opt$cohort
id_col <- opt$id_column # Vector of identifier columns 
pheno_fp=opt$pheno # AD exposure (phenotype of cohort)
MRS_fp=opt$mrs # AD MRS (predictor) 
covs_fp=opt$covs # Covariates 'GS_test_covs.txt'
outdir <- opt$outdir # File path of output directory

# sinking all output to a log file 

sink(paste0(outdir, cohort, "_MRS_AD_assoc.log"))

###############################################################################

# Read in covariates, phenotype and MRS files 

###############################################################################

MRS <- read.table(MRS_fp, header = T) # 'GS_AD_MRS.txt', file with IID and MRS
print(colnames(MRS))

# check that there is a MRS column in the file 

if('AD_MRS' %in% colnames(MRS) == FALSE){
  stop('No AD_MRS column in the MRS file')
} else {
  print('AD_MRS column in the MRS file')
}

# support phenotype .csv files or .txt files 

if (endsWith(pheno_fp, '.csv')){
  ad_pheno <- read.csv(pheno_fp, header = T)
} else if (endsWith(pheno_fp, '.txt')) {
  ad_pheno <- read.table(pheno_fp, header = T)
} else {
  stop('Unsupported phenotype file, please provide the phenotype as a .csv or .txt file')
}

# check that there is an antidep column in the file 

if('antidep' %in% colnames(ad_pheno) == FALSE){
  stop('No antidep column in the phenotype file')
} else {
  print('antidep column in the phenotype file')
}

# remove missing values (if any?)
ad_pheno <- ad_pheno %>% filter(!is.na(antidep)) 
print(head(ad_pheno))
# logging phenotype characteristics 
print(paste0('Read in the Antidepressant exposure phenotype for ', cohort, ' : Number of cases: ',
             nrow(ad_pheno %>% 
                    filter(antidep==1)), 
             'Number of controls: ',
             nrow(ad_pheno%>% 
                    filter(antidep==0))))

all_covs <- read.table(covs_fp, header = T) # 'GS_test_covs.txt', covariate file which has PCs in it

print(paste0("Covariates read in ", paste(colnames(all_covs %>% select(-all_of(id_col))), collapse = ", ")))

#merge the phenotype, MRS and the covariate file together 

MRS_covs <- merge(MRS, all_covs, by = id_col)
MRS_covs_pheno <- merge(MRS_covs, ad_pheno, by = id_col)

# logging phenotype characteristics after merging 

print(paste0('Read in the Antidepressant exposure phenotype for ', cohort, ' after merging with MRS and covariates: Number of cases: ',
             nrow(MRS_covs_pheno %>% 
                    filter(antidep==1)), 
             'Number of controls: ',
             nrow(MRS_covs_pheno%>% 
                    filter(antidep==0))))


###############################################################################
  
  # Generalised linear mixed model (GLMM)
  
###############################################################################

# Fit a logistic mixed effects model 
# logit link function and binomial family 
# Outcome - Antidepressant exposure phenotype 
# Predictor - Antidepressant MRS (from MRS_calc.R)
# Covariates - Age, sex, Monocyte cell proportions, lymphocyte cell proportions,
# AHRR methylation M values + top 10 Genetic principal components 

assoc_mod <- glmer(as.factor(antidep)~ scale(AD_MRS) + scale(age) + scale(Mono) + 
                     scale(lymphocytes) + scale(cg05575921)+
                      as.factor(sex_coded) + scale(C1) + scale(C2) + 
                     scale(C3) + scale(C4) + scale(C5) + scale(C6) +
                     scale(C7) + scale(C8) + scale(C9) + scale(C10) +
                     (1|Batch), data = MRS_covs_pheno, family = 'binomial')


# Extract the fixed effect estimates, standard errors and p-value 

print(assoc_mod)
print(summary(assoc_mod)$coefficients %>% as.data.frame())

assoc_coefs <- summary(assoc_mod)$coefficients %>% as.data.frame()
  
# save the coefficients 
print(paste0('Saving the full model coefficients to ', outdir, cohort, "_MRS_AD_coefficients.txt"))

write.table(assoc_coefs, paste0(outdir, cohort, "_MRS_AD_coefficients.txt"), row.names = F, quote = F)


###############################################################################

# Null model and McFaddens R-squared

###############################################################################

# Fit a model the same but without the MRS to calculate McFaddens R-squared  
# Fit a logistic mixed effects model 
# logit link function and binomial family 
# Outcome - Antidepressant exposure phenotype 
# Predictor- NULL
# Covariates - Age, sex, Monocyte cell proportions, lymphocyte cell proportions,
# AHRR methylation M values + top 10 Genetic principal components 

print('Fitting a Null model')

null_mod <- glmer(as.factor(antidep)~ scale(age) + scale(Mono) + scale(lymphocytes) + scale(cg05575921)+
                      as.factor(sex_coded) + scale(C1) + scale(C2) + 
                    scale(C3) + scale(C4) + scale(C5) + scale(C6) +
                    scale(C7) + scale(C8) + scale(C9) + scale(C10) +
                    (1|Batch), data = MRS_covs_pheno, family = 'binomial')


print('Null model summary')
print(null_mod)
print('Null model coefficients')
print(summary(null_mod)$coefficients %>% as.data.frame())
# calculate McFaddens R square 

mcf_r2 <- 1-logLik(assoc_mod)/logLik(null_mod)
print(paste0('McFaddens R2:', mcf_r2))

# create a table for results 
  
loglikelihoods <- data.frame(Cohort=cohort, 
                               loglik_MRS = logLik(assoc_mod), 
                               loglik_null = logLik(null_mod),
                               mcfad_R2 = mcf_r2)

print('Log likelihood table of the full and null model')
print(loglikelihoods)
# save R2 and log likelihood values 
  
print(paste0('Saving the log likelihood table to ',outdir, cohort, "_MRS_AD_logL.txt"))
write.table(loglikelihoods, paste0(outdir, cohort, "_MRS_AD_logL.txt"), row.names = F, quote = F)
  



