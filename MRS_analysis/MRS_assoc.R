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
library(tibble)
library(pROC)

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

# logging phenotype characteristics 
print(paste0('Read in the Antidepressant exposure phenotype for ', cohort, ' : Number of cases: ',
             nrow(ad_pheno %>% 
                    filter(antidep==1)), 
             'Number of controls: ',
             nrow(ad_pheno%>% 
                    filter(antidep==0))))

all_covs <- read.table(covs_fp, header = T) # 'GS_test_covs.txt', covariate file which has PCs in it

print(paste0("Covariates read in ", paste(colnames(all_covs %>% dplyr::select(-all_of(id_col))), collapse = ", ")))

#merge the phenotype, MRS and the covariate file together 

MRS_covs <- merge(MRS, all_covs, by = id_col)
MRS_covs_pheno <- merge(MRS_covs, ad_pheno, by = id_col)

# logging phenotype characteristics after merging 

print(paste0('Read in the Antidepressant exposure phenotype for ', cohort, ' after merging with MRS and covariates: Number of cases: ',
             nrow(MRS_covs_pheno %>% 
                    filter(antidep==1)), 
             ' \n Number of controls: ',
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
warnings()
print(assoc_mod)
print(summary(assoc_mod)$coefficients %>% as.data.frame())

assoc_coefs <- summary(assoc_mod)$coefficients %>% as.data.frame()
assoc_coefs <- rownames_to_column(assoc_coefs, var = "Coefficient")

# save the coefficients 
print(paste0('Saving the full model coefficients to ', outdir, cohort, "_MRS_AD_coefficients.txt"))

write.table(assoc_coefs, paste0(outdir, cohort, "_MRS_AD_coefficients.txt"), row.names = F, quote = F)

###############################################################################

# Calculating the AUC and ROC graph

###############################################################################

predicted_probs <- predict(assoc_mod, type = 'response')

# takes the true outcomes (MRS_covs_pheno$antidep)
# and the predicted probabilities for the 1 ('case') class
# returns false positive and true positive rates for different classification thresholds

roc_curve <- roc(MRS_covs_pheno$antidep, predicted_probs)
auc_value <- auc(roc_curve)

# save ROC curve object for plotting all cohorts together
print(paste0('Saving the ROC curve object for plotting all cohorts together to rds object: ', outdir, cohort, '_roc_curve.rds'))
saveRDS(roc_curve, paste0(outdir, cohort, '_roc_curve.rds'))

# ROC Graph 
print(paste0('Saving the ROC curve for the cohort alone to ', outdir, cohort, '_assoc_ROC_curve.pdf'))
cairo_pdf(file = paste0(outdir, cohort, '_assoc_ROC_curve.pdf'), width = 8, height = 6)
plot.roc(roc_curve, col = "blue", lwd =2, main = paste0('ROC Curve: ', cohort))
dev.off()

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

warnings()

print('Null model summary')
print(null_mod)

print('Null model coefficients')
print(summary(null_mod)$coefficients %>% as.data.frame())

# calculate McFaddens R square 

mcf_r2 <- 1-logLik(assoc_mod)/logLik(null_mod)
print(paste0('McFaddens R2:', mcf_r2))

# calculate Nagelkerke's R square 
# Nagelkerke's R square 
# calculate the likelihoods of the full and null models

lik_assoc_mod <- exp(as.numeric(logLik(assoc_mod)))
lik_null_mod <- exp(as.numeric(logLik(null_mod)))

# calculate cox and snell R2 (the numerator in the Nagelkerkes R2)

cox_snell_r2 <- 1-(lik_null_mod/lik_assoc_mod)^(2/nobs(assoc_mod))

# calculate nagelkerkes r2
nagel_denominator <- 1-(lik_null_mod)^(2/nobs(assoc_mod))
nagel_r2 <- cox_snell_r2/nagel_denominator

# create a table for results 
  
modelmetrics <- data.frame(Cohort=cohort, 
                               loglik_MRS = logLik(assoc_mod), 
                               loglik_null = logLik(null_mod),
                               lik_MRS = exp(logLik(assoc_mod)),
                               lik_null = exp(logLik(null_mod)),
                               mcfad_R2 = mcf_r2, 
                               cox_snell_r2 = cox_snell_r2,
                               nagelkerke_r2 = nagel_r2,
                               AUC = auc_value)



print('Model metric table')
print(modelmetrics)

# save R2 and model metric table 
  
print(paste0('Saving the model metric table to ',outdir, cohort, "_MRS_AD_logL.txt"))
write.table(modelmetrics, paste0(outdir, cohort, "_MRS_AD_modmetrics.txt"), row.names = F, quote = F)
  



