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
library(caret)

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
summary(MRS_covs_pheno)
MRS_covs_pheno <- na.omit(MRS_covs_pheno) #6NAs in the cell count columns 
summary(MRS_covs_pheno)
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

assoc_mod <- glmer(as.factor(antidep)~ scale(AD_MRS) + scale(age) + scale(neut) + 
                     scale(lymph) + scale(cg05575921)+ scale(mono) + scale(eosin) + 
                     scale(baso) + 
                      as.factor(sex_coded)  +
                     (1|array), data = MRS_covs_pheno, family = 'binomial')

## When get genetic PCs -- add these ! 
# + scale(C1) + scale(C2) + 
#scale(C3) + scale(C4) + scale(C5) + scale(C6) +
 # scale(C7) + scale(C8) + scale(C9) + scale(C10)


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

# Confusion matrix 

###############################################################################

binary_predictions <- ifelse(predicted_probs > 0.5, 1, 0)
table(MRS_covs_pheno$antidep, binary_predictions)
conf_matrix <- confusionMatrix(data = factor(binary_predictions, levels = c(0,1)), reference = factor(MRS_covs_pheno$antidep, levels = c(0,1)), positive = "1")
#https://stackoverflow.com/questions/23891140/r-how-to-visualize-confusion-matrix-using-the-caret-package
draw_confusion_matrix <- function(cm) {
  
  total <- sum(cm$table)
  res <- as.numeric(cm$table)
  
  # Generate color gradients. Palettes come from RColorBrewer.
  greenPalette <- c("#F7FCF5","#E5F5E0","#C7E9C0","#A1D99B","#74C476","#41AB5D","#238B45","#006D2C","#00441B")
  redPalette <- c("#FFF5F0","#FEE0D2","#FCBBA1","#FC9272","#FB6A4A","#EF3B2C","#CB181D","#A50F15","#67000D")
  getColor <- function (greenOrRed = "green", amount = 0) {
    if (amount == 0)
      return("#FFFFFF")
    palette <- greenPalette
    if (greenOrRed == "red")
      palette <- redPalette
    colorRampPalette(palette)(100)[10 + ceiling(90 * amount / total)]
  }
  
  # set the basic layout
  layout(matrix(c(1,1,2)))
  par(mar=c(2,2,2,2))
  plot(c(100, 345), c(300, 450), type = "n", xlab="", ylab="", xaxt='n', yaxt='n')
  title(paste0(cohort, ' : CONFUSION MATRIX'), cex.main=2)
  
  
  # create the matrix 
  classes = colnames(cm$table)
  rect(150, 430, 240, 370, col=getColor("green", res[1]))
  text(195, 435, classes[1], cex=1.2)
  rect(250, 430, 340, 370, col=getColor("red", res[3]))
  text(295, 435, classes[2], cex=1.2)
  text(125, 370, 'Predicted', cex=1.3, srt=90, font=2)
  text(245, 450, 'Actual', cex=1.3, font=2)
  rect(150, 305, 240, 365, col=getColor("red", res[2]))
  rect(250, 305, 340, 365, col=getColor("green", res[4]))
  text(140, 400, classes[1], cex=1.2, srt=90)
  text(140, 335, classes[2], cex=1.2, srt=90)
  
  # add in the cm results
  text(195, 400, res[1], cex=1.6, font=2, col='white')
  text(195, 335, res[2], cex=1.6, font=2, col='white')
  text(295, 400, res[3], cex=1.6, font=2, col='white')
  text(295, 335, res[4], cex=1.6, font=2, col='white')
  
  # add in the specifics 
  plot(c(100, 0), c(100, 0), type = "n", xlab="", ylab="", main = "DETAILS", xaxt='n', yaxt='n')
  text(10, 85, names(cm$byClass[1]), cex=1.2, font=2)
  text(10, 70, round(as.numeric(cm$byClass[1]), 3), cex=1.2)
  text(30, 85, names(cm$byClass[2]), cex=1.2, font=2)
  text(30, 70, round(as.numeric(cm$byClass[2]), 3), cex=1.2)
  text(50, 85, names(cm$byClass[5]), cex=1.2, font=2)
  text(50, 70, round(as.numeric(cm$byClass[5]), 3), cex=1.2)
  text(70, 85, names(cm$byClass[6]), cex=1.2, font=2)
  text(70, 70, round(as.numeric(cm$byClass[6]), 3), cex=1.2)
  text(90, 85, names(cm$byClass[7]), cex=1.2, font=2)
  text(90, 70, round(as.numeric(cm$byClass[7]), 3), cex=1.2)
  
  # add in the accuracy information 
  text(20, 35, names(cm$overall[1]), cex=1.5, font=2)
  text(20, 20, round(as.numeric(cm$overall[1]), 3), cex=1.4)
  text(50, 35, names(cm$overall[2]), cex=1.5, font=2)
  text(50, 20, round(as.numeric(cm$overall[2]), 3), cex=1.4)
  text(80, 35, names(cm$byClass[11]), cex=1.5, font=2)
  text(80, 20, round(as.numeric(cm$byClass[11]), 3), cex=1.4)
}

cairo_pdf(file = paste0(outdir, cohort, '_assoc_conf_matrix.pdf'), width = 8, height = 6)
draw_confusion_matrix(conf_matrix)
dev.off()

###############################################################################

# Precision Recall Curve 

###############################################################################

pr_curve <- pr.curve(MRS_covs_pheno$antidep, predicted_probs, curve = T)

cairo_pdf(file = paste0(outdir, cohort, '_assoc_precision_recall.pdf'), width = 8, height = 6)
plot(pr_curve, main= paste0(cohort, ' : Precision Recall Curve'), col = 'red')
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

null_mod <- glmer(as.factor(antidep)~ scale(age) + scale(neut) + 
                    scale(lymph) + scale(cg05575921)+ scale(mono) + scale(eosin) + 
                    scale(baso) + 
                    as.factor(sex_coded)  +
                    (1|array), data = MRS_covs_pheno, family = 'binomial')

## When get genetic PCs -- add these ! 
# + scale(C1) + scale(C2) + 
#scale(C3) + scale(C4) + scale(C5) + scale(C6) +
# scale(C7) + scale(C8) + scale(C9) + scale(C10)
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
  



