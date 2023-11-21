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

parse <- OptionParser()

# setting up options for the filepaths to the correct files
option_list <- list(
  make_option('--cohort', type='character', help="Cohort, ideally no spaces (for graphs and documentation)", action='store'),
  make_option('--DNAm', type='character', help="The filepath for processed DNAm file", action='store'),
  make_option('--id_column', type = 'character', default="IID", help = "Column names of identifier column in DNAm object", action = 'store'),
  make_option('--weights', type = 'character', default = 'big_lasso_450K_selfrep.txt', help = "Weights file provided for calculating MRS"),
  make_option('--pheno', type = 'character', help = 'File path to antidepressant exposure phenotype file'),
  make_option('--outdir', type = 'character', help = 'The filepath for output directory', action = 'store')
)

args = commandArgs(trailingOnly=TRUE)
opt <- parse_args(OptionParser(option_list=option_list), args=args)
cohort <- opt$cohort
DNAm_fp=opt$DNAm # DNAm file
weights_fp=opt$weights # Weights file from GS
id_col <- opt$id_column # Vector of identifier columns 
pheno_fp=opt$pheno
outdir <- paste0(cohort,"_", opt$outdir) # output file path

# sink output to a text file with the same name as outpath with a '.log' extension
sink(paste0(out_dir, cohort, "_AD_MRS.log"))
print(paste0('Calculating a MRS for ', cohort))
print(paste0('Read in the processed DNAm file from: ', DNAm_fp))
print(paste0('Read in the weights file from: ', weights_fp))
print(paste0('Output to be saved in: ', outdir))

###############################################################################

# Read in the files

###############################################################################

DNAm <- read_table(DNAm_fp) %>% as.data.frame()
weights <- read.table(weights_fp, header = F) 
colnames(weights) <- c('CpG', 'Weight')
weights <- weights[-1,]

###############################################################################

# Calculate the MRS 

###############################################################################
print('Converting DNAm to long format')
DNAm_long <- DNAm %>% 
  pivot_longer(-c(all_of(id_col)), 
               names_to = "CpG", 
               values_to= "Mval")
print('Merging with the weights')
DNAm_long <- merge(DNAm_long, weights, by = 'CpG')
print('Calculating the MRS')
MRS <- DNAm_long %>% group_by(!!sym(id_col)) %>%
  summarise(weighted_sum = sum(Weight*Mval, na.rm = T)) %>% as.data.frame()


###############################################################################

# MRS Distributions

###############################################################################

print('Plotting distributions')

# Distribution of the MRS across the DNAm cohort 

MRS_dist <- ggplot(MRS, aes(x = weighted_sum)) + 
  geom_histogram() + 
  theme_minimal() + 
  labs(x = 'Methylation Risk Score', y = 'Count')+
  ggtitle(cohort)

ggsave(paste0(out_dir, cohort, "_AD_MRS_overalldist.png"), MRS_dist, width = 8, height = 6, device='png', dpi=300)

# looking at Distribution in AD exposed and AD not exposed (violin plots)

if (endsWith(pheno_fp, '.csv')){
ad_pheno <- read.csv(pheno_fp, header = T)
} else if (endsWith(pheno_fp, '.txt')) {
  ad_pheno <- read.table(pheno_fp, header = T)
} else {
  stop('Unsupported phenotype file, please provide the phenotype as a .csv or .txt file')
}
ad_pheno <- ad_pheno %>% filter(!is.na(antidep)) # remove missing values if there are any 
ad_pheno_MRS <- merge(ad_pheno, MRS, by = 'IID')

MRS_pheno_dists <- ggplot(ad_pheno_MRS, aes(x = weighted_sum, fill = as.factor(antidep))) + 
  geom_histogram(alpha = 0.8) + 
  theme_minimal() + 
  labs(x = 'Methylation Risk Score', y = 'Count', fill = 'Self-reported AD use') +
  ggtitle(cohort)

ggsave(paste0(out_dir, cohort, "_AD_MRS_phenodist.png"), MRS_pheno_dists, width = 8, height = 6, device='png', dpi=300)

###############################################################################

# Saving the methylation risk score 

###############################################################################
outfile <- paste0(outdir, cohort, '_AD_MRS.txt')
print(paste0('Saving the methylation risk score to ', outfile))

colnames(MRS)[2] <- 'AD_MRS'
write.table(MRS, outfile, row.names = F, quote = F)


