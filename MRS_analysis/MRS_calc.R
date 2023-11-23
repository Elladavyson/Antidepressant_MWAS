###############################################################################

# Set up libraries and options/files

###############################################################################

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
  make_option('--weights', type = 'character', default = 'GS_AD_MRS_weights.txt', help = "Weights file provided for calculating MRS"),
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
outdir <- opt$outdir

# sink output to a text file with the same name as outpath with a '.log' extension
sink(paste0(outdir, cohort, "_AD_MRS.log"))
print(paste0('Calculating a MRS for ', cohort))
print(paste0('Read in the processed DNAm file from: ', DNAm_fp))
print(paste0('Read in the weights file from: ', weights_fp))
print(paste0('Output to be saved in: ', outdir))

###############################################################################

# Read in the files

###############################################################################

DNAm <- read_table(DNAm_fp) %>% as.data.frame()
print(paste0('The DNAm file has data for ', nrow(DNAm), ' participants and ', DNAm %>% select(starts_with("cg")) %>% ncol(), " CpG sites"))
weights <- read.table(weights_fp, header = T) 
print(paste0('The weights file has weights for ', nrow(weights), ' CpGs'))

if(DNAm %>% select(starts_with("cg")) %>% ncol() != nrow(weights)){
  print(paste0('Number of probes read in for DNAm file (', DNAm %>% select(starts_with("cg")) %>% ncol(),
  ') does not match the number of weights provided (', nrow(weights), ')'))
} else {
  print(paste0('Number of probes read in for DNAm matches the number of weights provided: n = ', nrow(weights)))
}
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

# Group by ID
# and then calculate a weighted sum of all the CpGs per IID

MRS <- DNAm_long %>% group_by(!!sym(id_col)) %>%
  summarise(weighted_sum = sum(Weight*Mval, na.rm = T)) %>% as.data.frame()

print(paste0('MRS calculated for ', nrow(MRS), ' people in the ', cohort, ' cohort'))

###############################################################################

# MRS Distributions

###############################################################################

print(paste0('Plotting distributions across the whole ', cohort, ' sample'))

# Distribution of the MRS across the DNAm cohort 

MRS_dist <- ggplot(MRS, aes(x = weighted_sum)) + 
  geom_histogram() + 
  theme_minimal() + 
  labs(x = 'Methylation Risk Score', y = 'Count')+
  ggtitle(cohort)

ggsave(paste0(outdir, cohort, "_AD_MRS_overalldist.png"), MRS_dist, width = 8, height = 6, device='png', dpi=300)

# looking at Distribution in AD exposed and AD not exposed (violin plots)

if (endsWith(pheno_fp, '.csv')){
ad_pheno <- read.csv(pheno_fp, header = T)
} else if (endsWith(pheno_fp, '.txt')) {
  ad_pheno <- read.table(pheno_fp, header = T)
} else {
  stop('Unsupported phenotype file, please provide the phenotype as a .csv or .txt file')
}

if('antidep' %in% colnames(ad_pheno) == FALSE){
  stop('No antidep column in the phenotype file')
} else {
  print('antidep column in the phenotype file')
}

ad_pheno <- ad_pheno %>% filter(!is.na(antidep)) # remove missing values if there are any 

print(paste0('Read in the Antidepressant exposure phenotype for ', cohort, ': Number of cases: ',
             nrow(ad_pheno %>% 
                    filter(antidep==1)), 
             ' Number of controls: ',
             nrow(ad_pheno%>% 
                    filter(antidep==0))))
ad_pheno_MRS <- merge(ad_pheno, MRS, by = 'IID')

print('Plotting MRS distributions for AD exposure cases and controls ')
MRS_pheno_dists <- ggplot(ad_pheno_MRS, aes(x = weighted_sum, fill = as.factor(antidep))) + 
  geom_histogram(alpha = 0.8) + 
  theme_minimal() + 
  labs(x = 'Methylation Risk Score', y = 'Count', fill = 'Self-reported AD use') +
  ggtitle(cohort)

ggsave(paste0(outdir, cohort, "_AD_MRS_phenodist.png"), MRS_pheno_dists, width = 8, height = 6, device='png', dpi=300)

###############################################################################

# Saving the methylation risk score 

###############################################################################
outfile <- paste0(outdir, cohort, '_AD_MRS.txt')
print(paste0('Saving the methylation risk score to ', outfile))

colnames(MRS)[2] <- 'AD_MRS'
write.table(MRS, outfile, row.names = F, quote = F)


