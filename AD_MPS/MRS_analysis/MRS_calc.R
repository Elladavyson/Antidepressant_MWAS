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
library(ggpubr)

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


###############################################################################

# Assess missingness of CpGs

###############################################################################

if(DNAm %>% select(starts_with("cg")) %>% ncol() != nrow(weights)){
  print(paste0('Number of probes read in for DNAm file (', DNAm %>% select(starts_with("cg")) %>% ncol(),
               ') does not match the number of weights provided (', nrow(weights), ')'))
  # If not all weights included, save the probes included in the MRS for the cohort 
  readr::write_lines(DNAm %>% select(starts_with("cg")) %>% colnames(), paste0(outdir,cohort, '_probesinMRS.txt'))
  print(paste0('Probes used in the MRS are saved in: ', outdir, cohort, '_probesinMRS.txt'))
} else {
  print(paste0('Number of probes read in for DNAm matches the number of weights provided: n = ', nrow(weights)))
}

missing_percentage <- DNAm %>% select(-!!sym(id_col)) %>%
  summarise_all(~ mean(is.na(.)) * 100) %>%
  gather(CpG, MissingPercentage) %>% arrange(desc(MissingPercentage))

print(paste0('The CpG with the highest level of missingness is: ', missing_percentage$CpG[1]))
print(paste0('There are ', nrow(missing_percentage %>% filter(MissingPercentage > 5)), ' CpGs with more than 5% missingness and ',
             nrow(missing_percentage %>% filter(MissingPercentage > 50)), ' with more than 50% missingness'))

print(missing_percentage)

missing_hist <- ggplot(missing_percentage, aes(x = MissingPercentage, fill = MissingPercentage < 50)) + 
  geom_histogram() + 
  theme_minimal() + 
  labs(x = '% of Missingness', y = 'Number of CpGs')+
  geom_vline(xintercept = 50, color = "red", linetype = 'dashed')+
  ggtitle(cohort)

missing_plot <- ggplot(missing_percentage %>% filter(MissingPercentage > 50),
                           aes(x = reorder(CpG, -MissingPercentage), y = MissingPercentage)) +
    geom_bar(stat='identity', fill = 'skyblue', color = 'black') + 
    labs(title = paste0(cohort, ': CpG missingness in MS'), x = "CpG", y = "% Missing") +
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

# Plot the % of missingness against the absolute value of the MRS coefficient 

missingness_weights <- merge(missing_percentage, weights, by = 'CpG')

missing_weights_plt <- ggplot(missingness_weights, aes(x = MissingPercentage, y = abs(Weight))) + 
  geom_point() + 
  theme_minimal() + 
  labs(x = '% of Missingness', y ='MS weight', 
       title = paste0(cohort, ': % missing vs MS weights for CpGs > 50% missingness'))

all_missing_plots <- ggarrange(missing_hist, missing_weights_plt, missing_plot, nrow = 3, ncol = 1)

ggsave(paste0(outdir, cohort, "_MRS_cpg_missingness.png"), missing_plot, width = 8, height = 6, device='png', dpi=300)
ggsave(paste0(outdir, cohort, "_MRS_cpg_missingness_hist.png"), missing_hist, width = 8, height = 6, device='png', dpi=300)
ggsave(paste0(outdir, cohort, "_MRS_cpg_missingness_weights.png"), missing_weights_plt, width = 8, height = 6, device='png', dpi=300)
ggsave(paste0(outdir, cohort, "_MRS_cpg_missingness_allplots.png"), all_missing_plots, width = 8, height = 6, device='png', dpi=300)

write.table(missing_percentage, paste0(outdir, cohort, '_cpg_missingness.txt'), row.names = F, quote = F)

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

# Distribution of the number of CpGs within the MRS (non-missing)

###############################################################################

# Showing as similar metric to the above missing plots
# But another way of looking at it 


num_cpgs <- DNAm_long %>%
  group_by(!!sym(id_col)) %>%
  summarise(cpgs_mrs = sum(!is.na(Mval)))

cpg_mrs_hist <- ggplot(num_cpgs, aes(x = cpgs_mrs)) + 
  geom_histogram() +
  theme_minimal() + 
  labs(x = 'Number of CpGs included in the MS', y = 'Frequency') + 
  ggtitle(paste0(cohort, ': Histogram of number of CpGs included in the MS'))

ggsave(paste0(outdir, cohort, "_indiv_numcpgs_MRS.png"), cpg_mrs_hist, width = 8, height = 6, device='png', dpi=300)

###############################################################################

# MRS Distributions

###############################################################################

print(paste0('Plotting distributions across the whole ', cohort, ' sample'))

# Distribution of the MRS across the DNAm cohort 

MRS_dist <- ggplot(MRS, aes(x = weighted_sum)) + 
  geom_histogram() + 
  theme_minimal() + 
  labs(x = 'Methylation Profile Score', y = 'Count')+
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
ad_pheno_MRS <- merge(ad_pheno, MRS, by = id_col)

print('Plotting MRS distributions for AD exposure cases and controls ')
MRS_pheno_dists <- ggplot(ad_pheno_MRS, aes(x = weighted_sum, fill = as.factor(antidep))) + 
  geom_histogram(alpha = 0.8) + 
  theme_minimal() + 
  labs(x = 'Methylation Profile Score', y = 'Count', fill = 'Self-reported AD use') +
  ggtitle(cohort)

ggsave(paste0(outdir, cohort, "_AD_MRS_phenodist.png"), MRS_pheno_dists, width = 8, height = 6, device='png', dpi=300)

###############################################################################

# Saving the methylation risk score 

###############################################################################
outfile <- paste0(outdir, cohort, '_AD_MRS.txt')
print(paste0('Saving the methylation risk score to ', outfile))

colnames(MRS)[2] <- 'AD_MRS'
write.table(MRS, outfile, row.names = F, quote = F)


