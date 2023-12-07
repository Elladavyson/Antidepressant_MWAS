###############################################################################

# Set up libraries and options/files

###############################################################################

.libPaths('/exports/igmm/eddie/GenScotDepression/users/edavyson/R/x86_64-pc-linux-gnu-library/4.1')
library(data.table)
library(dplyr)
library(optparse)
library(readr)
library(tools)
library(ggplot2)
library(tidyr)

parse <- OptionParser()

# setting up options for the filepaths to the correct files
option_list <- list(
  make_option('--cohort', type='character', help="Cohort, ideally no spaces (for graphs and documentation)", action='store'),
  make_option('--DNAm', type='character', help="The filepath for DNAm file", action='store'),
  make_option('--probes', type = 'character', help= "The filepath for the list of probes to be extracted (for MRS or CpG look-up)", action = 'store'),
  make_option('--id_column', type = 'character', default="IID", help = "Column names of identifier column", action = 'store'),
  make_option('--analysis', type = 'character', help = 'Name of analysis preprocessing is being performed for', action = 'store'),
  make_option('--outdir', type = 'character', help = 'The filepath for output directory', action = 'store')
)

args = commandArgs(trailingOnly=TRUE)
opt <- parse_args(OptionParser(option_list=option_list), args=args)

# setting up arguments from the options 
print('Setting up the options')
cohort <- opt$cohort
DNAm_filepath=opt$DNAm # DNAm file
probes_filepath=opt$probes # Probe list 
id_col <- opt$id_column # Vector of identifier columns 
analysis <- opt$analysis
out_dir <- opt$outdir

sink(paste0(out_dir, cohort, "_", analysis, "_DNAm_preproc.log"))

if (analysis == 'sig'){
  print('Preprocessing DNAm for the plotting the distributions of significant CpGs from the GS MWAS in an external cohort')
} else if (analysis == 'mrs') {
  print('Preprocessing DNAm for the calculation of a methylation risk score with weights from BIGLASSO in GS')
} else {
  stop('Please provide either sig or mrs to the analysis argument')
}
print(paste0('DNAm file from : ', DNAm_filepath))
print(paste0('List of probes from : ', probes_filepath))
print(paste0('ID column : ', id_col))
print(paste0('Output to be saved in : ', out_dir))


###############################################################################

# Read in the files

###############################################################################

if (endsWith(DNAm_filepath, ".rds")){
  DNAm <- readRDS(DNAm_filepath)
} else if (endsWith(DNAm_filepath, ".txt")){
  DNAm <- readr::read_table(DNAm_filepath) %>% as.data.frame()
} else {
  stop("Unsupported file format. Please provide a file with .rds or .txt extension")
}

probes <- readr::read_lines(probes_filepath)
print('Read in files')
print(paste0('Number of probes in DNAm file: ',DNAm %>% select(starts_with("cg")) %>% ncol()))

###############################################################################

# Filter the DNAm to just the CpGs used in the MRS (and identifiers)

###############################################################################
print('Filter to just the MRS CpGs')
DNAm_MRS <- DNAm %>% select(c(all_of(id_col), intersect(names(DNAm), probes)))

print(paste0('Filtered to ', DNAm_MRS %>% select(starts_with("cg")) %>% ncol(), ' CpGs')) 
rm(DNAm) # remove large methylation object 

###############################################################################

# Scale the DNAm levels 

###############################################################################
print('Scaling the CpG columns')

DNAm_MRS_std <- DNAm_MRS %>% mutate(across(-c(all_of(id_col)), scale))

# Plotting the distribution of unstandardised and standardised M vals
# for 3 randomly selected CpGs 

cpgs <- sample(setdiff(names(DNAm_MRS_std), id_col), 3)
DNAm_both <- rbind(
DNAm_MRS %>% 
  select(c(all_of(id_col), all_of(cpgs))) %>%
  pivot_longer(cols = -c(all_of(id_col)), 
               names_to = "CpG", 
               values_to = "Mval") %>% 
  as.data.frame() %>% 
  mutate(Values = 'Unstandardised'),
DNAm_MRS_std %>% 
  select(c(all_of(id_col), all_of(cpgs))) %>%
  pivot_longer(cols = -c(all_of(id_col)), 
               names_to = "CpG", 
               values_to = "Mval") %>%
  mutate(Values = 'Standardised')
)

DNAm_dists <- ggplot(DNAm_both, aes(x = Mval, fill = CpG)) +
  geom_histogram() + 
  facet_grid(Values~CpG) +
  ggtitle(paste0(cohort, ': Random sample of probes - standardisation'))

ggsave(filename=paste0(out_dir, cohort, "_", analysis, "_DNAm_preproc_std.png"),DNAm_dists, 
       width = 8, height = 6, device='png', dpi=300)


###############################################################################

# Write out the filtered and standardised DNAm data 

###############################################################################

outfile <- paste0(out_dir, cohort, "_", analysis, "_DNAm_preproc.txt")
write.table(DNAm_MRS_std, outfile, row.names = F, quote = F)
sink()
