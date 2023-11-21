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
  make_option('--DNAm', type='character', help="The filepath for DNAm file", action='store'),
  make_option('--probes', type = 'character', help= "The filepath for the list of probes for the MRS", action = 'store'),
  make_option('--id_column', type = 'character', default="IID", help = "Column names of identifier column", action = 'store'),
  make_option('--out', type = 'character', help = 'The filepath for output', action = 'store')
)

args = commandArgs(trailingOnly=TRUE)
opt <- parse_args(OptionParser(option_list=option_list), args=args)

# setting up arguments from the options 
print('Setting up the options')
DNAm_filepath=opt$DNAm # DNAm file
probes_filepath=opt$probes # Probe list 
id_col <- opt$id_column # Vector of identifier columns 
outfile <- opt$out

sink(paste0(file_path_sans_ext(outfile), ".log"))
print(paste0('DNAm file from : ', DNAm_filepath))
print(paste0('List of probes from : ', probes_filepath))
print(paste0('ID column : ', id_col))
print(paste0('Output to be saved : ', outfile))

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
print(paste0('Number of probes in DNAm file: ',ncol(DNAm)))

###############################################################################

# Filter the DNAm to just the CpGs used in the MRS (and identifiers)

###############################################################################
print('Filter to just the MRS CpGs')
DNAm_MRS <- DNAm %>% select(c(all_of(id_col), intersect(names(DNAm), probes)))

print(paste0('Filtered to ', ncol(DNAm_MRS)-1, ' CpGs')) # not counting the ID column
rm(DNAm) # remove large methylation object 

###############################################################################

# Scale the DNAm levels 

###############################################################################
print('Scaling the CpG columns')

DNAm_MRS_std <- DNAm_MRS %>% mutate(across(-c(all_of(id_col)), scale))

# pick a random CpG to plot the distribution?

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

DNAm_dists <- ggplot(DNAm_both, aes(x = Mval, fill = CpG)) + geom_histogram() + facet_grid(Values~CpG) + ggtitle('Random sample of probes - standardisation')
ggsave(filename=paste0(file_path_sans_ext(outfile), ".png"),DNAm_dists, width = 8, height = 6, device='png', dpi=300)


###############################################################################

# Write out the filtered and standardised DNAm data 

###############################################################################
print(paste0('Writing the processed DNAm file to ', outfile))

write.table(DNAm_MRS_std, outfile, row.names = F, quote = F)
sink()
