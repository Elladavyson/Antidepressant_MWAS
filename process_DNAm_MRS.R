###############################################################################

# Set up libraries and options/files

###############################################################################


library(data.table)
library(dplyr)
library(optparse)
library(readr)

parse <- OptionParser()

# setting up options for the filepaths to the correct files
option_list <- list(
  make_option('--DNAm', type='character', help="The filepath for DNAm file", action='store'),
  make_option('--probes', type = 'character', help= "The filepath for the list of probes for the MRS", action = 'store'),
  make_option('--id_columns', type = 'character', nargs=2, default=c("FID", "IID"), help = "Column names of identifier column(s)", action = 'store'),
  make_option('--out', type = 'character', help = 'The filepath for output', action = 'store')
)

args = commandArgs(trailingOnly=TRUE)
opt <- parse_args(OptionParser(option_list=option_list), args=args)

# setting up arguments from the options 

DNAm_filepath=opt$DNAm # DNAm file
probes_filepath=opt$probes # Probe list 
id_cols <- opt$id_columns # Vector of identifier columns 
outfile <- opt$out

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

###############################################################################

# Filter the DNAm to just the CpGs used in the MRS (and identifiers)

###############################################################################

DNAm_MRS <- DNAm %>% select(c(id_cols, intersect(names(DNAm), probes)))

###############################################################################

# Scale the DNAm levels 

###############################################################################

DNAm_MRS_std <- DNAm_MRS %>% mutate(across(-c(FID, IID), scale))

###############################################################################

# Write out the filtered and standardised DNAm data 

###############################################################################

write.table(DNA_MRS_std, outfile, row.names = F, quote = F)
