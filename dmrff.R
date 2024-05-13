.libPaths('edavyson/R/x86_64-pc-linux-gnu-library/4.1')
library(data.table)
library(tidyverse)
library(dmrff)
library(readr)
library(dplyr)
library(optparse)
library(officer)

parse <- OptionParser()

# setting up options to sopecify which MWAS results files we are running the DMR on
option_list <- list(
  make_option('--phenotype', type='character', help="The phenotype (for selecting the right MWAS summary statistics", action='store')
)

args = commandArgs(trailingOnly=TRUE)

# extracting the phenotype info (i.e "antidep_pheno1_clean_appt" or "selfrep_pheno3")

opt <- parse_args(OptionParser(option_list=option_list), args=args)
pheno_name=opt$phenotype # phenotype file 

# setting the outpath for the output files

outpath <- 'edavyson/antidep_project/MWAS_results/GRM_unadjusted_res_09_05/MOA/resid_pheno_ORM/06_10_update/DMR/'
f.out=paste0(outpath, pheno_name,'_dmr_res.tsv') #output location 

# setting up a log file for the script

logging <- function(x,x.f=gsub('.tsv','.log',f.out),init=F, timing = F){ #function for pasting things into a log file 
  if (timing) {
    x <- paste0(x, " - Time taken: ", Sys.time() - start_time)
  }
  
  x <- paste0(x, collapse = ',')
  x=x %>% paste0(.,collapse = ',')
  if(init==T){ # i.e first thing being logged 
    system(paste0('echo ',x,'> ',x.f)) # write over the file
  }else{ # not the first thing being logged 
    system(paste0('echo ',x,'>> ',x.f)) # append to the file
  }
}

logging('Setting the WD', init = T)

setwd('scratch/s2112198')

logging(paste0('WD: ', getwd()))

# read in the MWAS summary stats 

logging('Reading in the stats file')

stats <- read_table(paste0('edavyson/antidep_project/MWAS_results/GRM_unadjusted_res_09_05/MOA/resid_pheno_ORM/06_10_update/GRM_unadjusted_', pheno_name, '_MOA_ORM_residph_standard_06_10.moa')) %>%
  as.data.frame()

logging(head(stats))

logging(paste0('Dimensions of MWAS summary statistics', dim(stats)))

#Below we assume that you have already performed an epigenome-wide association analysis
#and have loaded your summary statistics in R as a data frame `stats`
#which has the following columns:
# - `estimate` (regression coefficient),
#- `se` (standard error of the coefficient),
#- `p.value`,
#- `chr` (chromosome of the CpG site),
#- `pos` (position of the CpG site on the chromosome).

stats <- stats %>% rename(
  estimate = b,
 p.value=p,
  chr=Chr,
  pos=bp
)

# reading in the methylation data (prepped for dmrff with dmrff_prep.R)
logging('Reading in the methylation file')
start_time <- Sys.time()
methylation <- readRDS('scratch/s2112198/mvals_dmrff_26_10.rds')
logging('Time reading in the methylation file', timing = TRUE)

# check that the rows in stats matches the rows in methylation

if (nrow(stats) != nrow(methylation)){
  logging('Number of CpGs in the summary statistics != CpGs in the methylation matrix')
  logging(paste0('The mismatched CpGs: ', stats$Probe[!(stats$Probe %in% row.names(methylation))]))
  logging('Filtering the CpGs to match those in the methylation matrix')
  stats <- stats %>% filter(Probe %in% row.names(methylation))
} else {
  logging('Number of CpGs in the summary statistics == CpGs in the methylation matrix')
}

logging('Performing the DMR tests using the dmrff package')

output <- capture.output({
  dmrs <- dmrff(
    estimate = stats$estimate,
    se = stats$se,
    p.value = stats$p.value,
    methylation = methylation,
    chr = stats$chr,
    pos = stats$pos,
    maxgap = 500,
    verbose = TRUE
  )
})

logging(output)
logging(head(dmrs))

# write out the dmr results 

write.table(dmrs, f.out, sep = '\t', row.names = F, quote = F)
