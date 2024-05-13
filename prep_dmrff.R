
# Preparation of the methylation data file #
# for the DMR analysis using the dmrff package #
# which requires the methylation data to be in a standard format #
# and read in along with MWAS summary statistics #
# essentially removing polymirphic and cross hybridising probes and AHRR which are removed in the MWAS analysis
# and subsetting to the people used in the MWAS (removing 16 people)- by merging with the oii_file.

.libPaths('edavyson/R/x86_64-pc-linux-gnu-library/4.1')
# set wd

setwd('scratch/s2112198')

# libraries
library(data.table)
library(dplyr)
library(tidyverse)

# read in updated mvals info

mvals <- readRDS('s2112198/GRM_uncorrected/GS_18869.rds')

# read in the polymorphic and cross-hybridising probes

probes.rm = 'methylation/STRADL/Mvalues/SNP_CH_probes'
ls.probe.rm = read.delim(probes.rm,header=F,stringsAsFactors = F)

# remove polymorphic and cross-hybridising probes and the AHRR probe

mvals <- mvals[!(rownames(mvals) %in% ls.probe.rm$V1 | rownames(mvals) == 'cg05575921'), ]

# people in oii file (in the MWAS)

oii_file <- read.table('mvals_GRM_uncorrected_29_06_OSCA_standard.oii', header = F)
sample_info <- read_table('scratch/s2112198/sample_info.txt')

sample_oii <- merge(sample_info, oii_file, by.x = 'IID', by.y = 'V2')

# keep only the columns which are in the merged sample_oii file (SentrixID)

mvals <- mvals[, colnames(mvals) %in% sample_oii$Sentrix_ID]

# cpgs are columns and rows are people
# write out the new rds for the dmrff analysis

write_rds(mvals, 'mvals_dmrff_26_10.rds')
