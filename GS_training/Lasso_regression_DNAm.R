### GS LASSO Regression ###

.libPaths('/exports/igmm/eddie/GenScotDepression/users/edavyson/R/x86_64-pc-linux-gnu-library/4.1')

# libraries

library(data.table)
library(dplyr)
library(readr)
library(biglasso)

# establish the start time of the script 

start_time <- Sys.time()

#--------------------------------------------------------------------------------------------

# Logging function 

#--------------------------------------------------------------------------------------------


log_path='/exports/eddie/scratch/s2112198/GS_train_DNAm_antidep_450K.log'

logging <- function(x,x.f = log_path,init=F, timing = F){ #function for pasting things into a log file 
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


#--------------------------------------------------------------------------------------------

# Reading in the GS DNAm Data (subsetting to overlapping 450K probes)

#--------------------------------------------------------------------------------------------

# saved as an efile
# osca --befile mvals_GRM_uncorrected_29_06_OSCA_standard --extract-probe 450_EPIC_overlap.txt --make-efile --out selfrep_stand_450K_overlap

#DNAm <- read_table('/exports/eddie/scratch/s2112198/test_probe_standard') %>% as.data.frame() # subset of probes for trouble shooting 
DNAm <- read_table('/exports/eddie/scratch/s2112198/selfrep_stand_450K_overlap') %>% as.data.frame() # the whole data object
logging('DNAm data read in', timing = T)

rownames(DNAm) = DNAm$IID
DNAm <- DNAm %>% select(-c(FID, IID))

#--------------------------------------------------------------------------------------------

# Reading in the residualised self-report phenotype 
# Residuals from a linear mixed model of self-report (GRM-residualised) phenotype regressed on MWAS covariates
# resid_antidep ~ scale(age)+scale(Mono)+scale(lymphocytes)+ scale(cg05575921)+ as.factor(sex_coded) + (1|Batch)

#--------------------------------------------------------------------------------------------

# from datastore location 
# /exports/igmm/datastore/GenScotDepression/users/edavyson/Antidep_methylation/MRS/selfrep_cov_residuals.txt'


logging('Reading in phenotype files')
resid_pheno <- read.table('/exports/eddie/scratch/s2112198/selfrep_cov_residuals.txt', header = T)

# make row names the IID's
rownames(resid_pheno) <- resid_pheno$IID


#--------------------------------------------------------------------------------------------

# Selecting the intersect between DNAm and phenotype data 

#--------------------------------------------------------------------------------------------

# Identify the people in the DNAm file which map those in the phenotype file 
# Should be all of them in the phenotype file 

a = intersect(row.names(DNAm), row.names(resid_pheno))
print(paste0('N=',length(a)))

# filter out the rows in the DNAm for the people in the phenotype file 

X <- DNAm[a,]

# filter out those in the phenotype file with DNAm data (should be all)

y <- resid_pheno[a,]

# Check ordering is correct
if (!all.equal(row.names(X), row.names(y))){
  stop('Methylation and phenotype data need to have matching IDs',call.=F)
}

# remove the large methylation object for memory 
rm(DNAm)

#--------------------------------------------------------------------------------------------

# Fitting Big LASSO

#--------------------------------------------------------------------------------------------

tmp.y=y[,'pheno_residuals']
X.bm <- as.big.matrix(X[!is.na(tmp.y),])
print('Big.matrix transformed')
y.input = tmp.y[!is.na(tmp.y)]
cvfit <- cv.biglasso(X.bm, y.input, seed = 1234, nfolds = 10, ncores = 1) 
print('Lasso finished')
lambda <- cvfit$lambda.min

#--------------------------------------------------------------------------------------------

# Saving Big LASSO results 

#--------------------------------------------------------------------------------------------

logging(cvfit)
logging(paste0('Lambda selected: ', lambda))

#writing a dataframe of all the non-zero coefficients 

out <- as.data.frame(cbind(rownames(coef(cvfit))[which(coef(cvfit) != 0)], coef(cvfit)[which(coef(cvfit) != 0)]))

write.table(out, file='/exports/eddie/scratch/s2112198/big_lasso_450K_selfrep.txt', quote=F, row.names=F, col.names=F)
print('Lasso results saved')

logging('Script done!', timing = TRUE)


#--------------------------------------------------------------------------------------------

# Post script processing (on EDDIE, documenting )

#--------------------------------------------------------------------------------------------

# Extracting a list of the probes 

#results <- read.table('big_lasso_450K_selfrep.txt', header = F)
#colnames(results) <- c('CpG', 'Weight')
#results <- results[-1,] # remove the intercept row (don't need?)
#write.table(results, 'GS_AD_MRS_weights.txt', row.names= F, quote = F)
#readr::write_lines(results$Weight, 'MRS_probes.txt')





