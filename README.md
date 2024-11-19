# Antidepressant_MWAS
Repository for the manuscript "Antidepressant Exposure and DNA Methylation: Insights from a Methylome-Wide Association Study", currently available on medrxiv (https://www.medrxiv.org/content/10.1101/2024.05.01.24306640v1)

## Phenotypes and covariates 

`Prescription_data.Rmd`: R Markdown for preprocessing and parsing the prescription records, calculating antidepressant treatment episodes and deriving prescription derived phenotypes (All and MDD-only). 

`Self_report.Rmd`: R Markdown for extracting the self-report questionnaire data to derive the self-report antidepressant exposure phenotypes (All and MDD only). 

`OSCA_covariates.R`: R script prepping the covariate data for OSCA and the MWAS analysis

`GCTA_BLUP.sh`: Shell script which submits the phenotypes to GCTA to be residualised against the GRM, using the BLUP tool. 

## Methylome-wide association study

`MOA_MWAS.sh`: Shell script for running the MOA MWAS on the antidepressant exposure phenotype and the given covariates.

## Differentially methylated regions analysis 

`prep_dmrff.R`: R script prepping the data for dmrff analysis (the correlation matrix) using the dmrff R package

`dmrff.R`: R script running the dmrff analysis on the MWAS summary statistics.

## Methylation Profile Score 

`LASSO_pheno_prep.R`: R script which regresses the GRM-residualised phenotype against the MWAS covariates before running the LASSO model. The extracted residuals are then inputted into the LASSO model. 

`LASSO_GS_MPS.R`: R script which runs the LASSO model on antidepressant exposure and DNAm, to identify the CpGs for the methylation profile score. 

*Further scripts, alongside the information provided to external cohorts is provided at another repository: https://github.com/Elladavyson/AD_MRS/*

## Results processing 

`MWAS_downstream`: The visualisation and comparison of MWAS results, alongside results from downstream analysis (DMRFF, FUMA, SynGO, MSigDB, Time correlations). 
