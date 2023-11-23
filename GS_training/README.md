# Training the MRS in Generation Scotland 

## Regressing the self-report AD phenotype on  covariates 

The MRS was trained using a LASSO model to select informative (CpGs) features for antidepressant exposure. 

`resid_pheno_GCTA.txt`: First the antidepressant phenotype was regressed on the GRM, using BLUP in GCTA . 

`pheno_covs_resid_lmer.R`: Then the residuals from this analysis were regressed on age, sex, lymphocyte cell proportions, monocyte cell proportions, AHRR methylation levels and (1| Batch). 

The residuals from this analysis were then taken forward into the LASSO model. 

## Running the LASSO model 

`Lasso_regression_DNAm.R`: The LASSO model `antidep_resid ~ DNAm` was run using the R package `biglasso`. 

 LAll non-zero features were then extracted with their corresponding weights. These are then saved as the weights file to be used in the calculation of MRS in external cohorts. 
