# AD_MRS
Repository for antidepressant MRS calculation in external cohorts. 


## GS training (training DNAm score in GS) *Scripts for reference* 

`resid_pheno_GCTA.txt`: Residualising AD exposure phenotype on the GRM using GCTA (BLUP).

`pheno_covs_resid_lmer.R`: Regressing the AD-GRM residuals on covariates including in MWAS (age, sex, WBC proportions, Batch and AHRR methylation).

`Lasso_regression_DNAm.R`: Running big lasso model of `antidepressant exposure (residuals) ~ DNAm` to get CpGs and corresponding weights for the AD MRS. 

## MRS analysis (Applying DNAm scores to external cohorts) *Scripts to run*

`process_DNAm_MRS.R` : Rscript for filtering and standardising DNAm data  DNAm data (in .rds/.txt format).

`process_DNAm_OSCA.sh`: Shell script for filtering and standardising DNAm data (in BOD format).

`MRS_calc.R`: Rscript to calculate antidepressant exposure MRS using GS weights.

`MRS_assoc.R`: Rscript to run GLMM of antidepressant exposure ~ MRS + covariates (age, sex, WBC proportions, Batch, AHRR methylation and genetic PCs).

*NB* The `MRS_assoc.R` script is a template for modelling the association between the phenotype and MRS, however we anticipate the use of other models to account for population and family structure of individual cohorts. In this case, please follow usual pipelines and return the coefficients and statistics of the model, McFadden's pseudo-R2 (more detail in the README.md in this folder), and documentation of the model used.

### Contact

If you have any questions, please contact me at s2112198@ed.ac.uk OR e.e.davyson@sms.ed.ac.uk
