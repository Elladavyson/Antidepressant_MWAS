# AD_MRS
Repository for antidepressant MRS calculation in external cohorts. 

## Training the MRS in Generation Scotland 

### Regressing the self-report AD phenotype on  covariates 

The MRS was trained using a LASSO model to select informative (CpGs) features for antidepressant exposure. First the antidepressant phenotype was regressed on the GRM, using BLUP in GCTA. Then the residuals from this analysis were regressed on age, sex, lymphocyte cell proportions, monocyte cell proportions, AHRR methylation levels and (1| Batch). The residuals from this analysis were then taken forward into the LASSO model. 

### Running the LASSO model 

The LASSO model `antidep_resid ~ DNAm + covars` was run using the R package `biglasso`. All non-zero features were then extracted with their corresponding weights. These are then saved as the weights file to be used in the calculation of MRS in external cohorts. 

## Calculating MRS in an external cohort 

### DNAm preprocessing 

*Optional*: As the project works with cohorts with EPIC and/or 450K array, we have restricted the MRS to just the overlapping CpGs on the 450K and EPIC arrays (n = 365, 912). A list of the CpGs are saved in the file "", and can be used to filter down the DNAm object (to reduce processing time) prior to the MRS calculation. 

Our MRS was trained on standardised DNAm levels: (X-mean)/standard deviation. Therefore we would ideally like the MRS to be calculated in external cohorts also on standardised methylation levels. 

Depending on how the DNAm data is stored, it can be standardised in different ways (example scripts for each).  

**R**: Using the function scale()
**OSCA**: Using the flag --std-probe

### Calculating MRS 

The R script `MRS_calc.R` allows you to calculate MRS for participants in the cohort given the following: 

Weights filepath: File path to the weights file provided from the training sample with Generation Scotland.
DNAm filepath: Filepath to the preprocessed DNAm file, each row is a participant and each column is a CpG site (excluding identifer columns)
Output filepath: Filepath for the output file to be saved (a table with ID identifer and MRS)

### Associational models 

Each cohort will have their own pipelines for associational models which can take into account relatedness/familial structure. Therefore, rather than a script there is detailed description of the predictor and covariates we would like included in the model. 
