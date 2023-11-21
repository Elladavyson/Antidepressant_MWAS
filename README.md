# AD_MRS
Repository for antidepressant MRS calculation in external cohorts. 

## Training the MRS in Generation Scotland 

### Regressing the self-report AD phenotype on  covariates 

The MRS was trained using a LASSO model to select informative (CpGs) features for antidepressant exposure. First the antidepressant phenotype was regressed on the GRM, using BLUP in GCTA. Then the residuals from this analysis were regressed on age, sex, lymphocyte cell proportions, monocyte cell proportions, AHRR methylation levels and (1| Batch). The residuals from this analysis were then taken forward into the LASSO model. 

### Running the LASSO model 

The LASSO model `antidep_resid ~ DNAm + covars` was run using the R package `biglasso`. All non-zero features were then extracted with their corresponding weights. These are then saved as the weights file to be used in the calculation of MRS in external cohorts. 

## Calculating MRS in an external cohort 

### Shared files needed 

`selfrep_lasso_probe_lst.txt` : Text file of the CpG 'cg' names used in calculating the MRS. 

### DNAm preprocessing 

The MRS consists of a weighted sum of 212 CpGs (trained in GS). The MRS was trained using *standardised DNAm levels* 
$$(X-mean)/standard deviation$$. Therefore we would like the MRS to be calculated also using standardised DNAm levels. 

**R** `process_DNAm_MRS.R` is a script which will read in the DNAm object (if saved as a .rds or .txt file) and will filter it to the CpGs in the MRS risk score and standardise (using `scale()`) the DNAm levels. Assumes that the DNAm object has rows as participants and columns as CpG names (alongside identifier column names). 

**OSCA**  `process_DNAm_OSCA.sh` is a template script to select and standardised DNAm object from BOD files (OSCA format). Change DNAm_filename and outfile_name accordingly. Assumes the DNAm object has rows as participants and columns as CpG names.


### Calculating MRS 

The R script `MRS_calc.R` allows you to calculate MRS for participants in the cohort given the following: 

Weights filepath: File path to the weights file provided from the training sample with Generation Scotland.
DNAm filepath: Filepath to the preprocessed DNAm file, each row is a participant and each column is a CpG site (excluding identifer columns)
Output filepath: Filepath for the output file to be saved (a table with ID identifer and MRS)

### Associational models 

Each cohort will have their own pipelines for associational models which can take into account relatedness/familial structure. Therefore, rather than a script there is detailed description of the predictor and covariates we would like included in the model. 
