# AD_MRS
Repository for antidepressant MRS calculation in external cohorts. 

## Training the MRS in Generation Scotland 

### Regressing the self-report AD phenotype on  covariates 

The MRS was trained using a LASSO model to select informative (CpGs) features for antidepressant exposure. 

`resid_pheno_GCTA.txt`: First the antidepressant phenotype was regressed on the GRM, using BLUP in GCTA . 

`pheno_covs_resid_lmer.R`: Then the residuals from this analysis were regressed on age, sex, lymphocyte cell proportions, monocyte cell proportions, AHRR methylation levels and (1| Batch). 

The residuals from this analysis were then taken forward into the LASSO model. 

### Running the LASSO model 

`Lasso_regression_DNAm`: The LASSO model `antidep_resid ~ DNAm + covars` was run using the R package `biglasso`. All non-zero features were then extracted with their corresponding weights. These are then saved as the weights file to be used in the calculation of MRS in external cohorts. 

## Calculating MRS in an external cohort 

### Shared files needed 

`selfrep_lasso_probe_lst.txt` : Text file of the CpG 'cg' names used in calculating the MRS. 
`big_lasso_450K_selfrep.txt`: Text file oif the CpGs and their weights 

### DNAm preprocessing 

The MRS consists of a weighted sum of 212 CpGs (trained in GS). The MRS was trained using *standardised DNAm levels* 
$$\(X-mean)\over Std $$

Therefore we would like the MRS to be calculated also using standardised DNAm levels. 

**R**: `process_DNAm_MRS.R` is a script which will read in the DNAm object (if saved as a .rds or .txt file) and will filter it to the CpGs in the MRS risk score and standardise (using `scale()`) the DNAm levels. Assumes that the DNAm object has rows as participants and columns as CpG names (alongside identifier column names). 

Arguments: 

*--cohort* : Cohort name, e.g 'GS' or 'LBC1936'

*--DNAm* : DNAm file path e.g '/Users/data/DNAm/DNAm.rds'

*--probes* : The file path for the list of probes used in the MRS e.g '/Users/data/DNAm/selfrep_lasso_probe_lst.txt'

*--id_column* : The name of the identifier column in the data e.g 'IID'

*--out_dir* : The directory where the results and graphs will be saved e.g /Users/data/DNAm/AD_MRS/

This script will produce a text file of the 212 standardised CpGs required for making the MRS, in a file `{cohort}_DNAm_preproc.txt`. As a sanity check, it also produces a graph of distributions from three randomly selected CpGs to overview their distributions pre and post scaling `{cohort}_DNAm_preproc_std.png`. 

**OSCA** : `process_DNAm_OSCA.sh` is a template script to select and standardised DNAm object from BOD files (OSCA format). Change DNAm_filename and outfile_name accordingly. Assumes the DNAm object has rows as participants and columns as CpG names.


### Calculating MRS 

The R script `MRS_calc.R` allows you to calculate MRS for participants. 

Arguments: 

*--cohort* : Cohort name, e.g 'GS' or 'LBC1936'

*--DNAm* : The file path for the preprocessed DNAm file from `process_DNAm_MRS.R`, e.g /Users/data/DNAm/AD_MRS/GS_DNAm_preproc.txt

*--id_column* : The column name of the identifier column

*--weights* : The file path for the MRS weights files provided by GS e.g /Users/data/DNAm/big_lasso_450K_selfrep.txt

*--pheno* : The file path to the antidepressant exposure phenotype file for your cohort. Script can read in either a .csv, or .txt file and should follow a format of `FID`(Optional), `IID` (Required), and `antidep`, with `antidep` being coded as 0 (no exposure) and 1 (exposure).  e.g /Users/data/phenos/AD_pheno.csv

*--outdir* : The directory where the results and graphs will be saved e.g  /Users/data/DNAm/AD_MRS/

This script will save a txt file of participants identifiers and their MRS, `{cohort}_AD_MRS.txt`. It should also save a histogram of the MRS distribution in everyone `{cohort}_AD_MRS_overalldist.png`, and another which has distributions for AD non-exposed and AD-exposed individuals in your cohort `{cohort}_AD_MRS_phenodist.png`. 

### Associational models 

Each cohort will have their own pipelines for associational models which can take into account relatedness/familial structure which may be unique to each cohort.

*Phenotype* : Antidepressant exposure phenotype 

*Predictor* : Antidepressant exposure MRS

*Covariates*: Age, Sex, Methylation levels at AHRR(cg05575921), lymphocyte cell proportions, monocyte cell proportions, Technical covariates included in your cohort (in GS, we often include (1|Batch)), Top 10 genetic PCs*

*Accounting for relatedness: Each cohort will have their own pipelines for best taking into account familial structures within their cohort, i.e for Twin Cohorts we anticipate perhaps the use of generalised estimation equations and `GEE` rather than a mixed linear model. Please do follow your own pipelines/versions for testing the association between the AD-exposure MRS and AD exposure, please include the same covariates and document the model/processes used.




If you have any questions, please contact me at s2112198@ed.ac.uk OR e.e.davyson@sms.ed.ac.uk
