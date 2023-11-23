
# Calculating MRS in an external cohorts

## Shared files needed 

`GS_AD_MRS_weights.txt`: Text file of the CpGs and their weights 

`MRS_probes.txt` : Text file of the CpG 'cg' names used in calculating the MRS. 

## DNAm preprocessing 

The MRS consists of a weighted sum of 212 CpGs (trained in GS). The MRS was trained using *standardised DNAm levels* 
$$X_{std}={(X-mean)\over Std}$$

Therefore we would like the MRS to be calculated also using standardised DNAm levels. 

**R**: 

`process_DNAm_MRS.R`: which will read in the DNAm object (if saved as a .rds or .txt file) and will filter it to the CpGs in the MRS risk score and standardise (using `scale()`) the DNAm levels. Assumes that the DNAm object has rows as participants and columns as CpG names (alongside identifier column names). 

Arguments: 

*--cohort* : Cohort name, e.g 'GS' or 'LBC1936'

*--DNAm* : DNAm file path e.g '/Users/data/DNAm/DNAm.rds'

*--probes* : The file path for the list of probes used in the MRS e.g '/Users/data/DNAm/selfrep_lasso_probe_lst.txt'

*--id_column* : The name of the identifier column in the data e.g 'IID'

*--outdir* : The directory where the results and graphs will be saved e.g /Users/data/DNAm/AD_MRS/

**Example** : Rscript process_DNAm_MRS.R --cohort GS --DNAm test_212_probes_unstd.rds --probes MRS_probes.txt --id_column IID --outdir /exports/eddie/scratch/s2112198/

**OSCA** : `process_DNAm_OSCA.sh` is a template script to select and standardised DNAm object from BOD files (OSCA format). Change DNAm_filename and outfile_name accordingly. Assumes the DNAm object has rows as participants and columns as CpG names.

### Output 

**R script** `process_DNAm_MRS.R`:

`{cohort}_DNAm_preproc.txt`: A text file of the standardised DNAm levels for the 212 CpGs required for making the MRS

`{cohort}_DNAm_preproc_std.png`: Distributions from three randomly selected CpGs to overview their distributions pre and post scaling (sanity check). 

`{cohort}_DNAm_preproc.log`: Log file

**OSCA** `process_DNAm_OSCA.sh`

`{cohort}_DNAm_preproc`: A text file of the standardised DNAm levels for the 212 CpGs required for making the MRS. 

## Calculating MRS 

The R script `MRS_calc.R`: will calculate MRS for participants. 

Arguments: 

*--cohort* : Cohort name, e.g 'GS' or 'LBC1936'

*--DNAm* : The file path for the preprocessed DNAm file from `process_DNAm_MRS.R`, e.g /Users/data/DNAm/AD_MRS/GS_DNAm_preproc.txt

*--id_column* : The column name of the identifier column (default == IID)

*--weights* : The file path for the MRS weights files provided by GS e.g /Users/data/DNAm/big_lasso_450K_selfrep.txt

*--pheno* : The file path to the antidepressant exposure phenotype file for your cohort. Script can read in either a .csv, or .txt file and should follow a format of `FID`(Optional), `IID` (Required), and `antidep`, with `antidep` being coded as 0 (no exposure) and 1 (exposure).  e.g /Users/data/phenos/AD_pheno.csv

*--outdir* : The directory where the results and graphs will be saved e.g  /Users/data/DNAm/AD_MRS/

**Example** : Rscript MRS_calc.R --cohort GS --DNAm GS_DNAm_preproc.txt --weights GS_AD_MRS_weights.txt --id_column IID --pheno selfrep_pheno3_methyl_03_05.csv --outdir /exports/eddie/scratch/s2112198/

### OutPut 

`{cohort}_AD_MRS.txt`: a txt file of participants identifiers and their MRS ('AD_MRS')

`{cohort}_AD_MRS_overalldist.png`: A histogram of the MRS distribution in everyone in the cohort

`{cohort}_AD_MRS_phenodist.png`: A histogram of distributions for AD non-exposed and AD-exposed individuals in the cohort

`{cohort}_AD_MRS.log`: Log file

## Associational models 

Key elements of the model: 

**Phenotype** : Antidepressant exposure phenotype 

**Predictor** : Antidepressant exposure MRS

**Key Covariates**: Age, Sex, Methylation levels at AHRR(cg05575921), lymphocyte cell proportions, monocyte cell proportions, Technical covariates as a random effect (in GS, we often include (1|Batch)), *If related cohort* Top 10 genetic PCs/Kinship matrix 

*A note on accounting for relatedness*: Each cohort will have their own pipelines to account for relatedness within their cohorts, i.e Twin Cohorts. In this case, we welcome the use of in-house pipelines and models for testing the association of Phenotype ~ MRS, including the key covariates detailed above. If not using the example above, please do document your model and analysis. 

### General Output 

We would like the coefficients of the model, alongside the standard errors, Z scores, P values and the McFaddens pseudo-R2 to assess the wellness of fit. To calculate McFaddens pseudo-R2, we need to model a *null* model, which is, the same model parameters as before but without the predictor of interest (**without the AD MRS**).

The McFaddens pseudo R2 is then calculated as a ratio of the loglikelihood of the full (including MRS) and null (excluding MRS) model: 

$$pseudo R^{2}=1-({LogLik(full model)\over LogLik(nullmodel)})$$

### Example Script 

`MRS_assoc.R` is a general script which tests the association of the MRS  with antidepressant exposure phenotype using a generalised mixed linear model (glmm).

*--cohort*: Cohort name, e.g 'GS' or 'LBC1936'

*--id_column*: The column name of the identifier column (default == IID)

*--mrs*: The filepath to the AD MRS file (made using MRS_calc.R,`{cohort}_AD_MRS.txt`, MRS column named `AD_MRS`)

Column names: 

**Identifier column(s)** = 'IID' ('FID')

**MRS column** = 'AD_MRS'

*--pheno*: The filepath to the AD phenotype file

Column names: 

**Identifier column(s)** = 'IID' ('FID')
**Phenotype column** = 'antidep', coded as **0**-**No AD exposure** and **1**- **AD exposure**

*--covs*: The filepath to the covariate file (format: Identifier columns and covariates to include in the model). 

Column names: 

**Age** = 'age'

**Sex** = 'sex_coded' (0/1)

**AHRR**= 'cg05575921'

**Monocyte cell proportions** = 'Mono'

**Lymphocyte cell proportions** = 'lymphocytes'

**Batch = 'Batch'**

**PCS = 'C1-10'**

*--outdir* : The directory where the results and graphs will be saved e.g  /Users/data/DNAm/AD_MRS/

**Example** : Rscript MRS_assoc.R --cohort GS --id_column IID --mrs GS_AD_MRS.txt --pheno selfrep_pheno3_methyl_03_05.csv --covs GS_test_covs_pcs.txt --id_column IID --outdir /exports/eddie/scratch/s2112198/

### Output 

`{cohort}_MRS_AD_coefficients.txt`: A table of the coefficients (Betas), Standard Errors, Z scores and P values from the association model results, extracted using `summary(assoc_mod)$coefficients %>% as.data.frame()`

`{cohort}_MRS_AD_logL.txt`: A table of the log likelihood of the model including a DNAm predictor (`assoc_mod`), and when not including the predictor (`null_mod`), and the McFaddans R2 `mcf_r2` calculated as `1-logLik(assoc_mod)/logLik(null_mod)`. 

`{cohort}_MRS_AD_assoc.log`: Log file 
