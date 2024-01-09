
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

*--probes* : The file path for the list of probes used in the MRS e.g '/Users/data/DNAm/MRS_probes.txt'

*--id_column* : The name of the identifier column in the data e.g 'IID'

*--analysis*: Either 'sig' or 'mrs', to denote which analysis the file is being formatted for (prevent overwriting).

*--outdir* : The directory where the results and graphs will be saved e.g /Users/data/DNAm/AD_MRS/

**Example** : Rscript process_DNAm_MRS.R --cohort GS --DNAm test_212_probes_unstd.rds --probes MRS_probes.txt --id_column IID --analysis mrs --outdir /exports/eddie/scratch/s2112198/

**OSCA** : `process_DNAm_OSCA.sh` is a template script to select and standardised DNAm object from BOD files (OSCA format). Change DNAm_filename and outfile_name accordingly. Assumes the DNAm object has rows as participants and columns as CpG names.

### Output 

**R script** `process_DNAm_MRS.R`:

`{cohort}_mrs_DNAm_preproc.txt`: A text file of the standardised DNAm levels for the 212 CpGs required for making the MRS

`{cohort}_mrs_DNAm_preproc_std.png`: Distributions from three randomly selected CpGs to overview their distributions pre and post scaling (sanity check). 

`{cohort}_mrs_DNAm_preproc.log`: Log file

**OSCA** `process_DNAm_OSCA.sh`

`{cohort}_DNAm_preproc`: A text file of the standardised DNAm levels for the 212 CpGs required for making the MRS. 

## Calculating MRS 

The R script `MRS_calc.R`: will calculate MRS for participants. 

Arguments: 

*--cohort* : Cohort name, e.g 'GS' or 'LBC1936'

*--DNAm* : The file path for the preprocessed DNAm file from `process_DNAm_MRS.R`, e.g /Users/data/DNAm/AD_MRS/GS_mrs_DNAm_preproc.txt

*--id_column* : The column name of the identifier column (default == IID)

*--weights* : The file path for the MRS weights files provided by GS e.g /Users/data/DNAm/GS_AD_MRS_weights.txt

*--pheno* : The file path to the antidepressant exposure phenotype file for your cohort. Script can read in either a .csv, or .txt file and should follow a format of `FID`(Optional), `IID` (Required), and `antidep`, with `antidep` being coded as 0 (no exposure) and 1 (exposure).  e.g /Users/data/phenos/AD_pheno.csv

*--outdir* : The directory where the results and graphs will be saved e.g  /Users/data/DNAm/AD_MRS/

**Example** : Rscript MRS_calc.R --cohort GS --DNAm GS_mrs_DNAm_preproc.txt --weights GS_AD_MRS_weights.txt --id_column IID --pheno selfrep_pheno3_methyl_03_05.csv --outdir /exports/eddie/scratch/s2112198/

### OutPut 

#### Data 

`{cohort}_AD_MRS.txt`: a txt file of participants identifiers and their MRS ('AD_MRS'), needed for the associational model (`MRS_assoc.R`). 

`{cohort}_cpg_missingness.txt`: A txt file of CpGs and the % of missingness within the sample 

#### Plots 

`{cohort}_AD_MRS_overalldist.png`: A histogram of the MRS distribution in everyone in the cohort

`{cohort}_AD_MRS_phenodist.png`: A histogram of distributions for AD non-exposed and AD-exposed individuals in the cohort

`{cohort}_MRS_cpg_missingness.png`: A bar plot of all the CpGs missing in more than 50% of individuals and the % of missingness.

`{cohort}_MRS_cpg_missingness_weights.png`: A scatter plot of the % of missingness and the weight of the CpGs in the MRS.

`{cohort}_MRS_cpg_missingness_hist.png`: A histogram of the % of missingness in the CpGs included in the MRS.

`{cohort}_MRS_cpg_missingness_all_plots.png`: All the missingness plots, `{cohort}_MRS_cpg_missingness.png`, `{cohort}_MRS_cpg_missingness_weights.png` and `{cohort}_MRS_cpg_missingness_hist.png` plotted together using `ggarrange()` from `ggpubr`. 

`{cohort}_indiv_numcpgs_MRS.png`: A histogram of the number of CpGs included within the MRS for individuals (another metric of missingness). 

#### Log file 

`{cohort}_AD_MRS.log`: Log file

## Associational models 

Key elements of the model: 

**Phenotype** : Antidepressant exposure phenotype 

**Predictor** : Antidepressant exposure MRS

**Key Covariates**: Age, Sex, Methylation levels at AHRR(cg05575921), lymphocyte cell proportions, monocyte cell proportions, Technical covariates as a random effect (in GS, we often include (1|Batch)), *If related cohort* Top 10 genetic PCs/Kinship matrix 

*A note on accounting for relatedness*: Each cohort will have their own pipelines to account for relatedness within their cohorts, i.e Twin Cohorts. In this case, we welcome the use of in-house pipelines and models for testing the association of Phenotype ~ MRS, including the key covariates detailed above. If not using the example above, please do document your model and analysis. 

### General Output 

We would like the coefficients of the model, alongside the standard errors, Z scores, P values. 
We will assess fit using the AUC, ROC curves and a few pseudo-R2 measures (McFaddens pseudo-R2 and Nagelkerke's R2).
To calculate McFaddens and Nagelkerke's pseudo-R2, we need to model a *null* model, which is, the same model parameters as before but without the predictor of interest (**without the AD MRS**).

The McFaddens pseudo R2 is then calculated as a ratio of the loglikelihood of the full (including MRS) and null (excluding MRS) model: 

$$ McFaddens pseudo R^{2}=1-({LogLik(full model) \over LogLik(null model)})$$

Nagelkerkes pseudo R2 is calculated as 

$$ Nagelkerkes pseudo R^{2} = ({1-({Likelihood(full model)\over Likelihood(null model)}^({2\over N})) \over (1- Likelihood(null model)^({2 \over N}))}) $$

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

**Example** : Rscript MRS_assoc.R --cohort GS --id_column IID --mrs GS_AD_MRS.txt --pheno selfrep_pheno3_methyl_03_05.csv --covs GS_test_covs_pcs.txt --outdir /exports/eddie/scratch/s2112198/

### Output 

#### Data

`{cohort}_MRS_AD_coefficients.txt`: A table of the coefficients (Betas), Standard Errors, Z scores and P values from the association model results, extracted using `summary(assoc_mod)$coefficients %>% as.data.frame()`

`{cohort}_MRS_AD_modemetrics.txt`: A table of the log likelihood of the model including a DNAm predictor (`assoc_mod`), and when not including the predictor (`null_mod`), McFaddans R2 `mcf_r2`, Cox and Snells R2 (the Numerator of Nagelkerke's R2), Nagelkerkes R2, and the AUC of the model.

`{cohort}_roc_curve.rds`: A roc_curve object from the `pROC` package for plotting all cohorts together.

#### Plots

`{cohort}_assoc_ROC_curve.pdf`: A ROC curve for the cohort. 

#### Logs 

`{cohort}_MRS_AD_assoc.log`: Log file 

## Demographic information 

`cohort_demographics.R` will generate a table of demographic information for your sample (useful for the manuscript and interpretating our results). Currently it formats information on age, sex, bmi, smoking, pack years, AD MRS (generated from `MRS_calc.R` and mdd status (if applicable). 

*--cohort*: Cohort name, e.g 'GS' or 'LBC1936'

*--id_column*: The column name of the identifier column (default == IID)

*--mrs*: The filepath to the AD MRS file (made using MRS_calc.R,`{cohort}_AD_MRS.txt`, MRS column named `AD_MRS`)

*--pheno*: The filepath to the AD phenotype file

Column names: 

**Identifier column(s)** = 'IID' ('FID')

**Phenotype column** = 'antidep', coded as **0**-**No AD exposure** and **1**- **AD exposure**

*--demo*: Filepath to the file containing demographic information

Column names:

**Age**: 'age', numeric()

**Sex**: 'sex_coded', numeric(), coded as 0 (Female) and 1 (Male)

**BMI**: 'bmi', numeric()

**Smoking status ('Have you ever smoked tobacco?')**: 'ever_smoke', numeric(), coded as 1 ('Current smoker'), 2 ('Former- < 12 months), 3('Former, > 12 months), 4 (Never smoked), 5 (NA) 

**Pack years**: 'pack_years', numeric()

**Lifetime MDD status**: 'mdd', numeric(), coded as 0 (Controls) and 1 (Cases)

**Example** : Rscript cohort_demographics.R --cohort GS --id_column IID --mrs GS_AD_MRS.txt --pheno selfrep_pheno3_methyl_03_05.csv --demo GS_demograph.txt --outdir /exports/eddie/scratch/s2112198/

### Output 

`{cohort}_demo_summary.txt`: A summary table of the demographics for those exposed to ADs and those not exposed

`{cohort}_MRS_demographics.txt`: A log file for the Rscript. 
