# Plotting the CpG distributions of significant CpGs from GS

## Files required 

List of probes significant in the GS MWAS: MWAS_signif_probes.txt
 
## Processing the DNAm 

`process_DNAm_MRS.R`: Can be used to extract and standardise the significant CpGs from a DNAm object. The exact same as for the MRS preprocessing, but instead extracting the significant MWAS CpGs rather than the probes for the MRS. 

Arguments: 

*--cohort* : Cohort name, e.g 'GS' or 'LBC1936'

*--DNAm* : DNAm file path e.g '/Users/data/DNAm/DNAm.rds'

*--probes* : The file path for the list of significant probes in the MWAS e.g '/Users/data/DNAm/MWAS_signif_probes.txt'

*--id_column* : The name of the identifier column in the data e.g 'IID'

*--analysis*: Either 'sig' or 'mrs', to denote which analysis the file is being formatted for (prevent overwriting).

*--outdir* : The directory where the results and graphs will be saved e.g /Users/data/DNAm/AD_MRS/

**Example** Rscript process_DNAm_MRS.R --cohort GS --DNAm test_212_probes_unstd.rds --probes MWAS_signif_probes.txt --id_column IID --analysis sig --outdir /exports/eddie/scratch/s2112198/

### Output 

`{cohort}_sig_DNAm_preproc.txt`: A text file of the standardised DNAm levels for the 212 CpGs required for making the MRS

`{cohort}_sig_DNAm_preproc_std.png`: Distributions from three randomly selected CpGs to overview their distributions pre and post scaling (sanity check). 

`{cohort}_sig_DNAm_preproc.log`: Log file



## Plotting the distributions 

`signif_cpg_ttest.R`: Will run t-tests on the CpG Mval distribution in cases and controls, and also plot the distributions of the probes. 

Arguments: 

*--cohort* : Cohort name, e.g 'GS' or 'LBC1936'

*--id_column* : 

*--probes_dat* : The file path for the standardised data table of significant CpGs e.g '/Users/data/DNAm/GS_DNAm_sig_preproc.txt'

*--pheno*: The file path for the antidepressant exposure phenotype (.txt or .csv)

*--outdir* : The directory where the results and graphs will be saved e.g /Users/data/DNAm/AD_MRS/

### Output 

**Graphs**

`{cohort}_all_dists.png`: Plots of the distributions of the probes in AD cases and controls

`{cohort}_violin_{cpg}.png`: Violin plots of a distribution, annotated with t-test results (one plot produced for each CpG)

`{cohort}_all_violins.png`: All Violin plots plotted together

**Tables**

`{cohort}_probe_ttest.tsv`: Results from the t-test in tabular format 

`{cohort}_signif_cpg_lookup.log`: Log file for script 

