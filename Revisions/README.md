# Antidepressant methylation revisions 
Code for the revisions of the antidepressant methylation paper

Insights from a Methylome-Wide Association Study of Antidepressant Exposure

## Demographics 

`demographic_statistics.R`: Getting demographics for the exposed and unexposed groups and running statistical tests on the differences between them
`mdd_group_psychiatric_assessment.R`: Getting more information on SCID diagnoses in the MDD group

## Sex stratified MWAS analysis (self-report)

`covar_qcovar_sex.R`: Generating the covariate files for sex stratified MWAS
`MWAS_sex`: Script running sex specific MOA MWAS models
`sex_OSCA_std`: Restandardising the methylation M-values for each sex tested
`sex_specific_effects.R`: Looking at the results of sex-stratified analyses, calculating Z-scores and P-values



