module load igmm/apps/osca/0.46 # or equivilent 

cohort="" #Input the cohort name , e.g "GS"
osca --befile DNAm_filename --extract-probe selfrep_lasso_probe_lst.txt --std-probe --make-efile --out ${cohort}_DNAm_preproc

