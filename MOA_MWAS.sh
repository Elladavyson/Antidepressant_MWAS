#!/bin/bash
#########################################################
#$ -N MOA_uncorrected_resid_ORM_06_10_standardised
#$ -l h_rt=48:00:00
#$ -l h_vmem=64G
#$ -cwd
#$ -pe sharedmem 8
#$ -e resid_pheno_GRMuncorrected_MOA_ORM_standardised_06_10
#$ -o resid_pheno_GRMuncorrected_MOA_ORM_standardised_06_10
#$ -M s2112198@ed.ac.uk
#$ -m baes

. /etc/profile.d/modules.sh
module load igmm/apps/R/4.1.0
module load igmm/apps/osca/0.46

SCRATCH='/exports/eddie/scratch/s2112198'

pheno=$1

OUTPUTDIR="edavyson/antidep_project/MWAS_results/GRM_unadjusted_res_09_05/MOA/resid_pheno_ORM/06_10_update"
RESIDDIR="edavyson/antidep_project/antidep_phenotypes/phenotypes/05_23_phenos/resid_phenos"

pheno_filename="${pheno}_nocolnames.pheno"

echo "Phenotype file: ${pheno_filename}"

output_filename="${OUTPUTDIR}/GRM_unadjusted_${pheno}_MOA_ORM_residph_standard_06_10"

echo "Output to be saved to file: ${output_filename}"

osca --moa --befile $SCRATCH/mvals_GRM_uncorrected_29_06_OSCA_standard --pheno $RESIDDIR/residualised_${pheno}_nocolnames.pheno --qcovar $SCRATCH/qcovs_11_07_GRMFID_nocolnames --covar $SCRATCH/covs_11_07_GRMFID_nocolnames --thread-num 8 --out ${output_filename}
