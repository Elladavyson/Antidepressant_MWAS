#!/bin/bash
#########################################################
#$ -N MOA_sex_segregated
#$ -l h_rt=48:00:00
#$ -l h_vmem=16G
#$ -cwd
#$ -pe sharedmem 8
#$ -e MOA_sex_segregated
#$ -o MOA_sex_segregated
#$ -M s2112198@ed.ac.uk
#$ -m baes

. /etc/profile.d/modules.sh
module load igmm/apps/osca/0.46

SCRATCH='/exports/eddie/scratch/s2112198'

pheno=$1
sex=$2

OUTPUTDIR="/exports/igmm/eddie/GenScotDepression/users/edavyson/antidep_project/revisions/output/MWAS_sex_segregated"
RESIDDIR="/exports/igmm/eddie/GenScotDepression/users/edavyson/antidep_project/antidep_phenotypes/phenotypes/05_23_phenos/resid_phenos"

pheno_filename="residualised_"${pheno}"_nocolnames.pheno"

echo "Phenotype file: ${pheno_filename}"

output_filename="${OUTPUTDIR}/"${pheno}"_GRMunadj_ORM_residph_"${sex}"_16_09"

echo "Output to be saved to file: ${output_filename}"

osca --moa --befile $SCRATCH/mvals_GRM_uncorrected_29_06_OSCA_standard --pheno $SCRATCH/residualised_${pheno}_${sex}_nocolnames.pheno --qcovar $SCRATCH/qcovs_11_07_GRMFID_${sex}_nocolnames --covar $SCRATCH/covs_11_07_GRMFID_${sex}_nocolnames --thread-num 8 --out ${output_filename}