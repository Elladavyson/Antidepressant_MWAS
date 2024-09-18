#!/bin/bash
#########################################################
#$ -N MOA_extra_probe
#$ -l h_rt=48:00:00
#$ -l h_vmem=16G
#$ -l rl9=false
#$ -cwd
#$ -e /exports/igmm/eddie/GenScotDepression/users/edavyson/antidep_project/revisions/output/
#$ -o /exports/igmm/eddie/GenScotDepression/users/edavyson/antidep_project/revisions/output/
#$ -M s2112198@ed.ac.uk
#$ -m baes

# First need to extract the probe from the OSCA object
#### nano extra_probe.txt
# cg08527546
# cp extra_probe.txt /exports/eddie/scratch/s2112198/
. /etc/profile.d/modules.sh
module load igmm/apps/osca/0.46
SCRATCH="/exports/eddie/scratch/s2112198"
OUTPUT="/exports/igmm/eddie/GenScotDepression/users/edavyson/antidep_project/revisions/output"
osca --befile ${SCRATCH}/mvals_GRM_uncorrected_29_06_OSCA_standard --extract-probe ${SCRATCH}/extra_probe.txt --make-efile --out ${OUTPUT}/selfrep_MDD_probe_17_09