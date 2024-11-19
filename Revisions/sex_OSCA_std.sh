# Create OSCA DNA methylation datasets which are standardised separately for females and males
# From the same OSCA DNA methylation object which was standardised (males + females together) for the main analysis
module load igmm/apps/osca/0.46
osca --befile mvals_GRM_uncorrected_29_06_OSCA_noAHRR --keep selfrep_females.list --std-probe --make-bod mvals_GRMuncor_female_standardised_23_09
osca --befile mvals_GRM_uncorrected_29_06_OSCA_noAHRR --keep selfrep_males.list --std-probe --make-bod mvals_GRMuncor_male_standardised_23_09