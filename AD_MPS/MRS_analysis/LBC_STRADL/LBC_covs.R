
################################################################################

## Creating LBC covariate files 

################################################################################

.libPaths('/exports/igmm/eddie/GenScotDepression/users/edavyson/R/x86_64-pc-linux-gnu-library/4.1')
library(readr)
library(data.table)
library(dplyr)
library(tidyr)

################################################################################

## Reading in the data 

################################################################################
setwd('/exports/eddie/scratch/s2112198/')
mvals_36 <- readRDS('LBC36_mvals_wave1.rds')
mvals_21 <- readRDS('LBC21_mvals_wave1.rds')
targets <- readRDS('targets_3489_bloodonly.rds')

# getting the smoking CpG

smoke_cg36 <- mvals_36 %>% select(ID_raw, cg05575921)
smoke_cg21 <- mvals_21 %>% select(ID_raw, cg05575921)

# Reading in the genetic PCs (TBC from Daniel)

################################################################################

## Plotting the distribution of the cell count variables 

################################################################################

targets_filter <- targets %>% filter(!(neut == 9999| lymph == 9999|mono == 9999|eosin == 9999| baso == 9999))
targets_filter_no_outliers <- targets_filter %>%
  filter(across(c(neut, lymph, mono, eosin, baso), ~ between(., quantile(., 0.25) - 1.5 * IQR(.), quantile(., 0.75) + 1.5 * IQR(.))))

targets_long <- targets_filter %>% 
  select(c(ID_raw, cohort, neut, lymph, mono, eosin, baso)) %>%
  pivot_longer(., c(neut, lymph, mono, eosin, baso), names_to = "Cell_count", values_to = "value")

targets_long_no_outliers <- targets_filter_no_outliers %>% 
  select(c(ID_raw, cohort, neut, lymph, mono, eosin, baso)) %>%
  pivot_longer(., c(neut, lymph, mono, eosin, baso), names_to = "Cell_count", values_to = "value")


cell_boxplots <- ggplot(targets_long, aes(x = Cell_count, y =value, fill = Cell_count)) + 
  geom_boxplot() + 
  scale_fill_viridis_d(alpha = 0.6) +
  theme_minimal() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 11)
  )+
  ggtitle("Box plots of cell counts in LBC1921 and LBC1936") + 
  xlab("") + facet_wrap(~cohort)

cell_boxplots_nooutliers <- ggplot(targets_long_no_outliers, aes(x = Cell_count, y =value, fill = Cell_count)) + 
  geom_boxplot() + 
  scale_fill_viridis_d(alpha = 0.6) +
  theme_minimal() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 11)
  )+
  ggtitle("Box plots of cell counts in LBC1921 and LBC1936") + 
  xlab("") + facet_wrap(~cohort)



cell_hists <- ggplot(targets_long, aes(x = value)) + geom_histogram() + facet_grid(cohort~Cell_count) + theme_minimal()

ggsave(filename = '/exports/eddie/scratch/s2112198/LBC_cellcount_hists.png', plot = cell_hists, width = 8, height = 6, device = "png", dpi = 300)
ggsave(filename = '/exports/eddie/scratch/s2112198/LBC_cellcount_boxplots.png', plot = cell_boxplots, width = 8, height = 6, device = "png", dpi = 300)
ggsave(filename = '/exports/eddie/scratch/s2112198/LBC_cellcount_boxplots_nooutliers.png', plot = cell_boxplots_nooutliers, width = 8, height = 6, device = "png", dpi = 300)

################################################################################

## Selecting columns from target file (age, sex, cell counts)

################################################################################

lbc36_covs <- targets %>% 
  filter(cohort == 'LBC36' & WAVE == 1) %>%
  select(ID_raw, age, neut, lymph, mono, eosin, baso, array, sex) %>%
  mutate(sex_coded = ifelse(sex == 'M', 1, 0))

lbc21_covs <- targets %>% 
  filter(cohort == 'LBC21' & WAVE == 1) %>%
  select(ID_raw, age, neut, lymph, mono, eosin, baso, array, sex) %>%
  mutate(sex_coded = ifelse(sex == 'M', 1, 0))


lbc36_covs_no_outliers <- targets_filter_no_outliers %>% 
  filter(cohort == 'LBC36' & WAVE == 1) %>%
  select(ID_raw, age, neut, lymph, mono, eosin, baso, array, sex) %>%
  mutate(sex_coded = ifelse(sex == 'M', 1, 0))

lbc21_covs_no_outliers <- targets_filter_no_outliers %>% 
  filter(cohort == 'LBC21' & WAVE == 1) %>%
  select(ID_raw, age, neut, lymph, mono, eosin, baso, array, sex) %>%
  mutate(sex_coded = ifelse(sex == 'M', 1, 0))


################################################################################

## Merging together

################################################################################

lbc36_covs <- merge(lbc36_covs, smoke_cg36, by.x = 'ID_raw')
lbc21_covs <- merge(lbc21_covs, smoke_cg21, by.x = 'ID_raw')
lbc36_covs_no_outliers <- merge(lbc36_covs_no_outliers, smoke_cg36, by.x = 'ID_raw')
lbc21_covs_no_outliers <- merge(lbc21_covs_no_outliers, smoke_cg21, by.x = 'ID_raw')

################################################################################

## Saving 

################################################################################

write.table(lbc36_covs, 'LBC1936_covs.txt', row.names = F, quote = F)
write.table(lbc21_covs, 'LBC1921_covs.txt', row.names = F, quote = F)
write.table(lbc36_covs_no_outliers, 'LBC1936_covs_nooutliers.txt', row.names = F, quote = F)
write.table(lbc21_covs_no_outliers, 'LBC1921_covs_nooutliers.txt', row.names = F, quote = F)
