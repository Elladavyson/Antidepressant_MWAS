# Investigating sex specific effects 
.libPaths(c(.libPaths(), '/exports/eddie3_homes_local/s2112198/R/x86_64-pc-linux-gnu-library/4.4'))
library(dplyr)
library(readr)
library(data.table)
library(ggplot2)
library(ggrepel)
library(tidyr)
library(ggpubr)
# Self-reported antidepressant exposure 
mwas_dir <- "/exports/igmm/eddie/GenScotDepression/users/edavyson/antidep_project/revisions/output/MWAS_sex_segregated/"
# Read in the MWAS results
male <- read.table(paste0(mwas_dir, "MOA_output/selfrep_pheno3_GRMunadj_ORM_residph_male_16_09.moa"), header = T)
female <- read.table(paste0(mwas_dir, "MOA_output/selfrep_pheno3_GRMunadj_ORM_residph_female_16_09.moa"), header = T)
# Formatting 
female[female$Probe == 'cg07023494', 'Chr'] <- 7
female[female$Probe == 'cg07023494', 'bp'] <- 158965466
female <- female %>% dplyr::rename(Name = Probe)
male[male$Probe == 'cg07023494', 'Chr'] <- 7
male[male$Probe == 'cg07023494', 'bp'] <- 158965466
male <- male %>% dplyr::rename(Name = Probe)
female <- female %>% select(-c(Gene, Orientation))
male <- male %>% select(-c(Gene, Orientation))

# MWAS plot function 
MWAS_manhattan <- function(OSCA_results, plottitle, annot_level = 9.42e-08, start_var=NULL, end_var = NULL){
  OSCA_results_format <- OSCA_results %>% 
  #compute chromosome size 
  group_by(Chr) %>%
  reframe(chr_len = max(bp)) %>%
  
  #Calculate the cumulative position of each chromosome 
  mutate(tot = cumsum(as.numeric(chr_len))-as.numeric(chr_len)) %>%
  dplyr::select(-chr_len) %>%
  
  #add this info to the original data set of sum stats 
  left_join(OSCA_results, ., by = c('Chr'='Chr')) %>%
  
  #add the cumulative position of each SNP 
  arrange(Chr,bp) %>% 
  mutate( BPcum = bp + tot)

axisdf = OSCA_results_format %>% group_by(Chr) %>% reframe(center=(max(BPcum) +min(BPcum))/2)
OSCA_results_signif <- OSCA_results_format %>% dplyr::filter(p < annot_level)
MWAS_plot <- ggplot(OSCA_results_format, aes(x = BPcum, y = -log10(p)))+geom_point(aes(color = as.factor(Chr)), alpha = ifelse(OSCA_results_format$p < annot_level, 1,0.8), size = 1.3) + 
  scale_color_manual(values = rep(c("darkgray", "coral"), 22)) + 
  #custom x axis 
  
  scale_x_continuous(label=axisdf$Chr, breaks= axisdf$center) + 
  scale_y_continuous(expand = c(0,0), limits = c(0,max(-log10(OSCA_results$p))+2)) +
  geom_text_repel(data = OSCA_results_signif, aes(label = Name))+
  geom_hline(yintercept = -log10(9.42e-08), linetype = 'dashed') +
  theme_bw() + theme(legend.position = "none", panel.border = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), plot.title= element_text(face = "bold", hjust = 0.5)) + 
  labs(x = 'Chromosome', title = plottitle)
return(MWAS_plot)
}

female_MWAS <- MWAS_manhattan(female, 'Self report antidepressant exposure: Females (n = 9710)')
male_MWAS <- MWAS_manhattan(male, 'Self report antidepressant exposure: Males (n = 6821)')
ggsave(filename = paste0(mwas_dir, 'selfrep_MWAS_female.png'), plot = female_MWAS, width = 8, height = 6, device='png', dpi=300)
ggsave(filename = paste0(mwas_dir, 'selfrep_MWAS_male.png'), plot = male_MWAS, width = 8, height = 6, device='png', dpi=300)

# Comparison plot 

comp_plot <- function(res1, res2, analysis1, analysis2, phenotype, probes){
   # add the analysis to the relevant column names for long format

  print('Appending analysis to column names')
  colnames(res1)[c(4,5,6)] <- paste0(colnames(res1)[c(4,5,6)], '_', analysis1)
  colnames(res2)[c(4,5,6)] <- paste0(colnames(res2)[c(4,5,6)], '_', analysis2)
  # merge them together 
  print('Merging the results together')
  results_mrg <- merge(res1, res2, by = c('Chr', 'Name', 'bp') )
print('Making long format')
# making into long format
results_long <- results_mrg %>% 
  pivot_longer(-c(Chr, Name, bp),
               names_to = c(".value", "Analysis"), 
               names_sep = "_") %>% as.data.frame()

# add standard errors 
print('Calculating standard errors')
results_long <- results_long %>% mutate(lowSE = b - se, highSE = b + se) %>% mutate(highCI = b + (1.96*se), lowCI = b - (1.96*se))

# subset the significant hits from the antidep_pheno1 (not subset to MDD)
print('Plotting')
res1_sig <- res1 %>% dplyr::filter(Name %in% probes)
num_signif <- nrow(res1_sig)
comp_plt <- ggplot(results_long %>% dplyr::filter(Name %in% res1_sig$Name), aes(
  x = b, y = reorder(Name, b), color = Analysis, shape = p < 0.05/num_signif)) + 
  geom_point(size = 2, position= position_dodge(width=0.5)) + 
  geom_errorbarh(aes(xmin = lowCI, xmax = highCI), height = 0, position= position_dodge(width=0.5)) + 
  theme_minimal() +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'black') + 
  labs(y = 'CpG', x = 'Standardised Effect Size', title = paste0(
         'Probes significant in the', phenotype), color = 'Phenotype')+guides(shape=guide_legend(title=paste0('p < 0.05/', num_signif)))

return(list(comp_plt, results_long))
}

sex_comp <- comp_plot(male, female, "Male", "Female", phenotype = "Self-report", probes = female %>% filter(p < 9.42e-08) %>% pull(Name))
sex_comp_plot <- sex_comp[[1]]
ggsave(filename = paste0(mwas_dir, 'selfrep_sex_comp_signif.png'), plot = sex_comp[[1]], width = 8, height = 6, device='png', dpi=300)

sex_comp_df <- sex_comp[[2]]
sex_comp_df_wide <- sex_comp_df %>% 
pivot_wider(names_from = "Analysis", 
values_from = c("b", "se", "p", "lowSE", "highSE", "highCI", "lowCI"), 
names_sep = "_")


sex_comp_df_wide <- sex_comp_df_wide %>% mutate(b_sexdiff = b_Male - b_Female)

# Correlation plotn 

beta_cor_sex <- ggplot(sex_comp_df_wide, aes(x = b_Male, y = b_Female)) + 
geom_point(alpha = 0.6)+geom_abline(slope = 1, intercept = 0, color = 'red', linetype = 'dashed') + 
stat_cor()+ theme_minimal() + 
theme(legend.position = 'right', plot.title= element_text(face = "bold", hjust = 0.5), axis.line = element_line(color = "black"), text = element_text(size = 12))+
ggtitle('Self-report MWAS: Sex stratified analysis') + xlim(-0.15, 0.25) + ylim(-0.15, 0.25) + labs(x = 'CpG Effects: Male (n = 6, 821)', y = 'CpG Effects: Female (n =  9, 710)')
# Calculate the Z scores of the effect estiamtes differences 
mean(sex_comp_df_wide$b_sexdiff)
sd(sex_comp_df$b_sexdiff)