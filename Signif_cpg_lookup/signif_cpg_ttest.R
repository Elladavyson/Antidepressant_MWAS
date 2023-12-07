################################################################################

# Set up libraries and options/files

################################################################################
.libPaths('/exports/igmm/eddie/GenScotDepression/users/edavyson/R/x86_64-pc-linux-gnu-library/4.1')
library(data.table)
library(dplyr)
library(ggplot2)
library(optparse)
library(tidyr)
library(ggpubr)

parse <- OptionParser()

# setting up options for the filepaths to the correct files
option_list <- list(
  make_option('--cohort', type='character', help="Cohort, ideally no spaces (for graphs and documentation)", action='store'),
  make_option('--id_column', type = 'character', default="IID", help = "Column names of identifier column in phenotype and covariate files", action = 'store'),
  make_option('--probes_dat', type = 'character', help = 'Methylation data file for the 7 significant CpGs', action = 'store'),
  make_option('--pheno', type = 'character', help = 'File path to antidepressant exposure phenotype file'),
  make_option('--outdir', type = 'character', help = 'The filepath for output directory', action = 'store')
)

args = commandArgs(trailingOnly=TRUE)
opt <- parse_args(OptionParser(option_list=option_list), args=args)

cohort <- opt$cohort
id_col <- opt$id_column # Vector of identifier columns 
probes_fp=opt$probes_dat # File path of methylation data file 
pheno_fp=opt$pheno # AD exposure (phenotype of cohort)
outdir <- opt$outdir # File path of output directory

################################################################################

# Read in data 

################################################################################

# standardised methylation levels (columns = ID, cgX, cgY, cgZ ....)

probes <- read.table(probes_fp, header = T)

# check the format (?)

## AD phenotype 
# support phenotype .csv files or .txt files 

if (endsWith(pheno_fp, '.csv')){
  ad_pheno <- read.csv(pheno_fp, header = T)
} else if (endsWith(pheno_fp, '.txt')) {
  ad_pheno <- read.table(pheno_fp, header = T)
} else {
  stop('Unsupported phenotype file, please provide the phenotype as a .csv or .txt file')
}

# check that there is an antidep column in the file 

if('antidep' %in% colnames(ad_pheno) == FALSE){
  stop('No antidep column in the phenotype file')
} else {
  print('antidep column in the phenotype file')
}

# remove missing values (if any?)
ad_pheno <- ad_pheno %>% filter(!is.na(antidep)) 
print(head(ad_pheno))
# logging phenotype characteristics 
print(paste0('Read in the Antidepressant exposure phenotype for ', cohort, ' : Number of cases: ',
             nrow(ad_pheno %>% 
                    filter(antidep==1)), 
             'Number of controls: ',
             nrow(ad_pheno%>% 
                    filter(antidep==0))))


################################################################################

# Merge data and visualise distributions

################################################################################

ad_pheno_meth <- merge(probes, ad_pheno, by = id_col)

ad_ph_meth_lg <- ad_pheno_meth %>% 
  pivot_longer(cols = starts_with("cg"), 
               names_to = c("CpG"),
               values_to = 'Mvals') %>% as.data.frame()

ad_ph_meth_lg <- ad_ph_meth_lg %>% mutate(antidep_name = ifelse(antidep == 0, 'Not AD exposed', 'AD exposed'))

# plot of the probe distributions in those exposed to ADs 
# and those not exposed 
# for deciding whether to use parametric/non-parametric test

dist_plots <- ggplot(ad_ph_meth_lg, aes(x = Mvals, fill = as.factor(antidep_name))) + 
  geom_histogram() + facet_grid(CpG~antidep_name)+  ylab('Count') + 
  xlab(paste0('Standardised M values')) + 
  theme_minimal() + theme(legend.position = 'none')

ggsave(filename= paste0(outdir, cohort, '_all_dists.png'), dist_plots, 
       width = 6, height = 10, device='png', dpi=300)

################################################################################

# Plotting functions 

################################################################################

data_summary <- function(x) {
  m <- mean(x)
  ymin <- m -sd(x)
  ymax <- m + sd(x)
  return(c(y=m, ymin = ymin, ymax = ymax))
}

probe_violins <- function(data, cpg) {
  print(cpg)
  # Welch's t-test of methylation levels in antidepressant exposed vs non exposed 
  cpg_ttest <- t.test(eval(as.symbol(cpg)) ~ antidep, data = data)
  cpg_tt_res <- data.table(CpG = cpg, t_stat = as.numeric(cpg_ttest$statistic),df = as.numeric(cpg_ttest$parameter), 
                                           p_val = as.numeric(cpg_ttest$p.value)) 
  plot <- ggplot(data, 
                 aes(x = as.factor(antidep), y = eval(as.symbol(cpg)),
                    color = as.factor(antidep), fill = as.factor(antidep)))+ 
    geom_violin(trim = FALSE, na.rm = TRUE, alpha= 0.3)+
    scale_color_brewer(palette = "Dark2", aesthetics= c("colour", "fill"),name = "Phenotype")+ 
    theme_classic()+
    geom_jitter(shape = 16, position=position_jitter(0.2), na.rm = TRUE)+
    stat_summary(fun.data =data_summary, shape = 23, color = "black", na.rm = TRUE) +
    labs(title = cpg) + 
    theme(legend.position= "none", text = element_text(size = 15),
    plot.title= element_text(face = "bold", hjust = 0.5))+
    scale_x_discrete(labels = c("Controls", "Cases")) + 
    ylim(min(data[,cpg]), max(data[,cpg] + 2*sd(data[,cpg])))+ 
    ylab('Standardised M vals') +
    labs(x = NULL)
    
    plot <- plot + annotate("text", x = 1.5, y =max(plot$data[,cpg]) +
                            2*sd(plot$data[,cpg]), label = paste0("p: ", signif(cpg_tt_res$p_val, 3)))+
      annotate("text", x = 1.5, y =max(plot$data[,cpg]) +
                 0.2*sd(plot$data[,cpg]), 
               label = ifelse(cpg_tt_res$p_val > 0.05, "NS", ifelse(cpg_tt_res$p_val < 0.01, "**", "*")))+
      geom_segment(aes(x = 1, xend = 2, y = max(plot$data[,cpg]) +
                   1*sd(plot$data[,cpg]), yend = max(plot$data[,cpg]) +
                     1*sd(plot$data[,cpg])), color = "black", linetype = "solid", size = 0.8)
  return(list(plot, cpg_tt_res))
}

################################################################################

# Plotting and saving probes 

################################################################################

violins <- list()
ttest_res <- data.frame()

for (i in 1:7) {
  cpg <- (ad_pheno_meth %>% select(starts_with("cg"))%>% colnames())[i]
  print(i)
  violin_func<- probe_violins(ad_pheno_meth, cpg)
  print(violin_func[1])
  print(violin_func[2])
  violins[[i]] <- violin_func[[1]]
  ggsave(filename = paste0(outdir, cohort, '_violin_', cpg, '.png'), plot = violin_func[[1]], 
         width = 7, height = 7, device='png', dpi=300)
  ttest_res <- rbind(ttest_res, violin_func[[2]])
}

all_violins <- ggarrange(plotlist = violins)
ggsave(filename = paste0(outdir, cohort, '_all_violins.png'), plot = all_violins, 
       width = 8, height = 10, device='png', dpi=300)

################################################################################

# Saving t-test results 

################################################################################

write.table(ttest_res, paste0(outdir, cohort, '_probe_ttest.tsv'), sep = '\t', row.names = F, quote = F)

