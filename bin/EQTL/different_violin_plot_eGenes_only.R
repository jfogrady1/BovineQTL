library(data.table)
library(tidyverse)
library(ggplot2)
library(UpSetR)
args = commandArgs(trailingOnly = T)
data_ALL <- fread("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/ALL.cis_qtl_fdr0.05.txt") %>% filter(is_eGene == T)
 data_CONTROL <- fread("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/CONTROL.cis_qtl_fdr0.05.txt")  %>% filter(is_eGene == T)
data_INFECTED <- fread("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/INFECTED.cis_qtl_fdr0.05.txt") %>% filter(is_eGene == T)
ALL_nominal <- fread("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/ALL.cis_qtl_pairs.ALL.unzipped.txt")
CONTROL_nominal <- fread("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/CONTROL.cis_qtl_pairs.ALL.unzipped.txt")
INFECTED_nominal <- fread("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/INFECTED.cis_qtl_pairs.ALL.unzipped.txt")

# Get rid of the replicated headings
ALL_nominal <- ALL_nominal %>% filter(phenotype_id != "phenotype_id")
CONTROL_nominal <- CONTROL_nominal %>% filter(phenotype_id != "phenotype_id")
INFECTED_nominal <- INFECTED_nominal %>% filter(phenotype_id != "phenotype_id")
# Plot gene wise threshold
threshold_plot_ALL <- data_ALL %>% select(1,19)
threshold_plot_CONTROL <- data_CONTROL %>% select(1,19)
threshold_plot_INFECTED <- data_INFECTED %>% select(1,19)

threshold_plot_ALL$Group <- "All"
threshold_plot_CONTROL$Group <- "Control"
threshold_plot_INFECTED$Group <- "Reactor"

threshold_plot <- rbind(threshold_plot_ALL, threshold_plot_CONTROL, threshold_plot_INFECTED)
ggplot(data = threshold_plot, aes(x = Group, y = -log10(pval_nominal_threshold), fill = Group)) + 
geom_violin(alpha = 0.5) + 
geom_boxplot(width = 0.07) + 
scale_fill_manual(values = c("#542788","#2166ac", "#b2182b")) + 
theme_bw() +
theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black"),
    axis.text.y = element_text(angle = 0, size = 15, colour = "black"),
    axis.title.y = element_text(size = 21, color = "black"),
    panel.grid.minor = element_blank(),
    legend.title = element_text(size = 15, color = "black"),
    legend.text = element_text(size = 15))
ggsave("Updated_violin_plot_thresholds.pdf", width = 10, height = 12, dpi = 600)
