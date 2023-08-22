library(UpSetR)
library(tidyverse)
library(ggridges)
library(ggplot2)
library(data.table)

args = commandArgs(trailing = T)

########################################################################################
########################################################################################
## 1. Upset plot of trans-eGenes
########################################################################################
########################################################################################

# Now onto upset plot of eGenes in reference groups
#args[1-3]
ALL_perm <- fread(args[1], header = F) %>% filter(V10 < 0.05) %>% select(1)
CONTROL_perm <- fread(args[2], header = F) %>% filter(V10 < 0.05) %>% select(1)
INFECTED_perm <- fread(args[3], header = F) %>% filter(V10 < 0.05) %>% select(1)
ALL_perm_genes <- unique(ALL_perm$V1)
CONTROL_perm_genes <- unique(CONTROL_perm$V1)
INFECTED_perm_genes <- unique(INFECTED_perm$V1)

listInput <- list(All = ALL_perm_genes, Control = CONTROL_perm$V1, Infected = INFECTED_perm$V1)

pdf(file=args[4], width = 12, height = 12, onefile=FALSE)
upset(fromList(listInput), order.by = "freq", sets.bar.color = c("#542788","#2166ac", "#b2182b"),
sets.x.label = "trans-eGenes", point.size = 4, line.size = 2,
mainbar.y.label = "trans-eGene intersections",
text.scale = 2.5, shade.alpha = 0.5)
dev.off()