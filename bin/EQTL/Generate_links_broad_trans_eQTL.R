
# Generate links for circos plot
library(data.table)
library(tidyverse)
args = commandArgs(trailingOnly = TRUE)
data <- fread(args[1]) %>% filter(V10 < 0.05)
# determine groups
data$intra = if_else(data$V2 == data$V5, TRUE, FALSE)
table(data$intra)
#FALSE  TRUE 
#4638   585

# get colour
data$linetype = ifelse(data$intra == TRUE, paste0("color=black_a1,thickness=2p"), paste0("color=vvdpred_a1,thickness=2p"))

# now get columns into right format
colnames(data) <- c("gene", "gene_chr", "gene_pos", "variant", "variant_chr", "variant_pos", "pval", "r2", "corr", "FDR", "intra", "linetype")
data$gene_pos_2 <- data$gene_pos +1
data$variant_pos_2 <- data$variant_pos + 1
data_plot <- data %>% select(gene_chr, gene_pos, gene_pos_2,variant_chr, variant_pos,variant_pos_2,linetype)
data_plot$gene_chr <- gsub("^", "btau", data_plot$gene_chr)
data_plot$variant_chr <- gsub("^", "btau", data_plot$variant_chr)
write.table(data_plot, args[2], sep = " ", row.names = F, col.names = F, quote = F)
