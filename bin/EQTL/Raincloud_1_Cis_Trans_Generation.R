# 1: Generating SNPs for LD calculation to cis SNPs of the same gene

# Arguments
# ${cis_eqtls} ${trans_eqtls} ${gene_loc} ${snp_loc}

## Packages
library(tidyverse)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggdist)
library(gghalves)
library(tidyquant)
args = commandArgs(trailingOnly=TRUE)

# load in data
# Arg[1] = cis data
# Arg[2] = trans data (both significant) = output from eQTL analysis (order and select both)
cis_data <- read.csv(args[1], sep = "\t")
trans_data <- read.csv(args[2], sep = "\t") # all data with most points - used also in TWAS so need to plot that


cis_data <- cis_data %>% filter(cis_data$FDR < 0.01)
trans_data <- trans_data %>% filter(trans_data$FDR < 0.01)

# ONLY consider trans-SNPs on different chromosomes
# need gene location information
# need snp location information also
# Arg[3] = gene location file used in normalisation procedure
gene_loc <- read.csv(args[3], sep = "\t") %>% select(1,2) # get gene location
colnames(gene_loc)[1] <- "gene"
colnames(gene_loc)[2] <- "gene_CHR"
trans_data <- left_join(trans_data, gene_loc) # join everything


# snp location
# Arg[4] = SNP location used in EQTL analysis
snp_loc <- read.csv(args[4], sep="\t") %>% select(1,2)
colnames(snp_loc)[1] <- "snps"
colnames(snp_loc)[2] <- "SNP_CHR"
trans_data <- left_join(trans_data, snp_loc)

dim(trans_data)
# Filter for SNPs which are on a different chromsome
trans_data <- trans_data %>% filter(gene_CHR != SNP_CHR) # find SNPs which arent on same CHR as gene
dim(trans_data) # putative eQTLs on different chromosomes
length(unique(trans_data$snps)) # numbers

# Have SNPs, now need the genes
# Get gene list from trans-eQTLs
gene_list <- unique(trans_data$gene)

# filter cis data for same genes
cis_data <- cis_data %>% filter(gene %in% gene_list) # note here we have significant trans SNPs not associated with cis-SNPs

# cross reference
gene_list <- gene_list[gene_list %in% cis_data$gene]
length(gene_list)
# filter trans data
trans_data <- trans_data %>% filter(gene %in% gene_list)

length(gene_list) # 503 genes with a trans association on different chromsome and also cis association

# generate list of SNPs to compare with
snp_list <- c(trans_data$snps, cis_data$snps)
snp_list <- as.data.frame(snp_list)
length(snp_list$snp_list) - length(unique(snp_list$snp_list)) # 4848 SNPs with cis and trans association
dim(snp_list)
write.csv(snp_list, "TransCisLDComparisonSNPS.txt", quote = F, row.names = F, col.names = F)
write.csv(trans_data$snps, "Trans_LD_comparisonSNPs.txt", quote = F, row.names = F, col.names = F) # will extract from column 1
write.csv(cis_data$snps, "Cis_LD_comparisonSNPs.txt", quote = F, row.names = F, col.names = F) # will extract from column 2 of outputted file
