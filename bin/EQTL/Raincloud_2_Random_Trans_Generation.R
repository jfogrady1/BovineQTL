# 2. Generating random SNPs on same chromosome with same AF for calculation of LD with cis SNPs

# Arguments
#${cis_eqtls} ${trans_eqtls} ${gene_loc} ${snp_loc} Trans.SNPs.AF.header.fixed.frq NonTrans.SNPs.AF.header.fixed.frq

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
# Arg[1] = cis data (ALL animals)
# Arg[2] = trans data (both significant) = output from eQTL analysis (order and select both)
cis_data <- read.csv(args[1], sep = "\t")
trans_data <- read.csv(args[2], sep = "\t") # all data with most points - used also in TWAS so need to plot that

cis_data <- cis_data %>% filter(cis_data$FDR < 0.01)
trans_data <- trans_data %>% filter(trans_data$FDR < 0.01)


# ONLY consider trans-SNPs on different chromosomes
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


###################################################################################################################################################################################
###################################################################################################################################################################################
################## Generation of Random trans SNPs and looking at LD between them and cis SNPs ####################################################################################
###################################################################################################################################################################################
###################################################################################################################################################################################


trans_data <- trans_data %>% select(1,7,8)
dim(trans_data) # need to generate 2627 (2616 unique)
length(unique(trans_data$SNP_CHR)) # from all 29 autosomes.

# 1. need to get allele frequency for current SNPs
# This command is run on the shell in nextflow
#system(paste("vcftools --gzvcf TransCisLDSNPs.vcf.gz --snps Trans_LD_comparisonSNPs.txt --freq --out Trans.SNPs.AF"))
#system(paste("tail -n+2 Trans.SNPs.AF.frq  > Trans.SNPs.AF.header.fixed.frq"))

# Read in the file produced by vcf tools
# Change the name of the columns
# Arg[5] <- the AF for Trans SNPs
AF = read.csv(args[5], sep = "\t", header = F)
colnames(AF) = c("CHROM", "POS_trans", "NALLELES", "NCHR", "REF:FREQ", "ALT:FREQ")

# Edit the trans_data to get it into long format(really to just get CHR and POSITION)
# Format columns to make sure they are the same for joining
trans_long = separate(trans_data, snps, into = c("CHR_trans", "POS_trans", "REF", "ALT"), sep = ":") 
trans_long$CHR_trans <- as.character(trans_long$CHR_trans)
trans_long$POS_trans <- as.character(trans_long$POS_trans)
AF$POS_trans <- as.character(AF$POS_trans)

# Join the trans SNPs with their Allele frequencies 
# Format the last 2 columns to extend out
# Select columns we want
trans_long = left_join(trans_long, AF)
trans_long = separate(trans_long, `REF:FREQ`, into = c("ref", "AF_ref"), sep = ":")
trans_long = separate(trans_long, `ALT:FREQ`, into = c("alt", "AF_alt"), sep = ":")
trans_long <- trans_long %>% select(1,2,3,4,11,13)
trans_long$Type = "Trans-eQTL"

# Now have AFs for trans SNPs which we will compare LD against top cis SNPs (done in LDR"FULL.interchrom.geno.ld)

# 2. Need to sample random SNPs from same chromosome and which have same allele frequency as each SNP
# Will subset out the trans SNPs which we have and compute the AF for all other SNPs
# Then will match SNPs with MAF on the same chromosome in a loop and obtain a corresponding SNP for each trans-eQTL
# Run on shell
#system(paste("vcftools --vcf ../../../data/SNP_data/imputation/Final_data_for_eqtl/Axiom.IMPUTED.FILTERED.simple.codes.EQTL.vcf --exclude Trans_LD_comparisonSNPs.txt --freq --out NonTRANS.SNPs.AF"))
#system(paste("tail -n+2 NonTRANS.SNPs.AF.frq  > NonTrans.SNPs.AF.header.fixed.frq")) # format the header

# Format the allele frequency file as above
# ARG[6] <- AF of NonTrans SNPs
AF = read.csv(args[6], sep = "\t", header = F)
colnames(AF) = c("CHR_trans", "POS_trans", "NALLELES", "NCHR", "REF:FREQ", "ALT:FREQ")
AF_long = separate(AF, `REF:FREQ`, into = c("ref", "AF_ref"), sep = ":")
AF_long = separate(AF_long, `ALT:FREQ`, into = c("alt", "AF_alt"), sep = ":") # AFs for ALL SNPs except the 2000 trans SNPs above


####################################################################
######### Generating RANDOM SNPS ###################################
####################################################################

# set up an empty data frame which will collect our random SNPs
# set colnames equal to trans long - this will make feasible to bind
Random_SNPs <- data.frame(matrix(ncol = 7, nrow = 0))
colnames(Random_SNPs) <- colnames(trans_long)

# set seed for reproducibility
set.seed(17897)

# for each row in our trans_long file (file which contains trans eQTLs on diff chromosomes)
for (row in 1:nrow(trans_long)) {
  chr = trans_long[row, "CHR_trans"] # get chromosome
  af = trans_long[row, "AF_ref"] # get reference allele freqeuncy
  
  # Filter our master AF file (excluding the trans SNPs) for SNPs with the same CHR and AF as our trans (in variables above)
  potent = AF_long %>% filter(CHR_trans == chr & AF_ref == af) %>% as.data.frame() # save potential matches as df
  
  # Select a random SNP from this filtered data set
  random_row = sample_n(potent, 1)
  # Call it a Random SNP
  random_row$Type = "Random"
  # Append it to our empty data frame set up before loop started
  Random_SNPs = rbind(Random_SNPs, random_row)
} 

# Reformat the Random SNP file
Random_SNPs_final = Random_SNPs %>% unite(., col = snps, CHR_trans, POS_trans, ref, alt, sep = ":")
# Select columns required
Random_SNPs_final <- Random_SNPs_final %>% select(1,6)




# Bind our random trans SNPs with cis SNPs
# Write to a file
snp_list <- c(Random_SNPs_final$snps, cis_data$snps) # from above
snp_list <- as.data.frame(snp_list)
write.csv(snp_list, "RandomTransCisLDComparisonSNPS.txt", quote = F, row.names = F, col.names = F)
write.csv(Random_SNPs_final$snps, "RandomTrans_LD_comparisonSNPs.txt", quote = F, row.names = F, col.names = F) # will extract from column 1
