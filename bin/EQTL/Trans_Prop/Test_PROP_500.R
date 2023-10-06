# Trans proportion statistics analysis on expressed transcription factors

# Read in the libraries
library(tidyverse)
library(dplyr)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
# Trans-eqtls
# These are removed if in high LD with cis-eQTLs and have FDR > 0.01
final_trans_FDR_01 <- read.csv("../../../work/cd/7fa9f208f71838895024288f03f321/Final_trans_FDR_01.txt", sep = ",")
# select unique SNPs for comparison
final_trans_FDR_01 <- final_trans_FDR_01[!duplicated(final_trans_FDR_01$snp),]

# Read in Transcription factors
TFs <- read.csv("../../../data/EQTL/Bos_taurus_TF.txt", sep = "\t") %>% dplyr::select(1,2,3,4,6)
CoF <- read.csv("../../../data/EQTL/Bos_taurus_Cof.txt", sep = "\t")
All <- rbind(TFs, CoF)



# Annotation file
# Filter for genes/TFs which have gene information available

annotation_file <- read.csv("../../../data/RNA_seq/Bovine_annotation_MF2.csv") %>% dplyr::select(1,2,3,6,7,8)
TFs <- right_join(All, annotation_file) # for Ratio analysis need location etc
TFs <- na.omit(TFs) # drop NAs -  were removed from our analysis during expression filtering


# Final data frame
# This will pick up number of genes close to our SNPs
Final_df <- data.frame(matrix(ncol = 13, nrow = 0))
colnames(Final_df) <- c("Distance",  "Total_Random", "Total_Trans", "Random_No_obs", "Random_Yes_obs", "Trans_No_obs","Trans_Yes_obs", "Random_No_exp", "Random_Yes_exp", "Trans_No_exp","Trans_Yes_exp", "Statistic", "P.value")


# Distance vector
# 10,000 - 1,000,000 bps
distance <- c(10000, 20000,30000,40000,50000,60000,70000,80000,90000,100000,200000,300000,400000,500000,600000,700000,800000,900000,1000000)



# AF for subsampling
# These AFs will be used to subsample a random set of snps
AF_trans <- read.csv("../../../work/dc/6598bc0de5eb5c613fa4018f4e6596/Trans.SNPs.AF.header.fixed.frq", sep = "\t", header = F)
AF_random <- read.csv("../../../work/dc/6598bc0de5eb5c613fa4018f4e6596/NonTrans.SNPs.AF.header.fixed.frq", sep = "\t", header = F)
AF <- rbind(AF_trans, AF_random)
colnames(AF) = c("CHR_trans", "POS_trans", "NALLELES", "NCHR", "REF:FREQ", "ALT:FREQ")
AF_long = separate(AF, `REF:FREQ`, into = c("ref", "AF_ref"), sep = ":")
AF_long = separate(AF_long, `ALT:FREQ`, into = c("alt", "AF_alt"), sep = ":") # AFs for ALL SNPs except the 2000 trans SNPs above
rm(AF)
rm(AF_trans)
rm(AF_random)

# SNP location
# Need this for our SNPs
snp_loc <- read.csv("../../../results/SNP_data/Imputation_Performance/eQTL_SNP_POS.txt", sep="\t") %>% dplyr::select(1,2, 3)
colnames(snp_loc)[1] <- "snps"
colnames(snp_loc)[2] <- "SNP_CHR"
colnames(snp_loc)[3] <- "SNP_POS"



# Gene location
gene_loc <- read.csv("../../../results/RNA-seq/Normalisation/gene_location.txt" ,sep = "\t") %>% dplyr::select(1,2,3,4) # get gene location
colnames(gene_loc)[1] <- "gene"
colnames(gene_loc)[2] <- "gene_CHR"
colnames(gene_loc)[3] <- "gene_START"
colnames(gene_loc)[4] <- "gene_END"

# Filtering for expressed genes
# Only consider TFs which are in our expressed genes
TFs <- TFs[TFs$Entrez_ID %in% gene_loc$gene,]
length(unique(TFs$Entrez_ID)) # FINAL total is 1660 expressed TFs / TF Cofactors
write.table(TFs, file = "Expressed_TFs.txt", row.names = F, col.names = T)


TFs$Chromosome = as.integer(TFs$Chromosome)



## The loop -------------------------------------------------------------------------------------------------------------------------
props_list <- list()

# set up an empty data frame which will collect our random SNPs in the loop
# set colnames equal to AF_long - this will make feasible to bind
Random_SNP_FINAL <- data.frame(matrix(ncol = 8, nrow = 0))
colnames(Random_SNP_FINAL) <- colnames(AF_long)


# Our groups here
group_1 <- rep(distance[1], 1000)
group_2 <- rep(distance[2], 1000)
group_3 <- rep(distance[3], 1000)
group_4 <- rep(distance[4], 1000)
group_5 <- rep(distance[5], 1000)
group_6 <- rep(distance[6], 1000)
group_7 <- rep(distance[7], 1000)
group_8 <- rep(distance[8], 1000)
group_9 <- rep(distance[9], 1000)
group_10 <- rep(distance[10], 1000)
group_11 <- rep(distance[11], 1000)
group_12 <- rep(distance[12], 1000)
group_13 <- rep(distance[13], 1000)
group_14 <- rep(distance[14], 1000)
group_15 <- rep(distance[15], 1000)
group_16 <- rep(distance[16], 1000)
group_17 <- rep(distance[17], 1000)
group_18 <- rep(distance[18], 1000)
group_19 <- rep(distance[19], 1000)
groups = c(group_1,
           group_2,
           group_3,
           group_4,
           group_5,
           group_6,
           group_7,
           group_8,
           group_9,
           group_10,
           group_11,
           group_12,
           group_13,
           group_14,
           group_15,
           group_16,
           group_17,
           group_18,
           group_19)


# Empty vector for each iteration
props = c()
props_trans = c()

# Loop for each distance
for (d in distance) {
  
  # change to numeric to ensure binding
  TFs$Chromosome = as.integer(TFs$Chromosome)
  
  # Number of trans SNPs (for calculating proportion)
  n_trans = length(final_trans_FDR_01$snp)
  
  # Merged to TFs by chromosome
  # Calculate distance to TFs
  Master_trans <- dplyr::left_join(final_trans_FDR_01, TFs, by = c("CHR_trans" = "Chromosome"))
  Master_trans$Distance_Start <- abs(as.numeric(Master_trans$POS_trans) - as.numeric(Master_trans$Start_location))
  Master_trans$Distance_End <- abs(as.numeric(Master_trans$POS_trans) - as.numeric(Master_trans$End_location))
  
  # Set up a column with distance of SNP to TFs
  Master_trans <- Master_trans %>% mutate(Distance = case_when(
    Master_trans$Distance_Start > Master_trans$Distance_End ~ Master_trans$Distance_End,
    Master_trans$Distance_End > Master_trans$Distance_Start ~ Master_trans$Distance_Start))
  
  # Get number of SNPs
  # Will use this to subset an equal number of SNPs (WITH REPLACEMENT) however many times
  
  
  Master_trans <- Master_trans %>% filter(Distance <= d) # filter for genes which are close to SNPs based on distance
  Master_trans <- Master_trans[order(Master_trans$snp),]
  Master_trans <- Master_trans %>% drop_na() # some genes do not have identifiers # doesnt matter as we are just pickng up TFs
  prop_trans <- length(unique(Master_trans$snp)) / n_trans
  props_trans <- append(props_trans, prop_trans)
  
  # set seed for reproducibility
  # Set up an empty dataframe for each permutation
  set.seed(1245)
  Random_final = replicate(1000, sample_n(AF_long, n_trans, replace = T), simplify = F)
  
  for (ele in 1:length(Random_final)) {
    
    # Reformat the Random SNP file
    Random_final[[ele]]$SNP_CHR <- as.numeric(Random_final[[ele]]$CHR_trans)
    Random_final[[ele]]$SNP_POS <- as.numeric(Random_final[[ele]]$POS_trans)
    Random_SNPs_final = Random_final[[ele]] %>% unite(., col = snps, CHR_trans, POS_trans, ref, alt, sep = ":")
    TFs$Chromosome <- as.numeric(TFs$Chromosome)
    
    
    # Join our random SNP set to TF by chromosome
    # Repeat process as above and calculate the distances
    Random_SNPs_final <- dplyr::left_join(Random_SNPs_final, TFs, by = c("SNP_CHR"="Chromosome"))
    Random_SNPs_final$Distance_Start <- abs(as.numeric(Random_SNPs_final$SNP_POS) - as.numeric(Random_SNPs_final$Start))
    Random_SNPs_final$Distance_End <- abs(as.numeric(Random_SNPs_final$SNP_POS) - as.numeric(Random_SNPs_final$End))
    Random_SNPs_final <- Random_SNPs_final %>% mutate(Distance = case_when(
      Random_SNPs_final$Distance_Start > Random_SNPs_final$Distance_End ~ Random_SNPs_final$Distance_End,
      Random_SNPs_final$Distance_End > Random_SNPs_final$Distance_Start ~ Random_SNPs_final$Distance_Start))
    
    # Calculate the proportion of SNPs close to a TF based on a distance
    Random_SNPs_final <- Random_SNPs_final %>% filter(Distance <= d) %>% drop_na() # filter for genes which are close to SNPs
    prop_ran <- length(unique(Random_SNPs_final$snps)) / n_trans
    print(c(d,ele))
    props <- append(props, prop_ran)
    
    # Save random SNPs for viewing
    Random_SNP_FINAL <- rbind(Random_SNP_FINAL, Random_SNPs_final)
  }
}

# Split up
test_split <- split(props, f=groups)

# Convert to a dataframe
df <- data.frame(matrix(unlist(test_split), nrow=length(test_split), byrow=TRUE),stringsAsFactors=FALSE)
# Bind in the distance
df <- cbind(distance, df)

# Convert to long format
df_long <- pivot_longer(df, cols = X1:X1000, names_to = "Iteration", values_to = "Prop_with_TFs")

# Getting the statistics
trans_df <- data.frame(distance, props_trans)
trans_df$Condition = "Trans-eQTL"
df_long_simple = df_long %>% dplyr::select(1,3)
df_long_simple$Condition = "Random"

MASTER_STAT <- df_long_simple %>% group_by(distance) %>% dplyr::summarize(Mean_NULL = mean(Prop_with_TFs), SD_NULL = sd(Prop_with_TFs))
MASTER_STAT$Obs <- props_trans
MASTER_STAT$Z <- (MASTER_STAT$Obs - MASTER_STAT$Mean_NULL) / MASTER_STAT$SD_NULL
MASTER_STAT$P_2_tailed <- (2*pnorm(q=MASTER_STAT$Z, lower.tail = F))
MASTER_STAT$P_1_tailed <-  pnorm(q=MASTER_STAT$Z, lower.tail = F)
MASTER_STAT$Bonferroni_2_tailed <- p.adjust(MASTER_STAT$P_2_tailed, method = "bonferroni")
MASTER_STAT$Bonferroni_1_tailed <- p.adjust(MASTER_STAT$P_1_tailed, method = "bonferroni")
MASTER_STAT$FDR_2_tailed <- p.adjust(MASTER_STAT$P_2_tailed, method = "BH")
MASTER_STAT$FDR_1_tailed <- p.adjust(MASTER_STAT$P_1_tailed, method = "BH")

# Write statistics
write.table(MASTER_STAT, file = "TF_PROP_STATS.csv", row.names = F, col.names = T)
colnames(trans_df) <- colnames(df_long_simple)


# Write the RAW STATS TO A TABLE
df_long_simple = rbind(trans_df, df_long_simple)
write.table(df_long_simple, file = "TF_Prop_RAW.txt", sep = "\t", row.names = F, col.names = T)
