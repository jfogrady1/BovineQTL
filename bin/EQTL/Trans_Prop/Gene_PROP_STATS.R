# Trans_ratio Statistic analysis
library(tidyverse)
library(dplyr)
library(ggplot2)
args = commandArgs(trailingOnly=TRUE)


# Trans-eqtls
# These are removed if in high LD with cis-eQTLs and have FDR > 0.01
final_trans_FDR_01 <- read.csv(args[1], sep = ",")
# select unique SNPs for comparison
final_trans_FDR_01 <- final_trans_FDR_01[!duplicated(final_trans_FDR_01$snp),]

# Read in Transcription factors
TFs <- read.csv(args[2], sep = "\t") %>% dplyr::select(1,2,3,4,6)
CoF <- read.csv(args[3], sep = "\t")
All <- rbind(TFs, CoF)



# Final data frame
Final_df <- data.frame(matrix(ncol = 13, nrow = 0))
colnames(Final_df) <- c("Distance",  "Total_Random", "Total_Trans", "Random_No_obs", "Random_Yes_obs", "Trans_No_obs","Trans_Yes_obs", "Random_No_exp", "Random_Yes_exp", "Trans_No_exp","Trans_Yes_exp", "Statistic", "P.value")


# Distance vector
# 10,000 - 1,000,000 bps
distance <- c(10000, 20000,30000,40000,50000,60000,70000,80000,90000,100000,200000,300000,400000,500000,600000,700000,800000,900000,1000000)

# AF for subsampling
# These AFs will be used to subsample a random set of snps
AF_trans <- read.csv(args[4], sep = "\t", header = F)
AF_random <- read.csv(args[5], sep = "\t", header = F)
AF <- rbind(AF_trans, AF_random)
colnames(AF) = c("CHR_trans", "POS_trans", "NALLELES", "NCHR", "REF:FREQ", "ALT:FREQ")
AF_long = separate(AF, `REF:FREQ`, into = c("ref", "AF_ref"), sep = ":")
AF_long = separate(AF_long, `ALT:FREQ`, into = c("alt", "AF_alt"), sep = ":") # AFs for ALL SNPs except the 2000 trans SNPs above
rm(AF)
rm(AF_trans)
rm(AF_random)


# SNP location
# Need this for our SNPs
snp_loc <- read.csv(args[6], sep="\t") %>% dplyr::select(1,2, 3)
colnames(snp_loc)[1] <- "snps"
colnames(snp_loc)[2] <- "SNP_CHR"
colnames(snp_loc)[3] <- "SNP_POS"



# Gene location
gene_loc <- read.csv(args[7], sep = "\t") %>% dplyr::select(1,2,3,4) # get gene location
colnames(gene_loc)[1] <- "gene"
colnames(gene_loc)[2] <- "gene_CHR"
colnames(gene_loc)[3] <- "gene_START"
colnames(gene_loc)[4] <- "gene_END"



gene_loc$gene_CHR <- as.numeric(gene_loc$gene_CHR)
#gene_loc <- gene_loc %>% filter(!(gene %in% All$Entrez_ID))
#dim(gene_loc)
## The loop -------------------------------------------------------------------------------------------------------------------------
props_list <- list()

# set up an empty data frame which will collect our random SNPs in the loop
# set colnames equal to AF_long - this will make feasible to bind
Random_SNP_FINAL <- data.frame(matrix(ncol = 8, nrow = 0))
colnames(Random_SNP_FINAL) <- colnames(AF_long)

group_1 <- rep(distance[1], 10000)
group_2 <- rep(distance[2], 10000)
group_3 <- rep(distance[3], 10000)
group_4 <- rep(distance[4], 10000)
group_5 <- rep(distance[5], 10000)
group_6 <- rep(distance[6], 10000)
group_7 <- rep(distance[7], 10000)
group_8 <- rep(distance[8], 10000)
group_9 <- rep(distance[9], 10000)
group_10 <- rep(distance[10], 10000)
group_11 <- rep(distance[11], 10000)
group_12 <- rep(distance[12], 10000)
group_13 <- rep(distance[13], 10000)
group_14 <- rep(distance[14], 10000)
group_15 <- rep(distance[15], 10000)
group_16 <- rep(distance[16], 10000)
group_17 <- rep(distance[17], 10000)
group_18 <- rep(distance[18], 10000)
group_19 <- rep(distance[19], 10000)
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

props = c()
props_trans = c()
for (d in distance) {
  n_trans = length(final_trans_FDR_01$snp)
  Master_trans <- dplyr::left_join(final_trans_FDR_01, gene_loc, by = c("CHR_trans" = "gene_CHR"))
  Master_trans$Distance_Start <- abs(as.numeric(Master_trans$POS_trans) - as.numeric(Master_trans$gene_START))
  Master_trans$Distance_End <- abs(as.numeric(Master_trans$POS_trans) - as.numeric(Master_trans$gene_END))
  Master_trans <- Master_trans %>% mutate(Distance = case_when(
    Master_trans$Distance_Start > Master_trans$Distance_End ~ Master_trans$Distance_End,
    Master_trans$Distance_End > Master_trans$Distance_Start ~ Master_trans$Distance_Start))
  # Get number of SNPs
  # Will use this to subset an equal number of SNPs however many times
  
  # finish with our observed set
  Master_trans <- Master_trans %>% filter(Distance <= d) # filter for genes which are close to SNPs based on distance
  Master_trans <- Master_trans[order(Master_trans$snp),]
  Master_trans <- Master_trans %>% drop_na() # some genes do not have identifiers # doesnt matter as we are just pickng up TFs
  prop_trans <- length(unique(Master_trans$snp)) / n_trans
  props_trans <- append(props_trans, prop_trans)
  
  head(Master_trans)
  set.seed(1245)
  
  Random_final = replicate(10000, sample_n(AF_long, n_trans, replace = T), simplify = F)
  for (ele in 1:length(Random_final)) {
    # Reformat the Random SNP file
    Random_final[[ele]]$SNP_CHR <- as.numeric(Random_final[[ele]]$CHR_trans)
    Random_final[[ele]]$SNP_POS <- as.numeric(Random_final[[ele]]$POS_trans)
    Random_SNPs_final = Random_final[[ele]] %>% unite(., col = snps, CHR_trans, POS_trans, ref, alt, sep = ":")
    Random_SNPs_final <- dplyr::left_join(Random_SNPs_final, gene_loc, by = c("SNP_CHR"= "gene_CHR"))
    Random_SNPs_final$Distance_Start <- abs(as.numeric(Random_SNPs_final$SNP_POS) - as.numeric(Random_SNPs_final$gene_START))
    Random_SNPs_final$Distance_End <- abs(as.numeric(Random_SNPs_final$SNP_POS) - as.numeric(Random_SNPs_final$gene_END))
    Random_SNPs_final <- Random_SNPs_final %>% mutate(Distance = case_when(
      Random_SNPs_final$Distance_Start > Random_SNPs_final$Distance_End ~ Random_SNPs_final$Distance_End,
      Random_SNPs_final$Distance_End > Random_SNPs_final$Distance_Start ~ Random_SNPs_final$Distance_Start))
    
    Random_SNPs_final <- Random_SNPs_final %>% filter(Distance <= d) %>% drop_na() # filter for genes which are close to SNPs
    prop_ran <- length(unique(Random_SNPs_final$snps)) / n_trans
    print(c(d, ele))
    props <- append(props, prop_ran)
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
df_long <- pivot_longer(df, cols = X1:X10000, names_to = "Iteration", values_to = "Prop_with_TFs")




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
write.table(MASTER_STAT, file = "Gene_PROP_STATS.csv", row.names = F, col.names = T)
colnames(trans_df) <- colnames(df_long_simple)
head(AF_long)

# Write the RAW STATS TO A TABLE
df_long_simple = rbind(trans_df, df_long_simple)
write.table(df_long_simple, file = "Gene_PROP_RAW.txt", sep = "\t", row.names = F, col.names = T)
