# Script to see if remaining trans-eQTLs reside close to TFs or Co TFs
library(data.table)
library(tidyverse)
library(vcfR)
library(ggridges)
# Read in inter and intra chromosomal trans-eQTLs
inter = fread("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/ALL_interchromosomal_LD_corrected.txt")
intra = fread("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/ALL_intrachromosomal_LD_corrected.txt")
intra$V11 <- "TRUE"

trans <- rbind(intra,inter)
colnames(trans) <- c("Gene", "gene_chr", "gene pos", "SNP", "SNP_chr", "SNP_pos", "Pvalue", "V8", "correlation", "FDR", "intrachromosomal")

dim(trans) # 3799 trans-eVariants associated with 49 genes
length(unique(trans$Gene))

# Now select top SNP
top_trans <- trans %>% group_by(Gene) %>% filter(Pvalue == min(Pvalue))
top_trans <- top_trans[order(top_trans$Pvalue, decreasing = F),]
top_trans <- top_trans[!(duplicated(top_trans$Pvalue)),]

# Remove ties
table(top_trans$Pvalue)

table(top_trans$intrachromosomal)
length(unique(trans_cis$Gene))


length(unique(top_trans$SNP))

# Now read in TFs and Co TFs
TF <- fread("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/TWAS/Bos_taurus_TF.txt") %>% select(Ensembl, Symbol)
coTF <- fread("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/TWAS/Bos_taurus_Cof.txt") %>% select(Ensembl, Symbol)
head(TF)


bed_file <- fread("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/ALL.expr_tmm_inv.bed.gz") %>% select(1:4)

TF <- TF %>% filter(Ensembl %in% bed_file$phenotype_id)
TFs <- left_join(TF, bed_file, by = c("Ensembl" = "phenotype_id"))
dim(TFs) # 973 expressed TFs/CoTFs

top_trans <- top_trans %>% select(4,5,6)
top_trans <- top_trans[,-1]
colnames(top_trans) <- c("SNP", "CHR", "SNP_POS")
head(top_trans)
colnames(TFs) <- c("Ensembl", "Symbol", "Chromosome", "POS", "POS2")
head(TFs)
# Final data frame
# This will pick up number of genes close to our SNPs
Final_df <- data.frame(matrix(ncol = 13, nrow = 0))
colnames(Final_df) <- c("Distance",  "Total_Random", "Total_Trans", "Random_No_obs", "Random_Yes_obs", "Trans_No_obs","Trans_Yes_obs", "Random_No_exp", "Random_Yes_exp", "Trans_No_exp","Trans_Yes_exp", "Statistic", "P.value")


# Distance vector
# 10,000 - 1,000,000 bps
distance <- c(10000, 20000,30000,40000,50000,60000,70000,80000,90000,100000,200000,300000,400000,500000,600000,700000,800000,900000,1000000)

# Load in vcf file of all SNPs
vcf_data <- read.vcfR("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/ALL_IMPUTED_UPDATED.vcf.gz")
snps_geno <- vcf_data@gt
head(vcf_data@fix)
snps_id <- vcf_data@fix %>% as.data.frame() %>% select(ID, CHROM, POS)
head(snps_id)
colnames(snps_geno) <- gsub("_.*", "", colnames(snps_geno))
snps_geno <- snps_geno[,-1]
snps_geno <- as.data.frame(snps_geno)

Random_SNP_FINAL <- data.frame(matrix(ncol = 8, nrow = 0))
colnames(Random_SNP_FINAL) <- colnames(snps_id)

# Our groups here
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

head(TFs)
# Empty vector for each iteration
props = c()
props_trans = c()
head(top_trans)
# Loop for each distance
for (d in distance) {
  
  # change to numeric to ensure binding
  TFs$Chromosome = as.numeric(TFs$Chromosome)
  top_trans$CHR = as.numeric(top_trans$CHR)
  # Number of trans SNPs (for calculating proportion)
  n_trans = 49
  print(n_trans)
  # Merged to TFs by chromosome
  # Calculate distance to TFs
  Master_trans <- dplyr::left_join(top_trans, TFs, by = c("CHR" = "Chromosome"))
  Master_trans$Distance_Start <- abs(as.numeric(Master_trans$SNP_POS) - as.numeric(Master_trans$POS))
  Master_trans$Distance_End <- abs(as.numeric(Master_trans$SNP_POS) - as.numeric(Master_trans$POS2))
  
  # Set up a column with distance of SNP to TFs
  Master_trans <- Master_trans %>% mutate(Distance = case_when(
    Master_trans$Distance_Start > Master_trans$Distance_End ~ Master_trans$Distance_End,
    Master_trans$Distance_End > Master_trans$Distance_Start ~ Master_trans$Distance_Start))
  
  # Get number of SNPs
  # Will use this to subset an equal number of SNPs (WITH REPLACEMENT) however many times
  length(unique(Master_trans$SNP))
  Master_trans <- Master_trans %>% filter(Distance <= d) # filter for genes which are close to SNPs based on distance
  Master_trans <- Master_trans[order(Master_trans$SNP),]
  Master_trans <- Master_trans %>% drop_na() # some genes do not have identifiers # doesnt matter as we are just pickng up TFs
  prop_trans <- length(unique(Master_trans$SNP)) / n_trans
  props_trans <- append(props_trans, prop_trans)

  # set seed for reproducibility
  # Set up an empty dataframe for each permutation
  set.seed(1285)
  Random_final = replicate(10000, sample_n(as.data.frame(snps_id), n_trans, replace = T), simplify = F)
    Random_final
    length(Random_final)
  head(snps_id)
  for (ele in 1:length(Random_final)) {
  
  # Reformat the Random SNP file
  Random_final[[ele]]$CHROM <- as.numeric(Random_final[[ele]]$CHROM)
  Random_final[[ele]]$POS <- as.numeric(Random_final[[ele]]$POS)
  TFs$Chromosome <- as.numeric(TFs$Chromosome)
  Random_SNPs_final <- Random_final[[ele]]
  
  # Join our random SNP set to TF by chromosome
  # Repeat process as above and calculate the distances
  Random_SNPs_final <- dplyr::left_join(Random_SNPs_final, TFs, by = c("CHROM"="Chromosome"))
  Random_SNPs_final$Distance_Start <- abs(as.numeric(Random_SNPs_final$POS.x) - as.numeric(Random_SNPs_final$POS.y))
  Random_SNPs_final$Distance_End <- abs(as.numeric(Random_SNPs_final$POS.x) - as.numeric(Random_SNPs_final$POS2))
  Random_SNPs_final <- Random_SNPs_final %>% mutate(Distance = case_when(
    Random_SNPs_final$Distance_Start > Random_SNPs_final$Distance_End ~ Random_SNPs_final$Distance_End,
    Random_SNPs_final$Distance_End > Random_SNPs_final$Distance_Start ~ Random_SNPs_final$Distance_Start))
  
  # Calculate the proportion of SNPs close to a TF based on a distance
  Random_SNPs_final <- Random_SNPs_final %>% filter(Distance <= d) %>% drop_na() # filter for genes which are close to SNPs
  prop_ran <- length(unique(Random_SNPs_final$ID)) / n_trans
  prop_ran
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
head(df)
# Convert to long format
df_long <- pivot_longer(df, cols = X1:X10000, names_to = "Iteration", values_to = "Prop_with_TFs")
head(df_long)
# Getting the statistics
trans_df <- data.frame(distance, props_trans)
trans_df$Condition = "Trans-eQTL"
df_long_simple = df_long %>% dplyr::select(1,3)
df_long_simple$Condition = "Random"
head(df_long_simple)

MASTER_STAT <- df_long_simple %>% group_by(distance) %>% dplyr::summarize(Mean_NULL = mean(Prop_with_TFs), SD_NULL = sd(Prop_with_TFs))
head(MASTER_STAT)
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



trans <-trans_df %>% dplyr::select(1,2)
trans$Condition = "Trans-eQTL"
colnames(trans)[2] <- "Ratio"
Random <- df_long_simple %>% filter(Condition == "Random") %>% dplyr::select(1,2)
Random$Condition <- "Random"
colnames(Random)[2] <- "Ratio"

data_plot <- rbind(trans,Random)
head(data_plot)
data_plot$distance <- data_plot$distance / 1000
distance = unique(data_plot$distance)
Plotting_simple <- data_plot %>% filter(distance >= 200) 


sum_plot <- data_plot %>% group_by(distance, Condition) %>% dplyr::summarize(Ratio = mean(Ratio))

table(data_plot$Condition)
# The Plot
ggplot(data = data_plot, aes(x = Ratio, y = distance, group = distance, colour = Condition)) +
  geom_point(size = 0.5, alpha = 1) +
  theme_bw() +
  labs(x = "Proportion",
       y = "Distance (000's bps)") +
  scale_color_manual(values = c("Trans-eQTL" = "#d7191c", "Random" = "#2c7bb6"), name = "") +
  scale_y_continuous(breaks= distance) +
  theme(axis.text.x = element_text(angle = 90)) + coord_flip() +
  geom_boxplot(data = Plotting_simple, outlier.color = NULL, outlier.alpha = 0, scale = 0.9, alpha = 0.2, fill = "#2c7bb6") + coord_flip() +
  geom_line(data = sum_plot[sum_plot$Condition == "Trans-eQTL",], col = "#d7191c", size = 0.5, alpha = 0.8, group = 1) +
  geom_line(data = sum_plot[sum_plot$Condition == "Random",], col = "#2c7bb6", size = 0.5, alpha = 0.8, group = 1) +
  xlim(0,1) +
  theme(axis.text.x = element_text(angle = 90, size = 12.5, colour = "black"),
        axis.text.y = element_text(angle = 0, size = 18, colour = "black"),
        legend.text = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 18, colour = "black"),
        axis.title.y = element_text(size = 18, color = "black"),
        legend.key.size = unit(3, "lines"))  # adjust legend key size)
