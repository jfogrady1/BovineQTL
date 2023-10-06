# Script to calculate intrachromosomal LD between trans-eQTLs and top cis-eQTLs of each gene

########################
# Load in the libraries
########################

# Load in library
library(SNPRelate)
library(data.table)
library(tidyverse)
library(devtools)
library("vcfR")
library(ggplot2)
library(ggridges)
args = commandArgs(trailingOnly = T)

#args[1] =  "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/ALL_IMPUTED_UPDATED.vcf.gz"
#args[2]  = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/ALL.cis_qtl_fdr0.05.txt"
#args[3] = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/ALL_trans_FDR_corrected.txt"
#args[4] = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/intra_chrom_trans_SNPs.txt"
#args[5] = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/transQTL/trans_intra.ld.gz"
#args[6] = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/transQTL/intra_transVar_LD_between_them.pdf"
#args[7] = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/intra_observed_cis_trans_SNPs.txt"
#args[8] = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/transQTL/trans_cis_intra_LD.geno.ld"
#args[9] = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/intra_expected_cis_trans_SNPs.txt"
#args[10] = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/transQTL/trans_cis_intra_LD_null.geno.ld"
#args[11] = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/transQTL/intra_trans_V_15_permuted.pdf"
#args[12] = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/transQTL/intra_trans_cis_permutation_results.txt"
#args[13]= "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/transQTL/intra_trans_cis_permutation_raw.txt"
#args[14] = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/ALL_intrachromosomal_LD_corrected.txt"
vcf <- args[1]
vcf_data <- read.vcfR(vcf)


# TOP cis-eQTLs
cis = fread(args[2])

# Trans-eQTL results
trans = fread(args[3], header = F) %>% filter(V10 < 0.05)
trans$intra = if_else(trans$V2 == trans$V5, TRUE, FALSE)
trans <- trans %>% filter(intra == TRUE)


# Now have trans SNPs which are on same chromosome to associated gene
length(unique(trans$V1))
# 26 eGenes were associated on same chromosome as trans-eQTLs
dim(trans) # 26 eGenes associated with 585 trans-evariants





################################################################
# Calcualte LD among trans-variants associated with same gene###
################################################################



# Now calculate LD matrix for individual SNPs associated with each of the 26 genes 

# Estimate LD between/among trans-eQTLs which are associated to gene on same chromsoome
write.table(unique(trans$V4), args[4], col.names = F, row.names = F, quote = F)
system(paste0("plink --cow --bfile /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/ALL_genotypes_tensor --keep-allele-order --r2 gz --ld-snp-list /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/intra_chrom_trans_SNPs.txt --out /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/transQTL/trans_intra"))



# Now bind this information to current information of trans-eQTLs
LD_trans = fread(args[5])
LD_trans = LD_trans %>% select(3,6,7)
intra_results = list()
head(LD_trans)
genes = unique(trans$V1)
genes
for (g in genes) {
    trans_gene = trans %>% filter(V1 == g) %>% left_join(., LD_trans, by = c("V4" = "SNP_A"))
    print(trans_gene)
    intra_results[[g]] <- trans_gene
}
intra_results
# Generate a dummy plot to visualise this for the 26 genes
plots_df = data.frame(matrix(ncol=2, nrow = 0))
colnames(plots_df) = c("phenotype", "r")
test <- intra_results$ENSBTAG00000022329 %>% as.data.frame()
head(test)
colnames(test)
for(ele in genes) {
    test_df = intra_results[[ele]] %>% as.data.frame()
    test_df <- test_df %>% filter(V4 != SNP_B)
    test_df <- test_df %>% select(V1, R2)
    test_df$R <- sqrt(test_df$R2)
    colnames(test_df) <- c("phenotype", "r")
    print(head(test_df))
    plots_df <- rbind(plots_df, test_df)
}

mean(plots_df$r)
sd(plots_df$r)
ggplot(data = plots_df, aes(x = phenotype, y = r, fill = phenotype)) + geom_boxplot() + coord_flip() +
  theme_bw() +
  labs(y = "Linkage disequilibrium (r) among trans-eVariants", x = "Intra-chromosomal trans-eGene") +
  theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black"),
        axis.text.y = element_text(angle = 0, size = 15, colour = "black"),
        axis.title.x = element_text(size = 18, colour = "black"),
        axis.title.y = element_text(size = 18, color = "black"),
        legend.position = "none") +
        scale_y_continuous(breaks = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))
ggsave(args[6], width = 12, height = 12, dpi = 600)
# Given high LD between SNPs associated with the same gene on the same chromosome in trans,
# We will just variant with the smallest p-value
final_intra_trans = data.frame(matrix(ncol = 13, nrow = 0))
colnames(final_intra_trans) <- colnames(intra_results$ENSBTAG00000022329)
for (g in genes) {
    trans_gene = trans %>% filter(V1 == g) %>% filter(V7 == min(V7))
    # remove SNPs if tie occurs
    trans_gene <- trans_gene[!duplicated(trans_gene$V7),]
    final_intra_trans = rbind(final_intra_trans, trans_gene)
}
dim(final_intra_trans) # 26 trans-evariants associated with 26 trans-eGenes


#################################################################################
# Calcualte LD among trans-variants and cis-variants associated with same gene###
#################################################################################


# 26 varaints and 26 eGenes
dim(final_intra_trans)
length(unique(final_intra_trans$V1))
# Now pull out top cis-eQTL from the cis eQTLs

# 25 cis-eGenes meaning one cis-eGene did not have any varaint in cis to be tested
cis <- cis %>% filter(phenotype_id %in% genes)

#24 cis-eGenes
cis <- cis %>% filter(is_eGene == TRUE)
ALL <- left_join(cis, final_intra_trans, by = c("phenotype_id" = "V1"))

# Get position of cis-eGene variant, the most significant one
# First get variants
cis_pos <- vcf_data@fix %>% as.data.frame() %>% select(3,2)
colnames(cis_pos)[2] <- "cis_pos"
ALL <- left_join(ALL, cis_pos, by = c("variant_id" = "ID"))


# V6 and cis_pos represeent the #trans_position and cis_position
# now calculate the distance between
ALL$distance = abs(as.numeric(ALL$V6) - as.numeric(ALL$cis_pos))
ALL <- ALL %>% select(1,7,32,24, 25, 26,33)
colnames(ALL) <- c("gene", "cis_variant", "cis_pos", "trans_variant", "CHR", "trans_pos", "distance")
observed <- ALL %>% select(cis_variant, trans_variant)
observed <- c(observed$cis_variant, observed$trans_variant)

# write to file and compute the LD between these varaints
write.table(observed, args[7], row.names = F, col.names = F, quote = F, sep = "\t")
system(paste0("vcftools --gzvcf /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/ALL_IMPUTED_UPDATED.vcf.gz --snps /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/intra_observed_cis_trans_SNPs.txt --geno-r2 --out /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/transQTL/trans_cis_intra_LD"))

# now get actual pairs
# Will use this for filtering the LD file
ALL$cis_var <- paste0(ALL$CHR, ":", ALL$cis_pos)
ALL$trans_var <- paste0(ALL$CHR, ":", ALL$trans_pos)
ALL$pairs = paste0(ALL$cis_var, "-", ALL$trans_var)
LD_trans = fread(args[8])
summary(LD_trans$R^2)

LD_trans$cis_pos = paste0(LD_trans$CHR, ":", LD_trans$POS1)
LD_trans$trans_pos = paste0(LD_trans$CHR, ":", LD_trans$POS2)
LD_trans$pairs <- paste0(LD_trans$cis_pos, "-", LD_trans$trans_pos)
LD_trans$pairs_reverse <- paste0(LD_trans$trans_pos, "-", LD_trans$cis_pos)


LD_trans <- LD_trans %>% filter(LD_trans$pairs %in% ALL$pairs | LD_trans$pairs_reverse %in% ALL$pairs)
LD_trans$r <- sqrt(LD_trans$`R^2`) #  = absolute correlation



sd(LD_trans$r)
# Now have LD between trans-eQTLs and cis-eQTLs of same gene on same chromosome



#############################################################
# Calcualte LD among null trans-variants and cis-variants ###
#############################################################


# Now need to sample 10000 SNPs which are > 5Mb but < +2SD away from mean of our set

# Set up the variables
min_distance = 5e6
max_distance = 5e6 + (2 * sd(ALL$distance))
max_distance
min_distance

set.seed(18794)
# now sample a 10 million SNPs with replacement to get a rough number of SNPs
vcf_master = vcf_data@fix %>% as.data.frame() %>% select(1,2,3)

# Sample the SNPs 500000
vcf_1 = vcf_master[sample(nrow(vcf_master), 3e7, replace = T),]
vcf_2 = vcf_master[sample(nrow(vcf_master), 3e7, replace = T),]

dim(vcf_1)
colnames(vcf_1) <- c("CHR1", "POS1", "VAR1")
colnames(vcf_2) <- c("CHR2", "POS2", "VAR2")


merged <- cbind(vcf_1, vcf_2)
rm(vcf_1, vcf_2)


# now order so we can calculate distance
merged <- merged[order(merged$CHR1, merged$CHR2),]

# Filter for same chromosome
merged <- merged %>% filter(CHR1 == CHR2)

# Calcualte distance
# filter for those which meet our criteria
merged$distance = abs(as.numeric(merged$POS1) - as.numeric(merged$POS2))
merged <- merged %>% filter(distance > min_distance & distance < max_distance)
head(merged)
dim(merged)

# Now sample 50,000
merged <- merged[sample(nrow(merged), 50000, replace = F),]

expected <- c(merged$VAR1, merged$VAR2)
head(expected)
write.table(expected, args[9], row.names = F, col.names = F, quote = F, sep = "\t")

if (!file.exists(args[10])) {
    system(paste0("vcftools --gzvcf /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/ALL_IMPUTED_UPDATED.vcf.gz --snps /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/intra_expected_cis_trans_SNPs.txt --geno-r2 --ld-window-bp ", max_distance," --ld-window-bp-min ",99999, " --out /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/transQTL/trans_cis_intra_LD_null"))
}

# here have set of samples which are

#############################################################
# Permute LD and calculate difference to observed 50000 #####
#############################################################


# Read in the file
# Calculate distance again
# Randomly sample 50,000 SNPs
data <- fread(args[10])
data$distance <- abs(as.numeric(data$POS1) - as.numeric(data$POS2))
data <- data %>% filter(distance > min_distance & distance < max_distance)
data$r <- sqrt(as.numeric(data$`R^2`))

data_s <- sample_n(data, 1000000)
dim(data_s)
# Now randmomly permute the median of 23 sets 100000 times
stats_df <- data.frame(matrix(ncol=2, nrow=0))
stats_temp <- data.frame(matrix(ncol=2, nrow=0))
colnames(stats_df) <- c("r", "permutation")
colnames(stats_temp) <- c("r", "permutation")
for (i in 1:10000) {
    set.seed(i)
    data_temp <- data_s[sample(nrow(data_s), 24, replace = T),]
    data_temp <- data_temp %>% select(r) 
    data_temp$permutation <- i
    stats_df <- rbind(stats_df, data_temp)

}

head(stats_df)
# Now onto permutation, what proportion of tests are not significantly different under the null
permutation_results = list()

for(i in 1:10000) {
    print(i)
    curr_round <- stats_df %>% filter(permutation == i)
    result <- wilcox.test(LD_trans$r, curr_round$r, alternative = "greater", exact = F)
    permutation_results[i] <- result$p.value

}
stats_df_p <- data.frame(matrix(ncol=2, nrow=0))
stats_temp <- data.frame(matrix(ncol=2, nrow=0))
colnames(stats_df_p) <- c("p", "permutation")
colnames(stats_temp) <- c("p", "permutation")
head(stats_temp)

for(i in 1:10000) {
    stats_temp[1,] <- c(permutation_results[i], i)
    stats_df_p <- rbind(stats_df_p, stats_temp)

}


#############################################################
# Get statistics and perform multiple testing correction ####
#############################################################
head(stats_df_p)
stats_df_p$padj <- p.adjust(as.numeric(stats_df_p$p), method = "BH")

# Plotting some distributions
num = sample(seq(1:10000), 20)
data_plot <- stats_df %>% filter(stats_df$permutation %in% num)
min_permute <- stats_df_p %>% filter(stats_df_p$padj == max(stats_df_p$padj)) %>% select(permutation) %>% as.vector()
min_plot <- stats_df %>% filter(stats_df$permutation %in% min_permute)
data_plot <- rbind(data_plot, min_plot)

data_plot$Type = "null"
data_plot$Type = "null"
LD_plot <- LD_trans %>% select(r)
LD_plot$permutation <- "observed"
LD_plot$Type = "observed"
LD_plot <- rbind(data_plot, LD_plot)


stats_df_p$signif = if_else(stats_df_p$padj < 0.05, TRUE, FALSE)
table(stats_df_p$signif)

0/100000  # permuted p-value adjustement 0.00161

ggplot(data = LD_plot, aes(x = permutation, y = r, fill = Type)) + geom_boxplot() +
  theme_bw() + coord_flip() +
  labs(y = "Linkage disequilibrium (r) among trans-eVariants", x = "Permutation set") +
  theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black"),
        axis.text.y = element_text(angle = 0, size = 15, colour = "black"),
        axis.title.x = element_text(size = 18, colour = "black"),
        axis.title.y = element_text(size = 18, color = "black"),
        legend.key.size = unit(1.2, 'cm'), #change legend key size
        legend.key.height = unit(1.2, 'cm'), #change legend key height
        legend.key.width = unit(1.2, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=12))

ggsave(args[11], width = 12, height = 12, dpi = 600)


write.table(stats_df_p, file = args[12], row.names = F, col.names = T, quote = F, sep ="\t")
write.table(stats_df, file = args[13], row.names = F, col.names = T, quote = F, sep ="\t")
dim(stats_df)

#############################################################
# Filter for eGenes with top SNP not in LD with cis-eQTL ####
#############################################################


# Final step is to filter for eGenes which were not in high LD with cis-eQTLs of the same gene

# Genes not associated in cis
Final_inter_trans_genes_no_cis <- final_intra_trans[!(V1 %in% ALL$gene),]

# Genes where top SNP is not in LD with cis-eQTL of same gene
trans = fread(args[3], header = F) %>% filter(V10 < 0.05)
trans$intra = if_else(trans$V2 == trans$V5, TRUE, FALSE)
trans <- trans %>% filter(intra == TRUE)
trans$trans_pos <- paste0(trans$V5, ":", trans$V6)

LD_trans = fread(args[8])
LD_trans$cis_pos = paste0(LD_trans$CHR, ":", LD_trans$POS1)
LD_trans$trans_pos = paste0(LD_trans$CHR, ":", LD_trans$POS2)
LD_trans$pairs <- paste0(LD_trans$cis_pos, "-", LD_trans$trans_pos)
LD_trans$pairs_reverse <- paste0(LD_trans$trans_pos, "-", LD_trans$cis_pos)
LD_trans <- LD_trans %>% filter(LD_trans$pairs %in% ALL$pairs | LD_trans$pairs_reverse %in% ALL$pairs)
LD_trans <- LD_trans[(as.numeric(LD_trans$`R^2`) < 0.01),]

# Perform filtering
trans <- trans %>% filter(trans_pos %in% LD_trans$cis_pos) # cis as it is the reverse pair
trans <- trans %>% select(-trans_pos)
trans_filtered <- rbind(Final_inter_trans_genes_no_cis, trans)

# Now collect the variants and filter
trans = fread(args[3], header = F) %>% filter(V10 < 0.05)
trans <- trans %>% filter(trans$V1 %in% trans_filtered$V1)
trans <- trans %>% filter(V10 < 0.01)

write.table(trans, args[14], col.names = F, row.names = F, quote = F, sep = "\t")
