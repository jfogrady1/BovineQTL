##########################################################################################
##########################################################################################
################## Script to identify if Trans-eQTLs are close to TFs/TF cofactors########
########################################################################################## 
##########################################################################################

library(tidyverse)
library(ggplot2)
library(tidyr)
library(ggdist)
library(gghalves)
library(tidyquant)
library(qqman)
library(ggsignif)
library(data.table)

args = commandArgs(trailingOnly=TRUE)
###################################################################################################################################################################
###################################################################################################################################################################
#################################### Get SNPs on different chromosomes and get their allele frequencies ###########################################################
###################################################################################################################################################################
###################################################################################################################################################################

# Get data - read in our trans-eQTLs which are deemed significant
# ALL TransFDR
# Args[1] ALL_TransFDR5.txt
trans_data <- read.csv(args[1], sep = "\t") %>% filter(FDR < 0.01)


# consider trans-SNPs on different chromosomes only
# need gene location information
# need snp location information also
# args[2] gene_loc (expressed genes)
gene_loc <- read.csv(args[2] ,sep = "\t") %>% dplyr::select(1,2)
colnames(gene_loc)[1] <- "gene"
colnames(gene_loc)[2] <- "gene_CHR"
trans_data <- left_join(trans_data, gene_loc) # join everything


# snp location
snp_loc <- read.csv(args[3] , sep="\t") %>% dplyr::select(1,2)
colnames(snp_loc)[1] <- "snps"
colnames(snp_loc)[2] <- "SNP_CHR"
trans_data <- left_join(trans_data, snp_loc)


# Filter for SNPs which are on a different chromsome
trans_data <- trans_data %>% filter(gene_CHR != SNP_CHR) # find SNPs which aren't on same CHR as gene

# How many trans SNPs on different chromosome do we have?
dim(trans_data) # This is where it is different to circos plot, here we have ~5857 SNPs which we need to sample

write.table(trans_data, file = "Trans_SNPS_diff_Chr_FDR1.txt", sep = ",", row.names = F, col.names = T) # take this file to statistical script
# FINAL FILTERING OF eQTLS
# get SNP id and reformat the column names
trans_data <- trans_data %>% separate(., col = "snps", into = c("CHR", "POS", "REF", "ALT"), sep = ":") 
trans_data <- trans_data %>% unite(., col = "TRANS_POS", "CHR", "POS", sep = ":")
colnames(trans_data)[4] <- "TRANS_GENE"
colnames(trans_data)[1] <- "FINAL_POS1" # trans positions

# get allele frequencies of the trans SNPs - note here we have the AFs for trans-SNPs with a cis-EQTL also
# args[4] Trans.SNPs.AF.header.fixed.frq
AF = read.csv(args[4], sep = "\t", header = F)
colnames(AF) = c("CHROM", "POS_Trans", "NALLELES", "NCHR", "REF:FREQ", "ALT:FREQ")
AF <- AF[-1,] # remove first row

dim(AF)
# Reformat the trans data so it can be jointed with AF data
trans_data = separate(trans_data, FINAL_POS1, into = c("CHR_trans", "POS_trans"), sep = ":") 
trans_data$CHR_trans <- as.character(trans_data$CHR_trans)
trans_data$POS_trans <- as.character(trans_data$POS_trans)
AF$POS_Trans <- as.character(AF$POS_Trans)
trans_data$CHR_trans <- as.numeric(trans_data$CHR_trans)


# Join the trans SNPs (which have a cis-eQTL to the same gene) with their Allele frequencies 
# Format the last 2 columns to extend out
# Select columns we want
trans_data = left_join(trans_data, AF, by = c("CHR_trans" = "CHROM",
                                              "POS_trans" = "POS_Trans"))
trans_data = separate(trans_data, `REF:FREQ`, into = c("ref", "AF_ref"), sep = ":")
trans_data = separate(trans_data, `ALT:FREQ`, into = c("alt", "AF_alt"), sep = ":")
trans_snps_no_cis <- trans_data %>% filter(., is.na(trans_data$AF_ref)) %>% dplyr::select(1,2,3,4) ## now have other trans SNPs that are on different chromosomes but not associated with a cis-eQTL of the same gene ##

# Now need to get the allele freqeuncy for other putative trans SNPs
# These are the trans_snps_no_cis data frame
# args[5] "NonTrans.SNPs.AF.header.fixed.frq"
AF2 <- read.csv(args[5], sep = "\t", header = F)
colnames(AF2) = c("CHROM", "POS_Trans", "NALLELES", "NCHR", "REF:FREQ", "ALT:FREQ")
AF2$POS_Trans <- as.character(AF2$POS_Trans)
AF2$CHROM = as.numeric(AF2$CHROM)
trans_snps_no_cis <- left_join(trans_snps_no_cis, AF2, by = c("CHR_trans" = "CHROM",
                                                              "POS_trans" = "POS_Trans"))
trans_snps_no_cis = separate(trans_snps_no_cis, `REF:FREQ`, into = c("ref", "AF_ref"), sep = ":")
trans_snps_no_cis = separate(trans_snps_no_cis, `ALT:FREQ`, into = c("alt", "AF_alt"), sep = ":")

# Merged trans snps and those trans snps which are in LD with cis -eQTLs
trans_data_final <- trans_data %>% filter(.,!is.na(NCHR)) # filter trans data for SNPs which are in LD with other cis-eTL snps
trans_data_final <- trans_data_final %>% dplyr::select(1,2,3,4,12,13,14,15,16,17)
trans_data_final <- rbind(trans_data_final, trans_snps_no_cis) # now have allele frequencies for all trans-SNPs
trans_data_final$SNP <- paste0(trans_data_final$CHR_trans, ":", trans_data_final$POS_trans)
dim(trans_data_final) # same as original data frame

# Now have AFs for all trans-SNPs


###############################################################################################
###############################################################################################
########### Removing trans-eQTLs in high LD with cis eQTLs ####################################
###############################################################################################
###############################################################################################


# Filtering for SNPs in LD with cis-eQTLs. Will use SNPS with LD > 0.01
# args[6] LDR2FULL.interchrom.geno.ld
LD_trans <- read.table(args[6], header = T) # interchromosomal LD  between trans and cis-eQTLs from VCF tools
LD_trans <- LD_trans %>% unite(., FINAL_POS1, "CHR1", "POS1", sep = ":" ) # get positions ready for filtering
LD_trans <- LD_trans %>% unite(., FINAL_POS2, "CHR2", "POS2", sep = ":" )
LD_trans_high_LD <- LD_trans %>% group_by(FINAL_POS1) %>% filter(R.2 > 0.01) # filter for SNPs above cut off (these SNPs have "high LD")
LD_trans_high_LD <- LD_trans_high_LD[order(LD_trans_high_LD$R.2),] # order the SNPs by high LD
LD_trans_high_LD <- LD_trans_high_LD[!duplicated(LD_trans_high_LD$FINAL_POS1),] # get unique SNPs in this list which have high LD


# Now need to remove these SNPs from our trans file
# These are SNPs which are on different chromosomes but which are not in LD with top cis-eQTLs, known as "high-LD" SNPs
# Will extract genes close to these SNPS that are remaining
trans_data_filtered <- trans_data_final %>% filter(!(trans_data_final$SNP %in% LD_trans_high_LD$FINAL_POS1)) # keep SNPs not in "high" LD - also keep SNPs that we do not have lD info for
dim(trans_data_filtered)
trans_data_filtered$snp_pos <- trans_data_filtered$SNP
trans_data_filtered <- trans_data_filtered %>% unite(., col = "snp", "SNP", "REF", "ALT", sep = ":")
trans_data_filtered <- trans_data_filtered %>% dplyr::select(1,2,3,6,7,8,9)

# Now have SNPs which are on different chromosomes and which have been filtered for LD and which have an FDR of 0.01
# Get the gene information etc for the remaining snps
trans_data <- read.csv(args[1], sep = "\t")
trans_data_filtered <- left_join(trans_data_filtered, trans_data, by = c("snp" = "snps")) # note, increase in rows because some SNPs are associated with multiple genes
trans_data_filtered$POS_trans <- as.numeric(trans_data_filtered$POS_trans)
trans_data_filtered = trans_data_filtered[order(trans_data_filtered$CHR_trans,trans_data_filtered$POS_trans),]
snp_final <- trans_data_filtered %>% dplyr::select(snp)
head(trans_data_filtered, 20)
write.table(trans_data_filtered, file = "Trans_SNPS_FDR1_nohighR2.txt", sep = ",", row.names = F, col.names = T) # take this file to statistical script
write.table(snp_final, file = "Final_Trans_SNPs_&_Genes_4_LD_Calc.txt", sep = "\t", col.names = F, row.names = F, quote = F)


###################################################################################################################################################################
###################################################################################################################################################################
########################################### Extracting middle SNP in locus of trans-eQTLs #########################################################################
###################################################################################################################################################################
###################################################################################################################################################################
trans_data_filtered <- trans_data_filtered %>% group_by(CHR_trans, gene) %>% mutate(distance_to_previous_SNP = c(0, diff(POS_trans))) %>% as.data.frame()
trans_data_filtered <- trans_data_filtered %>% group_by(CHR_trans, gene) %>% mutate(distance_to_next_SNP = c(diff(POS_trans), 0)) %>% as.data.frame()


trans_data_filtered <- trans_data_filtered %>% mutate(locus_1 = case_when(trans_data_filtered$distance_to_previous_SNP <= 1e6 & trans_data_filtered$distance_to_previous_SNP > 0 ~ "locus",
                                                                          trans_data_filtered$distance_to_previous_SNP == 0 ~ "single",
                                                                          trans_data_filtered$distance_to_previous_SNP > 1e6  ~ "distant"))
trans_data_filtered <- trans_data_filtered %>% mutate(locus_2 = case_when(trans_data_filtered$distance_to_next_SNP <= 1e6 & trans_data_filtered$distance_to_next_SNP > 0 ~ "locus",
                                                                        trans_data_filtered$distance_to_next_SNP == 0 ~ "single",
                                                                        trans_data_filtered$distance_to_next_SNP > 1e6  ~ "distant"
                                                                        ))



# locus 1 = previous SNP
trans_data_filtered <- trans_data_filtered %>% mutate(locus_3 = case_when(trans_data_filtered$locus_1 == "locus" & trans_data_filtered$locus_2 == "single" ~ "Cluster",
                                                                          trans_data_filtered$locus_1 == "locus" & trans_data_filtered$locus_2 == "distant" ~ "Cluster",
                                                                          trans_data_filtered$locus_1 == "distant" & trans_data_filtered$locus_2 == "locus" ~ "Cluster",
                                                                          trans_data_filtered$locus_1 == "single" & trans_data_filtered$locus_2 == "locus" ~ "Cluster",
                                                                          trans_data_filtered$locus_1 == "locus" & trans_data_filtered$locus_2 == "locus" ~ "Cluster",
                                                                          trans_data_filtered$locus_1 == "single" & trans_data_filtered$locus_2 == "single" ~ "single",
                                                                          trans_data_filtered$locus_1 == "distant" & trans_data_filtered$locus_2 == "single" ~ "distant",
                                                                          trans_data_filtered$locus_1 == "single" & trans_data_filtered$locus_2 == "distant" ~ "distant"
                                                                          ))
snps_single <- trans_data_filtered %>% filter(trans_data_filtered$locus_3 == "single")
distant_snps <- trans_data_filtered %>% filter(trans_data_filtered$locus_3 == "distant")
trans_data_filtered <- trans_data_filtered %>% filter(trans_data_filtered$locus_3 == "Cluster")

# Need to extract middle SNPs in a locus
# group by chromosome and gene (locus associated with the gene), count how many SNPs are there and then find the one in the middle
trans_data_extract <- trans_data_filtered %>% group_by(CHR_trans, gene) %>% dplyr::summarise(total_count=n(), extract_SNP = ceiling(n()/2), .groups = 'drop') #%>% order(gene) #%>% select(extract_SNP)


# Split data up into a list, we will cycle through this and at each interval, extract the middle SNP we want
#trans_data_filtered <- trans_data_filtered[order(trans_data_filtered$gene),]
trans_list <- trans_data_filtered %>% group_split(CHR_trans, gene)
middle <- as.vector(as.numeric(trans_data_extract$extract_SNP))
final_trans <- data.frame(matrix(ncol = 17, nrow = 0))
length(trans_list)
colnames(final_trans) = colnames(trans_data_filtered)



for (loc in 1:length(trans_list)){
  gene <- trans_list[[loc]]
  row = gene[middle[loc],]
  final_trans <- rbind(final_trans, row)
}

dim(final_trans)
final_trans <- rbind(final_trans, snps_single)
final_trans <- rbind(final_trans, distant_snps)
dim(final_trans)

## At this stage, have extracted SNPs in each locus
write.table(final_trans, file = "Final_trans_FDR_01.txt", sep = ",", row.names = F, col.names = T) # take this file to statistical script
