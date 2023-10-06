###################################################################################################
###################################################################################################
################# This script does 2 things                                                   #####
################# Prepares trans and cis SNPs for caluclaiton of LD                           #####
################# Prepares random set of 'trans' SNPS (same chr and AF) for calcualtion of LD #####
#################                                                                             #####
###################################################################################################
###################################################################################################



## Packages
library(tidyverse)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggdist)
library(gghalves)
library(tidyquant)

# load in data
# Arg[1] = cis data
# Arg[2] = trans data (both significant) = output from eQTL analysis (order and select both)
cis_data <- read.csv("../../results/EQTL/RAW/ALL_CisFDR5.txt", sep = "\t")
trans_data <- read.csv("../../results/EQTL/RAW/ALL_TransFDR5.txt", sep = "\t") # all data with most points - used also in TWAS so need to plot that




# ONLY consider trans-SNPs on different chromosomes
# need gene location information
# need snp location information also
# Arg[3] = gene location file used in normalisation procedure
gene_loc <- read.csv("../../results/RNA-seq/Normalisation/gene_location.txt", sep = "\t") %>% select(1,2) # get gene location
colnames(gene_loc)[1] <- "gene"
colnames(gene_loc)[2] <- "gene_CHR"
trans_data <- left_join(trans_data, gene_loc) # join everything


# snp location
# Arg[4] = SNP location used in EQTL analysis
snp_loc <- read.csv("../../results/SNP_data/Imputation_Performance/eQTL_SNP_POS.txt", sep="\t") %>% select(1,2)
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

write.csv(snp_list, "TransCisLDComparisonSNPS.txt", quote = F, row.names = F, col.names = F)
write.csv(trans_data$snps, "Trans_LD_comparisonSNPs.txt", quote = F, row.names = F, col.names = F) # will extract from column 1
write.csv(cis_data$snps, "Cis_LD_comparisonSNPs.txt", quote = F, row.names = F, col.names = F) # will extract from column 2 of outputted file

dim(snp_list)
# Subset these trans SNPs and calculate interchomosomal LD
# firt output will be subsetted (i.e. interchromosomal SNPs between both cis and trans SNPs)
system(paste("vcftools --gzvcf ../../results/SNP_data/Imputation_Performance/SNP_DATA_All.IMPUTED.FILTERED.RENAMED.vcf.gz --snps TransCisLDComparisonSNPS.txt --recode --recode-INFO-all --stdout | bgzip -c > TransCisLDSNPs.vcf.gz && tabix -p vcf TransCisLDSNPs.vcf.gz"))
system(paste("vcftools --gzvcf TransCisLDSNPs.vcf.gz  --interchrom-geno-r2 --out LDR2FULL")) # LD comparison with random SNPs

###################################################################################################################################################################################
###################################################################################################################################################################################
################## Generation of Random trans SNPs and looking at LD between them and cis SNPs ####################################################################################
###################################################################################################################################################################################
###################################################################################################################################################################################


trans_data <- trans_data %>% select(1,7,8)
dim(trans_data) # need to generate 2627 (2616 unique)
length(unique(trans_data$SNP_CHR)) # from all 29 autosomes.

# 1. need to get allele frequency for current SNPs
system(paste("vcftools --gzvcf TransCisLDSNPs.vcf.gz --snps Trans_LD_comparisonSNPs.txt --freq --out Trans.SNPs.AF"))
system(paste("tail -n+2 Trans.SNPs.AF.frq  > Trans.SNPs.AF.header.fixed.frq"))

# Read in the file produced by vcf tools
# Change the name of the columns
AF = read.csv("LD_work/Trans.SNPs.AF.header.fixed.frq", sep = "\t", header = F)
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

system(paste("vcftools --vcf ../../../data/SNP_data/imputation/Final_data_for_eqtl/Axiom.IMPUTED.FILTERED.simple.codes.EQTL.vcf --exclude Trans_LD_comparisonSNPs.txt --freq --out NonTRANS.SNPs.AF"))
system(paste("tail -n+2 NonTRANS.SNPs.AF.frq  > NonTrans.SNPs.AF.header.fixed.frq")) # format the header


# Format the allele frequency file as above

AF = read.csv("LD_work/NonTrans.SNPs.AF.header.fixed.frq", sep = "\t", header = F)
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

# set seed
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

# Prepare SNPs for vcftools
# Calcualte LD between Random SNPs and cis eqtls
system(paste("vcftools --vcf ../../../data/SNP_data/imputation/Final_data_for_eqtl/Axiom.IMPUTED.FILTERED.simple.codes.EQTL.vcf --snps RandomTransCisLDComparisonSNPS.txt --recode --recode-INFO-all --stdout | bgzip -c > RandomTransCisLDSNPs.vcf.gz && tabix -p vcf RandomTransCisLDSNPs.vcf.gz"))
system(paste("vcftools --gzvcf RandomTransCisLDSNPs.vcf.gz  --interchrom-geno-r2 --out LDR2RANDOMFULL"))


# Remove variables taking up space and which we no longer require
rm(snp_loc)
rm(potent)
rm(gene_loc)
rm(AF)
rm(LD)
rm(Random_SNPs_final)
rm(snp_list)
rm(AF_long)



#########################################################################################################
#########################################################################################################
#################################### Comparison of Groups ###############################################
#########################################################################################################
#########################################################################################################

# Getting LD between cis and trans SNPs
# read in all data
trans_snps = read.table("LD_work/Trans_LD_comparisonSNPs.txt", header = F)
ran_snps = read.table("LD_work/RandomTrans_LD_comparisonSNPs.txt", header = T, sep = ",") %>% select(1)
cis_snps = read.table("LD_work/Cis_LD_comparisonSNPs.txt", header = F)



# Format the DATA
trans_snps = trans_snps %>% separate(., V1, into = c("CHR", "POS", "REF", "ALT"), sep = ":")
trans_snps = trans_snps %>% unite(., col = "FINAL_POS1", "CHR", "POS", sep =":")
ran_snps = ran_snps %>% separate(., snps, into = c("CHR", "POS", "REF", "ALT"), sep = ":")
ran_snps = ran_snps %>% unite(., col = "FINAL_POS1", "CHR", "POS", sep =":")
cis_snps = cis_snps %>% separate(., V1, into = c("CHR", "POS", "REF", "ALT"), sep = ":")
cis_snps = cis_snps %>% unite(., col = "FINAL_POS2", "CHR", "POS", sep =":")


# Do it for our LD SNPs
LD_trans = read.table("LD_work/LDR2FULL.interchrom.geno.ld", header = T) # interchromosomal
LD_trans <- LD_trans %>% unite(., FINAL_POS1, "CHR1", "POS1", sep = ":" ) # get positions ready for filtering
LD_trans <- LD_trans %>% unite(., FINAL_POS2, "CHR2", "POS2", sep = ":" )
LD_trans_clean <- LD_trans[LD_trans$FINAL_POS1 %in% trans_snps$FINAL_POS1,]
LD_trans_clean <- LD_trans_clean[LD_trans_clean$FINAL_POS2 %in% cis_snps$FINAL_POS2,]
rm(LD_trans)


# DO it for the Random SNPS
LD_random = read.table("LD_work/LDR2RANDOMFULL.interchrom.geno.ld", header = T) # random interchromosomal
LD_random <- LD_random %>% unite(., FINAL_POS1, "CHR1", "POS1", sep = ":" ) # get positions ready for filtering
LD_random <- LD_random %>% unite(., FINAL_POS2, "CHR2", "POS2", sep = ":" )
LD_random_clean <- LD_random[LD_random$FINAL_POS1 %in% ran_snps$FINAL_POS1,]
LD_random_clean <- LD_random_clean[LD_random_clean$FINAL_POS2 %in% cis_snps$FINAL_POS2,]
rm(LD_random)


# Formatting
LD_trans_final <- left_join(LD_trans_clean, trans_snps, by = c("FINAL_POS1"))
LD_random_final <- left_join(LD_random_clean, ran_snps, by = c("FINAL_POS1"))
rm(LD_trans_clean)
rm(LD_random_clean)

LD_trans_final <- LD_trans_final %>% select(1,2,4)
LD_random_final <- LD_random_final %>% select(1,2,4)
LD_trans_final$Condition = "Trans-eQTL"
LD_random_final$Condition = "Random"
colnames(LD_trans_final)[3] = "R2"
colnames(LD_random_final)[1] = "FINAL_POS1"
colnames(LD_random_final)[3] = "R2"


#########################################################################################
##################### Filtering both for cis and trans SNPs of SAME GENE ################
#########################################################################################


# Need to get back in Trans genes and Cis genes
# Filter Trans SNPs with a cis-eQTL for the same genes
# Group by Trans SNP then select max R2 value for comparison.
# Use these SNPs to compare V LD of Random
cis_data <- read.csv("Results/AllCisFDR5.csv", sep = ",")
trans_data <- read.csv("Results/AllTransFDR5.csv", sep = ",")


# consider trans-SNPs on different chromosomes
# need gene location information
# need snp location information also
gene_loc <- read.csv("Inputs/GeneLocSymbolFullLength.txt") %>% select(1,2) # get gene location
colnames(gene_loc)[1] <- "gene"
colnames(gene_loc)[2] <- "gene_CHR"
trans_data <- left_join(trans_data, gene_loc) # join everything


# snp location
snp_loc <- read.csv("../MatrixQTL/MQTL_inputs/SNP_location.txt", sep=",") %>% select(1,2)
colnames(snp_loc)[1] <- "snps"
colnames(snp_loc)[2] <- "SNP_CHR"
trans_data <- left_join(trans_data, snp_loc)


# Filter for SNPs which are on a different chromsome
trans_data <- trans_data %>% filter(gene_CHR != SNP_CHR) # find SNPs which arent on same CHR as gene
dim(trans_data) # putative eQTLs on different chromosomes
length(unique(trans_data$gene)) # numbers

# Have SNPs, now need the genes
# Get gene list from trans-eQTLs
gene_list <- unique(trans_data$gene)

# filter cis data for same genes
cis_data <- cis_data %>% filter(gene %in% gene_list) # note here we have significant trans SNPs not associated with cis-SNPs

# cross reference
gene_list <- gene_list[gene_list %in% cis_data$gene]

# filter trans data
trans_data <- trans_data %>% filter(gene %in% gene_list)

length(gene_list) # 154 genes with a trans association on different chromsome and also cis association



# FINAL FILTERING OF eQTLS
trans_data <- trans_data %>% separate(., col = "snps", into = c("CHR", "POS", "REF", "ALT"), sep = ":")
trans_data <- trans_data %>% select(1,2,5)
trans_data <- trans_data %>% unite(., col = "TRANS_POS", "CHR", "POS", sep = ":")
colnames(trans_data)[2] <- "TRANS_GENE"
colnames(trans_data)[1] <- "FINAL_POS1" # trans positions

cis_data <- cis_data %>% separate(., col = "snps", into = c("CHR", "POS", "REF", "ALT"), sep = ":")

cis_data <- cis_data %>% select(1,2,5)
cis_data <- cis_data %>% unite(., col = "CIS_POS", "CHR", "POS", sep = ":")
colnames(cis_data)[2] <- "CIS_GENE"
colnames(cis_data)[1] <- "FINAL_POS2" # cis positions


LD_trans_final <- LD_trans_final %>% left_join(., trans_data)

LD_trans_final <- LD_trans_final %>% left_join(., cis_data)


# filter for TRANS-CIS of same GENE
LD_trans_final <- LD_trans_final %>% filter(TRANS_GENE == CIS_GENE)



###################################################################################################
######################## Random SNP-ciseQTL extraction ############################################
###################################################################################################

# NEED to filter RANDOM SNPs to get the same number of comparisons
# Cycle through each row in LD_trans
# count occurrence of each cis_SNP
# Subset that many rows with that SNP from the LD_random Final - i.e. that SNP is in the data X times same as above
# APPEND TO A FINAL TEST DATA SET
LD_RANDOM_FINAL <- data.frame(matrix(ncol = 6, nrow = 0)) # Final data frame
colnames(LD_RANDOM_FINAL) <- colnames(LD_trans_final) # same column names



set.seed(17897) # reproducibility
counter = 0 # to keep track


cis_pos = unique(LD_trans_final$FINAL_POS2) # SNPs we have to collect
for (row in 1:length(cis_pos)) { # cycle through all SNPs
  counter = counter + 1
  print(counter)
  cur_cis = cis_pos[row] # get cis Position
  count = as.numeric(length(which(LD_trans_final$FINAL_POS2 == cur_cis))) # see how many times that SNP is in our reference file
  potent = LD_random_final %>% filter(FINAL_POS2 == cur_cis) # filter for current pos
  random_row = sample_n(potent, count) # sample df for that position count times
  
  # Get columns in order
  # Will facilitate binding below
  random_row$TRANS_GENE = NA 
  random_row$CIS_GENE = NA
  
  # Bind to final
  LD_RANDOM_FINAL = rbind(LD_RANDOM_FINAL, random_row)
}

  
# Look at some plots
hist(LD_RANDOM_FINAL$R2)
hist(LD_trans_final$R2)

# some statistics
mean(LD_trans_final$R2)
mean(LD_RANDOM_FINAL$R2)

# Final data frame
Plotting = rbind(LD_trans_final, LD_RANDOM_FINAL)
Plotting$R <- sqrt(Plotting$R2) # square root for plotting

# Write to save it
write.table(Plotting, file = "Trans_Random_R2_plotting.txt", sep = "\t", row.names = F, quote = F, col.names = T)

# Dataframe for circos plot
Plotting %>% filter(Condition == "Trans-eQTL") %>% write.table(.,file = "Trans_R2_CircosPlotting.txt", sep = "\t", row.names = F, quote = F, col.names = T)

# remove variables
rm(snp_loc)
rm(LD_random_final_test)
rm(LD_trans_final_test)
rm(LD_RANDOM_FINAL)
rm(potent)
rm(cis_data)
rm(gene_loc)
rm(cis_snps)
rm(ran_snps)
rm(random_row)
rm(trans_data)
rm(trans_snps)

#########################################################################################################
#########################################################################################################
#################################### Plotting ###########################################################
#########################################################################################################
######################################################################################################### 




######################################
######### r2 #########################
######################################


# Plotting r2 overlapping
ggplot(data = Plotting, aes(x=R2, group=Condition, fill=Condition)) +
  geom_density(adjust=1.5, alpha = 0.4) +
  theme_ipsum()

# plotting r2 separately
ggplot(data = Plotting, aes(x=R2, group=Condition, fill=Condition)) +
  geom_density(adjust = 1) +
  theme_ipsum() +
  facet_wrap(~Condition) +
  theme(
    panel.spacing = unit(0.1, "lines"),
    axis.ticks.x=element_blank()
  )

######################################
######### r2 ##########################
######################################
ggplot(data = Plotting, aes(x=R, group=Condition, fill=Condition)) +
  geom_density(adjust=1.5, alpha = 0.5) +
  theme_ipsum()

ggplot(data = Plotting, aes(x=R, group=Condition, fill=Condition)) +
  geom_density(adjust=1.5) +
  theme_ipsum() +
  facet_wrap(~Condition) +
  theme(
    panel.spacing = unit(0.1, "lines"),
    axis.ticks.x=element_blank()
  )


ggplot(data = Plotting, aes(x=R, group=Condition, fill=Condition)) +
  geom_density(adjust=1.5, position="fill") +
  theme_ipsum()


wilcox.test(Plotting$R2 ~ Plotting$Condition, alternative = "less")
#interpret - R2 value for Trans-eQTLs is significantly higher than that of Random-SNPs
#data:  Plotting$R2 by Plotting$Condition
#W = 9.0313e+13, p-value < 2.2e-16
#alternative hypothesis: true location shift is not equal to 0


###################################################################################################
###################### Rain Cloud Plot ############################################################
###################################################################################################
     
ggplot(Plotting, aes(x = Condition, y = R2, fill = Condition)) + 
  ggdist::stat_halfeye(
    adjust = .5, 
    width = .6, 
    .width = 0, 
    justification = -.2, 
    point_colour = NA,
    alpha = 0.8
  ) + 
  geom_boxplot(
    width = .15, 
    outlier.colour = NA,
    alpha = 0.5,
    
  ) +
  #stat_boxplot(geom = 'errorbar', width = .15) +
    #ggdist::stat_dots(
     # side = "left",
      #justification = 1.1,
      #binwidth = .25
    #) +
  theme_tq() +
  labs(x = "Group",
    #y = parse(text='r^2'),
    fill = "Group"
  ) +
  coord_flip() +
  scale_fill_manual(values = c("Trans-eQTL" = "#d7191c", "Random" = "#2c7bb6"))
