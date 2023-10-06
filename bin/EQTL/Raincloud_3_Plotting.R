# 3. Getting everything together 

# Arguments
args = commandArgs(trailingOnly=TRUE)
# Trans_LD_comparisonSNPs.txt RandomTrans_LD_comparisonSNPs.txt Cis_LD_comparisonSNPs.txt LDR2FULL.interchrom.geno.ld LDR2RANDOMFULL.interchrom.geno.ld ${cis_snps} ${trans_snps} ${gene_loc} ${snp_loc}

#########################################################################################################
#########################################################################################################
#################################### Comparison of Groups ###############################################
#########################################################################################################
#########################################################################################################
library(tidyverse)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggdist)
library(gghalves)
library(tidyquant)
library(hrbrthemes)
library(ggsignif) 
hrbrthemes::import_roboto_condensed()
# Getting LD between cis and trans SNPs
# read in all data
#arg1 - arg3
trans_snps = read.table(args[1], header = T)
ran_snps = read.table(args[2], header = T, sep = ",") %>% select(1)
cis_snps = read.table(args[3], header = T)



# Format the DATA
trans_snps = trans_snps %>% separate(., x, into = c("CHR", "POS", "REF", "ALT"), sep = ":")
trans_snps = trans_snps %>% unite(., col = "FINAL_POS1", "CHR", "POS", sep =":")
ran_snps = ran_snps %>% separate(., x, into = c("CHR", "POS", "REF", "ALT"), sep = ":")
ran_snps = ran_snps %>% unite(., col = "FINAL_POS1", "CHR", "POS", sep =":")
cis_snps = cis_snps %>% separate(., x, into = c("CHR", "POS", "REF", "ALT"), sep = ":")
cis_snps = cis_snps %>% unite(., col = "FINAL_POS2", "CHR", "POS", sep =":")


# Do it for our LD SNPs
LD_trans = read.table(args[4], header = T) # interchromosomal
LD_trans <- LD_trans %>% unite(., FINAL_POS1, "CHR1", "POS1", sep = ":" ) # get positions ready for filtering
LD_trans <- LD_trans %>% unite(., FINAL_POS2, "CHR2", "POS2", sep = ":" )
LD_trans_clean <- LD_trans[LD_trans$FINAL_POS1 %in% trans_snps$FINAL_POS1,]
LD_trans_clean <- LD_trans_clean[LD_trans_clean$FINAL_POS2 %in% cis_snps$FINAL_POS2,]
rm(LD_trans)


# DO it for the Random SNPS
LD_random = read.table(args[5], header = T) # random interchromosomal
LD_random <- LD_random %>% unite(., FINAL_POS1, "CHR1", "POS1", sep = ":" ) # get positions ready for filtering
LD_random <- LD_random %>% unite(., FINAL_POS2, "CHR2", "POS2", sep = ":" )
LD_random_clean <- LD_random[LD_random$FINAL_POS1 %in% ran_snps$FINAL_POS1,]
LD_random_clean <- LD_random_clean[LD_random_clean$FINAL_POS2 %in% cis_snps$FINAL_POS2,]
rm(LD_random)

dim(LD_random_clean)
dim(LD_trans_clean)

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
cis_data <- read.csv(args[6], sep = "\t")
trans_data <- read.csv(args[7], sep = "\t") # all data with most points - used also in TWAS so need to plot that

cis_data <- cis_data %>% filter(cis_data$FDR < 0.01)
trans_data <- trans_data %>% filter(trans_data$FDR < 0.01)

# ONLY consider trans-SNPs on different chromosomes
# Arg[3] = gene location file used in normalisation procedure
gene_loc <- read.csv(args[8], sep = "\t") %>% select(1,2) # get gene location
colnames(gene_loc)[1] <- "gene"
colnames(gene_loc)[2] <- "gene_CHR"
trans_data <- left_join(trans_data, gene_loc) # join everything


# snp location
# Arg[4] = SNP location used in EQTL analysis
snp_loc <- read.csv(args[9], sep="\t") %>% select(1,2)
colnames(snp_loc)[1] <- "snps"
colnames(snp_loc)[2] <- "SNP_CHR"
trans_data <- left_join(trans_data, snp_loc)
trans_data <- left_join(trans_data, snp_loc)


# Filter for SNPs which are on a different chromsome
trans_data <- trans_data %>% filter(gene_CHR != SNP_CHR) # find SNPs which arent on same CHR as gene
dim(trans_data) # putative eQTLs on different chromosomes
length(unique(trans_data$gene)) # numbers (135 genes with a cis and trans association)

# Have SNPs, now need the genes
# Get gene list from trans-eQTLs
gene_list <- unique(trans_data$gene)

# filter cis data for same genes
cis_data <- cis_data %>% filter(gene %in% gene_list) # note here we have significant trans SNPs not associated with cis-SNPs

# cross reference
gene_list <- gene_list[gene_list %in% cis_data$gene]

# filter trans data
trans_data <- trans_data %>% filter(gene %in% gene_list)

length(gene_list) # 67 genes with a trans association on different chromosome and also cis association
# This is what we are comparing: Trans-eQTLs LD to Cis SNPs of the same gene



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

# Get in eQTL data
LD_trans_final <- LD_trans_final %>% left_join(., trans_data)
LD_trans_final <- LD_trans_final %>% left_join(., cis_data)


# filter for TRANS-CIS of same GENE
LD_trans_final <- LD_trans_final %>% filter(TRANS_GENE == CIS_GENE)


length(unique(LD_trans_final$FINAL_POS1)) # 622 SNPs in trans with an association in cis

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
colnames(LD_trans_final)


set.seed(17897) # reproducibility
counter = 0 # to keep track


cis_pos = unique(LD_trans_final$FINAL_POS2) # SNPs we have to collect
for (row in 1:length(cis_pos)) { # cycle through all SNPs
  counter = counter + 1
  cur_cis = cis_pos[row] # get cis Position
  count = as.numeric(length(which(LD_trans_final$FINAL_POS2 == cur_cis))) # see how many times that SNP is in our reference file
  potent = LD_random_final %>% filter(FINAL_POS2 == cur_cis) # filter for current pos
  random_row = sample_n(potent, count) # sample df for that position count times
  
  # Above - we are counting number of cis-eQTLs in our trans-eQTL file
  # using our Random LD file (for same number of trans-eQTL), we are extracting rows where the cis-eQTLs match
  # This will give us a file for comparison of LD between random SNPs and cis-eQTLs THAT IS THE SAME AS ABOVE
  # The only difference is the "trans" SNPs
  
  # Get columns in order
  # Will facilitate binding below
  random_row$TRANS_GENE = NA 
  random_row$CIS_GENE = NA
  
  # Bind to final
  LD_RANDOM_FINAL = rbind(LD_RANDOM_FINAL, random_row)
}


# Look at some plots
typeof(LD_RANDOM_FINAL$R2)

# some statistics
mean(LD_trans_final$R2)
mean(LD_RANDOM_FINAL$R2)

# Final data frame
Plotting = rbind(LD_trans_final, LD_RANDOM_FINAL)
Plotting$R <- sqrt(Plotting$R2) # square root for plotting

# Write to save it

write.table(Plotting, file = "Trans_Random_R2_plotting.txt", sep = "\t", row.names = F, quote = F, col.names = T)

# Dataframe for circos plot
# This file represents the R2 between trans-eQTLs and cis-eQTLs of the same gene
Plotting %>% filter(Condition == "Trans-eQTL") %>% write.table(.,file = "Trans_R2_CircosPlotting.txt", sep = "\t", row.names = F, quote = F, col.names = T)

# remove variables
rm(snp_loc)

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
Plotting <- read.csv("Trans_Random_R2_plotting.txt", sep = "\t")


######################################
######### r2 ##########################
######################################

#ggplot(data = Plotting, aes(x=R, group=Condition, fill=Condition)) +
 # geom_density(adjust=1.5) +
  #theme_ipsum() +
  #facet_wrap(~Condition) +
  #theme(
    #panel.spacing = unit(0.1, "lines"),
    #axis.ticks.x=element_blank()
  #)


#ggplot(data = Plotting, aes(x=R, group=Condition, fill=Condition)) +
 # geom_density(adjust=1.5, position="fill") +
  #theme_ipsum()


wilcox.test(Plotting$R2 ~ Plotting$Condition)
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
  stat_boxplot(geom = 'errorbar', width = .15) +
  theme_tq() +
  labs(x = "Group",
       y = expression(paste(italic("r")^2)),
       fill = "Group") + # Use atop() and italic() to italicize only the "Trans" part
  geom_signif(annotation = c(""), y_position = c(0.38), xmin = c(1), xmax = c(2), tip_length = c(0.5, 0.04)
  ) +
  annotate("text", x=1.5, y=0.39, label = expression(italic(P) < 2.2~"x"~10^-16), angle = 270, size = 6.5) +
  coord_flip() +
  scale_fill_manual(values = c("Trans-eQTL" = "#d7191c", "Random" = "#2c7bb6")) +
  theme(axis.text.x = element_text(angle = 0, size = 20, colour = "black"),
        axis.text.y = element_text(angle = 0, size = 20, colour = "black"),
        legend.text = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 25, colour = "black", face = "bold"),
        axis.title.y = element_text(size = 25, color = "black", face = "bold"),
        legend.position = "none")
ggsave("Raincloud_plot.png", width = 15, height = 8, dpi = 1000)
??annotate
