library("data.table")
library("tidyverse")

Gtex <- read.table("eqtl_study/eqtl_nextflow/data/EQTL/Blood_Trans_eQTL_FDR_0.05_filter.txt", header = T)

reference_inter <- read.table("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/ALL_interchromosomal_LD_corrected.txt")
reference_intra <- read.table("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/ALL_intrachromosomal_LD_corrected.txt") %>% mutate(V11 = "TRUE")
reference <- rbind(reference_inter, reference_intra)


Gtex$SNP <- gsub("_", ":",Gtex$SNP)
#######################################################################################
########################################################################################
## 1. AC methodology
########################################################################################
########################################################################################

#To calculate AC, we took the SNP-gene combination with the lowest p-value in the discovery cohort 
#for each gene if it was significant and matched these with the same SNP-gene pairs in the replication cohort if it was also significant. 
#We then determined the percentage of those significant SNP-gene pairs had the same allelic direction of effect compared to the discovery cohort. 

# Note, significant in both so will need the two files, permutation (gene level significance and nominal p value)
#args[1] <- "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/ALL.cis_qtl_fdr0.05.txt"
#args[2] <- "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/CONTROL.cis_qtl_fdr0.05.txt"
#args[3] <-  "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/INFECTED.cis_qtl_fdr0.05.txt"
#args[4] <- "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/EQTL/Blood.nominals.2rd.txt"
#args[5] <- "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/EQTL/Gtex_thresholds_determined.txt"
# Read in the data
# select top snp which were significant genome wide - note we are comparing the top eQTL in our dataset 

head(reference)
# match the SNPs together.

Gtex_AC <- Gtex %>% filter(SNP %in% reference$V4)

head(Gtex_AC)
reference_AC <- reference %>% filter(reference$V4 %in% Gtex_AC$SNP)
head(reference_AC)
