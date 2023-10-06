
####### Gene level eQTL analysis
args = commandArgs(trailingOnly=TRUE)

library(ggplot2)
library(tidyverse)
library(dplyr)
library(stringr)
library(MatrixEQTL)

# Already have SNP data ok, just need to sort out the expression (convert to symbol)
# And convert the gene_locations (symbols) and keep the start and end of gene



# Set paramaters such as model, name of genotype and expression file
useModel = modelLINEAR
expression_file_name <- args[1] # ARG 1
args1 <- args[1]
out <- unlist(strsplit(args1, "_"))[1]
covariate_file_name <- args[2] # Arg2
genotype_file_name <- args[3] # ARG3


# Sort out the SNP location
snp_location <- as.data.frame(read.csv(args[4], sep=",", header = T)) # arg4


# sort out gene_location
gene_location <- read.csv(args[5], sep = ",") # arg 5

# Set the p-value threshold
pvOutputThreshold_tra = 0.00005;# need this as we are considering local SNP gene pairs

# Define covariance matrix for error term

# Open and load in the data
snps = SlicedData$new();
snps$fileDelimiter = ",";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in pieces of 2,000 rows
snps$LoadFile( genotype_file_name );

gene = SlicedData$new();
gene$fileDelimiter = ",";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in pieces of 2,000 rows
gene$LoadFile( expression_file_name );

cvrt = SlicedData$new();
cvrt$fileDelimiter = ",";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
cvrt$fileSliceSize = 1;      # read file in pieces of 2,000 rows
cvrt$LoadFile( covariate_file_name );


me <- Matrix_eQTL_main(snps = snps, 
                       gene = gene,
                       cvrt = cvrt,
                       errorCovariance = numeric(),
                       output_file_name = paste0(out,"_Trans_RAW.txt"), 
                       pvOutputThreshold = 0.00001,
                       verbose = TRUE, 
                       output_file_name.cis = paste0(out,"_Cis_RAW.txt"), # ../TWAS/TWAS_outputs/Mediator_
                       pvOutputThreshold.cis = 0.05,
                       snpspos = snp_location, 
                       genepos = gene_location,
                       cisDist = 1000000,
                       pvalue.hist = "qqplot",
                       useModel = modelLINEAR,
                       min.pv.by.genesnp = FALSE ,
                       noFDRsaveMemory = FALSE)






# Get number of eGenes
length(unique(me$cis$eqtls$gene))
length(which(me$cis$eqtls$FDR < 0.05))


# Getting cis and trans eQTLs with significant FDR
cis_eqtl <- me$cis$eqtls[me$cis$eqtls$FDR < 0.05,]
cis_genes <- unique(cis_eqtl$gene)
trans_eqtl <- me$trans$eqtls[me$trans$eqtls$FDR < 0.05,]
trans_genes <- unique(trans_eqtl$gene)




write.table(cis_eqtl, file = paste0(out, "_CisFDR5.txt"), sep = "\t", row.names = F, quote = F)
write.table(trans_eqtl, file = paste0(out, "_TransFDR5.txt"), sep = "\t", row.names = F, quote = F)

#genes
write.table(cis_genes, file = paste0(out, "_CisGeneSIG.txt"), sep = "\t", row.names = F, quote = F)
write.table(trans_genes, file = paste0(out, "_TransGeneSIG.txt"), sep = "\t", row.names = F, quote = F)
