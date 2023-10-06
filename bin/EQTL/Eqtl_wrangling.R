############################################################################
############################################################################
####################### EQTL analysis ######################################
############################################################################
############################################################################



args = commandArgs(trailingOnly=TRUE)
# load in selected libraries
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(MatrixEQTL)
library(stringi)
library(parallel)
??parallel
#############################################################################
####################### Data input ##########################################
#############################################################################

#######################-------------------------------------------------------
## 1. Genotype data ###-------------------------------------------------------
#######################-------------------------------------------------------
args1 <- args[1]
out <- unlist(strsplit(args1, "_"))[1]
# Read in the data and select the column ID and the samples defined above
genotype <- read.table(args[3], sep = "\t", header = T) %>% select(-X)


head(genotype)
# Change the column name of column 1 to match matrixQTL example
colnames(genotype)[1] <- "id"
rownames(genotype) <- genotype$id


#########################-------------------------------------------------------
## 2. Expression data ###-------------------------------------------------------
#########################-------------------------------------------------------


# Read in the file and select the columns we need as above
expression <- read.csv(args[1] ,sep="\t")
head(expression)

# Change the column name of column 1 to id
colnames(expression)[1] <- "id"

# As above, need column 1 to be converted to numeric so replace GeneID:
# with nothing and convert to numeric


expression$id <- str_replace_all(expression$id, "GeneID:", "")
expression$id <- as.numeric(expression$id)


## NB. Below we will change the rows (gene id) in expression file to ensure
## That we only have genes with actual locations available
## This is to ensure that we can properly match cis-eQTLs.

##################################################-------------------------------------------------------
## 3. PEER data (contains measured covariates) ###-------------------------------------------------------
##################################################-------------------------------------------------------

# Read in the file and set as  a matrix
Peer_factors <- read.csv(args[2])

# Change the row names to equal the column name of the expression file
# This will be important when we transpose the file down below
colnames(Peer_factors) <- colnames(expression)
rownames(Peer_factors) <- 1:nrow(Peer_factors)
# Set up an empty vector equal to the number of **Current** columns
id <- c(1:length(colnames(Peer_factors)))

# Reset the id column variable to something more meaningful
Peer_factors[,"id"] <- c(1:length(rownames(Peer_factors)))

# Set to a character: this is the case in the example documentation
Peer_factors[,"id"] <- as.character(Peer_factors[,"id"])
Peer_factors
# Set the column names to numbers rather than V1 etc

# Write the data to the a new file
# Transpose the data, make sure there are no quotes --> essential

# Very important to note that there are 15 hidden factors and 7 measured covariates (Batch Infection and PC1-5)
# in the output from peer. Initially thought we would have to merge both peer file
# and covariate file but this is not the case



#############################-------------------------------------------------------
## 4. Gene location file ####-------------------------------------------------------
#############################-------------------------------------------------------
# At this point, run the TSS_determination.py script
# Note require an annotation file with coordinates of features on plus and minus strands
# This file is required to show the location (bp) and chr of gene features
# Get this information from master annotation file
# Read in the file

gene_data <- read.csv(args[5], sep = "\t", header = T)

# Create a variable with the same ids as expression ids (above)
ids <- as.numeric(expression$id)
ids <- expression$id



# Subset the genes in the gene data file: the final genes are ones which passed 
# expression filtering and normalisation
head(gene_data)
gene_data <- subset(gene_data, Geneid %in% ids)

# Check the dimensions, should be the same as expression dimensions
dim(gene_data)

# Convert to numeric to ensure it is read properly by matrix qtl
gene_data$Geneid <- as.numeric(gene_data$Geneid)

# Remove any NAs
gene_data <- gene_data %>%
  drop_na()



# Set the first column name as id to match expression
colnames(gene_data)[1] <- "id"

# ensure rownames are matched so MatrixQTL can match up
rownames(expression) <- expression$id

# Get rid of any characters so remove X and MT chr names
# These shouldn't have any impact on results as we do not have SNPs for them
gene_data["Chromosome"][gene_data["Chromosome"] == "X"] <- 30
gene_data["Chromosome"][gene_data["Chromosome"] == "MT"] <- 31



head(gene_data, 1000)
dim(gene_data)

#Convert to numeric matrix
gene_data <- (apply(gene_data, 2, as.numeric))




###########################-------------------------------------------------------
## 5. SNP position file ###-------------------------------------------------------
###########################-------------------------------------------------------


# Rename the columns as is in MatrixQTL documentation
snp_data <- read.table(args[4], sep="\t", header = F)
head(snp_data)
colnames(snp_data) <- c("snp", "chr", "pos")
snp_data <- as.data.frame(snp_data)



rownames(snp_data) <- snp_data[,1]
rownames(gene_data) <- gene_data[,1]
# Check dimensions of all to make sure they are ok
dim(genotype)
dim(expression)
dim(Peer_factors)
dim(gene_data)
dim(snp_data)

##########################################
## 7. Write all files to the directory ###
##########################################

# These files are then opened by MatrixQTL
write.table(genotype, file = paste0(out,"_GT_MatrixQTL_format.txt"), sep = ",", quote = F, row.names = F)
write.table(expression, file = paste0(out,"_EXP_MatrixQTL_format.txt"), sep = ",", quote = F, row.names = F)
write.table(Peer_factors, file = paste0(out,"_COV_MatrixQTL_format.txt"), sep = ",", quote = F, col.names = T, row.names = F)
write.table(snp_data, file = "SNP_loc_MatrixQTL_format.txt", sep = ",", row.names = FALSE, quote = F)
write.table(gene_data, file = "gene_loc_MatrixQTL_format.txt", sep = ",", row.names = FALSE, quote = F)
