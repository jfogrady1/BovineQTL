#Script to plot allele frequencies

# Use following as reference to see what we would be expecting
#https://www.researchgate.net/figure/Different-patterns-of-allele-frequencies-in-the-EAF-plot-These-different-patterns-have_fig2_261881356

library(tidyverse)
library(ggplot2)
library(tidyquant)
library(ggpointdensity)
library(viridis)

args = commandArgs(trailingOnly=TRUE)


# WGS file of EU samples
# Here specify the HF file if required
# ouput form vcftools -freq command
WGS_file <- read.table(args[1], sep ="\t", header = F, skip = 1)
colnames(WGS_file) <- c("CHR", "POS", "Allele_number", "Total_alleles", "REF", "ALT")
WGS_file <- WGS_file %>% separate(., REF, into = c("REF_allele", "REF_AF"), sep = ":")
WGS_file <- WGS_file %>% select (1,2,5,6)


# SNP file
Target_file <- read.table(args[2], sep = "\t", header = F, skip = 1)
colnames(Target_file) <- c("CHR", "POS", "Allele_number", "Total_alleles", "REF", "ALT")
Target_file <- Target_file %>% separate(., REF, into = c("T_REF_allele", "T_REF_AF"), sep = ":")
Target_file <- Target_file %>% select(1,2,5,6)


# Merge Files
# Use left_join as need to keep the SNPs in the WGS file which came from the intersection
# This will avoid NAs
merged = left_join(WGS_file, Target_file, by = c("CHR"="CHR", "POS" = "POS"))

merged$REF_AF <- as.numeric(merged$REF_AF)
merged$T_REF_AF <- as.numeric(merged$T_REF_AF)

## Correlation
# note 3 rows have missing values which are not included
# These are Insertions not picked up by genotyping platform

merged = drop_na(merged)
test = cor.test(merged$REF_AF, merged$T_REF_AF, method = "spearman", exact = F)
test$estimate
write.table(merged, "RAW_AF_Comparsion_Values.txt", row.names = F, col.names = T)
# Ggplot - 
ggplot(data = merged, aes(x = REF_AF, y = T_REF_AF)) +
  geom_point(size = 0.25, color = "#2c7bb6", alpha = 1) +
  labs(x = "Allele frequency in European WGS reference set", 
       y = "Allele frequency in target set") +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  scale_x_continuous(breaks = seq(0, 1, 0.2)) +
  annotate('segment', x = 0, xend = 1, y = 0, yend = 1,
           linewidth = 1,
           alpha = 0.8,
           color = "#d7191c") +
  theme_tq(base_size = 14) +
  theme(axis.title  = element_text(colour = "black"))

args[3]
ggsave(args[3], width=10, height=8, dpi = 600)

png(args[4], width = 10, height = 8, units = "in", res = 600)
smoothScatter(merged$REF_AF,merged$T_REF_AF,
              pch = 10, col = "red", nrpoints = 150, bandwidth = 0.04, ret.selection = F,
              nbin = 120, xlab = "Allele frequency in European WGS reference set",
              ylab = "Allele frequency in target set")

dev.off()


# base R
Points.2.remove <- smoothScatter(merged$REF_AF,merged$T_REF_AF,
                                 pch = 10, col = "red", nrpoints = 150, bandwidth = 0.04, ret.selection = T,
                                 nbin = 120, xlab = "Allele Frequency in European WGS Reference Set",
                                 ylab = "Allele Frequency in Target Set")

Points.2.remove
# Remove nrPoints identified above
SNPs.2.remove <- merged[Points.2.remove,] %>% select(CHR, POS, REF_AF,T_REF_AF)
SNPs.2.remove.positions <- merged[Points.2.remove,] %>% select(CHR, POS)

colnames(SNPs.2.remove) <- c("CHR", "POS", "WGS_REF_AF", "Target_REF_AF")
write.table(SNPs.2.remove, file = args[5], quote = F, col.names = T, row.names = F, sep = "\t")
write.table(SNPs.2.remove.positions, file = args[6], quote = F, col.names = F, row.names = F, sep = "\t")