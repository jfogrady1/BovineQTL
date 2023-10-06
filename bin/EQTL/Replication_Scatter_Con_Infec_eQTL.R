library(tidyverse)
library(dplyr)
#install.packages("ggpubr")
library(ggpubr)
library(tidyquant)
library(qvalue)
library(cowplot)
library(org.Hs.eg.db)
# Read in files 
# take the most significant SNP-gene pair in each data file and plot these together
# Arg1 TB_cis eQTLs
args = commandArgs(trailingOnly=TRUE)
TB_cis_raw <- as.data.frame(read.csv("../../results/EQTL/RAW/INFECTED_CisFDR5.txt", sep = "\t", header = T)) 


# Arg2 Control cis eQTLs
Con_cis_raw <- read.csv("../../results/EQTL/RAW/CONTROL_CisFDR5.txt", sep = "\t", header = T)

Con_cis_raw$type <- "Con"
TB_cis_raw$type <- "TB"


# Comparison of effect size of common-eQTL-gene associations
TB_cis_snp <- TB_cis_raw[TB_cis_raw$snps %in% Con_cis_raw$snps,]
Con_cis_snp <- Con_cis_raw[Con_cis_raw$snps %in% TB_cis_snp$snps,]
TB_cis_snp <- unite(TB_cis_snp, association, snps, gene, sep = ":")
Con_cis_snp <- unite(Con_cis_snp, association, snps, gene, sep = ":")
Con_cis_snp <- Con_cis_snp[Con_cis_snp$association %in% TB_cis_snp$association,]
TB_cis_snp <- TB_cis_snp[TB_cis_snp$association %in% Con_cis_snp$association,]
both_snp <- left_join(TB_cis_snp, Con_cis_snp, by = c("association"))
corr <- cor.test(both_snp$beta.x, both_snp$beta.y, method = "spearman", exact = F)
corr

# Plot
ggplot(data = both_snp, aes(x = beta.x, y = beta.y)) +
  geom_point(alpha = 0.8, size = 0.8) +
  theme_tq() +
  labs(x = "β value of cis-eQTLs in bTB-reactor group\n     (Matrix-eQTL slope)",
       y = "β value of cis-eQTLs in control group\n     (Matrix-eQTL slope)" ) +
  annotate("text", x=-1, y=2, label=paste0("rho = ",round(corr$estimate, 2)), size = 6.5) +
  annotate("text", x=-0.92, y=1.8, label = expression(italic(P) < 2.2~"x"~10^-16), size = 6.5) + # change this for publication
  geom_smooth(method = "lm", linetype = "dashed", col = "darkgrey", se = F) +
  theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black"),
        axis.text.y = element_text(angle = 0, size = 15, colour = "black"),
        legend.text = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 16, colour = "black"),
        axis.title.y = element_text(size = 16, color = "black"))

ggsave("Cis_replication.png", width=10, height=8, dpi=600)

# Comparison of top eQTL for common genes
head(Con_cis)
TB_cis <- TB_cis_raw[TB_cis_raw$gene %in% Con_cis_raw$gene,]
Con_cis <- Con_cis_raw[Con_cis_raw$gene %in% TB_cis_raw$gene,]
Con_cis <- Con_cis %>% dplyr::group_by(gene) %>% dplyr::summarise(FDR_con = min(FDR))
TB_cis <- TB_cis %>% dplyr::group_by(gene) %>% dplyr::summarise(FDR_TB = min(FDR))
df_both <- left_join(Con_cis, TB_cis)
head(Con_cis)

# Determine for each SNP-gene assocation, was it more significant in the Control or infected group
df_both$FDR_min <- if_else(df_both$FDR_con < df_both$FDR_TB, df_both$FDR_con, df_both$FDR_TB)


# Many SNPs will have same FDR so need to group by gene (to get the SNP gene pairs)
# Filter for the line which has the minimum FDR value of all the SNPs associated in cis with that gene
# Extract that first line for each gene
#df_both_final_really <- df_both_final %>% group_by(gene) %>% filter(FDR_min == min(FDR_min)) %>% slice(1)


df_both <- df_both %>% mutate(Group = case_when(df_both$FDR_TB < 0.005 & df_both$FDR_con >= 0.005 ~ "bTB reactor cis-eGene",
                                                df_both$FDR_con < 0.005 & df_both$FDR_TB >= 0.005 ~ "Control cis-eGene",
                                                df_both$FDR_con < 0.005 & df_both$FDR_TB < 0.005 ~ "General cis-eGene",
                                                df_both$FDR_con >= 0.005 & df_both$FDR_TB >= 0.005 ~ "Not significant"))

# Args3 = DE results
diff_expression <- read.csv("../../results/RNA-seq/DESEQ2/DE_Genes_Con_V_Infected.txt", sep = "\t")


df_both <- left_join(df_both, diff_expression, by = c("gene" = "Gene_ID"))
df_both$DE <- case_when(df_both$log2FoldChange >= 0 & df_both$padj < 0.05 ~ "DE up regulated",
                        df_both$log2FoldChange <= 0 & df_both$padj < 0.05 ~ "DE down regulated",
                        df_both$padj > 0.05 & df_both$Group == "bTB reactor cis-eGene" ~ "Not DE bTB reactor",
                        df_both$padj > 0.05 & df_both$Group == "Control cis-eGene" ~ "Not DE control")



tab <- table(df_both$Group)
ns <- as.character(tab[4])
btB <- as.character(tab[1])
cont <- as.character(tab[2])
General <- as.character(tab[3])
total <- as.numeric(ns) + as.numeric(btB) + as.numeric(cont) + as.numeric(General)
total
# 5209 genes common to bTB and reactor group which had a significant eQTL
5209 - 1119 - 2084 - 553
# Plot relationship between bTB and Control
ggplot(data = df_both, aes(x = -log10(FDR_TB), y = -log10(FDR_con), colour = Group)) + geom_point(size = 0.825, alpha = 1) +
  scale_color_manual(values = c("bTB reactor cis-eGene" = "#e08214", "Control cis-eGene" = "#8073ac", "General cis-eGene" = "darkgrey" , "Not significant" = "lightgrey")) +
  scale_x_continuous(limits = c(0.5,18), breaks = c(0,3,6,9,12,15,18)) +
  scale_y_continuous(limits = c(0.5,18), breaks = c(0,3,6,9,12,15,18)) +
  geom_vline(xintercept = -log10(0.005), col = "black", linetype="dashed", size = 0.5) + geom_hline(yintercept = -log10(0.005), col = "black", linetype = "dashed", size = 0.5) +
  labs(y=expression(-log[10](italic(P)[adj]~"Control")),
       x=expression(-log[10](italic(P)[adj]~"bTB reactor"))) +
  annotate("text", x=18, y=0.66, label=btB, col="#e08214", size = 4, fontface = "bold")  +
  annotate("text", x=0.66, y=18, label=cont, col = "#8073ac", size = 4, fontface = "bold" ) +
  annotate("text", x=18, y=18, label=General, col = "darkgrey", size = 4, fontface = "bold")+ 
  theme_bw() +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black"),
        axis.text.y = element_text(angle = 0, size = 15, colour = "black"),
        axis.title.x = element_text(size = 21, colour = "black", face = "bold"),
        axis.title.y = element_text(size = 21, color = "black", face = "bold"),
        panel.grid.minor = element_blank(),
        legend.title = element_text(size = 15, color = "black", face = "bold"),
        legend.text = element_text(size = 15)) +
  guides(color = guide_legend(override.aes = list(size = 4), title = "Group"))



ggsave("Cis_Con_Infec_comparison.png", width=10, height=8, dpi=600)

# Now highlight those genes which are DE between the two groups.
df_both_reduced <- df_both %>% filter(Group == "bTB reactor cis-eGene" | Group == "Control cis-eGene")
df_both_reduced <- df_both_reduced[order(df_both_reduced$diffexpressed, decreasing = T),]
df_both_reduced <- df_both_reduced %>% mutate(Alpha = if_else(DE == "DE up regulated" | DE == "DE down regulated", true = 1, false =0.1))

tab_btb <- df_both_reduced %>% filter(Group == "bTB reactor cis-eGene")
tab_bTB <- table(tab_btb$DE)
tab_bTB_down <- tab_bTB[1]
tab_bTB_up <- tab_bTB[2]

tab_con <- df_both_reduced %>% filter(Group == "Control cis-eGene")
tab_con <- table(tab_con$DE)
tab_con_down <- tab_con[1]
tab_con_up <- tab_con[2]

tab_con                                      
ggplot(data = df_both_reduced, aes(x = -log10(FDR_TB), y = -log10(FDR_con), colour = Group)) + geom_point(size = 2) +
  #scale_color_manual(values = c("bTB reactor cis-eGene" = "#e08214", "Control cis-eGene" = "#8073ac"), 1, 1, 0.2, 0.2) +
  scale_x_continuous(limits = c(0.5, 11),breaks = c(seq(1,11,1))) +
  scale_y_continuous(limits = c(0.5, 11),breaks = c(seq(1,11,1))) +
  geom_point(df_both_reduced, mapping = aes(x = -log10(FDR_TB), y = -log10(FDR_con), col=DE), alpha = df_both_reduced$Alpha) +
  scale_color_manual(values = c("DE up regulated" = alpha("#ca0020", 1), "DE down regulated" = alpha("#0571b0", 1), "Not DE bTB reactor" = alpha("#e08214", 0.15), "Not DE control" = alpha("#8073ac", 0.15), 
                                guide_legend(title="Differentially Expressed (DE)"))) +
  geom_vline(xintercept = -log10(0.005), col = "black", linetype="dashed", size = 0.5) + geom_hline(yintercept = -log10(0.005), col = "black", linetype = "dashed", size = 0.5) +
  xlab(bquote(~-Log[10]~italic(P)[adj]~bTB~reactor~group)) +
  ylab(bquote(~-Log[10]~italic(P)[adj]~Control~group)) +
  guides(color = guide_legend(override.aes = list(size = 2))) + 
  theme_test(base_size = 12) +
  annotate("text", x=10, y=2, label=paste(tab_bTB_up), col="#ca0020")  +
  annotate("text", x = 10.5, y = 2.05, label = "\u2191", size = 10, col = "#ca0020" ) +
  annotate("text", x = 10.5, y = 1.55, label = "↓", size = 10, col = "#0571b0" ) +
  annotate("text", x=10, y=1.5, label=tab_bTB_down, col = "#0571b0") +
  annotate("text", x= 0.95, y=10, label=paste(tab_con_up), col="#ca0020")  +
  annotate("text", x = 1.35, y = 10.1, label = "\u2191", size = 15, col = "#ca0020" ) +
  annotate("text", x = 2.05, y = 10.1, label = "↓", size = 15, col = "#0571b0" ) +
  annotate("text", x=1.65, y=10, label=tab_con_down, col = "#0571b0") +
  theme(axis.text.x = element_text(angle = 0, size = 20, colour = "black"),
        axis.text.y = element_text(angle = 0, size = 20, colour = "black"),
        legend.text = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 21, colour = "black", face = "bold"),
        axis.title.y = element_text(size = 21, color = "black", face = "bold"),
        legend.key.size = unit(3, "lines"))

ggsave("Cis_DE_overlayed.png", width=10, height=8, dpi=600)

bTB_genes <- df_both_reduced %>% filter(Group == "bTB reactor cis-eGene") %>% filter(DE != "Not DE bTB reactor")
Con_genes <- df_both_reduced %>% filter(Group == "Control cis-eGene") %>% filter(DE != "Not DE control")
dim(df_both)
df_both_final <- df_both %>% filter(Group == "bTB reactor cis-eGene" | Group == "Control cis-eGene" | Group == "General cis-eGene")
df_both_final <- df_both_final[order(df_both_final$diffexpressed, decreasing = T),]
df_both_final$DE <- if_else(df_both_final$Group == "General cis-eGene", true = "Not DE general", false = df_both_final$DE)
df_both_final$DE <- if_else(df_both_final$Group == "General cis-eGene" & df_both_final$padj < 0.05 & df_both_final$log2FoldChange > 0, true = "DE up regulated", false = df_both_final$DE)
df_both_final$DE <- if_else(df_both_final$Group == "General cis-eGene" & df_both_final$padj < 0.05 & df_both_final$log2FoldChange < 0, true = "DE down regulated", false = df_both_final$DE)
df_both_final <- df_both_final %>% mutate(Alpha = if_else(DE == "DE up regulated" | DE == "DE down regulated", true = 1, false =0.11))

tab_gen <- df_both_final %>% filter(Group == "General cis-eGene")
tab_gen <- table(tab_gen$DE)
tab_gen_down <- tab_gen[1]
tab_gen_up <- tab_gen[2]
Gen_genes <- df_both_final %>% filter(Group == "General cis-eGene") %>% filter(DE != "Not DE general")


ggplot(data = df_both_final, aes(x = -log10(FDR_TB), y = -log10(FDR_con), colour = Group)) + geom_point(size = 2) +
  #scale_color_manual(values = c("bTB reactor cis-eGene" = "#e08214", "Control cis-eGene" = "#8073ac"), 1, 1, 0.2, 0.2) +
  scale_x_continuous(limits = c(0.1, 19),breaks = c(seq(1,19,2))) +
  scale_y_continuous(limits = c(0.1, 19),breaks = c(seq(1,19,2))) +
  geom_point(df_both_final, mapping = aes(x = -log10(FDR_TB), y = -log10(FDR_con), col=DE), alpha = df_both_final$Alpha, size = 0.825) +
  scale_color_manual(values = c("DE up regulated" = alpha("#ca0020", 1), "DE down regulated" = alpha("#0571b0", 1), "Not DE bTB reactor" = alpha("#e08214", 0.25), "Not DE control" = alpha("#8073ac", 0.35), "Not DE general" = alpha("black", 0.25), 
                                guide_legend(title="Differentially Expressed (DE)"))) +
  geom_vline(xintercept = -log10(0.005), col = "black", linetype="dashed", size = 0.5) + geom_hline(yintercept = -log10(0.005), col = "black", linetype = "dashed", size = 0.5) +
  labs(y=expression(-log[10](italic(P)[adj]~"Control")),
       x=expression(-log[10](italic(P)[adj]~"bTB reactor"))) +
  theme_bw() +
  annotate("text", x=16.65, y=1.5, label=paste(tab_bTB_up), col="#ca0020",fontface = "bold", alpha = 0.8, size = 5)  +
  annotate("text", x = 17.4, y = 1.5, label = "\u2191", size = 10, col = "#ca0020",fontface = "bold", alpha = 0.8) +
  annotate("text", x = 19, y = 1.50, label = "↓", size = 10, col = "#0571b0",fontface = "bold", alpha = 0.8) +
  annotate("text", x=18.15, y=1.5, label=tab_bTB_down, col = "#0571b0", size = 5,fontface = "bold", alpha = 0.8) +
  annotate("text", x= 0.1, y=18.1, label=paste(tab_con_up), col="#ca0020", size = 5,fontface = "bold", alpha = 0.8)  +
  annotate("text", x = 0.9, y = 18.1, label = "\u2191", size = 10, col = "#ca0020",fontface = "bold", alpha = 0.8 ) +
  annotate("text", x = 2.05, y = 18.1, label = "↓", size = 10, col = "#0571b0", fontface = "bold", alpha = 0.8) + 
  annotate("text", x=1.45, y=18.1, label=tab_con_down, col = "#0571b0", size = 5,fontface = "bold", alpha = 0.8) +
  annotate("text", x=18.15, y=18.1, label=tab_gen_down, col = "#0571b0", size = 5,fontface = "bold", alpha = 0.8) +
  annotate("text", x = 19, y = 18.1, label = "↓", size = 10, col = "#0571b0",fontface = "bold", alpha = 0.8) +
  annotate("text", x=17.4, y=18.1, label = "\u2191", size = 10, col = "#ca0020",fontface = "bold", alpha = 0.8) +
  annotate("text", x=16.65, y=18.1, label = tab_gen_up, size = 5, col = "#ca0020",fontface = "bold", alpha = 0.8) +
  theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black"),
        axis.text.y = element_text(angle = 0, size = 15, colour = "black"),
        axis.title.x = element_text(size = 21, colour = "black", face = "bold"),
        axis.title.y = element_text(size = 21, color = "black", face = "bold"),
        panel.grid.minor = element_blank(),
        legend.title = element_text(size = 15, color = "black", face = "bold"),
        legend.text = element_text(size = 15)) +
  guides(color = guide_legend(override.aes = list(size = 4), title = "Gene Set"))


ggsave("Cis_DE_Gen_overlayed_FINAL.png", width=10, height=8, dpi=600)



bTB_genes <- bTB_genes[order(bTB_genes$padj, decreasing = F),]
Con_genes <- Con_genes[order(Con_genes$padj, decreasing = F),]
Gen_genes <- Gen_genes[order(Gen_genes$padj, decreasing = F),]

ALL_DE_cis_genes <- rbind(bTB_genes, Con_genes, Gen_genes)
ALL_DE_cis_genes <- ALL_DE_cis_genes[order(ALL_DE_cis_genes$padj, decreasing = F),]
write.table(bTB_genes, file = "bTB_DE_ciseGenes.txt", sep = "\t", col.names = T, row.names = F)
write.table(Con_genes, file = "Con_DE_ciseGenes.txt", sep = "\t", col.names = T, row.names = F)
write.table(Gen_genes, file = "Gen_DE_ciseGenes.txt", sep = "\t", col.names = T, row.names = F)
write.table(ALL_DE_cis_genes, file = "ALL_DE_ciseGenes.txt", sep = "\t", col.names = T, row.names = F)

