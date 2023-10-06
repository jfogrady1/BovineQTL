###############################################################
###                                                        ####
###  Figure - Circos plot of LD between cis and Trans SNPs ####
###                                                        ####  
###############################################################
args = commandArgs(trailingOnly=TRUE)

library(tidyverse)
library(plyr)
library(dplyr)
library(ComplexHeatmap)
library(devtools)
library(gridBase)
#devtools::install_github("jokergoo/circlize")
library(circlize)
#BiocManager::install("BSgenome.Btaurus.UCSC.bosTau9")
library("BSgenome.Btaurus.UCSC.bosTau9")




chromInfo = as.data.frame(read.chromInfo(species = "bosTau9")$df)

chrs <-  paste0("chr", c(1:29))

# Read in the LD data
LD <- read.table(args[1], header = T) %>% select(1,2,3,7)
LD <- as.data.frame(LD)
LD <- LD %>% separate(., col = "FINAL_POS1", into = c("chr1", "Pos1"), sep = ":")
LD <- LD %>% separate(., col = "FINAL_POS2", into = c("chr2", "Pos2"), sep = ":")

LD <- LD %>% mutate(COLOUR = case_when(
  R2 >= 0.05 & R2 < 0.10 ~ "#ffeda0" ,
  R2 >= 0.10 & R2 < 0.15 ~ "#feb24c" ,
  R2 >= 0.15 & R2 < 0.20 ~ "#f03b20" ,
  R2 >= 0.20 ~ "#000000"))


LD$chr1 <- sub("^","chr",LD$chr1)
LD$chr2 <- sub("^","chr",LD$chr2)
table(LD$COLOUR)
fontsize(20)
png(filename = "Trans_Cis_LD.png", width = 12, height = 12, units = "in", res = 800)
circos.par("start.degree" = 90)
df = data.frame(
  chr  = chromInfo$chr,
  start = chromInfo$start,
  end   = chromInfo$end)
df$annotation = "bovine"
df$year = "2023"
circos.initializeWithIdeogram(df, chromosome.index = chrs)
LD_plot <- LD %>% filter(!is.na(COLOUR))
LD_plot <- LD_plot[order(LD_plot$R2, decreasing = F),]

head(LD_plot)
for (row in 1:nrow(LD_plot)) {
  sector_index1 = LD_plot[row,"chr1"]
  point_1 = as.numeric(LD_plot[row,"Pos1"])
  sector_index2 = LD_plot[row,"chr2"]
  point_2 = as.numeric(LD_plot[row,"Pos2"])
  col_line = LD_plot[row,"COLOUR"]
  circos.link(sector.index1 = sector_index1, point1 = point_1, sector.index2 = sector_index2, point2 = point_2, col = col_line)
}

lgd_lines = Legend(at = c("0.05-0.10", "0.10-0.15", "0.15-0.20", "0.20+"), type = "lines", 
                   legend_gp = gpar(col = c("#ffeda0","#feb24c","#f03b20","#000000"), lwd = 8.5, fontsize = 80), title_gp = gpar(col = "black", fontsize = 30), title_position = "topleft", 
                   title = expression(italic(r)^2))

draw(lgd_lines, x = unit(10, "in"), y = unit(1.25, "in"), just = c("left", "bottom"))

dev.off()
circos.clear()
