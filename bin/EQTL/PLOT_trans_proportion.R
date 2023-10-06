## Plotting
library(ggridges)
library(tidyquant)
library(tidyverse)
library(dplyr)
args = commandArgs(trailingOnly=TRUE)
STATS <- read.table(file = "/home/workspace/jogrady/TF_PROP_STATS.csv", header = T)
Plotting <- read.table(file = "TF_Prop_RAW.txt", header = T)
trans <- Plotting %>% filter(Condition == "Trans-eQTL") %>% dplyr::select(1,2)
trans <- trans %>% group_by(distance) %>% dplyr::summarise(Ratio_Trans = mean(Prop_with_TFs))
trans$Condition = "Trans-eQTL"
colnames(trans)[2] <- "Ratio"
Random <- Plotting %>% filter(Condition == "Random") %>% dplyr::select(1,2)
Random$Condition <- "Random"
colnames(Random)[2] <- "Ratio"


random_point <- Random %>%
 filter(row_number() %% 100 == 1)

head(random_point)
dim(random_point)
data_plot <- rbind(trans,Random)
data_plot$distance <- data_plot$distance / 1000
distance = unique(data_plot$distance)
Plotting_simple <- data_plot %>% filter(distance >= 0) 


sum_plot <- data_plot %>% group_by(distance, Condition) %>% dplyr::summarize(Ratio = mean(Ratio))

# The Plot
ggplot(data = data_plot, aes(x = Ratio, y = distance, group = distance, colour = Condition)) +
  #geom_point(size = 0.5, alpha = 1) +
  theme_bw() +
  labs(x = "Proportion",
       y = "Distance (000's bps)") +
  scale_color_manual(values = c("Trans-eQTL" = "#d7191c", "Random" = "#2c7bb6"), name = "") +
  scale_y_continuous(breaks= distance) +
  theme(axis.text.x = element_text(angle = 90)) + coord_flip() +
  geom_boxplot(data = Plotting_simple, alpha = 0.5, outlier.colour = NULL, outlier.shape = NA, fill = "#2c7bb6") + coord_flip() +
  geom_line(data = sum_plot[sum_plot$Condition == "Trans-eQTL",], col = "#d7191c", linewidth = 0.5, alpha = 1, group = 1) +
  geom_line(data = sum_plot[sum_plot$Condition == "Random",], col = "#2c7bb6", linewidth = 0.5, alpha = 1, group = 1) +
  xlim(0,1) +
  theme(axis.text.x = element_text(angle = 90, size = 11.5, colour = "black"),
        axis.text.y = element_text(angle = 0, size = 18, colour = "black"),
        legend.text = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 25, colour = "black"),
        axis.title.y = element_text(size = 25, color = "black"),
        legend.position = "none",
        legend.key.size = unit(3, "lines"))  # adjust legend key size)
ggsave("Trans_TF_proportions.pdf", width = 15, height = 8, dpi = 600)