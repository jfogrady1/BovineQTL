## Plotting
library(ggridges)
library(tidyquant)
library(tidyverse)
library(dplyr)
args = commandArgs(trailingOnly=TRUE)
STATS <- read.table(file = args[1], header = T)
Plotting <- read.table(file = args[2], header = T)

trans <- Plotting %>% filter(Condition == "Trans-eQTL") %>% dplyr::select(1,2)
trans <- trans %>% group_by(distance) %>% dplyr::summarise(Ratio_Trans = mean(Prop_with_TFs))
trans$Condition = "Trans-eQTL"
colnames(trans)[2] <- "Ratio"
Random <- Plotting %>% filter(Condition == "Random") %>% dplyr::select(1,2)
Random$Condition <- "Random"
colnames(Random)[2] <- "Ratio"


data_plot <- rbind(trans,Random)
data_plot$distance <- data_plot$distance / 1000
distance = unique(data_plot$distance)
Plotting_simple <- data_plot %>% filter(distance >= 200) 


sum_plot <- data_plot %>% group_by(distance, Condition) %>% dplyr::summarize(Ratio = mean(Ratio))


# The Plot
ggplot(data = data_plot, aes(x = Ratio, y = distance, group = distance, colour = Condition)) +
 geom_point(size = 0.2, alpha = 1) +
theme_tq() +
labs(x = "Proportion",
     y = "Distance (000's bps)") +
scale_color_manual(values = c("Trans-eQTL" = "#d7191c", "Random" = "#2c7bb6")) +
scale_y_continuous(breaks= distance) +
theme(axis.text.x = element_text(angle = 90)) + coord_flip() +
geom_density_ridges(data = Plotting_simple, rel_min_height=0.0002, alpha = 0.2, quantile_lines = TRUE, quantiles = c(0.975), scale = 0.9, fill = "#2c7bb6") + coord_flip() +
geom_line(data = sum_plot[sum_plot$Condition == "Trans-eQTL",], col = "#d7191c", size = 0.5, alpha = 0.8, group = 1) +
geom_line(data = sum_plot[sum_plot$Condition == "Random",], col = "#2c7bb6", size = 0.5, alpha = 0.8, group = 1)

ggsave("Trans_Gene_proportions.tiff", width = 15, height = 8, dpi = 600)