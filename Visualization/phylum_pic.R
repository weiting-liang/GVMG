###stacked hist phylum SGB
#R 4.2.3
workdir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(workdir)
rm(list = ls())

library(reshape2)
library(ggplot2)

#data
data <- read.csv('bac_phylum_sgb.csv')

data2_long <- melt(data[, c(1:3)])
data2_long$phylum <- factor(data2_long$phylum, levels=c(data$phylum[order(data$all)]))
color <- c('#1F78B4', '#A6CEE2')


p <- ggplot()+
  geom_bar(aes(y=value, x=phylum, fill=variable), data=data2_long,
           stat='identity', position='stack')+
  geom_text(aes(x=phylum, y=all+10, label=uSGB_proportion_text), data=data, size=3)+
  labs(y='Number of Species', title='Novel Species by Phylum')+
  scale_fill_manual(values=color)+
  theme_bw()+
  coord_flip()+
  theme(panel.grid.major.y = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank())

ggsave('Phylum_usgb.pdf', p, width = 10, height = 6)




