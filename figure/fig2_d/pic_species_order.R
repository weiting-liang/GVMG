###liangweiting
###pic order and species proportion
###R 4.2.3
workdir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(workdir)
rm(list = ls())

#packagase
library(ggplot2)
library(tidyr)

#data
order_data <- read.csv('order_mags.txt', header=FALSE)
species_data <- read.csv('species_count.csv')

#order_data
order_df1 <- as.data.frame(table(order_data))
order_df1 <- order_df1[order(order_df1$Freq, decreasing=TRUE), ]
top5 <- order_df1[1:5, ]
top5$V1 <- as.character(top5$V1)
other_count <- sum(order_df1$Freq) - sum(top5$Freq)
top5[6, ] <- c('Others', other_count)
top5[, 'x'] <- 'Order'
top5$Freq <- as.numeric(top5$Freq)
top5[, 'proportion'] <- top5$Freq/sum(top5$Freq)*100


color <- c("#b2d183", "#999999", "#E1B378","#56B4E9", "#F0E442","#E69F00")

top5 <- top5[order(top5$proportion, decreasing=TRUE), ]
top5$V1 <- factor(top5$V1, levels=top5$V1)

colnames(top5) <- c('Order', 'Freq', 'x', 'proportion')

###species order
species_data$order <- factor(species_data$order, levels=top5$Order)
p3 <- ggplot(species_data, aes(x=reorder(species, mags_count), y=mags_count, fill=order))+
  geom_bar(stat='identity')+
  scale_fill_manual(values=color)+
  coord_flip()+
  xlab('species(assembly mags>100)')+
  ylab('number of MAGs')+
  theme_bw()
  
ggsave('species_order_pic.pdf', p3, width=10, height=8)
