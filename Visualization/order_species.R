###liangweiting
###pic order and species proportion
###R 4.2.3
workdir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(workdir)
rm(list = ls())

#packages
library(readxl)
library(ggplot2)

#data
data <- read_excel("species_100mags_order.xlsx")

#count top5 order
order_counts <- table(data$order)
sorted_order_counts <- sort(order_counts, decreasing = TRUE)
Top5_order <- names(sorted_order_counts)[1:5]

#change order input
data$order <- ifelse(data$order %in% Top5_order, data$order, "Others")

#pic
color <- c( "#56B4E9", "pink","#F0E442","#b2d183","#999999","#E69F00")
#species_data$order <- factor(species_data$order, levels=top5$Order)
#class(species_data$mags_count)
p3 <- ggplot(data, aes(x=reorder(species, mags_count), y=mags_count, fill=order))+
  geom_bar(stat='identity')+
  scale_fill_manual(values=color)+
  coord_flip()+
  xlab('species(assembly mags>100)')+
  ylab('number of MAGs')+
  theme_bw()

ggsave('order.pdf', p3, width=8, height=10)


