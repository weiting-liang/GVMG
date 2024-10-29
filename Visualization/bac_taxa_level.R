###R 4.0.5 
###liangweiting
workdir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(workdir)
rm(list = ls())

library('cowplot')
library('ggplot2')

data <- read.csv('bac_taxa_level.csv', header=T)

p1 <- ggplot(data, aes(x= reorder(level, total), y=total))+
  geom_bar(stat='identity', color='black', fill= 'blue', width=0.8, position=position_dodge(width=0.9))+
  coord_flip()+
  theme_bw()+
  scale_y_continuous(breaks = NULL, limits=c(0,1000))+
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text= element_text(size=12))+
  theme(panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank())+
  geom_text(aes(label=total), size=4, hjust=-0.5)+
  xlab('')+
  ylab('')+
  labs(title='Total no. of SGBs from vagina')

p2 <- ggplot(data, aes(x= reorder(level, total), y=new))+
  geom_bar(stat='identity', color='black', fill= 'purple', width=0.8, position=position_dodge(width=0.9))+
  coord_flip()+
  theme_bw()+
  scale_y_continuous(breaks = NULL, limits=c(0,1000))+
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank())+
  theme(panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank())+
  geom_text(aes(label=new), size=4, hjust=-0.5)+
  xlab('')+
  ylab('')+
  labs(title='no. of uSGBs')


p <- plot_grid(p1, p2, labels = "AUTO", ncol=2, align='v', byrow = TRUE)

ggsave('f2b.pdf', p, width = 7, height = 3)
