###R 4.0.5 liangweiting
workdir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(workdir)
rm(list = ls())

library('cowplot')
library('ggplot2')

data <- read.csv('taxa_level.csv', header=T)

p1 <- ggplot(data, aes(x= reorder(level, total), y=total))+
  geom_bar(stat='identity', color='black', fill= 'grey', width=0.8, position=position_dodge(width=0.9))+
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
  geom_bar(stat='identity', color='black', fill= '#E1B378', width=0.8, position=position_dodge(width=0.9))+
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

p3 <- ggplot(data, aes(x= reorder(level, total), y=archaea))+
  geom_bar(stat='identity', color='black', fill= '#E1B378', width=0.8, position=position_dodge(width=0.9))+
  coord_flip()+
  theme_bw()+
  scale_y_continuous(breaks = NULL, limits=c(0,1000))+
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank())+
  theme(panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank())+
  geom_text(aes(label=archaea), size=4, hjust=-0.5)+
  xlab('')+
  ylab('')+
  labs(title='no. of Archaeas')

p <- plot_grid(p1, p2, p3, labels = "AUTO", ncol=3, align='v', byrow = TRUE)

ggsave('taxa_number.pdf', p, width = 10, height = 5)
