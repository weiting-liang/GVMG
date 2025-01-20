#dingqiuxia
library(ecodist)
library(dplyr)

# Batch read genetic distance files and generate NMDS datasets.
path<-"./Patristic_distance_forNMDS/"
filenames<-dir(path)
filepath<-sapply(filenames,function(x){paste(path,x,sep="/")})

data<-lapply(filepath,function(x){read.table(x)})

for (i in 1:length(filenames)){
  temp<-data[[i]]
  ddat<-as.dist(temp)
  gv.nmds<-nmds(ddat, mindim = 1, maxdim = 2, nits = 10)
  save(gv.nmds,file=paste(filenames[i],"nmds.rda",sep="."))
}

library(ggplot2)
PatristicDist<-read.table("distance_forNMDS/SGB769_PatristicDistMatrix.txt",head=T,row.names=1,check.names = F)
load("SGB769_PatristicDistMatrix.txt.nmds.rda")
# choose the best two-dimensional solution to work with
gv.nmin <- min(gv.nmds, dims=2)

# rotate the configuration to maximize variance
gv.rot <- princomp(gv.nmin)$scores

phe<-read.table("Source.phe.txt",head=T,row.names = 1,sep="\t")
phe<-as.matrix(phe)
phe1<-as.data.frame(phe[rownames(PatristicDist),])
rownames(phe1)<-rownames(PatristicDist)
colnames(phe1)<-colnames(phe)
#phe1$Region_of_samples[is.na(phe1$Region_of_samples)]<-"Isolate"
table(phe1$Source_group)

dat<-data.frame()
dat <- cbind(gv.nmin,phe1$Source_group) 
#dat <- na.omit(dat)
colnames(dat)<-c("MDS1","MDS2","Source")
dat$Source<-as.factor(dat$Source)

Region_phylogenetic_dist_nmds<-"Bifidobacterium breve"
library(ggplot2)
library(dplyr)
pdf("Bifidobacterium breve.source.pdf",width=12,height = 10)
ggplot(dat,mapping=aes(x=MDS1,y=MDS2,color=Source))+geom_point(size=3.75)+
  labs(x="MDS1",
       y="MDS2",
       title=Region_phylogenetic_dist_nmds)+scale_color_manual(values = c( "#88cda4","#097969","#305faa","#99b0b1","#E6E6E6","#BFBFBF","#7F7F7F","#ad4716","#aa8736","#E6F598","#FDAE61","#FEE08B"))+
  theme(legend.text = element_text(size = 15, face = 'bold'),
        legend.title = element_text(size = 15, face = 'bold'),
        axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        panel.background = element_blank(), 
       panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),        
      axis.line = element_line(size = 0.8),      
      axis.ticks = element_line(size = 0.8),     
      panel.border = element_rect(color = "black", fill = NA, size = 0.8)       )
dev.off()
