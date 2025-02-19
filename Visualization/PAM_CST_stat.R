#haolilan
library(dplyr)
library(reshape2)
library(tibble)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(readxl)

###
setwd("E:/OneDrive - BGI Hong Kong Tech Co., Limited/projects/Merge_multi_cohorts/VGC2024/")
outdir <- paste0(getwd(),"/")
theme_set(theme_classic()+theme(axis.text.x = element_text(angle = 90, vjust=0.5,hjust = 1)))

## CST ####################################
### data #####
prof.scaled <- read.delim("Whole.relative_taxonomic_abundance.recal.txt")
cst.stat  <- read.table(paste0(outdir,"cluster.scale.C12.Sample.txt"),header = T)
cst.avg.top3  <- read.table(paste0(outdir,"cluster.scale.C12.TaxaTop3.txt"),header = T)
## cluster samplenum via cohort, new 
phen.cluster <- read.delim(paste0(outdir,"phen.cluster.txt"),header = T)
phen.cluster$Samp <- phen.cluster$Prof_ID_update

## CST: China & America #####
### prop.test #####
prefix <- "AmericaChina" 
phen.cluster$group <- phen.cluster$Country
phen.input <- phen.input%>%subset(group%in%c("China","USA"))%>%select(Samp,group)%>%.[!duplicated(.),]
phen.input$group <- factor(phen.input$group,levels = c("China","USA"))
table(phen.input$group)
##
data <- merge(cst.stat,phen.input,by="Samp",all.y = T)
data.m<- data%>%group_by(clustering,group)%>%summarise(count=n())%>%dcast(.,clustering~group)
data.m[is.na(data.m)]<-0
colnames(data.m)
colnames(data.m) <- c("rename","case","control")
## any>0
raw <- data.m%>%dplyr::select(rename,case,control)%>%
  subset(case>0|case>0)
raw[is.na(raw)] <- 0
### function prop.test ####
name <- raw[,1]%>%as.matrix(.)
data <- raw[,2:3]%>%as.matrix(.)
sum <- colSums(data)
res <- c()
for (i in 1:nrow(data)){
  ct <- prop.test(unlist(data[i,]),colSums(data))
  res <- c(res,name[i],(data[i,]),ct$statistic,ct$p.value,ct$estimate)
}
res <- matrix(res,byrow = T,ncol = 7)%>%data.frame()
colnames(res) <- c("rename","value1","value2","statistic","p","Prop1","Prop2")
res$p.adj <- p.adjust(res$p,method = "fdr")
res$Prop1<-as.numeric(res$Prop1) ## mean rate
res$Prop2<-as.numeric(res$Prop2) ## mean rate
###
write.table(res,paste0(outdir,prefix,".CST.prop.test.txt"),quote = F,sep = "\t")

#### plot #####
#res <- read.table(paste0(outdir,"/cluster/AmericaChina.CST.prop.test.txt"),header=T)
raw.data <- res #%>%subset(p.adj<0.05)
## add 0.1*min.nonzero
min.nonzero <-setdiff(c(raw.data$Prop1,raw.data$Prop2),0)%>%min()
raw.data$China.prop <- ifelse(raw.data$Prop1==0,0.1*min.nonzero,as.vector(raw.data$Prop1))
raw.data$America.prop <- ifelse(raw.data$Prop2==0,0.1*min.nonzero,as.vector(raw.data$Prop2))
raw.data$mag.fc <- raw.data$China.prop/raw.data$America.prop
## labels
raw.data$China.label <- paste0(round(raw.data$Prop1,4)*100,"%") 
raw.data$America.label <- paste0(round(raw.data$Prop2,4)*100,"%") 
##
data <- melt(raw.data,id.var="rename")%>%subset(variable%in%c("America.prop","China.prop"))%>%mutate(value=as.numeric(value))
data$labels <- paste0(round(data$value,4)*100,"%") 

#### order by avg.fc############
order.rename <- raw.data$rename[order(raw.data$America.prop)]%>%as.vector()
data$rename <- factor(data$rename,levels = order.rename)

p1 <- ggplot(data,aes(x=value,y=rename,fill=variable,label=labels))+
  geom_col(position = "dodge")+
  geom_text(size = 2,vjust = -1, position = position_dodge(.9),angle=-90,color = "grey40")+
  scale_fill_manual(values = c("#FDBF6F","#1F78B4"))+
  xlim(0,0.6)+
  labs(x="Vagitype Rate",y="",fill="")+theme_pubclean()

raw.data$rename <- factor(raw.data$rename,levels = order.rename)
raw.data$temp <- as.factor(1)
p3<-ggplot(raw.data,aes(x=temp,y=rename,fill=-log10(p.adj)))+
  geom_tile(color="grey60")+coord_equal()+
  scale_fill_gradientn(colors = brewer.pal(8,"Greens")[3:8])+
  geom_text(aes(label=round(p.adj,2)))+
  labs(x="-Log10(P adjust)",y="",fill="-Log10(P adjust)")+theme_minimal()

ggarrange(p1,p3,align = "hv",legend = "bottom",ncol = 2,widths = c(1,1))
ggsave(paste0(outdir,"CST.prop.pdf"),width = 8,heigh=8)

##
values=c(brewer.pal(12,"Set3"),brewer.pal(12,"Paired"),brewer.pal(8,"Set2"),brewer.pal(8,"Set1"),rep("grey",40))
cst.avg.top3$clustering <-factor(cst.avg.top3$clustering,levels = order.rename)
p4 <- ggplot(cst.avg.top3,aes(x=mean,y=clustering,fill=taxa))+
  geom_col(position = "dodge")+
  xlim(0,0.8)+
  geom_text(aes(label=taxa),size=2,vjust = 0.3,hjust= 0,  position = position_dodge(.9))+
  scale_fill_manual(values = values)+theme_minimal()

ggarrange(p1,p3,p4,align = "hv",legend = "bottom",ncol = 3,widths = c(2,1,2))
ggsave(paste0(outdir,"CST.meantop.prop.pdf"),width = 12,heigh=8)

## shannon - sgb #####
alpha.stat <- function(profile){
  library(vegan)
  dataT<-t(profile) ###rm the tax &tax.all 
  shannon.wiener=data.frame(values=diversity(dataT, "shannon"))
  simpson=data.frame(values=diversity(dataT,"simpson"))#????simpsonָ??
  s=data.frame(values=specnumber(dataT)) #?????????ۼ???S  richness
  pielou=data.frame(values=shannon.wiener/log(s)) #????Pielou???ȶ?ָ??
  ##
  alpha<-data.frame(Richness=s$values,
                    Shannon=shannon.wiener$values,
                    Pielou=pielou$values,
                    SeqID=rownames(dataT))
  return(alpha)
}
alpha.sgb <- alpha.stat(prof.scaled)
write.table(alpha.sgb,paste0(outdir,"/alpha.sgb.txt"),quote = F,row.names = F,sep = "\t")
##
col.cluster <- c("#901C42","#9F267E","#D42C7F","#BF8EBE","#F98400","#FF0000","#74A9CF" ,"#045A8D","#A7D38B","#44A15A","#006536", "#EDEB64")
names(col.cluster) <- c("CVT_G4","CVT_G6","CVT_G1","CVT_G8","CVT_G10","CVT_G7","CVT_L3","CVT_L2","CVT_G3","CVT_G5","CVT_G2","CVT_G9" )
order.rename <- c( "CVT_G10","CVT_G5",  "CVT_G4",  "CVT_G2",  "CVT_G1",  "CVT_G7",  "CVT_G8",
                   "CVT_G3",  "CVT_G6",  "CVT_G9",  "CVT_L3",  "CVT_L2" )
##
cst.alpha.sgb <- merge(cst.stat%>%rename(SeqID=Samp),
                       melt(alpha.sgb,id.vars = "SeqID",variable.name = "alpha"),by="SeqID",all = T)
cst.alpha.sgb$clustering <- factor(cst.alpha.sgb$clustering,levels = order.rename)
cst.alpha.sgb%>%subset(alpha!="Pielou")%>%
  ggplot(aes(clustering,value,color=clustering))+coord_flip()+
  geom_boxplot(aes(fill=clustering),alpha =0.6)+
  scale_fill_manual(values = col.cluster)+
  scale_color_manual(values = col.cluster)+
  facet_wrap(~alpha,scales = "free")
ggsave(paste0(outdir,"CST.alpha.pdf"),width = 10,heigh=8)
