#haolilan
library(ggplot2)
library(reshape2)
library(dplyr)
library(readxl)
library(tibble)
library(ggrepel)
library(stringr)
library(pheatmap)
library(ggpubr)
library(viridis)
library(RColorBrewer)
options(ggplot2.continuous.colour="viridis")
options(ggplot2.continuous.fill = "viridis")
options(
  ggplot2.discrete.colour = c(brewer.pal(8,"Set2"),brewer.pal(12,"Paired")),
  ggplot2.discrete.fill = c(brewer.pal(8,"Set2"),brewer.pal(12,"Paired"))
)
theme_set(theme_bw()+theme(panel.grid = element_blank()))
###
outdir   <- paste0(getwd(),"/")

### input ####
prof.votu <- read.delim("vOTU_relab_v2025.txt")
prof.votu <- prof.votu[grep("t__",prof.votu$clade_name,invert = T),]
rownames(prof.votu) <- NULL
prof.votu <- prof.votu%>%column_to_rownames("clade_name")
rownames(prof.votu) <- rownames(prof.votu)%>%gsub("^.*__","",.)

##### vOTUInfo,CC_vOTUInfo #######################
vOTUInfo <- read_excel("info_20250203.xlsx",sheet = "vOTU")
HpvInfo <- read_excel("info_20250203.xlsx",sheet = "hpv_type")
HpvInfo <- HpvInfo %>%  mutate(vOTU=vOTUInfo$vOTU[match(vOTU_ID,vOTUInfo$vOTU_ID)])
id.hrHpv <- paste0("Human papillomavirus ",c(52,16,58,56,51,66,39,31,18,59,33,45,35))
HpvInfo$`High-risk`<- ifelse(HpvInfo$types%in%id.hrHpv,"High-risk",NA)

vOTUInfo <- vOTUInfo %>% dplyr::select(taxa_prof,vOTU,`Number of viral sequences`,`Taxonomic classification (family level)`,`Host assignment (phylum level)` :Lifestyle,vOTU_ID,virus_type )
vOTUInfo$Host_if <- ifelse(vOTUInfo$`Host assignment (phylum level)`!="-","MircobeHost","NonMircobeHost")
write.table(vOTUInfo,"reshape.vOTUInfo.txt",sep = "\t",quote = F,row.names = F)
write.table(HpvInfo,"reshape.HpvInfo.txt",sep = "\t",quote = F,row.names = F)
##
vOTUInfo$Family_level <- vOTUInfo$`Taxonomic classification (family level)` ## protien
vOTUInfo$Family_level[vOTUInfo$Family_level=="-"] <- "Unclassified"
vOTUInfo$HPV <- HpvInfo$species[match(vOTUInfo$vOTU,HpvInfo$vOTU)]
vOTUInfo$`HPV type` <- HpvInfo$types[match(vOTUInfo$vOTU,HpvInfo$vOTU)]
vOTUInfo <- vOTUInfo %>%
  mutate(Genome_sampleID = gsub("_NODE.*$|_k.*$","",vOTU_ID)) #  source_sampleID == Genome_sampleID
##
re_SampleInfo <- read.delim("reshape.CC_SampleInfo.txt",header = T) 
setdiff(vOTUInfo$Genome_sampleID, re_SampleInfo$BioSample.ID)
intersect(vOTUInfo$Genome_sampleID, re_SampleInfo$BioSample.ID)
vOTUInfo <- vOTUInfo %>%
  mutate(Country = re_SampleInfo$Country[match(Genome_sampleID,re_SampleInfo$BioSample.ID)],
         Source  = re_SampleInfo$Source[match(Genome_sampleID,re_SampleInfo$BioSample.ID)]
  )
dim(vOTUInfo)
table(vOTUInfo$Country)%>%addmargins()

## Genome size ##################
vOTU.contig <- read_excel("info_20250203.xlsx",sheet = "ALL",skip = 2)
vOTU.contig$GenomeSize  <- vOTU.contig$...3

my_comparisons <- list(c("Medium-quality","Complete"),c("Medium-quality","High-quality" ),c("Complete","High-quality"))
ggplot(vOTU.contig,aes(x=`Genome quality`,y=GenomeSize/1000,fill=`Genome quality`))+
  geom_boxplot()+
  labs(y="Genome size (Kbp)",x="Genome quality")+
  stat_compare_means(comparisons=my_comparisons,method = "wilcox.test")
ggsave("vOTU Genome quality.pdf",width = 8,height = 6)

### HPV types #############
HpvInfo$types[HpvInfo$types=="?"]<-"New Type"
Hpvtypes <- HpvInfo%>%group_by(types,species)%>%summarise(Num_virus_seq=sum(`Number of viral sequences`),Num_vOTU=n())%>%arrange(-Num_virus_seq)
Hpvtypes$Num_virus_seq_rate <- prop.table(Hpvtypes$Num_virus_seq)%>%round(.,4)*100
Hpvtypes$types <- factor(Hpvtypes$types,levels = rev(unique(Hpvtypes$types)))
Hpvtypes$hr.risk <- HpvInfo$`High-risk`[match(Hpvtypes$types,HpvInfo$types)]
write.table(Hpvtypes,"HPV_type_contigs_summary.txt",quote = F, sep = "\t")
## 或者着色换成高风险和低风险。
ggplot(Hpvtypes,aes(types,Num_virus_seq,fill=hr.risk))+geom_col(width = 0.8)+
  geom_point(aes(size=Num_vOTU,color=hr.risk,y=-10))+
  coord_flip()+scale_size(range = c(2,4))
ggsave("HPV_type_contigs_summary.pdf",width = 6, height = 8)

## Host info ##################
vOTUInfo_host <- vOTUInfo%>%subset(`Host assignment (phylum level)`!="-" )
vOTUInfo_host$Num_HostPhylum <- vOTUInfo_host$`Host assignment (phylum level)`%>%str_count(.,pattern=",")+1
vOTUInfo_host$Num_HostSpecies <- vOTUInfo_host$`Host assignment (species level)`%>%str_count(.,pattern=",")+1
vOTUInfo_host <- vOTUInfo_host%>%select(vOTU,Num_HostPhylum,Num_HostSpecies, `Number of viral sequences`:virus_type)

## host phylum
temp <- str_split_fixed(vOTUInfo_host$`Host assignment (phylum level)`,", ",max(vOTUInfo_host$Num_HostPhylum))
vOTU.P <- cbind(vOTUInfo_host[,c("vOTU","Number of viral sequences")],temp)%>%melt(.,value.name = "Host.phylum",id.vars=c("vOTU", "Number of viral sequences"))%>%
  subset(Host.phylum!="")%>%select(-variable)%>%.[!duplicated(.),]
xtabs(~Host.phylum,vOTU.P)%>%data.frame()%>%arrange(-Freq)
vOTU.P$Host.phylum%>%unique()%>%length()

## host species
temp <- str_split_fixed(vOTUInfo_host$`Host assignment (species level)`,", ",max(vOTUInfo_host$Num_HostSpecies))
vOTU.S <- cbind(vOTUInfo_host[,c("vOTU","Number of viral sequences")],temp)%>%melt(.,value.name = "Host.species",id.vars=c("vOTU", "Number of viral sequences"))%>%
  subset(Host.species!="")%>%select(-variable)%>%.[!duplicated(.),]
vOTU.S$Host.species%>%unique()%>%length()

## host genus
vOTU.S$Host.genus <- gsub("_unknown","",vOTU.S$Host.species)%>%gsub(" .*$","",.)
taxa.lac <- vOTU.S$Host.genus[grep("actobacillus",vOTU.S$Host.genus)]%>%unique()
taxa.lac
vOTU.S$Host.genus <- ifelse(vOTU.S$Host.genus%in% taxa.lac, "Lactobacillus",as.vector(vOTU.S$Host.genus))
vOTU.G <- vOTU.S%>%group_by(Host.genus,vOTU)%>%summarise(N=n())%>%arrange(-N)
xtabs(~Host.genus,vOTU.G)%>%data.frame()%>%arrange(-Freq)%>%head()
xtabs(~Host.genus,vOTU.G)%>%data.frame()%>%arrange(-Freq)%>%head()%>%ggplot(.,aes(x=Host.genus,y=Freq))+geom_col(fill="#A6CEE3")
##
id.votu.bifi <- vOTU.G$vOTU[vOTU.G$Host.genus=="Bifidobacterium"]
id.votu.lac <- vOTU.G$vOTU[vOTU.G$Host.genus=="Lactobacillus"]

#### multiple host species #########
addmargins(xtabs(~Num_HostSpecies,vOTUInfo_host%>%subset(vOTU%in%id.votu.bifi)))
addmargins(xtabs(~Num_HostSpecies,vOTUInfo_host%>%subset(vOTU%in%id.votu.lac)))

#### intra Bifidobacterium vaginales ########
vOTU.bifi <- vOTU.S%>%subset(Host.genus%in% "Bifidobacterium")
vOTU.bifi$Host.species%>%unique()
bv <- setdiff(vOTU.bifi$Host.species%>%unique(),c("Bifidobacterium breve" ,"Bifidobacterium dentium",     "Bifidobacterium longum", "Bifidobacterium adolescentis","Bifidobacterium pseudocatenulatum",
                                            "Bifidobacterium bifidum","Bifidobacterium scardovii","Bifidobacterium_unknown"  ))
bv
vOTU.bv <- vOTU.S%>%subset(Host.species%in%bv)
xtabs(~vOTU,vOTU.bv)%>%table()%>%addmargins()

#### intra genome ########
addmargins(xtabs(~N,vOTU.G%>%subset(Host.genus=="Bifidobacterium")))
addmargins(xtabs(~N,vOTU.G%>%subset(Host.genus=="Lactobacillus")))
pb <- xtabs(~N,vOTU.G%>%subset(Host.genus=="Bifidobacterium"))%>%data.frame()%>%ggplot(.,aes(x=N,y=Freq,label=Freq))+geom_col(fill="#A6CEE3")+geom_text(vjust=0)
pl <- xtabs(~N,vOTU.G%>%subset(Host.genus=="Lactobacillus"))%>%data.frame()%>%ggplot(.,aes(x=N,y=Freq,label=Freq))+geom_col(fill="#FF7F00")+geom_text(vjust=0)
ggarrange(pl+labs(title = "Lactobacillus"),pb+labs(title = "Bifidobacterium"),widths = c(1,3))
ggsave("Host.intraspecies.stat.pdf",width = 10,height = 6)
#
votu.lac <- dcast(vOTU.S%>%subset(Host.genus=="Lactobacillus"),vOTU~Host.species)
votu.lac <- votu.lac %>% subset(vOTU%in%vOTU.G$vOTU[vOTU.G$N>1&vOTU.G$Host.genus=="Lactobacillus"])
rownames(votu.lac) <- NULL
votu.lac <- votu.lac%>%column_to_rownames("vOTU")
votu.lac[!is.na(votu.lac)]<-1
votu.lac[is.na(votu.lac)]<-0
votu.lac <- apply(votu.lac, c(1,2), as.numeric)
votu.lac <- votu.lac[rowSums(votu.lac)!=0,colSums(votu.lac)!=0]
ph.l <- pheatmap(votu.lac,color = colorRampPalette(c("#FFF5EB","#ff7f00"))(10),border_color = "#FFF5EB",main = "Lactobacillus N>1")

#
votu.bifi <- dcast(vOTU.S%>%subset(Host.genus=="Bifidobacterium"),vOTU~Host.species)
votu.bifi <- votu.bifi %>% subset(vOTU%in%vOTU.G$vOTU[vOTU.G$N>2&vOTU.G$Host.genus=="Bifidobacterium"])
rownames(votu.bifi) <- NULL
votu.bifi <- votu.bifi%>%column_to_rownames("vOTU")
votu.bifi[!is.na(votu.bifi)]<-1
votu.bifi[is.na(votu.bifi)]<-0
votu.bifi <- apply(votu.bifi, c(1,2), as.numeric)
votu.bifi <- votu.bifi[rowSums(votu.bifi)!=0,colSums(votu.bifi)!=0]
ph.b <- pheatmap(votu.bifi,color = colorRampPalette(c("#FFF5EB","#ff7f00"))(10),border_color = "#FFF5EB",main = "Bifidobacterium N>2")

ggarrange(ph.l$gtable,ph.b$gtable,widths = c(1.25,2),align = "hv")
ggsave("Host.intraspecies.hp.pdf",width = 12,height = 6)

## occurrence  #########
#### whole #####
prefix <- "whole_"
prof.input <- prof.votu
tax.stat  <- data.frame(
  avg = apply(prof.input, 1, mean),
  sd = apply(prof.input, 1, sd),
  occurence = apply(prof.input, 1, function(x)length(x[x>0.0001])/length(x)),
  maxtaxa = apply(prof.input, 1, max))%>%
  arrange(avg)
tax.stat$vOTU <- rownames(tax.stat)
tax.stat$taxa <- vOTUInfo$taxa_prof[match(tax.stat$vOTU,vOTUInfo$vOTU)]
tax.stat$phylum <- tax.stat$taxa %>% gsub("^.*p__","",.)%>%gsub(".c__.*$","",.)
# tax.stat$family <- tax.stat$taxa %>% gsub("^.*f__","",.)%>%gsub(".g__.*$","",.)
tax.stat$family <- vOTUInfo$Family_level[match(tax.stat$vOTU,vOTUInfo$vOTU)]
tax.stat$HPVtype <- vOTUInfo$`HPV type`[match(tax.stat$vOTU,vOTUInfo$vOTU)]%>%gsub("Human papillomavirus","HPV",.)
tax.stat$host.phylum  <- vOTUInfo$`Host assignment (phylum level)`[match(tax.stat$vOTU,vOTUInfo$vOTU)]
tax.stat$host.species <- vOTUInfo$`Host assignment (species level)`[match(tax.stat$vOTU,vOTUInfo$vOTU)]
tax.stat <- tax.stat %>% select(-taxa)
write.table(tax.stat,paste0(outdir,prefix,"vOTU.profile.S.average.txt"),quote = F,sep = "\t",row.names = F)
### 
ggplot(tax.stat,aes(log10(occurence)))+geom_histogram()
top <- tax.stat%>%subset(occurence>0.05)%>%arrange(-occurence)%>%.[,"vOTU"]
p1 <- ggplot(tax.stat%>%subset(occurence>0.05),aes(occurence,vOTU,fill=family))+geom_col(alpha=0.8)+scale_y_discrete(limits=top)+labs(subtitle = "occurence(rel_ab>10e-4)\n ")
p2 <- ggplot(tax.stat%>%subset(occurence>0.05),aes(occurence,vOTU,fill=host.phylum,label=host.species))+geom_col(alpha=0.8)+scale_y_discrete(limits=top)+    labs(subtitle = "occurence(rel_ab>10e-4)\n ")+geom_text(size = 2,vjust = 0,hjust=1, position = position_dodge(.9),angle=0,color = "grey20")
ggarrange(p1,p2,widths=c(1,2),nrow = 1,align = "hv",legend = "bottom")
ggsave(paste0(outdir,prefix,"vOTU_top_occurrence_0.05.pdf"),width = 16,height = 8)
 
top <- tax.stat%>%subset(occurence>0.05)%>%arrange(-occurence)%>%.[,"vOTU"]
tax.stat$vOTU.label <- gsub("^.*s__","",tax.stat$vOTU)
ggplot(tax.stat%>%subset(vOTU%in%top),aes(occurence,log10(avg/100),color=family,size=maxtaxa))+geom_point(alpha=0.8)+
  geom_text_repel(aes(label=vOTU.label),size=3)+
  labs(title = "Top of mean vOTU rel_ab")
ggsave(paste0(outdir,prefix,"vOTU_top_occurrence_avg.pdf"),width = 8,height = 8)


#### area #########
phen.input <- re_SampleInfo%>%
  subset(Prof_ID%in% colnames(prof.votu))%>%
  subset(Country %in% c("USA","China"))%>%
  mutate(group = Country)
table(phen.input$group)

for (a in unique(phen.input$group)){
  prefix <- paste0(a,"_")
  prof.input <- prof.votu %>%.[,phen.input$Prof_ID[phen.input$group==a]]%>%.[rowSums(.)!=0,colSums(.)!=0]
  tax.stat  <- data.frame(
    avg = apply(prof.input, 1, mean),
    sd = apply(prof.input, 1, sd),
    occurence = apply(prof.input, 1, function(x)length(x[x>0.01])/length(x)),
    maxtaxa = apply(prof.input, 1, max))%>%
    arrange(avg)
  tax.stat$vOTU <- rownames(tax.stat)
  tax.stat$taxa <- vOTUInfo$taxa_prof[match(tax.stat$vOTU,vOTUInfo$vOTU)]
  tax.stat$phylum <- tax.stat$taxa %>% gsub("^.*p__","",.)%>%gsub(".c__.*$","",.)
  tax.stat$family <- tax.stat$taxa %>% gsub("^.*f__","",.)%>%gsub(".g__.*$","",.)
  tax.stat$HPVtype <- vOTUInfo$`HPV type`[match(tax.stat$vOTU,vOTUInfo$vOTU)]%>%gsub("Human papillomavirus","HPV",.)
  tax.stat$host.phylum  <- vOTUInfo$`Host assignment (phylum level)`[match(tax.stat$vOTU,vOTUInfo$vOTU)]
  tax.stat$host.species <- vOTUInfo$`Host assignment (species level)`[match(tax.stat$vOTU,vOTUInfo$vOTU)]
  tax.stat <- tax.stat %>% select(-taxa)
  write.table(tax.stat,paste0(outdir,prefix,"vOTU.profile.S.average.txt"),quote = F,sep = "\t",row.names = F)
  ### 
  ggplot(tax.stat,aes(log10(occurence)))+geom_histogram()
  top <- tax.stat%>%subset(occurence>0.1)%>%arrange(-occurence)%>%.[,"vOTU"]
  p1 <- ggplot(tax.stat%>%subset(occurence>0.1),aes(occurence,vOTU,fill=family))+geom_col(alpha=0.8)+scale_y_discrete(limits=top)+    labs(subtitle = "occurence(rel_ab>10e-4)\n ")
  p2 <- ggplot(tax.stat%>%subset(occurence>0.1),aes(occurence,vOTU,fill=host.phylum,label=host.species))+geom_col(alpha=0.8)+scale_y_discrete(limits=top)+    labs(subtitle = "occurence(rel_ab>10e-4)\n ")+geom_text(size = 2,vjust = 0,hjust=1, position = position_dodge(.9),angle=0,color = "grey20")
  ggarrange(p1,p2,widths=c(1,2),nrow = 1,align = "hv",legend = "bottom")
  ggsave(paste0(outdir,prefix,"vOTU_top_occurrence.pdf"),width = 16,height = 8)
  ### Top 40
  top <- tax.stat%>%subset(occurence>0.1)%>%arrange(-avg)%>%.[,"vOTU"]
  tax.stat$vOTU.label <- gsub("^.*s__","",tax.stat$vOTU)
  ggplot(tax.stat%>%subset(vOTU%in%top),aes(occurence,log10(avg/100),color=family,size=maxtaxa))+geom_point(alpha=0.8)+
    geom_text_repel(aes(label=vOTU.label),size=3)+
    labs(title = "Top of mean vOTU rel_ab")
  ggsave(paste0(outdir,prefix,"vOTU_top_occurrence_avg.pdf"),width = 8,height = 8)
}


## SU-CCS 2018/2019/2021  ############
lm_2factors <- function(profile,phen,prefix.file){
  
  id.same <- intersect(colnames(profile),rownames(phen))
  profile <- profile%>% .[,id.same]%>%.[rowSums(.)!=0,colSums(.)!=0]
  phen    <- phen[id.same,] 
  
  identical(colnames(profile),rownames(phen))
  prof <- t(profile)
  ##
  p_species <- data.frame(1:ncol(prof))
  rownames(p_species) <- colnames(prof)
  p_species$names <- colnames(prof)
  
  coef_species <- data.frame(1:ncol(prof))
  rownames(coef_species) <- colnames(prof)
  coef_species$names <- colnames(prof)

  for (i in 1:length(phen)){
    groups <- phen[,i]
    prof.input <- prof[!is.na(groups),]%>%.[rowSums(.)!=0,colSums(.)!=0]
    
    groups <- groups[!is.na(groups)]
    groups <- as.factor(groups)
    
    pseudocount_overall <- min(prof.input[prof.input!=0], na.rm=T) *0.01
    
    temp <- data.frame(1:ncol(prof.input))
    rownames(temp) <- colnames(prof.input)
    temp$names <- colnames(prof.input)
    
    coef_temp <- data.frame(1:ncol(prof.input))
    rownames(coef_temp) <- colnames(prof.input)
    coef_temp$names <- colnames(prof.input)

    for (j in 1:ncol(prof.input)) {
      gmm <-
        glm( scale(log(prof.input[, j]+pseudocount_overall))  ~ groups)

      if (nrow(summary(gmm)$coefficient) < 2) {
        temp[j, 1] <- NA
        coef_temp[j, 1] <- NA
      } else {
        temp[j, 1] <- summary(gmm)$coefficient[2, 4]
        coef_temp[j, 1] <- summary(gmm)$coefficient[2, 1]
      }
    }
    colnames(temp) <- c(colnames(phen)[i], "names"  )
    colnames(coef_temp) <- c(colnames(phen)[i], "names"  )
    p_species <- merge(p_species, temp,by="names",all=T)
    coef_species <- merge(coef_species, coef_temp,by="names",all=T)
  }
  
  p_species <- p_species%>%column_to_rownames("names")
  p_species <- p_species[, -1]
  colnames(p_species) <- colnames(phen)
  
  coef_species <- coef_species%>%column_to_rownames("names")
  coef_species <- coef_species[, -1]
  colnames(coef_species) <-  colnames(phen)

  rm(gmm, temp, coef_temp)
  
  #glm.out <- data.frame(coef = unlist(coef_species[1,]),p.value = unlist(p_species[1,]),p.adj = p.adjust(p_species[1,], method = "fdr"))
  glm.out <- merge(coef_species%>%rownames_to_column("names")%>%melt(variable.name = "factor",value.name = "coefficient"),
                   p_species%>%rownames_to_column("names")%>%melt(variable.name = "factor",value.name = "pvalue"),
                   by=c("names","factor"),all=T)
  glm.out$p.adj <- p.adjust(glm.out$pvalue, method = "fdr")
  glm.out.0.05 <- glm.out %>% subset(pvalue < 0.05)
  
  write.csv(glm.out,paste0(outdir,"/lm.out.",prefix.file,".csv"),quote = F,row.names = T)
  write.csv(glm.out.0.05,paste0(outdir,"/lm.out.p.",prefix.file,".csv"),quote = F,row.names = T)
  return(glm.out)
}
##
prefix.file <- "SU-CCS"
phen.SEL<- phen.merge.profile %>%subset(Source%in%c("SU-CCS2018/2019 (Peacock)","SU-CCS2021 (Peacock)" ))%>%
  dplyr::select(SeqID,Menopause:VVC)
## 只保留数值型变量
ID <- phen.SEL[,"SeqID"]
phen.SEL <- phen.SEL[,sapply(phen.SEL, function(x)is.numeric(x))]
rownames(phen.SEL) <- NULL
rownames(phen.SEL) <- ID
##
phen.raw <- phen.SEL
profile <- prof.votu%>%.[apply(.,1,function(x)sum(x>0)>50),]%>%.[rowSums(.)!=0,colSums(.)!=0]
id.same <- intersect(colnames(profile),rownames(phen.raw))
profile <- profile%>% .[,id.same]%>%.[rowSums(.)!=0,colSums(.)!=0]
phen    <- phen.raw[id.same,]
identical(colnames(profile),rownames(phen))
res.lm<-lm_2factors(profile,phen,prefix.file)
res.lm.INFO<- merge(res.lm,vOTUInfo%>%select(taxa_prof:`HPV type`)%>%rename(names=vOTU),all.x = T)
write.table(res.lm.INFO,paste0(outdir,"/lm.out.p.Suzhou.vOTUInfo.txt"),quote = F,row.names = T,sep = "\t")

####plots##################
data.m <- merge(melt((profile/100)%>%rownames_to_column("names"),variable.name = "SeqID"),phen%>%rownames_to_column("SeqID"),all.x=T)
data.mean <- rbind(
          data.m %>% group_by(Menopause,names)%>%summarise(mean=mean(value))%>%
            dcast(.,names~Menopause,value.var = "mean")%>% mutate(factor="Menopause")%>%select(names,`0`,`1`,factor),
          data.m %>% group_by(HPV_test,names)%>%summarise(mean=mean(value))%>%
            dcast(.,names~HPV_test,value.var = "mean")%>% mutate(factor="HPV_test")%>%select(names,`0`,`1`,factor),
          data.m %>% group_by(BV,names)%>%summarise(mean=mean(value))%>%
            dcast(.,names~BV,value.var = "mean")%>% mutate(factor="BV")%>%select(names,`0`,`1`,factor))

data.draw<- merge(res.lm.INFO %>%subset(factor%in%c("BV","HPV_test","Menopause")),
                  data.mean,by=c("names","factor"),all.x = T)
pseudocount_overall <- 0.1*min(data.draw$`1`[data.draw$`1`!=0],data.draw$`0`[data.draw$`0`!=0])
pseudocount_overall
data.draw$log2fc <- log2((data.draw$`1`+pseudocount_overall)/(data.draw$`0`+pseudocount_overall))
data.draw$Host.genus <- data.draw$`Host assignment (species level)`%>%gsub(" .*$","",.)%>%gsub("_unknown","",.)
data.draw$Host.genus[grep(",",data.draw$Host.genus )] <- "Multiple"
data.draw$virus_type[is.na(data.draw$virus_type)] <- "Unclassified Family"
data.draw$factor <- factor(data.draw$factor,levels = c("HPV_test","Menopause","BV"))
##
data.draw$delabel <- NA
id.genus <- xtabs(~Host.genus,data.draw)%>%data.frame()%>%arrange(-Freq)%>%.[1:7,1]
id.genus <- c(as.vector(id.genus),"Metamycoplasma")
data.draw$delabel[!(data.draw$Host.genus%in%id.genus)]  <- "Other_MicrobesHost"
data.draw$delabel[data.draw$Host.genus%in%id.genus] <- data.draw$Host.genus[data.draw$Host.genus%in%id.genus]
data.draw$delabel[!is.na(data.draw$`HPV type`)] <- data.draw$`HPV type`[!is.na(data.draw$`HPV type`)] %>% gsub("Human papillomavirus .*$","HPV",.)
data.draw$delabel[data.draw$delabel=="-"]  <- "Unssigned_Host"

data.draw$delabel
data.draw$delabel <- factor(data.draw$delabel,levels = c("Bifidobacterium","Prevotella","Lactobacillus",
                                                         "Metamycoplasma","Fannyhessea","Dialister",       
                                                        "Multiple","Unssigned_Host","Other_MicrobesHost","HPV" ))
## coef
ggplot(data=data.draw%>%subset(p.adj<0.05), aes(x=coefficient, y=-log10(p.adj), col=delabel)) + #,shape=virus_type
  geom_point(aes(size=abs(coefficient)),position = position_jitterdodge(),size=1,alpha=0.8) +   
  geom_hline(yintercept = -log10(0.05),linetype=2,color="red")+
  geom_vline(xintercept = 0,linetype=2,color="grey30")+
  theme_light()+facet_wrap(~factor,scales = "free_y")+  
  theme(panel.grid = element_blank())+
  geom_text_repel(aes(label=names),size=1.5)
ggsave(paste0(outdir,"lm.out.coef.p0.05.pdf"),width = 14,height = 7)
# 
# ## 
# data.draw$abs_coef <- ifelse(data.draw$coefficient>0,"1","-1")%>%as.factor()
# temp <- data.draw%>%subset(p.adj<0.05)%>%group_by(factor,abs_coef,delabel)%>%summarise(N=n())
# temp$N<- ifelse(temp$abs_coef=="-1",-1*as.vector(temp$N),as.vector(temp$N))
# temp$delabel <- factor(temp$delabel,levels = c("Bifidobacterium","Prevotella","Lactobacillus",
#                                                          "Metamycoplasma","Fannyhessea","Dialister",       
#                                                          "Multiple","Unssigned_Host","Other_MicrobesHost","HPV" ))
# #temp <- temp %>%subset(!(delabel%in%c("Unssigned_Host","Other_MicrobesHost")))
# ggplot(data=temp, aes(x=N, y=delabel, fill=delabel)) + #,shape=virus_type
#   geom_col() +  xlim(-max(abs(temp$N)),max(abs(temp$N)))+
#   geom_vline(xintercept = 0,linetype=1,color="grey40")+
#   theme_light()+facet_wrap(~factor,scales = "fixed")+  
#   theme(panel.grid = element_blank())
# ggsave(paste0(outdir,"lm.out.coef.stat.p0.05.pdf"),width = 18,height = 4)




## Correlation between SGB & vOTUs #################
outdir <- paste0(getwd(),"/")
dis.lst <- list.dirs(path = paste0(outdir,"/sparcc/output_202502"))[-1]
dis.lst
##
SGBInfo <- read.delim("reshape.SGBProfInfo.n759.txt",header = T)
taxa.info <- data.frame(names=c(SGBInfo$taxa_prof,vOTUInfo$vOTU),taxa=c(SGBInfo$taxa,vOTUInfo$taxa_prof),rename = c(SGBInfo$rename,vOTUInfo$vOTU))
taxa.info$domain <- taxa.info$taxa%>%gsub("^d__","",.)%>%gsub(";.*$","",.)
taxa.info$domain[taxa.info$domain == "Unclassified"] <- "Viruses" 
taxa.info$phylum <- taxa.info$taxa%>%gsub("^.*p__","",.)%>%gsub("d__Fungi;p_","",.)%>%gsub(";.*$","",.)
taxa.info$family <- taxa.info$taxa%>%gsub("^.*f__","",.)%>%gsub(";.*$","",.)
taxa.info$genus <- taxa.info$taxa%>%gsub("^.*g__","",.)%>%gsub(";.*$","",.)
##
#taxa.info$rename <- taxa.info$names
#taxa.info$rename[str_detect(taxa.info$rename,"^SGB")] <- paste0(taxa.info$genus,"_",taxa.info$names)[str_detect(taxa.info$rename,"^SGB")] 
taxa.info$rename[str_detect(taxa.info$rename,"^vOTU")] <- paste0(taxa.info$family,"_",taxa.info$names)[str_detect(taxa.info$rename,"^vOTU")] 
#
taxa.lst<- data.frame(ID=taxa.info$names,genus=taxa.info$genus,rename=taxa.info$rename)%>%
  .[!duplicated(.),]
id.genus <- c("Bifidobacterium","Prevotella","Lactobacillus","Dialister" )
id.species <- taxa.lst$ID[taxa.lst$genus%in%id.genus]

##
taxa.votu <- vOTUInfo%>%select(vOTU:`Host assignment (species level)`,`HPV type`)
taxa.votu$Host.genus <- taxa.votu$`Host assignment (species level)`%>%gsub(" .*$","",.)%>%gsub("_unknown","",.)
taxa.votu$delabel <- taxa.votu$Host.genus
taxa.votu$delabel[grep(",",taxa.votu$delabel )] <- "Multiple"
taxa.votu$delabel[!(taxa.votu$Host.genus%in%id.genus)]  <- "Other_MicrobesHost"
taxa.votu$delabel[taxa.votu$Host.genus=="-"]  <- "Unssigned_Host"
taxa.votu$delabel[!is.na(taxa.votu$`HPV type`)] <- taxa.votu$`HPV type`[!is.na(taxa.votu$`HPV type`)] %>% gsub("Human papillomavirus .*$","HPV",.)
taxa.votu$delabel[taxa.votu$delabel=="?"]<-"Unidentified"
taxa.votu$rename <- taxa.info$rename[match(taxa.votu$vOTU,taxa.info$names)]

for (i in c(5,6)){
  dir <- dis.lst[i] #1,5,6
  prefix <- gsub(".*/","",dir)#%>%gsub(".merge_filter","",.)
  cor_prof <- read.table(paste0(dir,"/cor_sparcc.out.txt"))#data.frame(cor$r) 
  p_prof   <- read.table(paste0(dir,"/pvals_two_sided.txt"))#data.frame(cor$P)
  ##
  cor_prof <- cor_prof[rownames(cor_prof)%in%id.species,grep("vOTU",colnames(cor_prof))]
  p_prof <- p_prof[rownames(p_prof)%in%id.species,grep("vOTU",colnames(p_prof))]
  ## any < 0.05
  p_prof <- p_prof[apply(p_prof,1,function(x)any(x<0.05)),apply(p_prof,2,function(x)any(x<0.05))]
  cor_prof <- cor_prof[rownames(p_prof),colnames(p_prof)]
  ## rename
  colnames(cor_prof) <- taxa.info$rename[match(colnames(cor_prof),taxa.info$names)]
  rownames(cor_prof) <- taxa.info$rename[match(rownames(cor_prof),taxa.info$names)]
  colnames(p_prof) <- taxa.info$rename[match(colnames(p_prof),taxa.info$names)]
  rownames(p_prof) <- taxa.info$rename[match(rownames(p_prof),taxa.info$names)]
  
  ###
  p_sign <- ifelse(p_prof<0.05,"*","")
  ## annotation
  anno_col <- data.frame(Host.genus = taxa.votu$delabel[match(colnames(cor_prof),taxa.votu$rename)] )
  rownames(anno_col) <- colnames(cor_prof)
  anno_row <- data.frame(Genus = taxa.lst$genus[match(rownames(cor_prof),taxa.lst$rename)] )
  rownames(anno_row) <-  rownames(cor_prof)
  #
  anno_colors = list(Genus = c("#66C2A5", "#FC8D62", "#8DA0CB","#FFD92F"),
                     Host.genus = c("#66C2A5", "#FC8D62", "#8DA0CB","#FFD92F","#B3B3B3", "#A6CEE3","#1F78B4","grey30"))
  names(anno_colors$Genus) <- c("Bifidobacterium","Prevotella","Lactobacillus","Dialister" )
  names(anno_colors$Host.genus) <- c("Bifidobacterium","Prevotella","Lactobacillus","Dialister", "Unssigned_Host", "Other_MicrobesHost" ,"HPV", "Unidentified")
  
  ##
  paletteLength <- 50
  # myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
  myColor <- colorRampPalette(c("#998ec3", "white", "#f1a340"))(paletteLength)
  myBreaks <- c(seq(min(cor_prof), 0, length.out=ceiling(paletteLength/2) + 1), 
                seq(max(cor_prof)/paletteLength, max(cor_prof), length.out=floor(paletteLength/2)))
  
  # Plot the heatmap
  pdf(paste0(outdir,"Cohort_sgb_vOTU.",prefix,".pdf"),width = 12,height = 9)
  pheatmap(cor_prof, color=myColor, breaks=myBreaks,#cluster_rows = F,
           annotation_col = anno_col,
           annotation_row = anno_row,
           annotation_colors = anno_colors,
           display_numbers = p_sign,
           main = prefix)
  dev.off()
}
