#haolilan

library(readxl);library(dplyr);library(tibble);library(ggplot2);library(reshape2)
####
theme_set(theme_bw()+theme(panel.grid = element_blank()))

## DATA LOAD ######
prof.scaled <- read.delim("Whole.relative_taxonomic_abundance.recal.txt")
##
phen.merge.profile <- read.table("reshape.phen.merge.profile.txt",header = T,sep = "\t")
###

## BV diff #######
diff_ancom <- function(profile, groups, outdir, prefix.file,cutoff_feature=0.0001, filter_preval=0.01,filter_N=10){
  source("ancom.R")
  library(compositions);library(readxl);library(dplyr);library(tibble)
  ##
  groups     <- groups %>% subset(SeqID %in% colnames(profile))
  groups$V2  <- factor(groups$V2,levels = c("control","case"))
  prof.input <- profile %>% dplyr::select(groups$SeqID)%>% .[,groups$SeqID]%>% .[rowSums(.)!=0,]
  ###
  prevalence <- data.frame(x= apply(prof.input, 1, function(x){sum(x>cutoff_feature)/length(x)}),
                           sum= apply(prof.input, 1, function(x){sum(x>cutoff_feature)})) %>% subset(x > filter_preval & sum > filter_N)
  prof.input <- prof.input[rownames(prevalence),] #%>%.[rowSums(.)!=0,] ????????????
  # Check variables in the correct order
  identical(colnames(prof.input),groups$SeqID)
  dim(prof.input)
  levels(groups$V2)
  table(groups$V2)
  ## calculate mean by aroup
  tem  <- apply(prof.input,1,function(x)aggregate(unlist(x),by=list(groups$V2),mean))%>%data.frame()%>%.[,grep("Group.1",colnames(.),invert = T)]
  tem <- data.frame(t(tem))
  colnames(tem) <- levels(groups$V2)
  tem$enrich <- ifelse(tem$control>tem$case,"control",ifelse(tem$control<tem$case,"case","equal"))
  tem$taxa_id <- rownames(tem)%>%gsub("..$","",.)
  ##
  ancom_table <- feature_table_pre_process(prof.input, groups, "SeqID", group_var = NULL, out_cut = 0.05, zero_cut = 1, lib_cut = 0, neg_lb)
  saveRDS(ancom_table, paste(outdir, "/", prefix.file, "_ancom_table.rds", sep = ""))
  
  # Run tests with covariates (age + sex)
  ancom_out <- ANCOM(ancom_table$feature_table, ancom_table$meta_data, 
                     ancom_table$structure_zeros, main_var = "V2", 
                     p_adj_method = "BH",adj_formula=NULL) #adj_formula = "proforma_sex + age ")
  ancom_out$out <- merge(ancom_out$out,tem,by="taxa_id",all = T)
  ancom_out$out <- ancom_out$out[order(ancom_out$out$W, decreasing = T),]
  write.table(ancom_out$out, paste(outdir, "/", prefix.file, "_ANCOM.txt", sep = ""),
              col.names = T, row.names = F, sep = "\t", quote = F)
  saveRDS(ancom_out$fig, paste(outdir, "/", prefix.file, "_ANCOM.rds", sep = ""))
}

diff_prevalence <- function(profile, groups, outdir, prefix.file,cutoff_feature=0.0001, filter_preval=0.01,filter_N=10){
  library(compositions);library(readxl);library(dplyr);library(tibble)
  ##
  groups     <- groups %>% subset(SeqID %in% colnames(profile))
  groups$V2  <- factor(groups$V2,levels = c("control","case"))
  prof.input <- profile %>% dplyr::select(groups$SeqID)%>% .[,groups$SeqID]%>% .[rowSums(.)!=0,]
  ###
  prevalence <- data.frame(x= apply(prof.input, 1, function(x){sum(x>cutoff_feature)/length(x)}),
                           sum= apply(prof.input, 1, function(x){sum(x>cutoff_feature)})) %>% subset(x > filter_preval & sum > filter_N)
  prof.input <- prof.input[rownames(prevalence),] #%>%.[rowSums(.)!=0,] ????????????
  # Check variables in the correct order
  identical(colnames(prof.input),groups$SeqID)
  dim(prof.input)
  levels(groups$V2)
  table(groups$V2)
  ## calculate mean by aroup
  tem  <- apply(prof.input,1,function(x)aggregate(unlist(x),by=list(groups$V2),function(x)length(x[x>0])/length(x)))%>%data.frame()%>%.[,grep("Group.1",colnames(.),invert = T)]
  tem <- data.frame(t(tem))
  colnames(tem) <- levels(groups$V2)
  tem$enrich <- ifelse(tem$control>tem$case,"control",ifelse(tem$control<tem$case,"case","equal"))
  tem$taxa_id <- rownames(tem)%>%gsub("..$","",.)
  tem$NumControl <- table(groups$V2)[1]
  tem$NumCase <- table(groups$V2)[2]
  write.table(tem, paste(outdir, "/", prefix.file, "_prevalence.txt", sep = ""),
              col.names = T, row.names = F, sep = "\t", quote = F)
}

### bv xw ####
outdir     <- paste0(getwd(),"/data/")
groups     <- read.table("match.BV.xw.csv",header = T,sep = ",");colnames(groups)[2]<-"V2"
table(groups$V2)
## SGB 
profile <- prof.scaled   ### 0.0001
prefix.file  <- "BV.xw.sgb"
intersect(groups$SeqID,colnames(profile)) %>%length()
diff_prevalence(profile, groups, outdir, prefix.file,cutoff_feature=0.0001,filter_preval=0.01,filter_N=5)
diff_ancom(profile, groups, outdir, prefix.file,cutoff_feature=0.0001,filter_preval=0.01,filter_N=5)

### bv SU-CCS ####
outdir     <- paste0(getwd(),"/data/")
groups     <- read.table("match.BV.SU-CCS.csv",header = T,sep = ",");colnames(groups)[2]<-"V2"
table(groups$V2)
## SGB 
profile <- prof.scaled   ### 0.0001
prefix.file  <- "BV.SU-CCS.sgb"
intersect(groups$SeqID,colnames(profile)) %>%length()
diff_prevalence(profile, groups, outdir, prefix.file,cutoff_feature=0.0001,filter_preval=0.01,filter_N=5)
diff_ancom(profile, groups, outdir, prefix.file,cutoff_feature=0.0001,filter_preval=0.01,filter_N=5)

### bv BJ-GC ####
outdir     <- paste0(getwd(),"/data/")
groups     <- read.table("match.BV.BJ-GC.csv",header = T,sep = ",");colnames(groups)[2]<-"V2"
table(groups$V2)
## SGB 
profile <- prof.scaled   ### 0.0001
prefix.file  <- "BV.BJ-GC.sgb"
intersect(groups$SeqID,colnames(profile)) %>%length()
diff_prevalence(profile, groups, outdir, prefix.file,cutoff_feature=0.0001,filter_preval=0.01,filter_N=5)
diff_ancom(profile, groups, outdir, prefix.file,cutoff_feature=0.0001,filter_preval=0.01,filter_N=5)



## plot comparison multiple-cohort #####################################################################
SGBInfo <- read.delim("reshape.SGBProfInfo.n759.txt",header = T)
taxa.info <- data.frame(names=SGBInfo$taxa_prof,taxa=SGBInfo$taxa,rename = SGBInfo$rename)

###
outdir     <- ""
dis.lst <- dir(path = paste0(outdir),pattern=".sgb_ANCOM.txt")
dis.lst
##
data <- c()
for (i in 2:length(dis.lst)){
  dir <- dis.lst[i] #
  prefix <- gsub("_ANCOM.txt","",dir)#%>%gsub(".merge_filter","",.)
  temp <- read.table(paste0(outdir,dir),header = T)
  temp <- subset(temp,detected_0.7=="TRUE")
  #temp <- subset(temp,case>0.0001 | control >0.0001)
  temp <- dplyr::select(temp,-detected_0.9:-detected_0.6)
  temp$dataset <- prefix
  data <- rbind(data,temp)
}
stat <- data%>%group_by(taxa_id,enrich)%>%summarise(n=n())
stat <- dcast(stat,taxa_id~enrich,value.var = "n")%>%arrange(-case,-control)
write.table(stat,paste0(outdir,"BV.ancom.sign.stat.txt"),quote = F,sep = "\t",row.names = F)
write.table(data,paste0(outdir,"BV.ancom.sign.txt"),quote = F,sep = "\t",row.names = F)
### same direction #####
stat.draw <- stat[rowSums(stat[,2:3],na.rm = T)>1,]
stat.draw[is.na(stat.draw)] <- 0
stat.com <- stat.draw%>%subset(case==0 | control==0)%>%
  arrange(case,control)
###
data.m <- melt(data%>%subset(taxa_id%in%stat.com$taxa_id),id.var = c("taxa_id","enrich","W","dataset"), variable.name="group", value.name = "mean")
data.m <- arrange(data.m,dataset,enrich,W)
data.m$taxa_id <- factor(data.m$taxa_id, levels = unique(stat.com$taxa_id))
data.m$enrich <- factor(data.m$enrich,levels = c("case","control"))
data.m$group <- factor(data.m$group,levels = c("case","control"))
data.m$dataset <- factor(data.m$dataset ,levels = c( "BV.xw.sgb","BV.SU-CCS.sgb","BV.BJ-GC.sgb"))
data.m$rename <- taxa.info$rename[match(data.m$taxa_id,taxa.info$names)]
n.con <- nrow(stat.com[stat.com$control!=0,])

col.bv<- c("#6A3D9A","#CAB2D6" )
names(col.bv) <-c("case","control")
col.axis.x<-c()
col.axis.x[which(stat.com$control==stat.com$case)] <- "black"
col.axis.x[which(stat.com$control>stat.com$case)] <- "#CAB2D6"
col.axis.x[which(stat.com$control<stat.com$case)] <- "#6A3D9A"

ggplot(data.m,aes(x=log10(mean),y=rename))+
  geom_line()+
  geom_point(aes(color=group),alpha=0.8,position = position_dodge(0.01),size =4)+
  #geom_hline(yintercept = n.con +0.5, linetype = 2)+
  scale_fill_manual(values=col.bv)+
  scale_color_manual(values=col.bv)+
  labs(x="Log10(The Average Relative Abundance)",y="")+facet_grid(~dataset)+
  theme(axis.text.y = element_text(color = col.axis.x),
          panel.grid.major.y = element_line(color = "grey",linetype = 2))
ggsave(paste0(outdir,"BV.ancom.sign.com.pdf"),width = 8,height = 8)


### different direction #####
stat.diff <- stat.draw%>%subset(case!=0 & control!=0)%>%
  arrange(case,control)
###
data.m <- melt(data%>%subset(taxa_id%in%stat.diff$taxa_id),id.var = c("taxa_id","enrich","W","dataset"), variable.name="group", value.name = "mean")
data.m <- arrange(data.m,dataset,enrich,W)
data.m$taxa_id <- factor(data.m$taxa_id, levels = unique(stat.diff$taxa_id))
data.m$enrich <- factor(data.m$enrich,levels = c("case","control"))
data.m$group <- factor(data.m$group,levels = c("case","control"))
data.m$dataset <- factor(data.m$dataset ,levels = c( "BV.xw.sgb","BV.SU-CCS.sgb","BV.BJ-GC.sgb"))
data.m$rename <- taxa.info$rename[match(data.m$taxa_id,taxa.info$names)]

col.bv<- c("#6A3D9A","#CAB2D6" )
names(col.bv) <-c("case","control")

col.axis.x<-c()
col.axis.x[which(stat.diff$control==stat.diff$case)] <- "black"
col.axis.x[which(stat.diff$control>stat.diff$case)] <- "#CAB2D6"
col.axis.x[which(stat.diff$control<stat.diff$case)] <- "#6A3D9A"

ggplot(data.m,aes(x=log10(mean),y=rename))+
  geom_line()+
  geom_point(aes(color=group),alpha=0.8,position = position_dodge(0.01),size =4)+
  #geom_hline(yintercept = n.con +0.5, linetype = 2)+
  scale_fill_manual(values=col.bv)+
  scale_color_manual(values=col.bv)+
  labs(x="The Average Relative Abundance (log10)",y="")+facet_grid(~dataset)+
  theme(axis.text.y = element_text(color = col.axis.x),
        panel.grid.major.y = element_line(color = "grey",linetype = 2))
ggsave(paste0(outdir,"BV.ancom.sign.diff.pdf"),width = 8,height = 4)

## BV SGB cooccurrence ########
#### heatmap - only BV sample#############
SGBInfo <- read.delim("reshape.SGBProfInfo.n759.txt",header = T)
taxa.info <- data.frame(names=SGBInfo$taxa_prof,taxa=SGBInfo$taxa,rename = SGBInfo$rename)
#
  dir <- "output_202502/BV.all.BV.merge_filter"
  prefix <- gsub(".*/","",dir)%>%gsub(".BV.merge_filter","",.)
  cor_prof <- read.table(paste0(dir,"/cor_sparcc.out.txt"))#data.frame(cor$r) 
  p_prof   <- read.table(paste0(dir,"/pvals_two_sided.txt"))#data.frame(cor$P)
  #
  cor_prof <- cor_prof[grep("vOTU",rownames(cor_prof),invert = T),grep("vOTU",colnames(cor_prof),invert = T)]
  p_prof   <- p_prof[rownames(cor_prof),colnames(cor_prof)]
  
  temp<- cor_prof
  temp[abs(temp)<0.5]<-0
  temp <- temp[rowSums(temp)!=1,colSums(temp)!=1]
  cor_prof  <- cor_prof[rownames(temp),colnames(temp)]
  p_prof    <- p_prof[rownames(temp),colnames(temp)]
  
  ## any < 0.05
  p_prof <- p_prof[apply(p_prof,1,function(x)any(x<0.05)),apply(p_prof,2,function(x)any(x<0.05))]
  cor_prof <- cor_prof[rownames(p_prof),colnames(p_prof)]
  ## rename
  colnames(cor_prof) <- taxa.info$rename[match(colnames(cor_prof),taxa.info$names)]
  rownames(cor_prof) <- taxa.info$rename[match(rownames(cor_prof),taxa.info$names)]
  colnames(p_prof) <- taxa.info$rename[match(colnames(p_prof),taxa.info$names)]
  rownames(p_prof) <- taxa.info$rename[match(rownames(p_prof),taxa.info$names)]
  
  p_sign <- ifelse(p_prof<0.05,"*","")
  ##
  paletteLength <- 50
  myColor <- colorRampPalette(c("#998ec3", "white", "#f1a340"))(paletteLength)
  myBreaks <- c(seq(min(cor_prof), 0, length.out=ceiling(paletteLength/2) + 1), 
                seq(max(cor_prof)/paletteLength, max(cor_prof), length.out=floor(paletteLength/2)))
  
  # Plot the heatmap
  pdf(paste0(outdir,"BV_sgb",prefix,".pdf"),width = 12,height = 10)
  pheatmap(cor_prof, color=myColor, breaks=myBreaks,#cluster_rows = F,
           # annotation_col = anno_col,
           # annotation_row = anno_row,
           # annotation_colors = anno_colors,
           display_numbers = p_sign,
           main = prefix)
  dev.off()


## BV function cytolysin & sialiase####
prefix.data <- "KO"
data.input <- read.delim("KO_mags_count.tsv",header = T) 

## ko 
vf.ko <- read.delim("vf_v1.txt",header = T) 
vf.ko.list <- c(na.omit(vf.ko$ko[vf.ko$group=="Sialidase"]),na.omit(vf.ko$ko[vf.ko$group=="Cytolysin_/_Inerolysin_/_Vaginolysin_production"]))
vf.ko.list <- vf.ko.list[vf.ko.list%in%rownames(data.input)]

### all the genome ####      
## Genome_sampleID: 获取 Genome_ID 和 Genome_quality
MagInfo <- read.delim("reshape.MagInfo.txt",header = T)

MagInfo.input <- MagInfo ##不匹配人群为主
MagInfo.input$Genome_ID <- MagInfo.input$Genome_ID%>%gsub("-",".",.)
xtabs(~Genome_quality,MagInfo.input)

## 提取需要KO 和 genome
setdiff(colnames(data.input),MagInfo.input$Genome_ID)
setdiff(MagInfo.input$Genome_ID,colnames(data.input))
data <- data.input[vf.ko.list,MagInfo.input$Genome_ID%in%colnames(data.input)]
data.m <- data.frame(t(data))%>%rownames_to_column("Genome_ID")
data.m$BVBoth <- ifelse(data.m$K01186>0&data.m$K11031>0,"Both",ifelse(data.m$K01186>0,"Sialidase",ifelse(data.m$K11031>0,"Cytolysin","None")))
data.m <- merge(data.m,MagInfo.input,by="Genome_ID")%>%na.omit()
data.m$SGB_Species <- gsub(" ","_",data.m$SGB_Species)
xtabs(~SGB_Species+Genome_quality,data.m)

###  Genome_quality 选高质量和近乎完整的 
q.stat <- merge(data.m%>%group_by(SGB_Species,Genome_quality,BVBoth)%>%summarise(n=n()),
                data.m%>%group_by(SGB_Species)%>%summarise(N=n()))
q.stat$rate <- q.stat$n/q.stat$N
q.stat$Genome_quality2 <- ifelse(q.stat$Genome_quality=="near-complete","high-quality",as.vector(q.stat$Genome_quality))

### Selection ####
data.m <- data.m%>%subset(Genome_quality%in%c("high-quality","near-complete"))
data.stat <- merge(data.m%>%group_by(SGB_Species,BVBoth)%>%summarise(n=n()),
                   data.m%>%group_by(SGB_Species)%>%summarise(N=n()))
data.stat$rate <- data.stat$n/data.stat$N
data.stat$labels <- paste0(data.stat$n,"/",data.stat$N,"\n",round(data.stat$rate,2)*100,"%") 
data.stat <- data.stat %>% subset(N>10)
###
id.Both <- data.stat$SGB_Species[data.stat$BVBoth%in%c("Sialidase","Cytolysin","Both") & data.stat$rate>0.1]
data.bv.draw <- data.stat%>%subset(SGB_Species%in%id.Both)
data.bv.draw$rename <- taxa.info$rename[match(data.bv.draw$SGB_Species,taxa.info$names)]

## The two enzymes proportion in SGBs: at least 10% of any enzymes
temp <- dcast(data.bv.draw,rename~BVBoth,value.var = "rate")
temp[is.na(temp)] <- 0
p<-pheatmap::pheatmap(temp%>%column_to_rownames("rename"))
order.cluster <- p$tree_row$labels[p$tree_row$order]
##
data.bv.draw$rename <- factor(data.bv.draw$rename,levels = order.cluster)

## cutoff 
cf.bv <- 0.9
id.select <- data.bv.draw$rename[data.bv.draw$rate>cf.bv& data.bv.draw$BVBoth%in%c("Sialidase","Cytolysin","Both")]
data.select.withNone <- data.bv.draw%>%subset(rename%in%id.select)
prefix <- paste0("Prop_",cf.bv)
ggplot(data.select.withNone,aes(x=n,y=rename,fill=BVBoth,label=labels))+
  geom_col(position = "fill")+
  geom_text(size = 2,vjust = -1, position = position_fill(.5),angle=-90,color = "grey20")+
  scale_fill_manual(values = c( "#6A3D9A" ,"#FB9A99", "grey80","#CAB2D6"))+
  geom_vline(xintercept = c(1-cf.bv,cf.bv),linetype = "dotdash")+
  labs(title = "The proportion of 2 enzymes",
       subtitle = paste0("Each species: Num of MAGs > 10; any enzymes at least ",cf.bv),x="",y="",fill="")+theme_pubclean()
ggsave(paste0(outdir,prefix,".species.BV2enzy.pdf"),width = 8,height = 12)



