#haolilan
args <- commandArgs(T)
##  profile  prefix
###
library(vegan)
library(cluster)
library(factoextra)
library(RColorBrewer)
library(dplyr)
###
profile_relab <- read.delim(args[1],header=T)
prefix <- args[2]
##
data <- t(profile_relab)
dist <- vegdist(data, method = "bray")
write.table(as.matrix(dist),paste0(prefix,".dist.txt"),quote = F,sep = "\t")

dist.lower <- dist
dist.lower[upper.tri(dist.lower)] <- NA
dist.melt.raw  <- dist.lower %>% rownames_to_column("source") %>% melt(id = "source", variable.name ="variable") %>% subset(!(value %in% c(0, NA)))
write.table(dist.melt.raw,paste0(prefix,".dist.melt.txt"),quote = F,sep = "\t")


## set.seed(123)
set.seed(123)
fit1 <- fviz_nbclust(data, pam, method = "silhouette", diss = bray_curtis_dist, k.max = 15)+labs(subtitle = "Silhouette method")

pdf(paste0(prefix,".Optimal_cluster.pdf"),width = 8,height = 6,onefile = T)
fit1
dev.off()

### fit1 
which.best <- fit1$data$clusters[which.max(fit1$data$y)]
object <-pam(bray_curtis_dist,k=which.best)
res.cluster <- data.frame(Samp = names(object$clustering),clustering = object$clustering)
write.table(res.cluster,paste0(prefix,".k=",which.best,".clustering.txt"),quote=F,row.name=F)

## taxa.top.mean
groups <- res.cluster;names(groups) <- c("V1","V2")
feature<-c()
top.n <- 10
for (k in seq(which.best)) {
  d = profile_relab[,groups$V1[groups$V2==k]]
  feature.temp <- data.frame(taxa=rownames(d),mean=apply(d, 1, mean),cluster=k)%>%arrange(-mean)%>%.[1:top.n,]#rownames(d)[order(apply(d, 1, mean),decreasing=TRUE)[1:10]]
  feature <- rbind(feature,feature.temp)
  write.table(feature,paste0(prefix,".k=",which.best,".Feature_meanTOP.txt"),quote=F,row.name=F)
}

ggplot(feature,aes(cluster,mean,fill=taxa))+geom_col(position = "stack")+scale_fill_manual(values=c(brewer.pal(12,"Set3"),brewer.pal(12,"Paired"),brewer.pal(8,"Set2"),brewer.pal(8,"Set1"),rep("grey",40)))
ggsave(paste0(prefix,".k=",which.best,".Feature_meanTOP.pdf"),width = 12,height = 8)

