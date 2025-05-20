#haolilan
## R Package
library(vegan)
library(cluster)
library(factoextra)

## profile
profile_relab <- read.delim(args[1],header=T)
prefix <- args[2]

## dist
data <- t(profile_relab)
bray_curtis_dist <- vegdist(data, method = "bray")

## fit 
set.seed(123)
fit1 <- fviz_nbclust(data, pam, method = "silhouette", diss = bray_curtis_dist, k.max = 15)+labs(subtitle = "Silhouette method")
pdf(paste0(prefix,".Optimal_cluster.pdf"),width = 8,height = 6,onefile = T)
fit1
dev.off()
### best fit 
which.best <- fit1$data$clusters[which.max(fit1$data$y)]
object <-pam(bray_curtis_dist,k=which.best)
res.cluster <- data.frame(Samp = names(object$clustering),clustering = object$clustering)
write.table(res.cluster,paste0(prefix,".k=",which.best,".clustering.txt"),quote=F,row.name=F)
