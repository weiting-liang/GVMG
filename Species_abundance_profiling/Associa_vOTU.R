#haolilan
# Association between vOTUs in population and phenotypes (Generalized Linear Model)

## SU-CCS  ############
# function 
glm_2factors <- function(profile,phen,prefix.file){
  
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
prefix.file <- "Suzhou"
phen.SEL<- phen.merge.profile %>%subset(Source%in%c("SU-CCS2018/2019 (Peacock)","SU-CCS2021 (Peacock)" ))%>%
  dplyr::select(SeqID,Menopause:VVC)
## keep numeric variates
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
res.lm<-glm_2factors(profile,phen,prefix.file)
 
res.lm.INFO<- merge(res.lm,vOTUInfo%>%select(taxa_prof:`HPV type`)%>%rename(names=vOTU),all.x = T)
 
write.table(res.lm.INFO,paste0(outdir,"/lm.out.p.Suzhou.vOTUInfo.txt"),quote = F,row.names = T,sep = "\t")
