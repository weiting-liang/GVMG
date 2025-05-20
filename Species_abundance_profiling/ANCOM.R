#haolilan
# Bacterial vaginosis specific signatures in bacterial SGBs
## https://github.com/FrederickHuangLin/ANCOM/blob/master/programs/ancom.R

diff_ancom <- function(profile, groups, outdir, prefix.file,cutoff_feature=0.0001, filter_preval=0.01,filter_N=10){
  source("./ancom.R") 
  library(compositions);library(readxl);library(dplyr);library(tibble)
  ##
  groups     <- groups %>% subset(SeqID %in% colnames(profile))
  groups$V2  <- factor(groups$V2,levels = c("control","case"))
  prof.input <- profile %>% dplyr::select(groups$SeqID)%>% .[,groups$SeqID]%>% .[rowSums(.)!=0,]
  ###
  prevalence <- data.frame(x= apply(prof.input, 1, function(x){sum(x>cutoff_feature)/length(x)}),
                           sum= apply(prof.input, 1, function(x){sum(x>cutoff_feature)})) %>% subset(x > filter_preval & sum > filter_N)
  prof.input <- prof.input[rownames(prevalence),] 
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
 
outdir     <- "./diff"
groups     <- read.table("./diff/match.BV.xw.csv",header = T,sep = ",");colnames(groups)[2]<-"V2"
table(groups$V2)
## SGB 
profile <-  read.delim("./data/Whole.relative_taxonomic_abundance.recal.txt")
prefix.file  <- "BV.xw.sgb"
intersect(groups$SeqID,colnames(profile)) %>%length()
diff_ancom(profile, groups, outdir, prefix.file,cutoff_feature=0.0001,filter_preval=0.01,filter_N=5)
