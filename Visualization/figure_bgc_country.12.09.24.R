#jiezhuye
#####--------------------------------
###        read data
####_---------------------------------

	#  read  core BGC coverage profile
		library(data.table)
		bgco = fread("BGC_coreCoverage_trim.csv",sep=",")
		bgco = as.data.frame(bgco)
		rownames(bgco) = bgco[,1]
		bgco=bgco[,-1]
		colnames(bgco) = gsub(".rmhost$","",gsub(".corecov","",colnames(bgco)))
		bgco[1:4,1:4]
		dim(bgco)

	# core BGC coverage 
		bgco[1:4,1:4]
		dim(bgco)

	# core BGC PRKM
	#  read  core BGC RPKM profile
		bgcpf = fread("BGC_coreRPKM_trim.csv",sep=",")
		bgcpf = as.data.frame(bgcpf)
		rownames(bgcpf) = bgcpf[,1]
		bgcpf=bgcpf[,-1]
		colnames(bgcpf) =  gsub(".rmhost$","",gsub(".coreRPKM","",colnames(bgcpf)))

		bgcpf[1:4,1:4]
		dim(bgcpf)
		rownames(bgcpf)[119]


		bgcpf[1:4,1:4]
	#### RPKM to RA
		bgcpf_j = sweep(bgcpf,2,apply(bgcpf,2,sum),"/")
		bgcpf_j[1:4,1:4]

	#### coverage align to RA
		bgco = bgco[match(rownames(bgcpf_j),rownames(bgco)),match(colnames(bgcpf_j),colnames(bgco))]
		bgcpf_j = t(bgcpf_j)

	# get bgc annotation
		# get bgc annotation
		library(openxlsx)
		bgc_a = read.xlsx("GCs.rename.MAGs.annotate.xlsx",sheet="annotate")
		colnames(bgc_a)[colnames(bgc_a)=="GCF_Type"]="Type"
		colnames(bgc_a)[colnames(bgc_a)=="rename_GCF"]="rename"
		

	
	# read sgb file
		#sgb - mag table
		tmp = read.xlsx("MAGs_genomes_stat241104_euk.xlsx",sheet="pro_SGB_index")
		
		cnt = read.xlsx("SGBs_table.xlsx",sheet="Sheet1")
		cnt$profid = tmp[match(cnt$index,tmp$index),"SGB"]
		rownames(cnt) = cnt$index


		## sgb information
		cnt <- read.csv("D:/Work/VIMAG2022/sgb_table.csv",row.names=1,head=T,as.is=T)
		k50 = cnt[,"number_of_public_isolation"]/cnt[,"mags_number"]


		#culture, uSGB, kSGB
		sp =unlist(lapply(strsplit(cnt$gbtk_anotation,";",fixed=T),function(x){x[length(x)]}))
		cnt$k50 = k50
		cnt$type = ifelse(sp=="s__","uSGB","kSGB")	
		cnt=cnt[order(cnt$number_of_public_isolation + cnt$number_of_culture_in_this_study,decreasing=T),]
		cnt=cnt[order(cnt$k50,decreasing=T),]
		head(cnt)	

	# mags information
	mags = read.xlsx("MAGs_genomes_stat241104_euk.xlsx",sheet="pro",colNames=T,startRow=2)
	head(mags)
	colnames(mags)[colnames(mags)=="Species-level_genomic_bin_(95%_ANI)"]="SGBid"

	# sample conf

		library(openxlsx)
		conf = read.xlsx("5-STable_20241202.xlsx",sheet="Table2",startRow =2)
			head(conf)
		conf$prfid = conf$BioSample.ID
		conf$prfid[!conf$BioSample.ID%in%colnames(bgco)] =NA
		conf$prfid[conf$Prof_ID%in%colnames(bgco)] = conf$Prof_ID[conf$Prof_ID%in%colnames(bgco)]
		conf$prfid[gsub("_A$","",conf$BioSample.ID)%in%colnames(bgco)]=gsub("_A$","",conf$BioSample.ID)[gsub("_A$","",conf$BioSample.ID)%in%colnames(bgco)]
		conf$prfid[gsub("_B$","",conf$BioSample.ID)%in%colnames(bgco)]=gsub("_B$","",conf$BioSample.ID)[gsub("_B$","",conf$BioSample.ID)%in%colnames(bgco)]

		mchpv = conf$Country[match(colnames(bgco),conf$prfid)]
		mchpv[!mchpv%in%c("China","USA") ]=NA
		table(mchpv)

	# BGC- mag
		library(data.table)
		bgc_mag = fread("MAGs_BGC_count.csv",sep=",")
		bgc_mag = as.data.frame(bgc_mag)
		rownames(bgc_mag) = bgc_mag[,1]
		bgc_mag= bgc_mag[,-1]
		bgc_mag[1:4,1:4]

		bgc_mag$SGBid =mags[ match(rownames(bgc_mag),mags$Genome_ID),"SGBid"]
		bgc_mag$sgb= bgc_mag$SGBid
		bgc_mag$taxon= cnt[pmatch(bgc_mag$SGBid,cnt$index,duplicates.ok=T),"rename"]
		bgc_mag$MAGs = rownames(bgc_mag)
		bgc_mag$rename = bgc_mag$taxon

	### color

		fill <- c("#5F9EA0", "#E1B378","#E64B35B2", "#F0E442","#00798c", "#b2d183","#999999", "#E69F00", "#56B4E9",
				pals::parula(10)[c(1,3,6)],"#B40404", "#0B6121", "#FFBF00","#F39B7FB2")
	


	# s016 bgc
		s016_bgc = colnames(bgc_mag)[which(apply(bgc_mag[bgc_mag$sgb=="S016.metaspades.concoct_dastools.bin.8.fa",],2,function(x){ sum(as.numeric(x))>0 }))]


	# fa bgc
		fa_bgc = colnames(bgc_mag)[which(apply(bgc_mag[bgc_mag$sgb=="MG14.megahit.concoct_dastools.bin.27.fa"|bgc_mag$sgb=="CES0010525.megahit.metabat2.bin.13.fa",]
				,2,function(x){ sum(as.numeric(x))>0 }))]
	#bf bgc, no bgc?

		bf_bgc = colnames(bgc_mag)[which(apply(bgc_mag[bgc_mag$sgb=="SXH005271.metaspades.maxbin2_dastools.bin.2.fa",]
				,2,function(x){ sum(as.numeric(x))>0 }))]
	# s__Bifidobacterium_breve
		bb_bgc = colnames(bgc_mag)[which(apply(bgc_mag[bgc_mag$sgb=="GCF_000213865.1_ASM21386v1_genomic.fa",]
				,2,function(x){ sum(as.numeric(x))>0 }))]

	#-------------------------------
	# chisq test
	#------------------------------------
		cut=0.5
		bgco1=  t(bgco)[,apply(t(bgco),2,function(x){ sum(x>cut)>50 })]
		bgco1=  bgco1[,apply(bgco1,2,function(x){ sum(x>cut)<nrow(bgco1)-50 })]
		dim(bgco1)
		pv=c()

		pv = apply(bgco1>=cut,2,function(x){
				chisq.test(x,factor(mchpv))$p.value
			})
		qv = p.adjust(pv,method="BH")
		rest1 = aggregate(bgco1>=cut,by=list(mchpv=factor(mchpv)),mean,na.rm=T)
			rownames(rest1)=rest1[,1]
			rest1 = rest1[,-1]
			rest1 = t(rest1)
			#match(names(pv),rownames(rest1))
			rest1 = data.frame(rest1,pv,qv)
			rest1$freqratio = rest1$China/rest1$USA
			freqx= rest1
			table(freqx$qv<0.05, freqx$freqratio>1)


	
			r1 = freqx[freqx$freqratio>1 & freqx$qv<1e-10 & freqx$China>0.05,]
			r1_freq1 = freqx
			# china enrich
			r1 = freqx[which(freqx$freqratio>1 & freqx$qv<1e-10 & freqx$China>=0.05),]
			r1 = data.frame(r1, bgc_a[match(rownames(r1),bgc_a$rename),c("Type","rename")])
			r1_freq1 = data.frame(r1_freq1, bgc_a[match(rownames(r1_freq1),bgc_a$rename),c("Type","rename")])

			sum(is.na(match(rownames(bgc_mag),mags$Genome_ID)))
			extract_host = function(r1){
				bx =c()
				for(i in 1:nrow(r1)){
					## sgb count
					#bb = tail(sort(table(mags[match(rownames(bgc_mag)[bgc_mag[,rownames(r1)[i]]>0],mags$Genome_ID),"SGBid"]
					bb = tail(sort(table(bgc_mag[bgc_mag[,rownames(r1)[i]]>0,"SGBid"])),1)
					#bb = sort(table(bgc_mag[bgc_mag[,rownames(r1)[i]]>0,2]))
					bx=rbind(bx,data.frame(sgb=cnt[gsub(".fa","",names(bb)),"rename"],num=as.numeric(bb)))
				}
				bx
			}
			r1 = data.frame(r1, extract_host(r1))
			r1_freq1 = data.frame(r1_freq1, extract_host(r1_freq1))
			dim(r1_freq1)
			write.csv(r1_freq1,"bgc_country_freq.csv")

			# virgo enrich
			r2 = freqx[which(freqx$freqratio<1 & freqx$qv<1e-10 & freqx$USA>=0.05),]
			r2 = r2[order(r2$qv)[1:15],]
			r2 = data.frame(r2, bgc_a[match(rownames(r2),bgc_a$rename),c("Type","rename")])
			r2 = data.frame(r2, extract_host(r2))


			r1 = rbind(r1,r2)

		r1_freq = r1
		
	#---------------------
	#wilcox.test
	#--------------------
		dim(bgcpf_j)

		pv=c()
		for(i in 1:ncol(bgcpf_j)){
			pv=c(pv,wilcox.test(bgcpf_j[,i]~factor(mchpv))$p.value)
		}

		pv = apply(bgcpf_j,2,function(x){
			wilcox.test(x~factor(mchpv))$p.value

		})
		qv = p.adjust(pv,method="BH")
		rest1 = aggregate(log10(bgcpf_j+1e-10),by=list(mchpv=factor(mchpv)),mean,na.rm=T)

			rownames(rest1)=rest1[,1]
			rest1 = rest1[,-1]
			rest1 = t(rest1)
			rest1 = data.frame(rest1,pv,qv)
			rest1$freqratio = rest1$China/rest1$USA

			r1 = rest1[which(rest1$freqratio<1 & rest1$qv<1e-10 & rest1$China>=log10(1e-9) ),]
			r1 = data.frame(r1, bgc_a[match(rownames(r1),bgc_a$rename),c("Type","rename")])
			extract_host = function(r1){
				bx =c()
				for(i in 1:nrow(r1)){
					## sgb count
					#bb = tail(sort(table(mags[match(rownames(bgc_mag)[bgc_mag[,rownames(r1)[i]]>0],mags$Genome_ID),"SGBid"]
					bb = tail(sort(table(bgc_mag[bgc_mag[,rownames(r1)[i]]>0,"SGBid"])),1)
					#bb = sort(table(bgc_mag[bgc_mag[,rownames(r1)[i]]>0,2]))
					bx=rbind(bx,data.frame(sgb=cnt[gsub(".fa","",names(bb)),"rename"],num=as.numeric(bb)))
				}
				bx
			}
			r1 = data.frame(r1, extract_host(r1))

			write.csv(data.frame(rest1,extract_host(rest1)),"bgc_country_abun.csv")

			# virgo enrich
			r2 = rest1[which(rest1$freqratio>1 & rest1$qv<1e-10 & rest1$USA>=log10(1e-9)),]
			r2 = r2[order(r2$qv)[1:15],]
			r2 = data.frame(r2, bgc_a[match(rownames(r2),bgc_a$rename),c("Type","rename")])
			r2 = data.frame(r2, extract_host(r2))


			r1 = rbind(r1,r2)
			table(rest1$qv<0.05,rest1$freqratio<1)



		#-----------------------------------------------
		#------ output 
		#--------------------------------------------------
			bxlist = unique(c(rownames(r1),rownames(r1_freq)))

			bxlist = bxlist[!is.na(match(bxlist,rownames(r1_freq1)))]


			rr = aggregate(apply(bgcpf_j,2,rank),by=list(mchpv=factor(mchpv)),mean,na.rm=T)
			rownames(rr)=rr[,1]
			rr= rr[,-1]
			rr=data.frame( t(rr))
			rr$enrich = ifelse(rr[,1]>rr[,2],1,0)
			xx1 = rest1[match(bxlist,rownames(rest1)),]
			xx1$enrich =rr[match(bxlist,rownames(rest1)),"enrich"]
			xx = r1_freq1[match(bxlist,rownames(r1_freq1)),]
			table( sign( xx1$enrich==0) == sign(xx$freqratio<1))
			bxlist = bxlist[sign( xx1$enrich==0) == sign(xx$freqratio<1)]

			xx1 = xx1[match(bxlist,rownames(xx1)),]
			xx1=xx1[order(xx1$freqratio),]
			bxlist= rownames(xx1)


			tmp = melt(as.matrix(r1_freq1[match(bxlist,rownames(r1_freq1)),c("China","USA")]))
			tmp$Var1 = factor(tmp$Var1,levels=bxlist)
			tmp$ord = as.numeric(factor(tmp$Var1,levels=bxlist))
			tmp$rename = paste0( gsub("vag","",tmp$Var1), "-", bgc_a[ pmatch(tmp$Var1,bgc_a$rename,duplicates.ok=T),"Type"])

			library(ggplot2)
 			g1 = ggplot(data=tmp) + geom_bar(aes(reorder(rename,ord),value, fill=factor(Var2)),
                           stat="identity", position=position_dodge())+
				 labs(x="", y="Prevalence") +
  				ggtitle("")+
				 scale_fill_manual(values=fill)+ theme_bw() +
				  theme(panel.grid.major.y = element_blank(),
	 				axis.text.x = element_text(angle = 90, hjust = 1)) +
				  theme(axis.line = element_line(size=1, colour = "black"),
       					 panel.grid.major = element_line(colour = "#d3d3d3"), 
							panel.grid.minor = element_blank(),
							legend.position = "none" ,
       					 panel.border = element_blank(), panel.background = element_blank())+
					coord_flip()
		
	
			tmp = melt(as.matrix(bgcpf_j[,match(bxlist,colnames(bgcpf_j))]))
			tmp = data.frame(tmp,locate=mchpv[pmatch(tmp$Var1,rownames(bgcpf_j),duplicates.ok=T)])
			tmp$Var2 = factor(tmp$Var2,levels=bxlist)
			tmp$ord = as.numeric(factor(tmp$Var2,levels=bxlist))
			tmp$rename = paste0( gsub("vag","",tmp$Var2), "-", bgc_a[ pmatch(tmp$Var2,bgc_a$rename,duplicates.ok=T),"Type"])
			tmp$label = r1_freq1[pmatch(as.character(tmp$Var2),rownames(r1_freq1),duplicates.ok=T),"sgb"]

library(dplyr)
library(grid)

label_summary <- tmp %>%
  group_by(rename) %>%
  summarize(label = first(label), ord = first(ord))

		g2 = ggplot(na.omit(tmp), aes(reorder(rename,ord),log10(value+1e-10), fill=factor(locate))) +
  				geom_boxplot(outlier.shape = NA)+
				 labs(x="", y="Log10RA") +
  				ggtitle("")+
				 scale_fill_manual(values=fill)+ theme_bw() +
				  theme(panel.grid.major.y = element_blank(),
	 				axis.text.y = element_blank() ) +
				  theme(axis.line = element_line(size=1, colour = "black"),
       					 panel.grid.major = element_line(colour = "#d3d3d3"), 
							panel.grid.minor = element_blank(),
						strip.text = element_text(size = 8),
       					 panel.border = element_blank(), panel.background = element_blank())+
					coord_flip()+
			 geom_text(data = label_summary, aes(x = reorder(rename, ord), y = max(log10(0 + 1e-10)) + 0.5, label = label), 
            hjust = 0, vjust = 0, size = 3, inherit.aes = FALSE)

		pdf("bgc_country.pdf",  
						 width = 7, height = 5,onefile = F)
		egg::ggarrange(
  	g1,g2,
  	nrow = 1,ncol = 2,
 	 widths = c(1,1)
  ) 
		dev.off()








