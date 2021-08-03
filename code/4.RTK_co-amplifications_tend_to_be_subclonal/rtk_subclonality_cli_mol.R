
#####################Survival
setwd('/pub5/xiaoyun/Jobs/J22/CopyNumberClonalityProject/Results/Section3/Resources')
rtks <- read.csv(file='RTKs.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)

gold.set <- readRDS('/pub5/xiaoyun/Jobs/J22/CopyNumberClonalityProject/Resource/CuratedData/gold_set.rds')
gene.het <- readRDS(file='/pub5/xiaoyun/Jobs/J22/CopyNumberClonalityProject/Results/Section3/Resources/gene.het.5.0.rds')
gene.het <- gene.het[, paste(gold.set, '01', sep = '-')]

tcga.cli.data <- readRDS('/pub5/xiaoyun/Jobs/J22/CopyNumberClonalityProject/Resource/CuratedData/tcga_glioma_cli_mol.rds')
tcga.cli.data <- tcga.cli.data[gold.set, ]
rownames(tcga.cli.data) <- paste0(rownames(tcga.cli.data), '-01', sep='')


# 
RtkCliMolAnalysis <- function(het.mat, rtk, group.index, cli.data){
 
 rtk.het <- het.mat[rtk, ]

 rtk.clo.sam <- colnames(rtk.het[, colSums(rtk.het == 1) >= 1])
 rtk.sub.sam <- colnames(rtk.het[, colSums(rtk.het == 2) >= 1])
 
 if(group.index=='subclonaldominant'){
  rtk.clo.sam <- setdiff(rtk.clo.sam, rtk.sub.sam)

 }else if(group.index=='clonaldominant'){
  rtk.sub.sam <- setdiff(rtk.sub.sam, rtk.clo.sam)
 
 }else{
  com.sam <- intersect(rtk.clo.sam, rtk.sub.sam)
  rtk.clo.sam <- setdiff(rtk.clo.sam, com.sam)
  rtk.sub.sam <- setdiff(rtk.sub.sam, com.sam)
  
 }
 
 rtk.clo.mat <- data.frame(group = rep(c('rtk clonal amp', 'rtk subclonal amp'), 
  times=c(length(rtk.clo.sam), length(rtk.sub.sam))))
 rownames(rtk.clo.mat) <- c(rtk.clo.sam, rtk.sub.sam)

 cli.data <- merge(cli.data, rtk.clo.mat, by = 'row.names')
 
 return(cli.data)
}

rtk.cli.mol <- RtkCliMolAnalysis(gene.het, rtks$Approved.symbol, 'exclude', tcga.cli.data)

