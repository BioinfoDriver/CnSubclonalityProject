
# RNA_SEQ expression FPKM
# setwd('/pub5/xiaoyun/Jobs/J22/EvoClass2.0/Revise/de_expression/FPKM')
# samples <- read.csv(file = 'gdc_sample_sheet.2020-07-22.tsv', sep = '\t', header = TRUE, stringsAsFactors = FALSE)
# samples$filename <- paste(samples$File.ID, samples$File.Name, sep = "/")
# samples <- subset(samples, Sample.Type == 'Primary Tumor')
# samples <- samples[!duplicated(samples$Case.ID), ]

# library(edgeR)
# fpkm <- readDGE(samples$filename, columns=c(1, 2))
# colnames(fpkm) <- samples$Case.ID
# saveRDS(fpkm, file='/pub5/xiaoyun/Jobs/J22/CopyNumberClonalityProject/Results/Section3/Resources/pan_gli_exp_fpkm.rds')

gene.exp <- readRDS(file='/pub5/xiaoyun/Jobs/J22/CopyNumberClonalityProject/Results/Section3/Resources/pan_gli_exp_tmm.rds')

# Tumor proliferating gene
setwd('/pub5/xiaoyun/Jobs/J22/EvoClass2.0/Section4/Results/Proliferating')
pcna.sig <- read.csv(file = 'PCNA_Signature.txt', header = TRUE,sep = '\t',stringsAsFactors =FALSE)

library('Homo.sapiens')
pcna.sig <- select(Homo.sapiens, keys = pcna.sig$Entrez.ID, 
	columns = c("GENENAME", "SYMBOL", "ENSEMBL"), keytype='ENTREZID')

######## clonal rtk alt Vs subclonal rtk alt
setwd('/pub5/xiaoyun/Jobs/J22/CopyNumberClonalityProject/Results/Section3/Resources')
rtks <- read.csv(file='RTKs.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)

gold.set <- readRDS('/pub5/xiaoyun/Jobs/J22/CopyNumberClonalityProject/Resource/CuratedData/gold_set.rds')
gene.het <- readRDS(file='/pub5/xiaoyun/Jobs/J22/CopyNumberClonalityProject/Results/Section3/Resources/gene.het.rds')
gene.het <- gene.het[, paste(gold.set, '01', sep = '-')]

RtkCliMolAnalysis <- function(het.mat, rtk, group.index){
 
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


 return(rtk.clo.mat)
}

rtk.alt.sam <- RtkCliMolAnalysis(gene.het, rtks$Approved.symbol, 'exclude')

rtk.gene.exp <- gene.exp[, intersect(colnames(gene.exp), rownames(rtk.alt.sam))]
rtk.alt.sam <- rtk.alt.sam[colnames(rtk.gene.exp), , FALSE]
rownames(rtk.gene.exp) <- substr(rownames(rtk.gene.exp), 1, 15)

# proliferating score——median score
pcna.gene.exp <- rtk.gene.exp[intersect(pcna.sig$ENSEMBL, rownames(rtk.gene.exp)), ]
pcna.median.score <- apply(pcna.gene.exp, 2, median)	
rtk.alt.sam$pcna.median.score <- pcna.median.score[rownames(rtk.alt.sam)]


#  proliferating score——GSVA score
# procoding genes
setwd('/pub5/xiaoyun/Jobs/J22/CopyNumberClonalityProject/Results/Section3/Resources')
gene.info <- read.csv(file='gene_with_protein_product.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
gene.info <- gene.info[, c('hgnc_id', 'ensembl_gene_id')]

all.gene.exp <- rtk.gene.exp[intersect(rownames(rtk.gene.exp), gene.info$ensembl_gene_id), ]
library(GSVA)
pcna.gsva.score <- gsva(expr=all.gene.exp, 
 gset.idx.list=list(pcna=intersect(pcna.sig$ENSEMBL, rownames(all.gene.exp))), method="gsva")
rtk.alt.sam$pcna.gsva.score <- pcna.gsva.score[, rownames(rtk.alt.sam)]

pcna.ssgsea.score <- gsva(expr=all.gene.exp, 
 gset.idx.list=list(pcna=intersect(pcna.sig$ENSEMBL, rownames(all.gene.exp))), method="ssgsea")
rtk.alt.sam$pcna.ssgsea.score <- pcna.ssgsea.score[, rownames(rtk.alt.sam)]


######### Tumor proliferating score compare
tcga.cli.data <- readRDS('/pub5/xiaoyun/Jobs/J22/CopyNumberClonalityProject/Resource/CuratedData/tcga_glioma_cli_mol.rds')
idhwt.gbm <- subset(tcga.cli.data, histological_grade %in% c('G4') & 
 IDH_CODEL_SUBTYPE %in% c('IDHwt'))$bcr_patient_barcode
idhwt.gbm <- intersect(paste0(idhwt.gbm, '-01'), rownames(rtk.alt.sam))

wilcox.test(pcna.median.score~group, rtk.alt.sam[idhwt.gbm, ], alternative = 'two.sided')$p.value # 0.04407467
wilcox.test(pcna.ssgsea.score~group, rtk.alt.sam[idhwt.gbm, ], alternative = 'two.sided')$p.value # 0.03344954
wilcox.test(pcna.gsva.score~group, rtk.alt.sam[idhwt.gbm, ], alternative = 'two.sided')$p.value # 0.08611389


# plot 
library('ggplot2')
pcna.score.plot <- ggplot(data = rtk.alt.sam[idhwt.gbm, ], aes(group, pcna.median.score)) +
 geom_boxplot(aes(fill = group), outlier.shape = NA) + geom_jitter(shape=16, position=position_jitter(0.2)) + 
 xlab('RTKs subtype') + ylab('Proliferating index') + guides(fill = guide_legend(title= 'RTKs subtype')) + 
 theme(axis.text.x = element_text(angle=30, vjust=0.5)) + stat_compare_means(label.x.npc = 'left', label.y.npc = 'top')

ggsave(filename='rtk_median_score_com.pdf', plot=pcna.score.plot, 
 path='/pub5/xiaoyun/Jobs/J22/CopyNumberClonalityProject/Results/Section3/Results/RTKProliferation')



