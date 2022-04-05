
# Methylation subtype
rtks <- read.csv(file='/data/RTKs.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)

gold.set <- readRDS('/data/gold_set.rds')
gene.het <- readRDS(file='/data/gene.het.rds')
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
rtk.alt.sam <- tibble::rownames_to_column(rtk.alt.sam, var = "bcr_patient_barcode")

tcga.cli.data <- readRDS('/data/tcga_glioma_cli_mol.rds')
gbm.cli.data <- subset(tcga.cli.data, Integrated_Diagnoses %in% c('Glioblastoma,IDHwt'))
gbm.cli.data$bcr_patient_barcode <- paste0(gbm.cli.data$bcr_patient_barcode, '-01')

gbm.cli.data <- merge(gbm.cli.data, rtk.alt.sam, by='bcr_patient_barcode')
# chisq.test(table(gbm.cli.data[, c('SUPERVISED_DNA_METHYLATION_CLUSTER', 'group')]))$p.value # 0.5428189

meth.stat <- data.frame(Freq=c(c(30, 3, 56)/89, c(7, 2, 18)/27), group=rep(c('RTKs clonal', 'RTKs subclonal'), each=3), 
 subtype=rep(c('Classic-like', 'LGm6-GBM', 'Mesenchymal-like'), times=2))


# transcription subtype
idhwt.gbm.subtype <- read.csv(file='/data/OriginalData/IDHwt_GBM_Exp_Subgroup.txt', 
 header=TRUE, sep='\t')
idhwt.gbm.subtype$sampleId <- gsub('\\.', '-', idhwt.gbm.subtype$sampleId)


idhwt.gbm.subtype <- merge(idhwt.gbm.subtype, rtk.alt.sam, by.x='sampleId', by.y='bcr_patient_barcode')
idhwt.gbm.subtype <- subset(idhwt.gbm.subtype, sampleId %in% 
 paste0(subset(tcga.cli.data, Integrated_Diagnoses %in% "Glioblastoma,IDHwt")$bcr_patient_barcode, '-01'))

# chisq.test(table(idhwt.gbm.subtype[, c('Group', 'group')]))$p.value # 0.01229864
# wilcox.test(SIMS~group, idhwt.gbm.subtype, alternative = 'two.sided')$p.value # 0.1958251

tran.stat <- data.frame(Freq=c(c(38, 22, 16)/76, c(6, 8, 13)/27), group=rep(c('RTKs clonal', 'RTKs subclonal'), each=3), 
 subtype=rep(c('Classical', 'Mesenchymal', 'Proneural'), times=2))



# Plot 
p1 <- ggbarplot(meth.stat, "group", "Freq", fill = "subtype", color = "subtype", 
 ylab = FALSE, xlab = FALSE, label = TRUE, lab.col = "white", lab.pos = "in", lab.size = 3, lab.vjust=2, 
 legend = 'none', x.text.angle=30) + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(size=7))

p2 <- ggbarplot(tran.stat, "group", "Freq", fill = "subtype", color = "subtype", 
 ylab = FALSE, xlab = FALSE, label = TRUE, lab.col = "white", lab.pos = "in", lab.size = 3, lab.vjust=2, 
 legend = 'none', x.text.angle=30) + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(size=7))

p3 <- ggplot(data = idhwt.gbm.subtype, aes(group, SIMS)) +
 geom_boxplot(aes(fill = group), outlier.shape = NA) + geom_jitter(shape=16, position=position_jitter(0.2)) + 
 xlab('RTKs subtype') + ylab('Transcriptional heterogeneity') + guides(fill = guide_legend(title= 'RTKs subtype')) + 
 theme(axis.text.x = element_text(angle=30, vjust=0.5)) + stat_compare_means(label.x.npc = 'left', label.y.npc = 'top')


ggsave(ggarrange(p1, p2, p3, ncol=3, nrow=2, common.legend = TRUE), 
 file='/result/Section4/rtk_subtype_com.pdf')


