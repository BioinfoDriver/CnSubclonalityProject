

# load
gli.cn.alt.frac <- readRDS(file='/data/gli_glod_cn_alt_frac.rds')

gli.cn.alt.frac$Integrated_Diagnoses <- factor(gli.cn.alt.frac$Integrated_Diagnoses, 
 levels=c('Oligodendroglioma,IDHmut-codel', 'Astrocytoma,IDHmut', 'Glioblastoma,IDHwt'))

gli.cn.alt.frac <- dplyr::arrange(gli.cn.alt.frac, Integrated_Diagnoses, desc(subclo_cn_alt_frac))

# Genome alteration fraction
genome.alt.frac <- data.frame(sample=rep(x=gli.cn.alt.frac$bcr_patient_barcode, times=4),
 alt.frac=c(gli.cn.alt.frac$subclo_amp_genome_frac, gli.cn.alt.frac$clo_amp_genome_frac, 
 gli.cn.alt.frac$subclo_del_genome_frac, gli.cn.alt.frac$clo_del_genome_frac), 
 type=rep(x=c('Subclonal gain', 'Clonal gain', 'Subclonal loss', 'Clonal loss'), each=nrow(gli.cn.alt.frac)))

genome.alt.frac$sample <- factor(genome.alt.frac$sample, levels=gli.cn.alt.frac$bcr_patient_barcode)
genome.alt.frac$type <- factor(genome.alt.frac$type, levels=c('Subclonal gain', 'Subclonal loss', 'Clonal gain', 'Clonal loss'))


# Subclonal copy numbuer alteration proportion
cn.alt.prop <- data.frame(sample=rep(x=gli.cn.alt.frac$bcr_patient_barcode, times=4),
 alt.prop=c(gli.cn.alt.frac$subclo_amp_alt_frac, gli.cn.alt.frac$clo_amp_alt_frac, 
 gli.cn.alt.frac$subclo_del_alt_frac, gli.cn.alt.frac$clo_del_alt_frac), 
 type=rep(x=c('Subclonal gain', 'Clonal gain', 'Subclonal loss', 'Clonal loss'), each=nrow(gli.cn.alt.frac)))

cn.alt.prop$sample <- factor(cn.alt.prop$sample, levels=gli.cn.alt.frac$bcr_patient_barcode)
cn.alt.prop$type <- factor(cn.alt.prop$type, levels=c('Subclonal gain', 'Subclonal loss', 'Clonal gain', 'Clonal loss'))
 
# plot
library('ggpubr')

gen.alt.plot.part1 <- ggbarplot(genome.alt.frac, x="sample", y="alt.frac", fill="type", color=NA, width=1, 
 palette=c("#E59398", "#8CB6D2", "#D42527", "#2171A9"), 
 sort.by.groups=TRUE, xlab=FALSE, ylab='Genome alteration fraction', legend.title='Alteration type') + 
 theme(axis.text.x=element_blank()) + theme(axis.ticks.x=element_blank()) + coord_cartesian(ylim = c(0, 0.30))

gen.alt.plot.part2 <- ggbarplot(genome.alt.frac, x="sample", y="alt.frac", fill="type", color=NA, width=1, 
 palette=c("#E59398", "#8CB6D2", "#D42527", "#2171A9"), 
 sort.by.groups=TRUE, xlab=FALSE, ylab='Genome alteration fraction', legend.title='Alteration type') + 
 theme(axis.text.x=element_blank()) + theme(axis.ticks.x=element_blank()) + coord_cartesian(ylim = c(0.35, 0.90))

cn.alt.plot <- ggbarplot(cn.alt.prop, x="sample", y="alt.prop", fill="type", color=NA, width=1,
 palette=c("#E59398", "#8CB6D2", "#D42527", "#2171A9"), sort.by.groups=TRUE, xlab=FALSE, ylab='SCNAs percentage', 
 legend.title = 'Alteration type') + theme(axis.text.x=element_blank()) + theme(axis.ticks.x=element_blank())

alt.plot <- ggpubr::ggarrange(gen.alt.plot.part2, gen.alt.plot.part1, cn.alt.plot, heights=c(0.3, 0.4, 0.3), 
  legend="top", common.legend=TRUE, ncol=1, nrow=3, align="v") 

ggsave(alt.plot, file='/result/Section1/cn_alt_frac.pdf')



# Hatmapannotation
library('circlize')
library('RColorBrewer')
library('ComplexHeatmap')
anno.dat <- gli.cn.alt.frac[, c('Integrated_Diagnoses', 'age', 'histological_grade', 'CDKN2AB', 
 'EGFR', 'CHR_7_GAIN_CHR_10_LOSS', 'TERT_PROMOTER_STATUS', 'ATRX_STATUS', 'MGMT_PROMOTER_STATUS')]

anno.dat$CDKN2AB <- ifelse(anno.dat$CDKN2AB == 1, 'Deletion', 'Diploid')
anno.dat$EGFR <- ifelse(anno.dat$EGFR == 1, 'Amplification', 'Diploid')
anno.dat$CHR_7_GAIN_CHR_10_LOSS <- ifelse(anno.dat$CHR_7_GAIN_CHR_10_LOSS == 'Gain chr 7 & loss chr 10', 
 'Chr +7/-10', 'No combined CNA')


pdf('/result/Section1/subclo_genome_heatmap.pdf')
 heat.anno <- HeatmapAnnotation(df=anno.dat, 
 col=list(Integrated_Diagnoses=setNames(c("#00AFBB", "#E7B800", "#FC4E07"), 
   c('Oligodendroglioma,IDHmut-codel', 'Astrocytoma,IDHmut', 'Glioblastoma,IDHwt')), 
 age=colorRamp2(c(0, 80), c("white", "red")), 
 histological_grade=setNames(brewer.pal(8, 'BuPu')[c(3, 4, 5)], c('G2', 'G3', 'G4')),  
 CDKN2AB=setNames(brewer.pal(8, 'Set1')[2:3], c("Deletion", "Diploid")), 
 EGFR=setNames(brewer.pal(8, 'Dark2')[1:2], c("Amplification", "Diploid")), 
 CHR_7_GAIN_CHR_10_LOSS=setNames(brewer.pal(12, 'Paired')[c(3, 12)], c("Chr +7/-10", "No combined CNA")), 
 TERT_PROMOTER_STATUS=setNames(brewer.pal(8, 'Set3')[4:5], c("Mutant", "WT")),  
 ATRX_STATUS=setNames(brewer.pal(8, 'Set3')[6:7], c("Mutant", "WT")),
 MGMT_PROMOTER_STATUS=setNames(brewer.pal(10, 'Set3')[c(1, 10)], c("Methylated", "Unmethylated"))))		   	   
 
 zero_row_mat <- matrix(nrow = 0, ncol = nrow(anno.dat))
 Heatmap(zero_row_mat, top_annotation = heat.anno)
dev.off()



