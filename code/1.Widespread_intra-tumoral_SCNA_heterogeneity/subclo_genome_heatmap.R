

# load
gli.cn.alt.frac <- readRDS(file='/data/gli_glod_cn_alt_frac.rds')
gli.cn.alt.frac <- dplyr::arrange(gli.cn.alt.frac, IDH_CODEL_SUBTYPE, desc(subclo_cn_alt_frac))

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
anno.dat <- gli.cn.alt.frac[, c('IDH_CODEL_SUBTYPE', 'cancer_type', 'molecular_histological_type', 'TRANSCRIPTOME_SUBTYPE', 
 'histological_grade', 'gender', 'age', 'MGMT_PROMOTER_STATUS', 'TERT_PROMOTER_STATUS', 'ATRX_STATUS')]

pdf('/result/Section1/subclo_genome_heatmap.pdf')
 heat.anno <- HeatmapAnnotation(df=anno.dat, 
 col=list(cancer_type=setNames(brewer.pal(8, 'Greys')[c(4, 6)], c('LGG', 'GBM')),
 IDH_CODEL_SUBTYPE=setNames(c("#00AFBB", "#E7B800", "#FC4E07"), c("IDHwt", "IDHmut-non-codel", "IDHmut-codel")), 
 molecular_histological_type=setNames(brewer.pal(8, 'Set1')[2:4], c("Glioblastoma", "Astrocytoma", "Oligodendroglioma")), 
 TRANSCRIPTOME_SUBTYPE=setNames(brewer.pal(8, 'Dark2')[1:4], c("ME", "CL", "NE", "PN")), 
 histological_grade=setNames(brewer.pal(8, 'BuPu')[c(3, 4, 5)], c('G2', 'G3', 'G4')), 
 gender=setNames(brewer.pal(12, 'Paired')[c(3, 12)], c('MALE', 'FEMALE')), 
 age=colorRamp2(c(0, 80), c("white", "red")), 
 MGMT_PROMOTER_STATUS=setNames(brewer.pal(10, 'Set3')[c(1, 10)], c("Methylated", "Unmethylated")), 
 TERT_PROMOTER_STATUS=setNames(brewer.pal(8, 'Set3')[4:5], c("Mutant", "WT")), 
 ATRX_STATUS=setNames(brewer.pal(8, 'Set3')[6:7], c("Mutant", "WT"))))		   	   
 
 zero_row_mat <- matrix(nrow = 0, ncol = nrow(anno.dat))
 Heatmap(zero_row_mat, top_annotation = heat.anno)
dev.off()





