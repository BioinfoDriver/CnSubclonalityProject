
# load data
rtks <- read.csv(file='/data/RTKs.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)

gold.set <- readRDS('/data/gold_set.rds')
gene.het <- readRDS(file='/data/gene.het.rds')
gene.het <- gene.het[, paste(gold.set, '01', sep = '-')]

# The patients of RTK clonal or subclonal amplification
rtk.het <- gene.het[rtks$Approved.symbol, ]
rtk.clo.sam <- colnames(rtk.het[, colSums(rtk.het == 1) >= 1])
rtk.sub.sam <- colnames(rtk.het[, colSums(rtk.het == 2) >= 1])
 
com.sam <- intersect(rtk.clo.sam, rtk.sub.sam)
rtk.clo.sam <- setdiff(rtk.clo.sam, com.sam)
rtk.sub.sam <- setdiff(rtk.sub.sam, com.sam)
 
# ABSOLUTE-based tumour purity
abs.puri.ploi <- readRDS('/data/tcga_gli_puri_ploi.rds')
rtk.clonality.purity <- abs.puri.ploi[c(rtk.clo.sam, rtk.sub.sam), ]
rtk.clonality.purity$rtk_status <- rep(c('Clonal', 'Subclonal'), times=c(length(rtk.clo.sam), length(rtk.sub.sam)))

# purity comparison
wilcox.test(purity ~ rtk_status, rtk.clonality.purity, alternative = 'two.sided')$p.value # 0.833399


# load SCNA burden
gli.cn.alt.frac <- readRDS(file='/data/gli_glod_cn_alt_frac.rds')
rtk.clonality.scna <- gli.cn.alt.frac[c(rtk.clo.sam, rtk.sub.sam), ]
rtk.clonality.scna$rtk_status <- rep(c('Clonal', 'Subclonal'), times=c(length(rtk.clo.sam), length(rtk.sub.sam)))


# SCNA comparison
wilcox.test(non_neutral_genome_frac ~ rtk_status, rtk.clonality.scna, alternative = 'two.sided')$p.value # 0.03527821
wilcox.test(clo_genome_frac ~ rtk_status, rtk.clonality.scna, alternative = 'two.sided')$p.value # 6.235393e-05
wilcox.test(subclo_genome_frac ~ rtk_status, rtk.clonality.scna, alternative = 'two.sided')$p.value # 1.833714e-10


# plot 
library('ggplot2')
p1 <- ggplot(data = rtk.clonality.purity, aes(rtk_status, purity)) +
 geom_boxplot(aes(fill = rtk_status), outlier.shape = NA) + geom_jitter(shape=16, position=position_jitter(0.2)) + 
 xlab('RTKs subtype') + ylab('Tumor purity') + guides(fill = guide_legend(title= 'RTKs subtype')) + 
 theme(axis.text.x = element_text(angle=30, vjust=0.5)) + stat_compare_means(label.x.npc = 'left', label.y.npc = 'top')

p2 <- ggplot(data = rtk.clonality.scna, aes(rtk_status, non_neutral_genome_frac)) +
 geom_boxplot(aes(fill = rtk_status), outlier.shape = NA) + geom_jitter(shape=16, position=position_jitter(0.2)) + 
 xlab('RTKs subtype') + ylab('Burden of SCNAs') + guides(fill = guide_legend(title= 'RTKs subtype')) + 
 theme(axis.text.x = element_text(angle=30, vjust=0.5)) + stat_compare_means(label.x.npc = 'left', label.y.npc = 'top')

p3 <- ggplot(data = rtk.clonality.scna, aes(rtk_status, clo_genome_frac)) +
 geom_boxplot(aes(fill = rtk_status), outlier.shape = NA) + geom_jitter(shape=16, position=position_jitter(0.2)) + 
 xlab('RTKs subtype') + ylab('Burden of clonal SCNAs') + guides(fill = guide_legend(title= 'RTKs subtype')) + 
 theme(axis.text.x = element_text(angle=30, vjust=0.5)) + stat_compare_means(label.x.npc = 'left', label.y.npc = 'top')

p4 <- ggplot(data = rtk.clonality.scna, aes(rtk_status, subclo_genome_frac)) +
 geom_boxplot(aes(fill = rtk_status), outlier.shape = NA) + geom_jitter(shape=16, position=position_jitter(0.2)) + 
 xlab('RTKs subtype') + ylab('Burden of subclonal SCNAs') + guides(fill = guide_legend(title= 'RTKs subtype')) + 
 theme(axis.text.x = element_text(angle=30, vjust=0.5)) + stat_compare_means(label.x.npc = 'left', label.y.npc = 'top')


ggsave(ggarrange(p1, p2, p3, p4, ncol=3, nrow=2, common.legend = TRUE), 
 file='/result/Section4/rtk_clonality_purity_scna.pdf')

###############################Glioblastoma

# purity comparison
gbm.rtk.clonality.purity <- rtk.clonality.purity[rownames(subset(rtk.clonality.scna, Integrated_Diagnoses == 'Glioblastoma,IDHwt')), ]
wilcox.test(purity ~ rtk_status, gbm.rtk.clonality.purity, alternative = 'two.sided')$p.value # 0.4206504


# SCNA comparison
gbm.rtk.clonality.scna <- subset(rtk.clonality.scna, Integrated_Diagnoses == 'Glioblastoma,IDHwt')
wilcox.test(non_neutral_genome_frac ~ rtk_status, gbm.rtk.clonality.scna, alternative = 'two.sided')$p.value # 2.104845e-06
wilcox.test(clo_genome_frac ~ rtk_status, gbm.rtk.clonality.scna, alternative = 'two.sided')$p.value # 0.01383941
wilcox.test(subclo_genome_frac ~ rtk_status, gbm.rtk.clonality.scna, alternative = 'two.sided')$p.value # 2.330514e-09


# plot 
p1 <- ggplot(data = gbm.rtk.clonality.purity, aes(rtk_status, purity)) +
 geom_boxplot(aes(fill = rtk_status), outlier.shape = NA) + geom_jitter(shape=16, position=position_jitter(0.2)) + 
 xlab('RTKs subtype') + ylab('Tumor purity') + guides(fill = guide_legend(title= 'RTKs subtype')) + 
 theme(axis.text.x = element_text(angle=30, vjust=0.5)) + stat_compare_means(label.x.npc = 'left', label.y.npc = 'top')

p2 <- ggplot(data = gbm.rtk.clonality.scna, aes(rtk_status, non_neutral_genome_frac)) +
 geom_boxplot(aes(fill = rtk_status), outlier.shape = NA) + geom_jitter(shape=16, position=position_jitter(0.2)) + 
 xlab('RTKs subtype') + ylab('Burden of SCNAs') + guides(fill = guide_legend(title= 'RTKs subtype')) + 
 theme(axis.text.x = element_text(angle=30, vjust=0.5)) + stat_compare_means(label.x.npc = 'left', label.y.npc = 'top')

p3 <- ggplot(data = gbm.rtk.clonality.scna, aes(rtk_status, clo_genome_frac)) +
 geom_boxplot(aes(fill = rtk_status), outlier.shape = NA) + geom_jitter(shape=16, position=position_jitter(0.2)) + 
 xlab('RTKs subtype') + ylab('Burden of clonal SCNAs') + guides(fill = guide_legend(title= 'RTKs subtype')) + 
 theme(axis.text.x = element_text(angle=30, vjust=0.5)) + stat_compare_means(label.x.npc = 'left', label.y.npc = 'top')

p4 <- ggplot(data = gbm.rtk.clonality.scna, aes(rtk_status, subclo_genome_frac)) +
 geom_boxplot(aes(fill = rtk_status), outlier.shape = NA) + geom_jitter(shape=16, position=position_jitter(0.2)) + 
 xlab('RTKs subtype') + ylab('Burden of subclonal SCNAs') + guides(fill = guide_legend(title= 'RTKs subtype')) + 
 theme(axis.text.x = element_text(angle=30, vjust=0.5)) + stat_compare_means(label.x.npc = 'left', label.y.npc = 'top')


ggsave(ggarrange(p1, p2, p3, p4, ncol=3, nrow=2, common.legend = TRUE), 
 file='/result/Section4/gbm_rtk_clonality_purity_scna.pdf')



