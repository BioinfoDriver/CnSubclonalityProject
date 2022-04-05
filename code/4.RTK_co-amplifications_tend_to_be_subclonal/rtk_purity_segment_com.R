
##########################################Purity comparison
# load data
rtks <- read.csv(file='/data/RTKs.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)

gold.set <- readRDS('/data/gold_set.rds')
gene.het <- readRDS(file='/data/gene.het.rds')
gene.het <- gene.het[, paste(gold.set, '01', sep = '-')]

# The patients of RTK clonal or subclonal amplification
rtk.het <- gene.het[rtks$Approved.symbol, ]


rtk.clo.sam <- colnames(rtk.het[, colSums(rtk.het == 1) >= 1])
rtk.sub.sam <- colnames(rtk.het[, colSums(rtk.het == 2) >= 1])

clo.sam <- colnames(gene.het[, colSums(gene.het == 1) >= 1])
sub.sam <- colnames(gene.het[, colSums(gene.het == 2) >= 1])


rtk.alt.sam <- unique(c(rtk.clo.sam, rtk.sub.sam))
rtk.non.alt.sam <- setdiff(unique(c(clo.sam, sub.sam)), rtk.alt.sam)
# rtk.non.alt.sam <- setdiff(colnames(gene.het[, colSums(gene.het !=0) >= 1]), rtk.alt.sam)

# ABSOLUTE-based tumour purity
abs.puri.ploi <- readRDS('/data/tcga_gli_puri_ploi.rds')
rtk.purity <- abs.puri.ploi[c(rtk.alt.sam, rtk.non.alt.sam), ]
rtk.purity$rtk_status <- rep(c('RTKs', 'Non-RTKs'), times=c(length(rtk.alt.sam), length(rtk.non.alt.sam)))

# purity comparison
wilcox.test(purity ~ rtk_status, rtk.purity, alternative = 'two.sided')$p.value # 0.7829576


library('ggplot2')
rtk.purity.plot <- ggplot(data = rtk.purity, aes(rtk_status, purity)) +
 geom_boxplot(aes(fill = rtk_status), outlier.shape = NA) + geom_jitter(shape=16, position=position_jitter(0.2)) + 
 xlab('RTKs subtype') + ylab('Tumor purity') + guides(fill = guide_legend(title= 'RTKs subtype')) + 
 theme(axis.text.x = element_text(angle=30, vjust=0.5)) + stat_compare_means(label.x.npc = 'left', label.y.npc = 'top')

ggsave(filename='rtk_purity_com.pdf', plot=ggarrange(rtk.purity.plot, ncol=3, nrow=2, common.legend = TRUE), 
 path='/result/Section4')



########################################segment size
gene.region <- readRDS(file='/data/gene.region.rds')
gli.puri.ploi <- readRDS(file='/data/tcga_gli_puri_ploi.rds')
in.file <- '/data/OriginalData/tcga_glioma_abs_seg.txt'


library('dplyr')
abs.seg.call <- read.table(file = in.file, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
abs.seg.call <- subset(abs.seg.call, Sample %in% paste(gold.set, '01', sep = '-'))

abs.seg.call <- subset(abs.seg.call, !is.na(Chromosome) & Chromosome != 23)

isect.samples <- intersect(abs.seg.call$Sample, rownames(gli.puri.ploi))
gli.puri.ploi <- gli.puri.ploi[isect.samples, ]
abs.seg.call <- subset(abs.seg.call, Sample %in% isect.samples)
abs.seg.call$ploidy <- gli.puri.ploi[match(abs.seg.call$Sample, rownames(gli.puri.ploi)), "ploidy"]
abs.seg.call$high_amp <- ifelse(abs.seg.call$Modal_Total_CN >= round(abs.seg.call$ploidy) + 2, TRUE, FALSE)

high.amp.seg <- subset(abs.seg.call, high_amp)
rel.cols <- c("Sample", "Chromosome", "Start", "End", "Length", 'Subclonal_HSCN_a1', 'Subclonal_HSCN_a2')
high.amp.seg <- high.amp.seg[, rel.cols]
high.amp.seg$Chromosome <- paste0('chr', high.amp.seg$Chromosome)
high.amp.seg$subclonal <- as.integer(high.amp.seg$Subclonal_HSCN_a1 | high.amp.seg$Subclonal_HSCN_a2)
high.amp.seg$index <- seq(nrow(high.amp.seg))
high.amp.seg$label <- apply(high.amp.seg[, c("Sample", "Chromosome", "Start", "End")], 1, 
	function(x) paste0(x, collapse='_'))

library('GenomicRanges')
gr <- makeGRangesFromDataFrame(high.amp.seg, keep.extra.columns=TRUE) 

rtk.seg.size <- overlapRegions(gr, subset(gene.region, hgnc_symbol %in% rtks$Approved.symbol), 
 colA=c('Sample', 'Length', 'index', 'label', 'subclonal'), colB=c('hgnc_symbol'))

non.rtk.seg.size <- overlapRegions(gr, subset(gene.region, !(hgnc_symbol %in% rtks$Approved.symbol)), 
 colA=c('Sample', 'Length', 'index', 'label', 'subclonal'), colB=c('hgnc_symbol'))


# the sizes of segments
rtk.segment.size <- data.frame(sizes=c(rtk.seg.size$Length[!duplicated(rtk.seg.size$index)], 
 non.rtk.seg.size$Length[!duplicated(non.rtk.seg.size$index)]), rtk_status=rep(c('RTK', 'Non-RTKs'), 
 times=c(sum(!duplicated(rtk.seg.size$index)), sum(!duplicated(non.rtk.seg.size$index)))))

wilcox.test(sizes ~ rtk_status, rtk.segment.size, alternative = 'two.sided')$p.value # 1.807867e-34


# plot
rtk.size.plot <- ggplot(data = rtk.segment.size, aes(rtk_status, sizes)) +
 geom_boxplot(aes(fill = rtk_status), outlier.shape = NA) + geom_jitter(shape=16, position=position_jitter(0.2)) + 
 xlab('RTKs subtype') + ylab('Size of segment') + guides(fill = guide_legend(title= 'RTKs subtype')) + 
 theme(axis.text.x = element_text(angle=30, vjust=0.5)) + stat_compare_means(label.x.npc = 'left', label.y.npc = 'top')

ggsave(filename='rtk_size_com.pdf', plot=ggarrange(rtk.size.plot, ncol=3, nrow=2, common.legend = TRUE), 
 path='/result/Section4')


#################################permutation test
abs.diff <- abs(median(rtk.seg.size$Length[!duplicated(rtk.seg.size$index)])-non.rtk.seg.size$Length)
prob[abs.diff > quantile(abs.diff, 0.12)] <- 0.000001
prob[quantile(abs.diff, 0.09) <=abs.diff & abs.diff < quantile(abs.diff, 0.12)] <- 0.0001
prob[quantile(abs.diff, 0.06) <=abs.diff & abs.diff < quantile(abs.diff, 0.09)] <- 0.01
prob[quantile(abs.diff, 0.03) <=abs.diff & abs.diff < quantile(abs.diff, 0.06)] <- 10
prob[abs.diff < quantile(abs.diff, 0.03)] <- 10000

non.rtk.seg.size$prob <- prob

set.seed(100)
ran.subclo.perc <- sapply(1:10000, function(x){
 
 id <- sample(x=unique(non.rtk.seg.size$index), size=length(unique(rtk.seg.size$index)), 
  replace = FALSE, prob=non.rtk.seg.size$prob[!duplicated(non.rtk.seg.size$index)])
 tmp <- subset(non.rtk.seg.size, index %in% id)
 
 p.value <- wilcox.test(tmp$Length, rtk.seg.size$Length[!duplicated(rtk.seg.size$index)], alternative = 'two.sided')$p.value
 print(p.value)
 if(p.value > 0.05){
	
	tmp <- subset(non.rtk.seg.size, label %in% tmp$label)
	return(sum(tmp$subclonal)/nrow(tmp))
 
 }
 else{
	return(NA)
 }
})

# > table(na.omit(ran.subclo.perc)[1:1000]>0.454)

# FALSE  TRUE 
    # 9   991
	
col.pal <- c("#E59398", "#D42527")
stat <- data.frame(Freq=c(991, 9), type=c('More', 'Less'))
p1 <- ggbarplot(stat, "type", "Freq", fill = "type", color = "type", ylab = 'Frequency', xlab = FALSE, 
 label = TRUE, lab.col = "white", lab.pos = "out", lab.size = 1, lab.vjust=2, palette = col.pal, 
 legend = 'none', x.text.angle=30) + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(size=3))

ggsave(filename='rtk_random_clonality_com.pdf', plot=ggarrange(p1, ncol=3, nrow=2, common.legend = TRUE), 
 path='/result/Section4')

