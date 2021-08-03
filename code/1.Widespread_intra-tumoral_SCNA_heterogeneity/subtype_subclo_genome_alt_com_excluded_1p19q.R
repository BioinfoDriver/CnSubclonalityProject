# @args 
#  infile = TCGA_mastercalls.abs_segtabs.fixed.txt
#  ploi: the sample’s ploidy (rounded to the nearest integer)
#
# @returns: a tibble

# Excluding 1p and 19q
CNAltFraction <- function(infile, ploi)
{	
	library('dplyr')

	abs.seg.call <- read.table(file = infile, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
	abs.seg.call <- subset(abs.seg.call, !is.na(Chromosome) & Chromosome != 23)

	isect.samples <- intersect(abs.seg.call$Sample, rownames(ploi))
	ploi <- ploi[isect.samples, ]
	abs.seg.call <- subset(abs.seg.call, Sample %in% isect.samples)
	abs.seg.call$ploidy <- ploi[match(abs.seg.call$Sample, rownames(ploi)), "ploidy"]

	IsNeutralSeg <- function(x) as.integer(abs.seg.call$Modal_Total_CN == round(abs.seg.call$ploidy))
	abs.seg.call$neutral_seg <- IsNeutralSeg(abs.seg.call)

	abs.seg.call$amp_seg <- as.integer(abs.seg.call$Modal_Total_CN > round(abs.seg.call$ploidy))
	abs.seg.call$del_seg <- as.integer(abs.seg.call$Modal_Total_CN < round(abs.seg.call$ploidy))

	IsSubClonal <- function(x) as.integer(x$Subclonal_HSCN_a1 | x$Subclonal_HSCN_a2)
	abs.seg.call$subclo_seg <- IsSubClonal(abs.seg.call)

	
	abs.seg.call <- subset(abs.seg.call, !(Chromosome == 1 & End < 123035434))
	abs.seg.call <- subset(abs.seg.call, !(Chromosome == 19 & Start > 27681782))
	
	sub.genome.stat <- abs.seg.call %>% dplyr::select(Sample, Length, neutral_seg, amp_seg, del_seg, subclo_seg) %>% 
	 group_by(Sample) %>% mutate(
	 neutral_genome_frac = sum(Length[neutral_seg == 1])/sum(Length),
	 non_neutral_genome_frac = sum(Length[neutral_seg == 0])/sum(Length),
	 
	 subclo_genome_frac = sum(Length[neutral_seg == 0 & subclo_seg == 1])/sum(Length),
	 clo_genome_frac = sum(Length[neutral_seg == 0 & subclo_seg == 0])/sum(Length),

	 subclo_amp_genome_frac = sum(Length[amp_seg == 1 & subclo_seg == 1])/sum(Length),
	 clo_amp_genome_frac = sum(Length[amp_seg == 1 & subclo_seg == 0])/sum(Length),
	 subclo_del_genome_frac = sum(Length[del_seg == 1 & subclo_seg == 1])/sum(Length),
	 clo_del_genome_frac = sum(Length[del_seg == 1 & subclo_seg == 0])/sum(Length),

	 subclo_cn_alt_frac = sum(Length[neutral_seg == 0 & subclo_seg == 1])/sum(Length[neutral_seg == 0]),
	 clo_cn_alt_frac = sum(Length[neutral_seg == 0 & subclo_seg == 0])/sum(Length[neutral_seg == 0]),

	 subclo_amp_alt_frac = sum(Length[amp_seg == 1 & subclo_seg == 1])/sum(Length[neutral_seg == 0]),
	 clo_amp_alt_frac = sum(Length[amp_seg == 1 & subclo_seg == 0])/sum(Length[neutral_seg == 0]),
	 subclo_del_alt_frac = sum(Length[del_seg == 1 & subclo_seg == 1])/sum(Length[neutral_seg == 0]),
	 clo_del_alt_frac = sum(Length[del_seg == 1 & subclo_seg == 0])/sum(Length[neutral_seg == 0]))

	sub.genome.stat <- sub.genome.stat[!duplicated(sub.genome.stat$Sample), ] %>% 
	 dplyr::select(!c(Length, neutral_seg, amp_seg, del_seg, subclo_seg))
	sub.genome.stat <- as.data.frame(sub.genome.stat)
	rownames(sub.genome.stat) <- sub.genome.stat$Sample
	sub.genome.stat <- sub.genome.stat[, -1]
	
	return(sub.genome.stat)
}

gli.puri.ploi <- readRDS(file='/pub5/xiaoyun/Jobs/J22/CopyNumberClonalityProject/Resource/CuratedData/tcga_gli_puri_ploi.rds')
in.file <- '/pub5/xiaoyun/Jobs/J22/CopyNumberClonalityProject/Resource/CuratedData/tcga_glioma_abs_seg.txt'
gli.cn.alt.frac <- CNAltFraction(infile=in.file, ploi=gli.puri.ploi)


# 5 samples without SCNAs, and 3 samples with extreme SCNAs (>90%) thus were excluded
gli.cn.alt.frac <- subset(gli.cn.alt.frac, !(non_neutral_genome_frac > 0.9 | non_neutral_genome_frac == 0))


# gold set
gold.set <- readRDS( file='/pub5/xiaoyun/Jobs/J22/CopyNumberClonalityProject/Resource/CuratedData/gold_set.rds')
pan.glio.cli.mol.data <- readRDS(file='/pub5/xiaoyun/Jobs/J22/CopyNumberClonalityProject/Resource/CuratedData/tcga_glioma_cli_mol.rds')

gold.set <- intersect(rownames(gli.cn.alt.frac), paste0(gold.set, '-01'))
gli.cn.alt.frac <- cbind(gli.cn.alt.frac[gold.set, ], pan.glio.cli.mol.data[substr(gold.set, 1, 12), ])

# subclonal genome stat
gli.cn.alt.frac$IDH_CODEL_SUBTYPE <- factor(gli.cn.alt.frac$IDH_CODEL_SUBTYPE, 
 levels = c('IDHmut-codel', 'IDHmut-non-codel', 'IDHwt'))


# > aggregate(subclo_cn_alt_frac ~ IDH_CODEL_SUBTYPE, data = gli.cn.alt.frac, summary)
  # IDH_CODEL_SUBTYPE subclo_cn_alt_frac.Min. subclo_cn_alt_frac.1st Qu.
# 1      IDHmut-codel              0.00000000                 0.00000000
# 2  IDHmut-non-codel              0.00000000                 0.24447929
# 3             IDHwt              0.00000000                 0.02800200
  # subclo_cn_alt_frac.Median subclo_cn_alt_frac.Mean subclo_cn_alt_frac.3rd Qu.	subclo_cn_alt_frac.Max.
# 1                0.09789246              0.33132051                 0.71801804	1              1.00000000
# 2                0.47042759              0.48045137                 0.72510763	2              1.00000000
# 3                0.17853358              0.27701768                 0.41893041	3              1.00000000


########################################## with or without 1p/19q
cn.alt.frac <- readRDS(file='/pub5/xiaoyun/Jobs/J22/CopyNumberClonalityProject/Resource/CuratedData/gli_glod_cn_alt_frac.rds')
gli.cn.alt.frac$all_clo_genome_frac <- cn.alt.frac[rownames(gli.cn.alt.frac), 'clo_genome_frac']
gli.cn.alt.frac$all_clo_cn_alt_frac <- cn.alt.frac[rownames(gli.cn.alt.frac), 'clo_cn_alt_frac']


codel.cn.alt.frac <- subset(gli.cn.alt.frac, IDH_CODEL_SUBTYPE == 'IDHmut-codel')
wilcox.test(codel.cn.alt.frac$all_clo_genome_frac, codel.cn.alt.frac$clo_genome_frac, paired=TRUE)$p.value
# 3.759611e-26
wilcox.test(codel.cn.alt.frac$all_clo_cn_alt_frac, codel.cn.alt.frac$clo_cn_alt_frac, paired=TRUE)$p.value
# 1.915642e-11



# subtype associated with genome alt fraction
SubtypeSubcloStat <- function(dat, alt.feac, sub.feac, label){
 
 dat <- dat[!is.na(dat[, sub.feac]), ]
 dat[, grep(sub.feac, colnames(dat))] <- factor(dat[, sub.feac], levels = label)
 
 p.value <- kruskal.test(get(alt.feac) ~ get(sub.feac), data = dat)$p.value
 res <- dunn.test::dunn.test(dat[, alt.feac], dat[, sub.feac], method = 'bh')
 res <- dplyr::bind_rows(res[2:length(res)]) %>% mutate(var1_var2 = paste(alt.feac, sub.feac, sep = " - "))
 
 return(list(p.value, res))
}

SubtypeComByFeacture <- function(features, subtype, labels){
	
 res.list <- lapply(features, function(feature){
  res <- SubtypeSubcloStat(gli.cn.alt.frac, alt.feac=feature, sub.feac=subtype, label=labels)
  return(res)
 })
 pvalues <- setNames(sapply(res.list, function(x) x[[1]]), features)
 com.pvalues <- sapply(res.list, function(x) x[2])
 com.pvalues <- do.call(rbind, com.pvalues)

  return(list(pvalues, com.pvalues))
}


features <- c('non_neutral_genome_frac', 'subclo_genome_frac', 'clo_genome_frac', 'subclo_cn_alt_frac')
SubtypeComByFeacture(features, subtype='IDH_CODEL_SUBTYPE', labels=c('IDHwt', 'IDHmut-non-codel', 'IDHmut-codel'))
# p-value = 4.67e-69
# p-value = 9.29e-22
# p-value = 5.37e-64
# p-value = 4.83e-15
             # Z            P   P.adjusted                     comparisons	var1_var2
# 1   -6.8262941 4.356811e-12 4.356811e-12 IDHmut-codel - IDHmut-non-codel	non_neutral_genome_
# 2  -16.7396220 3.370746e-63 1.011224e-62            IDHmut-codel - IDHwt	non_neutral_genome_
# 3  -11.2889556 7.439712e-30 1.115957e-29        IDHmut-non-codel - IDHwt	non_neutral_genome_
# 4   -9.2666817 9.601351e-21 2.880405e-20 IDHmut-codel - IDHmut-non-codel	subclo_genome_frac
# 5   -8.6600771 2.357227e-18 3.535841e-18            IDHmut-codel - IDHwt	subclo_genome_frac
# 6    1.3976888 8.110327e-02 8.110327e-02        IDHmut-non-codel - IDHwt	subclo_genome_frac
# 7   -4.3829554 5.854005e-06 5.854005e-06 IDHmut-codel - IDHmut-non-codel	clo_genome_frac
# 8  -15.1984144 1.811421e-52 5.434264e-52            IDHmut-codel - IDHwt	clo_genome_frac
# 9  -12.5399665 2.255930e-36 3.383895e-36        IDHmut-non-codel - IDHwt	clo_genome_frac
# 10  -6.3439143 1.119997e-10 1.679996e-10 IDHmut-codel - IDHmut-non-codel	subclo_cn_alt_frac
# 11  -0.4725798 3.182565e-01 3.182565e-01            IDHmut-codel - IDHwt	subclo_cn_alt_frac
# 12   7.4442611 4.874427e-14 1.462328e-13        IDHmut-non-codel - IDHwt	subclo_cn_alt_frac




##################################################plot
setwd('/pub5/xiaoyun/Jobs/J22/CopyNumberClonalityProject/Results/Section1/Results/Exclude1P19qAnalysis')
library('ggpubr')
clo.frac.plot <- ggpaired(codel.cn.alt.frac, cond1="all_clo_genome_frac", cond2="clo_genome_frac", 
 color="condition", line.color="gray", palette="npg", xlab=FALSE, ylab="Clonal SCNAs Burden") + theme(aspect.ratio=1)

clo.perc.plot <- ggpaired(codel.cn.alt.frac, cond1="all_clo_cn_alt_frac", cond2="clo_cn_alt_frac", 
 color="condition", line.color="gray", palette="npg", xlab=FALSE, ylab="Clonal SCNAs Percentage") + theme(aspect.ratio=1)

ggsave(ggarrange(clo.frac.plot, clo.perc.plot, ncol=2, nrow=2), filename='codel_clo_genome_com.pdf')


# gli.cn.alt.frac <- subset(gli.cn.alt.frac, histological_grade %in% c('G2' ,'G3'))
comparison <- combn(x = c('IDHmut-codel', 'IDHmut-non-codel', 'IDHwt'), m = 2, simplify = FALSE)
shape.s <- c(21, 24, 22)[as.numeric(gli.cn.alt.frac$IDH_CODEL_SUBTYPE)]
color.s <- c("#FC4E07", "#E7B800", "#00AFBB")[as.numeric(gli.cn.alt.frac$IDH_CODEL_SUBTYPE)]
 
subclo.bur.plot <- ggboxplot(data=gli.cn.alt.frac, x='IDH_CODEL_SUBTYPE', y='subclo_genome_frac', xlab = FALSE, 
 ylab = 'Subclonal SCNAs Burden', legend = 'none', alpha = 1.0, add = 'jitter', palette=, 
 add.params=list(color=color.s, shape=shape.s, size = 0.8), outlier.shape=NA) + theme(aspect.ratio=1)

subclo.bur.plot <- subclo.bur.plot + stat_compare_means(comparisons = comparison) + 
 stat_compare_means(label.x.npc = 'left', label.y.npc = 'top')

clo.bur.plot <- ggboxplot(data=gli.cn.alt.frac, x='IDH_CODEL_SUBTYPE', y='clo_genome_frac', xlab = FALSE, 
 ylab = 'Clonal SCNAs Burden', legend = 'none', alpha = 1.0, add = 'jitter', palette=c("#FC4E07", "#E7B800", "#00AFBB"), 
 add.params=list(color=color.s, shape=shape.s, size = 0.8), outlier.shape=NA) + theme(aspect.ratio=1)

clo.bur.plot <- clo.bur.plot + stat_compare_means(comparisons = comparison) + 
 stat_compare_means(label.x.npc = 'left', label.y.npc = 'top')

scna.bur.plot <- ggboxplot(data=gli.cn.alt.frac, x='IDH_CODEL_SUBTYPE', y='non_neutral_genome_frac', xlab = FALSE, 
 ylab = 'SCNAs Burden', legend = 'none', alpha = 1.0, add = 'jitter', palette=c("#FC4E07", "#E7B800", "#00AFBB"), 
 add.params=list(color=color.s, shape=shape.s, size = 0.8), outlier.shape=NA) + theme(aspect.ratio=1)

scna.bur.plot <- scna.bur.plot + stat_compare_means(comparisons = comparison) + 
 stat_compare_means(label.x.npc = 'left', label.y.npc = 'top')

subclo.perc.plot <- ggboxplot(data=gli.cn.alt.frac, x='IDH_CODEL_SUBTYPE', y='subclo_cn_alt_frac', xlab = FALSE, 
 ylab = 'Subclonal SCNAs Percentage', legend = 'none', alpha = 1.0, add = 'jitter', palette=c("#FC4E07", "#E7B800", "#00AFBB"), 
 add.params=list(color=color.s, shape=shape.s, size = 0.8), outlier.shape=NA) + theme(aspect.ratio=1)

subclo.perc.plot <- subclo.perc.plot + stat_compare_means(comparisons = comparison) + 
 stat_compare_means(label.x.npc = 'left', label.y.npc = 'top')
  
ggsave(ggarrange(scna.bur.plot, clo.bur.plot, subclo.bur.plot, subclo.perc.plot,
 ncol=2, nrow=2), filename='excluded_1p19q_scna_burden.pdf')


