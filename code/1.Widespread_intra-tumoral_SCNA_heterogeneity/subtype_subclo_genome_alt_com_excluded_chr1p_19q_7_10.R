# @args 
#  infile = TCGA_mastercalls.abs_segtabs.fixed.txt
#  ploi: the sample’s ploidy (rounded to the nearest integer)
#
# @returns: a tibble

# Excluding chromosome 1p, 19q, 7 and 10 
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
	abs.seg.call <- subset(abs.seg.call, !(Chromosome == 7 | Chromosome == 10))
	
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

gli.puri.ploi <- readRDS(file='/data/tcga_gli_puri_ploi.rds')
in.file <- '/data/OriginalData/tcga_glioma_abs_seg.txt'
gli.cn.alt.frac <- CNAltFraction(infile=in.file, ploi=gli.puri.ploi)


# 10 samples without SCNAs, and 3 samples with extreme SCNAs (>90%) thus were excluded
gli.cn.alt.frac <- subset(gli.cn.alt.frac, !(non_neutral_genome_frac > 0.9 | non_neutral_genome_frac == 0))


# gold set
gold.set <- readRDS( file='/data/gold_set.rds')
pan.glio.cli.mol.data <- readRDS(file='/data/tcga_glioma_cli_mol.rds')

gold.set <- intersect(rownames(gli.cn.alt.frac), paste0(gold.set, '-01'))
gli.cn.alt.frac <- cbind(gli.cn.alt.frac[gold.set, ], pan.glio.cli.mol.data[substr(gold.set, 1, 12), ])

# subclonal genome stat
gli.cn.alt.frac$Integrated_Diagnoses <- factor(gli.cn.alt.frac$Integrated_Diagnoses, 
 levels = c('Oligodendroglioma,IDHmut-codel', 'Astrocytoma,IDHmut', 'Glioblastoma,IDHwt'))


# > aggregate(subclo_cn_alt_frac ~ Integrated_Diagnoses, data = gli.cn.alt.frac, summary)
            # Integrated_Diagnoses subclo_cn_alt_frac.Min. subclo_cn_alt_frac.1st Qu.
# 1 Oligodendroglioma,IDHmut-codel              0.00000000                 0.00000000
# 2             Astrocytoma,IDHmut              0.00000000                 0.17908652
# 3             Glioblastoma,IDHwt              0.00000000                 0.03621721
  # subclo_cn_alt_frac.Median subclo_cn_alt_frac.Mean subclo_cn_alt_frac.3rd Qu. subclo_cn_alt_frac.Max.
# 1                0.07091772              0.32884289                 0.71958943              1.00000000
# 2                0.50879726              0.48509075                 0.74873027              1.00000000
# 3                0.28879292              0.35262991                 0.56957415              1.00000000



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
SubtypeComByFeacture(features, subtype='Integrated_Diagnoses', 
 labels=c('Oligodendroglioma,IDHmut-codel', 'Astrocytoma,IDHmut', 'Glioblastoma,IDHwt'))
# p-value = 1.12e-21
# p-value = 7.71e-19
# p-value = 8.60e-20
# p-value = 4.08e-09
           # Z            P   P.adjusted                                         comparisons               var1_var2
# 1  -2.808449 2.489036e-03 2.489036e-03             Astrocytoma,IDHmut - Glioblastoma,IDHwt non_neutral_genome_frac
# 2   7.075503 7.445380e-13 1.116807e-12 Astrocytoma,IDHmut - Oligodendroglioma,IDHmut-codel non_neutral_genome_frac
# 3   9.797065 5.795196e-23 1.738559e-22 Glioblastoma,IDHwt - Oligodendroglioma,IDHmut-codel non_neutral_genome_frac
# 4   2.395923 8.289286e-03 8.289286e-03             Astrocytoma,IDHmut - Glioblastoma,IDHwt      subclo_genome_frac
# 5   8.930327 2.123751e-19 6.371253e-19 Astrocytoma,IDHmut - Oligodendroglioma,IDHmut-codel      subclo_genome_frac
# 6   7.393630 7.143701e-14 1.071555e-13 Glioblastoma,IDHwt - Oligodendroglioma,IDHmut-codel      subclo_genome_frac
# 7  -5.116918 1.552844e-07 2.329266e-07             Astrocytoma,IDHmut - Glioblastoma,IDHwt         clo_genome_frac
# 8   4.627916 1.846821e-06 1.846821e-06 Astrocytoma,IDHmut - Oligodendroglioma,IDHmut-codel         clo_genome_frac
# 9   9.153051 2.767396e-20 8.302187e-20 Glioblastoma,IDHwt - Oligodendroglioma,IDHmut-codel         clo_genome_frac
# 10  4.673842 1.478081e-06 2.217122e-06             Astrocytoma,IDHmut - Glioblastoma,IDHwt      subclo_cn_alt_frac
# 11  5.784808 3.629750e-09 1.088925e-08 Astrocytoma,IDHmut - Oligodendroglioma,IDHmut-codel      subclo_cn_alt_frac
# 12  2.176707 1.475123e-02 1.475123e-02 Glioblastoma,IDHwt - Oligodendroglioma,IDHmut-codel      subclo_cn_alt_frac
  
      

##################################################plot
setwd('/result/Section1/')
library('ggpubr')

comparison <- combn(x = c('Oligodendroglioma,IDHmut-codel', 'Astrocytoma,IDHmut', 'Glioblastoma,IDHwt'), m = 2, simplify = FALSE)
shape.s <- c(21, 24, 22)[as.numeric(gli.cn.alt.frac$Integrated_Diagnoses)]
color.s <- c("#FC4E07", "#E7B800", "#00AFBB")[as.numeric(gli.cn.alt.frac$Integrated_Diagnoses)]
 
subclo.bur.plot <- ggboxplot(data=gli.cn.alt.frac, x='Integrated_Diagnoses', y='subclo_genome_frac', xlab = FALSE, 
 ylab = 'Subclonal SCNAs Burden', legend = 'none', alpha = 1.0, add = 'jitter', palette=, 
 add.params=list(color=color.s, shape=shape.s, size = 0.8), outlier.shape=NA) + theme(aspect.ratio=1)

subclo.bur.plot <- subclo.bur.plot + stat_compare_means(comparisons = comparison) + 
 stat_compare_means(label.x.npc = 'left', label.y.npc = 'top')

clo.bur.plot <- ggboxplot(data=gli.cn.alt.frac, x='Integrated_Diagnoses', y='clo_genome_frac', xlab = FALSE, 
 ylab = 'Clonal SCNAs Burden', legend = 'none', alpha = 1.0, add = 'jitter', palette=c("#FC4E07", "#E7B800", "#00AFBB"), 
 add.params=list(color=color.s, shape=shape.s, size = 0.8), outlier.shape=NA) + theme(aspect.ratio=1)

clo.bur.plot <- clo.bur.plot + stat_compare_means(comparisons = comparison) + 
 stat_compare_means(label.x.npc = 'left', label.y.npc = 'top')

scna.bur.plot <- ggboxplot(data=gli.cn.alt.frac, x='Integrated_Diagnoses', y='non_neutral_genome_frac', xlab = FALSE, 
 ylab = 'SCNAs Burden', legend = 'none', alpha = 1.0, add = 'jitter', palette=c("#FC4E07", "#E7B800", "#00AFBB"), 
 add.params=list(color=color.s, shape=shape.s, size = 0.8), outlier.shape=NA) + theme(aspect.ratio=1)

scna.bur.plot <- scna.bur.plot + stat_compare_means(comparisons = comparison) + 
 stat_compare_means(label.x.npc = 'left', label.y.npc = 'top')

subclo.perc.plot <- ggboxplot(data=gli.cn.alt.frac, x='Integrated_Diagnoses', y='subclo_cn_alt_frac', xlab = FALSE, 
 ylab = 'Subclonal SCNAs Percentage', legend = 'none', alpha = 1.0, add = 'jitter', palette=c("#FC4E07", "#E7B800", "#00AFBB"), 
 add.params=list(color=color.s, shape=shape.s, size = 0.8), outlier.shape=NA) + theme(aspect.ratio=1)

subclo.perc.plot <- subclo.perc.plot + stat_compare_means(comparisons = comparison) + 
 stat_compare_means(label.x.npc = 'left', label.y.npc = 'top')
  
ggsave(ggarrange(scna.bur.plot, clo.bur.plot, subclo.bur.plot, subclo.perc.plot,
 ncol=2, nrow=2), filename='excluded_chr1p_19q_7_10_scna_burden.pdf')


