# @args 
#  infile = TCGA_mastercalls.abs_segtabs.fixed.txt
#  ploi: the sample’s ploidy (rounded to the nearest integer)
#
# @returns: a tibble

# Excluding chromosome 7 and chromosome 10 
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


# 3 samples without SCNAs, and 3 samples with extreme SCNAs (>90%) thus were excluded
gli.cn.alt.frac <- subset(gli.cn.alt.frac, !(non_neutral_genome_frac > 0.9 | non_neutral_genome_frac == 0))


# gold set
gold.set <- readRDS( file='/data/gold_set.rds')
pan.glio.cli.mol.data <- readRDS(file='/data/tcga_glioma_cli_mol.rds')

gold.set <- intersect(rownames(gli.cn.alt.frac), paste0(gold.set, '-01'))
gli.cn.alt.frac <- cbind(gli.cn.alt.frac[gold.set, ], pan.glio.cli.mol.data[substr(gold.set, 1, 12), ])



########################################## with or without 1p/19q
cn.alt.frac <- readRDS(file='/data/gli_glod_cn_alt_frac.rds')
gli.cn.alt.frac$all_clo_genome_frac <- cn.alt.frac[rownames(gli.cn.alt.frac), 'clo_genome_frac']
gli.cn.alt.frac$all_clo_cn_alt_frac <- cn.alt.frac[rownames(gli.cn.alt.frac), 'clo_cn_alt_frac']


codel.cn.alt.frac <- subset(gli.cn.alt.frac, Integrated_Diagnoses == 'Glioblastoma,IDHwt')
wilcox.test(codel.cn.alt.frac$all_clo_genome_frac, codel.cn.alt.frac$clo_genome_frac, paired=TRUE)$p.value
# 5.151425e-54
wilcox.test(codel.cn.alt.frac$all_clo_cn_alt_frac, codel.cn.alt.frac$clo_cn_alt_frac, paired=TRUE)$p.value
# 3.023056e-18


##################################################plot
setwd('/result/Section1')
library('ggpubr')
clo.frac.plot <- ggpaired(codel.cn.alt.frac, cond1="all_clo_genome_frac", cond2="clo_genome_frac", 
 color="condition", line.color="gray", palette="npg", xlab=FALSE, ylab="Clonal SCNAs Burden") + theme(aspect.ratio=1)

clo.perc.plot <- ggpaired(codel.cn.alt.frac, cond1="all_clo_cn_alt_frac", cond2="clo_cn_alt_frac", 
 color="condition", line.color="gray", palette="npg", xlab=FALSE, ylab="Clonal SCNAs Percentage") + theme(aspect.ratio=1)

ggsave(ggarrange(clo.frac.plot, clo.perc.plot, ncol=2, nrow=2), filename='chr7_10_clo_genome_com.pdf')


