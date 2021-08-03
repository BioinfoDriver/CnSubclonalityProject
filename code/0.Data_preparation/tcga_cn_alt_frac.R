# @args 
#  infile = TCGA_mastercalls.abs_segtabs.fixed.txt
#  ploi: the sample’s ploidy (rounded to the nearest integer)
#
# @returns: a tibble

CNAltFraction <- function(infile, ploi)
{	
	library('dplyr')

	abs.seg.call <- read.table(file = infile, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
	abs.seg.call <- subset(abs.seg.call, !is.na(Chromosome) & Chromosome != 23)

	isect.samples <- intersect(abs.seg.call$Sample, rownames(ploi))
	ploi <- ploi[isect.samples, ]
	abs.seg.call <- subset(abs.seg.call, Sample %in% isect.samples)
	abs.seg.call$ploidy <- ploi[match(abs.seg.call$Sample, rownames(ploi)), "ploidy"]

	IsNeutralSeg <- function(x) as.integer(x$Modal_Total_CN == round(x$ploidy))
	abs.seg.call$neutral_seg <- IsNeutralSeg(abs.seg.call)

	abs.seg.call$amp_seg <- as.integer(abs.seg.call$Modal_Total_CN > round(abs.seg.call$ploidy))
	abs.seg.call$del_seg <- as.integer(abs.seg.call$Modal_Total_CN < round(abs.seg.call$ploidy))

	IsSubClonal <- function(x) as.integer(x$Subclonal_HSCN_a1 | x$Subclonal_HSCN_a2)
	abs.seg.call$subclo_seg <- IsSubClonal(abs.seg.call)



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
in.file <- '/data/tcga_glioma_abs_seg.txt'
gli.cn.alt.frac <- CNAltFraction(infile=in.file, ploi=gli.puri.ploi)


# 2 samples without SCNAs, and 3 samples with extreme SCNAs (>90%) thus were excluded
gli.cn.alt.frac <- subset(gli.cn.alt.frac, !(non_neutral_genome_frac > 0.9 | non_neutral_genome_frac == 0))
saveRDS(gli.cn.alt.frac, file='/data/tcga_gli_cn_alt_frac.rds') 


