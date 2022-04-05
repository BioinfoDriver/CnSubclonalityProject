# @args 
#  infile = TCGA_mastercalls.abs_segtabs.fixed.txt
#  ploi: the sample’s ploidy (rounded to the nearest integer)
#
# @returns: a tibble
ClosubcloSegNum <- function(infile, ploi, ana.sams)
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


	# abs.seg.call <- subset(abs.seg.call, !(Chromosome == 1 & End < 123035434))
	# abs.seg.call <- subset(abs.seg.call, !(Chromosome == 19 & Start > 27681782))
	# abs.seg.call <- subset(abs.seg.call, !(Chromosome == 7 | Chromosome == 10))
	
	# filter
	abs.seg.call <- subset(abs.seg.call, Sample %in% paste0(ana.sams, '-01'))
	print(table(subset(abs.seg.call, del_seg == 1 | amp_seg == 1)$subclo_seg))
	
	# chr1p/19q deletion clonality
	p1q19.del.seg.call <- subset(abs.seg.call, (Chromosome == 1 & End < 123035434 & del_seg == 1) | 
	 (Chromosome == 19 & Start > 27681782 & del_seg == 1))
	
	print(table(p1q19.del.seg.call$subclo_seg))
	
	non.p1q19.del.seg.call <- subset(abs.seg.call, !((Chromosome == 1 & End < 123035434 & del_seg == 1) | 
	 (Chromosome == 19 & Start > 27681782 & del_seg ==1))) 
	non.p1q19.del.seg.call <- subset(non.p1q19.del.seg.call, del_seg == 1 | amp_seg == 1)
	
	print(table(non.p1q19.del.seg.call$subclo_seg))
	
	print(fisher.test(cbind(table(p1q19.del.seg.call$subclo_seg), 
	 table(non.p1q19.del.seg.call$subclo_seg)))$p.value)
	
	
	# chr7 amplication & chr10 deletion clonality
	chr7.amp.10.del.seg.call <- subset(abs.seg.call, (Chromosome == 7 & amp_seg == 1) | 
	 (Chromosome == 10 & del_seg == 1))
	
	print(table(chr7.amp.10.del.seg.call$subclo_seg))
	
	non.chr7.amp.10.del.seg.call <- subset(abs.seg.call, !((Chromosome == 7 & amp_seg == 1) | 
	 (Chromosome == 10 & del_seg ==1))) 
	non.chr7.amp.10.del.seg.call <- subset(non.chr7.amp.10.del.seg.call, del_seg == 1 | amp_seg == 1)
	
	print(table(non.chr7.amp.10.del.seg.call$subclo_seg))	
	
	print(fisher.test(cbind(table(chr7.amp.10.del.seg.call$subclo_seg), 
	 table(non.chr7.amp.10.del.seg.call$subclo_seg)))$p.value)
}


gli.puri.ploi <- readRDS(file='/data/tcga_gli_puri_ploi.rds')
in.file <- '/data/tcga_glioma_abs_seg.txt'

# gold set
gold.set <- readRDS(file='/data/gold_set.rds')
pan.glio.cli.mol.data <- readRDS(file='/data/tcga_glioma_cli_mol.rds')
pan.glio.cli.mol.data <- pan.glio.cli.mol.data[gold.set, ]

# all samples
ClosubcloSegNum(infile=in.file, ploi=gli.puri.ploi, 
 ana.sams=rownames(pan.glio.cli.mol.data))

# Astrocytoma,IDHmut
ClosubcloSegNum(infile=in.file, ploi=gli.puri.ploi, 
 ana.sams=rownames(subset(pan.glio.cli.mol.data, Integrated_Diagnoses == 'Astrocytoma,IDHmut')))

# Glioblastoma,IDHwt
ClosubcloSegNum(infile=in.file, ploi=gli.puri.ploi, 
 ana.sams=rownames(subset(pan.glio.cli.mol.data, Integrated_Diagnoses == 'Glioblastoma,IDHwt')))

# Oligodendroglioma,IDHmut-codel
ClosubcloSegNum(infile=in.file, ploi=gli.puri.ploi, 
 ana.sams=rownames(subset(pan.glio.cli.mol.data, Integrated_Diagnoses == 'Oligodendroglioma,IDHmut-codel')))


#########################################################Plot
{
library('ggpubr')

## chr1p/19q deletion
col.pal <- c("#E59398", "#D42527")
chr1p19q.all.stat <- data.frame(freq=c(0.821, 0.179, 0.689, 0.311), # p.value = 7.225271e-44
 group=rep(c('chr1p/19q del', 'other SCNAs'), each=2), clonality=rep(c('Clonal', 'Subclonal'), times=2))
chr1p19q.all.stat$clonality <- factor(chr1p19q.all.stat$clonality, levels=c('Subclonal', 'Clonal'))

chr1p19q.ast.stat <- data.frame(freq=c(0.563, 0.437, 0.611, 0.389), # p.value = 0.069572
 group=rep(c('chr1p/19q del', 'other SCNAs'), each=2), clonality=rep(c('Clonal', 'Subclonal'), times=2))
chr1p19q.ast.stat$clonality <- factor(chr1p19q.ast.stat$clonality, levels=c('Subclonal', 'Clonal'))

chr1p19q.gli.stat <- data.frame(freq=c(0.751, 0.249, 0.739, 0.261), # p.value = 0.5573164
 group=rep(c('chr1p/19q del', 'other SCNAs'), each=2), clonality=rep(c('Clonal', 'Subclonal'), times=2))
chr1p19q.gli.stat$clonality <- factor(chr1p19q.gli.stat$clonality, levels=c('Subclonal', 'Clonal'))

chr1p19q.oli.stat <- data.frame(freq=c(0.915, 0.085, 0.724, 0.276), # p.value = 2.609388e-50
 group=rep(c('chr1p/19q del', 'other SCNAs'), each=2), clonality=rep(c('Clonal', 'Subclonal'), times=2))
chr1p19q.oli.stat$clonality <- factor(chr1p19q.oli.stat$clonality, levels=c('Subclonal', 'Clonal'))


plot.chr1p19q.all <- ggbarplot(chr1p19q.all.stat, "group", "freq", fill = "clonality", color="clonality", 
 ylab='Propotion of alterations', xlab=FALSE, label=TRUE, lab.col="white", lab.pos="in", lab.size=2, 
 lab.vjust=2, palette=col.pal) + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(size=6))

plot.chr1p19q.ast <- ggbarplot(chr1p19q.ast.stat, "group", "freq", fill = "clonality", color="clonality", 
 ylab='Propotion of alterations', xlab=FALSE, label=TRUE, lab.col="white", lab.pos="in", lab.size=2, 
 lab.vjust=2, palette=col.pal) + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(size=6))

plot.chr1p19q.gli <- ggbarplot(chr1p19q.gli.stat , "group", "freq", fill = "clonality", color="clonality", 
 ylab='Propotion of alterations', xlab=FALSE, label=TRUE, lab.col="white", lab.pos="in", lab.size=2, 
 lab.vjust=2, palette=col.pal) + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(size=6))

plot.chr1p19q.oli <- ggbarplot(chr1p19q.oli.stat, "group", "freq", fill = "clonality", color="clonality", 
 ylab='Propotion of alterations', xlab=FALSE, label=TRUE, lab.col="white", lab.pos="in", lab.size=2, 
 lab.vjust=2, palette=col.pal) + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(size=6))


## entire chromosome 7 and loss of entire chromosome 10
chr7chr10.all.stat <- data.frame(freq=c(0.738, 0.262, 0.691, 0.309), # p.value = 4.456287e-11
 group=rep(c('chr7 amp & chr10 del', 'other SCNAs'), each=2), clonality=rep(c('Clonal', 'Subclonal'), times=2))
chr7chr10.all.stat$clonality <- factor(chr7chr10.all.stat$clonality, levels=c('Subclonal', 'Clonal'))

chr7chr10.ast.stat <- data.frame(freq=c(0.641, 0.359, 0.607, 0.393), # p.value = 0.03447926
 group=rep(c('chr7 amp & chr10 del', 'other SCNAs'), each=2), clonality=rep(c('Clonal', 'Subclonal'), times=2))
chr7chr10.ast.stat$clonality <- factor(chr7chr10.ast.stat$clonality, levels=c('Subclonal', 'Clonal'))

chr7chr10.gli.stat <- data.frame(freq=c(0.764, 0.236, 0.733, 0.267), # p.value = 0.0001477675
 group=rep(c('chr7 amp & chr10 del', 'other SCNAs'), each=2), clonality=rep(c('Clonal', 'Subclonal'), times=2))
chr7chr10.gli.stat$clonality <- factor(chr7chr10.gli.stat$clonality, levels=c('Subclonal', 'Clonal'))

chr7chr10.oli.stat <- data.frame(freq=c(0.766, 0.234, 0.794, 0.206), # p.value = 0.2693551
 group=rep(c('chr7 amp & chr10 del', 'other SCNAs'), each=2), clonality=rep(c('Clonal', 'Subclonal'), times=2))
chr7chr10.oli.stat$clonality <- factor(chr7chr10.oli.stat$clonality, levels=c('Subclonal', 'Clonal'))


plot.chr7chr10.all <- ggbarplot(chr7chr10.all.stat, "group", "freq", fill = "clonality", color="clonality", 
 ylab='Propotion of alterations', xlab=FALSE, label=TRUE, lab.col="white", lab.pos="in", lab.size=2, 
 lab.vjust=2, palette=col.pal) + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(size=6))

plot.chr7chr10.ast <- ggbarplot(chr7chr10.ast.stat, "group", "freq", fill = "clonality", color="clonality", 
 ylab='Propotion of alterations', xlab=FALSE, label=TRUE, lab.col="white", lab.pos="in", lab.size=2, 
 lab.vjust=2, palette=col.pal) + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(size=6))

plot.chr7chr10.gli <- ggbarplot(chr7chr10.gli.stat , "group", "freq", fill = "clonality", color="clonality", 
 ylab='Propotion of alterations', xlab=FALSE, label=TRUE, lab.col="white", lab.pos="in", lab.size=2, 
 lab.vjust=2, palette=col.pal) + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(size=6))

plot.chr7chr10.oli <- ggbarplot(chr7chr10.oli.stat, "group", "freq", fill = "clonality", color="clonality", 
 ylab='Propotion of alterations', xlab=FALSE, label=TRUE, lab.col="white", lab.pos="in", lab.size=2, 
 lab.vjust=2, palette=col.pal) + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(size=6))


ggsave(ggarrange(plot.chr1p19q.all, plot.chr1p19q.ast, plot.chr1p19q.gli, plot.chr1p19q.oli,
 plot.chr7chr10.all, plot.chr7chr10.ast, plot.chr7chr10.gli, plot.chr7chr10.oli, ncol=4, nrow=3, common.legend = TRUE), 
 file='result/Section1/arm_level_clonality.pdf')
}


# region
chr.region <- read.csv(file='/data/chromosome.band.hg19.txt', 
 header=TRUE, sep='\t', stringsAsFactors=FALSE)

chr.region <- subset(chr.region, !(X.chrom %in% c('chrX', 'chrY')))
chr.region$X.chrom <- gsub('chr', '', chr.region$X.chrom)
chr.region$arm <- substr(chr.region$name, 1, 1)
chr.region$index <- paste0(chr.region$X.chrom, chr.region$arm)

chr.region <- sapply(split(chr.region, factor(chr.region$index)), function(x){
 region <- c(x$X.chrom[1], x$chromStart[1],  x$chromEnd[length(x$chromEnd)], x$index[1])
 names(region) <- c('chrom', 'chromStart', 'chromEnd', 'arm')
 return(region)
})
chr.region <- as.data.frame(t(chr.region))

library(dplyr)
chr.region <- chr.region %>% mutate_at(1:3, as.numeric)
chr.region <- chr.region[order(chr.region$chrom), ]


ArmLevelclonality <- function(infile, ploi, ana.sams, arm.region)
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

	
	# filter
	abs.seg.call <- subset(abs.seg.call, Sample %in% paste0(ana.sams, '-01'))
		
	arm.level.clonality.stat <- lapply(seq(nrow(arm.region)), function(index){
	 each.arm <- arm.region[index, ]
	 if(substr(each.arm$arm, nchar(each.arm$arm)-1, nchar(each.arm$arm))=='p'){
	  arm.seg.call <- subset(abs.seg.call, Chromosome == each.arm$chrom & End <= each.arm$chromEnd)
	 
	 }else{
	  arm.seg.call <- subset(abs.seg.call, Chromosome == each.arm$chrom & Start >= each.arm$chromStart)
	 }
	
	 sub.genome.stat <- arm.seg.call %>% dplyr::select(Sample, Length, neutral_seg, amp_seg, del_seg, subclo_seg) %>% 
	  group_by(Sample) %>% mutate( 
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
	 
	 sub.genome.stat$arm <- each.arm$arm
	 return(sub.genome.stat)
	})
	
	arm.level.clonality.stat <- do.call(rbind, arm.level.clonality.stat)
	return(arm.level.clonality.stat)
}

arm.level.clonality <- ArmLevelclonality(infile=in.file, ploi=gli.puri.ploi, 
 ana.sams=rownames(pan.glio.cli.mol.data), arm.region=chr.region)


arm.clonal.amp <- as.data.frame(table(subset(arm.level.clonality, clo_amp_genome_frac>0.5)$arm))
arm.clonal.amp$cna_type <- 'Amp'
arm.clonal.amp$clonality <- 'Clonal'

arm.subclonal.amp <- as.data.frame(table(subset(arm.level.clonality, subclo_amp_genome_frac>0.5)$arm))
arm.subclonal.amp$cna_type <- 'Amp'
arm.subclonal.amp$clonality <- 'Subclonal'

arm.clonal.del <- as.data.frame(table(subset(arm.level.clonality, clo_del_genome_frac>0.5)$arm))
arm.clonal.del$cna_type <- 'Del'
arm.clonal.del$clonality <- 'Clonal'

arm.subclonal.del <- as.data.frame(table(subset(arm.level.clonality, subclo_del_genome_frac>0.5)$arm))
arm.subclonal.del$cna_type <- 'Del'
arm.subclonal.del$clonality <- 'Subclonal'


# arm水平显著拷贝数改变
# https://gdac.broadinstitute.org/runs/analyses__2016_01_28/reports/cancer/GBMLGG-TP/CopyNumber_Gistic2/nozzle.html

amp.list <- c('1p','1q','7p','7q','10p','19p','19q','20p','20q')
del.list <- c('1p','4p','4q','6p','6q','9p','9q','10p','10q','11p','11q','12q','13q','14q','15q','18p','18q','19q', '22q')



arm.het <- rbind(merge(subset(arm.clonal.amp, Var1 %in% amp.list), subset(arm.subclonal.amp, Var1 %in% amp.list), by='Var1'),
merge(subset(arm.clonal.del, Var1 %in% del.list), subset(arm.subclonal.del, Var1 %in% del.list), by='Var1'))

p.values <- lapply(seq(nrow(arm.het)), function(i){
		
	clo.n <- arm.het[i, 'Freq.x']
	alt.n <- arm.het[i, 'Freq.x'] + arm.het[i, 'Freq.y']
	
	clo.perc <- clo.n/alt.n
	p.value <- binom.test(x=clo.n, n=alt.n, p = 0.5, alternative = "two.sided", conf.level = 0.95)$p.value
	   
	return(c(clo.perc, p.value))
})

p.values <- do.call(rbind, p.values)
colnames(p.values) <- c('clo.alt.perc', 'p.value')
p.values <- as.data.frame(p.values)
p.values$fdr <- p.adjust(p=p.values$p.value, method = 'fdr')
p.values$arm <- arm.het$Var1
p.values$cna_type <- arm.het$cna_type.x

{
   # clo.alt.perc      p.value          fdr arm cna_type
# 1     0.5267857 6.367998e-01 6.877438e-01 19p      Amp
# 2     0.5789474 4.176922e-01 5.126222e-01 19q      Amp
# 3     0.6129032 2.810415e-01 3.993748e-01  1p      Amp
# 4     0.5789474 4.176922e-01 5.126222e-01  1q      Amp
# 5     0.7175573 6.717866e-07 3.023040e-06 20p      Amp
# 6     0.6111111 1.336742e-01 2.623743e-01 20q      Amp
# 7     0.7398844 1.405659e-19 9.488196e-19  7p      Amp
# 8     0.6337449 3.651282e-05 1.232308e-04  7q      Amp
# 9     0.8181818 4.675186e-35 1.262300e-33 10p      Del
# 10    0.7692308 5.026162e-17 2.714127e-16 10q      Del
# 11    0.3870968 2.810415e-01 3.993748e-01 11p      Del
# 12    0.3448276 1.360459e-01 2.623743e-01 11q      Del
# 13    0.4761905 8.776143e-01 9.062943e-01 12q      Del
# 14    0.6294118 9.187496e-04 2.480624e-03 13q      Del
# 15    0.6645161 5.120483e-05 1.536145e-04 14q      Del
# 16    0.4597701 5.202916e-01 6.107771e-01 15q      Del
# 17    0.5781250 2.604355e-01 3.993748e-01 18p      Del
# 18    0.5862069 2.370471e-01 3.993748e-01 18q      Del
# 19    0.7773438 1.292216e-19 9.488196e-19 19q      Del
# 20    0.9186992 3.089822e-23 4.171259e-22  1p      Del
# 21    0.5546875 2.504399e-01 3.993748e-01 22q      Del
# 22    0.5138889 9.062943e-01 9.062943e-01  4p      Del
# 23    0.5505618 3.965702e-01 5.126222e-01  4q      Del
# 24    0.6236559 2.201859e-02 5.404563e-02  6p      Del
# 25    0.7073171 4.848718e-06 1.870220e-05  6q      Del
# 26    0.6101695 1.174774e-01 2.623743e-01  9p      Del
# 27    0.4411765 6.075914e-01 6.835403e-01  9q      Del
}


#########################################################Plot
sig.arm.level.clonality <- rbind(subset(rbind(arm.clonal.amp, arm.subclonal.amp), Var1 %in% amp.list), 
 subset(rbind(arm.clonal.del, arm.subclonal.del), Var1 %in% del.list)) 

# Number
sig.arm.level.clonality <- subset(sig.arm.level.clonality, !(Var1 == '10p' & cna_type == 'Amp'))
sig.arm.level.clonality$Var1 <- paste(sig.arm.level.clonality$cna_type, '(', sig.arm.level.clonality$Var1, ')', sep='')
sig.arm.level.clonality$cna_type <- paste(sig.arm.level.clonality$cna_type, sig.arm.level.clonality$clonality, sep='_')

# Frequency
sig.arm.clo.freq <- arm.het
sig.arm.clo.freq$Freq.x <- sig.arm.clo.freq$Freq.x/(sig.arm.clo.freq$Freq.x + sig.arm.clo.freq$Freq.y)
sig.arm.clo.freq$Freq.y <- 1-sig.arm.clo.freq$Freq.x
tmp1 <- sig.arm.clo.freq[, 1:4]
colnames(tmp1) <- c('Var1', 'Freq', 'cna_type', 'clonality')
tmp2 <- sig.arm.clo.freq[, c(1, 5:7)]
colnames(tmp2) <- c('Var1', 'Freq', 'cna_type', 'clonality')
sig.arm.clo.freq <- rbind(tmp1, tmp2)

sig.arm.clo.freq$Var1 <- paste(sig.arm.clo.freq$cna_type, '(', sig.arm.clo.freq$Var1, ')', sep='')
sig.arm.clo.freq$cna_type <- paste(sig.arm.clo.freq$cna_type, sig.arm.clo.freq$clonality, sep='_')

# sort
lab.order <- c("Amp(1p)", "Del(1p)", "Amp(1q)", "Del(4p)", "Del(4q)", "Del(6p)", "Del(6q)",
"Amp(7p)", "Amp(7q)", "Del(9p)", "Del(9q)", "Del(10p)", "Del(10q)", "Del(11p)", "Del(11q)",
"Del(12q)", "Del(13q)", "Del(14q)", "Del(15q)", "Del(18p)", "Del(18q)",
 "Amp(19p)", "Amp(19q)", "Del(19q)", "Amp(20p)", "Amp(20q)", "Del(22q)")

sig.arm.level.clonality$Var1 <- factor(sig.arm.level.clonality$Var1, levels=lab.order)
sig.arm.clo.freq$Var1 <- factor(sig.arm.clo.freq$Var1, levels=lab.order)
sig.arm.level.clonality$cna_type <- factor(sig.arm.level.clonality$cna_type, 
 levels=c('Amp_Subclonal', 'Amp_Clonal', 'Del_Subclonal', 'Del_Clonal'))
sig.arm.clo.freq$cna_type <- factor(sig.arm.clo.freq$cna_type, 
 levels=c('Amp_Subclonal', 'Amp_Clonal', 'Del_Subclonal', 'Del_Clonal'))


col.pal <- c("#E59398", "#D42527", "#8CB6D2", "#2171A9")
p1 <- ggbarplot(sig.arm.level.clonality, "Var1", "Freq", fill = "cna_type", color = "cna_type", 
 ylab = FALSE, xlab = FALSE, label = TRUE, lab.col = "white", lab.pos = "in", lab.size = 1, lab.vjust=2, palette = col.pal, 
 legend = 'none', x.text.angle=30) + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(size=3))

p2 <- ggbarplot(sig.arm.clo.freq, "Var1", "Freq", fill = "cna_type", color = "cna_type", 
 ylab = FALSE, xlab = FALSE, label = TRUE, lab.col = "white", lab.pos = "in", lab.size = 1, lab.vjust=2, palette = col.pal, 
 legend = 'none', x.text.angle=30) + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(size=3))


ggsave(ggarrange(p1, p2, ncol=1, nrow=3, common.legend = TRUE), file='/result/Section1/arm_level_clonal_freq.pdf')


