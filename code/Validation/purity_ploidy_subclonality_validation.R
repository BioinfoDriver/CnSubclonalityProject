
# Purity and Ploidy Comparisons
glioma.puri.ploi <- readRDS(file='/data/pcawg.glioma.puri.ploi.rds')
glioma.cna.data <- readRDS(file='/data/pcawg.glioma.cna.data.rds')


# ABSOLUTE-based tumour purity
setwd('/data/OriginalData/')
abs.call.data <- read.csv(file='TCGA_mastercalls.abs_tables_JSedit.fixed.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
abs.puri.ploi <- subset(abs.call.data, substr(array, 1, 12) %in% glioma.puri.ploi$submitted_donor_id)
abs.puri.ploi <- subset(abs.puri.ploi, substr(array, 14, 15) == '01' & call.status == 'called')

abs.puri.ploi <- abs.puri.ploi[, c('array', 'purity', 'ploidy')]
abs.puri.ploi$array <- substr(abs.puri.ploi$array, 1, 12)
colnames(abs.puri.ploi) <- c('array', 'purity.ABSOLUTE', 'ploidy.ABSOLUTE')


pcawg.puri.ploi <- glioma.puri.ploi[, c('submitted_donor_id', 'purity', 'ploidy')]
colnames(pcawg.puri.ploi) <- c('donor_id', 'purity.PCAWG', 'ploidy.PCAWG')

abs.pcawg.puri.ploi <- merge(abs.puri.ploi, pcawg.puri.ploi, by.x='array', by.y='donor_id')

CorTest <- function(puri.ploi){
	puri.sta <- cor.test(puri.ploi[, 2], puri.ploi[, 4], method="spearman", exact=FALSE)
	ploi.sta <- cor.test(puri.ploi[, 3], puri.ploi[, 5], method="spearman", exact=FALSE)
	
	sta <- c(puri.sta$estimate, puri.sta$p.value, ploi.sta$estimate, ploi.sta$p.value)
	names(sta) <- c('purity cor', 'purity pvalue', 'ploidy cor', 'ploidy pvalue')

	return(sta)
}

# > CorTest(abs.pcawg.puri.ploi)
   # purity cor purity pvalue    ploidy cor ploidy pvalue 
 # 8.626471e-01  5.198813e-14  8.723716e-01  1.227740e-14 



# Subclonal genome fraction
donor_id <- rep(names(glioma.cna.data), lapply(glioma.cna.data, nrow))
glioma.cna.data <- do.call(rbind, glioma.cna.data)
glioma.cna.data$donor_id <- donor_id
glioma.cna.data$ploidy <- glioma.puri.ploi$ploidy[match(glioma.cna.data$donor_id, glioma.puri.ploi$submitted_donor_id)]
glioma.cna.data <- subset(glioma.cna.data, !is.na(battenberg_nMaj1_A))


SWGCNAltFraction <- function(abs.seg.call){	
	library('dplyr')

	abs.seg.call <- subset(abs.seg.call, chromosome != 'X')

	IsNeutralSeg <- function(dat){
	  
	 neutral <- sapply(seq(nrow(dat)), function(index){
	  x <- dat[index, , FALSE]
	 
	  if(is.na(x$battenberg_nMaj2_A)){
	   nt <- as.integer(x$battenberg_nMaj1_A + x$battenberg_nMin1_A == round(x$ploidy))
	  }else{
	   nt <- as.integer(x$battenberg_nMaj1_A + x$battenberg_nMin1_A == round(x$ploidy) &
	    x$battenberg_nMaj2_A + x$battenberg_nMin2_A == round(x$ploidy))
	  }
	 })
	 return(neutral)
	}
	abs.seg.call$neutral_seg <- IsNeutralSeg(abs.seg.call)

	IsSubClonal <- function(x) as.integer(!is.na(x$battenberg_nMaj2_A))
	abs.seg.call$subclo_seg <- IsSubClonal(abs.seg.call)
    
	abs.seg.call$Length <- abs.seg.call$start - abs.seg.call$end

	sub.genome.stat <- abs.seg.call %>% dplyr::select(donor_id, Length, neutral_seg, subclo_seg) %>% 
	 group_by(donor_id) %>% mutate(
	 neutral_genome_frac = sum(Length[neutral_seg == 1])/sum(Length),
	 non_neutral_genome_frac = sum(Length[neutral_seg == 0])/sum(Length),
	 
	 subclo_genome_frac = sum(Length[neutral_seg == 0 & subclo_seg == 1])/sum(Length),
	 clo_genome_frac = sum(Length[neutral_seg == 0 & subclo_seg == 0])/sum(Length),

	 subclo_cn_alt_frac = sum(Length[neutral_seg == 0 & subclo_seg == 1])/sum(Length[neutral_seg == 0]),
	 clo_cn_alt_frac = sum(Length[neutral_seg == 0 & subclo_seg == 0])/sum(Length[neutral_seg == 0]))

	sub.genome.stat <- sub.genome.stat[!duplicated(sub.genome.stat$donor_id), ] %>% 
	 dplyr::select(!c(Length, neutral_seg, subclo_seg))
	sub.genome.stat <- as.data.frame(sub.genome.stat)
	rownames(sub.genome.stat) <- sub.genome.stat$donor_id
	sub.genome.stat <- sub.genome.stat[, -1]
	
	return(sub.genome.stat)
}

pcawg.cn.alt.frac <- SWGCNAltFraction(glioma.cna.data)
rownames(pcawg.cn.alt.frac) <- paste0(rownames(pcawg.cn.alt.frac), '-01')

abs.cn.alt.frac <- readRDS(file='/data/gli_glod_cn_alt_frac.rds')

int.sect <- intersect(rownames(pcawg.cn.alt.frac), rownames(abs.cn.alt.frac))
pcawg.cn.alt.frac <- pcawg.cn.alt.frac[int.sect, c('subclo_genome_frac', 'clo_genome_frac')]
colnames(pcawg.cn.alt.frac) <- c('subclonality.WGS', 'clonality.WGS')

abs.cn.alt.frac <- abs.cn.alt.frac[int.sect, c('subclo_genome_frac', 'clo_genome_frac')]
colnames(abs.cn.alt.frac) <- c('subclonality.SNP', 'clonality.SNP')

abs.pcawg.scna.clonality <- cbind(abs.cn.alt.frac, pcawg.cn.alt.frac)


# cor.test(~subclonality.WGS+subclonality.SNP, data=abs.pcawg.scna.clonality, method="spearman")
# # p-value = 6.437e-05, rho = 0.6017432
# cor.test(~clonality.WGS+clonality.SNP, data=abs.pcawg.scna.clonality, method="spearman")
# # p-value = 5.763e-07, rho = 0.7391399


library('ggpubr')
pdf('/result/Validation/puri_ploi_subclonality_snp_was.pdf')

p1 <- ggscatter(data=abs.pcawg.puri.ploi, x='purity.ABSOLUTE', y='purity.PCAWG', 
 xlab='Purity(SNP-array)', ylab='Purity(WGS)', add = "reg.line", cor.coef=TRUE, 
 cor.coeff.args=list(method="spearman", label.x.npc="left", label.y.npc="top"), size=0.9, alpha=0.5) + 
 scale_x_continuous(breaks = seq(0, 1.0, 0.2), limit=c(0, 1.0)) + scale_y_continuous(breaks = seq(0, 1.0, 0.2))

p2 <- ggscatter(data=abs.pcawg.puri.ploi, x='ploidy.ABSOLUTE', y='ploidy.PCAWG', 
 xlab='Ploidy(SNP-array)', ylab='Ploidy(WGS)', add = "reg.line", cor.coef=TRUE, 
 cor.coeff.args=list(method="spearman", label.x.npc="left", label.y.npc="top"), size=0.9,  alpha=0.5) + 
 scale_x_continuous(breaks = 1:5, limit=c(1, 5)) + scale_y_continuous(breaks = 1:5, limit=c(1, 5))

p3 <- ggscatter(data=abs.pcawg.scna.clonality, x='subclonality.SNP', y='subclonality.WGS', 
 xlab='Subclonality(SNP-array)', ylab='Subclonality(WGS)', add = "reg.line", cor.coef=TRUE, 
 cor.coeff.args=list(method="spearman", label.x.npc="left", label.y.npc="top"), size=0.9, alpha=0.5) + 
 scale_x_continuous(breaks = seq(0, 0.8, 0.2), limit=c(0, 0.8)) + scale_y_continuous(breaks = seq(0, 0.8, 0.2), limit=c(0, 0.8))

p4 <- ggscatter(data=abs.pcawg.scna.clonality, x='clonality.SNP', y='clonality.WGS', 
 xlab='Clonality(SNP-array)', ylab='Clonality(WGS)', add = "reg.line", cor.coef=TRUE, 
 cor.coeff.args=list(method="spearman", label.x.npc="left", label.y.npc="top"), size=0.9, alpha=0.5) + 
 scale_x_continuous(breaks = seq(0, 0.6, 0.2), limit=c(0, 0.6)) + scale_y_continuous(breaks = seq(0, 0.6, 0.2), limit=c(0, 0.6))

ggarrange(p1, p2, p3, p4, ncol=2, nrow=2)
dev.off()





