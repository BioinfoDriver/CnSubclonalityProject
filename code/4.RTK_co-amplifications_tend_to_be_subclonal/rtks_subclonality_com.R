
rtks <- read.csv(file='/data/RTKs.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)

gold.set <- readRDS('/data/gold_set.rds')
gene.het <- readRDS(file='/data/gene.het.rds')
gene.het <- gene.het[, paste(gold.set, '01', sep = '-')]

tcga.cli.data <- readRDS('/data/tcga_glioma_cli_mol.rds')



RtkNonRtkAmpClonalityCom <- function(het.mat, rtk, cli.data, subtype, grade){
 
 cli.data <- subset(cli.data, IDH_CODEL_SUBTYPE %in% subtype & histological_grade %in% grade)
 rtk.het <- het.mat[rtk, intersect(colnames(het.mat), paste0(rownames(cli.data), '-01'))]
 
 print(table(rtk.het))
 
 non.rtk.het <- het.mat[setdiff(rownames(het.mat), rtk), intersect(colnames(het.mat), paste0(rownames(cli.data), '-01'))]
 print(table(non.rtk.het))
 
 mat <- c(table(rtk.het)['1'], table(rtk.het)['2'], table(non.rtk.het)['1'], table(non.rtk.het)['2'])
 mat <- matrix(data=mat, nrow=2, byrow=FALSE)
 print(mat)
 
 p.value <- fisher.test(mat)$p.value
 return(p.value)
}

RtkNonRtkAmpClonalityCom(gene.het, rtks$Approved.symbol, tcga.cli.data, 
	c('IDHmut-codel', 'IDHmut-non-codel', 'IDHwt'), c('G2', 'G3', 'G4', NA))
RtkNonRtkAmpClonalityCom(gene.het, rtks$Approved.symbol, tcga.cli.data, c('IDHwt'), c('G4'))



MulRtkUniRtkAmpClonalityCom <- function(het.mat, rtk, cli.data, subtype, grade){
 
 cli.data <- subset(cli.data, IDH_CODEL_SUBTYPE %in% subtype & histological_grade %in% grade)
 
 rtk.het <- het.mat[rtk, intersect(colnames(het.mat), paste0(rownames(cli.data), '-01'))]


 mul.rtk.het <- rtk.het[, colSums(rtk.het >= 1) >1 ]
 uni.rtk.het <- rtk.het[, colSums(rtk.het >= 1) ==1 ]

 print(table(mul.rtk.het))
 print(table(uni.rtk.het)) 
 
 mat <- c(table(mul.rtk.het)['1'], table(mul.rtk.het)['2'], table(uni.rtk.het)['1'], table(uni.rtk.het)['2'])
 mat <- matrix(data=mat, nrow=2, byrow=FALSE)
 print(mat)
 
 p.value <- fisher.test(mat)$p.value
 return(p.value)
}

MulRtkUniRtkAmpClonalityCom(gene.het, rtks$Approved.symbol, tcga.cli.data, 
	c('IDHmut-codel', 'IDHmut-non-codel', 'IDHwt'), c('G2', 'G3', 'G4', NA))
MulRtkUniRtkAmpClonalityCom(gene.het, rtks$Approved.symbol, tcga.cli.data, c('IDHwt'), c('G4'))



###############################################################################

RtkAmplificationFrequency <- function(het.mat, rtk, cli.data, subtype, grade){
 
 cli.data <- subset(cli.data, IDH_CODEL_SUBTYPE %in% subtype & histological_grade %in% grade)
 rtk.het <- het.mat[rtk, intersect(colnames(het.mat), paste0(rownames(cli.data), '-01'))]
 rtk.het[rtk.het < 0] <- 0
 
 rtk.alt.pat <- sum(colSums(rtk.het >= 1) >= 1)
 total.pat <- ncol(rtk.het)
 rtk.alt.fre <- rtk.alt.pat/total.pat
 
 rtk.alt.stat <- lapply(seq(nrow(rtk.het)), function(x){
	alt <- rtk.het[x, ]
	alt.pat <- sum(alt != 0)
	alt.clo.pat <- sum(alt == 1)
	alt.sub.pat <- sum(alt == 2)
	
	alt.fre <- alt.pat/total.pat
	clo.fre <- alt.clo.pat/total.pat
	sub.fre <- alt.sub.pat/total.pat
	sub.perc <- alt.sub.pat/alt.pat
	
	stat <- c(alt.pat, alt.clo.pat, alt.sub.pat, alt.fre, clo.fre, sub.fre, sub.perc)
	
	names(stat) <- c('alt.n', 'clo.alt.n', 'sub.alt.n', 'alt.freq', 'clo.alt.freq', 'sub.alt.freq', 'sub.alt.perc')
	return(stat)
 })
 rtk.alt.stat <- do.call(rbind, rtk.alt.stat)
 rownames(rtk.alt.stat) <- rownames(rtk.het)
 
 rtk.alt.stat <- as.data.frame(rtk.alt.stat)
 rtk.alt.stat <- rtk.alt.stat[order(rtk.alt.stat$alt.n, decreasing=TRUE), ]
 
 rtk.alt.stat$total.alt.n <- rtk.alt.pat
 rtk.alt.stat$total.alt.freq <- rtk.alt.fre
 rtk.alt.stat$n.patient <- total.pat
 
 return(rtk.alt.stat)
}

RtkAmplificationFrequency(gene.het, rtks$Approved.symbol, tcga.cli.data, 
	c('IDHmut-codel', 'IDHmut-non-codel', 'IDHwt'), c('G2', 'G3', 'G4', NA))
RtkAmplificationFrequency(gene.het, rtks$Approved.symbol, tcga.cli.data, c('IDHwt'), c('G4'))



#########################################################Plot
col.pal <- c("#D42527", "#E59398")
all.rtk.vs.nonrtk.stat <- data.frame(freq=c(0.552, 0.448, 0.466, 0.534), group=rep(c('RTK', 'Non-RTK'), each=2), 
	clonality=rep(c('High-level clonal amplification', 'High-level subclonal amplification'), times=2))

gbm.rtk.vs.nonrtk.stat <- data.frame(freq=c(0.551, 0.449, 0.480, 0.520), group=rep(c('RTK', 'Non-RTK'), each=2), 
	clonality=rep(c('High-level clonal amplification', 'High-level subclonal amplification'), times=2))

plot.rtkvsnonrtk.all <- ggbarplot(all.rtk.vs.nonrtk.stat, "group", "freq", fill = "clonality", color="clonality", 
 ylab='Propotion of amplifications', xlab=FALSE, label=TRUE, lab.col="white", lab.pos="in", lab.size=2, 
 lab.vjust=2, palette=col.pal) + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(size=6))

plot.rtkvsnonrtk.gbm <- ggbarplot(gbm.rtk.vs.nonrtk.stat, "group", "freq", fill = "clonality", color="clonality", 
 ylab='Propotion of amplifications', xlab=FALSE, label=TRUE, lab.col="white", lab.pos="in", lab.size=2, 
 lab.vjust=2, palette=col.pal) + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(size=6))



all.mul.uni.rtk.stat <- data.frame(freq=c(0.524, 0.476, 0.728, 0.272), group=rep(c('co-amp', 'non-co-amp'), each=2), 
	clonality=rep(c('High-level clonal amplification', 'High-level subclonal amplification'), times=2))

gbm.mul.uni.rtk.stat <- data.frame(freq=c(0.505, 0.495, 0.833, 0.167), group=rep(c('co-amp', 'non-co-amp'), each=2), 
	clonality=rep(c('High-level clonal amplification', 'High-level subclonal amplification'), times=2))

plot.mulvsunirtk.all <- ggbarplot(all.mul.uni.rtk.stat, "group", "freq", fill = "clonality", color="clonality", 
 ylab='Propotion of amplifications', xlab=FALSE, label=TRUE, lab.col="white", lab.pos="in", lab.size=2, 
 lab.vjust=2, palette=col.pal) + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(size=6))

plot.mulvsunirtk.gbm <- ggbarplot(gbm.mul.uni.rtk.stat, "group", "freq", fill = "clonality", color="clonality", 
 ylab='Propotion of amplifications', xlab=FALSE, label=TRUE, lab.col="white", lab.pos="in", lab.size=2, 
 lab.vjust=2, palette=col.pal) + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(size=6))



ggsave(ggarrange(plot.rtkvsnonrtk.all, plot.rtkvsnonrtk.gbm, 
 plot.mulvsunirtk.all, plot.mulvsunirtk.gbm, ncol=2, nrow=2, common.legend = TRUE), 
 file='result/Section4/rtk_clonality_compare.pdf')



