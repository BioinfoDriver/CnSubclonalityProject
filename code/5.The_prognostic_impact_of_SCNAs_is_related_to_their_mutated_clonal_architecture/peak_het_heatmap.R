
PeakHetPlot <- function(het.mat, mole.subtype, subtype.peaks, outfile){
	library(ggpubr)
	peaks <- rownames(het.mat)
	
	het.mat <- as.data.frame(t(het.mat))
	het.mat <- cbind(het.mat, mole.subtype[rownames(het.mat), 'IDH_CODEL_SUBTYPE', FALSE])
	het.mat <- data.table::melt(het.mat, id = 'IDH_CODEL_SUBTYPE', na.rm = TRUE)
	het.mat$label <- paste(het.mat$IDH_CODEL_SUBTYPE, het.mat$variable, sep = '_')


	het.mat <- split(het.mat, het.mat$label)
	het.mat <- lapply(het.mat, function(x){
		
		count <- dplyr::count(x, value)
		count <- cbind(subtype = x[, 1][1:nrow(count)], variable = x[, 2][1:nrow(count)], count)
	 
		count <- subset(count, value != 0)
		count$alt_freq <- as.numeric(sprintf("%0.2f", count$n/sum(count$n)))

		return(count)
	})
	het.mat <- do.call(rbind, het.mat)


	het.mat <- lapply(c('IDHmut-codel', 'IDHmut-non-codel', 'IDHwt'), function(subgroup){
		return(subset(het.mat, subtype==subgroup & variable %in% subtype.peaks[[subgroup]]))
	
	})
	het.mat <- do.call(rbind, het.mat)
	
	
	for(subgroup in c('IDHmut-codel', 'IDHmut-non-codel', 'IDHwt')){
		if(length(peaks) > length(subtype.peaks[[subgroup]])){
					
			mis.label <- setdiff(peaks, subtype.peaks[[subgroup]])
			mis.mat <- data.frame(subtype=subgroup, variable=mis.label, value=1, n=0, alt_freq=0)
			
			het.mat <- rbind(het.mat, mis.mat)
		}
	}
	

	het.mat$value <- factor(het.mat$value, levels = c(2, 1, -2, -1))
	col.pal <- c("#E59398", "#D42527", "#8CB6D2", "#2171A9")

	pdf(file=outfile, width=7, height=7)
	
	p1 <- ggbarplot(subset(het.mat, subtype=='IDHmut-codel'), "variable", "alt_freq", fill = "value", color = "value", 
	 ylab = FALSE, xlab = FALSE, label = TRUE, lab.col = "white", lab.pos = "in", lab.size = 1, lab.vjust=2, palette = col.pal, 
	 legend = 'none', x.text.angle=30) + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(size=3))
 
	p2 <- ggbarplot(subset(het.mat, subtype=='IDHmut-non-codel'), "variable", "alt_freq", fill = "value", color = "value", 
	 ylab = FALSE, xlab = FALSE, label = TRUE, lab.col = "white", lab.pos = "in", lab.size = 1, lab.vjust=2, palette = col.pal, 
	 legend = 'none', x.text.angle=30) + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(size=3))
 
	p3 <- ggbarplot(subset(het.mat, subtype=='IDHwt'), "variable", "alt_freq", fill = "value", color = "value", 
	 ylab = FALSE, xlab = FALSE, label = TRUE, lab.col = "white", lab.pos = "in", lab.size = 1, lab.vjust=2, palette = col.pal, 
	 legend = 'none', x.text.angle=30) + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(size=3))
 
	print(ggarrange(p1, p2, p3, heights=c(0.1, 0.1, 0.1, 0.7), ncol=1, nrow=4))
 

	p4 <- ggbarplot(subset(het.mat, subtype=='IDHmut-codel'), "variable", "n", fill = "value", color = "value", ylab = FALSE, 
	 xlab = FALSE, label = TRUE, lab.col = "white", lab.pos = "in", lab.size = 1, lab.vjust=2, palette = col.pal, 
	 legend = 'none', x.text.angle=30) + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(size=3))

	p5 <- ggbarplot(subset(het.mat, subtype=='IDHmut-non-codel'), "variable", "n", fill = "value", color = "value", ylab = FALSE, 
	 xlab = FALSE, label = TRUE, lab.col = "white", lab.pos = "in", lab.size = 1, lab.vjust=2, palette = col.pal, 
	 legend = 'none', x.text.angle=30) + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(size=3)) + 
	 scale_y_continuous(limits=c(0, 100))

	p6 <- ggbarplot(subset(het.mat, subtype=='IDHwt'), "variable", "n", fill = "value", color = "value", ylab = FALSE, 
	 xlab = FALSE, label = TRUE, lab.col = "white", lab.pos = "in", lab.size = 1, lab.vjust=2, palette = col.pal, 
	 legend = 'none', x.text.angle=30) + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(size=3))

	print(ggarrange(p4, p5, p6, heights=c(0.15, 0.15, 0.15, 0.55), ncol=1, nrow=4))
	
	dev.off()

}


gold.set <- readRDS('/data/gold_set.rds')
wgistic.het <- readRDS('/data/wgistic_het.rds')
gli.cli.mol.data <- readRDS(file='/data/tcga_glioma_cli_mol.rds')
subtype.peaks <- readRDS(file='/data/high_fre_alt_peaks.rds')


gli.cli.mol.data <- gli.cli.mol.data[gold.set, ]
rownames(gli.cli.mol.data) <- paste(rownames(gli.cli.mol.data), '01', sep = '-')
wgistic.het <- wgistic.het[, paste(gold.set, '01', sep = '-')]


PeakHetPlot(het.mat=wgistic.het, mole.subtype=gli.cli.mol.data, subtype.peaks, '/result/Section5/wpeak_het_heatmap.pdf')
