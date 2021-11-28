
library(survival)
ClonalitySurvPvalue <- function(clonality.data, surv.data, surv.type, grade, subtype){
	
	surv.data <- subset(surv.data, histological_grade %in% grade & Integrated_Diagnoses  %in% subtype)
	
	surv.pvalue <- lapply(colnames(clonality.data), function(peak){

		surv.data <- merge(surv.data, clonality.data[, peak, FALSE], by = 'row.names')
		surv.data <- surv.data[, c('Row.names', surv.type, paste0(surv.type, '_time'), peak, 
					  'age', 'gender', 'histological_grade')]
	
		colnames(surv.data) <- c('Patient_ID', 'event', 'time', 'sample.label', 'age', 'gender', 'grade')
		surv.data <- subset(surv.data, !(is.na(event) | is.na(time) | is.na(sample.label)))
		
		surv.data$sample.label[surv.data$sample.label == 0] <- 'wt'
		surv.data$sample.label[surv.data$sample.label %in% c(1, -1)] <- 'clonal'
		surv.data$sample.label[surv.data$sample.label %in% c(2, -2)] <- 'subclonal'
		
		# filter
		count <- table(surv.data$sample.label)
		count <- count[intersect(c('clonal', 'subclonal', 'wt'), names(count))]
		subtype <- names(count)[count >= 10]
		
		combn <- c('clonal wt', 'subclonal wt', 'clonal subclonal')
		pvalue.list <- list(c(NA, NA, NA), c(NA, NA, NA), c(NA, NA, NA))
		names(pvalue.list) <- combn
		
		if(length(subtype) >= 2){
			subtype.list <- combn(x = subtype, m = 2, simplify = FALSE)
			names(subtype.list) <- sapply(subtype.list, paste0, collapse = ' ')
		
			pvalue <- sapply(subtype.list, function(label){

				sur.data <- subset(surv.data, sample.label %in% label)
				if(all(label %in% c('wt', 'clonal')))
					sur.data$sample.label <- factor(sur.data$sample.label, levels = c('wt', 'clonal'))
				
				else if(all(label %in% c('wt', 'subclonal')))
					sur.data$sample.label <- factor(sur.data$sample.label, levels = c('wt', 'subclonal'))
				else
					sur.data$sample.label <- factor(sur.data$sample.label, levels = c('clonal', 'subclonal'))
				
				lr.test <- survdiff(Surv(time, event)~sample.label, data=sur.data)
				lr.p <- 1 - pchisq(lr.test$chisq, length(lr.test$n) - 1)
					
				univ.model <- coxph(Surv(time, event)~sample.label, data = sur.data)
				univ.p <- summary(univ.model)$coefficients[1, 5]
				
				sur.data <- subset(sur.data, !(is.na(age) | is.na(gender) | is.na(grade)))
				
				if(length(unique(sur.data$grade)) > 1)
					multiv.model <- coxph(Surv(time, event)~sample.label + age + gender + grade, 
					 data = sur.data, control = coxph.control(50))
				else
					multiv.model <- coxph(Surv(time, event)~sample.label + age + gender, 
					 data = sur.data, control = coxph.control(50))				
				 
				multiv.p <- summary(multiv.model)$coefficients[1, 5]
			
				return(c(lr.p, univ.p, multiv.p))
			}, simplify = FALSE)
			
			pvalue.list[names(pvalue)] <- pvalue
		}

		pvalue.list <- unlist(pvalue.list)
		return(pvalue.list)
	})
	surv.pvalue <- do.call(rbind, surv.pvalue)
	rownames(surv.pvalue) <- colnames(clonality.data)
	return(surv.pvalue)
}

CalSurvPvalue <- function(clonality.data, surv.data, grade, subtype){
	
	sur.os.pvalue <- ClonalitySurvPvalue(clonality.data=clonality.data, surv.data=surv.data, 
		surv.type='os', grade = grade, subtype = subtype)
	sur.dss.pvalue <- ClonalitySurvPvalue(clonality.data=clonality.data, surv.data=surv.data, 
		surv.type='dss', grade = grade, subtype = subtype)
	sur.pfi.pvalue <- ClonalitySurvPvalue(clonality.data=clonality.data, surv.data=surv.data, 
		surv.type='pfi', grade = grade, subtype = subtype)
	
	sur.pvalue <- cbind(sur.os.pvalue, sur.dss.pvalue, sur.pfi.pvalue)
	return(sur.pvalue)
}


wgistic.peak.het <- readRDS('/data/wgistic_het.rds')
tcga.cli.data <- readRDS('/data/tcga_glioma_cli_mol.rds')
gold.set <- readRDS('/data/gold_set.rds')


tcga.cli.data <- tcga.cli.data[gold.set, ]
rownames(tcga.cli.data) <- paste0(rownames(tcga.cli.data), '-01', sep='')
wgistic.peak.het <- as.data.frame(t(wgistic.peak.het[, paste(gold.set, '01', sep = '-')]))


setwd('/result/Section5')
codel.sur.pvalue <- CalSurvPvalue(wgistic.peak.het, tcga.cli.data, grade = c('G2', 'G3'), subtype = c('Oligodendroglioma,IDHmut-codel'))
noncodel.sur.pvalue <- CalSurvPvalue(wgistic.peak.het, tcga.cli.data, grade = c('G2', 'G3' ,'G4'), subtype = c('Astrocytoma,IDHmut'))
wt.sur.pvalue <- CalSurvPvalue(wgistic.peak.het, tcga.cli.data, grade = c('G4'), subtype = c('Glioblastoma,IDHwt'))

# write.table(codel.sur.pvalue, file='codel.sur.pvalue.txt', sep='\t', quote=FALSE, row.names=TRUE, col.names=TRUE)
# write.table(noncodel.sur.pvalue, file='noncodel.sur.pvalue.txt', sep='\t', quote=FALSE, row.names=TRUE, col.names=TRUE)
# write.table(wt.sur.pvalue, file='wt.sur.pvalue.txt', sep='\t', quote=FALSE, row.names=TRUE, col.names=TRUE)


all.sur.pvalue <- cbind(wt.sur.pvalue, noncodel.sur.pvalue, codel.sur.pvalue)
all.sur.pvalue <- all.sur.pvalue[apply(all.sur.pvalue, 1, function(x) any(x[seq(1, 81, 3)] < 0.05, na.rm=T)), ]
				

library('circlize')
library('ComplexHeatmap')
library('RColorBrewer')

anno.data <- data.frame(subtype = rep(x = c('Glioblastoma,IDHwt', 'Astrocytoma,IDHmut', 'Oligodendroglioma,IDHmut-codel'), each = 27), 
	surtype = rep.int(rep(x = c('OS', 'DSS', 'PFI'), each = 9), times = 3), 
	comtype = rep.int(rep(x = c('Clonal vs Wt', 'Subclonal vs Wt', 'Clonal vs Subclonal'), each = 3), times = 9), 
	pvaluetype = rep.int(x = c('lr.p', 'univ.p', 'mult.p'), times = 27))


col_fun = colorRamp2(c(0, 30), c("white", "red"))
mol.col <- c("#00AFBB", "#E7B800", "#FC4E07")
names(mol.col) <- c('Glioblastoma,IDHwt', 'Astrocytoma,IDHmut', 'Oligodendroglioma,IDHmut-codel')

sur.col <- brewer.pal(8, 'Set1')[c(2, 3, 4)] 
names(sur.col) <- c('OS', 'DSS', 'PFI')

com.col <- brewer.pal(8, 'Dark2')[1:3] 
names(com.col) <- c('Clonal vs Wt', 'Subclonal vs Wt', 'Clonal vs Subclonal')
   
pva.col <- brewer.pal(8, 'Dark2')[4:6]
names(pva.col) <- c('lr.p', 'univ.p', 'mult.p')


ha = HeatmapAnnotation(df = anno.data, col = list(subtype = mol.col, surtype = sur.col, comtype = com.col, pvaluetype = pva.col), 
show_legend = FALSE, show_annotation_name = FALSE)

pdf('result/Section5/pvalue_heatmap.pdf', width=14)
	Heatmap(matrix = -log2(all.sur.pvalue), cluster_columns = FALSE, cluster_rows = FALSE, column_split = rep(1:3, each = 27), 
	show_column_names = FALSE, row_names_side = 'left', row_names_gp = gpar(fontsize = 7), name = "p value", col = col_fun,
	cell_fun = function(j, i, x, y, width, height, fill) {
		grid.text(sprintf("%.3f", all.sur.pvalue[i, j]), x, y, gp = gpar(fontsize = 2.5))}, 
	column_title = NULL, top_annotation = ha, show_heatmap_legend = FALSE)
dev.off()

