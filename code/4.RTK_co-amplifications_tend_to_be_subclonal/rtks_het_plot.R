
############################################Receptor tyrosine kinases
setwd('/pub5/xiaoyun/Jobs/J22/CopyNumberClonalityProject/Results/Section3/Resources')
rtks <- read.csv(file='RTKs.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
rtks$HGNC.ID <- gsub(pattern='HGNC:', replacement='', rtks$HGNC.ID)


gold.set <- readRDS('/pub5/xiaoyun/Jobs/J22/CopyNumberClonalityProject/Resource/CuratedData/gold_set.rds')
gene.het <- readRDS(file='/pub5/xiaoyun/Jobs/J22/CopyNumberClonalityProject/Results/Section3/Resources/gene.het.5.0.rds')

wrtk.het <- gene.het[rtks$Approved.symbol, ]
wrtk.het <- wrtk.het[, paste(gold.set, '01', sep = '-')]
wrtk.het[wrtk.het < 0] <- 0
wrtk.het <- wrtk.het[, colSums(abs(wrtk.het))!=0]


# annotation
gli.cli.mol.data <- readRDS(file='/pub5/xiaoyun/Jobs/J22/CopyNumberClonalityProject/Resource/CuratedData/tcga_glioma_cli_mol.rds')
rownames(gli.cli.mol.data) <- paste(rownames(gli.cli.mol.data), '01', sep = '-')

# gli.cli.mol.data <- subset(gli.cli.mol.data, histological_grade %in% ('G4'))
com.samples <- intersect(colnames(wrtk.het), rownames(gli.cli.mol.data))
gli.cli.mol.data <- gli.cli.mol.data[com.samples, ]
wrtk.het <- wrtk.het[, com.samples]


# sort
wrtk.het <- wrtk.het[order(rowSums(wrtk.het != 0), decreasing=T), ]
flag <- paste(rep(letters, rep(10, length(letters))), 0:9, sep = "")	
order.index <- c()
for(i in 1:ncol(wrtk.het)){
 index <- which(wrtk.het[, i] != 0)
 order.index <- c(order.index, paste(flag[index], collapse = ""))
}
wrtk.het <- wrtk.het[, order(order.index)]



############ plot
library('ComplexHeatmap')
library('RColorBrewer')

mol.col <- c("#00AFBB", "#E7B800", "#FC4E07")
names(mol.col) <- c("IDHwt", "IDHmut-non-codel", "IDHmut-codel")

ha = HeatmapAnnotation(df = gli.cli.mol.data[, 'IDH_CODEL_SUBTYPE', FALSE], col = list(IDH_CODEL_SUBTYPE = mol.col), 
 show_legend = FALSE, show_annotation_name = FALSE)


pdf('/pub5/xiaoyun/Jobs/J22/CopyNumberClonalityProject/Results/Section3/Results/RTKClonalityPlot/all_rtk_plot.5.0.pdf')
 het.col <- c("#FC9272", "#CB181D", "#F7F7F7", "#084594", "#9ECAE1")
 names(het.col) <- c(2, 1, 0, -1, -2)
 Heatmap(matrix = wrtk.het, cluster_columns = FALSE, cluster_rows = FALSE, col = het.col, 
  column_split = seq(1, ncol(wrtk.het)), row_split = seq(1, nrow(wrtk.het)), 
  show_column_names = FALSE, row_names_side = 'left', top_annotation = ha, row_names_gp = gpar(fontsize = 8), 
  heatmap_legend_param = list(title = "Clonal status", at = c(2, 1, 0, -1, -2), column_gap=unit(0.5, "mm"), 
  labels = c("subclonal amplification", "clonal amplification", "wild type", 
  "clonal deletion", "subclonal deletion"), border = "black"), show_heatmap_legend = FALSE)
dev.off()
