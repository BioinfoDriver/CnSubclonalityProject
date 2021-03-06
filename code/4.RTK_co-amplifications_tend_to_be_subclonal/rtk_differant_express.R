
# RNA-SEQ raw counts
setwd('/data/OriginalData/counts')
sample.info <- read.csv(file = 'gdc_sample_sheet.2020-07-21.tsv', sep = '\t', header = TRUE, stringsAsFactors = FALSE)

sample.info <- subset(sample.info, Sample.Type == 'Primary Tumor')
sample.info <- sample.info[!duplicated(sample.info$Case.ID), ]
sample.info$File.Name <- paste(sample.info$File.ID, sample.info$File.Name, sep = "/")

sample.info <- sample.info[, c('Sample.ID', 'File.Name', 'Project.ID')]

library("DESeq2")
count <- DESeqDataSetFromHTSeqCount(sampleTable=sample.info, directory=".", design= ~ 1)
saveRDS(count, file='/data/pan_gli_exp_count.rds')


# procoding genes
gene.info <- read.csv(file='/data/OriginalData/gene_with_protein_product.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
gene.info <- gene.info[, c('symbol', 'ensembl_gene_id')]


######## clonal rtk alt Vs subclonal rtk alt
rtks <- read.csv(file='/data/RTKs.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
tcga.cli.data <- readRDS('/data/tcga_glioma_cli_mol.rds')

gold.set <- readRDS('/data/gold_set.rds')
gene.het <- readRDS(file='/data/gene.het.rds')
gene.het <- gene.het[, paste(gold.set, '01', sep = '-')]


RtkCliMolAnalysis <- function(het.mat, rtk, group.index, cli.data, grade, subtype){
 
 rtk.het <- het.mat[rtk, ]

 rtk.clo.sam <- colnames(rtk.het[, colSums(rtk.het == 1) >= 1])
 rtk.sub.sam <- colnames(rtk.het[, colSums(rtk.het == 2) >= 1])
 
 if(group.index=='subclonaldominant'){
  rtk.clo.sam <- setdiff(rtk.clo.sam, rtk.sub.sam)

 }else if(group.index=='clonaldominant'){
  rtk.sub.sam <- setdiff(rtk.sub.sam, rtk.clo.sam)
 
 }else{
  com.sam <- intersect(rtk.clo.sam, rtk.sub.sam)
  rtk.clo.sam <- setdiff(rtk.clo.sam, com.sam)
  rtk.sub.sam <- setdiff(rtk.sub.sam, com.sam)
  
 }
 
 rtk.clo.mat <- data.frame(group = rep(c('rtk.clonal.amp', 'rtk.subclonal.amp'), 
  times=c(length(rtk.clo.sam), length(rtk.sub.sam))))
 rownames(rtk.clo.mat) <- c(rtk.clo.sam, rtk.sub.sam)


 cli.data <- subset(cli.data, histological_grade %in% grade & Integrated_Diagnoses %in% subtype)
 rtk.clo.mat <- rtk.clo.mat[intersect(rownames(rtk.clo.mat), paste0(rownames(cli.data), '-01')), , FALSE]

 return(rtk.clo.mat)
}


# Differential expression analysis
exp.data <- counts(count)
rownames(exp.data) <- substr(rownames(exp.data), 1, 15)
colnames(exp.data) <- substr(colnames(exp.data), 1, 15)

keep <- rowSums(exp.data > 1) > ncol(exp.data)*0.5
exp.data <- exp.data[keep, ]

rtk.alt.sam <- RtkCliMolAnalysis(gene.het, rtks$Approved.symbol, 'exclude', tcga.cli.data, 'G4', 'Glioblastoma,IDHwt')
exp.data <- exp.data[intersect(rownames(exp.data), gene.info$ensembl_gene_id), 
	intersect(colnames(exp.data), rownames(rtk.alt.sam))]

rtk.alt.sam <- rtk.alt.sam[colnames(exp.data), , FALSE]
colnames(rtk.alt.sam) <- 'phenotype'
rtk.alt.sam$phenotype <- factor(rtk.alt.sam$phenotype, levels=c('rtk.clonal.amp', 'rtk.subclonal.amp'))
rtk.alt.sam$batch <- subset(tcga.cli.data, bcr_patient_barcode %in% substr(rownames(rtk.alt.sam), 1, 12))$cancer_type


RunDESeq2 <- function(countMatrix, pData){

	library(DESeq2);
	dds <- DESeqDataSetFromMatrix(countData=countMatrix, colData=pData, design = ~ batch + phenotype);
	
	dds <- DESeq(dds, parallel=T);
	dds <- replaceOutliersWithTrimmedMean(dds);

	res <- results(dds, alpha = 0.1, cooksCutoff=FALSE);
	res <- res[order(res$padj), ];
	
	return(res);
}

diff.exp <- RunDESeq2(countMatrix = exp.data, pData = rtk.alt.sam)
diff.exp$symbol <- gene.info$symbol[match(rownames(diff.exp), gene.info$ensembl_gene_id)]

saveRDS(diff.exp, file='/data/rtk_idhwt_gbm_diff_exp.rds')

setwd('/result/Section4/')
idhwt.gbm.diff.exp <- subset(diff.exp, abs(log2FoldChange) > 1 & padj < 0.05)
write.table(idhwt.gbm.diff.exp, "idhwt_gbm_diff_exp_gene.txt", sep="\t", quote=F, col.names=TRUE, row.names=TRUE)


