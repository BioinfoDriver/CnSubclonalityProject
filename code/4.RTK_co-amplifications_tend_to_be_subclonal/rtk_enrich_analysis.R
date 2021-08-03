
setwd('/pub5/xiaoyun/Jobs/J22/CopyNumberClonalityProject/Results/Section3/Results/RTKDifferentExpress')
diff.exp <- readRDS(file='rtk_idhwt_gbm_diff_exp.rds')


# preranked gene list for performing GSEA
PrerankedGSEAR <- function(diff.exp, file.name){
 
 diff.exp <- diff.exp[!is.na(diff.exp$pvalue), ]
 diff.exp$ranks <- sign(diff.exp$log2FoldChange) * -log10(diff.exp$pvalue)

 diff.ranks <- diff.exp[, c('symbol', 'ranks')]
 diff.ranks <- diff.ranks[!is.na(diff.ranks$symbol), ]
 
 write.table(diff.ranks, file=file.name, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
}

PrerankedGSEAR(diff.exp, "idhwt_gbm_diff_gene_list.rnk")


# https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/
# GO Enrichment Analysis of a gene set
GOEnrichmentAnalysis <- function(diff.exp, fc.cutoff, fdr.cutoff){
 library(org.Hs.eg.db)
 library(clusterProfiler)
 background <- rownames(diff.exp)
 diff.exp <- subset(diff.exp, abs(log2FoldChange) >= fc.cutoff & padj < fdr.cutoff)
 gene.list <- rownames(diff.exp)
 
 enrich <- enrichGO(gene=gene.list, OrgDb=org.Hs.eg.db, keyType="ENSEMBL", ont="BP", minGSSize=10, maxGSSize=500,
  pvalueCutoff=0.05, pAdjustMethod="fdr", qvalueCutoff=0.1, readable=TRUE, universe=background)

 enrich <- as.data.frame(enrich)
 return(enrich)
}

go.enrich <- GOEnrichmentAnalysis(diff.exp, fc.cutoff=1, fdr.cutoff=0.05)


# Gene Set Enrichment Analysis of Gene Ontology
GSEAofGO <- function(diff.exp){
 library(clusterProfiler)
 library(org.Hs.eg.db)
 gene.list <- diff.exp$log2FoldChange
 names(gene.list) <- rownames(diff.exp)
 gene.list = sort(gene.list, decreasing = TRUE)

 gse.res <- gseGO(geneList=gene.list, ont ="BP", OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", 
	minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05, verbose = TRUE, pAdjustMethod = "BH")
			 
 return(gse.res)
}

gsea.go.enrich <- GSEAofGO(diff.exp)



# Gene Set Enrichment Analysis of KEGG

GSEAofKEGG <- function(diff.exp){
 library(clusterProfiler)
 library(org.Hs.eg.db)
 
 gene.list <- diff.exp$log2FoldChange
 names(gene.list) <- rownames(diff.exp)
 gene.list = sort(gene.list, decreasing = TRUE)

 
 ENTREZ.ID <- mapIds(org.Hs.eg.db, key=names(gene.list), column="ENTREZID", keytype="ENSEMBL", multiVals="first")
 ENTREZ.ID <- na.omit(ENTREZ.ID)
 
 gene.list <- gene.list[names(ENTREZ.ID)]
 gene.list = sort(gene.list, decreasing = TRUE)

 gse.res <- gseKEGG(geneList = gene.list, organism = 'hsa', keyType = "ncbi-geneid", 
	minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05, pAdjustMethod = "BH")
 
 return(gse.res)
}

gsea.kegg.enrich <- GSEAofKEGG(diff.exp)


setwd('/pub5/xiaoyun/Jobs/J22/CopyNumberClonalityProject/Results/Section3/Results/RTKEnrichment/')

write.table(gsea.go.enrich[, c(1, 2, 5, 6, 7)], file='gsea.go.enrich.txt', sep='\t', 
	col.names=TRUE, row.names=FALSE, quote=FALSE)


require(DOSE)
pdf(file='gsea_go_enrich.pdf')
dotplot(gsea.go.enrich, showCategory=30, split=".sign", font.size=6) + facet_grid(.~.sign)
dev.off()



show.term <- c('GO:0042391', 'GO:0099003', 'GO:0044786', 'GO:0099504', 'GO:0071156', 'GO:0006260', 'GO:0006836', 
 'GO:0050803', 'GO:0050804', 'GO:0000819', 'GO:0007093', 'GO:0044843', 'GO:0050807', 'GO:0022604')



mitotic spindle organization,
regulation of synaptic vesicle exocytosis
cell cycle checkpoint signaling
regulation of neuronal synaptic plasticity
dendrite morphogenesis
exocytic process
mitotic cell cycle checkpoint
synaptic transmission, GABAergic
synaptic vesicle recycling
regulation of neurotransmitter receptor activity
sister chromatid segregation
synapse organization
regulation of inflammatory response
regulation of neurotransmitter levels
regulation of postsynaptic membrane neurotransmitter receptor levels
regulation of postsynapse organization
regulation of cell morphogenesis
regulation of synapse organization
neutral lipid catabolic process
modulation of chemical synaptic transmission
regulation of synapse structure or activity
neurotransmitter transport
potassium ion transmembrane transport
DNA replication
cell cycle G1/S phase transition
regulation of cell cycle arrest
synaptic vesicle cycle
cell cycle DNA replication
vesicle-mediated transport in synapse
regulation of membrane potential
regulation of ion transmembrane transport
