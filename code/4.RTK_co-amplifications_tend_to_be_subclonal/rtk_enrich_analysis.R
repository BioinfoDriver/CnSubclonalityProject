
# load data
diff.exp <- readRDS(file='/data/rtk_idhwt_gbm_diff_exp.rds')


# https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/
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


setwd('result/Section4/')
write.table(gsea.go.enrich[1:50, c(1, 2, 5, 6, 7)], file='gsea.go.enrich.txt', sep='\t', 
	col.names=TRUE, row.names=FALSE, quote=FALSE)

require(DOSE)
pdf(file='gsea_go_enrich.pdf')
dotplot(gsea.go.enrich, showCategory=30, split=".sign", font.size=6) + facet_grid(.~.sign)
dev.off()

