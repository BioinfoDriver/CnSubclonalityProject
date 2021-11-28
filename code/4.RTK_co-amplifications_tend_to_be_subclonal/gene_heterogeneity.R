
############################################HGNC protein-coding gene
setwd('/data/OriginalData/')
gene.info <- read.csv(file='gene_with_protein_product.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
gene.info <- gene.info[, c('hgnc_id', 'symbol', 'name', 'location', 'alias_symbol', 
 'gene_family', 'gene_family_id', 'entrez_id', 'ensembl_gene_id')]
gene.info$hgnc_id  <- gsub(pattern='HGNC:', replacement='', gene.info$hgnc_id)

library(biomaRt)
ensembl <- useEnsembl(biomart = 'genes', 
 dataset = 'hsapiens_gene_ensembl', host = 'https://feb2014.archive.ensembl.org', version = 75)

gene.region <- getBM(attributes=c('hgnc_id', 'hgnc_symbol', 'chromosome_name', 'start_position', 
 'end_position', 'strand', 'band'), filters=c('hgnc_id', 'hgnc_symbol', 'chromosome_name'), 
 values=list(gene.info$hgnc_id, gene.info$symbol, c(1:22, 'X', 'Y')), mart=ensembl)
gene.region <- gene.region[!duplicated(gene.region[, c('hgnc_id', 'hgnc_symbol')]), ]
 
library(GenomicRanges)
gene.region$strand <- '*'
gene.region$chromosome_name <- paste0('chr', gene.region$chromosome_name)
colnames(gene.region) <- c('hgnc_id', 'hgnc_symbol', 'chromosome_name', 'start', 'end', 'strand', 'band')
gene.region <- makeGRangesFromDataFrame(df=gene.region, keep.extra.columns=TRUE, ignore.strand=FALSE, seqinfo=NULL,
 seqnames.field="chromosome_name", start.field="start", end.field="end", strand.field="strand", starts.in.df.are.0based=FALSE)


saveRDS(gene.region, file='/data/gene.region.rds')


############################################

GetMatchedAbsoluteCalls <- function(absGRL, max.cn=Inf)
{
    library('RaggedExperiment')
	if(max.cn < Inf)
        absGRL <- endoapply(absGRL, 
            function(x) x[.totalCN(x) <= max.cn])

    RaggedExperiment(absGRL)
}


MaxScore <- function(scores, ranges, qranges) max(scores, na.rm=TRUE)
Wmean <- function(scores, ranges, qranges)
{
    isects <- pintersect(ranges, qranges)
    s <- sum(scores * width(isects)) / sum(width(isects))
    return(round(s))
}

QueryGisticSubAlt <- function(ra, query, 
    sum.method=c("any", "wmean"), sum.label, ext.range=0)
{
    sum.method <- match.arg(sum.method)
    sum.method <- ifelse(sum.method == "wmean", Wmean, MaxScore)

    library('RaggedExperiment')
	if(ext.range)
    {
        start(query) <- sapply(start(query),
            function(s) ifelse(s < ext.range, 1, s-ext.range))
        end(query) <- sapply(end(query), function(s) s + ext.range)
    }

    qa <- qreduceAssay(ra, query, 
                simplifyReduce=sum.method, i=sum.label, background=NA)
    return(qa)
}

abs.seg.call <- readRDS(file='/data/tcga_gli_seg_call_grl.rds')
gene.region <- readRDS(file='/data/gene.region.rds')

ra <- GetMatchedAbsoluteCalls(abs.seg.call)
wround.gene.subcl <- QueryGisticSubAlt(ra, query=gene.region, sum.label='score', sum.method="wmean")
wround.gene.alt <- QueryGisticSubAlt(ra, query=gene.region, sum.label='Modal_Total_CN', sum.method="wmean")


QualitativeGisticAlt <- function(alt, ploi){
	
	isect.samples <- intersect(colnames(alt), rownames(ploi))
	ploi <- ploi[isect.samples, ]
	alt <- alt[, isect.samples]
	
	alt.mat <- sapply(colnames(alt), function(samp){
		alt.vec <- alt[, samp]
		wildtype.index.1 <- which(is.na(alt.vec))
		
		amp.index <- which(alt.vec >= round(ploi[samp, 'ploidy']) + 2)
		del.index <- which(alt.vec <= round(ploi[samp, 'ploidy']) - 2)
		wildtype.index.2 <- which((alt.vec > round(ploi[samp, 'ploidy']) - 2) & (alt.vec < round(ploi[samp, 'ploidy']) + 2))
	
		alt.vec[wildtype.index.1] <- 0
		alt.vec[wildtype.index.2] <- 0
		alt.vec[amp.index] <- 1
		alt.vec[del.index] <- -1
		return(alt.vec)
	})
	
	return(alt.mat)
}


abs.puri.ploi <- readRDS(file='/data/tcga_gli_puri_ploi.rds')
wround.gene.alt <- QualitativeGisticAlt(wround.gene.alt, abs.puri.ploi)

rownames(wround.gene.alt) <- rownames(wround.gene.subcl) <- gene.region$hgnc_symbol


GisticHetPeak <- function(subcl, alt)
{

	samples <- intersect(colnames(subcl), colnames(alt))
	subcl <- subcl[, samples]
	alt <- alt[, samples]	

	# wild type -> 0; clonal gain -> 1; subclonal gain -> 2; 
	# wild type -> 0; clonal loss -> -1; subclonal loss -> -2; 

	het.matrix <- matrix(NA, nrow = nrow(alt), ncol = ncol(alt))
	rownames(het.matrix) <- rownames(alt)
	colnames(het.matrix) <- colnames(alt)
	
	
	het.matrix[alt == 0] <- 0
	
	het.matrix[alt == 1 & subcl == 0] <- 1
	het.matrix[alt == 1 & subcl == 1] <- 2
	
	het.matrix[alt == -1 & subcl == 0] <- -1
	het.matrix[alt == -1 & subcl == 1] <- -2
	
	return(het.matrix)
}

gene.het <- GisticHetPeak(wround.gene.subcl, wround.gene.alt)

saveRDS(gene.het, file='/data/gene.het.rds')




