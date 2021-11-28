

GetMatchedAbsoluteCalls <- function(absGRL, max.cn=Inf)
{
    library('RaggedExperiment')
	if(max.cn < Inf)
        absGRL <- endoapply(absGRL, 
            function(x) x[x$Modal_Total_CN <= max.cn])

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

abs.seg.call <- readRDS(file='/data/glass_seg_call_grl.rds')
gene.region <- readRDS(file='/data/gene.region.rds')

ra <- GetMatchedAbsoluteCalls(abs.seg.call, max.cn=14)
wround.gene.subcl <- QueryGisticSubAlt(ra, query=gene.region, sum.label='score', sum.method="wmean")
wround.gene.alt <- QueryGisticSubAlt(ra, query=gene.region, sum.label='Modal_Total_CN', sum.method="wmean")


QualitativeGisticAlt <- function(alt, ploi){
	
	isect.samples <- intersect(colnames(alt), rownames(ploi))
	ploi <- ploi[isect.samples, ]
	alt <- alt[, isect.samples]
	
	alt.mat <- sapply(colnames(alt), function(samp){
		alt.vec <- alt[, samp]
		wildtype.index.1 <- which(is.na(alt.vec))
		
		amp.index <- which(alt.vec >= 2*round(ploi[samp, 'ploidy']) + 1)
		del.index <- which(alt.vec <= round(ploi[samp, 'ploidy']) - 1)
		wildtype.index.2 <- which((alt.vec > round(ploi[samp, 'ploidy']) - 1) & (alt.vec < 2*round(ploi[samp, 'ploidy']) + 1))
	
		alt.vec[wildtype.index.1] <- 0
		alt.vec[wildtype.index.2] <- 0
		alt.vec[amp.index] <- 1
		alt.vec[del.index] <- -1
		return(alt.vec)
	})
	
	return(alt.mat)
}

setwd('/data/OriginalData/syn17038081')
abs.puri.ploi <- read.csv(file='variants_seqz_params.tsv', header=TRUE, sep='\t', stringsAsFactors=FALSE)
rownames(abs.puri.ploi) <- paste0(substr(abs.puri.ploi$pair_barcode, 1, 18), substr(abs.puri.ploi$pair_barcode, 22, 29))


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

saveRDS(gene.het, file='/data/glass_gene_het.rds')
