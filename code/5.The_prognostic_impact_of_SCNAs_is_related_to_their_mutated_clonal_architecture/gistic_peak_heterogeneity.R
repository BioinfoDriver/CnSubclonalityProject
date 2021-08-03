

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

abs.seg.call <- readRDS(file='/pub5/xiaoyun/Jobs/J22/CopyNumberClonalityProject/Resource/CuratedData/tcga_gli_seg_call_grl.rds')
gistic.peak <- readRDS(file='/pub5/xiaoyun/Jobs/J22/CopyNumberClonalityProject/Resource/gbmlgg_gistic_peak.rds')


ra <- GetMatchedAbsoluteCalls(abs.seg.call)
wround.gistic.subcl <- QueryGisticSubAlt(ra, query=rowRanges(gistic.peak), sum.label='score', sum.method="wmean")
mround.gistic.subcl <- QueryGisticSubAlt(ra, query=rowRanges(gistic.peak), sum.label='score', sum.method="any")

wround.gistic.alt <- QueryGisticSubAlt(ra, query=rowRanges(gistic.peak), sum.label='Modal_Total_CN', sum.method="wmean")
mround.gistic.alt <- QueryGisticSubAlt(ra, query=rowRanges(gistic.peak), sum.label='Modal_Total_CN', sum.method="any")



QualitativeGisticAlt <- function(alt, ploi){
	
	isect.samples <- intersect(colnames(alt), rownames(ploi))
	ploi <- ploi[isect.samples, ]
	alt <- alt[, isect.samples]
	
	alt.mat <- sapply(colnames(alt), function(samp){
		alt.vec <- alt[, samp]
		wildtype.index.1 <- which(is.na(alt.vec))
		

		amp.index <- which(alt.vec > round(ploi[samp, 'ploidy']))
		del.index <- which(alt.vec < round(ploi[samp, 'ploidy']))
		wildtype.index.2 <- which(alt.vec == round(ploi[samp, 'ploidy']))

		
		alt.vec[wildtype.index.1] <- 0
		alt.vec[wildtype.index.2] <- 0
		alt.vec[amp.index] <- 1
		alt.vec[del.index] <- -1
		return(alt.vec)
	})
	
	return(alt.mat)
}


abs.puri.ploi <- readRDS(file='/pub5/xiaoyun/Jobs/J22/CopyNumberClonalityProject/Resource/CuratedData/tcga_gli_puri_ploi.rds')

wround.gistic.alt <- QualitativeGisticAlt(wround.gistic.alt, abs.puri.ploi)
mround.gistic.alt <- QualitativeGisticAlt(mround.gistic.alt, abs.puri.ploi)



rownames(wround.gistic.alt) <- rownames(mround.gistic.alt) <- 
 paste0(substr(mcols(gistic.peak)$type, 1, 3), '_', mcols(gistic.peak)$band)


# retain dominant type
RetainDomType <- function(alt.mat){
	dom.alt.mat <- lapply(rownames(alt.mat), function(band){
		alt.vec <- alt.mat[band, ]
		type <- substr(band, 1, 3)
		
		if(type == 'Amp'){
			alt.vec[alt.vec == -1] <- 0
		
		}else{
			alt.vec[alt.vec == 1] <- 0
		
		}
	
		return(alt.vec)
	})
	dom.alt.mat <- do.call(rbind, dom.alt.mat)
	rownames(dom.alt.mat) <- rownames(alt.mat)
	
	return(dom.alt.mat)
}

wround.gistic.alt <- RetainDomType(wround.gistic.alt)
mround.gistic.alt <- RetainDomType(mround.gistic.alt)



GisticHetPeak <- function(gistic.subcl, gistic.alt)
{

	samples <- intersect(colnames(gistic.subcl), colnames(gistic.alt))
	gistic.subcl <- gistic.subcl[, samples]
	gistic.alt <- gistic.alt[, samples]	

	# wild type -> 0; clonal gain -> 1; subclonal gain -> 2; 
	# wild type -> 0; clonal loss -> -1; subclonal loss -> -2; 

	het.matrix <- matrix(NA, nrow = nrow(gistic.alt), ncol = ncol(gistic.alt))
	rownames(het.matrix) <- rownames(gistic.alt)
	colnames(het.matrix) <- colnames(gistic.alt)
	
	
	het.matrix[gistic.alt == 0] <- 0
	
	het.matrix[gistic.alt == 1 & gistic.subcl == 0] <- 1
	het.matrix[gistic.alt == 1 & gistic.subcl == 1] <- 2
	
	het.matrix[gistic.alt == -1 & gistic.subcl == 0] <- -1
	het.matrix[gistic.alt == -1 & gistic.subcl == 1] <- -2
	
	return(het.matrix)
}

wgistic.het <- GisticHetPeak(wround.gistic.subcl, wround.gistic.alt)
mgistic.het <- GisticHetPeak(mround.gistic.subcl, mround.gistic.alt)


saveRDS(wgistic.het, file='/pub5/xiaoyun/Jobs/J22/CopyNumberClonalityProject/Resource/CuratedData/wgistic_het.rds')
saveRDS(mgistic.het, file='/pub5/xiaoyun/Jobs/J22/CopyNumberClonalityProject/Resource/CuratedData/mgistic_het.rds')

saveRDS(wround.gistic.alt, file='/pub5/xiaoyun/Jobs/J22/CopyNumberClonalityProject/Resource/CuratedData/wgistic_alt.rds')
saveRDS(mround.gistic.alt, file='/pub5/xiaoyun/Jobs/J22/CopyNumberClonalityProject/Resource/CuratedData/mgistic_alt.rds')



