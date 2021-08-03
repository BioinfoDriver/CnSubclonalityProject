
## GISTIC
# outdir: a directory to store temporary files
# @returns: a RangedSummarizedExperiment
Gistic2RSE <- function(ctype=c("LGG", "GBM", "GBMLGG"), peak=c("wide", "narrow", "full"), outdir){
    ctype <- match.arg(ctype)
    peak <- match.arg(peak)
    rname <- paste("gistic", ctype, peak, sep="_")
    
    library('SummarizedExperiment')
    # download the tar
    BROAD.URL <- paste0("http://gdac.broadinstitute.org/",
                "runs/analyses__2016_01_28/data/ctype/20160128")
 
    GISTIC.FILE <- paste0("gdac.broadinstitute.org_ctype-TP.",
                    "CopyNumber_Gistic2.Level_4.2016012800.0.0.tar.gz")

    url <- file.path(BROAD.URL, GISTIC.FILE)
    url <- gsub("ctype", ctype, url)
    if(ctype == "LAML") url <- sub("TP", "TB", url)
    else if(ctype == "SKCM") url <- sub("TP", "TM", url)
    
    dest.file <- paste0(ctype, "_gistic2.tar.gz")
    dest.file <- file.path(outdir, dest.file)
    download.file(url, dest.file)
    files <- untar(dest.file, list=TRUE)
    
    # extract the lesions file
    basef <- "all_lesions.conf_99.txt"
    basef <- files[grepl(basef, files)]
    untar(dest.file, files=basef)

    # read the lesion file
    gistic <- read.delim(basef, as.is=TRUE)
    file.remove(dest.file)    
    unlink(dirname(basef), recursive=TRUE, force=TRUE)

    # transform to RangedSummarizedExperiment
    rel.rows <- grepl("Peak +[0-9]+$", gistic[,1])
    gistic <- gistic[rel.rows,]

    # (a) get the ranges from chosen peaks
    peak.col <- grep(peak, c("wide", "narrow", "full")) + 2
    ranges <- gistic[rel.rows,peak.col]
    ranges <- sub("\\(probes [0-9]+:[0-9]+\\) *$", "", ranges)   
    ranges <- as(ranges, "GRanges")
    ind <- seqnames(ranges) != "chrX"
    ranges <- ranges[ind]
    gistic <- gistic[as.logical(ind),]
    ind <- orderSeqlevels(seqlevels(ranges))
    seqlevels(ranges) <- seqlevels(ranges)[ind]
    ind <- order(ranges)
    ranges <- ranges[ind]
    gistic <- gistic[ind, ]
  
    # (b) get the peak type (amplification / deletion) 
    peak.type <- sapply(gistic[,1], 
        function(x) unlist(strsplit(x, " "))[1], USE.NAMES=FALSE)
    
	# () get descriptors of peak
	bands <- gistic[,2]
	bands <- trimws(bands, which = "right")
	
    # (d) create the SE
    rel.cols <- grepl("^TCGA", colnames(gistic))
    gistic <- gistic[rel.cols]
    gistic <- as.matrix(gistic)
    gisticSE <- SummarizedExperiment(assays=list(counts=gistic))
    rowRanges(gisticSE) <- ranges
    rowData(gisticSE)$type <- peak.type 
	rowData(gisticSE)$band <- bands 
	
    # ensure consistent naming with subtypes and absolute
    colnames(gisticSE)<- gsub("\\.", "-", colnames(gisticSE))
    colnames(gisticSE) <- TCGAutils::TCGAbarcode(colnames(gisticSE), sample=TRUE)
    colnames(gisticSE)<- sub("A$", "", colnames(gisticSE))

    return(gisticSE)
}

file.path <- '/data/OriginalData/'
gbmlgg.gistic.peak <- Gistic2RSE(ctype='GBMLGG', peak = 'wide', outdir = file.path)
saveRDS(gbmlgg.gistic.peak, file = '/data/OriginalData/gbmlgg_gistic_peak.rds')
