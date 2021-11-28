# infile = TCGA_mastercalls.abs_tables_JSedit.fixed.txt
# @returns: a matrix with sample IDs as rownames and columns 'purity', 'ploidy',
#           'Genome.doublings',  and 'Subclonal.genome.fraction'  
ReadPurityPloidy <- function(infile)
{
    dat <- read.delim(infile, as.is=TRUE)
    dat <- dat[, -c(2, 7, 8, 10)]
    rownames(dat) <- dat[, 1]
    dat <- dat[, -1]
	
	dat <- dat[!is.na(dat[, 1]), ]
    return(dat)
}

in.file <- '/data/OriginalData/tcga_glioma_purity_ploidy.txt'
gli.puri.ploi <- ReadPurityPloidy(infile=in.file)
saveRDS(gli.puri.ploi, file='/data/tcga_gli_puri_ploi.rds')


### processing raw data
# @args 
#  infile = TCGA_mastercalls.abs_segtabs.fixed.txt
#  ploi: the sample’s ploidy (rounded to the nearest integer)
#
# @returns: a GRangesList with per-sample ABSOLUTE calls
ReadAbsolute <- function(infile, ploi)
{
	library('GenomicRanges')
	
    df <-  read.delim(infile, as.is=TRUE)
    df <-  df[!is.na(df[, "Chromosome"]), ]
    df <-  df[df[, "Chromosome"] != 23, ]
    df[, "Chromosome"] <- paste0("chr", df[, "Chromosome"])

	isect.samples <- intersect(df$Sample, rownames(ploi))
	ploi <- ploi[isect.samples, ]
	df <- subset(df, Sample %in% isect.samples)
	
	df$ploidy <- ploi[match(df$Sample, rownames(ploi)), "ploidy"]
	
	df <- df[df$Modal_Total_CN != round(df$ploidy), ]
	
	IsSubClonal <- function(x) as.integer(x$Subclonal_HSCN_a1 | x$Subclonal_HSCN_a2)
	df$score <- IsSubClonal(df)

    rel.cols <- c("Sample", "Chromosome", "Start", "End", "Modal_Total_CN", "ploidy", "score")
    df <- df[, rel.cols]
    grl <- makeGRangesListFromDataFrame(df, 
        split.field="Sample", keep.extra.columns=TRUE) 
    
	return(grl)	
}

gli.puri.ploi <- readRDS('/data/tcga_gli_puri_ploi.rds')

in.file <- '/data/OriginalData/tcga_glioma_abs_seg.txt'

gli.seg.call.grl <- ReadAbsolute(infile=in.file, ploi=gli.puri.ploi)
saveRDS(gli.seg.call.grl, file='/data/tcga_gli_seg_call_grl.rds')
