
setwd('/data/OriginalData/syn17038081')

glass.seqz.seg <- read.csv(file='variants_seqz_seg.tsv', header=TRUE, sep='\t', stringsAsFactors=FALSE)
glass.seqz.ploidy <- read.csv(file='variants_seqz_params.tsv', header=TRUE, sep='\t', stringsAsFactors=FALSE)

glass.seqz.seg <- merge(glass.seqz.seg, glass.seqz.ploidy[, c('pair_barcode', 'ploidy')], by='pair_barcode')

glass.seqz.seg$cna_status <- ifelse(glass.seqz.seg$copy_number == round(glass.seqz.seg$ploidy), 'neutral', 
 ifelse(glass.seqz.seg$copy_number > round(glass.seqz.seg$ploidy), 'amplified', 'deleted'))

glass.seqz.seg <- subset(glass.seqz.seg, cna_status != 'neutral')

glass.seqz.seg$pair_barcode <- paste0(substr(glass.seqz.seg$pair_barcode, 1, 18), substr(glass.seqz.seg$pair_barcode, 22, 29))


gold.set <- read.csv(file='analysis_gold_set.tsv', header=TRUE, sep='\t', stringsAsFactors=FALSE)
library(dplyr)
gold.set <- gold.set %>% mutate(ROW_ID = NULL, ROW_VERSION = NULL)
gold.set$tumor_barcode_a <- paste0(substr(gold.set$tumor_pair_barcode, 1, 18), substr(gold.set$tumor_pair_barcode, 22, 29))

gold.set$tumor_barcode_b <- paste0(substr(gold.set$tumor_pair_barcode, 1, 12), 
 substr(gold.set$tumor_pair_barcode, 19, 21), substr(gold.set$tumor_pair_barcode, 16, 18), 
 substr(gold.set$tumor_pair_barcode, 22, 29))

# curated
{
# > gold.set$tumor_barcode_a[!(gold.set$tumor_barcode_a %in% glass.seqz.seg$pair_barcode)]
 # [1] "GLSS-19-0266-TP-01-01D-WXS" "TCGA-06-0152-TP-02-02D-WGS"
 # [3] "TCGA-06-0171-TP-02-02D-WGS" "TCGA-19-0957-TP-02-02D-WXS"
 # [5] "TCGA-DH-A669-TP-12-12D-WGS" "TCGA-DU-5870-TP-11-11D-WGS"
 # [7] "TCGA-DU-6397-TP-11-11D-WGS" "TCGA-DU-6404-TP-11-11D-WGS"
 # [9] "TCGA-DU-6407-TP-13-13D-WGS" "TCGA-DU-7304-TP-12-12D-WGS"
# [11] "TCGA-FG-5963-TP-11-11D-WXS" "TCGA-FG-5965-TP-11-11D-WGS"
# [13] "TCGA-FG-A4MT-TP-11-11D-WGS" "TCGA-TM-A7CF-TP-11-11D-WGS"
# [15] "TCGA-TQ-A7RK-TP-11-11D-WGS" "TCGA-TQ-A7RV-TP-21-21D-WGS"
# [17] "TCGA-TQ-A8XE-TP-11-11D-WGS"
# > 
# > gold.set$tumor_barcode_b[!(gold.set$tumor_barcode_b %in% glass.seqz.seg$pair_barcode)]
 # [1] "TCGA-06-0125-R1-01-11D-WGS" "TCGA-06-0152-R1-02-01D-WGS"
 # [3] "TCGA-06-0171-R1-02-11D-WGS" "TCGA-06-0211-R1-01-02D-WGS"
 # [5] "TCGA-06-0221-R1-01-11D-WGS" "TCGA-19-0957-R1-02-11D-WXS"
 # [7] "TCGA-DH-A669-R1-12-11D-WGS" "TCGA-DU-5870-R1-11-12D-WGS"
 # [9] "TCGA-DU-6397-R1-11-12D-WGS" "TCGA-DU-6404-R2-11-11D-WGS"
# [11] "TCGA-DU-6407-R2-13-11D-WGS" "TCGA-DU-7304-R1-12-12D-WGS"
# [13] "TCGA-FG-5963-R1-11-12D-WXS" "TCGA-FG-5965-R2-11-11D-WGS"
# [15] "TCGA-FG-A4MT-R1-11-11D-WGS" "TCGA-TM-A7CF-R1-11-11D-WGS"
# [17] "TCGA-TQ-A7RK-R2-11-11D-WGS" "TCGA-TQ-A7RV-R1-21-11D-WGS"
# [19] "TCGA-TQ-A8XE-R1-11-11D-WGS"

extra.patients <- c(gold.set$tumor_barcode_a[!(gold.set$tumor_barcode_a %in% glass.seqz.seg$pair_barcode)], 
 gold.set$tumor_barcode_b[!(gold.set$tumor_barcode_b %in% glass.seqz.seg$pair_barcode)])
extra.patients <- unique(substr(extra.patients, 1, 12))

seqz.seg.sample <- unique(glass.seqz.seg$pair_barcode)
# seqz.seg.sample[substr(seqz.seg.sample, 1, 12) %in% extra.patients]

gold.set.curated <- data.frame(tumor_pair_barcode=NA, case_barcode=NA,
tumor_barcode_a=c("GLSS-19-0266-TP-01-02D-WXS","TCGA-06-0125-TP-01-01D-WGS","TCGA-06-0152-TP-02-01D-WGS",
"TCGA-06-0171-TP-02-01D-WGS","TCGA-06-0211-TP-01-01D-WGS","TCGA-06-0221-TP-01-01D-WGS",
"TCGA-19-0957-TP-01-01D-WXS","TCGA-DH-A669-TP-12-01D-WGS","TCGA-DU-5870-TP-11-01D-WGS",
"TCGA-DU-6397-TP-11-01D-WGS","TCGA-DU-6404-TP-11-01D-WGS","TCGA-DU-6407-TP-13-01D-WGS",
"TCGA-DU-7304-TP-12-01D-WGS","TCGA-FG-5963-TP-11-01D-WXS","TCGA-FG-5965-TP-11-01D-WGS",
"TCGA-FG-A4MT-TP-11-01D-WGS","TCGA-TM-A7CF-TP-11-01D-WGS","TCGA-TQ-A7RK-TP-11-01D-WGS",
"TCGA-TQ-A7RV-TP-21-01D-WGS","TCGA-TQ-A8XE-TP-11-01D-WGS"), 
tumor_barcode_b=c("GLSS-19-0266-R1-01-01D-WXS","TCGA-06-0125-R1-11-01D-WGS","TCGA-06-0152-R1-01-01D-WGS",
"TCGA-06-0171-R1-11-01D-WGS","TCGA-06-0211-R1-02-01D-WGS","TCGA-06-0221-R1-11-01D-WGS",
"TCGA-19-0957-R1-11-01D-WXS","TCGA-DH-A669-R1-11-01D-WGS","TCGA-DU-5870-R1-12-01D-WGS",
"TCGA-DU-6397-R1-12-01D-WGS","TCGA-DU-6404-R1-21-01D-WGS","TCGA-DU-6407-R1-12-01D-WGS",
"TCGA-DU-7304-R1-12-01D-WGS","TCGA-FG-5963-R1-12-01D-WXS","TCGA-FG-5965-R1-11-01D-WGS",
"TCGA-FG-A4MT-R1-11-01D-WGS","TCGA-TM-A7CF-R1-11-01D-WGS","TCGA-TQ-A7RK-R1-11-01D-WGS",
"TCGA-TQ-A7RV-R1-11-01D-WGS","TCGA-TQ-A8XE-R1-11-01D-WGS"))

gold.set.curated$case_barcode <- substr(gold.set.curated$tumor_barcode_a, 1, 12)

}

gold.set <- rbind(gold.set, gold.set.curated)
gold.set <- subset(gold.set, tumor_barcode_a %in% glass.seqz.seg$pair_barcode & 
 tumor_barcode_b %in% glass.seqz.seg$pair_barcode)
 

library(regioneR)
glass.seqz.seg$chrom <- paste0('chr', glass.seqz.seg$chrom)
glass.seqz.seg$chrom[glass.seqz.seg$chrom == 'chr23'] <- 'chrX'

cna.subclonality <- lapply(seq(nrow(gold.set)), function(index){
 
 sam.ids <- gold.set[index, ]
 regionAAmp <- subset(glass.seqz.seg, pair_barcode == sam.ids[1, 3] & cna_status == 'amplified')[, 
	 c('chrom', 'start', 'end', 'pair_barcode', 'copy_number', 'ploidy')]
 regionADel <- subset(glass.seqz.seg, pair_barcode == sam.ids[1, 3] & cna_status == 'deleted')[, 
	 c('chrom', 'start', 'end', 'pair_barcode', 'copy_number', 'ploidy')]

 
 regionBAmp <- subset(glass.seqz.seg, pair_barcode == sam.ids[1, 4] & cna_status == 'amplified')[, 
	 c('chrom', 'start', 'end', 'pair_barcode', 'copy_number', 'ploidy')]
 regionBDel <- subset(glass.seqz.seg, pair_barcode == sam.ids[1, 4] & cna_status == 'deleted')[, 
	 c('chrom', 'start', 'end', 'pair_barcode', 'copy_number', 'ploidy')]

 
 resAmp <- overlapRegions(regionAAmp, regionBAmp, min.pctA=30, only.boolean=TRUE) # 
 resDel <- overlapRegions(regionADel, regionBDel, min.pctA=30, only.boolean=TRUE) # 
 
 regionAAmp$score <- ifelse(resAmp==TRUE, 0, 1)
 regionADel$score <- ifelse(resDel==TRUE, 0, 1)
 
 res <- rbind(regionAAmp, regionADel)
 res <- res[, c('pair_barcode', 'chrom', 'start', 'end', 'copy_number', 'ploidy', 'score')]
 colnames(res)<- c("Sample", "Chromosome", "Start", "End", "Modal_Total_CN", "ploidy", "score")
 
 return(res)
})

cna.subclonality <- do.call(rbind, cna.subclonality)

library('GenomicRanges')
glass.seg.call.grl <- makeGRangesListFromDataFrame(cna.subclonality, split.field="Sample", keep.extra.columns=TRUE) 
		
saveRDS(glass.seg.call.grl, file='/data/glass_seg_call_grl.rds')


 
