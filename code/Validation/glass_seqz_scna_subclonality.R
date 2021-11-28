

setwd('/data/OriginalData/syn17038081')
glass.seqz.seg <- read.csv(file='variants_seqz_seg.tsv', header=TRUE, sep='\t', stringsAsFactors=FALSE)
glass.seqz.ploidy <- read.csv(file='variants_seqz_params.tsv', header=TRUE, sep='\t', stringsAsFactors=FALSE)

glass.seqz.seg <- merge(glass.seqz.seg, glass.seqz.ploidy[, c('pair_barcode', 'ploidy')], by='pair_barcode')

glass.seqz.seg$cna_status <- ifelse(glass.seqz.seg$copy_number == round(glass.seqz.seg$ploidy), 'neutral', 
 ifelse(glass.seqz.seg$copy_number > round(glass.seqz.seg$ploidy), 'amplified', 'deleted'))

glass.seqz.seg$pair_barcode <- paste0(substr(glass.seqz.seg$pair_barcode, 1, 18), substr(glass.seqz.seg$pair_barcode, 22, 29))


gold.set <- read.csv(file='analysis_gold_set.tsv', header=TRUE, sep='\t', stringsAsFactors=FALSE)
library(dplyr)
gold.set <- gold.set %>% mutate(ROW_ID = NULL, ROW_VERSION = NULL)


gold.set$tumor_barcode_a <- paste0(substr(gold.set$tumor_pair_barcode, 1, 19), substr(gold.set$tumor_barcode_a, 17, 23))

gold.set$tumor_barcode_b <- paste0(substr(gold.set$tumor_barcode_b, 1, 15), 
 substr(gold.set$tumor_pair_barcode, 16, 18), substr(gold.set$tumor_barcode_b, 16, 23))

gold.set$seq_method <- substr(gold.set$tumor_barcode_a, 24, 26)

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
gold.set.curated$seq_method <- substr(gold.set.curated$tumor_barcode_a, 24, 26)

}

gold.set <- rbind(gold.set, gold.set.curated)
gold.set <- subset(gold.set, tumor_barcode_a %in% glass.seqz.seg$pair_barcode & 
 tumor_barcode_b %in% glass.seqz.seg$pair_barcode)
 

library(regioneR)
glass.seqz.seg$chrom <- paste0('chr', glass.seqz.seg$chrom)
glass.seqz.seg$chrom[glass.seqz.seg$chrom == 'chr23'] <- 'chrX'

subclo.cna.perc <- lapply(seq(nrow(gold.set)), function(index){
 
 sam.ids <- gold.set[index, ]
 regionA <- subset(glass.seqz.seg, pair_barcode == sam.ids[1, 3] & cna_status != 'neutral')[, 
	 c('chrom', 'start', 'end', 'pair_barcode', 'copy_number', 'cna_status')]
 colnames(regionA) <- paste0(colnames(regionA), 'A')
 
 regionB <- subset(glass.seqz.seg, pair_barcode == sam.ids[1, 4] & cna_status != 'neutral')[, 
	 c('chrom', 'start', 'end', 'pair_barcode', 'copy_number', 'cna_status')]
 colnames(regionB) <- paste0(colnames(regionB), 'B')
 
 res <- overlapRegions(regionA, regionB, colA=4:6, colB=4:6, type="any", min.bases=1, get.bases=TRUE)
 
 # subclonal scna percentage
 equal.size <- sum(subset(res, cna_statusA == cna_statusB)$ov.bases)
 unequal.size <- sum(subset(res, cna_statusA != cna_statusB)$ov.bases)
 
 regionA.size <- sum(regionA$endA - regionA$startA +1)
 regionB.size <- sum(regionB$endB - regionB$startB +1)
 
 subclo.cna.perc <- 1-equal.size/(regionA.size+regionB.size-equal.size-unequal.size)
 
 subclo.cna.perc.t <- 1-equal.size/regionA.size
 subclo.cna.perc.r <- 1-equal.size/regionB.size
 
 
 # clonal/subclonal scna burden
 regionA <- subset(glass.seqz.seg, pair_barcode == sam.ids[1, 3])[, c('chrom', 'start', 'end')]
 colnames(regionA) <- paste0(colnames(regionA), 'A')
 seg.length.t <- sum(regionA$endA - regionA$startA +1)
 
 
 regionB <- subset(glass.seqz.seg, pair_barcode == sam.ids[1, 4])[, c('chrom', 'start', 'end')]
 colnames(regionB) <- paste0(colnames(regionB), 'B') 
 seg.length.r <- sum(regionB$endB - regionB$startB +1)
 
 res <- overlapRegions(regionA, regionB, type="any", min.bases=1, get.bases=TRUE)
 
 
 clonal.cna.bur <- equal.size/(seg.length.t + seg.length.r - sum(res$ov.bases))
 subclonal.cna.bur <- (regionA.size+regionB.size-2*equal.size-unequal.size)/(seg.length.t + seg.length.r - sum(res$ov.bases))
 
 
 clonal.cna.bur.t <- equal.size/seg.length.t
 clonal.cna.bur.r <- equal.size/seg.length.r
 
 
 subclonal.cna.bur.t <- (regionA.size-equal.size)/seg.length.t
 subclonal.cna.bur.r <- (regionB.size-equal.size)/seg.length.r
 
 
 res <- c(subclo.cna.perc, subclo.cna.perc.t, subclo.cna.perc.r, clonal.cna.bur, subclonal.cna.bur, 
	clonal.cna.bur.t, clonal.cna.bur.r, subclonal.cna.bur.t, subclonal.cna.bur.r)
 res <- setNames(res, c('subclo.cna.perc', 'subclo.cna.perc.t', 'subclo.cna.perc.r', 'clonal.cna.bur', 'subclonal.cna.bur', 
	'clonal.cna.bur.t', 'clonal.cna.bur.r', 'subclonal.cna.bur.t', 'subclonal.cna.bur.r'))

 return(res)
})

subclo.cna.perc <- do.call(rbind, subclo.cna.perc)

gold.set <- data.frame(gold.set, subclo.cna.perc)
gold.set$tumor_barcode_a <- substr(gold.set$tumor_barcode_a, 1, 15)
gold.set$tumor_barcode_b <- substr(gold.set$tumor_barcode_b, 1, 15)

# Integrated diagnosis (combined tissue-based histological and molecular diagnosis)
clinical.data <- read.csv(file='clinical_surgeries.tsv', header=TRUE, sep='\t', stringsAsFactors=FALSE)
clinical.data <- clinical.data %>% mutate(ROW_ID = NULL, ROW_VERSION = NULL, case_barcode = NULL)

glass.cna.ith <- merge(gold.set, clinical.data, by.x='tumor_barcode_a', by.y='sample_barcode')

glass.cna.ith$Integrated_Diagnoses <- NA
glass.cna.ith$Integrated_Diagnoses[glass.cna.ith$idh_codel_subtype == 'IDHmut-codel' & 
 glass.cna.ith$grade %in% c('II', 'III')] <- 'Oligodendroglioma,IDHmut-codel'

glass.cna.ith$Integrated_Diagnoses[glass.cna.ith$idh_codel_subtype == 'IDHmut-noncodel' | 
 (glass.cna.ith$histology == 'Glioblastoma' & glass.cna.ith$idh_codel_subtype == 'IDHmut-codel')] <- 'Astrocytoma,IDHmut'

glass.cna.ith$Integrated_Diagnoses[glass.cna.ith$idh_codel_subtype == 'IDHwt' & 
 glass.cna.ith$histology == 'Glioblastoma'] <- 'Glioblastoma,IDHwt'


saveRDS(glass.cna.ith, file='/data/glass_seqz_cna_ith.rds')


