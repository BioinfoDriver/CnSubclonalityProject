
# whitelisted samples
setwd('/data/OriginalData')
sam.qua.anno <- read.csv(file='merged_sample_quality_annotations.tsv', sep='\t', header=TRUE, stringsAsFactors=FALSE)

gli.sam.snp.anno <- subset(sam.qua.anno, cancer.type %in% c('LGG', 'GBM') & platform=='Genome_Wide_SNP_6')
# In total, 1111 samples(GBM: 598; LGG: 513) 
gli.sam.snp.anno <- subset(gli.sam.snp.anno, substr(aliquot_barcode, 14, 15) == '01')

# samples were cleaned based on the quality metrics
gli.sam.snp.anno <- subset(gli.sam.snp.anno, patient_annotation == '') # 81 samples were excluded (GBM: 65; LGG: 16)
gli.sam.snp.anno <- subset(gli.sam.snp.anno, sample_annotation == '') # 60 samples were excluded (GBM: 60)
gli.sam.snp.anno <- subset(gli.sam.snp.anno, aliquot_annotation == '') # 22 samples were excluded (GBM: 22)

# excluded Whole Genome Amplification (WGA) samples
gli.sam.snp.anno <- subset(gli.sam.snp.anno, substr(aliquot_barcode, 20,  20) != 'G') # 21 samples were excluded (GBM: 21)

# excluded technical replicates
gli.sam.snp.anno <- subset(gli.sam.snp.anno, !duplicated(patient_barcode)) # 9 samples were excluded (GBM: 9)


# The Pan-Cancer clinical data
tcga.cli.data <- readRDS(file='/data/tcga_cli_data.rds') 
gli.cli.data <- subset(tcga.cli.data, cancer_type %in% c('GBM', 'LGG')) # In total, 1111 samples(GBM: 596; LGG: 515)  
gli.cli.data <- subset(gli.cli.data, histological_type %in% c('Astrocytoma', 'Oligoastrocytoma', 
 'Oligodendroglioma', 'Untreated primary (de novo) GBM')) # 51 samples were excluded (GBM: 51)

gli.sam.snp.anno <- subset(gli.sam.snp.anno, patient_barcode %in% rownames(gli.cli.data)) # 31 samples were excluded (GBM: 30, LGG: 1)
gli.cli.data <- gli.cli.data[gli.sam.snp.anno$patient_barcode, ] # In total, 887 samples(GBM: 391; LGG: 496) 



# Tissue Source Site Codes
tss.code <- read.csv(file='tissue_source_site_codes.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
lgg.code <- subset(tss.code, Study.Name == 'Brain Lower Grade Glioma')$TSS.Code
gbm.code <- subset(tss.code, Study.Name == 'Glioblastoma multiforme')$TSS.Code


# ABSOLUTE purity/ploidy file
abs.call.dat <- read.csv(file='TCGA_mastercalls.abs_tables_JSedit.fixed.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
gli.abs.call.dat <- subset(abs.call.dat, substr(array, 6, 7) %in% c(lgg.code, gbm.code))
gli.abs.call.dat <- subset(gli.abs.call.dat, substr(array, 14, 15) == '01') # In total, 1077 samples(GBM: 566; LGG: 511) 
gli.abs.call.dat <- subset(gli.abs.call.dat, solution=='new') # 200 samples were excluded (GBM: 191; LGG: 9) 
gli.abs.call.dat <- subset(gli.abs.call.dat, !(is.na(purity) |is.na(ploidy))) # 7 samples were excluded (GBM: 5; LGG: 2)
# In total, 796 samples(GBM: 310; LGG: 486); 74 samples were excluded (GBM: 60; LGG: 14)
gli.abs.call.dat <- subset(gli.abs.call.dat, array %in% substr(gli.sam.snp.anno$aliquot_barcode, 1, 15))


# ABSOLUTE SCNA calls
abs.seg.dat <- read.csv(file='TCGA_mastercalls.abs_segtabs.fixed.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
gli.abs.seg.dat <- subset(abs.seg.dat, substr(Sample, 6, 7) %in% c(lgg.code, gbm.code))
gli.abs.seg.dat <- subset(gli.abs.seg.dat, substr(Sample, 14, 15) == '01') # In total, 1079 samples(GBM: 567; LGG: 512) 
gli.abs.seg.dat <- subset(gli.abs.seg.dat, solution == 'new') # 202 samples were excluded (GBM: 192; LGG: 10)
gli.abs.seg.dat <- subset(gli.abs.seg.dat, !is.na(Modal_Total_CN)) # 7 samples were excluded (GBM: 5; LGG: 2)
# In total, 796 samples(GBM: 310; LGG: 486); 74 samples were excluded (GBM: 60; LGG: 14)
gli.abs.seg.dat <- subset(gli.abs.seg.dat, Sample %in% substr(gli.sam.snp.anno$aliquot_barcode, 1, 15))


# save
write.table(gli.abs.call.dat, file='tcga_glioma_purity_ploidy.txt', sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(gli.abs.seg.dat, file='tcga_glioma_abs_seg.txt', sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)


