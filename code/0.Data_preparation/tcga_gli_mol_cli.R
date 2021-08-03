
# load clinical characteristics
setwd('/pub5/xiaoyun/Jobs/J22/EvoClass2.0/OriginalData/PanCanAtlas')
pan.can.cli.data <- read.csv("clinical_PANCAN_patient_with_followup.tsv", header=TRUE, stringsAsFactors=FALSE, 
 sep='\t', na.strings = c("", "[Not Available]", "[Not Applicable]", '[Unknown]', 
 '[Discrepancy]', '[Not Evaluated]', 'Not listed in Medical Record'))
pan.glio.cli.data <- subset(pan.can.cli.data, acronym %in% c('LGG', 'GBM')) # 1111

# common clinical characteristics
analyzed.cli.fea <- c('karnofsky_performance_score', 'laterality', 'family_history_of_cancer', 'tumor_location', 
 'first_presenting_symptom')

pan.glio.cli.data <- pan.glio.cli.data[, c("bcr_patient_barcode", analyzed.cli.fea)]


# load clinical characteristics-(extent of surgical resection)
setwd('/pub5/xiaoyun/Jobs/J22/CopyNumberClonalityProject/OriginalData/GDC')
sample <- read.csv("sample.tsv", header=TRUE, stringsAsFactors=FALSE, sep='\t')

gbm.esr.data <- read.csv("nationwidechildrens.org_ssf_tumor_samples_gbm.txt", header=TRUE, stringsAsFactors=FALSE, 
 sep='\t', na.strings = c('[Not Available]')) 
gbm.esr.data <- gbm.esr.data[-c(1:2), c('bcr_patient_barcode', 'bcr_sample_uuid', 'tumor_sample_procurement_method')]

lgg.esr.data <- read.csv("nationwidechildrens.org_ssf_tumor_samples_lgg.txt", header=TRUE, stringsAsFactors=FALSE, 
 sep='\t', na.strings = c('[Not Available]')) 
lgg.esr.data <- lgg.esr.data[-c(1:2), c('bcr_patient_barcode', 'bcr_sample_uuid', 'tumor_sample_procurement_method')]

pan.glio.esr.data <- rbind(gbm.esr.data, lgg.esr.data)
pan.glio.esr.data <- cbind(pan.glio.esr.data, 
 sample[match(tolower(pan.glio.esr.data$bcr_sample_uuid), sample$sample_id), 'sample_submitter_id', FALSE])
pan.glio.esr.data <- subset(pan.glio.esr.data, substr(sample_submitter_id, 14, 15)=='01')
pan.glio.esr.data <- pan.glio.esr.data[, -c(2, 4)]

# load molecular features
setwd('/pub5/xiaoyun/Jobs/J22/CopyNumberClonalityProject/OriginalData/cBioPortal')
pan.glio.mol.data <- read.csv(file='data_clinical_sample.txt', sep='\t', header=TRUE, comment.char="#", 
 stringsAsFactors=FALSE, na.strings = "NA")

# salient molecular features
analyzed.mol.fea <- c('IDH_STATUS', 'TERT_PROMOTER_STATUS', 'TERT_EXPRESSION_STATUS', 'ATRX_STATUS', 
 'TELOMERE_MAINTENANCE', 'MGMT_PROMOTER_STATUS', 'CHR_7_GAIN_CHR_10_LOSS', 'CHR_19_20_CO_GAIN', 
 'IDH_1P19Q_SUBTYPE', 'IDH_CODEL_SUBTYPE', 'TRANSCRIPTOME_SUBTYPE', 'SUPERVISED_DNA_METHYLATION_CLUSTER')
  
pan.glio.mol.data <- pan.glio.mol.data[, c("PATIENT_ID", analyzed.mol.fea)]


# load cinical prognostic data
pan.can.pro.data <- readRDS(file='/pub5/xiaoyun/Jobs/J22/CopyNumberClonalityProject/Resource/tcga_cli_data.rds') 
pan.glio.pro.data <- subset(pan.can.pro.data, cancer_type %in% c('GBM', 'LGG')) # In total, 1111 samples(GBM: 596; LGG: 515)  
pan.glio.pro.data <- pan.glio.pro.data[, -c(5:6)]


## merge
pan.glio.cli.mol.data <- merge(pan.glio.cli.data, pan.glio.mol.data, by.x='bcr_patient_barcode', by.y='PATIENT_ID')
pan.glio.cli.mol.data <- cbind(pan.glio.cli.mol.data, pan.glio.pro.data[pan.glio.cli.mol.data$bcr_patient_barcode, ])
pan.glio.cli.mol.data <- merge(pan.glio.cli.mol.data, pan.glio.esr.data, by='bcr_patient_barcode', all.x=TRUE)


# Histological molecular subtypes
pan.glio.cli.mol.data$molecular_histological_type <- sapply(seq(nrow(pan.glio.cli.mol.data)), function(index){
	his.mol.data <- pan.glio.cli.mol.data[index, ]
	if(his.mol.data$cancer_type == 'GBM'){
		return('Glioblastoma')
	
	}else{
		if(is.na(his.mol.data$IDH_CODEL_SUBTYPE)){
			return(NA)
		
		}else if(his.mol.data$IDH_CODEL_SUBTYPE == 'IDHmut-codel'){
			return('Oligodendroglioma')
		
		}else{
			return('Astrocytoma')
		}
	}
})

pan.glio.cli.mol.data <- pan.glio.cli.mol.data[, c(1, 20:22, 24, 23, 34, 33, 2:6, 25:32, 7:19)]
rownames(pan.glio.cli.mol.data) <- pan.glio.cli.mol.data$bcr_patient_barcode
saveRDS(pan.glio.cli.mol.data, file='/pub5/xiaoyun/Jobs/J22/CopyNumberClonalityProject/Resource/CuratedData/tcga_glioma_cli_mol.rds')

