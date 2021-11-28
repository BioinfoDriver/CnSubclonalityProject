
# load clinical characteristics
setwd('/data/OriginalData')
pan.can.cli.data <- read.csv("clinical_PANCAN_patient_with_followup.tsv", header=TRUE, stringsAsFactors=FALSE, 
 sep='\t', na.strings = c("", "[Not Available]", "[Not Applicable]", '[Unknown]', 
 '[Discrepancy]', '[Not Evaluated]', 'Not listed in Medical Record'))
pan.glio.cli.data <- subset(pan.can.cli.data, acronym %in% c('LGG', 'GBM')) # 1111

# common clinical characteristics
analyzed.cli.fea <- c('karnofsky_performance_score', 'laterality', 'family_history_of_cancer', 'tumor_location', 
 'first_presenting_symptom')
pan.glio.cli.data <- pan.glio.cli.data[, c("bcr_patient_barcode", analyzed.cli.fea)]


# load clinical characteristics-(extent of surgical resection)
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
pan.glio.mol.data <- read.csv(file='data_clinical_sample.txt', sep='\t', header=TRUE, comment.char="#", 
 stringsAsFactors=FALSE, na.strings = "NA")

# salient molecular features
analyzed.mol.fea <- c('IDH_STATUS', 'TERT_PROMOTER_STATUS', 'TERT_EXPRESSION_STATUS', 'ATRX_STATUS', 
 'TELOMERE_MAINTENANCE', 'MGMT_PROMOTER_STATUS', 'CHR_7_GAIN_CHR_10_LOSS', 'CHR_19_20_CO_GAIN', 
 'IDH_1P19Q_SUBTYPE', 'IDH_CODEL_SUBTYPE', 'TRANSCRIPTOME_SUBTYPE', 'SUPERVISED_DNA_METHYLATION_CLUSTER')
  
pan.glio.mol.data <- pan.glio.mol.data[, c("PATIENT_ID", analyzed.mol.fea)]


# load cinical prognostic data
pan.can.pro.data <- readRDS(file='/data/tcga_cli_data.rds') 
pan.glio.pro.data <- subset(pan.can.pro.data, cancer_type %in% c('GBM', 'LGG')) # In total, 1111 samples(GBM: 596; LGG: 515)  
pan.glio.pro.data <- pan.glio.pro.data[, -c(5:6)]


# CDKN2A/B Deletion, EGFR amplification
gene.cnv.alt <- read.table(file = 'all_thresholded.by_genes.txt', sep = '\t', header = TRUE, stringsAsFactors = FALSE)
rownames(gene.cnv.alt) <- gene.cnv.alt$Gene.Symbol
gene.cnv.alt <- gene.cnv.alt[, -c(1:3)]
colnames(gene.cnv.alt) <- gsub('\\.', '-', substr(colnames(gene.cnv.alt), 1, 12))
gene.cnv.alt <- as.data.frame(t(gene.cnv.alt))
gene.cnv.alt$CDKN2AB <- ifelse((gene.cnv.alt$CDKN2A == -2) | (gene.cnv.alt$CDKN2B == -2), 1, 0)
gene.cnv.alt$EGFR <- ifelse(gene.cnv.alt$EGFR == 2, 1, 0)

## merge
pan.glio.cli.mol.data <- merge(pan.glio.cli.data, pan.glio.mol.data, by.x='bcr_patient_barcode', by.y='PATIENT_ID')
pan.glio.cli.mol.data <- cbind(pan.glio.cli.mol.data, pan.glio.pro.data[pan.glio.cli.mol.data$bcr_patient_barcode, ])
pan.glio.cli.mol.data <- merge(pan.glio.cli.mol.data, pan.glio.esr.data, by='bcr_patient_barcode', all.x=TRUE)

pan.glio.cli.mol.data$CDKN2AB <- gene.cnv.alt[pan.glio.cli.mol.data$bcr_patient_barcode, 'CDKN2AB']
pan.glio.cli.mol.data$EGFR <- gene.cnv.alt[pan.glio.cli.mol.data$bcr_patient_barcode, 'EGFR']

# Histological molecular subtypes WHO 2016
pan.glio.cli.mol.data$mol_his_type <- sapply(seq(nrow(pan.glio.cli.mol.data)), function(index){
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


pan.glio.cli.mol.data$MOLECULAR_SUBTYPE <- pan.glio.cli.mol.data$IDH_CODEL_SUBTYPE
pan.glio.cli.mol.data$MOLECULAR_SUBTYPE[pan.glio.cli.mol.data$MOLECULAR_SUBTYPE == 'IDHmut-non-codel'] <- 'IDHmut'
pan.glio.cli.mol.data$MOLECULAR_SUBTYPE[pan.glio.cli.mol.data$MOLECULAR_SUBTYPE == 'IDHmut-codel' & 
 pan.glio.cli.mol.data$histological_grade == 'G4'] <- 'IDHmut'

pan.glio.cli.mol.data$MOLECULAR_SUBTYPE[pan.glio.cli.mol.data$histological_grade == 'G4' & 
 pan.glio.cli.mol.data$IDH_STATUS == 'Mutant'] <- 'IDHmut'
pan.glio.cli.mol.data$MOLECULAR_SUBTYPE[pan.glio.cli.mol.data$histological_grade == 'G4' & 
 pan.glio.cli.mol.data$IDH_STATUS == 'WT'] <- 'IDHwt'


# WHO 2021
# glioblastomas will comprise only IDH wild-type tumors
gbm.mol.his.modified <- pan.glio.cli.mol.data$mol_his_type == 'Glioblastoma' & 
 pan.glio.cli.mol.data$MOLECULAR_SUBTYPE == 'IDHmut'
 
pan.glio.cli.mol.data$mol_his_type[gbm.mol.his.modified] <- 'Astrocytoma'

# IDH wildtype diffuse astrocytic tumors
astro.mol.his.modified <- (pan.glio.cli.mol.data$mol_his_type == 'Astrocytoma' & 
 pan.glio.cli.mol.data$MOLECULAR_SUBTYPE == 'IDHwt') & 
 (pan.glio.cli.mol.data$CHR_7_GAIN_CHR_10_LOSS == 'Gain chr 7 & loss chr 10' |
 pan.glio.cli.mol.data$TERT_PROMOTER_STATUS == 'Mutant' | pan.glio.cli.mol.data$EGFR == 1)

pan.glio.cli.mol.data$mol_his_type[astro.mol.his.modified] <- 'Glioblastoma'
pan.glio.cli.mol.data$histological_grade[astro.mol.his.modified] <- 'G4'

# IDH-mutant diffuse astrocytic tumors are graded as 2, 3, or 4
astro.grade.modified <- pan.glio.cli.mol.data$mol_his_type == 'Astrocytoma' & 
 pan.glio.cli.mol.data$MOLECULAR_SUBTYPE == 'IDHmut' & pan.glio.cli.mol.data$CDKN2AB == 1
pan.glio.cli.mol.data$histological_grade[astro.grade.modified] <- 'G4'

# Integrated diagnosis (combined tissue-based histological and molecular diagnosis)
pan.glio.cli.mol.data$Integrated_Diagnoses <- apply(pan.glio.cli.mol.data, 1, function(x){
 
 if(is.na(x['mol_his_type'])){
	return(paste(paste0(x['histological_type'], ',', 'NOS')))
 
 }else{
	return(paste0(x['mol_his_type'], ',', x['MOLECULAR_SUBTYPE']))
 }
})

pan.glio.cli.mol.data$Integrated_Diagnoses[pan.glio.cli.mol.data$Integrated_Diagnoses == 'Glioblastoma,NA'] <- 'Glioblastoma,NOS'
pan.glio.cli.mol.data$Integrated_Diagnoses[pan.glio.cli.mol.data$Integrated_Diagnoses == 'Astrocytoma,IDHwt'] <- 'Astrocytoma,NEC'
pan.glio.cli.mol.data$Integrated_Diagnoses[pan.glio.cli.mol.data$Integrated_Diagnoses == 'Astrocytoma,NEC' & 
 (is.na(pan.glio.cli.mol.data$TERT_PROMOTER_STATUS) | is.na(pan.glio.cli.mol.data$EGFR) | 
 is.na(pan.glio.cli.mol.data$CHR_7_GAIN_CHR_10_LOSS))] <- 'Astrocytoma,NOS'


# save
pan.glio.cli.mol.data <- pan.glio.cli.mol.data[, c(1, 19:24, 16, 36:38, 7:8, 34:35, 13, 15, 9:12, 14, 17:18, 2:6, 33, 25:32)]
rownames(pan.glio.cli.mol.data) <- pan.glio.cli.mol.data$bcr_patient_barcode

pan.glio.cli.mol.data$histological_grade[pan.glio.cli.mol.data$Integrated_Diagnoses == 'Astrocytoma,IDHmut' & 
 pan.glio.cli.mol.data$histological_grade %in% c('G2', 'G3') & is.na(pan.glio.cli.mol.data$CDKN2AB)] <- NA

saveRDS(pan.glio.cli.mol.data, file='/data/tcga_glioma_cli_mol.rds')

