
TCGAClinicalData <- function(infile)
{
	library('readxl')
	
	cli.data <- read_xlsx(path=infile, sheet = 1)
	cli.data <- as.data.frame(cli.data)
	rownames(cli.data) <- cli.data$bcr_patient_barcode

	cli.data <- cli.data[, c('type', 'age_at_initial_pathologic_diagnosis', 'gender', 'race', 
	 'ajcc_pathologic_tumor_stage', 'clinical_stage', 'histological_type', 'histological_grade', 
	 'OS', 'OS.time', 'DSS', 'DSS.time', 'DFI', 'DFI.time', 'PFI', 'PFI.time')]
	colnames(cli.data) <- c('cancer_type', 'age', 'gender', 'race', 'ajcc_stage', 'clinical_stage', 'histological_type', 
	 'histological_grade', 'os', 'os_time', 'dss', 'dss_time', 'dfi', 'dfi_time', 'pfi', 'pfi_time')

	# race: [Not Available], [Not Evaluated], [Unknown]
	cli.data$race[cli.data$race %in% c('[Not Available]', '[Not Evaluated]', '[Unknown]')] <- NA
	cli.data$race[cli.data$race %in% c('BLACK OR AFRICAN AMERICAN')] <- 'BLACK'
	cli.data$race[cli.data$race %in% c('AMERICAN INDIAN OR ALASKA NATIVE', 
	 'ASIAN', 'NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER')] <- 'OTHER'

	# ajcc_stage: [Discrepancy], [Not Available], [Not Applicable], [Unknown], 'Stage X', 'I/II NOS'
	cli.data$ajcc_stage[cli.data$ajcc_stage %in% c('[Discrepancy]', '[Not Available]', 
	 '[Not Applicable]', '[Unknown]', 'Stage X', 'I/II NOS')] <- NA
	cli.data$ajcc_stage[cli.data$ajcc_stage %in% c('IS', 'Stage IA', 'Stage IB')] <- 'Stage I'
	cli.data$ajcc_stage[cli.data$ajcc_stage %in% c('Stage IIA', 'Stage IIB', 'Stage IIC')] <- 'Stage II'
	cli.data$ajcc_stage[cli.data$ajcc_stage %in% c('Stage IIIA', 'Stage IIIB', 'Stage IIIC')] <- 'Stage III'
	cli.data$ajcc_stage[cli.data$ajcc_stage %in% c('Stage IVA', 'Stage IVB', 'Stage IVC')] <- 'Stage IV'

	# clinical_stage: [Discrepancy], [Not Applicable], [Not Available]
	cli.data$clinical_stage[cli.data$clinical_stage %in% c('[Discrepancy]', '[Not Applicable]', '[Not Available]')] <- NA
	cli.data$clinical_stage[cli.data$clinical_stage %in% c('I', 'Stage IA', 'Stage IA1', 'Stage IA2', 'Stage IB', 
	 'Stage IB1', 'Stage IB2', 'Stage IC', 'Stage IS')] <- 'Stage I'
	cli.data$clinical_stage[cli.data$clinical_stage %in% c('IIa', 'IIb', 'Stage IIA', 'Stage IIA1', 'Stage IIA2', 
	 'Stage IIB', 'Stage IIC')] <- 'Stage II'
	cli.data$clinical_stage[cli.data$clinical_stage %in% c('III', 'Stage IIIA', 'Stage IIIB', 'Stage IIIC', 
	 'Stage IIIC1', 'Stage IIIC2')] <- 'Stage III'
	cli.data$clinical_stage[cli.data$clinical_stage %in% c('IVa', 'IVb', 'Stage IVA', 'Stage IVB', 'Stage IVC')] <- 'Stage IV'

	# histological_grade: [Discrepancy], [Not Available], [Unknown], GX, GB
	cli.data$histological_grade[cli.data$histological_grade %in% c('[Discrepancy]', '[Not Available]', '[Unknown]', 'GX', 'GB')] <- NA
	cli.data$histological_grade[cli.data$histological_grade %in% c('Low Grade')] <- 'G1'
	cli.data$histological_grade[cli.data$histological_grade %in% c('High Grade')] <- 'G3'
	cli.data$histological_grade[cli.data$cancer_type == 'GBM'] <- 'G4'

	return(cli.data)
}

in.file <- '/data/OriginalData/TCGA-CDR-SupplementalTableS1.xlsx'
tcga.cli.data <- TCGAClinicalData(infile=in.file)
saveRDS(tcga.cli.data, file='/data/tcga_cli_data.rds') 


