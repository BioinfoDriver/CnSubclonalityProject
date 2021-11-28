

CliSamCountStat <- function(cli.features, cli.data){
  
  sam.count <- dim(cli.data)[1]
  cli.stat.result <- lapply(cli.features, function(sin.cli.fea){
	
	# Distinguish the variable as continuous or categorical
	cli.varibale <- cli.data[, sin.cli.fea]
    if(all(is.na(cli.varibale))){
      value <- NA
    }else{
      cli.varibale <- na.omit(cli.varibale)
	  cli.varibale <- as.character(cli.varibale)
      value <- all(!is.na(as.numeric(cli.varibale))) 
    }
	
	if(is.na(value)){
	cli.df <- as.matrix(data.frame(cli.feature = sin.cli.fea, cli.feature.stat = 'Unknown'))
	
	}else if(value == TRUE){
	# continuous, especially for age
		cli.varibale <- as.numeric(cli.varibale)
		cli.df <- as.matrix(data.frame(cli.feature = "Age", cli.feature.stat = paste(median(cli.varibale), '(',min(cli.varibale), '-', max(cli.varibale), ')', sep='')))
	
	}else if(value == FALSE){
	# categorical
	cli.table <- table(cli.data[, sin.cli.fea])
    cli.table.frac <- cli.table/sam.count
    unknown.count <- sam.count - sum(cli.table)
    unknown.frac <- 1-sum(cli.table.frac)
    cli.feature.count <- c(as.numeric(cli.table), unknown.count)
    cli.feature.frac <- c(round(cli.table.frac, 2)*100, round(unknown.frac, 2)*100)
    
    cli.fea.count.frac <- paste(cli.feature.count, '(', cli.feature.frac, ')', sep = '')
    
    cli.df <- as.matrix(data.frame(cli.feature = c(sin.cli.fea, names(cli.table), 'Unknown'), cli.feature.stat = c(sin.cli.fea, cli.fea.count.frac), stringsAsFactors = FALSE))

	}
	t(cli.df)
  })
  
  cli.stat.result
}

# Silver set
gli.cn.alt.frac <- readRDS(file='/data/tcga_gli_cn_alt_frac.rds')

silver.set <- rownames(gli.cn.alt.frac) 
pan.glio.cli.mol.data <- readRDS(file='/data/tcga_glioma_cli_mol.rds')
silver.set.cli.data <- pan.glio.cli.mol.data[substr(silver.set, 1, 12), ]

# Statistical analysis of clinical molecular characteristics on silver set
library('tidyverse')
silver.set.cli.data <- silver.set.cli.data %>%
  mutate(kps_group = as.character(cut_width(karnofsky_performance_score, width = 10, boundary = 10)),
		 kps_group = ifelse(kps_group %in% c('[20,30]', '(30,40]', '(40,50]', '(50,60]', '(60,70]'), "<=70", kps_group))

cli.features <- c("age", "gender", "race", "histological_grade", "tumor_sample_procurement_method", 
 "kps_group", "laterality", "family_history_of_cancer", "tumor_location", "first_presenting_symptom", 
 'TERT_PROMOTER_STATUS', 'CHR_7_GAIN_CHR_10_LOSS',  'CDKN2AB', 'EGFR', 'IDH_STATUS', 
 'IDH_1P19Q_SUBTYPE', 'IDH_CODEL_SUBTYPE', 'MOLECULAR_SUBTYPE', "mol_his_type",
 'ATRX_STATUS', 'TERT_EXPRESSION_STATUS','TELOMERE_MAINTENANCE', 'CHR_19_20_CO_GAIN', 
 'MGMT_PROMOTER_STATUS', 'TRANSCRIPTOME_SUBTYPE', 'SUPERVISED_DNA_METHYLATION_CLUSTER')


# silver.statistic <- CliSamCountStat(cli.features, silver.set.cli.data)

silver.set.cli.data$EGFR <- ifelse(silver.set.cli.data$EGFR == 1, 'YES', 'NO')
silver.set.cli.data$CDKN2AB <- ifelse(silver.set.cli.data$CDKN2AB == 1, 'YES', 'NO')

gold.set.cli.data <- subset(silver.set.cli.data, Integrated_Diagnoses %in% 
 c('Oligodendroglioma,IDHmut-codel', 'Astrocytoma,IDHmut', 'Glioblastoma,IDHwt'))
gold.statistic <- CliSamCountStat(cli.features, gold.set.cli.data)


# Statistical analysis of clinical molecular characteristics on each subgroups
subtype.name <- unique(gold.set.cli.data$Integrated_Diagnoses)
subtype.statistic <- lapply(subtype.name, function(per.type){
	per.subtype.characters <- subset(gold.set.cli.data, Integrated_Diagnoses == per.type)
	per.subtype.statistic <- CliSamCountStat(cli.features, per.subtype.characters)
	per.subtype.statistic
})
names(subtype.statistic) <- subtype.name

gold.set <- gold.set.cli.data$bcr_patient_barcode
saveRDS(gold.set, file='/data/gold_set.rds')

