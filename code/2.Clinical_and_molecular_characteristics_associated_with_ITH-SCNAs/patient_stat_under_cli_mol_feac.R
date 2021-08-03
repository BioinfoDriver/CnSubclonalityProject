# load data
gli.cn.alt.frac <- readRDS(file='/data/gli_glod_cn_alt_frac.rds')

library(dplyr)
library(ggplot2)
gli.cn.alt.frac <- gli.cn.alt.frac %>% mutate(
	kps_group = as.character(cut_width(karnofsky_performance_score, width = 10, boundary = 10)),
	kps_group = ifelse(kps_group %in% c('[20,30]', '(30,40]', '(40,50]', '(50,60]', '(60,70]'), "<=70", kps_group))

# > aggregate(age ~ IDH_CODEL_SUBTYPE, data=gli.cn.alt.frac, median)
  # IDH_CODEL_SUBTYPE age
# 1             IDHwt  61
# 2  IDHmut-non-codel  36
# 3      IDHmut-codel  44

# 
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

CliMolPatientSta <- function(dat, subtype, cancer.type, cli.mol.feacs){
 
 dat <- subset(dat, IDH_CODEL_SUBTYPE %in% subtype & cancer_type %in% cancer.type)
 dat <- dat %>% mutate(age_group = ifelse(age >= median(age),  'elder', 'younger'))
 
 stat.list <- lapply(cli.mol.feacs, function(cli.mol.feac){
  
  char.sta <- table(dat[, cli.mol.feac])
  char <- names(char.sta)[char.sta>5]
  dat.del <- subset(dat, get(cli.mol.feac) %in% char)
  
  stat <- CliSamCountStat(cli.mol.feac, dat.del)[[1]]
  return(stat)
 })
  
 return(stat.list)
}


codel.feacs <- c("age_group", "gender", "histological_grade", "tumor_sample_procurement_method", "kps_group", "laterality", 
 "family_history_of_cancer", "tumor_location", "first_presenting_symptom", 'TERT_EXPRESSION_STATUS', 'TRANSCRIPTOME_SUBTYPE')
 
non.codel.feacs <- c("age_group", "gender", "race", "histological_grade", "molecular_histological_type",
 "tumor_sample_procurement_method", "kps_group", "laterality", "family_history_of_cancer", "tumor_location", 
 "first_presenting_symptom", 'MGMT_PROMOTER_STATUS', 'ATRX_STATUS', "TERT_PROMOTER_STATUS", 
 'TERT_EXPRESSION_STATUS', 'TELOMERE_MAINTENANCE', 'TRANSCRIPTOME_SUBTYPE', 'SUPERVISED_DNA_METHYLATION_CLUSTER')
 
wt.feacs <- c("age_group", "gender", "race", "histological_grade", "molecular_histological_type",
 "tumor_sample_procurement_method", "kps_group", "laterality", "family_history_of_cancer", "tumor_location", 
 "first_presenting_symptom", 'CHR_7_GAIN_CHR_10_LOSS', 'CHR_19_20_CO_GAIN', 'MGMT_PROMOTER_STATUS', 
 'ATRX_STATUS', "TERT_PROMOTER_STATUS", 'TERT_EXPRESSION_STATUS','TELOMERE_MAINTENANCE', 
 'TRANSCRIPTOME_SUBTYPE', 'SUPERVISED_DNA_METHYLATION_CLUSTER')
 
CliMolPatientSta(gli.cn.alt.frac, 'IDHmut-codel', c('LGG', 'GBM'), codel.feacs)
CliMolPatientSta(gli.cn.alt.frac, 'IDHmut-non-codel', c('LGG', 'GBM'), non.codel.feacs)
CliMolPatientSta(gli.cn.alt.frac, 'IDHwt', c('LGG', 'GBM'), wt.feacs)


