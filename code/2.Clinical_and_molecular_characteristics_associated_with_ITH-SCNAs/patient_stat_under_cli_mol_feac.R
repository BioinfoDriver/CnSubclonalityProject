# load data
gli.cn.alt.frac <- readRDS(file='/data/gli_glod_cn_alt_frac.rds')

library(dplyr)
library(ggplot2)
gli.cn.alt.frac <- gli.cn.alt.frac %>% mutate(
	kps_group = as.character(cut_width(karnofsky_performance_score, width = 10, boundary = 10)),
	kps_group = ifelse(kps_group %in% c('[20,30]', '(30,40]', '(40,50]', '(50,60]', '(60,70]'), "<=70", kps_group))

gli.cn.alt.frac$CDKN2AB <- ifelse(gli.cn.alt.frac$CDKN2AB == 1, 'Deletion', 'Diploid')
gli.cn.alt.frac$EGFR <- ifelse(gli.cn.alt.frac$EGFR == 1, 'Amplification', 'Diploid')


# > aggregate(age ~ Integrated_Diagnoses, data=gli.cn.alt.frac, median)
            # Integrated_Diagnoses age
# 1             Astrocytoma,IDHmut  36
# 2             Glioblastoma,IDHwt  61
# 3 Oligodendroglioma,IDHmut-codel  44

# 样本数目统计
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

CliMolPatientSta <- function(dat, subtype, cli.mol.feacs, filter.num){
 
 dat <- subset(dat, Integrated_Diagnoses %in% subtype)
 dat <- dat %>% mutate(age_group = ifelse(age >= median(age),  'elder', 'younger'))
 
 stat.list <- lapply(cli.mol.feacs, function(cli.mol.feac){
  
  char.sta <- table(dat[, cli.mol.feac])
  char <- names(char.sta)[char.sta >= filter.num]
  dat.del <- subset(dat, get(cli.mol.feac) %in% char)
  
  stat <- CliSamCountStat(cli.mol.feac, dat.del)[[1]]
  return(stat)
 })
  
 return(stat.list)
}


codel.feacs <- c("age_group", "gender", "histological_grade", "tumor_sample_procurement_method", "kps_group", "laterality", 
 "family_history_of_cancer", "tumor_location", "first_presenting_symptom", 'TERT_EXPRESSION_STATUS', 'TRANSCRIPTOME_SUBTYPE')
 
non.codel.feacs <- c("age_group", "gender", "histological_grade", "tumor_sample_procurement_method", "kps_group", "laterality", 
 "family_history_of_cancer", "tumor_location", "first_presenting_symptom", 'CDKN2AB', 'ATRX_STATUS', 'TERT_EXPRESSION_STATUS', 
 'MGMT_PROMOTER_STATUS', 'TRANSCRIPTOME_SUBTYPE', 'SUPERVISED_DNA_METHYLATION_CLUSTER')
 
wt.feacs <- c("age_group", "gender", "race", "tumor_sample_procurement_method", "kps_group", "laterality", 
 "family_history_of_cancer", "tumor_location", "first_presenting_symptom", 'CHR_7_GAIN_CHR_10_LOSS', 'CDKN2AB', 'EGFR', 
 'TERT_EXPRESSION_STATUS', 'CHR_19_20_CO_GAIN', 'MGMT_PROMOTER_STATUS', 'TRANSCRIPTOME_SUBTYPE', 'SUPERVISED_DNA_METHYLATION_CLUSTER')
 
 
CliMolPatientSta(gli.cn.alt.frac, 'Oligodendroglioma,IDHmut-codel', codel.feacs, 10)
CliMolPatientSta(gli.cn.alt.frac, 'Astrocytoma,IDHmut', non.codel.feacs, 10)
CliMolPatientSta(gli.cn.alt.frac, 'Glioblastoma,IDHwt', wt.feacs, 10)


