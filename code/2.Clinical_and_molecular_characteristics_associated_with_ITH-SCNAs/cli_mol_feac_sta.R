
# load data
gli.cn.alt.frac <- readRDS(file='/data/gli_glod_cn_alt_frac.rds')

library(dplyr)
library(ggplot2)
gli.cn.alt.frac <- gli.cn.alt.frac %>% mutate(
	kps_group = as.character(cut_width(karnofsky_performance_score, width = 10, boundary = 10)),
	kps_group = ifelse(kps_group %in% c('[20,30]', '(30,40]', '(40,50]', '(50,60]', '(60,70]'), "<=70", kps_group))

gli.cn.alt.frac$CDKN2AB <- ifelse(gli.cn.alt.frac$CDKN2AB == 1, 'Deletion', 'Diploid')
gli.cn.alt.frac$EGFR <- ifelse(gli.cn.alt.frac$EGFR == 1, 'Amplification', 'Diploid')


# median, range
CliMolFeactureSta <- function(dat, subtype, features, cli.mol.feacs, filter.num){
 
 dat <- subset(dat, Integrated_Diagnoses %in% subtype)
 dat <- dat %>% mutate(age_group = ifelse(age >= median(age),  'elder', 'younger'))
 
 stat.list <- lapply(cli.mol.feacs, function(cli.mol.feac){
  
  char.sta <- table(dat[, cli.mol.feac])
  char <- names(char.sta)[char.sta >= filter.num]
  dat.del <- subset(dat, get(cli.mol.feac) %in% char)
  
  Merge <- function(x, y){

   by.by <- colnames(x)[1]
   z <- merge(x, y, by=by.by)

   return(z)
  }
  
  stat <- lapply(features, function(feac){
   
   st <- aggregate(get(feac) ~ get(cli.mol.feac), data=dat.del, summary)
   st <- cbind(st[, 1, FALSE], st[, 2][, c(2, 3, 5)])
   
   st[, 2] <- as.numeric(sprintf("%0.5f", st[, 2]))
   st[, 3] <- as.numeric(sprintf("%0.5f", st[, 3]))
   st[, 4] <- as.numeric(sprintf("%0.5f", st[, 4]))
   
   st$medianIQR <- paste(st[, 3], '(', st[, 2], ', ', st[, 4], ')', sep='')
   st <- st[, c(1, 5)]
   
   return(st)
  })
  
  stat <- Reduce(function(x, y) Merge(x, y), stat)
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
 

features <- c('clo_genome_frac', 'subclo_genome_frac')
CliMolFeactureSta(gli.cn.alt.frac, 'Oligodendroglioma,IDHmut-codel', features, codel.feacs, 10)
CliMolFeactureSta(gli.cn.alt.frac, 'Astrocytoma,IDHmut', features, non.codel.feacs, 10)
CliMolFeactureSta(gli.cn.alt.frac, 'Glioblastoma,IDHwt', features, wt.feacs, 10)

