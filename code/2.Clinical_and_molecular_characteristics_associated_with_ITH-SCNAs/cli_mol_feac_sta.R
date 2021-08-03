
# load data
gli.cn.alt.frac <- readRDS(file='/pub5/xiaoyun/Jobs/J22/CopyNumberClonalityProject/Resource/CuratedData/gli_glod_cn_alt_frac.rds')

library(dplyr)
library(ggplot2)
gli.cn.alt.frac <- gli.cn.alt.frac %>% mutate(
	kps_group = as.character(cut_width(karnofsky_performance_score, width = 10, boundary = 10)),
	kps_group = ifelse(kps_group %in% c('[20,30]', '(30,40]', '(40,50]', '(50,60]', '(60,70]'), "<=70", kps_group))


# median, range
CliMolFeactureSta <- function(dat, subtype, cancer.type, features, cli.mol.feacs){
 
 dat <- subset(dat, IDH_CODEL_SUBTYPE %in% subtype & cancer_type %in% cancer.type)
 dat <- dat %>% mutate(age_group = ifelse(age >= median(age),  'elder', 'younger'))
 
 stat.list <- lapply(cli.mol.feacs, function(cli.mol.feac){
  
  char.sta <- table(dat[, cli.mol.feac])
  char <- names(char.sta)[char.sta>5]
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
 
non.codel.feacs <- c("age_group", "gender", "race", "histological_grade", "molecular_histological_type",
 "tumor_sample_procurement_method", "kps_group", "laterality", "family_history_of_cancer", "tumor_location", 
 "first_presenting_symptom", 'MGMT_PROMOTER_STATUS', 'ATRX_STATUS', "TERT_PROMOTER_STATUS", 
 'TERT_EXPRESSION_STATUS', 'TELOMERE_MAINTENANCE', 'TRANSCRIPTOME_SUBTYPE', 'SUPERVISED_DNA_METHYLATION_CLUSTER')
 
wt.feacs <- c("age_group", "gender", "race", "histological_grade", "molecular_histological_type",
 "tumor_sample_procurement_method", "kps_group", "laterality", "family_history_of_cancer", "tumor_location", 
 "first_presenting_symptom", 'CHR_7_GAIN_CHR_10_LOSS', 'CHR_19_20_CO_GAIN', 'MGMT_PROMOTER_STATUS', 
 'ATRX_STATUS', "TERT_PROMOTER_STATUS", 'TERT_EXPRESSION_STATUS','TELOMERE_MAINTENANCE', 
 'TRANSCRIPTOME_SUBTYPE', 'SUPERVISED_DNA_METHYLATION_CLUSTER')
 

features <- c('clo_genome_frac', 'subclo_genome_frac')
CliMolFeactureSta(gli.cn.alt.frac, 'IDHmut-codel', c('LGG', 'GBM'), features, codel.feacs)
CliMolFeactureSta(gli.cn.alt.frac, 'IDHmut-non-codel', c('LGG', 'GBM'), features, non.codel.feacs)
CliMolFeactureSta(gli.cn.alt.frac, 'IDHwt', c('LGG', 'GBM'), features, wt.feacs)

