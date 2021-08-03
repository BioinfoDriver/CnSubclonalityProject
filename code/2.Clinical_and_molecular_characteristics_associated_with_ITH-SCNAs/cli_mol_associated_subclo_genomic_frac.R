
# load data
gli.cn.alt.frac <- readRDS(file='/pub5/xiaoyun/Jobs/J22/CopyNumberClonalityProject/Resource/CuratedData/gli_glod_cn_alt_frac.rds')

# age
library(dplyr)
library(ggplot2)
AgeSubcloGenomeFracCor <- function(dat, cancer.type){

 dat <- subset(dat, cancer_type %in% cancer.type)
 cor.pvalue <- lapply(split(dat, dat$IDH_CODEL_SUBTYPE), function(subclo.genome.frac){
  
  features <- c('subclo_genome_frac', 'clo_genome_frac')
 
  pvalues.list <- sapply(features, function(feac){
  res <- cor.test(~ age + get(feac), data = subclo.genome.frac, method = 'spearman', exact = FALSE)
  res <- c(res$estimate, res$p.value)
  names(res) <- c('cor', 'cor.p.value')
  
  subclo.genome.frac <- subclo.genome.frac %>% mutate(age_group_median = ifelse(age >= median(age),  'elder', 'younger'))
  median.p.value <- setNames(wilcox.test(get(feac) ~ age_group_median, subclo.genome.frac)$p.value, 'towgroup.pvalue')
  
  
  subclo.genome.frac <- subclo.genome.frac %>% 
   mutate(age_group = as.character(cut_width(age, width = 10, boundary = 10)), 
          age_group = ifelse(age_group %in% c("[10,20]", "(20,30]"), "<=30", age_group), 
          age_group = ifelse(age_group %in% c("(70,80]", "(80,90]"), ">70", age_group))
 
  group.p.value <- setNames(kruskal.test(get(feac) ~ age_group, subclo.genome.frac)$p.value, 'multigroup.pvalue')
  
  return(c(res, median.p.value, group.p.value))
  })
  
  return(pvalues.list)
 })
 
 return(cor.pvalue)
}

AgeSubcloGenomeFracCor(gli.cn.alt.frac, c('LGG', 'GBM'))


# other features
gli.cn.alt.frac <- gli.cn.alt.frac %>%
  mutate(kps_group = as.character(cut_width(karnofsky_performance_score, width = 10, boundary = 10)),
		 kps_group = ifelse(kps_group %in% c('[20,30]', '(30,40]', '(40,50]', '(50,60]', '(60,70]'), "<=70", kps_group))


CliMolAssoSubcloBurden <- function(dat, features, subtype, cancer.type, cli.mol.feacs){
 dat <- subset(dat, IDH_CODEL_SUBTYPE %in% subtype & cancer_type %in% cancer.type)
 
 feac.pvalues <- lapply(features, function(feac){
  pvalues <- lapply(cli.mol.feacs, function(cli.mol.feac){
   
   num <- table(dat[, cli.mol.feac])
   char <- names(num)[num>5]
   dat.del <- subset(dat, get(cli.mol.feac) %in% char)
  
   if(n_distinct(dat.del[, cli.mol.feac])<=1){
    p <- NA
	
   }else if(n_distinct(dat.del[, cli.mol.feac])>2){
     p <- kruskal.test(get(feac) ~ get(cli.mol.feac), dat.del)$p.value
	
   }else{
     p <- wilcox.test(get(feac) ~ get(cli.mol.feac), dat.del)$p.value
   
   }
   return(p)
  
  })
  
  
 names(pvalues) <- cli.mol.feacs
 return(pvalues)
 
 })
 feac.pvalues <- do.call(cbind, feac.pvalues)
 colnames(feac.pvalues) <- features
 return(feac.pvalues)
}

features <- c('subclo_genome_frac', 'clo_genome_frac')

codel.feacs <- c("gender", "race", "histological_grade", "tumor_sample_procurement_method", 
 "kps_group", "laterality", "family_history_of_cancer", "tumor_location", "first_presenting_symptom",  
 'ATRX_STATUS', "TERT_PROMOTER_STATUS", 'TERT_EXPRESSION_STATUS','TELOMERE_MAINTENANCE', 'TRANSCRIPTOME_SUBTYPE')
 
non.codel.feacs <- c("gender", "race", "histological_grade", "molecular_histological_type",
 "tumor_sample_procurement_method", "kps_group", "laterality", "family_history_of_cancer", "tumor_location", 
 "first_presenting_symptom", 'MGMT_PROMOTER_STATUS', 'ATRX_STATUS', "TERT_PROMOTER_STATUS", 
 'TERT_EXPRESSION_STATUS', 'TELOMERE_MAINTENANCE', 'TRANSCRIPTOME_SUBTYPE', 'SUPERVISED_DNA_METHYLATION_CLUSTER')
 
wt.feacs <- c("gender", "race", "histological_grade", "molecular_histological_type",
 "tumor_sample_procurement_method", "kps_group", "laterality", "family_history_of_cancer", "tumor_location", 
 "first_presenting_symptom", 'CHR_7_GAIN_CHR_10_LOSS', 'CHR_19_20_CO_GAIN', 'MGMT_PROMOTER_STATUS', 
 'ATRX_STATUS', "TERT_PROMOTER_STATUS", 'TERT_EXPRESSION_STATUS','TELOMERE_MAINTENANCE', 
 'TRANSCRIPTOME_SUBTYPE', 'SUPERVISED_DNA_METHYLATION_CLUSTER')
 
codel.pvalues <- CliMolAssoSubcloBurden(gli.cn.alt.frac, features, 'IDHmut-codel', c('LGG', 'GBM'), codel.feacs)
non.codel.pvalues <- CliMolAssoSubcloBurden(gli.cn.alt.frac, features, 'IDHmut-non-codel', c('LGG', 'GBM'), non.codel.feacs)
wt.pvalues <- CliMolAssoSubcloBurden(gli.cn.alt.frac, features, 'IDHwt', c('LGG', 'GBM'), wt.feacs)
