
# load data
gli.cn.alt.frac <- readRDS(file='/data/gli_glod_cn_alt_frac.rds')

# age
library(dplyr)
library(ggplot2)
AgeSubcloGenomeFracCor <- function(dat){

 cor.pvalue <- lapply(split(dat, dat$Integrated_Diagnoses), function(subclo.genome.frac){
  
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

age.pvalues <- AgeSubcloGenomeFracCor(gli.cn.alt.frac)


# other features
gli.cn.alt.frac <- gli.cn.alt.frac %>%
  mutate(kps_group = as.character(cut_width(karnofsky_performance_score, width = 10, boundary = 10)),
		 kps_group = ifelse(kps_group %in% c('[20,30]', '(30,40]', '(40,50]', '(50,60]', '(60,70]'), "<=70", kps_group))

gli.cn.alt.frac$CDKN2AB <- ifelse(gli.cn.alt.frac$CDKN2AB == 1, 'Deletion', 'Diploid')
gli.cn.alt.frac$EGFR <- ifelse(gli.cn.alt.frac$EGFR == 1, 'Amplification', 'Diploid')


CliMolAssoSubcloBurden <- function(dat, features, subtype, cli.mol.feacs, filter.num){
 dat <- subset(dat, Integrated_Diagnoses %in% subtype)
 
 feac.pvalues <- lapply(features, function(feac){
  
  pvalues <- sapply(cli.mol.feacs, function(cli.mol.feac){
   
   if(cli.mol.feac == 'age'){
  
   subclo.genome.frac <- dat %>% mutate(age_group_median = ifelse(age >= median(age),  'elder', 'younger'))
   p <- wilcox.test(get(feac) ~ age_group_median, subclo.genome.frac)$p.value   
	
   }else{
   
    num <- table(dat[, cli.mol.feac])
    char <- names(num)[num >= filter.num]
    dat.del <- subset(dat, get(cli.mol.feac) %in% char)
    
    if(n_distinct(dat.del[, cli.mol.feac])<=1){
     p <- NA
	 
    }else if(n_distinct(dat.del[, cli.mol.feac])>2){
      p <- kruskal.test(get(feac) ~ get(cli.mol.feac), dat.del)$p.value
	 
    }else{
      p <- wilcox.test(get(feac) ~ get(cli.mol.feac), dat.del)$p.value
    
    }
   }
   
   return(p)
  
  })
  
 qvalues <- p.adjust(pvalues, 'fdr') 
 pqvalues <- data.frame(pvalues, qvalues)
 rownames(pqvalues) <- cli.mol.feacs
 return(pqvalues)
 
 })
 
 feac.pvalues <- do.call(cbind, feac.pvalues)
 colnames(feac.pvalues) <- paste0(rep(features, each=2), c('_pvalue', '_qvalue'))
 return(feac.pvalues)
}

features <- c('subclo_genome_frac', 'clo_genome_frac')


codel.feacs <- c("age", "gender", "histological_grade", "tumor_sample_procurement_method", "kps_group", "laterality", 
 "family_history_of_cancer", "tumor_location", "first_presenting_symptom", 'TERT_EXPRESSION_STATUS', 'TRANSCRIPTOME_SUBTYPE')
 
non.codel.feacs <- c("age", "gender", "histological_grade", "tumor_sample_procurement_method", "kps_group", "laterality", 
 "family_history_of_cancer", "tumor_location", "first_presenting_symptom", 'CDKN2AB', 'ATRX_STATUS', 'TERT_EXPRESSION_STATUS', 
 'MGMT_PROMOTER_STATUS', 'TRANSCRIPTOME_SUBTYPE', 'SUPERVISED_DNA_METHYLATION_CLUSTER')
 
wt.feacs <- c("age", "gender", "race", "tumor_sample_procurement_method", "kps_group", "laterality", 
 "family_history_of_cancer", "tumor_location", "first_presenting_symptom", 'CHR_7_GAIN_CHR_10_LOSS', 'CDKN2AB', 'EGFR', 
 'TERT_EXPRESSION_STATUS', 'CHR_19_20_CO_GAIN', 'MGMT_PROMOTER_STATUS', 'TRANSCRIPTOME_SUBTYPE', 'SUPERVISED_DNA_METHYLATION_CLUSTER')
 

codel.pvalues <- CliMolAssoSubcloBurden(gli.cn.alt.frac, features, 'Oligodendroglioma,IDHmut-codel', codel.feacs, 10)
non.codel.pvalues <- CliMolAssoSubcloBurden(gli.cn.alt.frac, features, 'Astrocytoma,IDHmut', non.codel.feacs, 10)
wt.pvalues <- CliMolAssoSubcloBurden(gli.cn.alt.frac, features, 'Glioblastoma,IDHwt', wt.feacs, 10)

