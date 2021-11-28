
# load data
gli.cn.alt.frac <- readRDS(file='/data/gli_glod_cn_alt_frac.rds')
library('tidyverse')
gli.cn.alt.frac <- gli.cn.alt.frac %>% mutate(
	kps_group = as.character(cut_width(karnofsky_performance_score, width = 10, boundary = 10)),
	kps_group = ifelse(kps_group %in% c('[20,30]', '(30,40]', '(40,50]', '(50,60]', '(60,70]'), "<=70", kps_group))

gli.cn.alt.frac$CDKN2AB <- ifelse(gli.cn.alt.frac$CDKN2AB == 1, 'Deletion', 'Diploid')

# Cox
source('/code/Function/Cox.function.R')

SCNAsBurdenSurAnalysis <- function(sur.dat, feacs, subtype, grade, disperse=FALSE){
 
 sur.dat <- subset(sur.dat, Integrated_Diagnoses %in% subtype & histological_grade %in% grade)
 
 cox.pvalues <- lapply(feacs, function(feac){
  
  cox.pvalue <- lapply(c('os', 'dss', 'pfi'), function(x){
   
   if(subtype == 'Astrocytoma,IDHmut'){
    
	if(length(grade) == 1){
	 
	 dat <- sur.dat[, c('bcr_patient_barcode', x, paste0(x, '_time'), feac, 'age', 'gender', 
	  'MGMT_PROMOTER_STATUS')]
	}else if(length(grade) != 3){
	
	 dat <- sur.dat[, c('bcr_patient_barcode', x, paste0(x, '_time'), feac, 'age', 'gender', 
	  'histological_grade', 'MGMT_PROMOTER_STATUS')]
	
	 dat <- subset(dat, !is.na(histological_grade))	
	
	}else{
	 
	 dat <- sur.dat[, c('bcr_patient_barcode', x, paste0(x, '_time'), feac, 'age', 'gender', 
	  'CDKN2AB', 'histological_grade', 'MGMT_PROMOTER_STATUS')]
	
	 dat <- subset(dat, !is.na(histological_grade))
	}
	dat <- subset(dat, !is.na(MGMT_PROMOTER_STATUS))
	
   }else{
    dat <- sur.dat[, c('bcr_patient_barcode', x, paste0(x, '_time'), feac, 'age', 'gender', 
	 'MGMT_PROMOTER_STATUS', 'kps_group')] 
    dat <- subset(dat, !is.na(MGMT_PROMOTER_STATUS) & !is.na(kps_group))
   }

   
   print(nrow(dat))
   
   if(disperse){
    dat[, feac] <- ifelse(dat[, feac] > median(dat[, feac]), 'High ITH', 'Low ITH')
   }else{
    dat[, feac] <- cut_width(x=dat[, feac], width=0.05, labels=FALSE)
   }
    
   dat <- subset(dat, !(is.na(get(x)) | is.na(get(paste0(x, '_time')))))  
   sur.times <- dat[, paste0(x, '_time')]
   sur.status <- dat[, x]
     
   cox.res <- Cox.function(sur.times, sur.status, dat)
   index <- grep(feac, cox.res$variate)
   
   if(disperse){
    cox.res <- cox.res[index+1, , FALSE]
   }else{
    cox.res <- cox.res[index, , FALSE]
   }

   cox.res$variate <- paste(feac, x, sep=' ')
   
   return(cox.res)
  })
  
  cox.pvalue <- do.call(rbind, cox.pvalue)
  return(cox.pvalue)
 })

 cox.pvalues <- do.call(rbind, cox.pvalues)
 return(cox.pvalues)
}

features <- c('non_neutral_genome_frac', 'subclo_genome_frac', 'clo_genome_frac')

non.codel.cox <- SCNAsBurdenSurAnalysis(gli.cn.alt.frac, features, 'Astrocytoma,IDHmut', c('G2', 'G3', 'G4'))
wt.cox <- SCNAsBurdenSurAnalysis(gli.cn.alt.frac, features, 'Glioblastoma,IDHwt', c('G4'))


lgg.non.codel.cox <- SCNAsBurdenSurAnalysis(gli.cn.alt.frac, features, 'Astrocytoma,IDHmut', c('G2', 'G3'))
gbm.non.codel.cox <- SCNAsBurdenSurAnalysis(gli.cn.alt.frac, features, 'Astrocytoma,IDHmut', c('G4'))


# save
setwd('/result/Section3')
SCNAs.burden.cox <- do.call(rbind, list(non.codel.cox, lgg.non.codel.cox, gbm.non.codel.cox, wt.cox))
write.table(SCNAs.burden.cox, file='SCNAs_burden_cox.txt', sep='\t', row.names=FALSE, col.names=TRUE, quote=FALSE)

