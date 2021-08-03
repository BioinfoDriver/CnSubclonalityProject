
# load data
gli.cn.alt.frac <- readRDS(file='/data/gli_glod_cn_alt_frac.rds')

# ABSOLUTE-based tumour purity
abs.puri.ploi <- readRDS('/data/tcga_gli_puri_ploi.rds')
abs.puri.ploi <- subset(abs.puri.ploi, purity>=0.65)

gli.cn.alt.frac <- gli.cn.alt.frac[rownames(abs.puri.ploi), ]


# Cox
source('/code/Function/Cox.function.R')

SCNAsBurdenSurAnalysis <- function(sur.dat, feacs, subtype, grade, disperse=FALSE){
 
 sur.dat <- subset(sur.dat, IDH_CODEL_SUBTYPE %in% subtype & histological_grade %in% grade)
 
 cox.pvalues <- lapply(feacs, function(feac){
  
  cox.pvalue <- lapply(c('os', 'dss', 'pfi'), function(x){
   
   if(length(grade) > 1){
    dat <- sur.dat[, c('bcr_patient_barcode', x, paste0(x, '_time'), feac, 'age', 'gender', 'histological_grade')]
   }else{
    dat <- sur.dat[, c('bcr_patient_barcode', x, paste0(x, '_time'), feac, 'age', 'gender')]
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

non.codel.cox <- SCNAsBurdenSurAnalysis(gli.cn.alt.frac, features, 'IDHmut-non-codel', c('G2', 'G3', 'G4'))
wt.cox <- SCNAsBurdenSurAnalysis(gli.cn.alt.frac, features, 'IDHwt', c('G2', 'G3', 'G4'))


lgg.non.codel.cox <- SCNAsBurdenSurAnalysis(gli.cn.alt.frac, features, 'IDHmut-non-codel', c('G2', 'G3'))
lgg.wt.cox <- SCNAsBurdenSurAnalysis(gli.cn.alt.frac, features, 'IDHwt', c('G2', 'G3'))


gbm.non.codel.cox <- SCNAsBurdenSurAnalysis(gli.cn.alt.frac, features, 'IDHmut-non-codel', c('G4'))
gbm.wt.cox <- SCNAsBurdenSurAnalysis(gli.cn.alt.frac, features, 'IDHwt', c('G4'))


# save
setwd('/result/Section3/')

SCNAs.burden.cox <- do.call(rbind, list(non.codel.cox, lgg.non.codel.cox, 
 gbm.non.codel.cox, wt.cox, lgg.wt.cox, gbm.wt.cox))

write.table(SCNAs.burden.cox, file='high_purity_SCNAs_burden_cox.txt', sep='\t', row.names=FALSE, col.names=TRUE, quote=FALSE)


