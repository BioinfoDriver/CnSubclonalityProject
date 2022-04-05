
#####################Survival
rtks <- read.csv(file='/data/RTKs.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)

gold.set <- readRDS('/data/gold_set.rds')
gene.het <- readRDS(file='/data/gene.het.rds')
gene.het <- gene.het[, paste(gold.set, '01', sep = '-')]

tcga.cli.data <- readRDS('/data/tcga_glioma_cli_mol.rds')
tcga.cli.data <- tcga.cli.data[gold.set, ]

gli.cn.alt.frac <- readRDS(file='/data/tcga_gli_cn_alt_frac.rds') 
gli.cn.alt.frac <- gli.cn.alt.frac[, c('non_neutral_genome_frac', 'subclo_genome_frac', 'clo_genome_frac')]
colnames(gli.cn.alt.frac) <- c('SCNA_burden', 'clo_SCNA_burden', 'subclo_SCNA_burden')


# TP53 mutation
load('/data/PureGliomaMutData.RData')
gli.sig.mut <- subset(pure.glioma.mut.data, Patient %in% paste0(gold.set, '-01')) # 760
gli.sig.mut <- subset(gli.sig.mut, Hugo_Symbol %in% 'TP53')
gli.sig.mut <- subset(gli.sig.mut, !(Variant_Classification %in% c("3'UTR", "5'UTR", "Intron", "Silent")))
tp53.mut <- gold.set %in% substr(gli.sig.mut$Patient, 1, 12)
names(tp53.mut) <- paste0(gold.set, '-01')

gli.cn.alt.frac$TP53 <- tp53.mut[rownames(gli.cn.alt.frac)]

# library('tidyverse')
tcga.cli.data <- tcga.cli.data %>% mutate(
	kps_score = as.character(cut_width(karnofsky_performance_score, width = 10, boundary = 10)),
	kps_score = ifelse(kps_score %in% c('[20,30]', '(30,40]', '(40,50]', '(50,60]', '(60,70]', '(70,80]'), "<=80", kps_score)) 

rownames(tcga.cli.data) <- paste0(tcga.cli.data$bcr_patient_barcode, '-01', sep='')

tcga.cli.data <- cbind(tcga.cli.data, gli.cn.alt.frac[rownames(tcga.cli.data), ])


# Cox
source('/code/Function/Cox.function.R')
RtkMulCoxSurAnalysis <- function(het.mat, rtk, group.index, subtype, grade, cli.data, up.time){
 
 rtk.het <- het.mat[rtk, ]

 rtk.clo.sam <- colnames(rtk.het[, colSums(rtk.het == 1) >= 1])
 rtk.sub.sam <- colnames(rtk.het[, colSums(rtk.het == 2) >= 1])
 
 if(group.index=='subclonaldominant'){
  rtk.clo.sam <- setdiff(rtk.clo.sam, rtk.sub.sam)

 }else if(group.index=='clonaldominant'){
  rtk.sub.sam <- setdiff(rtk.sub.sam, rtk.clo.sam)
 
 }else{
  com.sam <- intersect(rtk.clo.sam, rtk.sub.sam)
  rtk.clo.sam <- setdiff(rtk.clo.sam, com.sam)
  rtk.sub.sam <- setdiff(rtk.sub.sam, com.sam)
  
 }
 
 # non.rtk.alt.sam <- setdiff(colnames(het.mat), c(rtk.clo.sam, rtk.sub.sam))
 
 # rtk.clo.mat <- data.frame(group = rep(c('rtk clonal amp', 'rtk subclonal amp', 'rtk wt'), 
  # times=c(length(rtk.clo.sam), length(rtk.sub.sam), length(non.rtk.alt.sam))))
 # rownames(rtk.clo.mat) <- c(rtk.clo.sam, rtk.sub.sam, non.rtk.alt.sam)

 rtk.clo.mat <- data.frame(group = rep(c('rtk clonal amp', 'rtk subclonal amp'), 
  times=c(length(rtk.clo.sam), length(rtk.sub.sam))))
 rownames(rtk.clo.mat) <- c(rtk.clo.sam, rtk.sub.sam)

 cli.data <- subset(cli.data, Integrated_Diagnoses %in% subtype & histological_grade %in% grade)
 cli.data <- merge(cli.data, rtk.clo.mat, by = 'row.names')
 
 # return(cli.data)
 # OS, DSS PFI
 mul.cox.res <- lapply(c('os', 'dss', 'pfi'), function(surv.type){
  
  surv.data <- cli.data[, c('Row.names', surv.type, paste0(surv.type, '_time'), 
   'age', 'gender', 'histological_grade', 'Integrated_Diagnoses', 'group', 'kps_score', 'TP53',
	'SCNA_burden', 'clo_SCNA_burden')]
  
  colnames(surv.data) <- c('Patient_ID', 'event', 'time', 'age', 'gender', 'grade', 
   'subtype', 'group', 'kps_score', 'TP53', 'SCNA_burden', 'clo_SCNA')
  
  print(nrow(surv.data))
  surv.data <- subset(surv.data, !is.na(grade))
  if (!is.null(up.time)) surv.data <- surv.data[surv.data$time <= up.time,]
  
  
  print(dim(surv.data))
  if(length(unique(surv.data$subtype)) > 1 & length(unique(surv.data$grade)) > 1){
   cox.res <- Cox.function(surv.data$time, surv.data$event, clinical.data=surv.data, clinical.variate=c(4:8, 10:12))
  }
  
  if(length(unique(surv.data$subtype)) == 1 & length(unique(surv.data$grade)) > 1){
   cox.res <- Cox.function(surv.data$time, surv.data$event, clinical.data=surv.data, clinical.variate=c(4, 5, 6, 8, 10:12))
  }
  
  if(length(unique(surv.data$subtype)) > 1 & length(unique(surv.data$grade)) == 1){
   cox.res <- Cox.function(surv.data$time, surv.data$event, clinical.data=surv.data, clinical.variate=c(4, 5, 7, 8, 10:12))
  }
  
  if(length(unique(surv.data$subtype)) == 1 & length(unique(surv.data$grade)) == 1){
   cox.res <- Cox.function(surv.data$time, surv.data$event, clinical.data=surv.data, clinical.variate=c(4, 5, 8:12))
  }  

  return(cox.res)
 })

 return(mul.cox.res)
}


RtkMulCoxSurAnalysis(gene.het, rtks$Approved.symbol, 'exclude', # 
 c('Glioblastoma,IDHwt'), c('G4'), subset(tcga.cli.data, !is.na(kps_score)), 1096)
