
glass.cna.ith <- readRDS('/data/glass_seqz_cna_ith.rds')

setwd('/data/OriginalData/syn17038081')
clinical.data <- read.csv(file='clinical_cases.tsv', header=TRUE, sep='\t', stringsAsFactors=FALSE)

glass.cna.ith <- merge(clinical.data, glass.cna.ith, by='case_barcode')
glass.cna.ith <- subset(glass.cna.ith, !is.na(Integrated_Diagnoses))
glass.cna.ith[glass.cna.ith=='NULL'] <- NA

glass.cna.ith$surgery_location[!( 
 glass.cna.ith$surgery_location %in% c('Temporal lobe', 'Parietal lobe', 'Frontal lobe'))] <- NA


# aggregate(case_age_diagnosis_years ~ Integrated_Diagnoses, data=glass.cna.ith, median)
glass.cna.ith$age_group <- sapply(seq(nrow(glass.cna.ith)), function(index){
 dat <- glass.cna.ith[index, , FALSE]
 if(dat$Integrated_Diagnoses == 'Oligodendroglioma,IDHmut-codel'){
  label <- ifelse(dat$case_age_diagnosis_years >=38, '>=38', '<38')
 }else if(dat$Integrated_Diagnoses == 'Astrocytoma,IDHmut'){
  label <- ifelse(dat$case_age_diagnosis_years >=33, '>=33', '<33')
 }else{
  label <- ifelse(dat$case_age_diagnosis_years >=56, '>=56', '<56')
 }
 return(label)
})


CliMolAssoSubcloPlot <- function(dat, features, subtype, cli.mol.feacs, out.path){
 library(ggpubr)
 library(ggplot2)
 dat <- subset(dat, Integrated_Diagnoses %in% subtype )
 
 for(cli.mol.feac in cli.mol.feacs){
   
  num <- table(dat[, cli.mol.feac])
  char <- names(num)[num>2]
  dat.del <- subset(dat, get(cli.mol.feac) %in% char)

  plot.list <- lapply(features, function(feac){
   
   box.plot <- ggboxplot(data = dat.del, x = cli.mol.feac, y = feac, color = cli.mol.feac, xlab = FALSE, 
   ylab = feac, legend = 'none', add = "jitter", shape = cli.mol.feac, outlier.shape = NA,
   alpha = 1.0, add.params=list(size = 0.8)) + theme(aspect.ratio = 1)
   
   
   comparison <- combn(x = na.omit(as.character(unique(dat.del[, cli.mol.feac]))), m = 2, simplify = FALSE)
   box.plot <- box.plot + stat_compare_means(comparisons = comparison) + 
   stat_compare_means(label.x.npc = 'left', label.y.npc = 'top')

 
   return(box.plot)
  })
  source('/code/Function/multiplot.r')
  
  pdf(file.path(out.path, paste0(cli.mol.feac, '.pdf')))
   multiplot(plotlist=plot.list, layout=matrix(seq(1, 4), ncol=2, byrow = TRUE))
  dev.off()
  
 }
}

features <- c('subclonal.cna.bur', 'clonal.cna.bur')
codel.feacs <- c("age_group",  "case_sex", "grade", "surgery_laterality", 'surgery_location', 'surgery_extent_of_resection')
codel.path <- '/result/Validation/cli_mol_asso/Codel' 
CliMolAssoSubcloPlot(glass.cna.ith, features, 'Oligodendroglioma,IDHmut-codel', codel.feacs, codel.path)


non.codel.feacs <- c("age_group", "case_sex", "grade", "surgery_laterality", 'surgery_location',
 'mgmt_methylation', 'surgery_extent_of_resection')
non.codel.path <- '/result/Validation/cli_mol_asso/Noncodel'
CliMolAssoSubcloPlot(glass.cna.ith, features, 'Astrocytoma,IDHmut', non.codel.feacs, non.codel.path)


wt.feacs <- c("age_group", "case_sex", 'mgmt_methylation', "surgery_laterality", 
 'surgery_location', 'surgery_extent_of_resection')
idh.wt.path <- '/result/Validation/cli_mol_asso/IDHwt' 
CliMolAssoSubcloPlot(glass.cna.ith, features, 'Glioblastoma,IDHwt', wt.feacs, idh.wt.path)
