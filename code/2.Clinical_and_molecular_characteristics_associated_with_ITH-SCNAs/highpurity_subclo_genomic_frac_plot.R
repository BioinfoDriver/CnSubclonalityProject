# load data
gli.cn.alt.frac <- readRDS(file='/data/gli_glod_cn_alt_frac.rds')

# ABSOLUTE-based tumour purity
abs.puri.ploi <- readRDS('/data/tcga_gli_puri_ploi.rds')
abs.puri.ploi <- subset(abs.puri.ploi, purity>=0.65)

gli.cn.alt.frac <- gli.cn.alt.frac[intersect(rownames(abs.puri.ploi), rownames(gli.cn.alt.frac)), ]


library(dplyr)
library(ggpubr)
library(ggplot2)
gli.cn.alt.frac <- gli.cn.alt.frac %>%
  mutate(age_group = as.character(cut_width(age, width = 10, boundary = 10)), 
         age_group = ifelse(age_group %in% c("[10,20]", "(20,30]"), "<=30", age_group), 
         age_group = ifelse(age_group %in% c("(70,80]", "(80,90]"), ">70", age_group), 
		 kps_group = as.character(cut_width(karnofsky_performance_score, width = 10, boundary = 10)),
		 kps_group = ifelse(kps_group %in% c('[20,30]', '(30,40]', '(40,50]', '(50,60]', '(60,70]'), "<=70", kps_group))

gli.cn.alt.frac$age_group <- factor(gli.cn.alt.frac$age_group, 
 levels=c('<=30', '(30,40]', '(40,50]', '(50,60]', '(60,70]', '>70'))
gli.cn.alt.frac$histological_grade <- factor(gli.cn.alt.frac$histological_grade, levels=c('G2', 'G3', 'G4'))
gli.cn.alt.frac$kps_group <- factor(gli.cn.alt.frac$kps_group, 
 levels=c('<=70', '(70,80]', '(80,90]', '(90,100]'))

gli.cn.alt.frac$CDKN2AB <- ifelse(gli.cn.alt.frac$CDKN2AB == 1, 'Deletion', 'Diploid')
gli.cn.alt.frac$EGFR <- ifelse(gli.cn.alt.frac$EGFR == 1, 'Amplification', 'Diploid')


CliMolAssoSubcloPlot <- function(dat, features, subtype, cli.mol.feacs, out.path){
 dat <- subset(dat, Integrated_Diagnoses %in% subtype )
 
 for(cli.mol.feac in cli.mol.feacs){
   
  num <- table(dat[, cli.mol.feac])
  char <- names(num)[num>5]
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

features <- c('subclo_genome_frac', 'clo_genome_frac')

codel.feacs <- c("age_group", "histological_grade", "laterality", 'TERT_EXPRESSION_STATUS')
codel.path <- '/result/Section2/Codel/HighPurity' 
CliMolAssoSubcloPlot(gli.cn.alt.frac, features, 'Oligodendroglioma,IDHmut-codel', codel.feacs, codel.path)


non.codel.feacs <- c("age_group", "gender", "histological_grade", "kps_group", "tumor_location", "first_presenting_symptom", 
 'CDKN2AB', 'TERT_EXPRESSION_STATUS', 'MGMT_PROMOTER_STATUS', 'SUPERVISED_DNA_METHYLATION_CLUSTER')
non.codel.path <- '/result/Section2/Noncodel/HighPurity'
CliMolAssoSubcloPlot(gli.cn.alt.frac, features, 'Astrocytoma,IDHmut', non.codel.feacs, non.codel.path)


wt.feacs <- c("age_group", 'CHR_7_GAIN_CHR_10_LOSS', 'CDKN2AB', 'EGFR', 'CHR_19_20_CO_GAIN',
 'TRANSCRIPTOME_SUBTYPE', 'SUPERVISED_DNA_METHYLATION_CLUSTER')
idh.wt.path <- '/result/Section2/IDHwt/HighPurity' 
CliMolAssoSubcloPlot(gli.cn.alt.frac, features, 'Glioblastoma,IDHwt', wt.feacs, idh.wt.path)




#############plot-grade-subclo genome alteration
grade.subclo.genome.frac <- ggboxplot(data=gli.cn.alt.frac, x='histological_grade', y='subclo_genome_frac', 
 legend='none', xlab=FALSE, ylab='Burden of subclonal SCNAs', color = "Integrated_Diagnoses", add='jitter', 
 shape = "Integrated_Diagnoses", palette = c("#00AFBB", "#E7B800", "#FC4E07"), size = 0.8, alpha = 1.0, 
 facet.by = "Integrated_Diagnoses", select=c('G2', 'G3', 'G4')) + theme(aspect.ratio = 1) + 
 stat_compare_means(comparisons = combn(x = c('G2', 'G3', 'G4'), m = 2, simplify = FALSE)) + 
 stat_compare_means(label.x.npc = 'left', label.y.npc = 'top')

grade.clo.genome.frac <- ggboxplot(data=gli.cn.alt.frac, x='histological_grade', y='clo_genome_frac', 
 legend='none', xlab=FALSE, ylab='Burden of clonal SCNAs', color = "Integrated_Diagnoses", add='jitter', 
 shape = "Integrated_Diagnoses", palette = c("#00AFBB", "#E7B800", "#FC4E07"), size = 0.8, alpha = 1.0, 
 facet.by = "Integrated_Diagnoses", select=c('G2', 'G3', 'G4')) + theme(aspect.ratio = 1) + 
 stat_compare_means(comparisons = combn(x = c('G2', 'G3', 'G4'), m = 2, simplify = FALSE)) + 
 stat_compare_means(label.x.npc = 'left', label.y.npc = 'top')

source('/code/Function/multiplot.r')
pdf(file='/result/Section2/highpurity_grade_subclo_genome_alt_frac.pdf')
 multiplot(plotlist=list(grade.subclo.genome.frac, grade.clo.genome.frac), cols=1)
dev.off()

