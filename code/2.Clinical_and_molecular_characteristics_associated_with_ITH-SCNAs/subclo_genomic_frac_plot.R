
# load data
gli.cn.alt.frac <- readRDS(file='/data/gli_glod_cn_alt_frac.rds')

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
codel.path <- '/result/Section2/Codel' 
CliMolAssoSubcloPlot(gli.cn.alt.frac, features, 'Oligodendroglioma,IDHmut-codel', codel.feacs, codel.path)


non.codel.feacs <- c("age_group", "gender", "histological_grade", "kps_group", "tumor_location", "first_presenting_symptom", 
 'CDKN2AB', 'TERT_EXPRESSION_STATUS', 'MGMT_PROMOTER_STATUS', 'SUPERVISED_DNA_METHYLATION_CLUSTER')
non.codel.path <- '/result/Section2/Noncodel'
CliMolAssoSubcloPlot(gli.cn.alt.frac, features, 'Astrocytoma,IDHmut', non.codel.feacs, non.codel.path)


wt.feacs <- c("age_group", 'CHR_7_GAIN_CHR_10_LOSS', 'CDKN2AB', 'EGFR', 'CHR_19_20_CO_GAIN',
 'TRANSCRIPTOME_SUBTYPE', 'SUPERVISED_DNA_METHYLATION_CLUSTER')
idh.wt.path <- '/result/Section2/IDHwt' 
CliMolAssoSubcloPlot(gli.cn.alt.frac, features, 'Glioblastoma,IDHwt', wt.feacs, idh.wt.path)



#############plot-age-subclo genome alteration
age.subclo.genome.frac.sca <- ggscatter(data=gli.cn.alt.frac, x='age', y='subclo_genome_frac', 
 legend='none', xlab='Age at diagnosis', ylab='Burden of subclonal SCNAs', add='reg.line', 
 conf.int=FALSE, cor.coef=TRUE, color = "Integrated_Diagnoses", shape = "Integrated_Diagnoses", 
 palette = c("#00AFBB", "#E7B800", "#FC4E07"), size = 0.8, alpha = 1.0, cor.method = 'spearman', 
 facet.by = "Integrated_Diagnoses") + theme(aspect.ratio = 1)

age.clo.genome.frac.sca <- ggscatter(data=gli.cn.alt.frac, x='age', y='clo_genome_frac', 
 legend='none', xlab='Age at diagnosis', ylab='Burden of clonal SCNAs', add='reg.line', 
 conf.int=FALSE, cor.coef=TRUE, color = "Integrated_Diagnoses", shape = "Integrated_Diagnoses", 
 palette = c("#00AFBB", "#E7B800", "#FC4E07"), size = 0.8, alpha = 1.0, cor.method = 'spearman', 
 facet.by = "Integrated_Diagnoses") + theme(aspect.ratio = 1)

source('/code/Function/multiplot.r')
pdf(file='/result/Section2/age_subclo_genome_alt_frac_corr.pdf')
 multiplot(plotlist=list(age.subclo.genome.frac.sca, age.clo.genome.frac.sca), cols=1)
dev.off()



# > aggregate(age ~ Integrated_Diagnoses, data=gli.cn.alt.frac, median)
            # Integrated_Diagnoses age
# 1             Astrocytoma,IDHmut  36
# 2             Glioblastoma,IDHwt  61
# 3 Oligodendroglioma,IDHmut-codel  44


#############plot-age(Discrete)-subclo genome alteration
gli.cn.alt.frac$age_median_group <- sapply(seq(nrow(gli.cn.alt.frac)), function(index){
 dat <- gli.cn.alt.frac[index, , FALSE]
 if(dat$Integrated_Diagnoses == 'Oligodendroglioma,IDHmut-codel'){
  label <- ifelse(dat$age >=44, '>=44', '<44')
 }else if(dat$Integrated_Diagnoses == 'Astrocytoma,IDHmut'){
  label <- ifelse(dat$age >=36, '>=36', '<36')
 }else{
  label <- ifelse(dat$age >=61, '>=61', '<61')
 }
 return(label)
})

age.subclo.genome.frac <- ggboxplot(data=gli.cn.alt.frac, x='age_median_group', y='subclo_genome_frac', legend='none', 
 xlab=FALSE, ylab='Burden of subclonal SCNAs', color = "Integrated_Diagnoses", add='jitter', shape = "Integrated_Diagnoses", 
 palette = c("#00AFBB", "#E7B800", "#FC4E07"), size = 0.8, alpha = 1.0, facet.by = "Integrated_Diagnoses") + 
 theme(aspect.ratio = 1) + stat_compare_means(label.x.npc = 'left', label.y.npc = 'top')

age.clo.genome.frac <- ggboxplot(data=gli.cn.alt.frac, x='age_median_group', y='clo_genome_frac', legend='none', 
 xlab=FALSE, ylab='Burden of clonal SCNAs', color = "Integrated_Diagnoses", add='jitter', shape = "Integrated_Diagnoses", 
 palette = c("#00AFBB", "#E7B800", "#FC4E07"), size = 0.8, alpha = 1.0, facet.by = "Integrated_Diagnoses") + 
 theme(aspect.ratio = 1) + stat_compare_means(label.x.npc = 'left', label.y.npc = 'top')

pdf(file='/result/Section2/age_subclo_genome_alt_frac.pdf')
 multiplot(plotlist=list(age.subclo.genome.frac, age.clo.genome.frac), cols=1)
dev.off()




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

pdf(file='/result/Section2/grade_subclo_genome_alt_frac.pdf')
 multiplot(plotlist=list(grade.subclo.genome.frac, grade.clo.genome.frac), cols=1)
dev.off()




