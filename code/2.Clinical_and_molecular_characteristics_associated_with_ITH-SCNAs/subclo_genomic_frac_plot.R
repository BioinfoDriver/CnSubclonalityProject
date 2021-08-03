
# load data
gli.cn.alt.frac <- readRDS(file='/pub5/xiaoyun/Jobs/J22/CopyNumberClonalityProject/Resource/CuratedData/gli_glod_cn_alt_frac.rds')

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


CliMolAssoSubcloPlot <- function(dat, features, subtype, cancer.type, cli.mol.feacs, out.path){
 dat <- subset(dat, IDH_CODEL_SUBTYPE %in% subtype & cancer_type %in% cancer.type)
 
 for(cli.mol.feac in cli.mol.feacs){
   
  num <- table(dat[, cli.mol.feac])
  char <- names(num)[num>5]
  dat.del <- subset(dat, get(cli.mol.feac) %in% char)

  plot.list <- lapply(features, function(feac){
   
   box.plot <- ggboxplot(data = dat.del, x = cli.mol.feac, y = feac, color = cli.mol.feac, xlab = FALSE, 
   ylab = feac, legend = 'none', add = "jitter", shape = cli.mol.feac, outlier.shape = NA,
   alpha = 1.0, add.params=list(size = 0.8)) + theme(aspect.ratio = 1)
   
   
   # comparison <- combn(x = na.omit(as.character(unique(dat.del[, cli.mol.feac]))), m = 2, simplify = FALSE)
   # box.plot <- box.plot + stat_compare_means(comparisons = comparison) + 
   # stat_compare_means(label.x.npc = 'left', label.y.npc = 'top')

 
   return(box.plot)
  })
  source('/pub5/xiaoyun/Jobs/J22/RScripts/LSY.RScripts/PlotFunction/multiplot.r')
  
  pdf(file.path(out.path, paste0(cli.mol.feac, '.pdf')))
   multiplot(plotlist=plot.list, layout=matrix(seq(1, 4), ncol=2, byrow = TRUE))
  dev.off()
  
 }
}

features <- c('subclo_genome_frac', 'clo_genome_frac')

codel.feacs <- c("age_group", "histological_grade", "laterality", 'TERT_EXPRESSION_STATUS')
codel.path <- '/pub5/xiaoyun/Jobs/J22/CopyNumberClonalityProject/Results/Section2/Results/Codel' 
CliMolAssoSubcloPlot(gli.cn.alt.frac, features, 'IDHmut-codel', c('LGG', 'GBM'), codel.feacs, codel.path)


non.codel.feacs <- c("age_group", "gender", "histological_grade", "molecular_histological_type", "kps_group", "tumor_location", 
 "first_presenting_symptom", 'MGMT_PROMOTER_STATUS', 'TERT_EXPRESSION_STATUS', 'SUPERVISED_DNA_METHYLATION_CLUSTER')
non.codel.path <- '/pub5/xiaoyun/Jobs/J22/CopyNumberClonalityProject/Results/Section2/Results/Noncodel'
CliMolAssoSubcloPlot(gli.cn.alt.frac, features, 'IDHmut-non-codel', c('LGG', 'GBM'), non.codel.feacs, non.codel.path)


wt.feacs <- c("age_group", "histological_grade", "molecular_histological_type", 'CHR_7_GAIN_CHR_10_LOSS', 
 'TELOMERE_MAINTENANCE', 'TRANSCRIPTOME_SUBTYPE', 'SUPERVISED_DNA_METHYLATION_CLUSTER')
idh.wt.path <- '/pub5/xiaoyun/Jobs/J22/CopyNumberClonalityProject/Results/Section2/Results/IDHwt' 
CliMolAssoSubcloPlot(gli.cn.alt.frac, features, 'IDHwt', c('LGG', 'GBM'), wt.feacs, idh.wt.path)



#############plot-age-subclo genome alteration
age.subclo.genome.frac.sca <- ggscatter(data=gli.cn.alt.frac, x='age', y='subclo_genome_frac', 
 legend='none', xlab='Age at diagnosis', ylab='Burden of subclonal SCNAs', add='reg.line', 
 conf.int=FALSE, cor.coef=TRUE, color = "IDH_CODEL_SUBTYPE", shape = "IDH_CODEL_SUBTYPE", 
 palette = c("#00AFBB", "#E7B800", "#FC4E07"), size = 0.8, alpha = 1.0, cor.method = 'spearman', 
 facet.by = "IDH_CODEL_SUBTYPE") + theme(aspect.ratio = 1)

age.clo.genome.frac.sca <- ggscatter(data=gli.cn.alt.frac, x='age', y='clo_genome_frac', 
 legend='none', xlab='Age at diagnosis', ylab='Burden of clonal SCNAs', add='reg.line', 
 conf.int=FALSE, cor.coef=TRUE, color = "IDH_CODEL_SUBTYPE", shape = "IDH_CODEL_SUBTYPE", 
 palette = c("#00AFBB", "#E7B800", "#FC4E07"), size = 0.8, alpha = 1.0, cor.method = 'spearman', 
 facet.by = "IDH_CODEL_SUBTYPE") + theme(aspect.ratio = 1)

source('/pub5/xiaoyun/Jobs/J22/RScripts/LSY.RScripts/PlotFunction/multiplot.r')
pdf(file='/pub5/xiaoyun/Jobs/J22/CopyNumberClonalityProject/Results/Section2/Results/age_subclo_genome_alt_frac_corr.pdf')
 multiplot(plotlist=list(age.subclo.genome.frac.sca, age.clo.genome.frac.sca), cols=1)
dev.off()



# > aggregate(age ~ IDH_CODEL_SUBTYPE, data=gli.cn.alt.frac, median)
  # IDH_CODEL_SUBTYPE age
# 1      IDHmut-codel  44
# 2  IDHmut-non-codel  36
# 3             IDHwt  61

#############plot-age(Discrete)-subclo genome alteration
gli.cn.alt.frac$age_median_group <- sapply(seq(nrow(gli.cn.alt.frac)), function(index){
 dat <- gli.cn.alt.frac[index, , FALSE]
 if(dat$IDH_CODEL_SUBTYPE == 'IDHmut-codel'){
  label <- ifelse(dat$age >=44, '>=44', '<44')
 }else if(dat$IDH_CODEL_SUBTYPE == 'IDHmut-non-codel'){
  label <- ifelse(dat$age >=36, '>=36', '<36')
 }else{
  label <- ifelse(dat$age >=61, '>=61', '<61')
 }
 return(label)
})

age.subclo.genome.frac <- ggboxplot(data=gli.cn.alt.frac, x='age_median_group', y='subclo_genome_frac', legend='none', 
 xlab=FALSE, ylab='Burden of subclonal SCNAs', color = "IDH_CODEL_SUBTYPE", add='jitter', shape = "IDH_CODEL_SUBTYPE", 
 palette = c("#00AFBB", "#E7B800", "#FC4E07"), size = 0.8, alpha = 1.0, facet.by = "IDH_CODEL_SUBTYPE") + 
 theme(aspect.ratio = 1) + stat_compare_means(label.x.npc = 'left', label.y.npc = 'top')

age.clo.genome.frac <- ggboxplot(data=gli.cn.alt.frac, x='age_median_group', y='clo_genome_frac', legend='none', 
 xlab=FALSE, ylab='Burden of clonal SCNAs', color = "IDH_CODEL_SUBTYPE", add='jitter', shape = "IDH_CODEL_SUBTYPE", 
 palette = c("#00AFBB", "#E7B800", "#FC4E07"), size = 0.8, alpha = 1.0, facet.by = "IDH_CODEL_SUBTYPE") + 
 theme(aspect.ratio = 1) + stat_compare_means(label.x.npc = 'left', label.y.npc = 'top')

source('/pub5/xiaoyun/Jobs/J22/RScripts/LSY.RScripts/PlotFunction/multiplot.r')
pdf(file='/pub5/xiaoyun/Jobs/J22/CopyNumberClonalityProject/Results/Section2/Results/age_subclo_genome_alt_frac.pdf')
 multiplot(plotlist=list(age.subclo.genome.frac, age.clo.genome.frac), cols=1)
dev.off()




#############plot-grade-subclo genome alteration
grade.subclo.genome.frac <- ggboxplot(data=gli.cn.alt.frac, x='histological_grade', y='subclo_genome_frac', 
 legend='none', xlab=FALSE, ylab='Burden of subclonal SCNAs', color = "IDH_CODEL_SUBTYPE", add='jitter', 
 shape = "IDH_CODEL_SUBTYPE", palette = c("#00AFBB", "#E7B800", "#FC4E07"), size = 0.8, alpha = 1.0, 
 facet.by = "IDH_CODEL_SUBTYPE", select=c('G2', 'G3', 'G4')) + theme(aspect.ratio = 1) + 
 stat_compare_means(comparisons = combn(x = c('G2', 'G3', 'G4'), m = 2, simplify = FALSE)) + 
 stat_compare_means(label.x.npc = 'left', label.y.npc = 'top')

grade.clo.genome.frac <- ggboxplot(data=gli.cn.alt.frac, x='histological_grade', y='clo_genome_frac', 
 legend='none', xlab=FALSE, ylab='Burden of clonal SCNAs', color = "IDH_CODEL_SUBTYPE", add='jitter', 
 shape = "IDH_CODEL_SUBTYPE", palette = c("#00AFBB", "#E7B800", "#FC4E07"), size = 0.8, alpha = 1.0, 
 facet.by = "IDH_CODEL_SUBTYPE", select=c('G2', 'G3', 'G4')) + theme(aspect.ratio = 1) + 
 stat_compare_means(comparisons = combn(x = c('G2', 'G3', 'G4'), m = 2, simplify = FALSE)) + 
 stat_compare_means(label.x.npc = 'left', label.y.npc = 'top')

source('/pub5/xiaoyun/Jobs/J22/RScripts/LSY.RScripts/PlotFunction/multiplot.r')
pdf(file='/pub5/xiaoyun/Jobs/J22/CopyNumberClonalityProject/Results/Section2/Results/grade_subclo_genome_alt_frac.pdf')
 multiplot(plotlist=list(grade.subclo.genome.frac, grade.clo.genome.frac), cols=1)
dev.off()




