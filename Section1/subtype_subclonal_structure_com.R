
# load data
gli.cn.alt.frac <- readRDS(file='/pub5/xiaoyun/Jobs/J22/CopyNumberClonalityProject/Resource/CuratedData/gli_glod_cn_alt_frac.rds')

# PhyloWGS subclonal structure
setwd('/pub5/xiaoyun/Jobs/J22/CopyNumberClonalityProject/Results/Section1/Resources/Franck_PlosGenetics_2018_TCGAITH/')
phylowgs.clo <- read.csv(file='TCGAPancancerITH.txt', header=TRUE, sep='\t', stringsAsFactors=FALSE)
rownames(phylowgs.clo) <- phylowgs.clo$sample_name

gli.cn.alt.frac <- merge(gli.cn.alt.frac, phylowgs.clo[, c('number.of.clones', 'Tree.score')], by='row.names', all.x=TRUE)
rownames(gli.cn.alt.frac) <- gli.cn.alt.frac$Row.names
gli.cn.alt.frac <- gli.cn.alt.frac[, -1]

# EXPANDS subclonal structure
library('readxl')
setwd('/pub5/xiaoyun/Jobs/J22/CopyNumberClonalityProject/Results/Section1/Resources/Andor_NatureMedicine_2016_TCGAITH')
expands.clo <- as.data.frame(read_xlsx(path='TCGAPancancerITH.xlsx', sheet=1, na="NA"))
rownames(expands.clo) <- paste0(expands.clo$"TCGA ID", '-01')

gli.cn.alt.frac <- merge(gli.cn.alt.frac, expands.clo[, 'CloneNumber(PurityNormalized)', FALSE], by='row.names', all.x=TRUE)
rownames(gli.cn.alt.frac) <- gli.cn.alt.frac$Row.names
gli.cn.alt.frac <- gli.cn.alt.frac[, -1]


# Pyclone subclonal structure
source('/pub5/xiaoyun/Jobs/J22/CopyNumberClonalityProject/Results/Section1/RScripts/PycloneFunctionModule.R')
# The number of detectable subclonal populations
file.path <- '/pub5/xiaoyun/Jobs/J22/EvoClass2.0/Section4/Results/PycloneITH/TCGA_LGG/'
lgg.pyclone.clo <- PyCloneClonalPopulationNumber(file.path)
file.path <- '/pub5/xiaoyun/Jobs/J22/EvoClass2.0/Section4/Results/PycloneITH/TCGA_GBM/'
gbm.pyclone.clo <- PyCloneClonalPopulationNumber(file.path)

pyclone.clo <- data.frame(Pyclone=c(lgg.pyclone.clo, gbm.pyclone.clo))
gli.cn.alt.frac <- merge(gli.cn.alt.frac, pyclone.clo, by='row.names', all.x=TRUE)



# plot
source('/pub5/xiaoyun/Jobs/J22/RScripts/LSY.RScripts/PlotFunction/multiplot.r')

ClonalPopulationBoxplot <- function(dat, subtype, feac, labels, cols){
 
 dat <- dat[!is.na(dat[, feac]), ]
 dat[, grep(subtype, colnames(dat))] <- factor(dat[, subtype], levels = labels)
 
 shape.s <- c(21, 24, 22, 25, 23, 3, 4)[as.numeric(dat[, subtype])]
 color.s <- cols[as.numeric(dat[, subtype])]

 
 box.plot <- ggboxplot(data = dat, x = subtype, y=feac, 
 xlab = FALSE, ylab = feac, legend = 'none', alpha = 1.0, add = 'jitter', 
 add.params=list(shape=shape.s, color = color.s, size = 0.8), outlier.shape = NA) + theme(aspect.ratio = 1)

 comparison <- combn(x = labels, m = 2, simplify = FALSE)
 box.plot <- box.plot + stat_compare_means(comparisons = comparison) + 
 stat_compare_means(label.x.npc = 'left', label.y.npc = 'top')
 
 return(box.plot)
}


FeacturePlotBySubtype <- function(dat, features, subtype, labels, cols, file.name){
 
 pdf(paste0('/pub5/xiaoyun/Jobs/J22/CopyNumberClonalityProject/Results/Section1/Results/ClonalPopulationCompare/', file.name))
 plot.list <- lapply(features, function(feature){
  p <- ClonalPopulationBoxplot(dat, subtype, feature, labels, cols)
  return(p)
 })

 multiplot(plotlist=plot.list, layout=matrix(seq(1, 4), ncol=2, byrow = TRUE))
 dev.off()
}

features <- c('number.of.clones', 'Tree.score', 'Pyclone', 'CloneNumber(PurityNormalized)')


FeacturePlotBySubtype(dat=gli.cn.alt.frac, 
 features, subtype='IDH_CODEL_SUBTYPE', labels=c('IDHwt', 'IDHmut-non-codel', 'IDHmut-codel'), 
 cols=c("#00AFBB", "#E7B800", "#FC4E07"), file.name='mole_subtype_subclo_popul.pdf')

FeacturePlotBySubtype(dat=subset(gli.cn.alt.frac, cancer_type=='LGG'), 
 features, subtype='IDH_CODEL_SUBTYPE', labels=c('IDHwt', 'IDHmut-non-codel', 'IDHmut-codel'), 
 cols=c("#00AFBB", "#E7B800", "#FC4E07"), file.name='lgg_mole_subtype_subclo_popul.pdf')

FeacturePlotBySubtype(dat=subset(gli.cn.alt.frac, cancer_type=='GBM'), 
 features, subtype='IDH_CODEL_SUBTYPE', labels=c('IDHwt', 'IDHmut-non-codel', 'IDHmut-codel'), 
 cols=c("#00AFBB", "#E7B800", "#FC4E07"), file.name='gbm_mole_subtype_subclo_popul.pdf')


FeacturePlotBySubtype(dat=gli.cn.alt.frac, 
 features, subtype='molecular_histological_type', labels=c('Glioblastoma', 'Astrocytoma', 'Oligodendroglioma'), 
 cols=c("#377EB8", "#4DAF4A", "#984EA3"), file.name='hist_subtype_subclo_popul.pdf')

FeacturePlotBySubtype(dat=subset(gli.cn.alt.frac, cancer_type=='LGG'), 
 features, subtype='molecular_histological_type', labels=c('Glioblastoma', 'Astrocytoma', 'Oligodendroglioma'), 
 cols=c("#377EB8", "#4DAF4A", "#984EA3"), file.name='lgg_hist_subtype_subclo_popul.pdf')

FeacturePlotBySubtype(dat=subset(gli.cn.alt.frac, cancer_type=='GBM'), 
 features, subtype='molecular_histological_type', labels=c('Glioblastoma', 'Astrocytoma', 'Oligodendroglioma'), 
 cols=c("#377EB8", "#4DAF4A", "#984EA3"), file.name='gbm_hist_subtype_subclo_popul.pdf')


