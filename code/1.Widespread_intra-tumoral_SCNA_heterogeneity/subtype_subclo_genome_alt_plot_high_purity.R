
# load data
gli.cn.alt.frac <- readRDS(file='/data/gli_glod_cn_alt_frac.rds')

# ABSOLUTE-based tumour purity
abs.puri.ploi <- readRDS('/data/tcga_gli_puri_ploi.rds')
abs.puri.ploi <- subset(abs.puri.ploi, purity>=0.65)

gli.cn.alt.frac <- gli.cn.alt.frac[rownames(abs.puri.ploi), ]


# plot
SubtypeSubcloBoxPlot <- function(dat, alt.feac, sub.feac, label, cols){
 
 dat <- dat[!is.na(dat[, sub.feac]), ]
 dat[, grep(sub.feac, colnames(dat))] <- factor(dat[, sub.feac], levels = label)
 
 library('ggpubr')
 
 shape.s <- c(21, 24, 22, 25, 23, 3, 4)[as.numeric(dat[, sub.feac])]
 color.s <- cols[as.numeric(dat[, sub.feac])]
 
 box.plot <- ggboxplot(data = dat, x = sub.feac, y=alt.feac, xlab = FALSE, ylab = alt.feac, 
 legend = 'none', alpha = 1.0, add = 'jitter',
 add.params=list(shape=shape.s, color = color.s, size = 0.8), outlier.shape = NA) + theme(aspect.ratio = 1)

 comparison <- combn(x = label, m = 2, simplify = FALSE)
 box.plot <- box.plot + stat_compare_means(comparisons = comparison) + 
 stat_compare_means(label.x.npc = 'left', label.y.npc = 'top')
 
 return(box.plot)
}


source('/code/Function/multiplot.r')

FeacturePlotBySubtype <- function(dat, features, subtype, labels, cols, out.path, file.name){
 
 pdf(paste0(out.path, file.name))
 plot.list <- lapply(features, function(feature){
  p <- SubtypeSubcloBoxPlot(dat, feature, subtype, labels, cols)
  return(p)
 })

 multiplot(plotlist=plot.list, layout=matrix(seq(1, 4), ncol=2, byrow = TRUE))
 dev.off()
}

out.path <- '/result/Section1/'
features <- c('non_neutral_genome_frac', 'subclo_genome_frac', 'clo_genome_frac', 'subclo_cn_alt_frac')


FeacturePlotBySubtype(dat=gli.cn.alt.frac, 
 features, subtype='IDH_CODEL_SUBTYPE', labels=c('IDHmut-codel', 'IDHmut-non-codel', 'IDHwt'), 
 cols=c("#FC4E07", "#E7B800", "#00AFBB"), out.path, file.name='high_purity_mole_subtype_subclo_alt_frac.pdf')

FeacturePlotBySubtype(dat=subset(gli.cn.alt.frac, cancer_type=='LGG'), 
 features, subtype='IDH_CODEL_SUBTYPE', labels=c('IDHmut-codel', 'IDHmut-non-codel', 'IDHwt'), 
 cols=c("#FC4E07", "#E7B800", "#00AFBB"), out.path, file.name='high_purity_lgg_mole_subtype_subclo_alt_frac.pdf')

FeacturePlotBySubtype(dat=subset(gli.cn.alt.frac, cancer_type=='GBM'), 
 features, subtype='IDH_CODEL_SUBTYPE', labels=c('IDHmut-codel', 'IDHmut-non-codel', 'IDHwt'), 
 cols=c("#FC4E07", "#E7B800", "#00AFBB"), out.path, file.name='high_purity_gbm_mole_subtype_subclo_alt_frac.pdf')
