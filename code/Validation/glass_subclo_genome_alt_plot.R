
# load data
glass.cna.ith <- readRDS('/data/glass_seqz_cna_ith.rds')

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


out.path <- '/result/Validation/'
t.features <- c('subclo.cna.perc.t', 'clonal.cna.bur.t', 'subclonal.cna.bur.t')

FeacturePlotBySubtype(dat=glass.cna.ith, t.features, 
 subtype='Integrated_Diagnoses', labels=c('Oligodendroglioma,IDHmut-codel', 'Astrocytoma,IDHmut', 'Glioblastoma,IDHwt'), 
 cols=c("#FC4E07", "#E7B800", "#00AFBB"), out.path, file.name='t_mole_subtype_subclo_alt_frac.pdf')


# kruskal.test(subclo.cna.perc.t ~ Integrated_Diagnoses, glass.cna.ith)$p.value
# 0.0001299936
# dplyr::bind_rows(dunn.test::dunn.test(glass.cna.ith$subclo.cna.perc.t, glass.cna.ith$Integrated_Diagnoses, method = 'bh')[2:5])
      # Z        P P.adjusted comparisons                                        
  # <dbl>    <dbl>      <dbl> <chr>                                              
# 1  3.52 0.000218   0.000327 Astrocytoma,IDHmut - Glioblastoma,IDHwt            
# 2  3.61 0.000150   0.000451 Astrocytoma,IDHmut - Oligodendroglioma,IDHmut-codel
# 3  1.44 0.0742     0.0742   Glioblastoma,IDHwt - Oligodendroglioma,IDHmut-codel



