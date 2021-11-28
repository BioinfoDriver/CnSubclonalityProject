
# load data
gli.cn.alt.frac <- readRDS(file='/data/gli_glod_cn_alt_frac.rds')

# ABSOLUTE-based tumour purity
abs.puri.ploi <- readRDS('/data/tcga_gli_puri_ploi.rds')
abs.puri.ploi <- subset(abs.puri.ploi, purity>=0.65)

gli.cn.alt.frac <- gli.cn.alt.frac[rownames(abs.puri.ploi), ]


# subtype associated with genome alt fraction
SubtypeSubcloStat <- function(dat, alt.feac, sub.feac, label){
 
 dat <- dat[!is.na(dat[, sub.feac]), ]
 dat[, grep(sub.feac, colnames(dat))] <- factor(dat[, sub.feac], levels = label)
 
 p.value <- kruskal.test(get(alt.feac) ~ get(sub.feac), data = dat)$p.value
 res <- dunn.test::dunn.test(dat[, alt.feac], dat[, sub.feac], method = 'bh')
 res <- dplyr::bind_rows(res[2:length(res)]) %>% mutate(var1_var2 = paste(alt.feac, sub.feac, sep = " - "))
 
 return(list(p.value, res))
}


SubtypeComByFeacture <- function(input.dat, features, subtype, labels){
	
 res.list <- lapply(features, function(feature){
  res <- SubtypeSubcloStat(dat=input.dat, alt.feac=feature, sub.feac=subtype, label=labels)
  return(res)
 })
 pvalues <- setNames(sapply(res.list, function(x) x[[1]]), features)
 com.pvalues <- sapply(res.list, function(x) x[2])
 com.pvalues <- do.call(rbind, com.pvalues)

  return(list(pvalues, com.pvalues))
}


features <- c('non_neutral_genome_frac', 'subclo_genome_frac', 'clo_genome_frac', 'subclo_cn_alt_frac')

SubtypeComByFeacture(gli.cn.alt.frac, features, subtype='Integrated_Diagnoses', 
 labels=c('Oligodendroglioma,IDHmut-codel', 'Astrocytoma,IDHmut', 'Glioblastoma,IDHwt'))
 

# p-value = 8.76e-46
# p-value = 1.20e-13
# p-value = 5.75e-38
# p-value = 7.20e-18

             # Z            P   P.adjusted                                         comparisons               var1_var2
# 1  -11.9486940 3.297674e-33 9.893021e-33             Astrocytoma,IDHmut - Glioblastoma,IDHwt non_neutral_genome_frac
# 2    0.8230647 2.052356e-01 2.052356e-01 Astrocytoma,IDHmut - Oligodendroglioma,IDHmut-codel non_neutral_genome_frac
# 3   11.8976748 6.084536e-33 9.126805e-33 Glioblastoma,IDHwt - Oligodendroglioma,IDHmut-codel non_neutral_genome_frac
# 4   -0.2056264 4.185414e-01 4.185414e-01             Astrocytoma,IDHmut - Glioblastoma,IDHwt      subclo_genome_frac
# 5    6.5270060 3.354867e-11 5.032301e-11 Astrocytoma,IDHmut - Oligodendroglioma,IDHmut-codel      subclo_genome_frac
# 6    7.2390983 2.258388e-13 6.775165e-13 Glioblastoma,IDHwt - Oligodendroglioma,IDHmut-codel      subclo_genome_frac
# 7  -12.6261552 7.575083e-37 2.272525e-36             Astrocytoma,IDHmut - Glioblastoma,IDHwt         clo_genome_frac
# 8   -3.4527176 2.774848e-04 2.774848e-04 Astrocytoma,IDHmut - Oligodendroglioma,IDHmut-codel         clo_genome_frac
# 9    7.9036803 1.353934e-15 2.030902e-15 Glioblastoma,IDHwt - Oligodendroglioma,IDHmut-codel         clo_genome_frac
# 10   5.4728639 2.214101e-08 3.321151e-08             Astrocytoma,IDHmut - Glioblastoma,IDHwt      subclo_cn_alt_frac
# 11   8.8122477 6.133061e-19 1.839918e-18 Astrocytoma,IDHmut - Oligodendroglioma,IDHmut-codel      subclo_cn_alt_frac
# 12   4.4755579 3.810602e-06 3.810602e-06 Glioblastoma,IDHwt - Oligodendroglioma,IDHmut-codel      subclo_cn_alt_frac



