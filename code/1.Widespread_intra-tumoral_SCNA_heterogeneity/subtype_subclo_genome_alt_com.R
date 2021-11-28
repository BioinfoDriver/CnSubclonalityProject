
library(dplyr)
# load data
gli.cn.alt.frac <- readRDS(file='/data/gli_glod_cn_alt_frac.rds')

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
  res <- SubtypeSubcloStat(input.dat, alt.feac=feature, sub.feac=subtype, label=labels)
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
# p-value = 1.35e-56
# p-value = 6.46e-19
# p-value = 5.94e-53
# p-value = 1.19e-28
            # Z            P   P.adjusted                                         comparisons               var1_var2
# 1  -13.277362 1.566175e-40 4.698525e-40             Astrocytoma,IDHmut - Glioblastoma,IDHwt non_neutral_genome_frac
# 2    1.802318 3.574771e-02 3.574771e-02 Astrocytoma,IDHmut - Oligodendroglioma,IDHmut-codel non_neutral_genome_frac
# 3   13.179785 5.735865e-40 8.603797e-40 Glioblastoma,IDHwt - Oligodendroglioma,IDHmut-codel non_neutral_genome_frac
# 4    1.196388 1.157726e-01 1.157726e-01             Astrocytoma,IDHmut - Glioblastoma,IDHwt      subclo_genome_frac
# 5    8.591265 4.300860e-18 1.290258e-17 Astrocytoma,IDHmut - Oligodendroglioma,IDHmut-codel      subclo_genome_frac
# 6    8.033300 4.744277e-16 7.116416e-16 Glioblastoma,IDHwt - Oligodendroglioma,IDHmut-codel      subclo_genome_frac
# 7  -15.194546 1.921579e-52 5.764736e-52             Astrocytoma,IDHmut - Glioblastoma,IDHwt         clo_genome_frac
# 8   -4.309137 8.194627e-06 8.194627e-06 Astrocytoma,IDHmut - Oligodendroglioma,IDHmut-codel         clo_genome_frac
# 9    8.371088 2.854433e-17 4.281650e-17 Glioblastoma,IDHwt - Oligodendroglioma,IDHmut-codel         clo_genome_frac
# 10   7.717575 5.928189e-15 8.892284e-15             Astrocytoma,IDHmut - Glioblastoma,IDHwt      subclo_cn_alt_frac
# 11  10.942747 3.599250e-28 1.079775e-27 Astrocytoma,IDHmut - Oligodendroglioma,IDHmut-codel      subclo_cn_alt_frac
# 12   4.969489 3.356479e-07 3.356479e-07 Glioblastoma,IDHwt - Oligodendroglioma,IDHmut-codel      subclo_cn_alt_frac



