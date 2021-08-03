
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

SubtypeComByFeacture(gli.cn.alt.frac, features, subtype='IDH_CODEL_SUBTYPE', labels=c('IDHwt', 'IDHmut-non-codel', 'IDHmut-codel'))
# p-value = 1.71e-44
# p-value = 1.90e-13
# p-value = 1.59e-37
# p-value = 5.93e-18
             # Z            P   P.adjusted                     comparisons	var1_var2
# 1   -0.7980388 2.124240e-01 2.124240e-01 IDHmut-codel - IDHmut-non-codel	non_neutral_genome_frac
# 2  -11.6780000 8.256762e-32 1.238514e-31            IDHmut-codel - IDHwt	non_neutral_genome_frac
# 3  -11.7472253 3.648792e-32 1.094638e-31        IDHmut-non-codel - IDHwt	non_neutral_genome_frac
# 4   -6.5044926 3.897802e-11 5.846704e-11 IDHmut-codel - IDHmut-non-codel	subclo_genome_frac
# 5   -7.1777200 3.544178e-13 1.063253e-12            IDHmut-codel - IDHwt	subclo_genome_frac
# 6   -0.1155942 4.539871e-01 4.539871e-01        IDHmut-non-codel - IDHwt	subclo_genome_frac
# 7    3.4543729 2.757871e-04 2.757871e-04 IDHmut-codel - IDHmut-non-codel	clo_genome_frac
# 8   -7.7857174 3.465937e-15 5.198906e-15            IDHmut-codel - IDHwt	clo_genome_frac
# 9  -12.5412649 2.219268e-36 6.657804e-36        IDHmut-non-codel - IDHwt	clo_genome_frac
# 10  -8.8124176 6.123769e-19 1.837131e-18 IDHmut-codel - IDHmut-non-codel	subclo_cn_alt_frac
# 11  -4.3890839 5.691457e-06 5.691457e-06            IDHmut-codel - IDHwt	subclo_cn_alt_frac
# 12   5.6411859 8.444147e-09 1.266622e-08        IDHmut-non-codel - IDHwt	subclo_cn_alt_frac

SubtypeComByFeacture(subset(gli.cn.alt.frac, cancer_type == 'LGG'), features, 
 subtype='IDH_CODEL_SUBTYPE', labels=c('IDHwt', 'IDHmut-non-codel', 'IDHmut-codel'))
# p-value = 1.35921230184229e-11
# p-value = 4.86187792225113e-10
# p-value = 2.32943572856334e-17
# p-value = 5.40992688272558e-17

            # Z            P   P.adjusted                     comparisons var1_var2
# 1   0.1459325 4.419873e-01 4.419873e-01 IDHmut-codel - IDHmut-non-codel non_neutral_genome_frac
# 2  -6.4652599 5.056242e-11 7.584362e-11            IDHmut-codel - IDHwt non_neutral_genome_frac
# 3  -6.7146192 9.427898e-12 2.828369e-11        IDHmut-non-codel - IDHwt non_neutral_genome_frac
# 4  -6.5428842 3.017177e-11 9.051532e-11 IDHmut-codel - IDHmut-non-codel subclo_genome_frac
# 5  -2.2715897 1.155565e-02 1.155565e-02            IDHmut-codel - IDHwt subclo_genome_frac
# 6   2.4526772 7.089876e-03 1.063481e-02        IDHmut-non-codel - IDHwt subclo_genome_frac
# 7   5.0894911 1.795129e-07 2.692693e-07 IDHmut-codel - IDHmut-non-codel clo_genome_frac
# 8  -4.5961238 2.152115e-06 2.152115e-06            IDHmut-codel - IDHwt clo_genome_frac
# 9  -8.4115856 2.022523e-17 6.067570e-17        IDHmut-non-codel - IDHwt clo_genome_frac
# 10 -8.2654199 6.960195e-17 2.088059e-16 IDHmut-codel - IDHmut-non-codel subclo_cn_alt_frac
# 11 -0.8301967 2.032138e-01 2.032138e-01            IDHmut-codel - IDHwt subclo_cn_alt_frac
# 12  5.1828885 1.092378e-07 1.638567e-07        IDHmut-non-codel - IDHwt subclo_cn_alt_frac

