
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

SubtypeComByFeacture(gli.cn.alt.frac, features, subtype='IDH_CODEL_SUBTYPE', labels=c('IDHwt', 'IDHmut-non-codel', 'IDHmut-codel'))
# p-value = 1.46e-52
# p-value = 1.35e-18
# p-value = 2.15e-50
# p-value = 6.37e-29
            # Z            P   P.adjusted                     comparisons	var1_var2
# 1   -1.772663 3.814229e-02 3.814229e-02 IDHmut-codel - IDHmut-non-codel	non_neutral_genome_frac
# 2  -12.650761 5.539814e-37 8.309721e-37            IDHmut-codel - IDHwt	non_neutral_genome_frac
# 3  -12.714421 2.458682e-37 7.376047e-37        IDHmut-non-codel - IDHwt	non_neutral_genome_frac
# 4   -8.583403 4.605365e-18 1.381609e-17 IDHmut-codel - IDHmut-non-codel	subclo_genome_frac
# 5   -7.895358 1.447414e-15 2.171121e-15            IDHmut-codel - IDHwt	subclo_genome_frac
# 6    1.443094 7.449706e-02 7.449706e-02        IDHmut-non-codel - IDHwt	subclo_genome_frac
# 7    4.287447 9.036907e-06 9.036907e-06 IDHmut-codel - IDHmut-non-codel	clo_genome_frac
# 8   -7.993960 6.533606e-16 9.800409e-16            IDHmut-codel - IDHwt	clo_genome_frac
# 9  -14.816689 5.714038e-50 1.714211e-49        IDHmut-non-codel - IDHwt	clo_genome_frac
# 10 -10.969556 2.676745e-28 8.030234e-28 IDHmut-codel - IDHmut-non-codel	subclo_cn_alt_frac
# 11  -4.979999 3.179229e-07 3.179229e-07            IDHmut-codel - IDHwt	subclo_cn_alt_frac
# 12   7.878348 1.658694e-15 2.488042e-15        IDHmut-non-codel - IDHwt	subclo_cn_alt_frac


SubtypeComByFeacture(subset(gli.cn.alt.frac, cancer_type == 'LGG'), features, 
 subtype='IDH_CODEL_SUBTYPE', labels=c('IDHwt', 'IDHmut-non-codel', 'IDHmut-codel'))
# p-value = 3.65e-14
# p-value = 1.31e-16
# p-value = 4.69e-25
# p-value = 1.32e-26
             # Z            P   P.adjusted                     comparisons	var1_var2
# 1   -0.7199845 2.357673e-01 2.357673e-01 IDHmut-codel - IDHmut-non-codel	non_neutral_genome_frac
# 2   -7.2577502 1.967908e-13 5.903724e-13            IDHmut-codel - IDHwt	non_neutral_genome_frac
# 3   -7.1797558 3.491803e-13 5.237704e-13        IDHmut-non-codel - IDHwt	non_neutral_genome_frac
# 4   -8.5366008 6.911406e-18 2.073422e-17 IDHmut-codel - IDHmut-non-codel	subclo_genome_frac
# 5   -3.5388218 2.009586e-04 3.014378e-04            IDHmut-codel - IDHwt	subclo_genome_frac
# 6    3.2543450 5.682708e-04 5.682708e-04        IDHmut-non-codel - IDHwt	subclo_genome_frac
# 7    6.1185858 4.720473e-10 7.080710e-10 IDHmut-codel - IDHmut-non-codel	clo_genome_frac
# 8   -4.7636795 9.504718e-07 9.504718e-07            IDHmut-codel - IDHwt	clo_genome_frac
# 9  -10.1518515 1.625678e-24 4.877035e-24        IDHmut-non-codel - IDHwt	clo_genome_frac
# 10 -10.3743375 1.622275e-25 4.866826e-25 IDHmut-codel - IDHmut-non-codel	subclo_cn_alt_frac
# 11  -1.8574954 3.162035e-02 3.162035e-02            IDHmut-codel - IDHwt	subclo_cn_alt_frac
# 12   6.5718445 2.484788e-11 3.727182e-11        IDHmut-non-codel - IDHwt	subclo_cn_alt_frac
