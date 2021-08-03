

DatStas <- function(cn.alt.frac, charac){
	
	dat.range <- aggregate(get(charac) ~ IDH_CODEL_SUBTYPE, data = cn.alt.frac, range)
	dat.mean <- aggregate(get(charac) ~ IDH_CODEL_SUBTYPE, data = cn.alt.frac, mean)
	dat.median <- aggregate(get(charac) ~ IDH_CODEL_SUBTYPE, data = cn.alt.frac, median)	
	
	mer.dat <- Reduce(function(x ,y) merge(x, y, by="IDH_CODEL_SUBTYPE"), list(dat.mean, dat.median, dat.range))
	mer.dat <- as.matrix(mer.dat)
	colnames(mer.dat) <- c('Subtype', 'Mean', 'Median', 'Min', 'Max')
	return(mer.dat)
}

gli.cn.alt.frac <- readRDS(file='/pub5/xiaoyun/Jobs/J22/CopyNumberClonalityProject/Resource/CuratedData/gli_glod_cn_alt_frac.rds')
gli.cn.alt.frac$IDH_CODEL_SUBTYPE <- factor(gli.cn.alt.frac$IDH_CODEL_SUBTYPE, 
 levels = c('IDHwt', 'IDHmut-non-codel', 'IDHmut-codel'))


# > summary(gli.cn.alt.frac$subclo_cn_alt_frac)
   # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00000 0.02705 0.25358 0.32602 0.53793 1.00000 


# > aggregate(clo_genome_frac ~ IDH_CODEL_SUBTYPE, data = gli.cn.alt.frac, summary)
  # IDH_CODEL_SUBTYPE clo_genome_frac.Min. clo_genome_frac.1st Qu.
# 1             IDHwt           0.00000000              0.09340608
# 2  IDHmut-non-codel           0.00000000              0.02071675
# 3      IDHmut-codel           0.00000000              0.05502464
  # clo_genome_frac.Median clo_genome_frac.Mean clo_genome_frac.3rd Qu.
# 1             0.14128243           0.15340080              0.18455951
# 2             0.04512781           0.06364297              0.06850856
# 3             0.05759575           0.08327542              0.08999306
  # clo_genome_frac.Max.
# 1           0.73475835
# 2           0.70775844
# 3           0.66338165


# > aggregate(subclo_genome_frac ~ IDH_CODEL_SUBTYPE, data = gli.cn.alt.frac, summary)
  # IDH_CODEL_SUBTYPE subclo_genome_frac.Min. subclo_genome_frac.1st Qu.
# 1             IDHwt            0.000000e+00               5.836879e-03
# 2  IDHmut-non-codel            0.000000e+00               1.528533e-02
# 3      IDHmut-codel            0.000000e+00               1.030554e-05
  # subclo_genome_frac.Median subclo_genome_frac.Mean subclo_genome_frac.3rd Qu.
# 1              3.313274e-02            7.101964e-02               1.015599e-01
# 2              4.086024e-02            6.554140e-02               8.352779e-02
# 3              5.974410e-04            2.433794e-02               3.232603e-02
  # subclo_genome_frac.Max.
# 1            6.067040e-01
# 2            5.355606e-01
# 3            3.782040e-01


# > aggregate(subclo_cn_alt_frac ~ IDH_CODEL_SUBTYPE, data = gli.cn.alt.frac, summary)
  # IDH_CODEL_SUBTYPE subclo_cn_alt_frac.Min. subclo_cn_alt_frac.1st Qu.
# 1             IDHwt            0.0000000000               0.0359915196
# 2  IDHmut-non-codel            0.0000000000               0.2455280269
# 3      IDHmut-codel            0.0000000000               0.0001033351
  # subclo_cn_alt_frac.Median subclo_cn_alt_frac.Mean subclo_cn_alt_frac.3rd Qu.
# 1              0.1824400297            0.2801004902               0.4259015533
# 2              0.4683398717            0.4829832520               0.7251076293
# 3              0.0090223726            0.1741072986               0.3448660897
  # subclo_cn_alt_frac.Max.
# 1            1.0000000000
# 2            1.0000000000
# 3            1.0000000000

