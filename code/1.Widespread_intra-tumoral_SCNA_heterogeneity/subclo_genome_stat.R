

DatStas <- function(cn.alt.frac, charac){
	
	dat.range <- aggregate(get(charac) ~ Integrated_Diagnoses, data = cn.alt.frac, range)
	dat.mean <- aggregate(get(charac) ~ Integrated_Diagnoses, data = cn.alt.frac, mean)
	dat.median <- aggregate(get(charac) ~ Integrated_Diagnoses, data = cn.alt.frac, median)	
	
	mer.dat <- Reduce(function(x ,y) merge(x, y, by="Integrated_Diagnoses"), list(dat.mean, dat.median, dat.range))
	mer.dat <- as.matrix(mer.dat)
	colnames(mer.dat) <- c('Subtype', 'Mean', 'Median', 'Min', 'Max')
	return(mer.dat)
}

gli.cn.alt.frac <- readRDS(file='/data/gli_glod_cn_alt_frac.rds')
gli.cn.alt.frac$Integrated_Diagnoses <- factor(gli.cn.alt.frac$Integrated_Diagnoses, 
 levels = c('Oligodendroglioma,IDHmut-codel', 'Astrocytoma,IDHmut', 'Glioblastoma,IDHwt'))


# > summary(gli.cn.alt.frac$subclo_cn_alt_frac)
   # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00000 0.02705 0.25451 0.32793 0.54267 1.00000 



# > aggregate(clo_genome_frac ~ Integrated_Diagnoses, data = gli.cn.alt.frac, summary)
            # Integrated_Diagnoses clo_genome_frac.Min. clo_genome_frac.1st Qu.
# 1 Oligodendroglioma,IDHmut-codel           0.00000000              0.05502464
# 2             Astrocytoma,IDHmut           0.00000000              0.02071675
# 3             Glioblastoma,IDHwt           0.00000000              0.10249448
  # clo_genome_frac.Median clo_genome_frac.Mean clo_genome_frac.3rd Qu.
# 1             0.05759575           0.08327542              0.08999306
# 2             0.04512781           0.06364297              0.06850856
# 3             0.14245336           0.15536972              0.18455951
  # clo_genome_frac.Max.
# 1           0.66338165
# 2           0.70775844
# 3           0.73475835


# > aggregate(subclo_genome_frac ~ Integrated_Diagnoses, data = gli.cn.alt.frac, summary)
            # Integrated_Diagnoses subclo_genome_frac.Min.
# 1 Oligodendroglioma,IDHmut-codel            0.000000e+00
# 2             Astrocytoma,IDHmut            0.000000e+00
# 3             Glioblastoma,IDHwt            0.000000e+00
  # subclo_genome_frac.1st Qu. subclo_genome_frac.Median subclo_genome_frac.Mean
# 1               1.030554e-05              5.974410e-04            2.433794e-02
# 2               1.528533e-02              4.086024e-02            6.554140e-02
# 3               6.463717e-03              3.560426e-02            7.232122e-02
  # subclo_genome_frac.3rd Qu. subclo_genome_frac.Max.
# 1               3.232603e-02            3.782040e-01
# 2               8.352779e-02            5.355606e-01
# 3               1.041065e-01            6.067040e-01


# > aggregate(subclo_cn_alt_frac ~ Integrated_Diagnoses, data = gli.cn.alt.frac, summary)
            # Integrated_Diagnoses subclo_cn_alt_frac.Min.
# 1 Oligodendroglioma,IDHmut-codel            0.0000000000
# 2             Astrocytoma,IDHmut            0.0000000000
# 3             Glioblastoma,IDHwt            0.0000000000
  # subclo_cn_alt_frac.1st Qu. subclo_cn_alt_frac.Median subclo_cn_alt_frac.Mean
# 1               0.0001033351              0.0090223726            0.1741072986
# 2               0.2455280269              0.4683398717            0.4829832520
# 3               0.0359915196              0.1811365481            0.2816177899
  # subclo_cn_alt_frac.3rd Qu. subclo_cn_alt_frac.Max.
# 1               0.3448660897            1.0000000000
# 2               0.7251076293            1.0000000000
# 3               0.4338482319            1.0000000000


