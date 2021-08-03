
gistic.alt <- readRDS(file='/pub5/xiaoyun/Jobs/J22/CopyNumberClonalityProject/Resource/CuratedData/wgistic_alt.rds')


gold.set <- readRDS('/pub5/xiaoyun/Jobs/J22/CopyNumberClonalityProject/Resource/CuratedData/gold_set.rds')
gli.cli.mol.data <- readRDS(file='/pub5/xiaoyun/Jobs/J22/CopyNumberClonalityProject/Resource/CuratedData/tcga_glioma_cli_mol.rds')

gli.cli.mol.data <- gli.cli.mol.data[gold.set, ]
rownames(gli.cli.mol.data) <- paste(rownames(gli.cli.mol.data), '01', sep = '-')

# samples in each subtype
subtype.samples <- lapply(c('IDHmut-codel', 'IDHmut-non-codel', 'IDHwt'), function(subtype){
 
 return(rownames(subset(gli.cli.mol.data, IDH_CODEL_SUBTYPE==subtype)))
 
})


# Peaks with high alteration frequency
subtype.peaks <- lapply(subtype.samples, function(samples){
 
 alt.n <- rowSums(gistic.alt[, samples])
 alt.cut <- ncol(gistic.alt[, samples])*0.05
 return(names(alt.n)[abs(alt.n) >= alt.cut])
})

names(subtype.peaks) <- c('IDHmut-codel', 'IDHmut-non-codel', 'IDHwt')

saveRDS(subtype.peaks, file='/pub5/xiaoyun/Jobs/J22/CopyNumberClonalityProject/Results/Section4/Resources/high_fre_alt_peaks.rds')
