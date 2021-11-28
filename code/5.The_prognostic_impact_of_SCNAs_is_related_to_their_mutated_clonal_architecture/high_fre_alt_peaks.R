
gistic.alt <- readRDS(file='/data/wgistic_alt.rds')


gold.set <- readRDS('/data/gold_set.rds')
gli.cli.mol.data <- readRDS(file='/data/tcga_glioma_cli_mol.rds')

gli.cli.mol.data <- gli.cli.mol.data[gold.set, ]
rownames(gli.cli.mol.data) <- paste(rownames(gli.cli.mol.data), '01', sep = '-')

# samples in each subtype
subtype.samples <- lapply(c('Oligodendroglioma,IDHmut-codel', 'Astrocytoma,IDHmut', 'Glioblastoma,IDHwt'), function(subtype){
 
 return(rownames(subset(gli.cli.mol.data, Integrated_Diagnoses==subtype)))
 
})


# Peaks with high alteration frequency
subtype.peaks <- lapply(subtype.samples, function(samples){
 
 alt.n <- rowSums(gistic.alt[, samples])
 alt.cut <- ncol(gistic.alt[, samples])*0.05
 return(names(alt.n)[abs(alt.n) >= alt.cut])
})

names(subtype.peaks) <- c('Oligodendroglioma,IDHmut-codel', 'Astrocytoma,IDHmut', 'Glioblastoma,IDHwt')

saveRDS(subtype.peaks, file='/data/high_fre_alt_peaks.rds')
