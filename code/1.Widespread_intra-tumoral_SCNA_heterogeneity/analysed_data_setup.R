
# gold set
gold.set <- readRDS( file='/data/gold_set.rds')

pan.glio.cli.mol.data <- readRDS(file='/data/tcga_glioma_cli_mol.rds')
gli.cn.alt.frac <- readRDS(file='/data/tcga_gli_cn_alt_frac.rds') 

gli.cn.alt.frac <- cbind(gli.cn.alt.frac[paste0(gold.set, '-01'), ], pan.glio.cli.mol.data[gold.set, ])
saveRDS(gli.cn.alt.frac, file='/data/gli_glod_cn_alt_frac.rds')



