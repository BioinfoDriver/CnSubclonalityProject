
# gold set
gold.set <- readRDS( file='/pub5/xiaoyun/Jobs/J22/CopyNumberClonalityProject/Resource/CuratedData/gold_set.rds')

pan.glio.cli.mol.data <- readRDS(file='/pub5/xiaoyun/Jobs/J22/CopyNumberClonalityProject/Resource/CuratedData/tcga_glioma_cli_mol.rds')
gli.cn.alt.frac <- readRDS(file='/pub5/xiaoyun/Jobs/J22/CopyNumberClonalityProject/Resource/CuratedData/tcga_gli_cn_alt_frac.rds') 

gli.cn.alt.frac <- cbind(gli.cn.alt.frac[paste0(gold.set, '-01'), ], pan.glio.cli.mol.data[gold.set, ])
saveRDS(gli.cn.alt.frac, file='/pub5/xiaoyun/Jobs/J22/CopyNumberClonalityProject/Resource/CuratedData/gli_glod_cn_alt_frac.rds')



