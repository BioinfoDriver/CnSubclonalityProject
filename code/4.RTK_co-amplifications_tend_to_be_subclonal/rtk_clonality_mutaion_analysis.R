
# load mutation
load('/data/PureGliomaMutData.RData')

pan.glioma.cancer.gene <- read.csv(file='/data/OriginalData/cancer_genes_multiple.txt', 
 header=TRUE, sep='\t', stringsAsFactors=FALSE)
driver.genes <- pan.glioma.cancer.gene$Cell_2016_PanGlioma

gold.set <- readRDS('/data/gold_set.rds')

# Filter
gli.sig.mut <- subset(pure.glioma.mut.data, Patient %in% paste0(gold.set, '-01'))
gli.sig.mut <- subset(gli.sig.mut, Hugo_Symbol %in% driver.genes)
gli.sig.mut <- subset(gli.sig.mut, !(Variant_Classification %in% c("3'UTR", "5'UTR", "Intron", "Silent")))


# mutated in more than 2% samples
mut.freq <- table(gli.sig.mut$Hugo_Symbol)/760
driver.genes <- names(mut.freq)[mut.freq >= 0.02]


sig.mut.mat <- sapply(unique(gli.sig.mut$Patient), function(patient){
 tmp <- subset(gli.sig.mut, Patient == patient)$Hugo_Symbol
 mut.index <- as.numeric(driver.genes %in% tmp)
 
 return(mut.index)
})

sig.mut.mat <- as.data.frame(t(sig.mut.mat))
colnames(sig.mut.mat) <- driver.genes



# RTK clonality subtype
rtks <- read.csv(file='/data/RTKs.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)

gold.set <- readRDS('/data/gold_set.rds')
gene.het <- readRDS(file='/data/gene.het.rds')
gene.het <- gene.het[, paste(gold.set, '01', sep = '-')]

tcga.cli.data <- readRDS('/data/tcga_glioma_cli_mol.rds')
tcga.cli.data <- tcga.cli.data[gold.set, ]


rtk.het <- gene.het[rtks$Approved.symbol, ]
rtk.clo.sam <- colnames(rtk.het[, colSums(rtk.het == 1) >= 1])
rtk.sub.sam <- colnames(rtk.het[, colSums(rtk.het == 2) >= 1])

com.sam <- intersect(rtk.clo.sam, rtk.sub.sam)
rtk.clo.sam <- setdiff(rtk.clo.sam, com.sam)
rtk.sub.sam <- setdiff(rtk.sub.sam, com.sam)


rtk.group <- data.frame(group = rep(c('clonal', 'subclonal'), times=c(length(rtk.clo.sam), length(rtk.sub.sam))))
rownames(rtk.group) <- c(rtk.clo.sam, rtk.sub.sam)


# merge
# length(intersect(pure.glioma.mut.data$Patient, rownames(rtk.group)))
# 206
sig.mut.mat <- merge(sig.mut.mat, rtk.group, by='row.names', all.y=TRUE)
sig.mut.mat <- tibble::column_to_rownames(sig.mut.mat, var = "Row.names")
sig.mut.mat[is.na(sig.mut.mat)] <- 0
sig.mut.mat$IDH <- sig.mut.mat$IDH1 + sig.mut.mat$IDH2
sig.mut.mat <- subset(sig.mut.mat, select = -c(IDH1, IDH2))

################## all patients
# sapply(setdiff(colnames(sig.mut.mat), 'group'), function(driver.gene){
	# table(sig.mut.mat[, c('group', driver.gene)])
# })

p.values <- sapply(setdiff(colnames(sig.mut.mat), 'group'), function(driver.gene){
   
	num <- table(sig.mut.mat[, driver.gene])
	char <- names(num)[num >= 3]
	dat.del <- subset(sig.mut.mat, get(driver.gene) %in% char)
	
	if(dplyr::n_distinct(dat.del[, driver.gene])<=1){
	 p <- NA
	
	}else{
	 stat <- table(dat.del[, c('group', driver.gene)])
	 p <- fisher.test(stat, alternative = "two.sided")$p.value
	
	}
	
	return(p)

})

p.adjust(p.values, method='fdr')


################## Glioblastoma,IDHwt
gbm.idhwt <- paste0(rownames(subset(tcga.cli.data, Integrated_Diagnoses == 'Glioblastoma,IDHwt')), '-01')
gbm.idhwt.mut.mat <- sig.mut.mat[intersect(rownames(sig.mut.mat), gbm.idhwt), ]

# sapply(setdiff(colnames(sig.mut.mat), 'group'), function(driver.gene){
	# table(gbm.idhwt.mut.mat[, c('group', driver.gene)])
# })


gbm.idhwt.p.values <- sapply(setdiff(colnames(sig.mut.mat), 'group'), function(driver.gene){
   
	num <- table(gbm.idhwt.mut.mat[, driver.gene])
	char <- names(num)[num >= 3]
	dat.del <- subset(gbm.idhwt.mut.mat, get(driver.gene) %in% char)
	
	if(dplyr::n_distinct(dat.del[, driver.gene])<=1){
	 p <- NA
	
	}else{
	 stat <- table(dat.del[, c('group', driver.gene)])
	 p <- fisher.test(stat, alternative = "two.sided")$p.value
	
	}
	
	return(p)

})

p.adjust(gbm.idhwt.p.values, method='fdr')


################## Glioblastoma,IDHwt
as.idhmut <- paste0(rownames(subset(tcga.cli.data, Integrated_Diagnoses == 'Astrocytoma,IDHmut')), '-01')
as.idhmut.mut.mat <- sig.mut.mat[intersect(rownames(sig.mut.mat), as.idhmut), ]

# sapply(setdiff(colnames(sig.mut.mat), 'group'), function(driver.gene){
	# table(as.idhmut.mut.mat[, c('group', driver.gene)])
# })


as.idhmut.p.values <- sapply(setdiff(colnames(sig.mut.mat), 'group'), function(driver.gene){
   
	num <- table(as.idhmut.mut.mat[, driver.gene])
	char <- names(num)[num >= 3]
	dat.del <- subset(as.idhmut.mut.mat, get(driver.gene) %in% char)
	
	if(dplyr::n_distinct(dat.del[, driver.gene])<=1){
	 p <- NA
	
	}else{
	 stat <- table(dat.del[, c('group', driver.gene)])
	 p <- fisher.test(stat, alternative = "two.sided")$p.value
	
	}
	
	return(p)

})

p.adjust(as.idhmut.p.values, method='fdr')


