
# load data
load('/data/PureGliomaMutData.RData')
load('/data/DDRGene.RData')
ddr.genes <- ddr.gene$Gene.Symbol

driver.genes <- read.csv(file='/data/OriginalData/cancer_genes_multiple.txt', 
 header=TRUE, sep='\t', stringsAsFactors=FALSE)
driver.genes <- driver.genes$Cell_2016_PanGlioma

gold.set <- readRDS('/data/gold_set.rds')
gli.cn.alt.frac <- readRDS(file='/data/gli_glod_cn_alt_frac.rds')
gli.cn.alt.frac <- gli.cn.alt.frac[, c('Integrated_Diagnoses', 'subclo_genome_frac','clo_genome_frac')]

# high purity
# abs.puri.ploi <- readRDS('/data/tcga_gli_puri_ploi.rds')
# abs.puri.ploi <- subset(abs.puri.ploi, purity>=0.65)
# gli.cn.alt.frac <- gli.cn.alt.frac[intersect(rownames(abs.puri.ploi), rownames(gli.cn.alt.frac)), ]


# Filter
gli.sig.mut <- subset(pure.glioma.mut.data, Patient %in% paste0(gold.set, '-01'))
gli.sig.mut <- subset(gli.sig.mut, Hugo_Symbol %in% unique(c(ddr.genes, driver.genes)))
gli.sig.mut <- subset(gli.sig.mut, !(Variant_Classification %in% c("3'Flank", "3'UTR", "5'UTR", "Intron", "Silent")))


# mutated in more than 2% samples
mut.freq <- table(gli.sig.mut$Hugo_Symbol)/760
high.freq.genes <- names(mut.freq)[mut.freq >= 0.02]

sig.mut.mat <- sapply(unique(gli.sig.mut$Patient), function(patient){
 tmp <- subset(gli.sig.mut, Patient == patient)$Hugo_Symbol
 mut.index <- as.numeric(high.freq.genes %in% tmp)
 
 return(mut.index)
})

sig.mut.mat <- as.data.frame(t(sig.mut.mat))
colnames(sig.mut.mat) <- high.freq.genes
sig.mut.mat$IDH <- sig.mut.mat$IDH1 + sig.mut.mat$IDH2
sig.mut.mat <- subset(sig.mut.mat, select = -c(IDH1, IDH2))


# associated with clonal and subclonal scna
gli.cn.alt.frac <- merge(gli.cn.alt.frac, sig.mut.mat, by='row.names', all.y=TRUE)
gli.cn.alt.frac <- tibble::column_to_rownames(gli.cn.alt.frac, var = "Row.names")
gli.cn.alt.frac[is.na(gli.cn.alt.frac)] <- 0


stat <- lapply(c('Glioblastoma,IDHwt', 'Astrocytoma,IDHmut', 'Oligodendroglioma,IDHmut-codel'), function(subtype){
 
 sub.stat <- subset(gli.cn.alt.frac, Integrated_Diagnoses == subtype)
 feac.pvalues <- lapply(c('subclo_genome_frac', 'clo_genome_frac'), function(feac){
  
  pvalues <- sapply(names(sig.mut.mat), function(mol.feac){
   
  
   num <- table(sub.stat[, mol.feac])
   char <- names(num)[num >= 5]
   dat.del <- subset(sub.stat, get(mol.feac) %in% char)
    
   if(dplyr::n_distinct(dat.del[, mol.feac])<=1){
     p <- NA
	 
   }else{
     p <- wilcox.test(get(feac) ~ get(mol.feac), dat.del)$p.value
    
    }
   
   return(p)
  
  })
  
 qvalues <- p.adjust(pvalues, 'fdr') 
 pqvalues <- data.frame(pvalues, qvalues)
 rownames(pqvalues) <- names(sig.mut.mat)
 return(pqvalues)
 
 })
 
 feac.pvalues <- do.call(cbind, feac.pvalues)
 colnames(feac.pvalues) <- paste0(rep(c('subclo', 'clonal'), each=2), c('_pvalue', '_qvalue'))
 
 return(feac.pvalues)
})


################## plot
gli.cn.alt.frac$TP53 <- factor(gli.cn.alt.frac$TP53)
library('ggplot2')
library('ggpubr')
# p1 <- ggplot(data = subset(gli.cn.alt.frac, Integrated_Diagnoses == 'Glioblastoma,IDHwt'), aes(TP53, subclo_genome_frac)) +
 # geom_boxplot(aes(shape = TP53), outlier.shape = NA) + geom_jitter(shape=16, position=position_jitter(0.2)) + 
 # xlab('Subtype') + ylab('Burden of subclonal SCNAs') + guides(shape = guide_legend(title= 'Subtype')) + 
 # theme(axis.text.x = element_text(angle=30, vjust=0.5)) + stat_compare_means(label.x.npc = 'left', label.y.npc = 'top')

p1 <- ggboxplot(data = subset(gli.cn.alt.frac, Integrated_Diagnoses == 'Glioblastoma,IDHwt'), 
 x = 'TP53', y = 'subclo_genome_frac', color = 'TP53', xlab = FALSE, ylab = 'Burden of subclonal SCNAs', legend = 'none', 
 add = "jitter", shape = 'TP53', outlier.shape = NA, alpha = 1.0, add.params=list(size = 0.8)) + 
 stat_compare_means(label.x.npc = 'left', label.y.npc = 'top') + theme(aspect.ratio = 1)


# p2 <- ggplot(data = subset(gli.cn.alt.frac, Integrated_Diagnoses == 'Astrocytoma,IDHmut'), aes(TP53, subclo_genome_frac)) +
 # geom_boxplot(aes(shape = TP53), outlier.shape = NA) + geom_jitter(shape=16, position=position_jitter(0.2)) + 
 # xlab('Subtype') + ylab('Burden of subclonal SCNAs') + guides(shape = guide_legend(title= 'Subtype')) + 
 # theme(axis.text.x = element_text(angle=30, vjust=0.5)) + stat_compare_means(label.x.npc = 'left', label.y.npc = 'top')


p2 <- ggboxplot(data = subset(gli.cn.alt.frac, Integrated_Diagnoses == 'Astrocytoma,IDHmut'), 
 x = 'TP53', y = 'subclo_genome_frac', color = 'TP53', xlab = FALSE, ylab = 'Burden of subclonal SCNAs', legend = 'none', 
 add = "jitter", shape = 'TP53', outlier.shape = NA, alpha = 1.0, add.params=list(size = 0.8)) + 
 stat_compare_means(label.x.npc = 'left', label.y.npc = 'top') + theme(aspect.ratio = 1)


# p3 <- ggplot(data = subset(gli.cn.alt.frac, Integrated_Diagnoses == 'Oligodendroglioma,IDHmut-codel'), aes(TP53, subclo_genome_frac)) +
 # geom_boxplot(aes(shape = TP53), outlier.shape = NA) + geom_jitter(shape=16, position=position_jitter(0.2)) + 
 # xlab('Subtype') + ylab('Burden of subclonal SCNAs') + guides(shape = guide_legend(title= 'Subtype')) + 
 # theme(axis.text.x = element_text(angle=30, vjust=0.5)) + stat_compare_means(label.x.npc = 'left', label.y.npc = 'top')

ggsave(ggarrange(p1, p2, ncol=2, nrow=2, common.legend = TRUE), file='/result/Section2/TP53_subclo_scna_com.pdf')





