
# ITH of SCNAs
glass.cna.ith <- readRDS('/data/glass_seqz_cna_ith.rds')

# glass.cna.ith$purity_a <- glass.purity$cellularity[match(paste(glass.cna.ith$tumor_barcode_a, 
 # glass.cna.ith$seq_method, sep='-'), glass.purity$tumor_barcode)]

# glass.cna.ith$purity_b <- glass.purity$cellularity[match(paste(glass.cna.ith$tumor_barcode_b, 
 # glass.cna.ith$seq_method, sep='-'), glass.purity$tumor_barcode)]
 
# Statistics
# table(glass.cna.ith[, c('treatment_tmz', 'treatment_radiotherapy')])
             # treatment_radiotherapy
# treatment_tmz FALSE NULL TRUE
        # FALSE    30    0   28
        # NULL     10   21   10
        # TRUE     26    0  118


# 2021 WHO subtype
glass.cna.ith <- subset(glass.cna.ith, !is.na(Integrated_Diagnoses)) # 219

# treatment-naive primary and matched post-treatment first-recurrence tumour samples
glass.cna.ith <- subset(glass.cna.ith, substr(tumor_barcode_a, 14, 15) == 'TP' & 
  substr(tumor_barcode_b, 14, 15) == 'R1') #190

# without treatment information
glass.cna.ith <- subset(glass.cna.ith, treatment_tmz != 'NULL' | treatment_radiotherapy != 'NULL') #188

# Statistics
# table(glass.cna.ith[, c('treatment_tmz', 'treatment_radiotherapy')])
             # treatment_radiotherapy
# treatment_tmz FALSE TRUE
        # FALSE    26   24
        # NULL      8    9
        # TRUE     21  100


# table(glass.cna.ith[, c('treatment_tmz', 'Integrated_Diagnoses')])
             # Integrated_Diagnoses
# treatment_tmz Astrocytoma,IDHmut Glioblastoma,IDHwt
        # FALSE                 26                 12
        # NULL                   7                  6
        # TRUE                  22                 91
             # Integrated_Diagnoses
# treatment_tmz Oligodendroglioma,IDHmut-codel
        # FALSE                             12
        # NULL                               4
        # TRUE                               8
# table(glass.cna.ith[, c('treatment_radiotherapy', 'Integrated_Diagnoses')])
                      # Integrated_Diagnoses
# treatment_radiotherapy Astrocytoma,IDHmut Glioblastoma,IDHwt
                 # FALSE                 33                  6
                 # TRUE                  22                103
                      # Integrated_Diagnoses
# treatment_radiotherapy Oligodendroglioma,IDHmut-codel
                 # FALSE                             16
                 # TRUE                               8



glass.cna.ith$cna.bur.r <- glass.cna.ith$clonal.cna.bur.r + glass.cna.ith$subclonal.cna.bur.r


features <- c("cna.bur.r", "subclo.cna.perc.r", "clonal.cna.bur.r" , "subclonal.cna.bur.r")			 

ScnaTreatAsso <- function(feacs, subtype, treatment){
	if(subtype!='all'){
		cna.ith <- subset(glass.cna.ith, Integrated_Diagnoses == subtype)
		cna.ith <- subset(cna.ith, get(treatment) != 'NULL')
	}else{
		cna.ith <- subset(glass.cna.ith, get(treatment) != 'NULL')
	}
	

	p.values <- sapply(feacs, function(feacture){
		p.value <- wilcox.test(get(feacture) ~ get(treatment), data=cna.ith)$p.value
		return(p.value)
	})
	
	return(p.values)
}


# ScnaTreatAsso(features, 'all', 'treatment_radiotherapy')
# ScnaTreatAsso(features, 'all', 'treatment_tmz')

ScnaTreatAsso(features, 'Astrocytoma,IDHmut', 'treatment_radiotherapy')
ScnaTreatAsso(features, 'Astrocytoma,IDHmut', 'treatment_tmz')

ScnaTreatAsso(features, 'Glioblastoma,IDHwt', 'treatment_radiotherapy')
ScnaTreatAsso(features, 'Glioblastoma,IDHwt', 'treatment_tmz')

ScnaTreatAsso(features, 'Oligodendroglioma,IDHmut-codel', 'treatment_radiotherapy')
ScnaTreatAsso(features, 'Oligodendroglioma,IDHmut-codel', 'treatment_tmz')

########### plot
setwd('/result/Section2/')
library('ggpubr')
scna.plot <- ggboxplot(subset(glass.cna.ith, treatment_tmz != 'NULL'), 'treatment_tmz', 'cna.bur.r', 
 xlab = FALSE, ylab = 'SCNA burden', facet.by='Integrated_Diagnoses', add = "jitter") + 
 stat_compare_means(label.x.npc = 'left', label.y.npc = 'middle')

subclo.perc.plot <- ggboxplot(subset(glass.cna.ith, treatment_tmz != 'NULL'), 'treatment_tmz', 'subclo.cna.perc.r', 
 xlab = FALSE, ylab = 'Subclonal SCNA percentage', facet.by='Integrated_Diagnoses', add = "jitter") + 
 stat_compare_means(label.x.npc = 'left', label.y.npc = 'middle')

clo.bur.plot <- ggboxplot(subset(glass.cna.ith, treatment_tmz != 'NULL'), 'treatment_tmz', 'clonal.cna.bur.r', 
 xlab = FALSE, ylab = 'Clonal SCNA burden', facet.by='Integrated_Diagnoses', add = "jitter") + 
 stat_compare_means(label.x.npc = 'left', label.y.npc = 'middle')

subclo.bur.plot <- ggboxplot(subset(glass.cna.ith, treatment_tmz != 'NULL'), 'treatment_tmz', 'subclonal.cna.bur.r', 
 xlab = FALSE, ylab = 'Subclonal SCNA burden', facet.by='Integrated_Diagnoses', add = "jitter") + 
 stat_compare_means(label.x.npc = 'left', label.y.npc = 'middle')


source('/code/Function/multiplot.r')
pdf(file='r_tmz_scna_asso_by_subtype.pdf')
 multiplot(plotlist=list(scna.plot, subclo.perc.plot, clo.bur.plot, subclo.bur.plot), cols=1)
dev.off()


########### plot
scna.plot <- ggboxplot(subset(glass.cna.ith, treatment_radiotherapy != 'NULL'), 'treatment_radiotherapy', 'cna.bur.r', 
 xlab = FALSE, ylab = 'SCNA burden', add = "jitter") + 
 stat_compare_means(label.x.npc = 'left', label.y.npc = 'top')

subclo.perc.plot <- ggboxplot(subset(glass.cna.ith, treatment_radiotherapy != 'NULL'), 'treatment_radiotherapy', 'subclo.cna.perc.r', 
 xlab = FALSE, ylab = 'Subclonal SCNA percentage', add = "jitter") + 
 stat_compare_means(label.x.npc = 'left', label.y.npc = 'top')

clo.bur.plot <- ggboxplot(subset(glass.cna.ith, treatment_radiotherapy != 'NULL'), 'treatment_radiotherapy', 'clonal.cna.bur.r', 
 xlab = FALSE, ylab = 'Clonal SCNA burden', add = "jitter") + 
 stat_compare_means(label.x.npc = 'left', label.y.npc = 'top')

subclo.bur.plot <- ggboxplot(subset(glass.cna.ith, treatment_radiotherapy != 'NULL'), 'treatment_radiotherapy', 'subclonal.cna.bur.r', 
 xlab = FALSE, ylab = 'Subclonal SCNA burden', add = "jitter") + 
 stat_compare_means(label.x.npc = 'left', label.y.npc = 'top')

ggsave(ggarrange(scna.plot, subclo.perc.plot, clo.bur.plot, subclo.bur.plot, 
 nrow = 2, ncol = 2, common.legend = TRUE), file='r_radio_scna_asso.pdf')
