
# CN concordance with WGS

# SCNA-SNP
glioma.puri.ploi <- readRDS(file='/data/pcawg.glioma.puri.ploi.rds')
glioma.cna.data <- readRDS(file='/data/pcawg.glioma.cna.data.rds')

donor_id <- rep(names(glioma.cna.data), lapply(glioma.cna.data, nrow))
glioma.cna.data <- do.call(rbind, glioma.cna.data)
glioma.cna.data$donor_id <- donor_id
glioma.cna.data$ploidy <- glioma.puri.ploi$ploidy[match(glioma.cna.data$donor_id, glioma.puri.ploi$submitted_donor_id)]

glioma.cna.data <- subset(glioma.cna.data, !is.na(absolute_broad_major_cn))
glioma.cna.data$total_cn <- glioma.cna.data$absolute_broad_major_cn + glioma.cna.data$absolute_broad_minor_cn

glioma.cna.data <- subset(glioma.cna.data, total_cn != round(ploidy) & chromosome != 'X')
glioma.cna.data$chromosome  <- paste0("chr", glioma.cna.data$chromosome)
glioma.cna.data$donor_id <- paste0(glioma.cna.data$donor_id, '-01')


glioma.cna.data <- glioma.cna.data[, c('donor_id', 'chromosome', 'start', 'end', 'total_cn')]
colnames(glioma.cna.data) <- c("Sample", "Chromosome", "Start", "End", "Modal_Total_CN")

library('GenomicRanges')
pcawg.call.grl <- makeGRangesListFromDataFrame(glioma.cna.data, 
        split.field="Sample", keep.extra.columns=TRUE) 

# SCNA-SNP
abs.call.grl <- readRDS(file='/data/tcga_gli_seg_call_grl.rds')


library(RaggedExperiment)
isect.samples <- intersect(names(pcawg.call.grl), names(abs.call.grl))
pcawg.call.grl<- pcawg.call.grl[isect.samples]
pcawg.calls.ra <- RaggedExperiment(pcawg.call.grl) 

abs.call.grl<- abs.call.grl[isect.samples]
abs.calls.ra <- RaggedExperiment(abs.call.grl) 

# GISTIC2 region
gistic.peak <- readRDS(file='/data/OriginalData/gbmlgg_gistic_peak.rds')
gistic <- rowRanges(gistic.peak)


Weightedmean <- function(score, range, qrange){
    w <- width(range)
    s <- sum(score * w) / sum(w)
    return(round(s))
}

pcawg.rassay <- qreduceAssay(pcawg.calls.ra, query=gistic, simplifyReduce=Weightedmean, background=2)
abs.rassay <- qreduceAssay(abs.calls.ra, query=gistic, simplifyReduce=Weightedmean, background=2)


# weighted mean
fract.equal <- vapply(seq_along(isect.samples), 
 function(i) mean(pcawg.rassay[,i] == abs.rassay[,i]), numeric(1))

fract.diff <- vapply(seq_along(isect.samples), 
 function(i) mean((pcawg.rassay[,i] == abs.rassay[,i]) | 
				 (pcawg.rassay[,i] == abs.rassay[,i] + 1) |
				 (pcawg.rassay[,i] == abs.rassay[,i] - 1)), numeric(1))

plot.data <- data.frame(type=rep(c('Equality', 'Difference'), each=length(fract.equal)), value=c(fract.equal, fract.diff))
p <- ggboxplot(data=plot.data, x='type', y='value', add='jitter', add.params = list(size=0.9,  alpha=0.5), 
 xlab = FALSE, ylab="Fraction of concordant", outlier.shape = NA, fill=c("#0072B2", "#D55E00"), ylim=c(0.70, 1.0))
p <- ggarrange(p, nrow=2, ncol=2)
ggsave(p, file='/result/Validation/CNA_concordance.pdf')


