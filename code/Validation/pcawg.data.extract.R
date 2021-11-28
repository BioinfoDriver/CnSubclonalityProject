
file.path <- '/data/OriginalData/PCAWG/CNA'
setwd(file.path)

cna.data <- sapply(list.files(file.path), function(f){
 cna <- read.csv(file=f, header=TRUE, sep='\t', stringsAsFactors=FALSE)
}, simplify = FALSE, USE.NAMES = TRUE)
names(cna.data) <- sapply(strsplit(names(cna.data), split='\\.'), function(x) x[1])


setwd('/data/OriginalData/PCAWG')
pcawg.puri.ploi <- read.table(file='consensus.20170218.purity.ploidy.txt', sep='\t', header=T, stringsAsFactors=F)
pcawg.sample <- read.table(file='pcawg_sample_sheet.tsv', sep='\t', header=T, stringsAsFactors=F)
pcawg.sample <- subset(pcawg.sample, aliquot_id %in% names(cna.data))

library('readxl')
pcawg.cli.data <- read_xlsx(path='pcawg_donor_clinical_August2016_v9 (2).xlsx', sheet=1)
pcawg.cli.data <- merge(pcawg.cli.data, pcawg.sample[, c('aliquot_id', 'icgc_donor_id', 'dcc_specimen_type')], by='icgc_donor_id')

# filter
glioma.cli.data <- subset(pcawg.cli.data, project_code %in% c('GBM-US', 'LGG-US') & dcc_specimen_type == 'Primary tumour - solid tissue')
glioma.cna.data <- cna.data[glioma.cli.data$aliquot_id]
glioma.puri.ploi <- subset(pcawg.puri.ploi, samplename %in% glioma.cli.data$aliquot_id)

glioma.puri.ploi$submitted_donor_id <- 
 glioma.cli.data$submitted_donor_id[match(glioma.puri.ploi$samplename, glioma.cli.data$aliquot_id)]
names(glioma.cna.data) <- glioma.cli.data$submitted_donor_id[match(names(glioma.cna.data), glioma.cli.data$aliquot_id)]


saveRDS(glioma.puri.ploi, file='/data/pcawg.glioma.puri.ploi.rds')
saveRDS(glioma.cna.data, file='/data/pcawg.glioma.cna.data.rds')

