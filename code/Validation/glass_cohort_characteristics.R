
glass.cna.ith <- readRDS(file='/data/glass_seqz_cna_ith.rds')
glass.cna.ith <- subset(glass.cna.ith, !is.na(Integrated_Diagnoses))
glass.cna.ith[glass.cna.ith=='NULL'] <- NA

setwd('/data/OriginalData/syn17038081')
clinical.data <- read.csv(file='clinical_cases.tsv', header=TRUE, sep='\t', stringsAsFactors=FALSE)
clinical.data <- merge(clinical.data, glass.cna.ith, by='case_barcode')

clinical.data$surgery_location[!(is.na(clinical.data$surgery_location) | 
 clinical.data$surgery_location %in% c('Temporal lobe', 'Parietal lobe', 'Frontal lobe'))] <- 'Other'


CliSamCountStat <- function(cli.features, cli.data){
  
  sam.count <- dim(cli.data)[1]
  cli.stat.result <- lapply(cli.features, function(sin.cli.fea){
	
	# Distinguish the variable as continuous or categorical
	cli.varibale <- cli.data[, sin.cli.fea]
    if(all(is.na(cli.varibale))){
      value <- NA
    }else{
      cli.varibale <- na.omit(cli.varibale)
	  cli.varibale <- as.character(cli.varibale)
      value <- all(!is.na(as.numeric(cli.varibale))) 
    }
	
	if(is.na(value)){
	cli.df <- as.matrix(data.frame(cli.feature = sin.cli.fea, cli.feature.stat = 'Unknown'))
	
	}else if(value == TRUE){
	# continuous, especially for age
		cli.varibale <- as.numeric(cli.varibale)
		cli.df <- as.matrix(data.frame(cli.feature = "Age", cli.feature.stat = paste(median(cli.varibale), '(',min(cli.varibale), '-', max(cli.varibale), ')', sep='')))
	
	}else if(value == FALSE){
	# categorical
	cli.table <- table(cli.data[, sin.cli.fea])
    cli.table.frac <- cli.table/sam.count
    unknown.count <- sam.count - sum(cli.table)
    unknown.frac <- 1-sum(cli.table.frac)
    cli.feature.count <- c(as.numeric(cli.table), unknown.count)
    cli.feature.frac <- c(round(cli.table.frac, 2)*100, round(unknown.frac, 2)*100)
    
    cli.fea.count.frac <- paste(cli.feature.count, '(', cli.feature.frac, ')', sep = '')
    
    cli.df <- as.matrix(data.frame(cli.feature = c(sin.cli.fea, names(cli.table), 'Unknown'), cli.feature.stat = c(sin.cli.fea, cli.fea.count.frac), stringsAsFactors = FALSE))

	}
	t(cli.df)
  })
  
  cli.stat.result
}

library('tidyverse')
cli.features <- c("case_sex", "case_age_diagnosis_years", "grade", 
 'surgery_extent_of_resection', 'surgery_laterality', 'surgery_location', "mgmt_methylation")

all.statistic <- CliSamCountStat(cli.features, clinical.data)


# Statistical analysis of clinical molecular characteristics on each subgroups
subtype.name <- unique(clinical.data$Integrated_Diagnoses)
subtype.statistic <- lapply(subtype.name, function(per.type){
	per.subtype.characters <- subset(clinical.data, Integrated_Diagnoses == per.type)
	per.subtype.statistic <- CliSamCountStat(cli.features, per.subtype.characters)
	per.subtype.statistic
})
names(subtype.statistic) <- subtype.name
