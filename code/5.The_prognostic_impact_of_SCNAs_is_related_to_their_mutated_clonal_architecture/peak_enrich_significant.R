

gold.set <- readRDS('/data/gold_set.rds')
wgistic.het <- readRDS('/data/wgistic_het.rds')
gli.cli.mol.data <- readRDS(file='/data/tcga_glioma_cli_mol.rds')
subtype.peaks <- readRDS(file='/data/high_fre_alt_peaks.rds')


gli.cli.mol.data <- gli.cli.mol.data[gold.set, ]
rownames(gli.cli.mol.data) <- paste(rownames(gli.cli.mol.data), '01', sep = '-')
wgistic.het <- wgistic.het[, paste(gold.set, '01', sep = '-')]

# samples in each subtype
subtype.samples <- sapply(c('IDHmut-codel', 'IDHmut-non-codel', 'IDHwt'), function(subtype){
 
 return(rownames(subset(gli.cli.mol.data, IDH_CODEL_SUBTYPE==subtype)))
 
}, USE.NAMES=TRUE)

# percentage of clonal alteration
peak.clo.alt.pro <- sapply(c('IDHmut-codel', 'IDHmut-non-codel', 'IDHwt'), function(subtype){
	
	het.mat <- wgistic.het[subtype.peaks[[subtype]], subtype.samples[[subtype]]]
	
	clo.pro <- rowSums(abs(het.mat)==1)/rowSums(het.mat != 0)
	
	return(clo.pro)
}, USE.NAMES=TRUE)

# Permutation
ran.peak.clo.alt.pro <- sapply(c('IDHmut-codel', 'IDHmut-non-codel', 'IDHwt'), function(subtype){
	
	het.mat <- wgistic.het[subtype.peaks[[subtype]], subtype.samples[[subtype]]]
		
	ran.clo.pro <- sapply(rownames(het.mat), function(i){
		
		clo.pro <- c()
		for(j in seq(100000)){
			
			ran.het <- sample(x=1:2, size=sum(het.mat[i, ] !=0), replace=TRUE)
			clo.pro <- c(clo.pro, sum(ran.het==1)/length(ran.het))
		}
			
		return(clo.pro)
	}, simplify = FALSE)
	
	return(ran.clo.pro)
}, USE.NAMES=TRUE)


p.values.list <- sapply(c('IDHmut-codel', 'IDHmut-non-codel', 'IDHwt'), function(subtype){
	peak.clo.pro <- peak.clo.alt.pro[[subtype]]
	
	p.values <- sapply(names(peak.clo.pro), function(peak){
	
		ran.peak.clo.pro <- ran.peak.clo.alt.pro[[subtype]][[peak]]
		p.value <- sum(peak.clo.pro[peak]>=ran.peak.clo.pro)/length(ran.peak.clo.pro)
		
		return(p.value)
	})
	return(p.values)
})


# The probability of clonal alteration 
clo.alt.pvalue <- sapply(c('IDHmut-codel', 'IDHmut-non-codel', 'IDHwt'), function(subtype){
	
	het.mat <- wgistic.het[subtype.peaks[[subtype]], subtype.samples[[subtype]]]
	p.values <- sapply(rownames(het.mat), function(i){
			
		clo.n <- sum(abs(het.mat[i, ])==1)
		alt.n <- sum(abs(het.mat[i, ])!=0)
		p.value <- binom.test(x=clo.n, n=alt.n, p = 0.5, alternative = "two.sided", conf.level = 0.95)$p.value
		   
		return(p.value)
	})
	
	return(p.values)
	
})


clo.alt.pvalue.adjust <- sapply(clo.alt.pvalue, function(x) p.adjust(p=x, method = 'fdr'))


