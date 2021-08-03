# load data
gli.cn.alt.frac <- readRDS(file='/pub5/xiaoyun/Jobs/J22/CopyNumberClonalityProject/Resource/CuratedData/gli_glod_cn_alt_frac.rds')

# Cox
source('/pub5/xiaoyun/Jobs/J22/EvoClass2.0/Section3/RScripts/Cox.function.R')

SCNAsBurdenSurCox <- function(sur.dat, feacs, subtype, grade){
 sur.dat <- subset(sur.dat, IDH_CODEL_SUBTYPE %in% subtype & histological_grade %in% grade)
 
 cox.pvalues <- lapply(feacs, function(feac){
  
  cox.pvalue <- lapply(c('os', 'dss', 'pfi'), function(x){
   
   if(length(grade) > 1){
    dat <- sur.dat[, c('bcr_patient_barcode', x, paste0(x, '_time'), feac, 'age', 'gender', 'histological_grade')]
   }else{
    dat <- sur.dat[, c('bcr_patient_barcode', x, paste0(x, '_time'), feac, 'age', 'gender')]
   }

   dat[, feac] <- cut_width(x=dat[, feac], width=0.05, labels=FALSE)
   
   dat <- subset(dat, !(is.na(get(x)) | is.na(get(paste0(x, '_time')))))   
   sur.times <- dat[, paste0(x, '_time')]
   sur.status <- dat[, x]
   
   cox.res <- Cox.function(sur.times, sur.status, dat)
   return(cox.res)
  })
  return(cox.pvalue)
 })
 return(cox.pvalues)
}

# Function: 根据多因素cox回归结果绘制森林图
###################################################################

Forestplot <- function(cox.result){
 require(forestplot)
  
 # 数据重组织
 cox.result <- subset(cox.result, variate %in% c("age", "clo_genome_frac") |cox.result[, 2] != " ")
 rownames(cox.result) <- c(cox.result[1:2, 1], cox.result[-c(1, 2), 2])
 cox.result <- cox.result[c('age', 'clo_genome_frac', 'MALE VS FEMALE', 
  'G3 VS G2', 'G4 VS G2'), c('multiv HR (95% CI for HR)',  'multiv p value')]
 
 # 构建需要输入的参数信息
 varibles <- c("Varibles", "Age at dignosis", "Burden of clonal SCNAs", "Gender", "Male VS Female", 
  "Histological grade", "G3 VS G2", "G4 VS G2") 
			   
 info.index <- c(2, 3, 5, 7, 8)
 # 用于绘图的HR信息
 HR <- rep(NA, length(varibles))
 multiCox.HR <- as.character(cox.result[, 1])
 HR[info.index] <- as.numeric(sapply(multiCox.HR, function(x) strsplit(x," \\(")[[1]][1]))
 lowerCI <- rep(NA, length(varibles))
 lowerCI[info.index] <- as.numeric(sapply(multiCox.HR, function(x) strsplit(strsplit(x," \\(")[[1]][2], "-")[[1]][1]))
 upperCI <- rep(NA, length(varibles))
 upperCI[info.index] <- as.numeric(sapply(multiCox.HR, function(x) strsplit(strsplit(x,"-")[[1]][2], "\\)")[[1]][1]))
 # 文字信息
 HR.info <- rep(NA, length(varibles))
 HR.info[info.index] <- paste(round(HR[info.index],2)," (", round(lowerCI[info.index],2), "-", round(upperCI[info.index],2),")",sep="")
 multiCox.P <- as.character(cox.result[,2])
 P.value <- rep(NA, length(varibles))
 P.value[info.index] <- multiCox.P
 tabletext <- cbind(varibles, HR.info, P.value)
 colnames(tabletext) <- c("Varible", "HR(95%CI)", "P")
 
 # 绘制森林图
 

 forest.plot <- forestplot(tabletext[,1:3],
            HR, lowerCI, upperCI, 
            txt_gp = fpTxtGp(ticks = gpar(cex = 0.3), xlab = gpar(cex = 0.5), cex = 0.5),
            graph.pos = 2, # 森林图为位置
            zero = 1, # 设置无效线
            boxsize = .1, # 规定方块的大小一致
            col = fpColors(box = "royalblue",line = "black"), 
            # vertices = TRUE, # 给线增加端点
            clip = c(.1, 8), # 超出范围的区间增加箭头
            xlab ="HR",
            line.margin = unit(5,"mm"),
            lineheight = unit(6,"mm"),
            xticks = (c(0, 0.5, 1:8)),
            graphwidth = unit(30,"mm"))

 return(forest.plot)
}


features <- c('clo_genome_frac')
multi.cox.res <- SCNAsBurdenSurCox(gli.cn.alt.frac, features, 'IDHmut-non-codel', c('G2', 'G3', 'G4'))

setwd('/pub5/xiaoyun/Jobs/J22/CopyNumberClonalityProject/Results/Section1/Results/Survival')
ggsave(plot = ggarrange(plotlist=lapply(multi.cox.res[[1]], Forestplot), ncol=2, nrow=2), 
 filename='non_codel_g234_forestplot.pdf')



