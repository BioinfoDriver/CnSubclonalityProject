
#####################Survival
rtks <- read.csv(file='/data/RTKs.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)

gold.set <- readRDS('/data/gold_set.rds')
gene.het <- readRDS(file='/data/gene.het.rds')
gene.het <- gene.het[, paste(gold.set, '01', sep = '-')]

tcga.cli.data <- readRDS('/data/tcga_glioma_cli_mol.rds')
tcga.cli.data <- tcga.cli.data[gold.set, ]
rownames(tcga.cli.data) <- paste0(rownames(tcga.cli.data), '-01', sep='')

# Survival plot
plot.surv <- function(clinical.data, upper.time = NULL, xscale = 1, xlab = "Time", median.time = TRUE, 
                      surv.median.line = "none", HR = FALSE, risk.table = TRUE, pval = TRUE, 
                      conf.int = FALSE, main = NULL, ylab = "Survival probability", break.x = 30.4375) {
  
  #载入相关R包
  require(survival)
  require(survminer)
  require(RColorBrewer)
  require(gridExtra)
  
	#确定事件类型和时间的单位
  # survival.event <- survival.event[1];
  # unit.xlabel <- unit.xlabel[1];
  
  #如果设置upper.time，则去除生存时间超过upper.time的样本
  if (!is.null(upper.time)) clinical.data <- clinical.data[clinical.data$time <= upper.time,]
  
  # #选择日期格式 
  # xSL <- data.frame(xScale=c(1,7,30,365.25),xLab=c("Days","Weeks","Months","Years"), stringsAsFactors=FALSE)
  # switch(unit.xlabel, year={xScale <- 365.25;}, month={xScale <- 30;}, week={xScale <- 7;}, day={xScale <- 1})
  # xLab <- xSL[which(xSL[,1]==xScale),2];
  
  #构造颜色
	if (!is.factor(clinical.data$sample.label)) 
    clinical.data$sample.label <- as.factor(clinical.data$sample.label)
  
  t.name <- levels(clinical.data$sample.label)

	if (length(t.name) > 6) stop("样本分组>6，超过函数接受范围")
    colors <- c("#808080","#EA4335","#4285F4","#FBBC05","#34A853","#000000") # 顺序：灰，红，蓝，黄，绿，黑
    t.col <- colors[1:length(t.name)]
    
  # 构造生存对象
  km.curves <- survfit(Surv(time, event)~sample.label, data=clinical.data)
    
    # 计算HR值和95%CI
	if (length(t.name) == 2) {
	  if (HR) {
        cox.obj <- coxph(Surv(time, event)~sample.label, data=clinical.data)
        tmp <- summary(cox.obj)
        HRs <- round(tmp$coefficients[ ,2], digits = 2)
        HR.confint.lower <- round(tmp$conf.int[,"lower .95"], 2)
        HR.confint.upper <- round(tmp$conf.int[,"upper .95"], 2)
        HRs <- paste0(HRs, " (", HR.confint.lower, "-", HR.confint.upper, ")")	  
	  }
  }

  # 构造生存图像中图例显示文字
  legend.content <- substr(names(km.curves$strata), start = 14, stop = 1000)
  
	# x轴刻度单位转换
  if (is.numeric(xscale) | (xscale %in% c("d_m", "d_y", "m_d", "m_y", "y_d", "y_m"))) {
    xscale = xscale
  } else {
    stop('xscale should be numeric or one of c("d_m", "d_y", "m_d", "m_y", "y_d", "y_m").')
  }
    
	# 隐函数：转换生存时间单位
  .format_xticklabels <- function(labels, xscale){
    # 1 year = 365.25 days
    # 1 month = 365.25/12 = 30.4375 days
    if (is.numeric(xscale)) xtrans <- 1/xscale
    else
      xtrans <- switch(xscale,
                       d_m = 12/365.25,
                       d_y = 1/365.25,
                       m_d = 365.25/12,
                       m_y = 1/12,
                       y_d = 365.25,
                       y_m = 12,
                       1
      )
    round(labels*xtrans,2)
  }

  # 在图中添加中位生存时间及其95%CI,放在副标题位置
  subtitle <- NULL
  if (median.time) {
    if (is.numeric(xscale)) {
      median.km.obj = km.curves
    } else if (xscale %in% c("d_m", "d_y", "m_d", "m_y", "y_d", "y_m")) {
      clinical.data$time <- .format_xticklabels(labels = clinical.data$time, xscale = xscale)
      median.km.obj <- survfit(Surv(time, event)~sample.label, data=clinical.data)
    }
    survival.time.info <- NULL
    survival.time.info <- rbind(survival.time.info, summary(median.km.obj)$table)
    median.survival <- round(survival.time.info[!duplicated(survival.time.info[,7:9]),7:9], digits = 2) # 注意：这里取得的置信区间上界可能为NA
    if (length(levels(clinical.data$sample.label)) == 1) {
      tmp1 <- levels(clinical.data$sample.label)
    } else {
      tmp1 <- do.call(rbind,strsplit(rownames(summary(median.km.obj)$table), split = "="))[,2]
    }
    tmp2 <- paste(median.survival[,1], "(", median.survival[,2], "-", median.survival[,3], ")")
    subtitle <- paste(tmp1, tmp2, sep = ":", collapse = "\n")
  }
    
  # ggsurvplot绘制生存图像
	ggsurv <- ggsurvplot(km.curves,               # survfit object with calculated statistics.
						 data = clinical.data,             # data used to fit survival curves.
						 palette = t.col,
						 
						 #图的主题构架
						 risk.table = risk.table,       # show risk table.
						 pval = pval,             # show p-value of log-rank test.
						 surv.median.line = surv.median.line,  # add the median survival pointer.
						 title = main,     #主标题名字
						 subtitle = subtitle, #副标题
						 font.main = 15,       #主标题字体大小              
						 xlab = xlab,   # customize X axis label.
						 ylab = ylab,   # customize Y axis label
						 xscale = xscale,
						 break.x.by = break.x, # m=30.4375, y=365.25


						 #图例设置
						 legend.title = "", #图例标题，一般不用，设为空
						 legend.labs = legend.content, #图例文字描述
						 legend = c(0.8,0.9), #图例的位置，取值在【0,1】之间
						 font.legend = 9,     #图例字体大小
						 
						 #risk table设置
						 tables.theme = theme_cleantable(),#table主题
						 risk.table.title = "No. at risk:",#table标题
						 risk.table.y.text.col = T, # 使用颜色代替Y轴文字
						 risk.table.y.text = FALSE, # Y轴不使用文字注释
						 tables.height = 0.15,      # table的高度
						 risk.table.fontsize = 3    # risk table内文字的大小
						)
	# HR标注的位置限定
  if (length(t.name) == 2) {
	  if (HR) 
      ggsurv$plot <- ggsurv$plot + ggplot2::annotate("text", x = max(km.curves$time)/13, 
                                                     y = 0.15, size = 5, label = paste("HR=", HRs))
  } 
	#图的标题居中
  ggsurv$plot <- ggsurv$plot + theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(size = 10), 
								                     plot.margin = unit(c(5.5, 5.5, 5.5, 50), "points"))
  #表的标题
  ggsurv$table <- ggsurv$table + theme(plot.title = element_text(hjust = -0.04),
                                       plot.margin = unit(c(5.5, 5.5, 5.5, 50), "points"))

  # 判断分类的类数，如果只有两类，就不必计算两两之间的log rank p值
  if(length(t.name) > 2) {
  # 计算pairwise的log rank的p值
  res <- pairwise_survdiff(Surv(time, event)~sample.label, data=clinical.data, p.adjust.method='none');
  pairwise.pvalue <- res$p.value;
  pairwise.pvalue[which(pairwise.pvalue < 0.001)] <- "<0.001";
  pairwise.pvalue[is.na(pairwise.pvalue)] <- "-"
  
  pairwise.pvalue <- apply(pairwise.pvalue, 2, substr, start = 1, stop = 6)
  # 添加表格
  tt <- ttheme_minimal(core = list(fg_params = list(col = "black"),bg_params = list(fill = NA, col = "black")),
                       colhead = list(fg_params = list(col = NA),bg_params = list(fill = t.col, col = "black")),
                       rowhead = list(fg_params = list(col = NA, hjust = 1),bg_params = list(fill = c("white",t.col[-1]), col = "black"))
                      )
  pairwise.table <- tableGrob(pairwise.pvalue, theme = tt)
  ggsurv <- ggarrange(ggarrange(ggsurv$plot, ggsurv$table, nrow=2, heights=c(2,0.5)),
                      pairwise.table, nrow=2, heights = c(2,0.5),
                      labels = c("","pairwise comparisons"),
                      hjust = 0, font.label = list(size = 15, face = "plain"))
  }

  ggsurv	
}


RtkLogRankSurAnalysis <- function(het.mat, rtk, cli.data, group.index, subtype, grade, up.time, file.name, break.x){
 
 rtk.het <- het.mat[rtk, ]

 rtk.clo.sam <- colnames(rtk.het[, colSums(rtk.het == 1) >= 1])
 rtk.sub.sam <- colnames(rtk.het[, colSums(rtk.het == 2) >= 1])
 
 if(group.index=='subclonaldominant'){
  rtk.clo.sam <- setdiff(rtk.clo.sam, rtk.sub.sam)

 }else if(group.index=='clonaldominant'){
  rtk.sub.sam <- setdiff(rtk.sub.sam, rtk.clo.sam)
 
 }else{
  com.sam <- intersect(rtk.clo.sam, rtk.sub.sam)
  rtk.clo.sam <- setdiff(rtk.clo.sam, com.sam)
  rtk.sub.sam <- setdiff(rtk.sub.sam, com.sam)
  
 }

 
 non.rtk.alt.sam <- setdiff(colnames(het.mat), c(rtk.clo.sam, rtk.sub.sam))
 
 # rtk.clo.mat <- data.frame(subtype = rep(c('rtk clonal amp', 'rtk subclonal amp', 'rtk wt'), 
  # times=c(length(rtk.clo.sam), length(rtk.sub.sam), length(non.rtk.alt.sam))))
 # rownames(rtk.clo.mat) <- c(rtk.clo.sam, rtk.sub.sam, non.rtk.alt.sam)
 
 rtk.clo.mat <- data.frame(subtype = rep(c('rtk clonal amp', 'rtk subclonal amp'), 
  times=c(length(rtk.clo.sam), length(rtk.sub.sam))))
 rownames(rtk.clo.mat) <- c(rtk.clo.sam, rtk.sub.sam) 
 
 cli.data <- subset(cli.data, IDH_CODEL_SUBTYPE %in% subtype & histological_grade %in% grade)
 cli.data <- merge(cli.data, rtk.clo.mat, by = 'row.names')
 
 # OS, DSS PFI
 sur.plot.list <- lapply(c('os', 'dss', 'pfi'), function(surv.type){
  
  surv.data <- cli.data[, c('Row.names', surv.type, paste0(surv.type, '_time'), 'subtype')]
  colnames(surv.data) <- c('Patient_ID', 'event', 'time', 'sample.label')
  surv.data <- subset(surv.data, !(is.na(event) | is.na(time) | is.na(sample.label)))

  surv.plot <- try(plot.surv(clinical.data = surv.data, xscale = 'd_m', xlab = 'Time in months', 
  risk.table = TRUE, median.time = TRUE, surv.median.line = 'hv', upper.time = up.time, break.x = break.x)) 
  
  return(surv.plot)
 })

 # plot 
 pdf(file.name)
  for(i in seq(length(sur.plot.list))){
   print(sur.plot.list[i])
  }
 dev.off()
 
}

setwd('/result/Section4')

RtkLogRankSurAnalysis(gene.het, rtks$Approved.symbol, tcga.cli.data, 'exclude', 
 c('IDHwt'), c('G4'), 1096, 'rtk_sur_idhwt_gbm.pdf', 30.4375*2)
RtkLogRankSurAnalysis(gene.het, rtks$Approved.symbol, tcga.cli.data, 'exclude', 
 c('IDHmut-non-codel'), c('G2', 'G3'), 3653, 'rtk_sur_idhwt_lgg.pdf', 30.4375*12)


# Cox
source('/code/Function/Cox.function.R')
RtkMulCoxSurAnalysis <- function(het.mat, rtk, group.index, subtype, grade, cli.data, up.time){
 
 rtk.het <- het.mat[rtk, ]

 rtk.clo.sam <- colnames(rtk.het[, colSums(rtk.het == 1) >= 1])
 rtk.sub.sam <- colnames(rtk.het[, colSums(rtk.het == 2) >= 1])
 
 if(group.index=='subclonaldominant'){
  rtk.clo.sam <- setdiff(rtk.clo.sam, rtk.sub.sam)

 }else if(group.index=='clonaldominant'){
  rtk.sub.sam <- setdiff(rtk.sub.sam, rtk.clo.sam)
 
 }else{
  com.sam <- intersect(rtk.clo.sam, rtk.sub.sam)
  rtk.clo.sam <- setdiff(rtk.clo.sam, com.sam)
  rtk.sub.sam <- setdiff(rtk.sub.sam, com.sam)
  
 }
 
 # non.rtk.alt.sam <- setdiff(colnames(het.mat), c(rtk.clo.sam, rtk.sub.sam))
 
 # rtk.clo.mat <- data.frame(group = rep(c('rtk clonal amp', 'rtk subclonal amp', 'rtk wt'), 
  # times=c(length(rtk.clo.sam), length(rtk.sub.sam), length(non.rtk.alt.sam))))
 # rownames(rtk.clo.mat) <- c(rtk.clo.sam, rtk.sub.sam, non.rtk.alt.sam)

 rtk.clo.mat <- data.frame(group = rep(c('rtk clonal amp', 'rtk subclonal amp'), 
  times=c(length(rtk.clo.sam), length(rtk.sub.sam))))
 rownames(rtk.clo.mat) <- c(rtk.clo.sam, rtk.sub.sam)

 cli.data <- subset(cli.data, IDH_CODEL_SUBTYPE %in% subtype & histological_grade %in% grade)
 cli.data <- merge(cli.data, rtk.clo.mat, by = 'row.names')
 
 # OS, DSS PFI
 mul.cox.res <- lapply(c('os', 'dss', 'pfi'), function(surv.type){
  
  surv.data <- cli.data[, c('Row.names', surv.type, paste0(surv.type, '_time'), 
   'age', 'gender', 'histological_grade', 'IDH_CODEL_SUBTYPE', 'group')]
  colnames(surv.data) <- c('Patient_ID', 'event', 'time', 'age', 'gender', 'grade', 'subtype', 'group')
  
  surv.data <- subset(surv.data, !is.na(grade))
  if (!is.null(up.time)) surv.data <- surv.data[surv.data$time <= up.time,]
  
  
  print(dim(surv.data))
  if(length(unique(surv.data$subtype)) > 1 & length(unique(surv.data$grade)) > 1){
   cox.res <- Cox.function(surv.data$time, surv.data$event, clinical.data=surv.data, clinical.variate=c(4:8))
  }
  
  if(length(unique(surv.data$subtype)) == 1 & length(unique(surv.data$grade)) > 1){
   cox.res <- Cox.function(surv.data$time, surv.data$event, clinical.data=surv.data, clinical.variate=c(4, 5, 6, 8))
  }
  
  if(length(unique(surv.data$subtype)) > 1 & length(unique(surv.data$grade)) == 1){
   cox.res <- Cox.function(surv.data$time, surv.data$event, clinical.data=surv.data, clinical.variate=c(4, 5, 7, 8))
  }
  
  if(length(unique(surv.data$subtype)) == 1 & length(unique(surv.data$grade)) == 1){
   cox.res <- Cox.function(surv.data$time, surv.data$event, clinical.data=surv.data, clinical.variate=c(4, 5, 8))
  }  

  return(cox.res)
 })

 return(mul.cox.res)
}


RtkMulCoxSurAnalysis(gene.het, rtks$Approved.symbol, 'exclude', 
 c('IDHwt'), c('G4'), tcga.cli.data, 1096)

RtkMulCoxSurAnalysis(gene.het, rtks$Approved.symbol, 'exclude', 
 c('IDHwt'), c('G2', 'G3'), tcga.cli.data, 3653)


