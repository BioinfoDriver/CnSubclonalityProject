
#  calculate the C-index for a cox model
#'@description 计算生存模型的C-index或者两个模型比较的p值
#'@param clinical.data: 临床变量信息，至少包含model中的变量信息
#'@param time: 数值型向量，患者对应的生存时间
#'@param event: 患者对应的生存状态，通常0=alive, 1=dead
#'@param models：列表，每一个元素是一个字符型向量，包含一个model所有变量对应的列名,至少包含两个元素
#'@param diff: 逻辑值，diff = FALSE不计算两个模型比较的p值，只返回对应的C-index值；diff = TRUE返回模型对应的C-index值及模型比较的p值
#'@return 返回一个data.frame，包含C-index值，置信区间，比较的p值（当diff = T时）
cIndex <- function(clinical.data, time, event, models, diff = TRUE){
  
  options(stringsAsFactors = FALSE)
  
  # load survival, Hmisc package
  suppressPackageStartupMessages(require(survival))
  suppressPackageStartupMessages(require(Hmisc))
  
  ##-------------------------------
  # CoxphObject：create a coxph object
  CoxphObject <- function(data, time, event, variables){
    formula <- as.formula(paste("Surv(time, event)~", paste(variables, collapse=" + ")))
    coxph.object <-coxph(formula, data = data)
    
    return(coxph.object)
  }
  
  ##-------------------------------
  # 输出的Cindex和置信区间
  c.ci <- function(rcorrobj){
    CIndex <- round(rcorrobj['C Index'], digits = 4)
    se     <- rcorrobj['S.D.']/2
    Lower <- round(CIndex - 1.96*se, digits = 4)
    Upper <- round(CIndex + 1.96*se, digits = 4)
    result <- c(CIndex, Lower, Upper)
    names(result) <- c("C-Index", "Lower", "Upper")
    
    return(result)
  }
  
  # 计算每个模型的C-index值和置信区间
  coxph.list <- lapply(models, function(x){CoxphObject(data = clinical.data, time = time, event = event, variables = x)})
  pred.models.coxph <- lapply(coxph.list, function(x){predict(x, type = "risk")})
  models.result <- lapply(pred.models.coxph, function(x){rcorr.cens(-x, Surv(time = time, event = event))})
  
  models.filter.result <- lapply(models.result, function(x){c.ci(rcorrobj = x)})
  
  # 是否进行C-index1的比较
  if (diff == FALSE) {
    result <- do.call(rbind, models.filter.result)
    conf.interval <- paste(result[, 'Lower'], result[, 'Upper'], sep = "-")
    result <- data.frame(result[, 'C-Index'], conf.interval)
    colnames(result) <- c('C-index', '95%CI')
    
    return(result)
  } else {
    # 计算比较的p值，都是和第一个模型比较
    compare.cindex <- lapply(pred.models.coxph[-1], function(x){rcorrp.cens(pred.models.coxph[[1]], x, Surv(time = time, event = event))})
    p.value <- c("-", unlist(lapply(compare.cindex, function(x)(round(2*(1 - pnorm(x[['Dxy']] / x[['S.D.']])), digits=4))))) 
    
    result <- do.call(rbind, models.filter.result)
    conf.interval <- paste(result[, 'Lower'], result[, 'Upper'], sep = "-")
    result <- data.frame(result[, 'C-Index'], conf.interval, p.value)
    colnames(result) <- c('C-index', '95%CI', 'p-value')
    
    return(result)
  }
}  


CIndexBySubtype <- function(cli.data, subtype, grade, models, calulate.p = FALSE){

 cli.dat <- subset(cli.data, Integrated_Diagnoses %in% subtype & histological_grade %in% grade) 
 
 cli.dat$subclo_genome_frac <- cut_width(x=cli.dat$subclo_genome_frac, width=0.05, labels=FALSE)
 cli.dat$non_neutral_genome_frac <- cut_width(x=cli.dat$non_neutral_genome_frac, width=0.05, labels=FALSE)
 cli.dat$clo_genome_frac <- cut_width(x=cli.dat$clo_genome_frac, width=0.05, labels=FALSE)
 
 
 mul.c.index <- lapply(c('os', 'dss', 'pfi'), function(sur.type){
  
  if(sur.type == 'os'){
   sur.dat <- subset(cli.dat, !(is.na(os) | is.na(os_time)))
   c.index <- cIndex(sur.dat, sur.dat$os_time, sur.dat$os, models, diff=calulate.p)
  
  }else if(sur.type == 'dss'){
   sur.dat <- subset(cli.dat, !(is.na(dss) | is.na(dss_time)))
   c.index <- cIndex(sur.dat, sur.dat$dss_time, sur.dat$dss, models, diff=calulate.p)
  
  }else{
   sur.dat <- subset(cli.dat, !(is.na(pfi) | is.na(pfi_time)))
   c.index <- cIndex(sur.dat, sur.dat$pfi_time, sur.dat$pfi, models, diff=calulate.p)
  
  }
  return(c.index)
 })
 
 mul.c.index <- do.call(rbind, mul.c.index)
 return(mul.c.index)
}

# load data
gli.cn.alt.frac <- readRDS(file='/data/gli_glod_cn_alt_frac.rds')
gli.cn.alt.frac <- subset(gli.cn.alt.frac, !is.na(histological_grade))


models <- list(c('subclo_genome_frac'), c('non_neutral_genome_frac'), c('clo_genome_frac'), 
 c('age', 'histological_grade'), 
 c('age', 'histological_grade', 'subclo_genome_frac'),
 c('age', 'histological_grade', 'non_neutral_genome_frac'), 
 c('age', 'histological_grade', 'clo_genome_frac'))

cindex.res <- CIndexBySubtype(gli.cn.alt.frac, 'Astrocytoma,IDHmut', c('G2', 'G3', 'G4'), models)


cindex.res <- data.frame(sur_type=rep(c('os', 'dss', 'pfi'), each=7), 
 models = rep(c('subclonal SCNAs', 'SCNAs', 'clonal SCNAs', 'age + CNS WHO grade', 
 'age + CNS WHO grade + subclonal SCNAs', 'age + CNS WHO grade + SCNAs', 'age + CNS WHO grade + clonal SCNAs'), times=3),
 cindex = cindex.res[, 'C-index'], ci_low=do.call(rbind, strsplit(cindex.res[, '95%CI'], split='-'))[, 1], 
 ci_high=do.call(rbind, strsplit(cindex.res[, '95%CI'], split='-'))[, 2])
library('dplyr')
cindex.res <- cindex.res %>% mutate(sur_type=factor(sur_type, levels=c('os', 'dss', 'pfi')), 
 models=factor(models, levels=c('subclonal SCNAs', 'SCNAs', 'clonal SCNAs', 'age + CNS WHO grade', 
 'age + CNS WHO grade + subclonal SCNAs', 'age + CNS WHO grade + SCNAs', 'age + CNS WHO grade + clonal SCNAs')),
 ci_low=as.numeric(ci_low), ci_high=as.numeric(ci_high))


cindex.plot <- ggplot(data=cindex.res, aes(models, cindex)) + geom_bar(aes(fill=models), stat="identity") + 
 geom_errorbar(aes(ymax=ci_high, ymin=ci_low), width=0.2) + facet_grid(. ~ sur_type) + 
 theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position="none")

setwd('/result/Section3/')
ggsave(cindex.plot, file='non_codel_g234_cindex.pdf')



# more compare
compare.models <- list(c('non_neutral_genome_frac'), c('clo_genome_frac'))
# > CIndexBySubtype(gli.cn.alt.frac, 'Astrocytoma,IDHmut', c('G2', 'G3', 'G4'), compare.models, TRUE)
    # C-index                               95%CI            p-value
# 1 0.6911318 0.619780056379083-0.762483521288976                  -
# 2 0.6615078 0.589337210155293-0.733678361436048             1.1419
# 3 0.6932715 0.614496276896878-0.772046646536996                  -
# 4 0.6783063 0.602427571793403-0.754184957208917             1.0699
# 5 0.5958181 0.532593999424523-0.659042248790046                  -
# 6 0.6087526 0.547666937638359-0.669838220273084 0.0475462439635823
# > CIndexBySubtype(gli.cn.alt.frac, 'Astrocytoma,IDHmut', c('G2', 'G3'), compare.models, TRUE)
    # C-index                               95%CI            p-value
# 1 0.6534477 0.554766023092232-0.752129296664824                  -
# 2 0.6554126 0.566072782483686-0.744752512264438 0.0857627189341428
# 3 0.6511949 0.541728506979297-0.760661198903056                  -
# 4 0.6748621  0.586962228012622-0.76276203669326  0.077751153005357
# 5 0.5706468 0.499332279980024-0.641961374717467                  -
# 6 0.5983153    0.530020664607-0.666609930572538 0.0106623873680873


compare.models <- list(
 c('age', 'histological_grade', 'non_neutral_genome_frac'), c('age', 'histological_grade', 'clo_genome_frac'))
# > CIndexBySubtype(gli.cn.alt.frac, 'Astrocytoma,IDHmut', c('G2', 'G3', 'G4'), compare.models, TRUE)
    # C-index                               95%CI             p-value
# 1 0.7487748 0.681233663632132-0.816315836465888                   -
# 2 0.7635758  0.69561910350302-0.831532435410918  0.0103364942216158
# 3 0.7526493  0.67785183013311-0.827446821119298                   -
# 4 0.7731214  0.70036882679735-0.845873947769124  0.0154082036595775
# 5 0.6365765  0.56953675421499-0.703616253455579                   -
# 6 0.6436415  0.57807377864933-0.709209224984098 0.00935412955079373
# > CIndexBySubtype(gli.cn.alt.frac, 'Astrocytoma,IDHmut', c('G2', 'G3'), compare.models, TRUE)
    # C-index                               95%CI            p-value
# 1 0.6602358 0.568212459054752-0.752259137944177                  -
# 2 0.6732762 0.579663563992205-0.766888776129267  0.044478523584889
# 3 0.6440717 0.545166079403666-0.742977302949276                  -
# 4 0.6769301 0.577607929899445-0.776252364218202 0.0336298706823268
# 5 0.5906911 0.517137509457677-0.664244684344586                  -
# 6 0.6082145  0.534309331517456-0.68211959125677 0.0508143496749791


