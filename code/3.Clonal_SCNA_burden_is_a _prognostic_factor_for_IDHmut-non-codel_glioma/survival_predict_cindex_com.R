
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

 cli.dat <- subset(cli.data, IDH_CODEL_SUBTYPE %in% subtype & histological_grade %in% grade) 
 
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
gli.cn.alt.frac <- readRDS(file='/pub5/xiaoyun/Jobs/J22/CopyNumberClonalityProject/Resource/CuratedData/gli_glod_cn_alt_frac.rds')

models <- list(c('subclo_genome_frac'), c('non_neutral_genome_frac'), c('clo_genome_frac'), 
 c('age', 'histological_grade'), 
 c('age', 'histological_grade', 'subclo_genome_frac'),
 c('age', 'histological_grade', 'non_neutral_genome_frac'), 
 c('age', 'histological_grade', 'clo_genome_frac'))

cindex.res <- CIndexBySubtype(gli.cn.alt.frac, 'IDHmut-non-codel', c('G2', 'G3', 'G4'), models)


cindex.res <- data.frame(sur_type=rep(c('os', 'dss', 'pfi'), each=7), 
 models = rep(c('subclonal SCNAs', 'SCNAs', 'clonal SCNAs', 'age + grade', 
 'age + grade + subclonal SCNAs', 'age + grade + SCNAs', 'age + grade + clonal SCNAs'), times=3),
 cindex = cindex.res[, 'C-index'], ci_low=do.call(rbind, strsplit(cindex.res[, '95%CI'], split='-'))[, 1], 
 ci_high=do.call(rbind, strsplit(cindex.res[, '95%CI'], split='-'))[, 2])
library('dplyr')
cindex.res <- cindex.res %>% mutate(sur_type=factor(sur_type, levels=c('os', 'dss', 'pfi')), 
 models=factor(models, levels=c('subclonal SCNAs', 'SCNAs', 'clonal SCNAs', 'age + grade', 
 'age + grade + subclonal SCNAs', 'age + grade + SCNAs', 'age + grade + clonal SCNAs')),
 ci_low=as.numeric(ci_low), ci_high=as.numeric(ci_high))



cindex.plot <- ggplot(data=cindex.res, aes(models, cindex)) + geom_bar(aes(fill=models), stat="identity") + 
 geom_errorbar(aes(ymax=ci_high, ymin=ci_low), width=0.2) + facet_grid(. ~ sur_type) + 
 theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position="none")

setwd('/pub5/xiaoyun/Jobs/J22/CopyNumberClonalityProject/Results/Section1/Results/Survival')
ggsave(cindex.plot, file='non_codel_g234_cindex.pdf')



# more compare
compare.models <- list(c('non_neutral_genome_frac'), c('clo_genome_frac'))
# > CIndexBySubtype(gli.cn.alt.frac, 'IDHmut-non-codel', c('G2', 'G3', 'G4'), compare.models, TRUE)
    # C-index                               95%CI            p-value
# 1 0.6911318 0.619780056379083-0.762483521288976                  -
# 2 0.6615078 0.589337210155293-0.733678361436048             1.1419
# 3 0.6932715 0.614496276896878-0.772046646536996                  -
# 4 0.6783063 0.602427571793403-0.754184957208917             1.0699
# 5 0.5958181 0.532593999424523-0.659042248790046                  -
# 6 0.6087526 0.547666937638359-0.669838220273084 0.0475462439635823
# > CIndexBySubtype(gli.cn.alt.frac, 'IDHmut-non-codel', c('G2', 'G3'), compare.models, TRUE)
    # C-index                               95%CI             p-value
# 1 0.6667554 0.582341353027693-0.751169426849869                   -
# 2 0.6704818 0.593681155556757-0.747282379178404   0.126344193602323
# 3 0.6674781  0.576577039991239-0.75837914306621                   -
# 4 0.6864654 0.610832990118104-0.762097876483649   0.127714175964158
# 5 0.5782904  0.511054431438113-0.64552635546023                   -
# 6 0.6045536 0.539914302235751-0.669192947614459 0.00552078277664014


compare.models <- list(
 c('age', 'histological_grade', 'non_neutral_genome_frac'), c('age', 'histological_grade', 'clo_genome_frac'))
# > CIndexBySubtype(gli.cn.alt.frac, 'IDHmut-non-codel', c('G2', 'G3', 'G4'), compare.models, TRUE)
    # C-index                               95%CI             p-value
# 1 0.7334789  0.66644401373823-0.800513829026677                   -
# 2 0.7425940 0.677798747101335-0.807389251379485 0.00714265844518258
# 3 0.7324826 0.661970855855276-0.802994341360501                   -
# 4 0.7452436 0.677700061503339-0.812787177475779  0.0126676492347511
# 5 0.6078003  0.537697613595697-0.67790308470616                   -
# 6 0.6235915 0.555630401347693-0.691552585479795  0.0512025433077923
# > CIndexBySubtype(gli.cn.alt.frac, 'IDHmut-non-codel', c('G2', 'G3'), compare.models, TRUE)
    # C-index                               95%CI            p-value
# 1 0.6807293  0.603235550781429-0.75822305981213                  -
# 2 0.6907107 0.613150546955028-0.768270799864242 0.0477360258329893
# 3 0.6799740 0.599689692037884-0.760258376770945                  -
# 4 0.6975008 0.617113294793133-0.777888328056591 0.0444884185771397
# 5 0.5779409 0.507000259904445-0.648881505623815                  -
# 6 0.5917715 0.521428498602013-0.662114541142345 0.0206992269586812


