setwd('/home/PJLAB/liangbilin/share/ICIs/RadClin/')

## 载入R包
library(rms)
library(survival)

############################## PFS ############################## 

## 读取LIHC数据
dat <- read.table("./results/merge/merge_data_PFS.tsv", header=TRUE, sep = '\t')
dd=datadist(dat)
options(datadist="dd")

## 构建COX比例风险模型
f2 <- psm(Surv(PFS, PFS_status) ~ Rad_score_pulmonary + Rad_score_enhancement + PLR + Bone_metastasis	, data =  dat, dist='lognormal')
med <- Quantile(f2) # 计算中位生存时间
surv <- Survival(f2) # 构建生存概率函数

## 绘制COX回归中位生存时间的Nomogram图
nom <- nomogram(f2, fun=function(x) med(lp=x),funlabel="PFS Time (month)")

## save fig
pdf(file = "PFS_normo.pdf",width =9, height = 6 ,family="GB1")
plot(nom)
dev.off()



############################## OS ##############################

## 读取LIHC数据
dat <- read.table("./results/merge/merge_data_OS.tsv", header=TRUE, sep = '\t')
dd=datadist(dat)
options(datadist="dd")

## 构建COX比例风险模型
f2 <- psm(Surv(OS, OS_status) ~ Rad_score_pulmonary + Liver_metastasis + Histologic_type, data =  dat, dist='lognormal')
med <- Quantile(f2) # 计算中位生存时间
surv <- Survival(f2) # 构建生存概率函数

## 绘制COX回归中位生存时间的Nomogram图
nom <- nomogram(f2, fun=function(x) med(lp=x),funlabel="OS Time (month)")

## save fig
pdf(file = "OS_normo.pdf",width =9, height = 6 ,family="GB1")
plot(nom)
dev.off()

