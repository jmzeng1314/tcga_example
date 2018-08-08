## 
### ---------------
###
### Create: Jianming Zeng
### Date: 2018-08-10 17:07:49
### Email: jmzeng1314@163.com
### Blog: http://www.bio-info-trainee.com/
### Forum:  http://www.biotrainee.com/thread-1376-1-1.html
### CAFS/SUSTC/Eli Lilly/University of Macau
### Update Log: 2018-08-10  First version
###
### ---------------

### https://github.com/jmzeng1314/GEO/blob/master/GSE11121/step5-surivival.R
# 进行COX比例风险模型。
# HR = 1: 无影响
# HR < 1: 风险降低
# HR > 1: 风险增高

rm(list=ls())
library(survival)
library(survminer)
load(file = 'TCGA-KIRC-miRNA-example.Rdata')
group_list=ifelse(substr(colnames(expr),14,15)=='01','tumor','normal')
table(group_list)
load(file='survival_input.Rdata')
head(phe)
exprSet[1:4,1:4]
## 挑选感兴趣的基因构建coxph模型 
# 2015-TCGA-ccRCC-5-miRNAs-signatures
# Integrated genomic analysis identifies subclasses and prognosis signatures of kidney cancer
# miR-21,miR-143,miR-10b,miR-192,miR-183
# hsa-mir-21,hsa-mir-143,hsa-mir-10b,hsa-mir-192,hsa-mir-183
e=t(exprSet[c('hsa-mir-21','hsa-mir-143','hsa-mir-10b','hsa-mir-192','hsa-mir-183'),])
e=log2(e)
colnames(e)=c('miR21','miR143','miR10b','miR192','miR183')
dat=cbind(phe,e)

dat$gender=factor(dat$gender)
dat$stage=factor(dat$stage)

colnames(dat)
s=Surv(time, event) ~ gender + stage + age + miR21+miR143+miR10b+miR192+miR183
s=Surv(time, event) ~ miR21+miR143+miR10b+miR192+miR183
model <- coxph(s, data = dat )
summary(model,data=dat)
options(scipen=1)
ggforest(model, data =dat, 
         main = "Hazard ratio", 
         cpositions = c(0.10, 0.22, 0.4), 
         fontsize = 1.0, 
         refLabel = "1", noDigits = 4)


fp <- predict(model)
summary(model,data=dat)
library(Hmisc)
options(scipen=200)
with(dat,rcorr.cens(fp,Surv(time, event)  ))

# 用于计算生存分析中的COX模型预测值与真实之间的区分度（discrimination），也称为Harrell's concordanceindex。
# C-index在0.5-1之间。0.5为完全不一致,说明该模型没有预测作用
# 1为完全一致,说明该模型预测结果与实际完全一致。


# 若要找到最佳模型，我们可以进行变量选择，可以采用逐步回归法进行分析







