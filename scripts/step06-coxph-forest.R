### ---------------
###
### Create: Jianming Zeng
### Date: 2019-04-02 21:59:01
### Email: jmzeng1314@163.com
### Blog: http://www.bio-info-trainee.com/
### Forum:  http://www.biotrainee.com/thread-1376-1-1.html
### CAFS/SUSTC/Eli Lilly/University of Macau
### Update Log:   2019-04-02  second version
###
### ---------------

### https://github.com/jmzeng1314/GEO/blob/master/GSE11121/step5-surivival.R

rm(list=ls())
options(stringsAsFactors = F)

Rdata_dir='../Rdata/'
Figure_dir='../figures/'
# 加载上一步从RTCGA.miRNASeq包里面提取miRNA表达矩阵和对应的样本临床信息。
load( file = 
        file.path(Rdata_dir,'TCGA-KIRC-miRNA-example.Rdata')
)
dim(expr)
dim(meta)
# 可以看到是 537个病人，但是有593个样本，每个样本有 552个miRNA信息。
# 当然，这个数据集可以下载原始测序数据进行重新比对，可以拿到更多的miRNA信息

# 这里需要解析TCGA数据库的ID规律，来判断样本归类问题。
group_list=ifelse(as.numeric(substr(colnames(expr),14,15)) < 10,'tumor','normal')

table(group_list)
exprSet=na.omit(expr)

table(group_list)
exprSet=na.omit(expr)
dim(exprSet)
load(  file = 
         file.path(Rdata_dir,'TCGA-KIRC-miRNA-survival_input.Rdata')
)
dim(exprSet) ## remove the nomral
head(phe)
exprSet[1:4,1:4]
head(colnames(exprSet))
head(phe$ID)
## 必须保证生存资料和表达矩阵，两者一致
all(substring(colnames(exprSet),1,12)==phe$ID)

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







