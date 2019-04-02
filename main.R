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


setwd('scripts/')

### ---------------
###
### install the packages 
###
### ---------------
#source('step00-install-packages.R')
# we don't need to run it by source.

### ---------------
###
### Get miRNA expression matrix from RTCGA package 
###
### ---------------
source('step01-getData-from-RTCGA.R') 
dim(expr)
dim(meta)
# 可以看到是 537个病人，但是有593个样本，每个样本有 552个miRNA信息。
# 当然，这个数据集可以下载原始测序数据进行重新比对，可以拿到更多的miRNA信息

# 还有一些其它获取表达矩阵的方式，这里不演示。

### ---------------
###
### Do DEG for miRNA expression matrix by DESeq2,edgeR,limma
###
### ---------------
source('step02-DEG-3-packages.R')

### ---------------
###
###  
###
### ---------------
source('step03-batch-logRank.R')
dim(expr)
dim(meta)


### ---------------
###
###  
###
### ---------------
source('step04-batch-coxp.R')
 

### ---------------
###
### 
###
### ---------------
source('step05-lasso.R') 


### ---------------
###
###  
###
### ---------------
source('step06-coxph-forest.R') 

### ---------------
###
### Get miRNA expression matrix from RTCGA package 
###
### ---------------
source('step07-risk-score-distribution.R')
dim(expr)
dim(meta)

### ---------------
###
### Get miRNA expression matrix from RTCGA package 
###
### ---------------
source('step08-Random-foreast.R')
dim(expr)
dim(meta)

### ---------------
###
### Get miRNA expression matrix from RTCGA package 
###
### ---------------
source('step09-miRNA-downstream.R')
dim(expr)
dim(meta)

### ---------------
###
### Get miRNA expression matrix from RTCGA package 
###
### ---------------
source('step10-maftools.R')
dim(expr)
dim(meta)

### ---------------
###
### Get miRNA expression matrix from RTCGA package 
###
### ---------------
source('step11-boxplot.R')
dim(expr)
dim(meta)

### ---------------
###
### Get miRNA expression matrix from RTCGA package 
###
### ---------------
source('step12-correlation.R')
dim(expr)
dim(meta)

### ---------------
###
### Get miRNA expression matrix from RTCGA package 
###
### ---------------
source('step13-split-cohort.R')
dim(expr)
dim(meta)

### ---------------
###
### Get miRNA expression matrix from RTCGA package 
###
### ---------------
source('step14-timeROC.R')
dim(expr)
dim(meta)

### ---------------
###
### Get miRNA expression matrix from RTCGA package 
###
### ---------------
source('step15-choose_lncRNA.R')
dim(expr)
dim(meta)

### ---------------
###
### Get miRNA expression matrix from RTCGA package 
###
### ---------------
source('step16-clinical-tables.R')
dim(expr)
dim(meta)


### ---------------
###
### Get miRNA expression matrix from RTCGA package 
###
### ---------------
source('step17-mutation-signatures.R')
dim(expr)
dim(meta)

### ---------------
###
### Get miRNA expression matrix from RTCGA package 
###
### ---------------
source('step18-SVM.R')
dim(expr)
dim(meta)

### ---------------
###
### Get miRNA expression matrix from RTCGA package 
###
### ---------------
source('step17-others.R')
dim(expr)
dim(meta)


