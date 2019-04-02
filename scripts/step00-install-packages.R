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

## 强调，不是所有的R包都需要安装成功的。
## 强调，不是所有的R包都需要安装成功的。
## 强调，不是所有的R包都需要安装成功的。
## 强调，不是所有的R包都需要安装成功的。
## 失败就失败，大不了从头再来，卸载R语言，从新开始。

## 强调，中国大陆的粉丝务必注意下载镜像。
## 强调，管是什么电脑，都请务必安装好R及Rstudio哦
# 所有的软件都安装在c盘哦，然后系统用户名最好是不要用中文，写代码最怕中文字符串哦！
# 生信0基础第一步，下载R和Rstudio并且安装在自己的电脑上面。官网链接是 
# - R: https://mirrors.tuna.tsinghua.edu.cn/CRAN/
# - RStudio：https://www.rstudio.com/products/rstudio/download/#download 
# 如果你的网络不好，可以从我整理的网盘下载，链接：https://share.weiyun.com/5hW6VAA  密码：3fuhrm



rm(list = ls()) 
#清空当前工作空间变量  
options()$repos  
#查看当前工作空间默认的下载包路径
options()$BioC_mirror 
#查看使用BioCManager下载包的默认路径
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") 
# 指定使用BioCManager下载的路径
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")) 
# 指定使用install.packages下载包的路径
options()$repos 
options()$BioC_mirror


# https://bioconductor.org/packages/release/bioc/html/GEOquery.html
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager") 
#判断是否存在BiocManger包，没有的话下载该包

#BiocManager::install("KEGG.db",ask = F,update = F)
#BiocManager::install(c("GSEABase","GSVA","clusterProfiler" ),ask = F,update = F)
#BiocManager::install(c("GEOquery","limma","impute" ),ask = F,update = F)
#BiocManager::install(c("genefu","org.Hs.eg.db","hgu133plus2.db" ),ask = F,update = F)

#判断是否存在这些包，不存在的话安装这些包
if(!require("KEGG.db")) BiocManager::install("KEGG.db",ask = F,update = F)
if(!require("GSEABase")) BiocManager::install("GSEABase",ask = F,update = F)
if(!require("GSVA")) BiocManager::install("GSVA",ask = F,update = F)
if(!require("clusterProfiler")) BiocManager::install("clusterProfiler",ask = F,update = F)
if(!require("GEOquery")) BiocManager::install("GEOquery",ask = F,update = F)
if(!require("limma")) BiocManager::install("limma",ask = F,update = F)
if(!require("impute")) BiocManager::install("impute",ask = F,update = F)
if(!require("genefu")) BiocManager::install("genefu",ask = F,update = F)
if(!require("org.Hs.eg.db")) BiocManager::install("org.Hs.eg.db",ask = F,update = F)
if(!require("hgu133plus2.db")) BiocManager::install("hgu133plus2.db",ask = F,update = F)
if(!require("ConsensusClusterPlus")) BiocManager::install("ConsensusClusterPlus",ask = F,update = F)

BiocManager::install(c('airway','DESeq2','edgeR','limma'),
                     ask = F,update = F)

if(! require("maftools")) BiocManager::install("maftools",ask = F,update = F)
if(! require("genefilter")) BiocManager::install("genefilter",ask = F,update = F)

if(! require("CLL")) BiocManager::install("CLL",ask = F,update = F)
if(! require("org.Hs.eg.db")) BiocManager::install('org.Hs.eg.db',ask = F,update = F)

if(! require("maftools")) BiocManager::install("maftools",ask = F,update = F)
if(! require("RTCGA")) BiocManager::install("RTCGA",ask = F,update = F)
if(! require("RTCGA.clinical")) BiocManager::install("RTCGA.clinical",ask = F,update = F)
# https://bioconductor.org/packages/3.6/data/experiment/src/contrib/RTCGA.clinical_20151101.8.0.tar.gzn
if(! require("RTCGA.miRNASeq")) BiocManager::install("RTCGA.miRNASeq",ask = F,update = F)


# 如果万一是 R3.4 版本之前的，请使用下面代码下载该项目所需要的包
# source("https://bioconductor.org/BiocManager::install.R") 
# library('BiocInstaller') 
# options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") 
# BiocInstaller::BiocManager::install("GEOquery")
# BiocInstaller::BiocManager::install(c("limma"))
# BiocInstaller::BiocManager::install(c("impute"))

options()$repos

#install.packages('WGCNA')
#install.packages(c("FactoMineR", "factoextra"))
#install.packages(c("ggplot2", "pheatmap","ggpubr"))

#判断是否存在这些包，不存在的话安装这些包
if(!require("WGCNA")) install.packages("WGCNA",update = F,ask = F)
if(!require("FactoMineR")) install.packages("FactoMineR",update = F,ask = F)
if(!require("factoextra")) install.packages("factoextra",update = F,ask = F)
if(!require("ggplot2")) install.packages("ggplot2",update = F,ask = F)
if(!require("pheatmap")) install.packages("pheatmap",update = F,ask = F)
if(!require("ggpubr")) install.packages("ggpubr",update = F,ask = F)
if(!require("glmnet")) install.packages("glmnet",update = F,ask = F)
if(!require("randomForest")) install.packages("randomForest",update = F,ask = F)



library("FactoMineR")
library("factoextra")

library(GSEABase)
library(GSVA)
library(clusterProfiler)
library(genefu)
library(ggplot2)
library(ggpubr)
library(hgu133plus2.db)
library(limma)
library(org.Hs.eg.db)
library(pheatmap)

install.packages(c("devtools","reshape2","pheatmap",
                   "ggplot2","ggfortify","stringr",
                   "survival","survminer","lars",
                   "glmnet","timeROC","ggpubr",
                   "randomForest","ROCR","Hmisc",
                   "caret","ggstatsplot","tableone", 
                   "devtools","reshape2","randomForest"))

library(devtools) 

