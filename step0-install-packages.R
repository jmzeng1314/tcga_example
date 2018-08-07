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

# https://github.com/jmzeng1314/biotrainee 

options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") 
source("http://bioconductor.org/biocLite.R") 

packs = c("devtools", "reshape2", "ggplot2", "pheatmap", "ggfortify", "stringr", "survival",
          "survminer", "lars", "glmnet", "timeROC", "ggpubr", "randomForest", "ROCR", "genefilter",
          "Hmisc", "caret", "airway","DESeq2","edgeR","limma", "CLL", "org.Hs.eg.db", "maftools")

if(! require(pacman)) install.packages("pacman", dependencies = TRUE)
pacman::p_load(packs, dependencies=TRUE, character.only = TRUE)

# check
pacman::p_loaded(packs, character.only = TRUE)
all(pacman::p_loaded(packs, character.only = TRUE))
