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
if(! require("devtools")) install.packages("devtools")
if(! require("reshape2")) install.packages("reshape2")
if(! require("ggplot2")) install.packages("ggplot2")
if(! require("pheatmap")) install.packages("pheatmap")
if(! require("ggfortify")) install.packages("ggfortify")
if(! require("stringr")) install.packages("stringr")
if(! require("survival")) install.packages("survival")
if(! require("survminer")) install.packages("survminer")
if(! require("lars")) install.packages("lars")
if(! require("glmnet")) install.packages("glmnet")

if(! require("timeROC")) install.packages("timeROC")
if(! require("ggpubr")) install.packages("ggpubr")

if(! require("randomForest")) install.packages("randomForest")
if(! require("ROCR")) install.packages("ROCR")
if(! require("genefilter")) install.packages("genefilter")
if(! require("Hmisc")) install.packages("Hmisc")
 
if(! require("caret")) install.packages("caret")
if(! require("genefilter")) install.packages("genefilter")
if(! require("Hmisc")) install.packages("Hmisc")


library(devtools) 
source("http://bioconductor.org/biocLite.R") 
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") 
library(BiocInstaller)
biocLite(c('airway','DESeq2','edgeR','limma'))

if(! require("CLL")) biocLite("CLL")
if(! require("org.Hs.eg.db")) biocLite(org.Hs.eg.db)
if(! require("maftools")) biocLite("maftools")







