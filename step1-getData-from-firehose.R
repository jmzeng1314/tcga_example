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
rm(list=ls())
options(stringsAsFactors = F)
kirc=read.table('broad-firehose/KIRC/expr.txt',header = T,sep='\t',fill=T)
dim(kirc)
kirc[1:4,1:4]
dim(na.omit(kirc))
ffpe=read.table('broad-firehose/KIRC-FFPE/expr.txt',header = T,sep='\t',fill=T)
dim(ffpe)
ffpe[1:4,1:4]
dim(na.omit(ffpe))
clinical=read.table('broad-firehose/clinical/KIRC.merged_only_clinical_clin_format.txt',
                    header = T,sep='\t',fill=T)

clinical[1:4,1:4]


