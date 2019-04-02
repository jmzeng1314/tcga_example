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
### Update Log: 2018-10-10  second version
###
### ---------------
rm(list=ls())
options(stringsAsFactors = F)
d='../Rdata/KIRC-broad-firehose/' 
## too many NA in the miRNA expression matrix from XENA

kirc=read.table(file.path(d,'KIRC/expr.txt'),header = T,sep='\t',fill=T)
dim(kirc)
kirc[1:4,1:4]
dim(na.omit(kirc))
expr=as.data.frame(expr[seq(1,nrow(expr),by=3),3:ncol(expr)])
kirc=kirc[-1,]
rownames(kirc)=kirc[,1]
kirc=kirc[,-1]
kirc=kirc[,seq(1,ncol(kirc),by=3)]
dim(kirc)
kirc[1:4,1:4]


ffpe=read.table(file.path(d,'KIRC-FFPE/expr.txt'),header = T,sep='\t',fill=T)
dim(ffpe)
ffpe[1:4,1:4]
dim(na.omit(ffpe))
clinical=read.table(file.path(d,'clinical/KIRC.merged_only_clinical_clin_format.txt'),
                    header = T,sep='\t',fill=T)

clinical[1:4,1:4]


