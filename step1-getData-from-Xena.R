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

# TCGA-XENA数据库 打包：https://share.weiyun.com/56URQ3a
# TCGA-GDC-somatic  链接：https://share.weiyun.com/5fx40jk 密码：7yenp9 

rm(list=ls())
options(stringsAsFactors = F)
if(file.exists('UCSC-XENA/miRNA_GA_gene.gz')){
  
  miRNA_GA=read.table('UCSC-XENA/miRNA_GA_gene.gz',header = T,sep = '\t')
  dim(miRNA_GA)
  miRNA_GA[1:4,1:4]
  dim(na.omit(miRNA_GA))
  miRNA_HiSeq=read.table('UCSC-XENA/miRNA_HiSeq_gene.gz',header = T,sep = '\t')
  dim(miRNA_HiSeq)
  miRNA_HiSeq[1:4,1:4]
  dim(na.omit(miRNA_HiSeq))
}


