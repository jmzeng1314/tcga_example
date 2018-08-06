## 
### ---------------
###
### Create: Jianming Zeng
### Date: 2018-08-010 17:07:49
### Email: jmzeng1314@163.com
### Blog: http://www.bio-info-trainee.com/
### Forum:  http://www.biotrainee.com/thread-1376-1-1.html
### CAFS/SUSTC/Eli Lilly/University of Macau
### Update Log: 2018-08-010  First version
###
### ---------------


rm(list=ls())
# Load the packages required to read XML files.
library("XML")
library("methods")
dir='/home/jmzeng/project/tcga/clinical/'
cl = lapply(list.files(path = dir ,pattern='*.xml$',recursive=T)
            , function(x){
              
              result <- xmlParse(file = file.path(dir,x)) 
              rootnode <- xmlRoot(result)  
              xmldataframe <- xmlToDataFrame( rootnode[2] ) 
              return(t(xmldataframe))
            })

cl_df <- t(do.call(cbind,cl))
save(cl_df,file = 'GDC_TCGA_LUAD_clinical_df.Rdata')

dir='/home/jmzeng/project/tcga/miRNAseq/'

mi = lapply(list.files(path = dir ,pattern='*.mirnas.quantification.txt$',recursive=T)
            , function(x){
              
              result <- read.table(file = file.path(dir,x),sep = '\t',header = T)[,1:2]  
              return( result )
            })

mi_df <- t(do.call(cbind,mi))
dim(mi_df)
mi_df[1:4,1:4]
colnames(mi_df)=mi_df[1,]
mi_df=mi_df[seq(2,nrow(mi_df),by=2),]
mi_df[1:4,1:4]


