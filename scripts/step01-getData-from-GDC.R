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
### Update Log: 2018-10-10  second version
###
### ---------------

# https://gdc.cancer.gov/access-data/gdc-data-transfer-tool
# mkdir -p ~/biosoft/gdc_client
# cd ~/biosoft/gdc_client/
# wget https://gdc.cancer.gov/system/files/authenticated%20user/0/gdc-client_v1.3.0_OSX_x64.zip 
# wget https://gdc.cancer.gov/system/files/authenticated%20user/0/gdc-client_v1.3.0_Ubuntu14.04_x64.zip

## ./gdc-client --help
## ./gdc-client download --help 
#  mkdir clinical
# ./gdc-client download -m ~/Documents/Nutstore/github/TCGA-KIRC-miRNA-example/GDC/gdc_manifest.2018-08-05-clinical.txt -d clinical
#  Successfully downloaded: 522
# mkdir miRNAseq
# ./gdc-client download -m ~/Documents/Nutstore/github/TCGA-KIRC-miRNA-example/GDC/gdc_manifest.2018-08-05-LUAD-miRNA-seq.txt -d miRNAseq
#  Successfully downloaded: 567
#  或者直接从微云下载：链接：https://share.weiyun.com/5XsyuzH 密码：68pm7e 

rm(list=ls())
options(stringsAsFactors = F)

m='GDC/gdc_manifest.2018-08-05-clinical.txt'
x1=read.table(m,header = T) 
m='GDC/gdc_manifest.2018-08-05-LUAD-miRNA-seq.txt'
x2=read.table(m,header = T) 

# Load the packages required to read XML files.
library("XML")
library("methods")
dir='/Users/jmzeng/biosoft/gdc_client/miRNAseq/'
if(dir.exists(dir)){
  
  all_fiels=list.files(path = dir ,pattern='*.xml$',recursive=T)
  all_fiels
  cl = lapply(all_fiels
              , function(x){
                #x=all_fiels[1]
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
  
  
}
