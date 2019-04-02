# https://portal.gdc.cancer.gov/projects/TCGA-LUSC
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
# https://github.com/jmzeng1314/tcga_example/blob/master/step1-getData-from-GDC.R

rm(list=ls())
# Load the packages required to read XML files.
library("XML")
library("methods")
dir='/Users/jimmy/data/tcga/lusc/clinical/'
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
   
  dir='/Users/jimmy/data/tcga/lusc/miRNA/'
  fls=list.files(path = dir ,pattern='*.mirnas.quantification.txt$',recursive=T)
  fls
  mi = lapply(fls
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

# Load the package required to read JSON files.
library("rjson")

# Give the input file name to the function.
result <- fromJSON(file = "tcga-lusc-miRNA.json")

# Print the result.
print(result)
fls=unlist(lapply(result,function(x){x[[1]]}))
cid=unlist(lapply(result,function(x){x[[6]][[1]][[2]]}))
id2fls=data.frame(cid=cid,fls=fls)

save(id2fls,mi_df,fls,cl_df,file = 'tcga_lusc_gdc_miRNA-clinical.Rdata')
rm(list = ls())
load(file = 'tcga_lusc_gdc_miRNA-clinical.Rdata')





