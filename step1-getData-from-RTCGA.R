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

if(F){ 
  options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
  
  source("http://bioconductor.org/biocLite.R") 
  options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")  
  if(! require('RTCGA')){
  BiocInstaller::biocLite('RTCGA',
                          ask = F,suppressUpdates = T)
  }
  if(! require('RTCGA.miRNASeq')){
  BiocInstaller::biocLite('RTCGA.miRNASeq',
                          ask = F,suppressUpdates = T)
  }
    if(! require('RTCGA.clinical')){
  BiocInstaller::biocLite('RTCGA.clinical',
                          ask = F,suppressUpdates = T)
    }
}

if(F){
  library(RTCGA.miRNASeq) 
  s=rownames(KIRC.miRNASeq)[seq(1,nrow(KIRC.miRNASeq),by=3)]
  expr <- expressionsTCGA(KIRC.miRNASeq)
  dim(expr)
  expr[1:40,1:4]
  expr=as.data.frame(expr[seq(1,nrow(expr),by=3),3:ncol(expr)])
  mi=colnames(expr)
  expr=apply(expr,1,as.numeric) 
  colnames(expr)=s
  rownames(expr)=mi
  expr[1:4,1:4]
  expr=na.omit(expr)
  expr=expr[apply(expr, 1,function(x){sum(x>1)>10}),]
  
  library(RTCGA.clinical) 
  meta <- KIRC.clinical
  tmp=as.data.frame(colnames(meta))
  meta[(grepl('patient.bcr_patient_barcode',colnames(meta)))]
  meta[(grepl('patient.days_to_last_followup',colnames(meta)))]
  meta[(grepl('patient.days_to_death',colnames(meta)))]
  meta[(grepl('patient.vital_status',colnames(meta)))]
  ## patient.race  # patient.age_at_initial_pathologic_diagnosis # patient.gender 
  # patient.stage_event.clinical_stage
  meta=as.data.frame(meta[c('patient.bcr_patient_barcode','patient.vital_status',
                            'patient.days_to_death','patient.days_to_last_followup',
                            'patient.race',
                            'patient.age_at_initial_pathologic_diagnosis',
                            'patient.gender' ,
                           'patient.stage_event.pathologic_stage')])
  #meta[(grepl('patient.stage_event.pathologic_stage',colnames(meta)))]
  
  save(expr,meta,file = 'TCGA-KIRC-miRNA-example.Rdata')
}
rm(list=ls())
load(file = 'TCGA-KIRC-miRNA-example.Rdata')
dim(expr)
dim(meta)