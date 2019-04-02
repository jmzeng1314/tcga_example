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
rm(list=ls())
options(stringsAsFactors = F)
# 注意，并不是说使用 RTCGA.miRNASeq包的数据是最佳选择，只是因为这个演示起来最方便。
# 因为GDC官网下载数据具有一定门槛，也不是每个人都必须学会的。
getwd()
Rdata_dir='../Rdata/'
Figure_dir='../figures/'

# 如果开启下面代码，就会从RTCGA.miRNASeq包里面提取miRNA表达矩阵和对应的样本临床信息。 
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
  ## 每次运行代码，就会重新生成文件。
  save(expr,meta,
       file = file.path(Rdata_dir,'TCGA-KIRC-miRNA-example.Rdata')
         )
}

## 我们已经运行了上面被关闭的代码，而且保存了miRNA表达矩阵和对应的样本临床信息
# 现在直接加载即可。
load( file = 
        file.path(Rdata_dir,'TCGA-KIRC-miRNA-example.Rdata')
)
dim(expr)
dim(meta)
# 可以看到是 537个病人，但是有593个样本，每个样本有 552个miRNA信息。
# 当然，这个数据集可以下载原始测序数据进行重新比对，可以拿到更多的miRNA信息

