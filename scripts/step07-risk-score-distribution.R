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

# https://mp.weixin.qq.com/s/J-vaFq1Vv-zR1LC0Wop1zw
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5051943/figure/djw200-F2/
# https://www.ncbi.nlm.nih.gov/pubmed/24893932
# Sixteen of the 111 most significantly altered miRNAs were associated with OS across different clinical subclasses of the TCGA-derived LUAD cohort. 
# A linear prognostic model of eight miRNAs (miR-31, miR-196b, miR-766, miR-519a-1, miR-375, miR-187, miR-331 and miR-101-1) was constructed and weighted by the importance scores from the supervised principal component method to divide patients into high- and low-risk groups. 
# Patients assigned to the high-risk group exhibited poor OS compared with patients in the low-risk group (hazard ratio [HR]=1.99, P <0.001). 
# The eight-miRNA signature is an independent prognostic marker of OS of LUAD patients and demonstrates good performance for predicting 5-year OS (Area Under the respective ROC Curves [AUC] = 0.626, P = 0.003), especially for non-smokers (AUC = 0.686, P = 0.023).

rm(list=ls())
options(stringsAsFactors = F)

Rdata_dir='../Rdata/'
Figure_dir='../figures/'


library(survival)
library(survminer)

# 这里举例的文章不一样，所以不再使用前面步骤的数据。

# 而是对TCGA-LUAD-miRNA重新处理，正好跟前面的步骤对应学习。

if(F){
  library(RTCGA.miRNASeq) 
  ??RTCGA.miRNASeq
  s=rownames(LUAD.miRNASeq)[seq(1,nrow(LUAD.miRNASeq),by=3)]
  expr <- expressionsTCGA(LUAD.miRNASeq)
  dim(expr)
  expr[1:40,1:4]
  expr=as.data.frame(expr[seq(1,nrow(expr),by=3),3:ncol(expr)])
  mi=colnames(expr)
  expr=apply(expr,1,as.numeric) 
  colnames(expr)=s
  rownames(expr)=mi
  expr[1:4,1:4]
  expr=na.omit(expr)
  dim(expr)
  expr=expr[apply(expr, 1,function(x){sum(x>1)>10}),]
  dim(expr)
  
  library(RTCGA.clinical) 
  ??RTCGA.clinical
  meta <- LUAD.clinical
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
  
  save(expr,meta,file = file.path(Rdata_dir,'TCGA-LUAD-miRNA-example.Rdata'))
} 
load(file = file.path(Rdata_dir,'TCGA-LUAD-miRNA-example.Rdata'))

group_list=ifelse(as.numeric(substr(colnames(expr),14,15)) < 10,'tumor','normal')
table(group_list)

if(F){
  exprSet=na.omit(expr)
  exprSet=exprSet[,group_list=='tumor']
  
  head(meta)
  colnames(meta)
  meta[,3][is.na(meta[,3])]=0
  meta[,4][is.na(meta[,4])]=0
  meta$days=as.numeric(meta[,3])+as.numeric(meta[,4])
  meta=meta[,c(1:2,5:9)]
  colnames(meta)
  colnames(meta)=c('ID','event','race','age','gender','stage',"days") 
  library(survival)
  library(survminer)
  meta$event=ifelse(meta$event=='alive',0,1)
  meta$age=as.numeric(meta$age)
  library(stringr) 
  meta$stage=str_split(meta$stage,' ',simplify = T)[,2]
  table(  meta$stage)
  boxplot(meta$age)
  meta$age_group=ifelse(meta$age>median(meta$age,na.rm = T),'older','younger')
  table(meta$race)
  meta$time=meta$days/30
  phe=meta
  
  head(phe)
  phe$ID=toupper(phe$ID) 
  phe=phe[match(substr(colnames(exprSet),1,12),phe$ID),]
  head(phe)
  exprSet[1:4,1:4]
  
  save(exprSet,phe,file=file.path(Rdata_dir,'TCGA-LUAD-survival_input.Rdata'))
}

load(file=file.path(Rdata_dir,'TCGA-LUAD-survival_input.Rdata'))
head(phe)
exprSet[1:4,1:4]
dim(exprSet)

# 这个时候是515个病人的673个miRNA表达矩阵。

## 挑选感兴趣的基因构建coxph模型 

# miR-31, miR-196b, miR-766, miR-519a-1, miR-375, miR-187, miR-331 and miR-101-1
# hsa-mir-31,hsa-mir-196b,hsa-mir-766,hsa-mir-519a-1,hsa-mir-375,hsa-mir-187,hsa-mir-331,hsa-mir-101-1
e=t(exprSet[c('hsa-mir-31','hsa-mir-196b','hsa-mir-766','hsa-mir-519a-1','hsa-mir-375','hsa-mir-187','hsa-mir-331','hsa-mir-101-1'),])
e=log2(e+1)
colnames(e)=c('miR31','miR196b','miR766','miR519a1','miR375','miR187','miR331','miR101')
dat=cbind(phe,e)

dat$gender=factor(dat$gender)
dat$stage=factor(dat$stage)

colnames(dat) 
s=Surv(time, event) ~ miR31+miR196b+miR766+miR519a1+miR375+miR187+miR331+miR101
model <- coxph(s, data = dat )
summary(model,data=dat)
options(scipen=1)
ggforest(model, data =dat, 
         main = "Hazard ratio", 
         cpositions = c(0.10, 0.22, 0.4), 
         fontsize = 1.0, 
         refLabel = "1", noDigits = 4)
 

# 对于生存数据预后模型评价很多采用C-index ，但c-index展示没有roc曲线图来的直观
new_dat=dat

fp <- predict(model,new_dat,type="risk");boxplot(fp)
fp <- predict(model,new_dat,type="expected") ;boxplot(fp)
plot(fp,phe$days)
fp <- predict(model,new_dat) ;boxplot(fp)
basehaz(model) 
library(Hmisc)
options(scipen=200)
with(new_dat,rcorr.cens(fp,Surv(time, event)  ))
# http://dni-institute.in/blogs/cox-regression-interpret-result-and-predict/

library(cowplot)
library(pheatmap)
# https://cran.r-project.org/web/packages/cowplot/vignettes/introduction.html
fp_dat=data.frame(s=1:length(fp),v=as.numeric(sort(fp )))
sur_dat=data.frame(s=1:length(fp),
                   t=phe[names(sort(fp )),'time'] ,
                   e=phe[names(sort(fp )),'event']  ) 
sur_dat$e=ifelse(sur_dat$e==0,'alive','death')
exp_dat=new_dat[names(sort(fp )),10:17]
plot.point=ggplot(fp_dat,aes(x=s,y=v))+geom_point();print(plot.point)
plot.sur=ggplot(sur_dat,aes(x=s,y=t))+geom_point(aes(col=e));print(plot.sur)

mycolors <- colorRampPalette(c("black", "green", "red"), bias = 1.2)(100)
tmp=t(scale(exp_dat))
tmp[tmp > 1] = 1
tmp[tmp < -1] = -1
plot.h=pheatmap(tmp,col= mycolors,show_colnames = F,cluster_cols = T)
plot.h=pheatmap(tmp,col= mycolors,show_colnames = F,cluster_cols = F)
plot_grid(plot.point, plot.sur, plot.h$gtable,
          labels = c("A", "B","C"),
          align = 'v',ncol = 1)

