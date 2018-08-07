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

### https://github.com/jmzeng1314/GEO/blob/master/GSE11121/step5-surivival.R

rm(list=ls())
library(survival)
library(survminer)
load(file = 'TCGA-KIRC-miRNA-example.Rdata')
group_list=ifelse(substr(colnames(expr),14,15)=='01','tumor','normal')
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
  # R里面实现生存分析非常简单！
  
  # 用my.surv <- surv(OS_MONTHS,OS_STATUS=='DECEASED')构建生存曲线。
  # 用kmfit2 <- survfit(my.surv~TUMOR_STAGE_2009)来做某一个因子的KM生存曲线。
  # 用 survdiff(my.surv~type, data=dat)来看看这个因子的不同水平是否有显著差异，其中默认用是的logrank test 方法。
  # 用coxph(Surv(time, status) ~ ph.ecog + tt(age), data=lung) 来检测自己感兴趣的因子是否受其它因子(age,gender等等)的影响。
  
  library(survival)
  library(survminer)
  meta$event=ifelse(meta$event=='alive',0,1)
  meta$age=as.numeric(meta$age)
  library(stringr) 
  meta$stage=str_split(meta$stage,' ',simplify = T)[,2]
  table(  meta$stage)
  boxplot(meta$age)
  meta$age_group=ifelse(meta$age>median(meta$age),'older','younger')
  table(meta$race)
  meta$time=meta$days/30
  phe=meta
  
  head(phe)
  phe$ID=toupper(phe$ID) 
  phe=phe[match(substr(colnames(exprSet),1,12),phe$ID),]
  head(phe)
  exprSet[1:4,1:4]
  
  save(exprSet,phe,file='survival_input.Rdata')
}

load(file='survival_input.Rdata')
head(phe)
exprSet[1:4,1:4]
# 利用ggsurvplot快速绘制漂亮的生存曲线图
sfit <- survfit(Surv(time, event)~gender, data=phe)
sfit
summary(sfit)
ggsurvplot(sfit, conf.int=F, pval=TRUE)
ggsurvplot(sfit,palette = c("#E7B800", "#2E9FDF"),
           risk.table =TRUE,pval =TRUE,
           conf.int =TRUE,xlab ="Time in months", 
           ggtheme =theme_light(), 
           ncensor.plot = TRUE)
## 多个 ggsurvplots作图生存曲线代码合并代码公布
sfit1=survfit(Surv(time, event)~gender, data=phe)
sfit2=survfit(Surv(time, event)~age_group, data=phe)
splots <- list()
splots[[1]] <- ggsurvplot(sfit1,pval =TRUE, data = phe, risk.table = TRUE)
splots[[2]] <- ggsurvplot(sfit2,pval =TRUE, data = phe, risk.table = TRUE)
# Arrange multiple ggsurvplots and print the output
arrange_ggsurvplots(splots, print = TRUE,  ncol = 2, nrow = 1, risk.table.height = 0.4)



## 挑选感兴趣的基因做差异分析
# 2015-TCGA-ccRCC-5-miRNAs-signatures
# Integrated genomic analysis identifies subclasses and prognosis signatures of kidney cancer
# miR-21,miR-143,miR-10b,miR-192,miR-183
tmp=as.data.frame(rownames(exprSet))
g='hsa-mir-21'
g='hsa-mir-143'
phe$gene=ifelse(exprSet[g,]>median(exprSet[g,]),'high','low')
table(phe$gene)
ggsurvplot(survfit(Surv(time, event)~gene, data=phe), conf.int=F, pval=TRUE)


## 批量生存分析 使用  logrank test 方法
mySurv=with(phe,Surv(time, event))
log_rank_p <- apply(exprSet , 1 , function(gene){
  # gene=exprSet[1,]
  phe$group=ifelse(gene>median(gene),'high','low')  
  data.survdiff=survdiff(mySurv~group,data=phe)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  return(p.val)
})
log_rank_p=sort(log_rank_p)
head(log_rank_p)
boxplot(log_rank_p)  
table(log_rank_p<0.01)
log_rank_p[log_rank_p<0.01]

library(pheatmap)
choose_gene=names(log_rank_p[log_rank_p<0.01])
choose_matrix=expr[choose_gene,]
choose_matrix[1:4,1:4]
choose_matrix=t(scale(t(log2(choose_matrix+1)))) 
## http://www.bio-info-trainee.com/1980.html
annotation_col = data.frame( group_list=group_list  )
rownames(annotation_col)=colnames(expr)
pheatmap(choose_matrix,show_colnames = F,annotation_col = annotation_col,
         filename = 'logRank_genes.heatmap.png')

library(ggfortify)
df=as.data.frame(t(choose_matrix))
df$group=group_list
png('logRank_genes.pca.png',res=120)
autoplot(prcomp( df[,1:(ncol(df)-1)] ), data=df,colour = 'group')+theme_bw()
dev.off()





