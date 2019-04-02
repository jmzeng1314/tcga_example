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

### https://github.com/jmzeng1314/GEO/blob/master/GSE11121/step5-surivival.R

rm(list=ls())
options(stringsAsFactors = F)

Rdata_dir='../Rdata/'
Figure_dir='../figures/'
# 加载上一步从RTCGA.miRNASeq包里面提取miRNA表达矩阵和对应的样本临床信息。
load( file = 
        file.path(Rdata_dir,'TCGA-KIRC-miRNA-example.Rdata')
)
dim(expr)
dim(meta)
# 可以看到是 537个病人，但是有593个样本，每个样本有 552个miRNA信息。
# 当然，这个数据集可以下载原始测序数据进行重新比对，可以拿到更多的miRNA信息

# 这里需要解析TCGA数据库的ID规律，来判断样本归类问题。
group_list=ifelse(as.numeric(substr(colnames(expr),14,15)) < 10,'tumor','normal')

table(group_list)
exprSet=na.omit(expr)

library(survival)
library(survminer)

# 这里做生存分析，已经不需要正常样本的表达矩阵了，所以需要过滤。
# 而且临床信息，有需要进行整理。
### survival analysis only for patients with tumor.
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
  
  save(exprSet,phe,
       file = 
         file.path(Rdata_dir,'TCGA-KIRC-miRNA-survival_input.Rdata')
  )
}
# 上面被关闭的代码，就是在整理临床信息和生存分析的表达矩阵。
# 整理好的数据，直接加载即可
load(  file = 
         file.path(Rdata_dir,'TCGA-KIRC-miRNA-survival_input.Rdata')
)
head(phe)
exprSet[1:4,1:4]


library(survival)
library(survminer)

## 批量生存分析 使用 coxph 回归方法
# http://www.sthda.com/english/wiki/cox-proportional-hazards-model
colnames(phe)
mySurv=with(phe,Surv(time, event))
# 对五百多个miRNA基因进行批量运行cox，需要一点点时间。
# 如果是mRNA-seq的表达矩阵， 通常耗时更长。
# 注意，如果是某些基因表达量为恒定，比如在所有样本为0，这个代码会爆仓
# 需要去除这样的基因，没有分析的必要性。

cox_results <-apply(exprSet , 1 , function(gene){
  # gene= exprSet[1,]
  group=ifelse(gene>median(gene),'high','low') 
  survival_dat <- data.frame(group=group,stage=phe$stage,age=phe$age,
                             gender=phe$gender,
                             stringsAsFactors = F)
  m=coxph(mySurv ~ gender + age + stage+ group, data =  survival_dat)
  
  beta <- coef(m)
  se <- sqrt(diag(vcov(m)))
  HR <- exp(beta)
  HRse <- HR * se
  
  #summary(m)
  tmp <- round(cbind(coef = beta, se = se, z = beta/se, p = 1 - pchisq((beta/se)^2, 1),
                     HR = HR, HRse = HRse,
                     HRz = (HR - 1) / HRse, HRp = 1 - pchisq(((HR - 1)/HRse)^2, 1),
                     HRCILL = exp(beta - qnorm(.975, 0, 1) * se),
                     HRCIUL = exp(beta + qnorm(.975, 0, 1) * se)), 3)
  return(tmp['grouplow',])
  
})
cox_results=t(cox_results)
table(cox_results[,4]<0.05)
cox_results[cox_results[,4]<0.05,]

## 批量生存分析 使用  logrank test 方法
mySurv=with(phe,Surv(time, event))
log_rank_p <- apply(exprSet , 1 , function(gene){
  # gene=exprSet[1,]
  phe$group=ifelse(gene>median(gene),'high','low')  
  data.survdiff=survdiff(mySurv~group,data=phe)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  return(p.val)
})

require("VennDiagram")
VENN.LIST=list(cox=rownames(cox_results[cox_results[,4]<0.05,]),
             log=names(log_rank_p[log_rank_p<0.05]))
venn.plot <- venn.diagram(VENN.LIST , NULL, 
                          fill=c("darkmagenta", "darkblue"), 
                          alpha=c(0.5,0.5), cex = 2, 
                          cat.fontface=4,  
                          main="overlap of coxph and log-rank test")
# 可以看到两种生存分析挑出来的基因重合度尚可。
grid.draw(venn.plot) 

save(log_rank_p,cox_results ,
     file = 
       file.path(Rdata_dir,'TCGA-KIRC-miRNA-survival_results.Rdata')
)
 
library(pheatmap)
choose_gene=rownames(cox_results[cox_results[,4]<0.05,])
choose_matrix=expr[choose_gene,]
choose_matrix[1:4,1:4] 
n=t(scale(t(log2(choose_matrix+1))))  #scale()函数去中心化和标准化
#对每个探针的表达量进行去中心化和标准化
n[n>2]=2 #矩阵n中归一化后，大于2的项，赋值使之等于2（相当于设置了一个上限）
n[n< -2]= -2 #小于-2的项，赋值使之等于-2（相当于设置了一个下限）
n[1:4,1:4]

## http://www.bio-info-trainee.com/1980.html
annotation_col = data.frame( group_list=group_list  )
rownames(annotation_col)=colnames(expr)

pheatmap(n,show_colnames = F,annotation_col = annotation_col,
         filename = 'cox_genes.heatmap.png')

library(ggfortify)
df=as.data.frame(t(choose_matrix))
df$group=group_list
png('cox_genes.pca.png',res=120)
autoplot(prcomp( df[,1:(ncol(df)-1)] ), data=df,colour = 'group')+theme_bw()
dev.off()

## 也可以尝试其它主成分分析的R包，视频就不继续没完没了的讲解了。


library("FactoMineR")
library("factoextra")  
## 这里的PCA分析，被该R包包装成一个简单的函数，复杂的原理后面讲解。
dat.pca <- PCA(t(choose_matrix), graph = FALSE) #'-'表示“非”
fviz_pca_ind(dat.pca,repel =T,
             geom.ind = "point", # show points only (nbut not "text")只显示点不显示文本
             col.ind =  group_list, # color by groups 颜色组
             # palette = c("#00AFBB", "#E7B800"),
             addEllipses = TRUE, # Concentration ellipses 集中成椭圆
             legend.title = "Groups"
)









