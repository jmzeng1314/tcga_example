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

## 绘图集大成者。
dev.off()
rm(list=ls())
library(survival)
library(survminer)
# 前面使用maftools解析了突变数据，就可以知道各个样本是否有某基因的突变。
load(file = '../Rdata/TCGA_KIRC_mut.Rdata')
load(file = '../Rdata/TCGA-KIRC-miRNA-example.Rdata')
group_list=ifelse(as.numeric(substr(colnames(expr),14,15)) < 10,'tumor','normal')

table(group_list)
load(file='../Rdata/survival_input.Rdata')
## 挑选感兴趣的miRNA来画表达差异的boxplot
# 2015-TCGA-ccRCC-5-miRNAs-signatures
# Integrated genomic analysis identifies subclasses and prognosis signatures of kidney cancer
# miR-21,miR-143,miR-10b,miR-192,miR-183
# hsa-mir-21,hsa-mir-143,hsa-mir-10b,hsa-mir-192,hsa-mir-183
dim(exprSet)
dim(expr)
expr[1:4,1:4]
head(phe)
head(mut)
dat=data.frame(gene=log2(exprSet['hsa-mir-10b',]+1),
               stage=phe$stage)
save(dat,file = 'data_for_boxplot.Rdata')
head(dat)
boxplot(dat$gene~dat$stage)
if(require('ggpubr')){
  library(ggpubr)
  # google search : ggpubr boxplot add p-value
  # http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/76-add-p-values-and-significance-levels-to-ggplots/
  p <- ggboxplot(dat, x = "stage", y = "gene",
                 color = "stage", palette = "jco",
                 add = "jitter")
  #  Add p-value
  p + stat_compare_means()
}
if(require('ggstatsplot')){
  library(ggstatsplot)
  ggbetweenstats(data = dat, x = stage,  y = gene)
}

# Change method
# Compute the analysis of variance
res.aov <- aov(gene ~ stage, data = dat)
# Summary of the analysis
summary(res.aov)
TukeyHSD(res.aov)


tail(sort(table(mut$Hugo_Symbol)))
VHL_mut=substr(as.character(
  as.data.frame( mut[mut$Hugo_Symbol=='VHL','Tumor_Sample_Barcode'])[,1] ),
  1,12)

library(dplyr)
  mut  %>% 
  filter(Hugo_Symbol=='VHL')  %>%   
  as.data.frame()  %>% 
  pull(Tumor_Sample_Barcode)   %>%  
  as.character()   %>%   
  substr(1,12)
  

dat=data.frame(gene=log2(exprSet['hsa-mir-10b',]),
               mut= substr(colnames(exprSet),1,12) %in% VHL_mut)
head(dat)
if(require('ggpubr')){
  library(ggpubr)
  # google search : ggpubr boxplot add p-value
  # http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/76-add-p-values-and-significance-levels-to-ggplots/
  p <- ggboxplot(dat, x = "mut", y = "gene",
                 color = "mut", palette = "jco",
                 add = "jitter")
  #  Add p-value
  p + stat_compare_means(method = "t.test")
}

if(require('ggstatsplot')){
  library(ggstatsplot)
  ggbetweenstats(data = dat, x = mut,  y = gene)
}


if(require('ggplot2')){
  library(ggplot2)
  ggplot(dat,aes(x=mut,y=gene))+
    geom_boxplot()+geom_jitter()+geom_violin()+
    theme_bw()
}








