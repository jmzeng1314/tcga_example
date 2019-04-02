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

## 绘图集大成者
dev.off()
rm(list=ls())
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
dat=data.frame(gene1=log2(exprSet['hsa-mir-10b',]+1),
               gene2=log2(exprSet['hsa-mir-143',]+1),
               stage=phe$stage)
save(dat,file = 'for_scatter.Rdata')
library(ggpubr)
# google search : ggpubr boxplot add p-value
# http://www.sthda.com/english/rpkgs/ggpubr/reference/stat_cor.html 
dat$stage=as.factor(dat$stage)
sp <- ggscatter(dat, x = "gene1", y = "gene2",
                add = "reg.line",  # Add regressin line 
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE # Add confidence interval
) 
# Add correlation coefficient
sp
sp + stat_cor(method = "pearson", label.x = 15, label.y = 20)

# Color by groups and facet
#::::::::::::::::::::::::::::::::::::::::::::::::::::
sp <- ggscatter( dat, x = "gene1", y = "gene2",
                color = "stage", palette = "jco",
                add = "reg.line", conf.int = TRUE)
sp + stat_cor(aes(color = stage),label.x = 15 )


