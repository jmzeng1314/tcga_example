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

table(group_list)
exprSet=na.omit(expr)
dim(exprSet)
load(  file = 
         file.path(Rdata_dir,'TCGA-KIRC-miRNA-survival_input.Rdata')
)
dim(exprSet) ## remove the nomral
head(phe)
exprSet[1:4,1:4]
head(colnames(exprSet))
head(phe$ID)
## 必须保证生存资料和表达矩阵，两者一致
all(substring(colnames(exprSet),1,12)==phe$ID)
 

library(randomForest)
library(ROCR)
library(genefilter)
library(Hmisc)
head(colnames(exprSet))
head(phe$ID)
## 必须保证生存资料和表达矩阵，两者一致
all(substring(colnames(exprSet),1,12)==phe$ID)

x=t(log2(exprSet+1))
y=phe$event
table(y)

tmp = as.vector(table(y))
num_classes = length(tmp)
min_size = tmp[order(tmp,decreasing=FALSE)[1]]
sampsizes = rep(min_size,num_classes)
tmp_rf=file.path(Rdata_dir,'TCGA_KIRC_miRNA_rf_output.Rdata')

if(!file.exists(tmp_rf)){
  rf_output=randomForest(x=x, y=y,importance = TRUE, ntree = 10001, proximity=TRUE )
  save(rf_output,file = tmp_rf)
}
load(file = tmp_rf)
rf_output
str(rf_output)


rf_importances=importance(rf_output, scale=FALSE)
head(rf_importances)
varImpPlot(rf_output, type=2, n.var=30, scale=FALSE, 
           main="Variable Importance (Gini) for top 30 predictors",cex =.7)
target_labels=as.vector(y)
MDSplot(rf_output, y, k=2, xlab="", ylab="", 
        pch=target_labels, palette=c("red", "blue"), main="MDS plot")

print(rf_output) 
choose_gene=rownames(tail(rf_importances[order(rf_importances[,2]),],50))

library(pheatmap) 
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
         filename = 'rf_genes.heatmap.png')

library(ggfortify)
df=as.data.frame(t(choose_matrix))
df$group=group_list
png('rf_genes.pca.png',res=120)
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









