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

### https://github.com/jmzeng1314/ML

rm(list=ls())
library(survival)
library(survminer)
load(file = 'TCGA-KIRC-miRNA-example.Rdata')
group_list=ifelse(substr(colnames(expr),14,15)=='01','tumor','normal')
table(group_list)
load(file='survival_input.Rdata')
head(phe)
exprSet[1:4,1:4]


# https://rpubs.com/xuefliang/153247 
# https://rpubs.com/kaz_yos/survival-auc
# https://cran.r-project.org/web/views/Survival.html

# 诊断试验中ROC曲线的绘制，一般金标准都是二分类变量，比如有病与无病。
# ROC曲线以假阳性率（1-speficicity）作为横坐标，以真阳性率（sensitivity）作为纵坐标。
# 曲线上的点代表不同临界点所对应的灵敏度和特异度对子。
# 通常将ROC曲线左上角那一点定为最佳临界点，此点的Youden指数（sensitivity-(1-specificity)）最大。
# 如果金标准是生存分析资料（生存时间overall survival time与生存状态 status）
# 需要通过R软件的survivalROC包来介绍生存资料的ROC分析，即时间依赖的ROC分析。


library(timeROC)
new_dat$fp=fp
with(new_dat,
     ROC <<- timeROC(T=time,#结局时间 
                     delta=event,#生存结局 
                     marker=fp,#预测变量 
                     cause=1,#阳性结局赋值，比如死亡与否
                     weighting="marginal",#权重计算方法，marginal是默认值，采用km计算删失分布
                     times=c(60,100),#时间点，选取10年和20年生存率
                     ROC = TRUE,
                     iid = TRUE)
)
# 画roc曲线
plot(ROC,time=60,col = "blue",add =FALSE)
#time是时间点，col是线条颜色，add指是否添加在上一张图中
plot(ROC,time=100,col = "red",add = F)

ROC$AUC
confint(ROC)





