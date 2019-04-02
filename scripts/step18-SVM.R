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

### https://github.com/jmzeng1314/ML

rm(list=ls())
load(file = '../Rdata/TCGA_KIRC_mut.Rdata')
load(file = '../Rdata/TCGA-KIRC-miRNA-example.Rdata')
group_list=ifelse(as.numeric(substr(colnames(expr),14,15)) < 10,'tumor','normal')

table(group_list)
load(file='../Rdata/survival_input.Rdata')

head(phe)
exprSet[1:4,1:4]

dim(expr)
set.seed(1234567890)
#sample(1:10,3)
k=sample(1:593,300)

t_exp=expr[,k] ##
v_exp=expr[,-k]## 

table(ifelse(substr(colnames(t_exp),14,15)=='01','tumor','normal'))
table(ifelse(substr(colnames(v_exp),14,15)=='01','tumor','normal'))

t_tumor=t_exp[,substr(colnames(t_exp),14,15)=='01']
v_tumor=t_exp[,substr(colnames(v_exp),14,15)=='01']

t_phe= phe[match(substr(colnames(t_tumor),1,12),phe$ID),]
v_phe=phe[match(substr(colnames(v_tumor),1,12),phe$ID),]
table(t_phe$stage)
table(v_phe$stage)


if(F){
  ## 切割数据 
  library(caret)
  set.seed(123456789)
  sam<- createDataPartition(phe$event, p = .5,list = FALSE)
  train <- phe[sam,]
  test <- phe[-sam,]
  #查看两组一些临床参数切割比例
  prop.table(table(train$stage))
  prop.table(table(test$stage)) 
}

head(t_phe)
t_tumor[1:4,1:4]
library(lars) 
library(glmnet) 
x=t(log2(t_tumor+1))
y=t_phe$event 

# model_lasso <- glmnet(x, y, family="binomial", nlambda=50, alpha=1)
# rf_output=randomForest(x=x, y=y,importance = TRUE, ntree = 10001, proximity=TRUE )
data=cbind(x,y)
colnames(data)[ncol(data)]='event'
data=as.data.frame(data)
data$event=as.factor(data$event)
library(e1071)
model = svm(formula = event ~ .,  # 这里的待预测的变量event是二分类变量，生与死。
            data = data,kernel = "linear")

summary(model) 
pred = predict(model, data)
table(pred,data$event)


x=t(log2(v_tumor+1))
y=v_phe$event 
data=cbind(x,y)
colnames(data)[ncol(data)]='event'
data=as.data.frame(data)
data$event=as.factor(data$event)
pred = predict(model, data)
table(pred,data$event)



