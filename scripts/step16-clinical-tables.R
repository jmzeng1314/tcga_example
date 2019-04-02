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

Rdata_dir='../Rdata/'
Figure_dir='../figures/'
f=file.path(Rdata_dir,'TCGA-LUAD-phe_clinical_tables.Rdata')


if(!file.exists(f)){
  library(RTCGA.clinical) 
  ??RTCGA.clinical
  meta <- LUAD.clinical
  tmp=as.data.frame(colnames(meta))
  meta[(grepl('patient.bcr_patient_barcode',colnames(meta)))]
  meta[(grepl('patient.days_to_last_followup',colnames(meta)))]
  meta[(grepl('patient.days_to_death',colnames(meta)))]
  meta[(grepl('patient.vital_status',colnames(meta)))]
  head(meta[(grepl('tnm',colnames(meta)))])
  tnm=meta[(grepl('tnm',colnames(meta)))]
  tnm=tnm[,4:6]
  ## patient.race  # patient.age_at_initial_pathologic_diagnosis # patient.gender 
  # patient.stage_event.clinical_stage
  clin=as.data.frame(meta[c('patient.bcr_patient_barcode','patient.vital_status',
                            'patient.days_to_death','patient.days_to_last_followup',
                            'patient.race',
                            'patient.age_at_initial_pathologic_diagnosis',
                            'patient.gender' ,
                            'patient.stage_event.pathologic_stage')])
  
  clin=cbind(clin,tnm)
  
  colnames(clin)
  clin[,3][is.na(clin[,3])]=0
  clin[,4][is.na(clin[,4])]=0
  clin$days=as.numeric(clin[,3])+as.numeric(clin[,4])
  clin=clin[,c(1:2,5:12)]
  colnames(clin)
  colnames(clin)=c('barcode','vital_status','race','age','gender','stage','m','n','t',"days") 
  library(survival)
  library(survminer)
  clin$event=ifelse(clin$vital_status=='alive',0,1)
  clin$age=as.numeric(clin$age)
  library(stringr) 
  clin$stage=str_split(clin$stage,' ',simplify = T)[,2]
  table(  clin$stage)
  boxplot(clin$age)
  clin$age_group=ifelse(clin$age>median(clin$age,na.rm = T),'older','younger')
  table(clin$race)
  clin$time=clin$days/30
  phe=clin
  
  head(phe)
  phe$barcode=toupper(phe$barcode)  
  
  save(phe,file=f)
}
load(file=f)


clinical_info=phe
head(clinical_info) 
# https://cran.r-project.org/web/packages/tableone/vignettes/introduction.html
## tableone package itself
library(tableone)
## survival pcakge for Mayo clinical_infoic's PBC data
library(survival)
data(pbc)
head(pbc)
CreateTableOne(data = pbc)

#首先对需要观测的临床特质值进行重新编码 
clinical_info$age<-as.numeric(clinical_info$age)
clinical_info$AGE<-factor(ifelse(clinical_info$age>60,'>60','<=60'),ordered = T)
clinical_info$gender<-factor(toupper(clinical_info$gender),levels=c("MALE", "FEMALE"),ordered = T)
clinical_info$stage<-factor(toupper(clinical_info$stage),ordered = T)
clinical_info$t<-factor(clinical_info$t,ordered = T)
clinical_info$n<-factor(clinical_info$n,ordered = T)
clinical_info$m<-factor(clinical_info$m,ordered = T) 
clinical_info$vital_status<-factor(toupper(clinical_info$vital_status),ordered = T)

clinical_info$race<-factor(clinical_info$race,ordered = T) 
# 去除不需要的临床信息
clinical_info=clinical_info[,-c(1,10:12)]
dput(names(clinical_info))
## Vector of variables to summarize
myVars <- dput(names(clinical_info))
## Vector of categorical variables that need transformation
catVars <- myVars[c(1,2,4:8,10)]
## ------------------------------------------------------------------------
##三线表类型之一  切割数据 
library(caret)
set.seed(123456789)
sam<- createDataPartition(clinical_info$vital_status, p = .5,list = FALSE)
train <- clinical_info[sam,]
test <- clinical_info[-sam,]
#查看两组一些临床参数切割比例
prop.table(table(train$stage))
prop.table(table(test$stage))
#添加分组
train$group<-'training datasets'
test$group<-'testing datasets'
clinical_info<-rbind(train,test)
clinical_info$group<-factor(clinical_info$group)
##生成三线表
vars <-colnames(clinical_info)[c(2:9,12,14,15)]
library(tableone)
## 最重要的三线表通常是以训练集和数据集来区分：group
tb_group<-CreateTableOne(vars = myVars, strata = c("group"), data = clinical_info,
                          factorVars = catVars) 
tab1<-print(tb_group, nonnormal = c('age','time'),
            exact = c(myVars,'AGE'), smd = TRUE)
summary(tab1)
tab_out<-print(tb_group, catDigits = 1, contDigits = 2, pDigits = 3,
           quote = FALSE, missing = T, explain = TRUE, printToggle = TRUE,
           test = TRUE, smd = T, noSpaces = FALSE, padColnames = FALSE,
           varLabels = FALSE, format = c("fp", "f", "p", "pf")[1],
           showAllLevels = FALSE, cramVars = NULL, dropEqual = FALSE,
           exact = NULL, nonnormal = NULL, minMax = FALSE)
## Save to a CSV file
write.csv(tab_out, file = "TCGA-LUAD-phe_clinical_tables1.csv")
# https://www.reed.edu/data-at-reed/resources/R/creating_html_tables.html
if(require(xtable)){
  library(xtable)
  xtable(tab_out)
}
if(require(DT)){
  file = "TCGA-LUAD-phe_clinical_tables1.html"
  y <- DT::datatable(tab_out,escape = F,rownames=F)
  DT::saveWidget(y,file)
}


# 如果要生成word表格，如果安装成功了ReporteRs和rJava
library(dplyr)
library(ReporteRs)
if(require('ReporteRs')){
  docx() %>% 
    addFlexTable(tab1 %>%
                   FlexTable( header.cell.props = cellProperties( background.color =  "#003366" ),
                              header.text.props = textBold( color = "white" ),
                              add.rownames = TRUE ) %>%
                   setZebraStyle( odd = "#DDDDDD", even = "#FFFFFF" ) ) %>%
    writeDoc( file = "clinical_infotable1.docx" )
}

### 其次也需要一个全部病人信息的临床三线表
tb_all<-CreateTableOne(vars = myVars,
                          data = clinical_info,
                          factorVars = catVars) 
tab2<-print(tb_all)
#生成临床全部信息表格
if(require('ReporteRs')){
  docx() %>% 
    addFlexTable(tab2 %>%
                   FlexTable( header.cell.props = cellProperties( background.color =  "#003366" ),
                              header.text.props = textBold( color = "white" ),
                              add.rownames = TRUE ) %>%
                   setZebraStyle( odd = "#DDDDDD", even = "#FFFFFF" ) ) %>%
    writeDoc( file = "clinical_infotable2.docx" )
}
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
##三线表类型之2 
#以性别和肿瘤分期为例
tb_gender<-CreateTableOne(vars = myVars, strata = c("gender"), data = clinical_info,
                         factorVars = catVars) 
tb_stage<-CreateTableOne(vars = myVars, strata = c("stage"), data = clinical_info,
                          factorVars = catVars)  
#打印表格
tab3<-print(tb_gender, nonnormal = c('age','time'),
            exact = c(myVars,'AGE'), smd = TRUE)
tab4<-print(tb_stage, nonnormal = c('age','time'),
            exact = c(myVars,'AGE'), smd = TRUE)
#生成word表格
if(require('ReporteRs')){
  docx() %>% 
    addFlexTable(tab3%>%
                   FlexTable( header.cell.props = cellProperties( background.color =  "#003366" ),
                              header.text.props = textBold( color = "white" ),
                              add.rownames = TRUE ) %>%
                   setZebraStyle( odd = "#DDDDDD", even = "#FFFFFF" ) ) %>%
    writeDoc( file = "clinical_infotable3.docx" )
  docx() %>% 
    addFlexTable(tab4 %>%
                   FlexTable( header.cell.props = cellProperties( background.color =  "#003366" ),
                              header.text.props = textBold( color = "white" ),
                              add.rownames = TRUE ) %>%
                   setZebraStyle( odd = "#DDDDDD", even = "#FFFFFF" ) ) %>%
    writeDoc( file = "clinical_infotable4.docx" )
}
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
