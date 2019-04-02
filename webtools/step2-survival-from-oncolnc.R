rm(list=ls())
options(stringsAsFactors = F)
# http://www.oncolnc.org/kaplan/?lower=50&upper=50&cancer=LGG&gene_id=93663&raw=ARHGAP18&species=mRNA
a=read.table('LGG_93663_ARHGAP18_50_50.csv',header = T,sep = ',',fill = T)
colnames(a)
dat=a
library(ggstatsplot)
ggbetweenstats(data =dat, x = Group,  y = Expression)
library(ggplot2)
library(survival)
library(survminer) 
table(dat$Status)
dat$Status=ifelse(dat$Status=='Dead',1,0)
sfit <- survfit(Surv(Days, Status)~Group, data=dat)
sfit
summary(sfit)
ggsurvplot(sfit, conf.int=F, pval=TRUE)
ggsave('survival_ARHGAP18_in_LGG.png')
ggsurvplot(sfit,palette = c("#E7B800", "#2E9FDF"),
           risk.table =TRUE,pval =TRUE,
           conf.int =TRUE,xlab ="Time in months", 
           ggtheme =theme_light(), 
           ncensor.plot = TRUE)
ggsave('survival_ARHGAP18_in_LGG.png')
