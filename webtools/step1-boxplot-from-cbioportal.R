rm(list=ls())
options(stringsAsFactors = F)
# Ovarian Serous Cystadenocarcinoma (TCGA, Nature 2011)
# All Complete Tumors  (316 samples) / 1 Genes
# Gene Set / Pathway is altered in 28 (8.9%) of queried samples
a=read.table('plot-data-ARHGAP18-TCGA-OV-cbioportal.txt',header = T,sep = '\t',fill = T)
colnames(a)=c('id','stage','gene','mut')
dat=a
library(ggstatsplot)
ggbetweenstats(data =dat, x = stage,  y = gene)
library(ggplot2)
ggsave('plot-again-ARHGAP18-TCGA-OV-cbioportal.png')

# Change method
# Compute the analysis of variance
res.aov <- aov(gene ~ stage, data = dat)
# Summary of the analysis
summary(res.aov)
TukeyHSD(res.aov)
