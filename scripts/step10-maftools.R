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
require(maftools) 
options(stringsAsFactors = F) 
# 理论上可以针对任何癌症的突变数据结果进行可视化和适当的统计。
laml = read.maf(maf = '../GDC/TCGA.KIRC.mutect.somatic.maf.gz')
laml 
project='TCGA_KIRC'

laml@data=laml@data[!grepl('^MT-',laml@data$Hugo_Symbol),]
#laml@data=laml@data[!grepl('^MUC',laml@data$Hugo_Symbol),]
laml@data$t_vaf = (laml@data$t_alt_count/laml@data$t_depth)
unique(laml@data$Tumor_Sample_Barcode)
getSampleSummary(laml) 
getGeneSummary(laml) 
getFields(laml)  

png(paste0('plotmafSummary_',project,'.png'),res = 150,width = 1080,height = 1080)
plotmafSummary(maf = laml, rmOutlier = TRUE,showBarcodes = T,
               addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()

## ---- fig.align='left',fig.height=5,fig.width=10, fig.align='left'-------
#We will draw oncoplots for top ten mutated genes.
png(paste0('oncoplot_top30_',project,'.png'),res = 150,width = 1080,height = 1080)
oncoplot(maf = laml, top = 30, fontSize = 12 ,showTumorSampleBarcodes = F )
dev.off()

table(laml@data$Tumor_Sample_Barcode)
 

if(F){
  oncoplot(maf = laml, top = 15, fontSize = 12,
           clinicalFeatures =c( "subtype" ) ,sortByAnnotation= T )
}

png(paste0('TMB_',project,'.png'),res = 150,width = 1080,height = 1080)
laml.mutload = tcgaCompare(maf = laml, cohortName = project)
dev.off()

png(paste0('Vaf_',project,'.png'),res = 150,width = 1080,height = 1080)
plotVaf(maf = laml ,top = 20)
dev.off()

if(F){
  
  dir.create(paste0('vaf_clust_',project ))
  table(laml@data$Tumor_Sample_Barcode)
  lapply(unique(laml@data$Tumor_Sample_Barcode), function(x){
    png(paste0('vaf_clust_',project,'/',x,'_vaf_clust.png'),res=120,width = 1080,height = 1080)
    het = inferHeterogeneity(maf = laml, tsb = x, vafCol = 't_vaf')
    print(het$clusterMeans) 
    plotClusters(clusters = het)
    dev.off()
  }) 
  
}


getFields(laml)  
laml@data$t_vaf = (laml@data$t_alt_count/laml@data$t_depth)
mut=laml@data[,c("Hugo_Symbol","Chromosome","Start_Position","Tumor_Sample_Barcode","t_vaf")]
mut$pos=paste(mut$Chromosome,mut$Start_Position,sep=':')

tail(sort(table(mut$Hugo_Symbol)))

save(mut,file = '../Rdata/TCGA_KIRC_mut.Rdata')












