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

laml = read.maf(maf = '../GDC/TCGA.KIRC.mutect.somatic.maf.gz')
laml 
project='TCGA_KIRC'

suppressPackageStartupMessages(library("deconstructSigs"))
suppressPackageStartupMessages(library("BSgenome"))
suppressPackageStartupMessages(library("BSgenome.Hsapiens.UCSC.hg38"))
options(stringsAsFactors = F)
mut=laml@data
head(mut)
mut=mut[mut$Variant_Type=='SNP',]
a=mut[,c(16,5,6,12,13)]
colnames(a)=c( "Sample","chr", "pos","ref",  "alt")

a$Sample=as.character(a$Sample)

plot(table(a$Sample),las=2)
sigs.input <- mut.to.sigs.input(mut.ref = a, 
                                sample.id = "Sample", 
                                chr = "chr", 
                                pos = "pos", 
                                ref = "ref", 
                                alt = "alt",
                                bsg = BSgenome.Hsapiens.UCSC.hg38)
class(sigs.input)
barplot(as.numeric(sigs.input[1,]))
head(t(sigs.input))
w=lapply(unique(a$Sample)[1:10], function(i){
  ## signatures.cosmic signatures.nature2013
  sample_1 = whichSignatures(tumor.ref = sigs.input[,], 
                             signatures.ref = signatures.cosmic, 
                             sample.id =  i, 
                             contexts.needed = TRUE,
                             tri.counts.method = 'default')
  print(i)
  return(sample_1$weights)
})
w=do.call(rbind,w)
library(pheatmap)
pheatmap(t(w),cluster_rows = F,cluster_cols = T)
pheatmap(w,cluster_rows = T,cluster_cols = F)

# Determine the signatures contributing to the two example samples
lapply(unique(a$Sample), function(i){
  ## signatures.cosmic signatures.nature2013
  sample_1 = whichSignatures(tumor.ref = sigs.input, 
                             signatures.ref = signatures.cosmic, 
                             sample.id =  i, 
                             contexts.needed = TRUE,
                             tri.counts.method = 'default')
  pdf(paste0(i,'.sig.pdf'))
  plotSignatures(sample_1, sub = i)
  dev.off()
})


