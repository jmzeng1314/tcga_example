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

# shell代码大家可以不要学，有点难度。

rm(list=ls())
## for human gencode.v25.annotation.gtf
# cat gencode.v25.annotation.gtf|perl -alne '{next unless $F[2] eq "gene";print}'|grep -w HAVANA |\
# cut -f 1,4,5,9| cut -d";" -f 1,2,4|sed 's/gene_id//g'|sed 's/gene_type//g'|sed 's/gene_name//g'|\
# sed 's/;//g'| sed 's/\"//g'|perl -alne '{/(ENSG\d+)/;print "$1\t$_"}' >human.gene.positions
if(F){
  a=read.table('data/human.gene.positions')[,c(2:4,1,6,7)]
  colnames(a)=c('chr','start','end','ensembl','type','symbol')
  length(unique(a$symbol))
  length(unique(a$ensembl))
  head(a)
  human_geneInfo_genecode_v25=a 
}
load('../human_geneInfo_genecode_v25.rda')
b=human_geneInfo_genecode_v25
head(b)
tail(sort(table(b$type)))
lincRNA_genes=b[b$type %in% c('lincRNA','antisense'),'ensembl']
a=read.table('../Rdata/TCGA-UVM.htseq_counts.tsv.gz',header = T,sep = '\t')
dim(a)
a[1:4,1:4]
library(stringr)
a$Ensembl_ID=str_split(a$Ensembl_ID,'[.]',simplify = T)[,1]

lincRNA=a[a$Ensembl_ID %in% lincRNA_genes,]







