install.packages("rjson")
# Load the package required to read JSON files.
library("rjson")

# Give the input file name to the function.
result <- fromJSON(file = "files.2018-11-14.json")

# Print the result.
print(result)

f=unlist(lapply(result, function(x){x$file_name}))
case_id=unlist(lapply(result, function(x){x$cases[[1]]$case_id}))
meta=data.frame(f=f,case_id=case_id)
a=read.table('clinical.tsv',sep = '\t',header = T)[,1:2]
meta=merge(meta,a,by='case_id')

fs=list.files(pattern = '*.htseq.counts.gz',recursive = T)
fs
tmp = lapply(fs, function(x){
  read.table(x) ### fread: data.table 
})
tmp=do.call(cbind,tmp)
rownames(tmp)=tmp[,1]
tmp=tmp[,seq(2,ncol(tmp),by = 2)]
fn=unlist(lapply(strsplit(fs,'/'), function(x)x[3]))
fn
colnames(tmp)=meta[match(fn,meta$f),3]
brca_male_gdc=tmp
write.csv(brca_male_gdc,'brca_male.csv')






