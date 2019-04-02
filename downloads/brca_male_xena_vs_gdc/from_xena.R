a=read.table('TCGA-BRCA.htseq_counts.tsv.gz',
             header = T,sep = '\t')
a[1:4,1:4]
rownames(a)=a[,1]

brca_male_xena=a[,match(gsub('-','.',colnames(brca_male_gdc)),
                        substring(colnames(a),1,12))]
brca_male_xena=2^brca_male_xena -1 
save(brca_male_xena,brca_male_gdc,
     file = 'brca_male.Rdata')

