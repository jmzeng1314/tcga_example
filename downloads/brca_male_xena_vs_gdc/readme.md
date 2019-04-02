# 比较一下从TCGA的GDC官网下载的RNA-seq的counts值和UCSC的XENA下载的区别

### 首先根据视频在TCGA的GDC官网下载所有男性乳腺癌患者的RNA-seq的counts值文件

如下:

```
.
├── [  96]  0077ddf1-c0ed-485e-87b2-a7f3a7133a86
│   └── [251K]  dfe06359-c16a-4959-bb04-298f18474a17.htseq.counts.gz
├── [  96]  00a9fc98-9f19-4a73-ae62-048fd7923b66
│   └── [247K]  e30e547e-eeff-4a1f-829b-b6ec9e79f02a.htseq.counts.gz
├── [  96]  0c00a19b-435a-4744-b142-feb64722abc1
│   └── [251K]  d998e259-4ea6-428d-8c47-4bc4656dbaa7.htseq.counts.gz
├── [  96]  50a55f2d-2edc-47eb-b6de-16d467b11d04
│   └── [249K]  0cb461ea-1a40-4f29-97d8-9337a74bc6c7.htseq.counts.gz
├── [  96]  6ec1d9a1-6c71-481c-ae53-a09c1c895029
│   └── [244K]  fcdaa827-e263-46e0-8a83-206fa5a009ba.htseq.counts.gz
├── [  96]  70702a9c-a74f-4703-9dd5-24a124ebc2d5
│   └── [248K]  f6b9d070-6d11-40f3-9ca0-a20779091825.htseq.counts.gz
├── [  96]  7a937bcc-decb-41e6-88e3-c50f21c34773
│   └── [240K]  6068fbcf-7836-41f6-8cd2-4cffde58d00a.htseq.counts.gz
├── [  96]  866cabeb-a740-4ba3-a680-dcface1fbafc
│   └── [248K]  b40a318f-fc61-49e5-ab41-c2e41510a12c.htseq.counts.gz
├── [  96]  ac333fd4-b9ac-4861-b994-152c765068b3
│   └── [248K]  1414d172-710a-4fb6-af82-f8b0c923b259.htseq.counts.gz
├── [  96]  c646e624-136a-4a5a-a11d-a47c3a398751
│   └── [243K]  06a8ba43-c1e0-4de1-9073-4b16c95245ca.htseq.counts.gz
├── [  96]  c7741a5b-70b8-41ad-87ea-10270d2c747d
│   └── [249K]  1a33a98d-2df8-4c2d-8042-5d2add5ce54c.htseq.counts.gz
├── [  96]  c9f65b52-9aae-4105-a53c-7018abebde66
│   └── [246K]  f5aa7410-7a2e-4267-ba3a-2a6a391f6b45.htseq.counts.gz
└── [  96]  e7e698aa-f2ee-4ef0-9f25-3a30a20f1d75
    └── [244K]  9838a921-f1fc-434a-abe3-d396e32c32b5.htseq.counts.gz

13 directories, 13 files

```

可以看到是13个样本的表达矩阵，而且是gz压缩的，每一个病人的表达矩阵一个文件夹所以可以简单使用R代码合并一下：

```r
fs=list.files(pattern = '*.htseq.counts.gz',recursive = T)
fs
tmp = lapply(fs, function(x){
  read.table(x) ### fread: data.table 
})
tmp=do.call(cbind,tmp)
```

这样合并的表达矩阵是缺乏样本ID的，没办法去跟其它数据对应上，其实是因为没有同步下载json文件，临床信息文件。

直接看视频吧。

### 然后看视频下载UCSC的XENA的所有男性乳腺癌患者的RNA-seq的counts值文件

这个时候下载的就是一个全部的表达矩阵，其实是可以根据上面TCGA的GDC官网挑选到ID来进行筛选。

```r
a=read.table('TCGA-BRCA.htseq_counts.tsv.gz',
             header = T,sep = '\t')
a[1:4,1:4]
rownames(a)=a[,1]

brca_male_xena=a[,match(gsub('-','.',colnames(brca_male_gdc)),
                        substring(colnames(a),1,12))]
brca_male_xena=2^brca_male_xena -1 
```

最后得到的两个表达矩阵可以看到是一模一样的。

所以，两种方法得到了印证。

