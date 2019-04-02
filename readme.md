### TCGA实战大全

首先需要自行根据我在生信技能树平台发布的系列教程来了解TCGA基础知识，需要至少14个小时的持续学习，目录见：[TCGA基础知识传送门](tgca-introduction.md)   

如果需要视频讲解，欢迎购买我的网易云课程：https://study.163.com/course/introduction/1006067243.htm (如无必要，请勿购买，谢谢理解)

### TCGA数据的探索最基本的就是3个需求：

- 根据各种指标(某基因突变与否，肿瘤分期)把样本分组来比较感兴趣基因的表现（表达，突变，甲基化）情况。
- 使用统计学方法看某个感兴趣基因的**重要性**，比如生存分析，差异分析等等。
- 看某两个感兴趣基因的相关性，调控或者其它。

### KIRC的miRNA实战

首先需要了解TCGA计划中的KIRC这个癌症背景知识，见[PPT](PPT-TCGA-ccRCC.pdf)

然后需要通读我们本次实战所**需要复现的文章**：[Integrated genomic analysis identifies subclasses and prognosis signatures of kidney cancer.](https://www.ncbi.nlm.nih.gov/pubmed/25826081)  该文章并没有任何特殊之处，纯粹是举个例子，这样类似的文章**多达3000篇。**

通过文章我们了解到了实现一个TCGA数据挖掘的**基本步骤**：

- 下载对应的TCGA数据，主要是根据癌症种类选择**6种数据**，比如KIRC的clinical和miRNA数据，这里有**8个数据中心**供选择。
- 把病人队列分成训**练集和测试集**，然后可能需要在GEO数据库也同步查找可供挖掘数据
- 然后**走一波统计分析**，比如差异分析，生存分析，lasso回归，随机森林等等找到目标基因集
- 接着**一波可视化说明**找到的基因集具有明显的意义，包括森林图，热图，火山图等等
- 对最后的基因集计算得到预测风险的公式，还有可视化展现风险因子关联情况。





### TCGA高阶分析

主要是针对TCGA的全部类型数据，包括：

- DNA Sequencing（包括全基因组和全外显子组的maf格式somatic突变数据）
- miRNA Sequencing （表达矩阵）
- Protein Expression（表达矩阵）
- mRNA Sequencing（测序的表达矩阵）
- Total RNA Sequencing（表达矩阵）
- Array-based Expression（芯片的表达矩阵）
- DNA Methylation （25/450/850K的甲基化芯片或者WGBS）
- Copy Number（主要是SNP6.0芯片，还有测序后计算的拷贝数变异情况）

首先可以使用maftools等工具来可视化全基因组和全外显子组的maf格式somatic突变数据，代码是：





### 网页工具大全

多不胜数，简单列举如下：

- TCGA官方文章列表 <https://tcga-data.nci.nih.gov/docs/publications/>
- 涉及到的平台及产出的数据： <https://tcga-data.nci.nih.gov/docs/publications/tcga/platformdesign.html>
- 数据下载权限控制：<https://tcga-data.nci.nih.gov/docs/publications/tcga/accesstiers.html>
- 官网数据存放中心：<https://portal.gdc.cancer.gov/>
- <http://www.cbioportal.org/index.do>
- <https://xenabrowser.net/datapages/>
- <https://gdac.broadinstitute.org/>
- <http://www.oncolnc.org/>
- <http://gepia.cancer-pku.cn/detail.php?gene=ERBB2>
- <http://ibl.mdanderson.org/tanric/_design/basic/index.html>
- TCIA Collections  <http://www.cancerimagingarchive.net/>
- Human Protein Atlas (HPA)  <https://www.proteinatlas.org/pathology>
- <https://dcc.icgc.org/>
- <https://dcc.icgc.org/repositories>

- GEO2R
- <http://mexpress.be/>
- 基因甲基化和表达数据库MethHC: <http://methhc.mbc.nctu.edu.tw/php/index.php>
- lncRNA功能研究神器：TANRIC数据库 <http://ibl.mdanderson.org/tanric/_design/basic/index.html>
- TCGA可视化网站GEPIA ： <http://gepia.cancer-pku.cn/detail.php#>
- 免疫 The  Cancer  Immunome  Atlas ：<https://tcia.at/home>
- 生存分析：<http://www.oncolnc.org/>

重点不是介绍这些网页工具的用法，如果真正理解了TCGA计划的前因后果以及数据规律，就很容易明白网页工具的设计逻辑，更重要的是可以合理利用网页工具，在它们的基础上面使用R语言做定制化的深度分析。

