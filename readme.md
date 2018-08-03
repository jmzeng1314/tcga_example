# 一个TCGA数据挖掘实例

### 首先需要了解TCGA计划背景知识及数据下载方法汇总

TCGA的28篇教程往期目录如下：

[使用R语言的cgdsr包获取TCGA数据](http://mp.weixin.qq.com/s?__biz=MzAxMDkxODM1Ng==&mid=2247486492&idx=1&sn=3a7251244377fdd4b2a3aa5c8cd1131a&chksm=9b484ca7ac3fc5b1a21202cf25ff15a8eec434424aa3e48787129fa6f5e66ebe57ffcb631772&scene=21#wechat_redirect)  （cBioPortal）

[TCGA的28篇教程- 使用R语言的RTCGA包获取TCGA数据](http://mp.weixin.qq.com/s?__biz=MzAxMDkxODM1Ng==&mid=2247486585&idx=1&sn=3035f6420904aad2c8161b362cdeb472&chksm=9b484cc2ac3fc5d479fc5bce3d68d4666b763652a21a55b281aad8c0c4df9b56b4d3b353cc4c&scene=21#wechat_redirect) （离线打包版本, FireBrowse）

[TCGA的28篇教程- 使用R语言的RTCGAToolbox包获取TCGA数据](http://mp.weixin.qq.com/s?__biz=MzAxMDkxODM1Ng==&mid=2247486728&idx=1&sn=3990dff5efccedc060443b7f3af3b6ee&chksm=9b484db3ac3fc4a51ee34ba578280d89ea6159a48d9dec1c7ecd5dbe208c53a1b36c439de75d&scene=21#wechat_redirect) （Broad Institute FireBrowse portal）

[ TCGA的28篇教程-  批量下载TCGA所有数据](http://mp.weixin.qq.com/s?__biz=MzAxMDkxODM1Ng==&mid=2247486746&idx=1&sn=b7c5ad7eff8cffb3620756f5feaff587&chksm=9b484da1ac3fc4b741a6e3b59ba1bf668a11e21eb610f1a1d4582e33d429c67c14e6659c0771&scene=21#wechat_redirect) （ UCSC的 XENA）

因为TCGA实在是一个跨时代的癌症研究项目，不能在下载这个基础环境耽误太多的功夫，下载渠道再多，也只需要一个好用的即可！

我以前在生信技能树论坛也写过 TCGA数据下载合集：

[以前下载TCGA数据非常简单，都在一个远程电脑里面](http://www.biotrainee.com/thread-820-1-1.html)

[现在下载TCGA数据也是非常方便，首先是GDC网站及客户端](http://www.biotrainee.com/thread-821-1-1.html) (基于mysql数据库进行条件过滤)

[现在下载TCGA数据也是非常方便，然后是firehose网站及客户端](http://www.biotrainee.com/thread-822-1-2.html) (基于远程电脑文件夹及文件名过滤)

[现在下载TCGA数据也是非常方便，接着是cgdsR和cbioportal](http://www.biotrainee.com/thread-824-1-3.html) (基于TCGA大文章分篇下载)

[现在下载TCGA数据也是非常方便，倒数第二个是Synapse](http://www.biotrainee.com/thread-825-1-3.html) (基于作者整理TCGA数据上传)

[现在下载TCGA数据也是非常方便，最后是各种杂七杂八的工具](http://www.biotrainee.com/thread-826-1-3.html) 

 下面就简单罗列几个还算是比较流行的TCGA下载器吧：

#### GDC官方下载工具

GDC给出了一系列的用户友好的选择框，你只需要根据条条框框来选择就可以下载到自己想要的数据，而不需要去几百个文件夹里面漫无目的的查找了。 <https://gdc-portal.nci.nih.gov/legacy-archive/search/f>  根据自定义搜索过滤条件拿到了 mainfest 文件就可以啦。

GDC客户端的说明书是：[https://docs.gdc.cancer.gov/Data ... de/Getting_Started/](https://docs.gdc.cancer.gov/Data_Transfer_Tool/Users_Guide/Getting_Started/) 
傻瓜式软件，非常简单！
[https://docs.gdc.cancer.gov/Data ... ownload_and_Upload/](https://docs.gdc.cancer.gov/Data_Transfer_Tool/Users_Guide/Data_Download_and_Upload/)
一般人只需要根据你搜索过滤得到的mainfest进行GDC下载数据即可，下载下来的文件，是每个样本一个文件夹，需要合并，需要了解为什么用XML来存储信息

### 实例文章



### 统计学背景

需要掌握的知识较多，包括：

- PCA分析
- LASSO回归
- COXPH