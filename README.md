## Ensembl_pairwise_alignment_20220927

利用ENSEMBL网站提供的API批量进行pairwise_alignment寻找某个序列的在其他物种的同源序列

使用:

`get_pairwise_alignment_result(human_reigon = "8:130248522-130248608:1",output_path = "./chr8_130248522-130248608/")`



## **GEO_data_analysis**

函数实现对GEO数据进行下载(表达矩阵临床信息平台注释文件)ID转换PCA差异分析绘制火山图热图



## WGCNA

实现WGCNA分析

`start.R` 初始化参数

`wgcna_function.R` WGCNA每个步骤写成函数

`run_WGCNA.R` 运行





## 富集分析

`run_go.R`

`run_kegg.R`

`run_reactome.R`

批量进行富集分析生成xlsx表格和气泡图等





## 获取基因注释信息

`get_gene_go_anno.R`

`get_gene_anno.R`