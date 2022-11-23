library(tidyverse)
library(openxlsx)


getwd()


# 准备表达数据
exp_data <- read.table("./ACC/ACC.txt",header = T,row.names = 1)
exp_data
dim(exp_data)

# 准备trait数据
trait <- read.table("./ACC/design.txt",header = T) %>% t() %>% .[-1,] %>% as.data.frame()
colnames(trait) <- "group"
group <- trait$group
trait <- model.matrix(~0 + group) %>% as.data.frame %>% set_names(levels(as.factor(group)) )





# 参数设置
#-------------------------------------------------------------------------------


#运行WGCNA脚本的路径
run_WGCNA_Rscript_path <- "./WGCNA/run.R"


# 是否添加基因注释信息
gene_anno <- NULL       # 不添加设置为NULL
gene_colnames = ""
gene_type_colname = ""

# 输出路径
out_path <- "./"

# 数据标准化 数据转换 方法
method = "vst"

# 去除离群样本
sampletree_cut_height <- NULL

# 过滤 基因
filter_gene_method = "mad"
keep_gene_num = 5000
filter_threshold = NULL

# 网络构建参数设置
cor_type = "pearson"
net_type = "signed"
TOM_type = "signed"
max_blockSize = 30000
mergeCutHeight = 0.25

#核心基因筛选
GS_threshold = 0.6
MM_threshold = 0.8
#-------------------------------------------------------------------------------

# 运行WGCNA
source(run_WGCNA_Rscript_path)





