


# @description: 对基因表达矩阵进行id转换 转换后重复的gene name 通过 avereps 取均值

# @parameters: exp_data 基因表达矩阵
# @parameters: ids 两列的数据框 一列为转换的ID 一列为gene name
# @parameters: join_by 要转换的ID的列名

# @return: exp_data_convert_id

# @example: exp_data_last <- id2symbol(exp_data = exp_dat,ids = ids,join_by = "ensembl_gene_id")


id2symbol <- function(exp_data,ids,join_by){
  
  library(dplyr)
  library(limma)
  
  exp_data[,join_by] <- rownames(exp_data)
  
  exp_data<- 
    exp_data %>% left_join(ids) %>% dplyr::select(-.data[[join_by]])
  
  print(exp_data[,ncol(exp_data)] %>% table() %>% table() )
  
  
  gene_symbol <- exp_data[,ncol(exp_data),drop = T]
  
  
  exp_data<- 
    exp_data %>% dplyr::select(-ncol(exp_data)) %>% as.matrix()
  
  rownames(exp_data) <- gene_symbol
  
  
  
  
  exp_data_convert_id <- avereps(exp_data) %>% as.data.frame()
  
  return(exp_data_convert_id)
  
}