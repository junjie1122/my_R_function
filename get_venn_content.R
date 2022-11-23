

# @description: 获取各个集合的内容

# @parameters: gene_list 基因向量构成的list,list要有name

# @return: all_intersect a list 包含所有可能的集合情况和集合内容

get_venn_content <- function(gene_list){
  
  
  library(VennDiagram)
  library(dplyr)
  all_intersect <- list()
  list_nums <- length(gene_list)
  df <- get.venn.partitions(gene_list)
  
  for ( i in 1:nrow(df)) {
    
    set <- colnames(df)[1:list_nums][unlist(df[i,1:4,drop = F])]
    set_name <- paste(set,collapse = "+")
    
    if (df[i,"..values.."] %>% unname() %>% unlist() %>% length() == 0) {
      
      #all_intersect[[set_name]] <- NULL
      next
    }
    all_intersect[[set_name]] <- df[i,"..values.."] %>% unname() %>% unlist()
  }
  
  
  return(all_intersect)
  

}
