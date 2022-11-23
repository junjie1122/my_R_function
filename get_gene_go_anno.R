


# @description: 获得基因的GO注释

# @parameters: species 指定物种 可选择 "human" or "mouse"
# @parameters: gene_type 基因类型 可选择 "ENSEMBL" or "SYMBOL"
# @parameters: format 输出格式 可选择 "longer" or "wider" longer每个GO trem占据一行,wider相同基因的所有GO trem压缩为一行 

# @return: 返回指定 format的GO annotate

# @example: 
            # gene = ene_list[[1]][1:5]
            # species = "mouse"
            # gene_type = "ENSEMBL"


get_gene_go_anno <- function(gene,species = "human",gene_type = "ENSEMBL",format = "wider"){
  
  library(biomaRt)
  library(stringr)
  library(dplyr)
  
  if (!exists("ensembl")) {
    
    if (species == "mouse") {
      ensembl <- useEnsembl(biomart = "genes")
      ensembl <- useDataset(dataset = "mmusculus_gene_ensembl", mart = ensembl)
      
      print("assign ensmbl to .GlobalEnv")
      ensembl <<- ensembl
    }
    
    
    if (species == "human") {
      ensembl <- useEnsembl(biomart = "genes")
      ensembl <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)
      
      print("assign ensmbl to .GlobalEnv")
      ensembl <<- ensembl
    }
    
  }
 
  gene_type <- switch (gene_type,
    "ENSEMBL" = 'ensembl_gene_id',
    "SYMBOL" = 'external_gene_name'
  )
  
  gene_go_anno <- getBM(attributes=c("ensembl_gene_id", "external_gene_name","gene_biotype",
                                     "description","go_id","namespace_1003","name_1006","definition_1006"), 
                        filters= gene_type, value =gene, mart = ensembl)
  
  
  
  gene_go_anno_format2 <- 
  gene_go_anno %>% 
    filter(go_id != "") %>%                         # 去掉空字符串,没有GO annotate
    mutate(GO_annotation = paste0(namespace_1003,": ",name_1006," (",go_id,")") ) %>% 
    select(ensembl_gene_id,external_gene_name,GO_annotation) %>% 
    group_by(ensembl_gene_id) %>% 
    summarise(GO_annotation = str_c(GO_annotation, collapse = ";;"))
  
  
  
  if (format == "wider") {
    
    return( gene_go_anno_format2)
  }
  
  
  if (format == "longer") {
    return( gene_go_anno)
  }
  
  

}