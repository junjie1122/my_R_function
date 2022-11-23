

# @description: 获得基因的biomaRt注释

# @parameters: species 指定物种 可选择 "human" or "mouse"
# @parameters: gene_type 基因类型 可选择 "ENSEMBL" or "SYMBOL"
# @parameters: version 指定ENSEMBL数据库版本  hg19 ensembl_version 为75  
# @parameters: host 指定其他版本的数据库 比如 hg19 "https://feb2014.archive.ensembl.org"
#                                               参考https://asia.ensembl.org/Help/ArchiveList  
#                                                   https://support.bioconductor.org/p/9140298/



# @parameters: add_attr 输入要增加的attributes的条目 向量形式
# @parameters: del_attr 输入要去除的attributes的条目 向量形式


# @return: 返回 gene_anno

# @example get_gene_go_anno(gene = gene)
#          gene_anno <- get_gene_go_anno(gene = rownames(exp_dat),species = "human",host = "https://feb2014.archive.ensembl.org",gene_type = "ENSEMBL")




get_gene_anno <- function(gene,species = "human",gene_type = "ENSEMBL",host = NULL,ensembl_version = NULL,add_attr = c(),del_attr = c()){
  
  library(biomaRt)
  library(stringr)
  library(dplyr)
  
  
  dataset <- switch (species,
    "human" = "hsapiens_gene_ensembl",
    "mouse" = "mmusculus_gene_ensembl"
  )
  
  
  if (!exists("ensembl")) {
    
      
      if (is.null(host)) {
        ensembl <- useEnsembl(biomart = "genes")
      }else{
        ensembl <- useEnsembl(biomart = "genes",host = host,version = ensembl_version)
      }
      
      
      ensembl <- useDataset(dataset = dataset, mart = ensembl)
      
      print("assign ensmbl to .GlobalEnv")
      ensembl <<- ensembl

  }
  
  
  gene_type <- switch (gene_type,
                       "ENSEMBL" = 'ensembl_gene_id',
                       "SYMBOL" = 'external_gene_name'
  )
  
  
  init_attr <- c('ensembl_gene_id',"external_gene_name","gene_biotype",
                 "chromosome_name","start_position","end_position", "description")
  
  
  attr <- c(init_attr,add_attr)
  
  attr <- attr[ !attr %in% del_attr]
  
  
  
  # hg19的的版本的ensembl数据库，attributes 和 filter 有些区别
  
  if (!is.null(ensembl_version)) {
    
    if (ensembl_version == 75) {
      
      attr <- attr %>% str_replace("external_gene_name","external_gene_id")
      
      gene_type <- switch (gene_type,
                           "ENSEMBL" = "ensembl_gene_id",
                           "SYMBOL" = "external_gene_id"
      )
      
    }
    
  }
  

  
  
  gene_anno <- getBM(attributes= attr, 
                     filters= gene_type, 
                     value =gene, 
                     mart = ensembl
  )
  
  # if (ensembl_version == 75) {
  # gene_anno <- gene_anno %>% rename("external_gene_name" = "external_gene_id")
  # }
  
  
  
  return(gene_anno)

}
