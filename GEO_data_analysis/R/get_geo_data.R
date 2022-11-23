

# @description: 下载GEO数据，用于下游分析

# @parameters: gse_id gse号
# @parameters: type array or rna-seq  rna-seq则下载supplement files里面的表达矩阵

# @parameters: title_matched_pattern 获取样本信息的title列的内容的匹配模式 为NULL不生成分组信息
#              比如 "PBMCs_Control|PBMCs_IS"
# @parameters: GPL_ids_cols id的列名 第一个为探针id的列名 第二三个为可能的GENE SYMBOL的列名
# @parameters: out_path 指定输出路径 如果为NULL 则不输出csv文件
# @parameters: local_series_matrix 为NULL使用在线下载数据，否则提供本地的series_matrix文件


# @return:
          # exp_dat = out_put_exp_dat,           # 表达矩阵
          # ids = ids,                           # ids
          # meta = meta,                         # 样本信息
          # GPL_anno_data = GPL_anno_data,       # 基因注释信息
          # group_data = group_data              # 分组信息


# @example:
          # gse_id <- "GSE22255"
          # title_matched_pattern <- "PBMCs_Control|PBMCs_IS"
          # out_path = "./output/

get_geo_data <- function(gse_id,type = "array",local_series_matrix = NULL,
                         title_matched_pattern = NULL,
                         GPL_ids_cols = c("ID","Gene symbol","Gene Symbol"),
                         out_path = "./output/"
                         
                         ){
  
  library(GEOquery)
  library(tidyverse)
  library(stringr)
  library(limma)
  options('download.file.method.GEOquery' = 'libcurl' )
  
  
  
  
  
  if (!is.null(out_path)) {
    if (!dir.exists(paste0(out_path))) {
      dir.create(paste0(out_path,"/get_geo_data_out/"),recursive = T)
      dir.create(paste0(out_path,"/RData"))
    }
  }
  
  
  
  
  if (!is.null(local_series_matrix)) {
    
    gset <- getGEO(filename = local_series_matrix )
    
    dat <- exprs(gset) %>% as.data.frame()
    meta <- gset@phenoData@data
    colnames(meta)=gsub("[:_]ch1", "", colnames(meta))
    
    GPL_anno_data <- gset@featureData@data
    
    
    # 如果提供title列的获取分组的正则匹配模式 则生成分组信息 group_data
    if (!is.null(title_matched_pattern)) {
      
      group_data <- meta %>% select(title) %>% 
        mutate(Group = str_extract_all(title,title_matched_pattern,simplify = T)[,1]) %>% select(-title)
    }else{
      group_data <- meta %>% select(title) %>% mutate(across(.fns = ~ gsub("[1-9]", "",.x)) )
    }
    
  }
  
  else{
    

    if (type == "rna-seq") {
      
      gset <- getGEO(gse_id, GSEMatrix =TRUE, AnnotGPL=F,getGPL = F)
      a = gset[[1]]
      
      #获取样本全部信息
      pd = pData(a)
      
      
      # 选择比较有用的信息，重命名列名
      meta = pd   # %>% dplyr::select(title, dplyr::ends_with("ch1"))
      colnames(meta)=gsub("[:_]ch1", "", colnames(meta))
      head(meta[,1:3])
      
      # 如果提供title列的获取分组的正则匹配模式 则生成分组信息 group_data
      if (!is.null(title_matched_pattern)) {
        
        group_data <- pd %>% select(title) %>% 
          mutate(Group = str_extract_all(title,title_matched_pattern,simplify = T)[,1]) %>% select(-title)
      }else{
        group_data <- pd %>% select(title) %>% mutate(across(.fns = ~ gsub("[1-9]", "",.x)) )
      }
      
      
      # 获取 supplementary_file 链接 并下载
      if (!is.null(out_path)) {
        download.file(url = gset[[1]]@experimentData@other$supplementary_file,
                      destfile = paste0(out_path,"/",gse_id,"_exp_data.txt.gz"),method = "auto")
        

        # 解压得到.txt的表达矩阵
        R.utils::gunzip(paste0(out_path,"/",gse_id,"_exp_data.txt.gz"),remove=F)
        
        
        exp_dat <- read_delim(paste0(out_path,"/",gse_id,"_exp_data.txt"), 
                                                      delim = "\t", escape_double = FALSE, 
                                                      trim_ws = TRUE)
        
      }
      

      exp_dat = exp_dat                    
      meta = meta                          
      group_data = group_data   
      GPL_anno_data = NULL   
      ids = NULL 
      
      
      if (!is.null(out_path)) {
        
        write.csv(exp_dat,file = paste0(out_path,"/get_geo_data_out/","exp_data.csv"),row.names = T)
        write.csv(meta,file = paste0(out_path,"/get_geo_data_out/","meta_data.csv"),row.names = T)
        write.csv(group_data,file = paste0(out_path,"/get_geo_data_out/","_group_data.csv"),row.names = T)
      }
      
      

    }
    
    
    
    if (type == "array") {
      
      
      # 下载
      gset <- getGEO(gse_id, GSEMatrix =TRUE, AnnotGPL=T,getGPL = T)
      class(gset)
      
      a = gset[[1]]                        # 第一个元素里面放的表达矩阵 
      dat = exprs(a) %>% as.data.frame()   # 将对象a 转化为表达矩阵， 基因ID是探针名
      head(dat)
      
      #获取样本全部信息
      pd = pData(a)
      
      
      # 选择比较有用的信息，重命名列名
      meta = pd   # %>% dplyr::select(title, dplyr::ends_with("ch1"))
      colnames(meta)=gsub("[:_]ch1", "", colnames(meta))
      head(meta[,1:3])
      
      
      # 如果提供title列的获取分组的正则匹配模式 则生成分组信息 group_data
      if (!is.null(title_matched_pattern)) {
        
        group_data <- pd %>% select(title) %>% 
          mutate(Group = str_extract_all(title,title_matched_pattern,simplify = T)[,1]) %>% select(-title)
      }else{
        group_data <- pd %>% select(title) %>% mutate(across(.fns = ~ gsub("[1-9]", "",.x)) )
      }
      
      
      GPL_name <- gset[[1]]@annotation
      
      # GPL列名的描述信息 
      gset[[1]]@featureData@varMetadata %>% head()
      
      
      GPL_anno_data <- gset[[1]]@featureData@data
      dim(GPL_anno_data)


    
      #只要ID列和Gene symbol列，并过滤 Gene symbol 为 "",NA,"---" 的行
      
      match_id_cols <- colnames(GPL_anno_data)[colnames(GPL_anno_data) %in% GPL_ids_cols]
      
      ids <- GPL_anno_data[,match_id_cols] %>% `colnames<-`(c("prob_id","gene_symbol")) %>% 
        filter(!gene_symbol %in% c("","---")) %>% 
        filter(!is.na(gene_symbol))
      
      head(ids)
      dim(ids)
      
  
      # 由于过滤了探针id 为了保证 dat 和 ids 的探针 都是共有的 进行匹配
      # 以防止后续 left_join 时出现 dat有 ids没有的探针 出现NA
      
      ids <-  ids[ids$prob_id %in% rownames(dat),]   #选取和dat行名 相匹配的行
      dim(ids)
      
      dat <-  dat[ids$prob_id,]                      #选取表达矩阵dat中的rowname与探针ID一样的行
      dim(dat)
    
    
      if (F) {
        
        # gene_symbol存在一个探针对应多个基因的情况，取第一个基因  (废弃)
        
        ids$gene_symbol <- trimws(str_split(ids$gene_symbol,'///',simplify = T)[,1])
        
        #id转换
        ids <-  ids[ids$prob_id %in% rownames(dat),]   #选取和dat行名 相匹配的行
        dim(ids)
        
        dat <-  dat[ids$prob_id,]                      #选取表达矩阵dat中的rowname与探针ID一样的行
        dim(dat)
        head(dat)
        
        
        dat <- cbind(dat,ids)  %>% dplyr::select(-prob_id) %>% dplyr::select(gene_symbol,everything())
        head(dat)
        
      }
      
      # 一个探针对应多个基因 保留所有对应情况
      ids <- ids %>% separate_rows(gene_symbol,sep = "\\/\\/\\/")
      dim(ids)
      
      
      # 连接 dat 和 ids
      tmp_dat <- dat
      tmp_dat$prob_id <- rownames(tmp_dat)
      
      tmp_dat <- left_join(tmp_dat,ids) %>% select( -c(prob_id) )
      
      gene_symbol <- tmp_dat$gene_symbol
      
      tmp_dat <- tmp_dat %>% select(-gene_symbol) %>% as.matrix() %>% modify(as.numeric)
      rownames(tmp_dat) <- gene_symbol
      dim(tmp_dat)
      
      
      #limma::avereps 函数 将重复gene_symbol 去重取均值
      out_put_exp_dat <- avereps(tmp_dat) %>% as.data.frame()
      head(out_put_exp_dat)
      dim(out_put_exp_dat)
    
    
    
      if (!is.null(out_path)) {
        
        write.csv(out_put_exp_dat,file = paste0(out_path,"/get_geo_data_out/","exp_data.csv"),row.names = T)
        write.csv(meta,file = paste0(out_path,"/get_geo_data_out/","meta_data.csv"),row.names = T)
        write.csv(group_data,file = paste0(out_path,"/get_geo_data_out/","_group_data.csv"),row.names = T)
        write.csv(ids,file = paste0(out_path,"/","ids.csv"),row.names = F)
        write.csv(GPL_anno_data,file = paste0(out_path,"/get_geo_data_out/","GPL_anno_data.csv"),row.names = T)
        
      }
    
      
      
    }
    

  
  }
  
  
  
  get_geo_out <- list(exp_dat = out_put_exp_dat,   # 表达矩阵
                      ids = ids,                           # ids
                      meta = meta,                         # 样本信息
                      GPL_anno_data = GPL_anno_data,       # 基因注释信息
                      group_data = group_data              # 分组信息
                      )
  
  
  
  
  if (!is.null(out_path)) {
    save(get_geo_out,file = paste0(out_path,"/RData/get_geo_out.RData"))
  }
  
  
  invisible(get_geo_out)
  
}






# 补充信息
##(3)下载Supplementary file
# gset$GSE22255_series_matrix.txt.gz@experimentData@other$supplementary_file

# gset[[1]]@experimentData@other$supplementary_file

if (F) {
  
  test_data <- data.frame(a = "id1",b = 100,c = 50)
  test_ids <- data.frame(a = c("id1","id1","id1"),sym = c("KK","jj","pp"))
  left_join( test_data,test_ids)
  
}


















