

# @description: 使用limma做差异分析

# @parameters: data 基因表达矩阵
# @parameters: data_type "RNA-seq" 或者 "array" 
# @parameters: coldata 样本信息 
# @parameters: ref_group 指定控制组的名称
# @parameters: my_formula 公式 比如 "~ gender + Group" 要研究的变量放最后面
# @parameters: use_filterByExpr True-使用filterByExpr函数过低表达, Flase不过滤(预先过滤好的情况下) 
# @parameters: use_adjust_p abs_logFC_cut p_cut 筛选差异基因的方式和阈值
# @parameters: out_path 指定输出路径 如果为NULL 则不输出差异分析结果的csv文件和limma_voom_Mean-variance_trand.pdf


# return: limma_DEG_result 


# 参考: https://www.jianshu.com/p/4a5508a83b3d
#       http://www.360doc.com/content/21/0714/12/76149697_986501444.shtml


 

run_limma <- function(data,data_type,coldata,ref_group,my_formula,
                      use_filterByExpr = T,
                      out_path = NULL,
                      use_adjust_p = F,abs_logFC_cut = 1,p_cut = 0.05){
  
  library(limma)
  library(dplyr)
  library(edgeR)
  
  if (!is.null(out_path)) {
    if (!dir.exists(paste0(out_path,"/run_limma_out/"))) {
      dir.create(paste0(out_path,"/run_limma_out/"),recursive = T)
    }
  }
  
  
  all_var <- my_formula %>% str_split(" ") %>% unlist() %>% .[!. %in% c("~","+")]
  
  
  
  if (length(all_var) > 1) {
    for (var in all_var) {
      assign(var,value = coldata[,var,drop = T] )
      
    }
  }

  
  coldata <- as.data.frame(coldata)
  # formal 要研究的变量最后面 提取最后面的变量
  study_col <- my_formula %>% str_remove_all(".* ")
  
  print(coldata)
  print(study_col)
  # 对研究的变量重新指定参考因子 控制组为参考因子
  assign(study_col,value = coldata[,study_col,drop = T] %>% unlist() %>% as.factor() %>% relevel(ref = ref_group) )
  
  
  
  
  design <- model.matrix(as.formula(my_formula))  
  rownames(design)=colnames(data)
  head(design)
  
  
  if (data_type == "array") {
    
    fit <- lmFit(data,design)                                             #为每个基因拟合线性模型
    fit2 <- eBayes(fit)
    colnames(coef(fit2))
    limmaDEG<- topTable(fit2,coef = ncol(design),n = Inf)                 #从线性模型拟合中提取排名靠前的基因的表格   
    limmaDEG <- limmaDEG[order(limmaDEG$adj.P.Val,decreasing = F),]
  }
  
  else if(data_type == "RNA-seq"){
    
    
    dge <- DGEList(counts=data)
    
    # 过滤低表达
    if (use_filterByExpr) {
      keep.exprs <- filterByExpr(dge,design = design)
      dge <- dge[keep.exprs,, keep.lib.sizes=FALSE]
    }

    
    # 计算文库因子
    dge <- calcNormFactors(dge,method = "TMM")
    
    
    if (!is.null(out_path)) {
      
      pdf(file = paste0(out_path,"/run_limma_out/","limma_voom_Mean-variance_trand.pdf"),width = 6,height = 5)
      v <- voom(dge,design,plot=TRUE)
      dev.off()
    }else{
      v <- voom(dge,design,plot=TRUE)
    }
    
    
    fit <- lmFit(v, design)
    fit <- eBayes(fit,)
    
    limmaDEG <- topTable(fit, coef=ncol(design),n = Inf)
    limmaDEG <- limmaDEG[order(limmaDEG$adj.P.Val,decreasing = F),]
    
  }
  else{
    stop("must be array or RNA-seq")
    
  }

  
  
  ## 筛选差异基因
  
  if (use_adjust_p == F) {
    
    limmaDEG$Change=ifelse(limmaDEG$P.Value>p_cut,'Not',                                        
                     ifelse( limmaDEG$logFC>=abs_logFC_cut,'Up',                               
                             ifelse( limmaDEG$logFC <= -abs_logFC_cut,'Down','Not') )  )    
    
  }else{
    
    limmaDEG$Change=ifelse(limmaDEG$adj.P.Val>0.05,'Not',                           
                     ifelse( limmaDEG$logFC>=abs_logFC_cut,'Up',                              
                             ifelse( limmaDEG$logFC <= -abs_logFC_cut,'Down','Not') )  )
    
  }
  
  

  print(table(limmaDEG$Change))
  
  limmaDEG_res <- limmaDEG_res
  
  if (!is.null(out_path)) {
    write.csv(limmaDEG_res,file = paste0(out_path,"/run_limma_out/","limma_DEG_result.csv"),row.names = T)
  }
 
  
  if (!is.null(out_path)) {
   save(limmaDEG_res,file = paste0(out_path,"/RData/","limmaDEG_res.RData") )
  }
  
  return(limmaDEG_res)
}











