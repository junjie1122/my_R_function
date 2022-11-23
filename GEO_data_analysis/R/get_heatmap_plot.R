


# @description: 差异基因热图

# @parameters: exp_data 表达矩阵
# @parameters: my_deg_gene 指定差异基因向量
# @parameters: deg_res limma差异分析的结果 ，my_deg_gene为NULL时使用deg_res的差异基因
# @parameters: top 是否选择pvalue最小的前几个基因 为NULL则全部
# @parameters: anno_col group_data
# @parameters: my_cols 颜色


# get_heatmap_plot(exp_data = exp_data,deg_res = limma_res,anno_col = group_data)



get_heatmap_plot <- function(exp_data,my_deg_gene = NULL,deg_res,top = NULL,anno_col,my_cols = c("green","red"),out_path,width = 6,height = 6){
  
  library(dplyr)
  library(pheatmap)
  library(scales)
  
  
  
  if (!is.null(out_path)) {
    if (!dir.exists(paste0(out_path))) {
      dir.create(paste0(out_path),recursive = T)
    }
  }
  
  
  
  
  if (is.null(my_deg_gene) | !is.null(deg_res)) {
    
    if (is.null(top)) {
      select_gene <- deg_res %>% filter(Change != "Not") %>% rownames()
      select_expr <- exp_data[select_gene,]
      
    }else{
      select_gene <- deg_res %>% filter(Change != "Not") %>% head(top) %>% rownames()
      select_expr <- exp_data[select_gene,]
      
    }

  }
  
  
  if (!is.null(my_deg_gene)) {
    elect_expr <- exp_data[my_deg_gene,]
  }
  
  
  
  
  
  anno_col <- group_data
  
  scale_data <- t(scale(t(select_expr)))
  
  
  value2 <- 1:20
  rr <- range(value2)
  svals <- (value2-rr[1])/diff(rr)
  f <- colorRamp(my_cols)
  colors <- rgb(f(svals)/255)
  
  
  
  pheatmap( scale_data ,fontsize = 10,border_color = NA,color =colors,show_colnames = F,
            annotation_col = anno_col,
            show_rownames = F,
            width = width,
            height = height,
            filename = paste0(out_path,"/heatmap_plot_DEG.pdf")
  
  
  )
  

  
}




