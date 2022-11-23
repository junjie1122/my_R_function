


# @description: 对差异分析结果做差异火山图

# @parameters: deg_res 差异分析结果
# @parameters: x y x轴logFC y轴为pvalue or adj.pvalue
# @parameters: use_adjust_p T or F 是否使用adj.pvalue来筛选差异基因
# @parameters: color_colname 新增的差异基因的状态的列 默认Change
# @parameters: out_path 为NULL时不生成文件
# @parameters: point_color 点的颜色
# @parameters: show_top_gene_name T or F 是否显示 top的基因的标签 
# @parameters: top 显示top前几 根据pvalue or adj.pvalue 作为top的标准
# @parameters: distance 标签距离最大 or 最小的 logFC的距离
# @parameters: expand x轴拓展的大小


# return: p

# deg_res = limma_res
# x = "logFC"
# y =  "P.Value"
# color_colname = "Change"
# abs_logFC_cut = 0.2



get_volcano_plot <- function(deg_res,x, y, use_adjust_p,color_colname = "Change",abs_logFC_cut = 1,p_cut = 0.05,
                             out_path = NULL,
                             point_color = c("#2874C5", "grey", "#f87669"),
                             show_top_gene_name = T,top = 10,point_size = 1,alpha = 0.8,
                             distance = 1,
                             width = 8,height = 8,
                             prefix = "",
                             expand = c(0.3,0.3),
                             nudge_y = 0
                         ){
  
  

  if (!is.null(out_path)) {
    if (!dir.exists(paste0(out_path))) {
      dir.create(paste0(out_path),recursive = T)
    }
  }
  
  
  
  library(tidyverse)
  library(ggrepel)
  
  
  palette =  point_color
  
  
  if (use_adjust_p == F) {
    
    
    deg_res[,color_colname] <- ifelse(deg_res$P.Value>p_cut,'Not',                                        
                                   ifelse( deg_res$logFC>=abs_logFC_cut,'Up',                               
                                           ifelse( deg_res$logFC <= -abs_logFC_cut,'Down','Not') )  )    
    
  }else{
    
    deg_res[,color_colname] <- ifelse(deg_res$adj.P.Val>0.05,'Not',                           
                                     ifelse( deg_res$logFC>=abs_logFC_cut,'Up',                              
                                             ifelse( deg_res$logFC <= -abs_logFC_cut,'Down','Not') )  )
  }
  
  
  
  
  
  
  
  p <- ggplot(data = deg_res, aes(x = .data[[x]], y = -log10(.data[[y]]))) +
    geom_point(aes(color = .data[[color_colname]]),size = point_size,alpha = alpha) +
    
    geom_vline(xintercept = abs_logFC_cut,color = palette[3], linetype="longdash",size = 0.4) +
    geom_vline(xintercept = -abs_logFC_cut,color = palette[1], linetype="longdash",size = 0.4) +
    
    geom_hline(yintercept = -log10(p_cut) ,color = palette[2],linetype="longdash", size = 0.4) +
    
    
    scale_color_manual(
      values = palette,
      guide = guide_legend(override.aes = list(label = ""))
    )+
    theme_bw()+
    theme(
      
      axis.text.x = element_text(size = 15,      
                                 color = "black", 
                                 angle = 0), 
      axis.text.y = element_text(size = 15,  
                                 color = "black",
                                 angle = 0),
      axis.title.x = element_text(size = 18, 
                                  color = "black",
                                  angle = 0),
      axis.title.y = element_text(size = 18, 
                                  color = "black",
                                  angle = 90),
      legend.title = element_text(color="black",
                                  size=18),
      
      legend.text = element_text(color="black", 
                                 size = 15),
      
      plot.margin = unit(c(1,1,1,1),"lines"),
      legend.key.size = unit(0.3, "inches"),
      
      panel.grid.major = element_blank(),     # 主要网格
      panel.grid.minor = element_blank(),      # 次要网格
      
    )
  
  
  
  
  
  if (show_top_gene_name == T) {
    
    
    deg_res$gene_name <- rownames(deg_res)
    
    Up_label_data <- deg_res[deg_res[,color_colname] == "Up",] %>% head(top)
    
    Down_label_data <- deg_res[deg_res[,color_colname] == "Down",] %>% head(top)
    
    
    max_logFC <- deg_res$logFC %>% max()
    min_logFC <- deg_res$logFC %>% min()
    
    
    
    p1 <- p +
      scale_x_continuous(expand = expand) +
      scale_y_continuous(expand = c(0.1,0.1)) +
      geom_label_repel(
        
        data = Up_label_data,
        size = 3,
        aes(label = gene_name),
        nudge_y      = nudge_y ,
        direction    = "y",
        hjust        = 0,
        segment.color = "gray50",
        segment.size = 0.1,
        segment.linetype = 1,
        max.overlaps = 10,
        max.iter = 1000000,
        max.time = 10,
        nudge_x =  max_logFC + distance - Up_label_data[,x],
        min.segment.length = 0,
        # fontface = "bold"
      )+
      geom_label_repel(
        
        data = Down_label_data,
        size = 3,
        aes(label = gene_name),
        nudge_y      = nudge_y ,
        direction    = "y",
        hjust        = 0,
        segment.size = 0.1,
        segment.color = "gray50",
        segment.linetype = 1,
        max.overlaps = 10,
        max.iter = 1000000,
        max.time = 10,
        nudge_x =  min_logFC - distance - Down_label_data[,x],
        min.segment.length = 0,
        # fontface = "bold"
      ) 
    
  }
  
  
  if (!is.null(out_path)) {
    
    pdf(paste0(out_path,"/",prefix,"_volcano_plot",".pdf"),width = width,height = height)
    print(p1)
    dev.off()
    
    write.csv(deg_res,file = paste0(out_path,"/",prefix,"_volcano_plot_data",".csv"))
    
  }

  res <- list(plot = p,plot_data = deg_res)
  
  return(p1)
  
}
