



# @description: PCA图

# exp_data
# coldata = group_data
# col_colname = "Group"
# get_pca_plot(exp_data = exp_data,coldata = group_data,col_colname = "Group",ellipse = T,out_path = "./output/")


get_pca_plot <- function(exp_data,coldata,col_colname,shape_colname = NULL,ellipse = F,
                     PC_X = "PC1",PC_Y = "PC2",point_size = 5,out_path,
                     width = 8,height = 8){
  
  
  library(dplyr)
  library(ggplot2)
  
  
  if (!is.null(out_path)) {
    if (!dir.exists(paste0(out_path))) {
      dir.create(paste0(out_path),recursive = T)
    }
  }
  
  
  
  
  pca <- prcomp(t(exp_data))
  percentVar <- round(100*(pca$sdev^2 / sum( pca$sdev^2 )),1)
  names(percentVar) <- paste0("PC",1:length(percentVar))
  
  
  d <- data.frame(PC1 = pca$x[,1],PC2 = pca$x[,2],PC3 = pca$x[,3],PC4 = pca$x[,4])
  d <- cbind(d,coldata)
  
  

  p <- ggplot(d, aes(.data[[PC_X]],.data[[PC_Y]])) +
    geom_point(aes(#shape = .data[[shape_colname]] , 
      colour = .data[[col_colname]]),size = point_size) +

    xlab(paste0(PC_X,", Variance: ", percentVar[PC_X], "%")) +
    ylab(paste0(PC_Y,", Variance: ", percentVar[PC_Y], "%")) +
    scale_color_manual(values = c("darkorange2", "dodgerblue4","darkgreen","darkyellow"))+
    theme_bw()+
    theme(
      
      axis.text.x = element_text(size = 15,  # 修改X轴上字体大小，
                                 color = "black", # 颜色
                                 #face = "bold", 
                                 angle = 0), #角度
      axis.text.y = element_text(size = 15,  
                                 color = "black",
                                 #face = "bold", 
                                 angle = 0),
      axis.title.x = element_text(size = 18, 
                                  color = "black",
                                  #face = "bold",
                                  angle = 0),
      axis.title.y = element_text(size = 18, 
                                  color = "black",
                                  #face = "bold", 
                                  angle = 90),
      legend.title = element_text(color="black", # 修改图例的标题
                                  size=18, 
                                  #face="bold"
      ),
      legend.text = element_text(color="black", # 设置图例标签文字
                                 size = 15,
                                 #face = "bold"
      ),
      plot.margin = unit(c(1,1,1,1),"lines"),
      legend.key.size = unit(0.3, "inches"),
      # panel.grid.major = element_blank(),     # 主要网格
      # panel.grid.minor = element_blank(),      # 次要网格
    )
  
  
  if (isTRUE(ellipse)) {
    p <- p+
      stat_ellipse(aes(color = .data[[col_colname]]),        ##添加置信区间
                   geom = "polygon",
                   alpha = 0,size = 0.8,
                   linetype = 2)
  }
  

  
  if (!is.null(out_path)) {
  pdf(file = paste0(out_path,"PCA_plot_",PC_X,"_",PC_Y,".pdf"),width = width,height = height)
  print(p)
  dev.off()
  }
  

  return(p)
  
  
}
