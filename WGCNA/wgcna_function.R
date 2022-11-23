library(openxlsx)
library(stringr)
library(DESeq2)
# library(edgeR)
library(patchwork)
library(WGCNA)
library(ggpubr)
library(biomaRt)
library(tidyverse)
library(dplyr)
library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("mutate", "dplyr")


old_max_print <- getOption("max.print")
options(max.print=50)
enableWGCNAThreads(nThreads = round(0.3*parallel::detectCores()) ) 
cor <- WGCNA::cor


my_theme <- function(){
  
  theme(
    
    axis.text.x = element_text(size = 15,
                               color = "black",
                               angle = 0), 
    axis.text.y = element_text(size = 15,  
                               color = "black",
                               angle = 0),
    
    axis.title.x = element_text(size = 18,
                                color = "black",
                                face = "bold",
                                angle = 0),
    
    axis.title.y = element_text(size = 18, 
                                color = "black",
                                face = "bold", 
                                angle = 90),
    
    legend.title = element_text(color="black",
                                size=18, 
                                face="bold"),
    
    legend.text = element_text(color="black", 
                               size = 15),
    plot.margin = unit(c(1,1,1,1),"lines"),
    legend.key.size = unit(0.3, "inches"),
    strip.text = element_text(size = 15)
  )

}




# 数据标准化 数据转换
data_normalization <- function(data,method = NULL){
  
  if (is.null(method) ) {
    data <- data
  }else if (method == "vst") {
    
    data <- varianceStabilizingTransformation(as.matrix(data) ,blind = T)
    
  }else if (method == "cpm") {
    
    data <- edgeR::cpm(data)
    
  }else if (method == "logcpm") {
    
    data <- log10(edgeR::cpm(data)+1)
    
  }else if (method == "log") {

    data <- log10(data+1)  
    
  }else{
    stop(paste0("not method :",method))
  }
  
  data <- as.data.frame(data)
  print("assign:data")
  assign("data",value = as.data.frame(data),envir = .GlobalEnv) 
  return(data)
  
}



# 过滤基因  并做检查异常基因和样本
gene_filter <- function(data,filter_gene_method,keep_gene_num = 5000,filter_threshold = NULL,out_path){
  
  dir.create(paste0(out_path,"/plot"),recursive = T)
  dir.create(paste0(out_path,"/RData"),recursive = T)
  

  if (is.null(filter_threshold)){

    # 按 指定的基因数目过滤
    if (keep_gene_num > 1) {
      if (filter_gene_method == "mad") {
        data <- data %>% mutate(mad = pmap_dbl(data,~ mad(c(...)))) %>% 
          arrange(desc(mad)) %>% .[1:keep_gene_num,] %>% 
          dplyr::select(-(mad))
      }else if (filter_gene_method == "var") {
        data <- data %>% mutate(var = pmap_dbl(data,~ var(c(...)))) %>% 
          arrange(desc(var)) %>% .[1:keep_gene_num,] %>% 
          dplyr::select(-(var))
      }
    }
    
    # 按 百分比 过滤
    if (keep_gene_num <=  1) {
      if (filter_gene_method == "mad") {
        data <- data %>% mutate(mad = pmap_dbl(data,~ mad(c(...)))) %>% 
          slice_max(mad, n = round(nrow(data)*keep_gene_num),with_ties = F) %>%
          dplyr::select(-(mad))
      }else if (filter_gene_method == "var") {
        data <- data %>% mutate(var = pmap_dbl(data,~ var(c(...)))) %>% 
          slice_max(var, n = round(nrow(data)*keep_gene_num),with_ties = F) %>% 
          dplyr::select(-(var))
      }
    }

  }
  # 按照指定阈值过滤
  else {

    data <- data %>% mutate(mad = pmap_dbl(data,~ mad(c(...)))) %>%
      filter(mad > filter_threshold) %>% dplyr::select(-(mad)) 
  }
  
  

  write.csv(data,file = paste0(out_path,"/WGCNA_input_data.csv"))
  data <- as.data.frame(t(data))
  dim(data)
  

  #检查异常基因和样本，
  gsg = goodSamplesGenes(data,verbose = 3)
  if (!gsg$allOK){
    
    if (sum(!gsg$goodGenes)>0)
      printFlush(paste("Removing genes:",
                       paste(names(data)[!gsg$goodGenes], collapse = ",")));
    if (sum(!gsg$goodSamples)>0)
      printFlush(paste("Removing samples:",
                       paste(rownames(data)[!gsg$goodSamples], collapse = ",")));
    # Remove the offending genes and samples from the data:
    data = data[gsg$goodSamples, gsg$goodGenes]
    
  }


  
  n_vlincRNA <- colnames(data) %>% str_detect("chr.*\\|") %>% sum()
  n_no_vlincRNA <- colnames(data) %>% str_detect("ENSG.*") %>% sum()

  

  print("assign: data, nGenes, nSamples,n_vlincRNA,n_no_vlincRNA")
  assign("data",value = as.data.frame(data),envir = .GlobalEnv)
  assign("nGenes",value = ncol(data),envir = .GlobalEnv)
  assign("nSamples",value = nrow(data),envir = .GlobalEnv)


  assign("n_vlincRNA",value = n_vlincRNA,envir = .GlobalEnv)
  assign("n_no_vlincRNA",value = n_no_vlincRNA,envir = .GlobalEnv)


  save(data,file = paste0(out_path,"/RData/data.RData"))
  save(nGenes,file = paste0(out_path,"/RData/nGenes.RData"))
  save(nSamples,file = paste0(out_path,"/RData/nSamples.RData"))
  
  return(data)

}



# 检查样本聚类
sample_clust <- function(data,trait,abHeight = "",width = 9,height = 6,out_path){
  

  
  
  sampleTree <- hclust(dist(data),method = "average") 
  traitColors = numbers2colors(trait, signed = FALSE)
  colnames(traitColors) <- colnames(trait)
  
  pdf(paste0(out_path,"/plot/Sample_dendrogram_and_trait_heatmap.pdf"),width = width,height = height)
  plotDendroAndColors(sampleTree, traitColors,abHeight = abHeight,
                      groupLabels = colnames(trait), 
                      main = "Sample dendrogram and trait heatmap") 
  
  dev.off()
  
  save(trait,file = paste0(out_path,"/RData/trait.RData"))
  
}



# remove outline sample
remove_outline_sample <- function(data,trait,cutHeight = 100000,width,height){
  
  sampleTree <- hclust(dist(data),method = "average") 
  
  clust = cutreeStatic(sampleTree, cutHeight = cutHeight,minSize = 10)
  print(table(clust))
  
  
  #保留没有离群的样本
  keepSamples = (clust==1)
  data = data[keepSamples,]
  trait <- trait[keepSamples,,drop = F]
  
  nGenes = ncol(data)
  nSamples = nrow(data)
  
  
  sampleTree <- hclust(dist(data),method = "average") 
  traitColors = numbers2colors(trait, signed = FALSE); 
  
  pdf(paste0(out_path,"/plot/Sample_dendrogram_and_trait_heatmap_remove_outline_sample.pdf"),width = width,height = height)
  plotDendroAndColors(sampleTree, traitColors, 
                      groupLabels = colnames(trait), 
                      main = "Sample dendrogram and trait heatmap") 
  dev.off()
  
  
  assign("data",value = data,envir = .GlobalEnv)
  assign("trait",value = trait,envir = .GlobalEnv)
  assign("nGenes",value = ncol(data),envir = .GlobalEnv)
  assign("nSamples",value = nrow(data),envir = .GlobalEnv)
  
  
  save(data,file = paste0(out_path,"/RData/data.RData"))
  save(trait,file = paste0(out_path,"/RData/trait.RData"))
  save(nGenes,file = paste0(out_path,"/RData/nGenes.RData"))
  save(nSamples,file = paste0(out_path,"/RData/nSamples.RData"))
  
  # return(data)
  
}




# sample_check
sample_check <- function(data,abHeight){

  data <- data_normalization(data = data,method = method)
  dim(data)


  data <- gene_filter(data = data,filter_gene_method = filter_gene_method,keep_gene_num = keep_gene_num,filter_threshold = filter_threshold)
  dim(data)


  sample_clust(data = data,trait = trait,width = 9,height = 6,abHeight = abHeight,out_path = out_path)

}





# 软阈值选择
power_choose <- function(data,cor_type,net_type,max_blockSize,width = 9,height = 6,out_path){
  
  powers = c(c(1:10), c(seq(12,30,2)))

  if(cor_type == "bicor"){

    sft <- pickSoftThreshold(blockSize = max_blockSize,data, powerVector = powers, verbose = 3,networkType =net_type,corFnc="bicor",corOptions = list(maxPOutliers= 0.05, robustY=FALSE)) 

  }else{

    sft <- pickSoftThreshold(blockSize = max_blockSize,data, powerVector = powers, verbose = 3,networkType =net_type) 

  }


  write.csv(sft,file = paste0(out_path,"/sft.csv"))
  
  # sft$powerEstimate 给出合适的power
  power <- sft$powerEstimate 
  
  
  # 第一张图 power和R^2关系图，第二张图 power和平均连接度关系
  pdf(paste0(out_path,"/plot/power_R2_Connectivity_relation.pdf"),width = 9,height = 6)
  par(mfrow = c(1,2));
  cex1 = 0.8;
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"))
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")
  
  abline(h=0.85,col="red")
  
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  dev.off()
  
  
  # 如果 sft$powerEstimate  为NA
  # 也就是没有合适的power使得R2大于0.85则使用下面的power
  type = net_type
  if (is.na(power)){
    power = ifelse(nSamples<20, ifelse(type == "unsigned", 9, 18),
                   ifelse(nSamples<30, ifelse(type == "unsigned", 8, 16),
                          ifelse(nSamples<40, ifelse(type == "unsigned", 7, 14),
                                 ifelse(type == "unsigned", 6, 12))
                   )
    )
  }
  
  
  print("assign: sft, power")
  print(power)
  assign("sft",value = sft,envir = .GlobalEnv)
  assign("power",value = power,envir = .GlobalEnv)

  save(power,file = paste0(out_path,"/RData/power.RData"))
  return(power)

}





# 网络构建
net_build <- function(data,power,max_blockSize,cor_type,net_type,TOM_type,mergeCutHeight,out_path){
  
  cor <- WGCNA::cor

  if (cor_type == "bicor") {

    net = blockwiseModules(data, power = power, maxBlockSize = max_blockSize,        
                         networkType = net_type,
                         TOMType = TOM_type,
                         corType = "bicor",
                         maxPOutliers = 0.05 ,robustY = FALSE,                              
                         minModuleSize =30,          
                         reassignThreshold = 0, mergeCutHeight = mergeCutHeight,   #mergeCutHeight是模块融合的阈值
                         numericLabels = TRUE, pamRespectsDendro = FALSE,
                         saveTOMs = TRUE,
                         saveTOMFileBase = paste0(out_path,"/RData/TOM"),
                         verbose = 3)
     
  }else{

    net = blockwiseModules(data, power = power, maxBlockSize = max_blockSize,        
                         networkType = net_type,
                         TOMType = TOM_type,                            
                         minModuleSize =30,          
                         reassignThreshold = 0, mergeCutHeight = mergeCutHeight,   #mergeCutHeight是模块融合的阈值
                         numericLabels = TRUE, pamRespectsDendro = FALSE,
                         saveTOMs = TRUE,
                         saveTOMFileBase = paste0(out_path,"/RData/TOM"),
                         verbose = 3)

  }

  
  # print("模块数量和基因数量:")
  print(table(net$colors))
  
  moduleColors <- labels2colors(net$colors)
  
  # 统计模块基因数量并保存
  table(moduleColors) %>% as_tibble() %>% 
    arrange(desc(n)) %>% 
    set_names(c("module","gene_count") ) %>% 
    write.csv(file = paste0(out_path,"/module_gene_count.csv"),row.names = F)
  
  #计算模块MES
  MEs0 = moduleEigengenes(data, moduleColors)$eigengenes          
  MEs = orderMEs(MEs0)                             
  head(MEs)
  
  
  print("assign: net, MEs, moduleColors")
  assign("net",value = net,envir = .GlobalEnv)
  assign("MEs",value = MEs,envir = .GlobalEnv)
  assign("moduleColors",value = moduleColors,envir = .GlobalEnv)
  
  
  net_build_out <- list(
    net = net,
    MEs = MEs,
    moduleColors = moduleColors
  )
  

  save(net,MEs,moduleColors,file = paste0(out_path,"/RData/net_build_out.RData"))
  # save(net_build_out,file = paste0(out_path,"/RData/net_build_out.RData"))
  return(net_build_out)

}


# 基因树状图
dendrograms_plot <- function(net,moduleColors,width = 9,height =6,out_path){
  
  pdf(file = paste0(out_path,"/plot/plotDendroAndColors.pdf"),width =width,height = height)
  plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],"Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  dev.off()
  
}



# 模块特征基因热图
mes_heatmap_plot <- function(data,MEs,moduleColors,width =9,height = 6,out_path){ 
  
  dir.create(paste0(out_path,"/plot/MEs_heatmap_plot"),recursive = T)
  
  for (col in colnames(MEs)) {
    
    ME=MEs[,col]
    ME_color <- str_remove_all(col,pattern = "ME")
    
    ## heatmap data
    h <- t(scale(data[,moduleColors == ME_color]))
    head(h)
    x = reshape2::melt(data = h,id.var = rownames(h))
    
    ## barplot data
    b <- data.frame(Sample = factor(rownames(data),levels = rownames(data)) ,
                    Eigengene = ME)
    
    p1 <- ggplot(x, aes(Var2, Var1)) +
      geom_tile(aes(fill = value),colour = alpha("white",0)) +
      scale_fill_gradient2(low = "green",mid = "black",high = "red")+
      theme_classic()+
      theme(legend.position = "none",
            axis.title = element_blank(),
            axis.text.y = element_blank(),
            axis.text.x = element_text(angle = 90,vjust = 0.5,size = 15),
            axis.line.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line.x = element_blank())
    
    ## barplot
    p2 <- ggplot(b,mapping  = aes(x = Sample,y = Eigengene)) + geom_bar(stat='identity',fill = ME_color)+
      theme_classic()+
      theme(axis.title = element_blank(),
            axis.line = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank())
    p3 <-  p1/p2
    
    
    pdf(paste0(out_path,"/plot/MEs_heatmap_plot/",col,".pdf"),width = width,height = height)
    print(p3)
    dev.off()
    
  }
}




# 模块特征基因关联trait
mes_trait_relation <- function(data,MEs,trait,width = 6,height = 9,mar = c(8,8,4,3),out_path){
  
  #计算模块MEs和trait相关性 进而选出出相关的模块
  moduleTraitCor = cor(MEs,trait, use = "p")           
  nSamples = nrow(data)
  nGenes = ncol(data)
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)   ##计算p值
  
  
  moduleTraitCor_Pvalue <- cbind(moduleTraitCor,moduleTraitPvalue)
  colnames(moduleTraitCor_Pvalue) <- c(paste0(colnames(moduleTraitCor),"_Cor"),paste0(colnames(moduleTraitCor),"_Pvalue"))
  head(moduleTraitCor_Pvalue)
  write.csv(moduleTraitCor_Pvalue,file = paste0(out_path,"/moduleTraitCor_Pvalue.csv"),row.names = T)
  
  
  # 根据模块数量设置高度
  height = 5+nrow(moduleTraitCor)*0.22
  
  
  textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPvalue, 1), ")", sep = "")
  dim(textMatrix) = dim(moduleTraitCor)
  
  pdf(file = paste0(out_path,"/plot/Module-trait_relationships.pdf"),width = width,height = height)
  par(mar = mar)   ## 设置图区变距大小 依次 下左上右
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = colnames(trait),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = FALSE,yLabelsPosition = "left",
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.8,
                 #textAdj = c(0.5, 0.5),
                 font.lab.x = 0.8,
                 font.lab.y = 0.8,
                 zlim = c(-1,1),
                 main = paste("Module-trait relationships"))
  dev.off()
  
  print("assign: moduleTraitCor, moduleTraitPvalue")
  assign("moduleTraitCor",value = moduleTraitCor,envir = .GlobalEnv)
  assign("moduleTraitPvalue",value = moduleTraitPvalue,envir = .GlobalEnv)
  
} 



# 筛选出和trait显著相关的模块
get_significent_mod <- function(moduleTraitPvalue){
  
  total_mod <- rownames(moduleTraitPvalue)
  significent_mod <- total_mod [rowSums(moduleTraitPvalue < 0.05) >= 1]
  
  #grey模块的基因是聚类是那些被剔除的基因的集合，无意义 不然后面会报错去掉
  significent_mod <- significent_mod[significent_mod != "MEgrey"]
  significent_mod<- substr(significent_mod,3,20)
  
  print("assign: significent_mod")
  assign("significent_mod",value = significent_mod,envir = .GlobalEnv)
  

  save(significent_mod,file = paste0(out_path,"/RData/significent_mod.RData"))
  return(significent_mod)
}


# 保存和trait显著相关的模块的基因
# merge_by 表达矩阵 和基因注释 合并的列名
# group_by 基因种类的名
save_significent_mod_gene <- function(data,MEs,moduleColors,significent_mod,gene_anno = NULL,merge_by = "ensembl_gene_id",group_by ="gene_biotype"  ,out_path){
  
  dir.create(paste0(out_path,"/significent_module_genes/txt"),recursive = T)
  dir.create(paste0(out_path,"/significent_module_genes/csv"),recursive = T)
  dir.create(paste0(out_path,"/significent_module_genes/summary"),recursive = T)
  # modNames = labels2colors(names(MEs))
  
  mod_color <- str_remove_all(names(MEs),pattern = "ME")
  t_data <- t(data) %>% as.data.frame()
  list_module_gene <- list()
  
  for (i in 1:length(significent_mod)) {
    
    module = significent_mod[i]                                                       
    moduleGenes = moduleColors == module
    Genes <- colnames(data)[moduleGenes]
    
    list_module_gene[[module]] <- Genes
    
    module_gene_exp_data <- t_data[Genes,]

    if (!is.null(gene_anno)) {

      module_gene_exp_data[,merge_by] <- rownames(module_gene_exp_data)
      module_gene_exp_data <- as.data.frame(module_gene_exp_data)
      module_gene_exp_data <- left_join(module_gene_exp_data,gene_annotate_with_description,by = merge_by)
      head(module_gene_exp_data)
      rownames(module_gene_exp_data) <- module_gene_exp_data[,merge_by,drop = T]
      
      # 模块基因类型分类
      module_gene_exp_data_summary <- module_gene_exp_data  %>% group_by(.data[[gene_biotype]])  %>% summarise(gene_count = n())
      total_gene <- sum(module_gene_exp_data_summary$gene_count)
      module_gene_exp_data_summary <- module_gene_exp_data_summary %>% 
          mutate(percent = gene_count/total_gene,total_gene = total_gene) %>% 
          arrange(desc(gene_count))

      write.csv(module_gene_exp_data_summary,file = paste0(out_path, "/significent_module_genes/summary/",module,"_genes_type_count.csv"), row.names = F )

    }


    write.table(x = Genes,file = paste0(out_path,"/significent_module_genes/txt/",module,"_genes.txt"),
                quote = F,row.names = F,col.names = F)
    
    write.csv(module_gene_exp_data,file = paste0(out_path,"/significent_module_genes/csv/",module,"_genes.csv"),
              row.names = T )
    
  }
  
  print("assign: list_module_gene")

  assign("list_module_gene",value = list_module_gene,envir = .GlobalEnv)
  save(list_module_gene,file = paste0(out_path,"/RData/list_module_gene.RData"))

  return(list_module_gene)
  
}




# 计算基因和模块特征基因MEs相关性 即 MM
get_MM <- function(data,MEs,moduleColors,out_path){
  
  nSamples <- nrow(data)
  tmp_MEs <- MEs %>% set_names(str_remove_all(colnames(.),"ME"))
  
  geneModuleMembership = as.data.frame(cor(data, tmp_MEs, use = "p")) %>% 
    rename_with(~paste0("MM.",.))
  
  geneModuleMembership_Pvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples)) %>% 
    rename_with(~paste0("p.",.))
  
  
  MM_table <- cbind(geneModuleMembership,geneModuleMembership_Pvalue) %>% mutate(module = moduleColors,.before = 1)
  head(MM_table)
  
  # out_path为 NULL 则不执行{}内代码 
  if (!is.null(out_path)) {
    write.csv(MM_table,file = paste0(out_path,"/MM_pMM_table.csv"),row.names = T)

    print("assign: geneModuleMembership, geneModuleMembership_Pvalue")
    assign("geneModuleMembership",value = geneModuleMembership,envir = .GlobalEnv)
    assign("geneModuleMembership_Pvalue",value = geneModuleMembership_Pvalue,envir = .GlobalEnv)
  }

  
  get_MM_out <- list(
    geneModuleMembership = geneModuleMembership,
    geneModuleMembership_Pvalue = geneModuleMembership_Pvalue
  )
  
  return(get_MM_out)
}



get_KME <- function(data,MEs,output){
  
  dataKME=signedKME(data, MEs)

  if (!is.null(out_path)) {
    write.csv(dataKME,file = paste0(out_path,"/KME_table.csv"),row.names = T)
    print("assign: dataKME")
    assign("dataKME",value = dataKME,envir = .GlobalEnv)    
    
  }
  return(dataKME)

}





# 计算基因和trait的相关性 即GS  在get_gs_mm_plot函数内被调用
get_GS <- function(data,trait,trait_col,out_path){
  
  
  if (!dir.exists(paste0(out_path,"/GS_table/"))) {
    dir.create(paste0(out_path,"/GS_table/"))
  }
  
  nSamples <- nrow(data)
  
  trait_select <- trait[,trait_col,drop = F]
  
  geneTraitSignificance = as.data.frame(cor(data, trait_select, use = "p")) %>% 
    set_names(paste("GS.", trait_col , sep=""))
  head(geneTraitSignificance)
  
  
  geneTraitSignificance_pvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples)) %>% 
    set_names(paste("p.GS.", trait_col, sep=""))

  head(geneTraitSignificance_pvalue)
  
  GS_table <- bind_cols(geneTraitSignificance,geneTraitSignificance_pvalue)
  
  write.csv(GS_table,file = paste0(out_path,"/GS_table/","gene_",trait_col,"_GS.csv"),row.names = T)

  # print("assign: geneTraitSignificance, geneTraitSignificance_pvalue")
  assign("geneTraitSignificance",value = geneTraitSignificance,envir = .GlobalEnv)
  assign("geneTraitSignificance_pvalue",value = geneTraitSignificance_pvalue,envir = .GlobalEnv)
  
  get_GS_out <- list(
    
    geneTraitSignificance = geneTraitSignificance,
    geneTraitSignificance_pvalue = geneTraitSignificance_pvalue
  )
  
  return(get_GS_out)
}


# GS MM scatter plot
get_gs_mm_plot <- function(data,MEs,trait,moduleColors,significent_mod,moduleTraitPvalue,width = 9,height = 6,out_path){
  
  if (!dir.exists(paste0(out_path,"/plot/GS_MM_plot"))) {
    dir.create(paste0(out_path,"/plot/GS_MM_plot"))
    dir.create(paste0(out_path,"/plot/GS_MM_plot/plot_data"),recursive = T)
  }
  
  get_MM_out <- get_MM(data = data,MEs = MEs,moduleColors = moduleColors,out_path = NULL)
  geneModuleMembership <- get_MM_out$geneModuleMembership

  #循环模块
  for ( module in significent_mod ) {
    
    moduleGenes = moduleColors==module
    #循环trait
    for (trait_col in colnames(trait)) {
      #模块与性状pvalue>0.05 则跳过进入下一循环
      if (moduleTraitPvalue[rownames(moduleTraitPvalue) == paste0("ME",module),trait_col ] > 0.05 ) {
        next
      }
      
      get_GS_out <- get_GS(data = data,trait = trait,trait_col = trait_col,out_path = out_path)
      
      geneTraitSignificance <- get_GS_out$geneTraitSignificance
      
      
      gs_mm_plot_data <- bind_cols(geneModuleMembership[,paste0("MM.",module),drop = F],geneTraitSignificance) %>% 
        modify_if(is.numeric,abs) %>% 
        filter(moduleGenes) %>% 
        arrange(desc(.data[[paste0("MM.",module)]]),desc(.data[[colnames(geneTraitSignificance)]]) )
      

      write.csv(gs_mm_plot_data,file = paste0(out_path,"/plot/GS_MM_plot/plot_data/",
                                              module,"_",trait_col,".csv"),row.names = T)
      
      
      pdf(file = paste0(out_path,"/plot/GS_MM_plot/cor_",paste0(module,"_",trait_col,".pdf")),width = width,height = height)
      verboseScatterplot(abs(geneModuleMembership[moduleGenes,paste0("MM.",module)]),         
                         abs(geneTraitSignificance[moduleGenes, 1]),
                         xlab = paste("Module Membership in", module, "module"),
                         ylab = "Gene significance",
                         main = paste("Module membership vs. gene significance\n"),
                         cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module) 
      
      dev.off()

    }
    
  }
 
}




merge_GS_MM_table <- function(out_path){
  
  
  MM_table <- read.csv(paste0(out_path,"/MM_pMM_table.csv"),row.names = 1)
  
  GS_file <- list.files( paste0(out_path,"/GS_table"),full.names = T )
  
  GS_table_all <- map_dfc(GS_file,read.csv,row.names = 1)
  
  merge_GS_MM_table <- bind_cols(MM_table,GS_table_all)
  
  write.csv(merge_GS_MM_table,file = paste0(out_path,"/merge_GS_MM_table.csv"))
  
}




# 加载TOM矩阵
get_tom <- function(net){
  load(net$TOMFiles)
  TOM <- as.matrix(TOM)
  assign("TOM",value = TOM,envir = .GlobalEnv)
  # return(TOM)
}


# TOM热图
get_tom_heatmap <- function(net,TOM,nselect,out_path,width = 9){
  
  geneTree = net$dendrograms[[1]]
  moduleColors = labels2colors(net$colors)
  # load(net$TOMFiles)
  # TOM <- as.matrix(TOM)
  
  #TOM矩阵转为距离矩阵
  dissTOM = 1 - TOM
  plotTOM = dissTOM ^ 7
  diag(plotTOM) = NA
  
  
  #随机选择部分基因 全部基因图太大，容易卡住
  set.seed(10)
  nGenes <- ncol(data)
  select <- sample(nGenes,size = nselect)
  selectTOM <- dissTOM[select,select]
  selectTOM[1:4,1:4]
  
  
  #重新聚类
  select_tree <- hclust(as.dist(selectTOM),method = "average")
  selectColors <- moduleColors[select]
  selectColors %>%head()
  
  plot_select_TOM = selectTOM ^ 7
  diag(plot_select_TOM) <- NA

  
  pdf(file = paste0(out_path,"/plot/Network_heatmap_plot_selected_genes.pdf"),width = width)
  TOMplot(plot_select_TOM,select_tree,selectColors,main = "Network heatmap plot, selected genes",col=gplots::colorpanel(250,'red',"orange",'lemonchiffon'))
  dev.off()
  
}


# MES聚类和热图
get_mes_cor_heatmap <- function(net,out_path,width = 9,height = 6){
  
  
  MEs_col <- net$MEs
  colnames(MEs_col) <- paste0("ME", labels2colors(
    as.numeric(str_replace_all(colnames(MEs_col),"ME",""))))
  MEs_col <- orderMEs(MEs_col)
  head(MEs_col)
  
  
  pdf(paste0(out_path,"/plot/Eigengene heatmap.pdf"),width = width,height = height)
  
  plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap",
                        marDendro = c(3,3,2,4),
                        marHeatmap = c(3,4,2,2), plotDendrograms = T, 
                        xLabelsAngle = 90)
  
  dev.off()
  
}


# 计算adj matrix
get_adj_matrix <- function(data,power,cor_type = NULL,net_type,out_path){
  
  if (cor_type == "bicor") {

     adj_matrix <- adjacency(datExpr= data, power = power,type = net_type,corFnc = "bicor",corOptions = list(maxPOutliers= 0.05, robustY=FALSE) )
  }else{

    adj_matrix = adjacency(datExpr= data, power = power,type = net_type)
  }
  adj_matrix[1:3,1:3]

  save(adj_matrix,file = paste0(out_path,"/RData/adj_matrix.RData") )
  return(adj_matrix)
}


# 计算连接度 from adj_matrix 矩阵 
# merge_by 同上 基因id的列名
get_connectivity <- function(adj_matrix,moduleColors,gene_anno = NULL,merge_by,out_path){
  
  dir.create(paste0(out_path,"/intramodularconnectivity/"),recursive = T)

  connect <- intramodularConnectivity(adj_matrix, colors = moduleColors)
  head(connect)
  connect[,"module"] <- moduleColors

  if (!is.null(gene_anno)) {

    connect[,merge_by] <- rownames(connect)
    connect <- connect  %>% group_by(module) %>% arrange(desc(connect$kWithin))
    connect <- left_join(connect, gene_anno, by = merge_by)
    head(connect)
  }
  

  
  connect %>% group_by(module) %>% 
    group_walk(~ write.csv(.x, file = paste0(out_path,"./intramodularconnectivity/",.y$module, "_connectivity.csv") ,row.names = F)  )

  write.csv(connect,file = paste0(out_path,"intramodularconnectivity.csv"),row.names = T)

  
  save(connect,file = paste0(out_path,"/RData/connect.RData"))

  print(head(connect))
  return(connect)
  
}




# 信息总结
info_summary <- function() {

  if (is.null(filter_threshold)) {
    
    info = paste0("标准化 ",method,"\n",
                  "过滤方法  ",filter_gene_method,"靠前的top",keep_gene_num,"\n",
                  "样本数  ",nSamples,"\n",
                  "基因数  ",nGenes,"\n",

                  "相关系数计算方法   ",cor_type,"\n",
                  "power软阈值   ",power,"\n",
                  "网络类型  ",net_type,"\n",
                  "模块融合阈值  ",mergeCutHeight,"\n",
                  "核心基因筛选阈值  ","GS_threshold: ",GS_threshold,"\t","MM_threshold: ",MM_threshold,"\n",
                  "输出路径  ",out_path,"\n"
                     )

    print(info)
    write.table(info,quote = F,row.names = F,col.names = F,file = paste0(out_path,"/info_summary.txt"),)
  }
  else{

    info = paste0(   "标准化  ",method,"\n",
                     "过滤方法  ",filter_gene_method,"大于",filter_threshold,"\n",
                     "样本数  ",nSamples,"\n",
                     "基因数  ",nGenes,"\n",

                     "相关系数计算方法   ",cor_type,"\n",
                     "power软阈值   ",power,"\n",
                     "网络类型  ",net_type,"\n",
                     "模块融合阈值  ",mergeCutHeight,"\n",
                     "核心基因筛选阈值  ","GS_threshold: ",GS_threshold, "\t" ,"MM_threshold: ",MM_threshold,"\n",
                     "输出路径  ",out_path,"\n"
                     )


    print(info)

    write.table(info,quote = F,row.names = F,col.names = F,file = paste0(out_path,"/info_summary.txt"))

  }



}




# 导出模块基因保存为cytoscape的输入格式
module_gene_to_cyto <- function(TOM,moduleColors,significent_mod,out_path){
  
  dir.create(paste0(out_path,"/Cytoscape/module_gene_to_Cytoscape"),recursive = T)
  
  for (module in significent_mod) {
    
    genenames = colnames(data)
    inModule = (moduleColors == module)
    modgene = genenames[inModule]
    
    modTOM = TOM[inModule, inModule]
    dim(modTOM)
    dimnames(modTOM) = list(modgene, modgene)
    
    # Export the network into edge and node list files Cytoscape can read
    cyt = exportNetworkToCytoscape(modTOM,
                                   edgeFile = paste(out_path,"/Cytoscape/module_gene_to_Cytoscape/CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
                                   # nodeFile = paste(out_path,"module_gene_to_Cytoscape/CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
                                   weighted = TRUE,
                                   threshold = 0.02,
                                   nodeNames = modgene,                               
                                   # nodeAttr = moduleColors[inModule]
                                   
                                   )
    
    
  }
  
}




# 核心基因的筛选
# 方法1
# 根据GS 和 MM  
get_hub_gene_gsmm <- function(data,MEs,trait,moduleColors,
                              significent_mod,moduleTraitPvalue,
                              GS_threshold,MM_threshold,out_path){
  
  
  dir.create(paste0(out_path,"hubgene_GS",GS_threshold,"_MM",MM_threshold,"/"))
  dir.create(paste0(out_path,"hubgene_GS",GS_threshold,"_MM",MM_threshold,"/","txt"))
  dir.create(paste0(out_path,"hubgene_GS",GS_threshold,"_MM",MM_threshold,"/","csv"))
  dir.create(paste0(out_path,"hubgene_GS",GS_threshold,"_MM",MM_threshold,"/","plot"))
  
  list_gs_mm_hubgene <- list()
  
  get_MM_out <- get_MM(data = data,MEs = MEs,moduleColors = moduleColors,out_path = out_path)
  geneModuleMembership <- get_MM_out$geneModuleMembership
  
  #循环模块
  for ( module in significent_mod ) {
    
    moduleGenes = moduleColors==module
    #循环trait
    for (trait_col in colnames(trait)) {
      #模块与性状pvalue>0.05 则跳过进入下一循环
      if (moduleTraitPvalue[rownames(moduleTraitPvalue) == paste0("ME",module),trait_col ] > 0.05 ) {
        next
      }
      
      get_GS_out <- get_GS(data = data,trait = trait,trait_col = trait_col,out_path = out_path)
      
      geneTraitSignificance <- get_GS_out$geneTraitSignificance
      
      # index <- abs(geneModuleMembership[,paste0("MM",module)])>MM_threshold &
      #   abs(geneTraitSignificance[, 1])>GS_threshold &
      #   moduleGenes
      # hubgene <- colnames(data)[index]
      
      hubgene_gsmm <- bind_cols(geneModuleMembership[,paste0("MM.",module),drop = F],geneTraitSignificance) %>% 
        modify_if(is.numeric,abs) %>% 
        filter(moduleGenes) %>% 
        filter(.data[[paste0("MM.",module)]]>MM_threshold,.data[[colnames(geneTraitSignificance)]] > GS_threshold) %>% 
        arrange(desc(.data[[paste0("MM.",module)]]),desc(.data[[colnames(geneTraitSignificance)]]) )
      
      
      hubgene_gsmm <- 
        hubgene_gsmm %>% mutate(ensembl_gene_id = rownames(hubgene_gsmm)) %>% left_join(gene_anno) %>% column_to_rownames("ensembl_gene_id")
      
      
      
      
      write.csv(hubgene_gsmm,file = paste0(out_path,"hubgene_GS",GS_threshold,"_MM",MM_threshold,"/","csv/",module,"_",trait_col,".csv"),row.names = T)
      write.table(rownames(hubgene_gsmm),
                  file = paste0(out_path,"hubgene_GS",GS_threshold,"_MM",MM_threshold,"/","txt/",module,"_",trait_col,".txt"),quote = F,row.names = F,col.names = F)
      
      
      list_name <- paste0(module,"_",trait_col)
      list_gs_mm_hubgene[[list_name]] <- rownames(hubgene_gsmm)
      
      
      
      pdf(file = paste0(out_path,"hubgene_GS",GS_threshold,"_MM",MM_threshold,"/","plot/",paste0(module,"_",trait_col,".pdf")),width = 9,height = 6)
      verboseScatterplot(abs(geneModuleMembership[moduleGenes,paste0("MM.",module)]),         
                         abs(geneTraitSignificance[moduleGenes, 1]),
                         xlab = paste("Module Membership in", module, "module"),
                         ylab = "Gene significance",
                         main = paste("Module membership vs. gene significance\n"),
                         cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module) 
      
      abline(v = MM_threshold, col = "red", lty = 2,lwd = 2)
      abline(h = GS_threshold, col = "red", lty = 2,lwd = 2)
      dev.off()
      
    }
    
  }
  
  
  assign(paste0("list_gs_mm_hubgene","_GS",GS_threshold,"_MM",MM_threshold),value = list_gs_mm_hubgene,envir = .GlobalEnv)
  
  save(list_gs_mm_hubgene,file = paste0(out_path,"RData/list_gs_mm_hubgene","_GS",GS_threshold,"_MM",MM_threshold,".RData"))
  
  
  return(list_gs_mm_hubgene)
  
}






#(多余的)
get_hub_gene_kme <- function(data,MEs,moduleColors,top_num,
                              significent_mod,
                              out_path){
  
  
  dir.create(paste0(out_path,"hubgene_KME"))
  dir.create(paste0(out_path,"hubgene_KME/txt"))
  dir.create(paste0(out_path,"hubgene_KME/csv"))
  
  
  list_KME_hubgene <- list()
  
  dataKME <- get_KME(data = data,MEs = MEs,output = NULL)
  
  
  #循环模块
  for ( module in significent_mod ) {
    
    moduleGenes = moduleColors==module
    
    
    module_KME_abs <- dataKME[moduleGenes,paste0("kME",module),drop = F] %>% 
      modify_at(1,abs) %>% 
      arrange(desc(.data[[paste0("kME",module)]]))
      
    
    module_KME_abs_top <- module_KME_abs %>% 
      .[1:top_num,,drop = F]
    
    
    list_KME_hubgene[[module]] <- rownames(module_KME_abs_top)
    
    write.csv(module_KME_abs_top,file = paste0(out_path,"hubgene_KME/csv/","top",top_num,"_",module,".csv"),row.names = T)
    write.table(rownames(module_KME_abs_top),
                file = paste0(out_path,"hubgene_KME/txt/","top",top_num,"_",module,".txt"),quote = F,row.names = F,col.names = F)
    
    
  }
  
  save(list_KME_hubgene,file = paste0(out_path,"RData/list_KME_hubgene.RData"))
  return(list_KME_hubgene)
}
    
  

# 方法2(有待改进，暂时别用)
# 根据 GS 和 模块内连接度
get_hub_gene_gs_connect <- function(connect,data,trait,moduleColors,
                                    significent_mod,moduleTraitPvalue,ntop = 10,out_path,
                                    width = 9,height = 6){
  
  dir.create(paste0(out_path,"hubgene_GS_connect"))
  dir.create(paste0(out_path,"hubgene_GS_connect/txt"))
  dir.create(paste0(out_path,"hubgene_GS_connect/csv"))
  dir.create(paste0(out_path,"hubgene_GS_connect/plot"))
  dir.create(paste0(out_path,"hubgene_GS_connect/plot_data"))
  
  
  list_gs_conn_hubgene <- list()
    
  for ( module in significent_mod ) {
    
    moduleGenes = moduleColors==module
    #循环trait
    for (trait_col in colnames(trait)) {
      #模块与性状pvalue>0.05 则跳过进入下一循环
      if (moduleTraitPvalue[rownames(moduleTraitPvalue) == paste0("ME",module),trait_col ] > 0.05 ) {
        next
      }
      
      get_GS_out <- get_GS(data = data,trait = trait,trait_col = trait_col,out_path = out_path)
      geneTraitSignificance <- get_GS_out$geneTraitSignificance

      
      gs_conn_data <- bind_cols(connect[moduleGenes,"kWithin",drop = F],geneTraitSignificance[moduleGenes,,drop = F]) %>% 
        modify_if(is.numeric,abs) %>% 
        arrange(desc(kWithin) )
      
      write.csv(gs_conn_data,file = paste0(out_path,"hubgene_GS_connect/plot_data/",module,"_",trait_col,".csv"),row.names = T)
      
      
      hubgene_gs_conn <- gs_conn_data %>% filter(.data[[colnames(geneTraitSignificance)]] > GS_threshold) %>% 
        slice_max(order_by = kWithin,n = ntop)
      
      write.csv(hubgene_gs_conn,file = paste0(out_path,"hubgene_GS_connect/csv/",module,"_",trait_col,".csv"),row.names = T)
      
      
      write.table(rownames(hubgene_gs_conn),
                  file = paste0(out_path,"hubgene_GS_connect/txt/",module,"_",trait_col,".txt"),quote = F,row.names = F,col.names = F)
      
      
      list_name <- paste0(module,"_",trait_col)
      list_gs_conn_hubgene[[list_name]] <-  rownames(hubgene_gs_conn)
      
      pdf(paste0(out_path,"hubgene_GS_connect/plot/",module,"_",trait_col,".pdf"),width = width,height = height)
      verboseScatterplot(connect$kWithin[moduleGenes],
                         abs(geneTraitSignificance[moduleGenes,]), col=module,
                         main=module,
                         xlab = "Connectivity", ylab = "Gene Significance" , abline = TRUE
                         )
      
      dev.off()

    }
  
  }
  
  assign("list_gs_conn_hubgene",value = list_gs_conn_hubgene,envir = .GlobalEnv)
  save(list_gs_conn_hubgene,file = paste0(out_path,"RData/list_gs_conn_hubgene.RData"))
  
  return(list_gs_conn_hubgene)

}



# 将核心基因导出为cytoscape的输入格式
hub_gene_to_cyto <- function(data,TOM,hubgene_list,get_hubgene_method,out_path){
  
  dir.create(paste0(out_path,"Cytoscape/",get_hubgene_method,"_hubgene_to_Cytoscape"),recursive =T)
  
  
  dimnames(TOM) <- list(colnames(data),colnames(data))
  
  for (i in 1:length(hubgene_list)) {
    
    hubgene_list_name <- names(hubgene_list)[i]
    hubgene_id <- hubgene_list[[i]]
    
    hubgene_TOM <- TOM[hubgene_id,hubgene_id]
    
    
    cyt_top = exportNetworkToCytoscape(hubgene_TOM,
                                       edgeFile = paste(out_path,"Cytoscape/",get_hubgene_method,"_hubgene_to_Cytoscape/","CytoscapeInput-edges-",hubgene_list_name , "_hubgene.txt", sep=""),
                                       # nodeFile = paste(out_path,get_hubgene_method,"_hubgene_to_Cytoscape/","CytoscapeInput-nodes-", hubgene_list_name , "_hubgene.txt", sep=""),
                                       weighted = TRUE,
                                       threshold = 0.02,
                                       nodeNames = hubgene_id,                               
                                       # #altNodeNames = modGenes,
                                       # nodeAttr = moduleColors[FilterGenes]
                                       
                                       )
    
    
  }

}


