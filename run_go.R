
# version 20220921

# @description: 批量做GO富集分析

# @parameters: genelist genelist 要有names
# @parameters: plot_type 生成图片类型 有 dot tree directed_acyclic
# @parameters: simplify simplify = T 时对生成的dotplot 减少显示相似term

# @parameters: keyType 输入的基因类型
# @parameters: OrgDb 数据库选择 such as org.Hs.eg.db, org.Mm.eg.db 
# @parameters: use_semData nCluster treeplot的需要的参数
# @parameters: outputpath 输出路径



run_go <- function(genelist,use_go_RData = F,go_level = NULL,plot_type,simplify = F,use_semData = F,nCluster = 5,outputpath,xlsxfilename = "",keyType ="ENSEMBL",OrgDb = "org.Hs.eg.db")
{ 
  
  library(conflicted)
  library(openxlsx)
  library(clusterProfiler)
  library(stringr)
  library(enrichplot)
  library(GOSemSim)
  library(ggnewscale)         # for tree
  library(readr)              # write_csv()
  library(ggplot2)
  # library(topGO)
  conflict_prefer("simplify", "clusterProfiler")
  
  if (!dir.exists(outputpath)) {
    dir.create(outputpath,recursive = T)
  }
  
  
  if (!use_go_RData) {
    
    get_go_result(genelist = genelist,xlsxfilename = xlsxfilename,outputpath = outputpath,OrgDb = OrgDb,
                  keyType = keyType)
    
  }
  
  
  
  for (name in names(genelist)) {
    
    load( paste0(outputpath,"/go_enrich_RData/",name,"_enrich_go.RData") )
    
     
    if (plot_type == "dot" & is.null(go_level) ) {
      go_dot_plot(go_enrich = go_enrich,name = name,outputpath = outputpath,simplify = simplify)
      
      go_dot_plot_005_01(go_enrich = go_enrich,name = name,outputpath = outputpath)
      
     
    }
    
    
    if (plot_type == "dot" & (!is.null(go_level)) ) {
      go_dot_plot_level(go_enrich = go_enrich,name = name,outputpath = outputpath,level = go_level)
    }
    
    
    if (plot_type == "directed_acyclic") {
      
      go_directed_acyclic_plot(go_enrich = go_enrich,name = name,outputpath = outputpath)

    }
    
    
    if (plot_type == "tree") {
      
      go_tree_plot(go_enrich = go_enrich,name = name,use_semData = use_semData , outputpath = outputpath,nCluster = nCluster,OrgDb = OrgDb)
      
      go_tree_plot_005_01(go_enrich = go_enrich,name = name,use_semData = use_semData, outputpath = outputpath,nCluster = nCluster,OrgDb = OrgDb)
      
    }
    
    
    if (plot_type == "all") {
      
      go_dot_plot(go_enrich = go_enrich,name = name,outputpath = outputpath)
      
      go_dot_plot_level(go_enrich = go_enrich,name = name,outputpath = outputpath,level = go_level)
      
      go_directed_acyclic_plot(go_enrich = go_enrich,name = name,outputpath = outputpath)
      
      go_tree_plot(go_enrich = go_enrich,name = name,outputpath = outputpath)
      
    }
    

  }
  

}




#==================================




get_go_result <-  function(genelist = genelist,outputpath = outputpath,xlsxfilename = xlsxfilename,keyType =keyType,OrgDb){
  
  
  if (!dir.exists(paste0(outputpath,"/go_enrich_RData/"))) {
    dir.create(paste0(outputpath,"/go_enrich_RData/"),recursive = T)
  }
  
  
  wb <- createWorkbook()
  wb2 <- createWorkbook()
  
  
  for (i in 1:length(genelist)) {
    
    name <- names(genelist[i])
    gene <- as.data.frame(genelist[i])
    gene <- gene[,1]
    
    print(c(i,name))
    
    if(str_detect(OrgDb,"Hs") & keyType == "ENSEMBL"){
      gene <- grep("^ENSG",gene,value = T)
    }

    
    
    
    if (keyType != "SYMBOL") {
      
      go_enrich <- enrichGO(gene          = gene,
                           OrgDb         = OrgDb,
                           keyType       = keyType,
                           ont           = "all",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 1,
                           qvalueCutoff  = 1,
                           readable = T
      )
      
    }else{
      
      go_enrich<- enrichGO(gene          = gene,
                           OrgDb         = OrgDb,
                           keyType       = keyType,
                           ont           = "all",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 1,
                           qvalueCutoff  = 1
      )      
    }
    
    
    save(go_enrich,file = paste0(outputpath,"/go_enrich_RData/",name,"_enrich_go.RData") )
    
    
    if (is.null(go_enrich)) {
      addWorksheet(wb,sheetName = name)
      addWorksheet(wb2,sheetName = name)
      
      next
    }
    
    
    go_result <- go_enrich@result
    
    addWorksheet(wb,sheetName = name)
    writeData(wb,i,go_result,colNames = T,rowNames = T)
    
    filter_go_all <- go_enrich[go_enrich@result$p.adjust<0.05,asis=T]
    addWorksheet(wb2,sheetName = name)
    writeData(wb2,i,filter_go_all,colNames = T,rowNames = T)
    
    saveWorkbook(wb,file = paste0(outputpath,"/GO_enrichment_result_all_",xlsxfilename,".xlsx"),overwrite = T)
    saveWorkbook(wb2,file = paste0(outputpath,"/GO_enrichment_result_padj005_",xlsxfilename,".xlsx"),overwrite = T)
  
  }
  
}



go_dot_plot <- function(go_enrich,name,outputpath,simplify){
  
  
  simp<- switch (simplify,
                 T = "_simplify",
                 F = ""
  )
  
  
  if (!dir.exists(paste0(outputpath,"/dotplot",simp))) {
    dir.create(paste0(outputpath,"/dotplot",simp),recursive = T)
  }
  
  
  
 
  
  filter_go_BP <- go_enrich[go_enrich@result$p.adjust<0.05 & go_enrich@result$ONTOLOGY == "BP",asis=T]
  filter_go_MF <- go_enrich[go_enrich@result$p.adjust<0.05 & go_enrich@result$ONTOLOGY == "MF",asis=T]
  filter_go_CC <- go_enrich[go_enrich@result$p.adjust<0.05 & go_enrich@result$ONTOLOGY == "CC",asis=T]

  
  
  if (length(filter_go_BP@result$ID) >= 1 ) {
    
    
    if (simplify == T) {
      filter_go_BP <- simplify(filter_go_BP)
    }
    
    
    p <- dotplot(filter_go_BP,showCategory = length(filter_go_BP@result$ID),x = "Count") +scale_size(range = c(4,8)) +scale_y_discrete(labels= function(x) str_wrap(x,width = 45))
    
    # go_dotplot_list[[paste0(name,"_bp")]] = p
    
    pdf(file =  paste0(outputpath,"/dotplot",simp,"/BP_dotplot_",name,".pdf"),width = 9,height = 5+(length(filter_go_BP@result$ID)*0.2))
    print(p)
    dev.off()
  }
  
  
  if (length(filter_go_MF@result$ID) >= 1 ) {
    
    
    
    if (simplify == T) {
      filter_go_MF <- simplify(filter_go_MF)
    }
    

    
    
    p <- dotplot(filter_go_MF,showCategory = length(filter_go_MF@result$ID),x = "Count") +scale_size(range = c(4,8)) +scale_y_discrete(labels= function(x) str_wrap(x,width = 45))
    
    # go_dotplot_list[[paste0(name,"_mf")]] = p
    
    pdf(file =  paste0(outputpath,"/dotplot",simp,"/MF_dotplot_",name,".pdf"),width = 9,height = 5+(length(filter_go_MF@result$ID)*0.2))
    print(p)
    dev.off()
  }
  

  if (length(filter_go_CC@result$ID) >= 1 ) {
    
    
    if (simplify == T) {
      filter_go_CC <- simplify(filter_go_CC)
    }
    

    
    p <- dotplot(filter_go_CC,showCategory = length(filter_go_CC@result$ID),x = "Count") +scale_size(range = c(4,8)) +scale_y_discrete(labels= function(x) str_wrap(x,width = 45))
    
    # go_dotplot_list[[paste0(name,"_cc")]] = p
    
    pdf(file =  paste0(outputpath,"/dotplot",simp,"/CC_dotplot_",name,".pdf"),width = 9,height = 5+(length(filter_go_CC@result$ID)*0.2))
    print(p)
    dev.off()
  }
  
  
  
  if (simplify == T) {
    
    if (!dir.exists(paste0(outputpath,"dotplot",simp,"/A-csv/"))) {
      dir.create(paste0(outputpath,"dotplot",simp,"/A-csv/"),recursive = T)
    }
    
    bind_rows(filter_go_BP@result,filter_go_MF@result,filter_go_CC@result) %>%
      write.csv(file = paste0(outputpath,"dotplot",simp,"/A-csv/",name,".csv"),row.names = F)
  }
  

  
  
}



go_dot_plot_level <- function(go_enrich,name,level,outputpath){
  
  
  
  if (!dir.exists(paste0(outputpath,"/dotplot_level_",level))) {
    dir.create(paste0(outputpath,"/dotplot_level_",level),recursive = T)
  }
  
  if (!dir.exists(paste0(outputpath,"/dotplot_level_",level,"/A-csv"))) {
    dir.create(paste0(outputpath,"/dotplot_level_",level,"/A-csv"),recursive = T)
  }
  
  
  
  
  filter_go_BP <- go_enrich[go_enrich@result$p.adjust<0.05 & go_enrich@result$ONTOLOGY == "BP",asis=T]
  filter_go_MF <- go_enrich[go_enrich@result$p.adjust<0.05 & go_enrich@result$ONTOLOGY == "MF",asis=T]
  filter_go_CC <- go_enrich[go_enrich@result$p.adjust<0.05 & go_enrich@result$ONTOLOGY == "CC",asis=T]
  
  filter_go_BP@ontology <- "BP"
  filter_go_MF@ontology <- "MF"
  filter_go_CC@ontology <- "CC"
  
  
  filter_go_BP <- gofilter(filter_go_BP, level=level)
  filter_go_MF <- gofilter(filter_go_MF, level=level)
  filter_go_CC <- gofilter(filter_go_CC, level=level)
  
  
  bind_rows(filter_go_BP@result,filter_go_MF@result,filter_go_CC@result) %>%
    write.csv(file = paste0(outputpath,"dotplot_level_",level,"/A-csv/",name,".csv"),row.names = F)
    
  
  if (length(filter_go_BP@result$ID) >= 1 ) {
    
    p <- dotplot(filter_go_BP,showCategory = length(filter_go_BP@result$ID),x = "Count") +scale_size(range = c(4,8)) +scale_y_discrete(labels= function(x) str_wrap(x,width = 45))
    
    # go_dotplot_list[[paste0(name,"_bp")]] = p
    
    pdf(file =  paste0(outputpath,"/dotplot_level_",level,"/BP_dotplot_",name,".pdf"),width = 9,height = 5+(length(filter_go_BP@result$ID)*0.2))
    print(p)
    dev.off()
  }
  
  
  if (length(filter_go_MF@result$ID) >= 1 ) {
    
    p <- dotplot(filter_go_MF,showCategory = length(filter_go_MF@result$ID),x = "Count") +scale_size(range = c(4,8)) +scale_y_discrete(labels= function(x) str_wrap(x,width = 45))
    
    # go_dotplot_list[[paste0(name,"_mf")]] = p
    
    pdf(file =  paste0(outputpath,"/dotplot_level_",level,"/MF_dotplot_",name,".pdf"),width = 9,height = 5+(length(filter_go_MF@result$ID)*0.2))
    print(p)
    dev.off()
  }
  
  
  if (length(filter_go_CC@result$ID) >= 1 ) {
    
    p <- dotplot(filter_go_CC,showCategory = length(filter_go_CC@result$ID),x = "Count") +scale_size(range = c(4,8)) +scale_y_discrete(labels= function(x) str_wrap(x,width = 45))
    
    # go_dotplot_list[[paste0(name,"_cc")]] = p
    
    pdf(file =  paste0(outputpath,"/dotplot_level_",level,"/CC_dotplot_",name,".pdf"),width = 9,height = 5+(length(filter_go_CC@result$ID)*0.2))
    print(p)
    dev.off()
    
  }
  
  
  
  
  
  
}



go_directed_acyclic_plot <- function(go_enrich,name,outputpath){
  
  
  
  if (!dir.exists(paste0(outputpath,"/directed_acyclic_plot"))) {
    dir.create(paste0(outputpath,"/directed_acyclic_plot"),recursive = T)
  }
  
  
  
  filter_go_BP <- go_enrich[go_enrich@result$p.adjust<0.05 & go_enrich@result$ONTOLOGY == "BP",asis=T]
  filter_go_MF <- go_enrich[go_enrich@result$p.adjust<0.05 & go_enrich@result$ONTOLOGY == "MF",asis=T]
  filter_go_CC <- go_enrich[go_enrich@result$p.adjust<0.05 & go_enrich@result$ONTOLOGY == "CC",asis=T]
  
  
  filter_go_BP@ontology <- "BP"
  filter_go_MF@ontology <- "MF"
  filter_go_CC@ontology <- "CC"
  
  if (length(filter_go_BP@result$ID) >= 1 ) {
    
    pdf(file =  paste0(outputpath,"/directed_acyclic_plot","/directed_acyclic_BP_",name,".pdf"),width = 9,height = 9)
    plotGOgraph(filter_go_BP)
    dev.off()
  }
  
  
  if (length(filter_go_MF@result$ID) >= 1 ) {
    
    pdf(file =  paste0(outputpath,"/directed_acyclic_plot","/directed_acyclic_MF_",name,".pdf"),width = 9,height = 9)
    plotGOgraph(filter_go_MF)
    dev.off()
  }
  
  
  if (length(filter_go_CC@result$ID) >= 1 ) {
    
    pdf(file =  paste0(outputpath,"/directed_acyclic_plot","/directed_acyclic_CC_",name,".pdf"),width = 9,height = 9)
    plotGOgraph(filter_go_CC)
    dev.off()
    
  }

  
}



go_tree_plot <- function(go_enrich,name,use_semData,nCluster = 5,outputpath,OrgDb){
  
  
  # if (!dir.exists(paste0(outputpath,"/tree_plot"))) {
  #   dir.create(paste0(outputpath,"/tree_plot"),recursive = T)
  # }
  
  if (use_semData) {
    
    dir_end = "_use_semData"
    dir.create(paste0(outputpath,"/tree_plot","_nCluster_",nCluster,dir_end),recursive = T)
    
  }else{
    dir_end = ""
    dir.create(paste0(outputpath,"/tree_plot","_nCluster_",nCluster,dir_end),recursive = T)
  }
  
  
  
  
  filter_go_BP <- go_enrich[go_enrich@result$p.adjust<0.05 & go_enrich@result$ONTOLOGY == "BP",asis=T]
  filter_go_MF <- go_enrich[go_enrich@result$p.adjust<0.05 & go_enrich@result$ONTOLOGY == "MF",asis=T]
  filter_go_CC <- go_enrich[go_enrich@result$p.adjust<0.05 & go_enrich@result$ONTOLOGY == "CC",asis=T]
  
  
  filter_go_BP@ontology <- "BP"
  filter_go_MF@ontology <- "MF"
  filter_go_CC@ontology <- "CC"
  
  if (use_semData) {
    
    if (!exists("hsGO_BP",)) {
      hsGO_BP<<- godata(OrgDb, ont="BP")
    }
    if (!exists("hsGO_MF")) {
      hsGO_MF <<- godata(OrgDb, ont="MF")
    }
    if (!exists("hsGO_CC")) {
      hsGO_CC <<- godata(OrgDb, ont="CC")
    }

  }

  
  
  if (length(filter_go_BP@result$ID) >= 1 ) {
    
    
    if (use_semData) {
      x <- pairwise_termsim(filter_go_BP,semData = hsGO_BP,method = "Wang")
    }else{
      x <- pairwise_termsim(filter_go_BP)
    }
    

    tryCatch(expr = {p <- treeplot(x,showCategory = length(filter_go_BP@result$ID),nCluster = nCluster)+coord_cartesian(clip = "off")
    pdf(file =  paste0(outputpath,"/tree_plot","_nCluster_",nCluster,dir_end,"/tree_plot_BP_",name,".pdf"),width = 12,height = 5+(length(filter_go_BP@result$ID)*0.2))
    print(p)
    dev.off()
                    },
             error = function(e){
               print(e)
             })
    

  }
  
  
  if (length(filter_go_MF@result$ID) >= 1 ) {
    
    
    if (use_semData) {
      x <- pairwise_termsim(filter_go_MF,semData = hsGO_MF,method = "Wang")
    }else{
      x <- pairwise_termsim(filter_go_MF)
    }


    tryCatch(expr = {
      p <- treeplot(x,showCategory = length(filter_go_MF@result$ID),nCluster = nCluster)+coord_cartesian(clip = "off")
      pdf(file =  paste0(outputpath,"/tree_plot","_nCluster_",nCluster,dir_end,"/tree_plot_MF_",name,".pdf"),width = 12,height = 5+(length(filter_go_MF@result$ID)*0.2))
      print(p)
      dev.off()
      
    },
             error = function(e){
               print(e)
             })
    

  }
  
  
  if (length(filter_go_CC@result$ID) >= 1 ) {
    

    if (use_semData) {
      x <- pairwise_termsim(filter_go_CC,semData = hsGO_CC,method = "Wang")
    }else{
      x <- pairwise_termsim(filter_go_CC)
    }
    
    
    tryCatch(expr = {
      
      p <- treeplot(x,showCategory = length(filter_go_CC@result$ID),nCluster = nCluster)+coord_cartesian(clip = "off")
      pdf(file =  paste0(outputpath,"/tree_plot","_nCluster_",nCluster,dir_end,"/tree_plot_CC_",name,".pdf"),width = 12,height = 5+(length(filter_go_CC@result$ID)*0.2))
      print(p)
      dev.off()
      
      
    },
             error = function(e){
               print(e)
             })
    

    
  }
  
  
  
  
  
}



go_dot_plot_005_01 <- function(go_enrich,name,outputpath,OrgDb){
  
  
  
  if (!dir.exists(paste0(outputpath,"/dotplot_005-01"))) {
    dir.create(paste0(outputpath,"/dotplot_005-01"),recursive = T)
  }
  
  
  
  if (file.exists( paste0(outputpath,"/dotplot_005-01/A-dotplot_",name,".csv"))) {

    file.remove( paste0(outputpath,"/dotplot_005-01/A-dotplot_",name,".csv") )

  }
  
  
  
  filter_go_BP <- go_enrich[go_enrich@result$p.adjust<0.1 & go_enrich@result$ONTOLOGY == "BP",asis=T]
  filter_go_MF <- go_enrich[go_enrich@result$p.adjust<0.1 & go_enrich@result$ONTOLOGY == "MF",asis=T]
  filter_go_CC <- go_enrich[go_enrich@result$p.adjust<0.1 & go_enrich@result$ONTOLOGY == "CC",asis=T]
  
  
  
  if (length(filter_go_BP@result$ID) >= 1 & min(filter_go_BP@result$p.adjust) >0.05 ) {
    
    p <- dotplot(filter_go_BP,showCategory = length(filter_go_BP@result$ID),x = "Count") +scale_size(range = c(4,8)) +scale_y_discrete(labels= function(x) str_wrap(x,width = 45))
    
    # go_dotplot_list[[paste0(name,"_bp")]] = p
    
    
    write_csv(as.data.frame(filter_go_BP@result),file = paste0(outputpath,"/dotplot_005-01/A-dotplot_",name,".csv"),append = T)
    
    
    pdf(file =  paste0(outputpath,"/dotplot_005-01/BP_dotplot_",name,".pdf"),width = 9,height = 5+(length(filter_go_BP@result$ID)*0.2))
    print(p)
    dev.off()
    
    
  }
  
  
  if (length(filter_go_MF@result$ID) >= 1 & min(filter_go_MF@result$p.adjust) >0.05 ) {
    
    p <- dotplot(filter_go_MF,showCategory = length(filter_go_MF@result$ID),x = "Count") +scale_size(range = c(4,8)) +scale_y_discrete(labels= function(x) str_wrap(x,width = 45))
    
    # go_dotplot_list[[paste0(name,"_mf")]] = p
    
    
    write_csv(as.data.frame(filter_go_MF@result),file = paste0(outputpath,"/dotplot_005-01/A-dotplot_",name,".csv"),append = T)
    
    pdf(file =  paste0(outputpath,"/dotplot_005-01/MF_dotplot_",name,".pdf"),width = 9,height = 5+(length(filter_go_MF@result$ID)*0.2))
    print(p)
    dev.off()
  }
  
  
  if (length(filter_go_CC@result$ID) >= 1 & min(filter_go_BP@result$p.adjust) >0.05 ) {
    
    p <- dotplot(filter_go_CC,showCategory = length(filter_go_CC@result$ID),x = "Count") +scale_size(range = c(4,8)) +scale_y_discrete(labels= function(x) str_wrap(x,width = 45))
    
    # go_dotplot_list[[paste0(name,"_cc")]] = p
    
    
    write_csv(as.data.frame(filter_go_CC@result),file = paste0(outputpath,"/dotplot_005-01/A-dotplot_",name,".csv"),append = T)
    
    
    pdf(file =  paste0(outputpath,"/dotplot_005-01/CC_dotplot_",name,".pdf"),width = 9,height = 5+(length(filter_go_CC@result$ID)*0.2))
    print(p)
    dev.off()
  }
  
}



go_tree_plot_005_01 <- function(go_enrich,name,use_semData,nCluster = 5,outputpath,OrgDb){
  
  
  # if (!dir.exists(paste0(outputpath,"/tree_plot_005-01"))) {
  #   dir.create(paste0(outputpath,"/tree_plot_005-01"),recursive = T)
  # }
  
  
  if (use_semData) {
    
    dir_end = "_use_semData"
    dir.create(paste0(outputpath,"/tree_plot_005-01","_nCluster_",nCluster,dir_end),recursive = T)
    
  }else{
    dir_end = ""
    dir.create(paste0(outputpath,"/tree_plot_005-01","_nCluster_",nCluster,dir_end),recursive = T)
  }
  
  
  
  
  filter_go_BP <- go_enrich[go_enrich@result$p.adjust<0.1 & go_enrich@result$ONTOLOGY == "BP",asis=T]
  filter_go_MF <- go_enrich[go_enrich@result$p.adjust<0.1 & go_enrich@result$ONTOLOGY == "MF",asis=T]
  filter_go_CC <- go_enrich[go_enrich@result$p.adjust<0.1 & go_enrich@result$ONTOLOGY == "CC",asis=T]
  
  
  filter_go_BP@ontology <- "BP"
  filter_go_MF@ontology <- "MF"
  filter_go_CC@ontology <- "CC"
  
  
  if (use_semData) {
    
    if (!exists("hsGO_BP",)) {
      hsGO_BP<<- godata(OrgDb, ont="BP")
    }
    if (!exists("hsGO_MF")) {
      hsGO_MF <<- godata(OrgDb, ont="MF")
    }
    if (!exists("hsGO_CC")) {
      hsGO_CC <<- godata(OrgDb, ont="CC")
    }

  }  

  
  
  if (length(filter_go_BP@result$ID) >= 1 & min(filter_go_BP@result$p.adjust) >0.05 ) {
    
    
    if (use_semData) {
      x <- pairwise_termsim(filter_go_BP,semData = hsGO_BP,method = "Wang")
    }else{
      x <- pairwise_termsim(filter_go_BP)
    }
    
    
    tryCatch(expr = {p <- treeplot(x,showCategory = length(filter_go_BP@result$ID),nCluster = nCluster)+coord_cartesian(clip = "off")
    pdf(file =  paste0(outputpath,"/tree_plot_005-01","_nCluster_",nCluster,dir_end,"/tree_plot_BP_",name,".pdf"),width = 12,height = 5+(length(filter_go_BP@result$ID)*0.2))
    print(p)
    dev.off()
    },
    error = function(e){
      e
    })
    
    
  }
  
  
  if (length(filter_go_MF@result$ID) >= 1 & min(filter_go_MF@result$p.adjust) >0.05 ) {
    
    
    if (use_semData) {
      x <- pairwise_termsim(filter_go_MF,semData = hsGO_MF,method = "Wang")
    }else{
      x <- pairwise_termsim(filter_go_MF)
    }
    
    
    tryCatch(expr = {
      p <- treeplot(x,showCategory = length(filter_go_MF@result$ID),nCluster = nCluster)+coord_cartesian(clip = "off")
      pdf(file =  paste0(outputpath,"/tree_plot_005-01","_nCluster_",nCluster,dir_end,"/tree_plot_MF_",name,".pdf"),width = 12,height = 5+(length(filter_go_MF@result$ID)*0.2))
      print(p)
      dev.off()
      
    },
    error = function(e){
      e
    })
    
    
  }
  
  
  if (length(filter_go_CC@result$ID) >= 1 & min(filter_go_CC@result$p.adjust) >0.05 ) {
    
    
    if (use_semData) {
      x <- pairwise_termsim(filter_go_CC,semData = hsGO_CC,method = "Wang")
    }else{
      x <- pairwise_termsim(filter_go_CC)
    }
    
    
    tryCatch(expr = {
      
      p <- treeplot(x,showCategory = length(filter_go_CC@result$ID),nCluster = nCluster)+coord_cartesian(clip = "off")
      pdf(file =  paste0(outputpath,"/tree_plot_005-01","_nCluster_",nCluster,dir_end,"/tree_plot_CC_",name,".pdf"),width = 12,height = 5+(length(filter_go_CC@result$ID)*0.2))
      print(p)
      dev.off()
      
      
    },
    error = function(e){
      e
    })
    
    
    
  }
  
  
  
  
  
}


