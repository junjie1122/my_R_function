
# version 20221009

# @description: 批量做Reactome富集分析

# @parameters: genelist genelist 要有names
# @parameters: keyType 输入的基因类型 
# @parameters: organism 物种选择 one of "human", "rat", "mouse"
# @parameters: showCategory plot展示的富集条数数量, 为NULL 展示所有adj.p小于0.05的

# @parameters: outputpath 输出路径



run_reactome <- function(genelist,xlsxfilename = "",keyType,organism,showCategory = NULL,outputpath)
{ 

  library(ReactomePA)
  library(openxlsx)
  library(clusterProfiler)
  library(stringr)
  
  
  if (!file.exists(outputpath)) {
    dir.create(outputpath,recursive = T)
  }
  
  
  if (!dir.exists(paste0(outputpath,"/reactome_enrich_RData/"))) {
    dir.create(paste0(outputpath,"/reactome_enrich_RData/"),recursive = T)
  }
  
  dir.create(paste0(outputpath,"/dotplot"))
  
  
  wb <- createWorkbook()
  wb2 <- createWorkbook()
  
  # reactome_dotplot_list <- list()


  for (i in 1:length(genelist)) {
    
    
    
    name <- names(genelist[i])
    
    print(c(i,name))
    
    gene <- as.data.frame(genelist[i])
    gene <- gene[,1]

    # 排除 vlincRNA
    if(organism == "human" & keyType == "ENSEMBL"){
      gene <- grep("^ENSG",gene,value = T)
    }
    
    
    
    OrgDb = switch (organism,
                    'human' = "org.Hs.eg.db",
                    'mouse' = "org.Mm.eg.db",
                    "Rat" = "org.Rn.eg.db"
    )
    
    
   
    # 基因ID转换
    if (keyType != "ENTREZID") {
      
      # 转换0个gene_id会报错 用tryCatch 报错时跳过 "ENTREZ_id is NULL" 赋值给 possibleError
      possibleError <- tryCatch(expr = { entrez <- bitr(gene,fromType = keyType,toType = "ENTREZID",OrgDb = OrgDb)},
                                error = function(e){
                                  "ENTREZ_id is NULL"
                                })
      
      if (length(possibleError) == 1) {
        addWorksheet(wb,sheetName = name)
        addWorksheet(wb2,sheetName = name)
        print(str_c(name,":ENTREZ_id is NULL"))
        next
      }
      entrez_gene <- entrez$ENTREZID
    }
    # keyType == "ENTREZID" 不用转换
    else{
      entrez_gene <- gene
    }
    
    
    if (keyType != "SYMBOL") {
      enrich_reactome <-enrichPathway(gene = entrez_gene,organism = organism,pvalueCutoff = 1,qvalueCutoff = 1,readable = T)
    }else{
      enrich_reactome <-enrichPathway(gene = entrez_gene,organism = organism,pvalueCutoff = 1,qvalueCutoff = 1)
    }
    
    # 保存富集分析原始结果
    save(enrich_reactome,file = paste0(outputpath,"/reactome_enrich_RData/",name,"_enrich_reactome.RData") )
    
    
    if (is.null(enrich_reactome)) {
      addWorksheet(wb,sheetName = name)
      addWorksheet(wb2,sheetName = name)
      next
    }
    
    
    enrich_reactome_result <- enrich_reactome@result
  
    
    addWorksheet(wb,sheetName = name)
    writeData(wb,i,enrich_reactome_result,colNames = T,rowNames = T,)
    

    filter_enrich_reactome <- enrich_reactome[enrich_reactome@result$p.adjust<0.05,asis=T] %>%as.data.frame()
    
    
    addWorksheet(wb2,sheetName = name)
    writeData(wb2,i,filter_enrich_reactome,colNames = T,rowNames = T)
    
    
    filter_enrich_reactome <- enrich_reactome[enrich_reactome@result$p.adjust<0.05,asis=T]

    if (length(filter_enrich_reactome@result$ID) >= 1 ) {
      
      

      showCategory <- ifelse(is.null(showCategory),length(filter_enrich_reactome@result$ID),showCategory)
      
      print(showCategory)

      p <- dotplot(filter_enrich_reactome,showCategory = showCategory,x = "Count") +scale_size(range = c(4,8)) +scale_y_discrete(labels= function(x) str_wrap(x,width = 45))

      pdf(file =  paste0(outputpath,"/dotplot/Reactome_dotplot_",name,".pdf"),width = 9,height = 5+showCategory*0.2)
      print(p)
      dev.off()

    }
    
  }

  # save(reactome_dotplot_list,file = paste0(outputpath,"/reactome_dotplot_list.RData"))
  saveWorkbook(wb,file = paste0(outputpath,"/reactome_enrichment_result_all_",xlsxfilename,".xlsx"),overwrite = T)
  saveWorkbook(wb2,file = paste0(outputpath,"/reactome_enrichment_result_padj005_",xlsxfilename,".xlsx"),overwrite = T)
}
