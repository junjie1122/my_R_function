

# @description: 批量做KEGG富集分析

# @parameters: genelist genelist 要有names
# @parameters: keyType 输入的基因类型
# @parameters: organism 物种选择 默认人 "hsa",	supported organism listed in 'https://www.genome.jp/kegg/catalog/org_list.html'
# @parameters: showCategory plot展示的富集条数数量, 为NULL 展示所有adj.p小于0.05的
# @parameters: outputpath 输出路径




run_kegg <- function(genelist,xlsxfilename = "",keyType ="SYMBOL",organism = "hsa",showCategory = NULL ,outputpath,showCategorynum = 15){ 
  
  
  R.utils::setOption( "clusterProfiler.download.method",'auto' )
  
  library(openxlsx)
  library(clusterProfiler)
  library(stringr)
  
  
  if (!file.exists(outputpath)) {
    dir.create(outputpath,recursive = T)
  }
  
  if (!dir.exists(paste0(outputpath,"/kegg_enrich_RData/"))) {
    dir.create(paste0(outputpath,"/kegg_enrich_RData/"),recursive = T)
  }
  
  dir.create(paste0(outputpath,"/dotplot"),recursive = T)
  
  
  wb <- createWorkbook()
  wb2 <- createWorkbook()
  

  
  for (i in 1:length(genelist)) {
    
    
    name <- names(genelist[i])
    
    print(c(i,name))
    
    gene <- as.data.frame(genelist[i])
    gene <- gene[,1]
    
    # 排除 vlincRNA
    if(organism == "hsa" & keyType == "ENSEMBL"){
      gene <- grep("^ENSG",gene,value = T)
    }
    
    
    
    if (keyType != "ENTREZID") {

	OrgDb = switch (organism,
                    'hsa' = "org.Hs.eg.db"
                     )





    
      ## 转换0个gene_id会报错 用tryCatch 报错时跳过 "ENTREZ_id is NULL" 赋值给 possibleError
      possibleError <- tryCatch(expr = { entrez <- bitr(gene,fromType = keyType,toType = "ENTREZID",OrgDb = "org.Hs.eg.db")},
                                error = function(e){
                                  "ENTREZ_id is NULL"
                                })
      
      if (all(possibleError %in%  "ENTREZ_id is NULL") ) {
        addWorksheet(wb,sheetName = name)
        addWorksheet(wb2,sheetName = name)
        print(str_c(name,":ENTREZ_id is NULL"))
        next
      }
    
      entrez_gene <- entrez$ENTREZID  
    }
    
    else{entrez_gene <- gene}
    
    
    
    
    enrich_kegg <- enrichKEGG(
      gene = entrez_gene,
      keyType = "kegg",
      organism  = organism,
      pvalueCutoff  = 1,
      pAdjustMethod  = "BH",
      qvalueCutoff  = 1,
    )
    
    
    
    if(organism == "hsa"){
      enrich_kegg <- setReadable(enrich_kegg,OrgDb = "org.Hs.eg.db",keyType = "ENTREZID")
    }

    
    # 保存富集分析原始结果
    save(enrich_kegg,file = paste0(outputpath,"/kegg_enrich_RData/",name,"_kegg_enrich.RData") )
    
    
    
    if (is.null(enrich_kegg)) {
      addWorksheet(wb,sheetName = name)
      addWorksheet(wb2,sheetName = name)
      
      next
    }
    
    
    kegg_result <- enrich_kegg@result
    
    
    addWorksheet(wb,sheetName = name)
    writeData(wb,i,kegg_result,colNames = T,rowNames = T,)
    
    
    
    filter_enrich_kegg <- enrich_kegg[enrich_kegg@result$p.adjust<0.05,asis=T] %>%as.data.frame()
    addWorksheet(wb2,sheetName = name)
    writeData(wb2,i,filter_enrich_kegg,colNames = T,rowNames = T)
    

    filter_enrich_kegg <- enrich_kegg[enrich_kegg@result$p.adjust<0.05,asis=T]
    
    if (length(filter_enrich_kegg@result$ID) >= 1 ) {
      
      
      showCategory <- ifelse(is.null(showCategory),length(filter_enrich_kegg@result$ID),showCategory)
      
      p <- dotplot(filter_enrich_kegg,showCategory = showCategory,x = "Count") +scale_size(range = c(4,8)) +scale_y_discrete(labels= function(x) str_wrap(x,width = 45))
      
      pdf(file =  paste0(outputpath,"/dotplot/kegg_dotplot_",name,".pdf"),width = 9,height = 5+showCategory*0.2 )
      print(p)
      dev.off()
    }
    
    
  }
  
  saveWorkbook(wb,file = paste0(outputpath,"/kegg_enrichment_result_all_",xlsxfilename,".xlsx"),overwrite = T)
  saveWorkbook(wb2,file = paste0(outputpath,"/kegg_enrichment_result_padj005_",xlsxfilename,".xlsx"),overwrite = T)
}
#####################################################################################################################  



