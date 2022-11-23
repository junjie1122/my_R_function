


# 获取当前run脚本的路径  wgcna_function脚本跟run脚本放一个文件夹
run_WGCNA_Rscript_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
source(paste0(run_WGCNA_Rscript_dir,"/wgcna_function.R"))


data <- data_normalization(data = exp_data,method = method) 
dim(data)


data <- gene_filter(data = data,filter_gene_method = filter_gene_method,keep_gene_num = keep_gene_num,filter_threshold = filter_threshold,out_path = out_path)
dim(data)




if (is.null(sampletree_cut_height)) {
  sample_clust(data = data,trait = trait,width = 9,height = 6,out_path = out_path)
}else{
  sample_clust(data = data,trait = trait,width = 9,height = 6,abHeight = sampletree_cut_height,out_path = out_path)
  remove_outline_sample(data = data,trait = trait,cutHeight = sampletree_cut_height,width = 9,height = 6)
}



power <- power_choose(data = data,cor_type = cor_type,net_type = net_type,max_blockSize = max_blockSize,width = 9,height = 6,out_path = out_path)
power



net_build_out <- net_build(data = data,power = power,max_blockSize = max_blockSize,cor_type = cor_type,
                            net_type = net_type,TOM_type = TOM_type,mergeCutHeight = mergeCutHeight,out_path = out_path) 



dendrograms_plot(net = net,moduleColors = moduleColors,width = 9,height = 6,out_path = out_path)


mes_heatmap_plot(data = data,MEs = MEs,moduleColors = moduleColors,
                 width = 9,height = 7,out_path = out_path)



mes_trait_relation(data = data,MEs = MEs,trait = trait,width = 6,height = 9,out_path = out_path)



significent_mod <- get_significent_mod(moduleTraitPvalue = moduleTraitPvalue)
significent_mod




list_module_gene <- save_significent_mod_gene(data = data,MEs = MEs,moduleColors = moduleColors,gene_anno = gene_anno,merge_by = gene_colnames,group_by =gene_type_colname,
                                              significent_mod = significent_mod,out_path = out_path)



get_MM_out <- get_MM(data = data,MEs = MEs,moduleColors = moduleColors,out_path = out_path)




get_gs_mm_plot(data = data,MEs = MEs,trait = trait,moduleColors = moduleColors,width = 9,height = 6,
               significent_mod = significent_mod,moduleTraitPvalue = moduleTraitPvalue,out_path = out_path)



TOM <- get_tom(net = net)

get_tom_heatmap(net = net,TOM = TOM,nselect = 1200,out_path = out_path,width = 9)


get_mes_cor_heatmap(net = net,out_path = out_path,width = 9,height = 6)



module_gene_to_cyto(TOM = TOM,moduleColors = moduleColors,
                    significent_mod = significent_mod,out_path = out_path)


merge_GS_MM_table(out_path = out_path)


info_summary()



adj_matrix <- get_adj_matrix(data = data,power = power,cor_type = cor_type,net_type = net_type,out_path = out_path)



connect <- get_connectivity(adj_matrix = adj_matrix,gene_anno = gene_anno,merge_by = gene_colnames,
                            moduleColors = moduleColors,out_path = out_path)



get_hub_gene_gsmm(data = data,MEs = MEs,trait = trait,moduleColors = moduleColors,significent_mod = significent_mod,
                  moduleTraitPvalue = moduleTraitPvalue,GS_threshold = GS_threshold,MM_threshold = MM_threshold,out_path = out_path)



rmarkdown::render('E:/R_function/WGCNA/wgcna_info_summary.qmd', output_file = paste0(out_path,"/WGCNA_info_summary") )





# GO reactome 富集分析
source("./R/run_go.R")
source("./R/run_reactome.R")

run_go(genelist = list_module_gene,outputpath = paste0(out_path,"GO/"),keyType = "ENSEMBL")     
run_reactome(genelist = list_module_gene,outputpath = paste0(out_path,"reactome/") , keyType = "ENSEMBL")




# get_hub_gene_gs_connect(connect = connect,data = data,trait = trait,moduleColors = moduleColors,
#                         significent_mod = significent_mod,moduleTraitPvalue = moduleTraitPvalue,ntop = 10)


# hub_gene_to_cyto(data = data,TOM = TOM,hubgene_list = list_gs_mm_hubgene,get_hubgene_method = "GS_MM",out_path = out_path)


# hub_gene_to_cyto(data = data,TOM = TOM,hubgene_list = list_gs_conn_hubgene,get_hubgene_method = "GS_connect",out_path = out_path)


# 富集分析
if (F) {

    # GO 富集分析
    source("./R/run_go.R")
    run_go(genelist = list_module_gene,outputpath = paste0(out_path,"GO/module_gene"),keyType = "SYMBOL")
    run_go(genelist = list_gs_conn_hubgene,outputpath = paste0(out_path,"GO/gs_connect_hubgene"),keyType = "SYMBOL")
    run_go(genelist = list_gs_mm_hubgene,outputpath = paste0(out_path,"GO/gs_mm_hubgene"),keyType = "SYMBOL")



    # KEGG富集分许
    source("./R/run_kegg.R")
    run_kegg(genelist = list_module_gene,outputpath = paste0(out_path,"KEGG/module_gene"),formgeneType = "SYMBOL")
    run_kegg(genelist = list_gs_conn_hubgene,outputpath = paste0(out_path,"KEGG/gs_connect_hubgene_gene"),formgeneType = "SYMBOL")
    run_kegg(genelist = list_gs_mm_hubgene,outputpath = paste0(out_path,"KEGG/gs_mm_hubgene_gene"),formgeneType = "SYMBOL")

}
