

# @description: 把list写入到xlsx的不用sheet

# @parameters: input 输入的list
# @parameters: output_file 输出文件路径
# @parameters: names 是否提供额外的sheet names，默认为NULL，使用list的names

# @return: 没有return

# @example: write_list2xslx(input = gene_list,output_file = "output_test/test_write_list2xlsx.xlsx")


write_list2xslx <- function(input,output_file,names = NULL,set_col_name = T,set_row_name = T){
  
  library(openxlsx)
  # input = gene_list
  # output = "output_test"

  path <- output_file %>% str_extract("(.*/){1}")
  
  if (!dir.exists(path)) {
    dir.create(path,recursive = T)
  }
  
  
  
  wb <- createWorkbook()
  
  if (is.null(names)) {
    names <- names(input)
  }
  
  for (i in 1:length(input)) {
    
    addWorksheet(wb,sheetName = names[i])

    writeData(wb,names[i],input[[i]],colNames = set_col_name,rowNames = set_row_name)
  }
  saveWorkbook(wb,file = output_file,overwrite = T)
}



