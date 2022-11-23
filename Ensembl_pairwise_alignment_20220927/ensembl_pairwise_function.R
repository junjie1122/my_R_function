

# @description: ENSEMBL数据库 pairwise alignment 找同源 

# @parameters: human_reigon eg: "8:130248522-130248608:1"
# @parameters: output_path 输出路径 

# @return: 绑定 align_fasta_all  align_reigon_all  到全局环境 

# @example: get_pairwise_alignment_result(human_reigon = "8:130248522-130248608:1",output_path = "./chr8_130248522-130248608/")


# 运行函数
get_pairwise_alignment_result <- function(human_reigon, output_path){
  
  
  library(httr)
  library(jsonlite)
  library(xml2)
  library(tidyverse)
  library(Biostrings)
  library(openxlsx)
  library(plyr)
  
  
  if (!dir.exists(output_path)) {
    dir.create(output_path,recursive = T)
  }
  
  
  
  name_table <<- read.xlsx("E:/R_function/Ensembl_pairwise_alignment_20220927/ensembl_pairwise_alignment_species_20220710.xlsx")
  
  
  align_fasta_all <- data.frame()
  align_reigon_all <- data.frame()
  error_message_all <- data.frame()
  
  for (i in 1:nrow(name_table)) {
    
    
    alig_dat <- pairwise_alignment_function(human_reigon = human_reigon,scientific_name = name_table[i,2],name_table = name_table)
    
    
    if (length(alig_dat) == 2) {
      
      align_fasta_all <- rbind.fill(align_fasta_all,alig_dat[[1]])
      align_reigon_all <- rbind.fill(align_reigon_all,alig_dat[[2]])
      assign("align_reigon_all",value = align_reigon_all,envir = .GlobalEnv)
      
      
      
    }else if(length(alig_dat) > 3){
      
      alig_result_g3 <- map_dfr(alig_dat,as.data.frame)
      alig_result_g3 <-  alig_result_g3 %>% select(.,species,seq_region,start,end,strand,seq, description) %>% 
        mutate(species = if_else(species != "homo_sapiens",name_table[i,1],"homo_sapiens"))
      
      
      
      write.xlsx(alig_result_g3,file = paste0(output_path,"reigon数多于3的比对结果_",name_table[i,1],".xlsx"))
      
      
    }else{
      
      error_message_all <- rbind.fill(error_message_all,as.data.frame(alig_dat[[1]]))
    }
    
    
  }
  
  
  
  
  align_fasta_all <- fix_fa(fa = align_fasta_all)
  assign("align_fasta_all",value = align_fasta_all,envir = .GlobalEnv)
  
  
  save_seqs <- c(seq(3,nrow(align_fasta_all),4),seq(4,nrow(align_fasta_all),4)) %>% sort()
  align_fasta_all_remove_human <- align_fasta_all[save_seqs,,drop = F]
  
  mult_reigon <- align_reigon_all[!is.na(align_reigon_all$seq_region_2),]
  
  
  
  # 得到mult_reigon_species_full_seq
  get_mult_reigon_species_full_seq(mult_reigon = mult_reigon,output_path = output_path,name_table = name_table)
  
  
  
  
  write.table(align_fasta_all,file = paste0(output_path,"pairwise_alignment_fasta.txt"),quote = F,row.names = F,col.names = F)
  
  write.table(align_fasta_all_remove_human,file = paste0(output_path,"pairwise_alignment_fasta_remove_human.txt"),quote = F,row.names = F,col.names = F)
  
  write.xlsx(align_reigon_all,file = paste0(output_path,"pairwise_alignment_reigon_info.xlsx"))
  
  write.table(error_message_all,file = paste0(output_path,"pairwise_alignment_error_info.txt"),quote = F,row.names = F,col.names = F)
  
  write.xlsx(mult_reigon ,file = paste0(output_path,"pairwise_alignment_mult_reigon_info.xlsx"))
  

  
}




reigon_to_seq <- function(reigon,scientific_name,name_table){
  
  df <- data.frame()
  server <- "https://rest.ensembl.org"
  ext <- paste0("/sequence/region/",scientific_name,"/",reigon,"?")
  
  
  # tryCatch(r <- RETRY(paste(server, ext, sep = ""), content_type("text/plain"),times = 6),
  #          error = function(e){
  #            r <- RETRY(paste(server, ext, sep = ""), content_type("text/plain"),times = 6)
  #          })
  # 
  # r <- GET(paste(server, ext, sep = ""), content_type("text/plain"))
  
  
  r <- RETRY(verb = "GET",paste(server, ext, sep = ""), content_type("text/plain"),times = 6)
  

  
  if (r$status_code == 400) {
    
    print(paste0(scientific_name,": fail to get seq"))
  }else{
    
    
    stop_for_status(r)
    
    #print(content(r))
    
    df[1,1] <- paste0(">",scientific_name)
    df[2,1] <- content(r)
    
    df[1,1] <- paste0(">",name_table[name_table$scientific.name == str_remove(df[1,1],">"),1] )
    
    return(df)
    
  }
  
  # return(list(Scientific_name,reigon,content(r)) )
  # return(content(r) )
  
}



pairwise_alignment_function <- function(human_reigon,scientific_name,name_table){
  
  
  human_reigon_start <- human_reigon %>% str_extract("(:.*:)") %>% str_remove_all("-.*") %>% str_remove(":") %>% as.numeric()
  human_reigon_end <- human_reigon %>% str_extract("(:.*:)") %>% str_remove_all(".*-") %>% str_remove(":") %>% as.numeric()
  
  human_seq_length <- human_reigon_end - human_reigon_start +1
  
  strand <- human_reigon %>% str_remove_all(".*:") %>% as.numeric()
  
  output_align_fasta <- data.frame()
  output_align_reigon <- data.frame()
  
  # reigon_col_sort <- c("common_name","scientific_name","seq_region_1","start_1","end_1","strand_1","seq_region_2",
  #                                 "start_2","end_2","strand_2","alig_seq_2","seq_region_3","start_3","end_3","strand_3","alig_seq_3")
  
  
  output_mult_alig_fasta <-  data.frame()
  
  
  server <- "https://rest.ensembl.org"
  ext <- paste0("/alignment/region/Homo_sapiens/",human_reigon,"?species_set=Homo_sapiens;species_set=",scientific_name,";method=LASTZ_NET")
  
  # r <- GET(paste(server, ext, sep = ""), content_type("application/json"))
  
  
  
  tryCatch(expr = {r <- RETRY(verb = "GET",paste(server, ext, sep = ""),content_type("application/json"),times = 6)
                  },
           error = function(e){
             print("try again")
             tryCatch(r <- RETRY(verb = "GET",paste(server, ext, sep = ""),content_type("application/json"),times = 6),
                      error = function(e){
                        print("try again again ")
                        tryCatch(r <- RETRY(verb = "GET",paste(server, ext, sep = ""),content_type("application/json"),times = 6),
                                 error = function(e){
                                   print("try again again again")
                                   r <- RETRY(verb = "GET",paste(server, ext, sep = ""),content_type("application/json"),times = 6)
                                   
                                 })
                       
                        
                      })
             
           })
  
  # r <- RETRY(verb = "GET",paste(server, ext, sep = ""),content_type("application/json"),times = 6)
  
  
  
  if (names(head(fromJSON(toJSON(content(r)))) )[1] == "error") {
    
    print(paste0("occur error: ",name_table[name_table$Scientific.name == scientific_name,"Common.name"],"   ",content(r)$error) )
    
    
    error_message <- paste0(name_table[name_table$Scientific.name == scientific_name,"Common.name"],"\t",content(r)$error)
    
    return(error_message)
    
  }else{
    
    stop_for_status(r)
    
    a <- head(fromJSON(toJSON(content(r))))
    
    result <- a$alignments 
    # print(result)
    
    
    #当只匹配一个区域的时候
    if (length(result) == 1) {
      
      
      output_align_reigon[1,"common_name"] <- name_table[name_table$Scientific.name == scientific_name,"Common.name"]
      output_align_reigon[1,"scientific_name"] <- scientific_name
      
      output_align_reigon[1,"seq_region_1"] <- result[[1]][2,"seq_region"] %>% unlist()
      
      output_align_reigon[1,"start_1"] <- result[[1]][2,"start"]  %>% unlist()
      output_align_reigon[1,"end_1"] <- result[[1]][2,"end"] %>% unlist()
      output_align_reigon[1,"strand_1"] <- result[[1]][2,"strand"] %>% unlist()
      
      
      #当找到的homologous sequence没有.....的时候 human_seq是完整的 
      if (result[[1]][1,"start"] ==  human_reigon_start & result[[1]][1,"end"] ==  human_reigon_end ) {
        
        human_seq <- result[[1]][1,"seq"] %>% unlist()
        ailgn_seq <- result[[1]][2,"seq"] %>% unlist()
        
        output_align_reigon[1,"alig_seq_1"] <- ailgn_seq
        
        
        #当找到的homologous sequence有..... human_seq只能显示部分,...对应的human_seq是没有的，start end的位置也是去除...后的位置
        #end不相等时 dot在左边
      }else if (result[[1]][1,"end"] != human_reigon_end ){
        
        human_seq <- result[[1]][1,"seq"] %>% unlist()
        ailgn_seq <- result[[1]][2,"seq"] %>% unlist()
        
        human_seq_wild <- reigon_to_seq(reigon = human_reigon,scientific_name = "homo_sapiens",name_table = name_table )
        human_seq_wild <- human_seq_wild[2,1]
        #human_seq_wild_complement <- complement(DNAString(human_seq_wild[2,1])) %>% toString()
        
        
        human_seq_cut <- human_seq %>% str_remove_all("-")  #防止有- locate报错
        
        start_site <- str_locate(human_seq_wild, human_seq_cut) %>% .[,"start"]
        
        # print(start_site)
        
        
        dot <- rep(".", as.numeric(start_site - 1) ) %>% str_c(collapse = "")
        ailgn_seq_add_dot <- paste0(dot,ailgn_seq)
        
        
        human_seq <- human_seq_wild
        ailgn_seq <- ailgn_seq_add_dot
        output_align_reigon[1,"alig_seq_1"] <- ailgn_seq
        
        #start不相等时 dot在右边
      }else if (result[[1]][1,"start"] !=  human_reigon_start) {
        
        human_seq <- result[[1]][1,"seq"] %>% unlist()
        ailgn_seq <- result[[1]][2,"seq"] %>% unlist()
        
        human_seq_wild <- reigon_to_seq(reigon = human_reigon,scientific_name = "homo_sapiens",name_table = name_table)
        human_seq_wild <- human_seq_wild[2,1]
        #human_seq_wild_complement <- complement(DNAString(human_seq_wild[2,1])) %>% toString()
        
        human_seq_cut <- human_seq %>% str_remove_all("-")
        
        end_site <- str_locate(human_seq_wild, human_seq_cut) %>% .[,"end"]
        
        
        dot <- rep(".",nchar(human_seq_wild)-end_site ) %>% str_c(collapse = "")
        ailgn_seq_add_dot <- paste0(ailgn_seq,dot)
        
        human_seq <- human_seq_wild
        ailgn_seq <- ailgn_seq_add_dot
        output_align_reigon[1,"alig_seq_1"] <- ailgn_seq
        
        #  start end 都不相等时 dot在左右边
      }else{
        
        human_seq <- result[[1]][1,"seq"] %>% unlist()
        ailgn_seq <- result[[1]][2,"seq"] %>% unlist()
        
        human_seq_wild <- reigon_to_seq(reigon = human_reigon,scientific_name = "homo_sapiens",name_table = name_table )
        human_seq_wild <- human_seq_wild[2,1]
        
        human_seq_cut <- human_seq %>% str_remove_all("-")
        
        start_site <- str_locate(human_seq_wild, human_seq_cut) %>% .[,"start"]
        print(start_site)
        
        dot_s <- rep(".", as.numeric(start_site-1) ) %>% str_c(collapse = "")
        
        end_site <- str_locate(human_seq_wild, human_seq_cut) %>% .[,"end"]
        dot_e <- rep(".",nchar(human_seq_wild)-end_site ) %>% str_c(collapse = "")
        
        ailgn_seq_add_dot <- paste0(dot_s,ailgn_seq,dot_e)
        
        human_seq <- human_seq_wild
        ailgn_seq <- ailgn_seq_add_dot
        output_align_reigon[1,"alig_seq_1"] <- ailgn_seq
        
      }
      
      
      
      
      
    }else if(length(result) == 2){
      
      
      if (strand == -1) {
        
        output_align_reigon[1,"common_name"] <- name_table[name_table$Scientific.name == scientific_name,"Common.name"]
        output_align_reigon[1,"scientific_name"] <- scientific_name
        
        output_align_reigon[1,"seq_region_1"] <- result[[2]][2,"seq_region"] %>% unlist()
        output_align_reigon[1,"start_1"] <- result[[2]][2,"start"]  %>% unlist()
        output_align_reigon[1,"end_1"] <- result[[2]][2,"end"] %>% unlist()
        output_align_reigon[1,"alig_seq_1"] <- result[[2]][2,"seq"] %>% unlist()
        output_align_reigon[1,"strand_1"] <- result[[2]][2,"strand"] %>% unlist()
        
        
        
        output_align_reigon[1,"seq_region_2"] <- result[[1]][2,"seq_region"] %>% unlist()
        output_align_reigon[1,"start_2"] <- result[[1]][2,"start"]  %>% unlist()
        output_align_reigon[1,"end_2"] <- result[[1]][2,"end"] %>% unlist()
        output_align_reigon[1,"alig_seq_2"] <- result[[1]][2,"seq"] %>% unlist()
        output_align_reigon[1,"strand_2"] <- result[[1]][2,"strand"] %>% unlist()
        
      }else{
        
        output_align_reigon[1,"common_name"] <- name_table[name_table$Scientific.name == scientific_name,"Common.name"]
        output_align_reigon[1,"scientific_name"] <- scientific_name
        
        output_align_reigon[1,"seq_region_1"] <- result[[1]][2,"seq_region"] %>% unlist()
        output_align_reigon[1,"start_1"] <- result[[1]][2,"start"]  %>% unlist()
        output_align_reigon[1,"end_1"] <- result[[1]][2,"end"] %>% unlist()
        output_align_reigon[1,"alig_seq_1"] <- result[[1]][2,"seq"] %>% unlist()
        output_align_reigon[1,"strand_1"] <- result[[1]][2,"strand"] %>% unlist()
        
        
        
        output_align_reigon[1,"seq_region_2"] <- result[[2]][2,"seq_region"] %>% unlist()
        output_align_reigon[1,"start_2"] <- result[[2]][2,"start"]  %>% unlist()
        output_align_reigon[1,"end_2"] <- result[[2]][2,"end"] %>% unlist()
        output_align_reigon[1,"alig_seq_2"] <- result[[2]][2,"seq"] %>% unlist()
        output_align_reigon[1,"strand_2"] <- result[[2]][2,"strand"] %>% unlist()
        
      }
      
      
      
      # 最后一个条件说明 human序列没有被裁剪，不存在....
      if (result[[1]][1,"start"] ==  human_reigon_start & result[[2]][1,"end"] ==  human_reigon_end & (as.numeric(result[[1]][1,"end"])+1 == result[[2]][1,"start"])  ) {
        
        
        human_seq <- if_else(strand == -1,paste0(unlist(result[[2]][1,"seq"]),unlist(result[[1]][1,"seq"]) ), 
                             paste0(unlist(result[[1]][1,"seq"]),unlist(result[[2]][1,"seq"]) ))
        
        
        ailgn_seq <- if_else(strand == -1,paste0(unlist(result[[2]][2,"seq"]),unlist(result[[1]][2,"seq"]) ), 
                             paste0(unlist(result[[1]][2,"seq"]),unlist(result[[2]][2,"seq"]) ))
        
        
        # mult_alig_human_seq <- ailgn_seq
        # output_mult_alig_fasta[1,1] <- ">human"
        # output_mult_alig_fasta[1,1] <- human_seq
        # output_mult_alig_fasta[1,1] <- paste0(">",name_table[name_table$Scientific.name == scientific_name,"Common.name"] )
        # output_mult_alig_fasta[1,1] <- mult_alig_human_seq
        
        
        
        
        #若 存在.... 也就是 网站里结果有.... 但 API代码得到的结果没...需要手动补上....
      }else {
        
        # print("error: multiple reigon == 2,exit ....")
        # stop()
        
      
        if (strand == -1) {
          
          human_seq_wild <- reigon_to_seq(reigon = human_reigon,scientific_name = "homo_sapiens",name_table = name_table )
          human_seq_wild <- human_seq_wild[2,1]
          
          alig_seq_1 <- unlist(result[[2]][2,"seq"])
          alig_seq_2 <- unlist(result[[1]][2,"seq"])
          
          
          start_site_1 <- str_locate(human_seq_wild,unlist(result[[2]][1,"seq"]) %>% str_remove_all("-")    ) %>% .[,"start"]
          end_site_1 <- str_locate(human_seq_wild,unlist(result[[2]][1,"seq"]) %>% str_remove_all("-")   ) %>% .[,"end"]
          
          start_site_2 <- str_locate(human_seq_wild,unlist(result[[1]][1,"seq"])  %>% str_remove_all("-")  ) %>% .[,"start"]
          end_site_2 <- str_locate(human_seq_wild,unlist(result[[1]][1,"seq"]) %>% str_remove_all("-")  ) %>% .[,"end"]
          
          
          dot_1 <- rep(".",start_site_1 - 1) %>% str_c(collapse = "")
          dot_2 <- rep(".",start_site_2 - end_site_1 - 1) %>% str_c(collapse = "")
          dot_3 <- rep(".",human_seq_length-end_site_2) %>% str_c(collapse = "")
          
          
          ailgn_seq_add_dot <- paste0(dot_1,alig_seq_1,dot_2,alig_seq_2,dot_3)
          
          # 想要的
          human_seq <- human_seq_wild
          ailgn_seq <- ailgn_seq_add_dot
          
        }else{
          
          human_seq_wild <- reigon_to_seq(reigon = human_reigon,scientific_name = "homo_sapiens",name_table = name_table )
          human_seq_wild <- human_seq_wild[2,1]
          
          alig_seq_1 <- unlist(result[[1]][2,"seq"])
          alig_seq_2 <- unlist(result[[2]][2,"seq"])
          
          start_site_1 <- str_locate(human_seq_wild,unlist(result[[1]][1,"seq"])%>% str_remove_all("-")  ) %>% .[,"start"]
          end_site_1 <- str_locate(human_seq_wild,unlist(result[[1]][1,"seq"])%>% str_remove_all("-")  ) %>% .[,"end"]
          
          start_site_2 <- str_locate(human_seq_wild,unlist(result[[2]][1,"seq"])%>% str_remove_all("-")  ) %>% .[,"start"]
          end_site_2 <- str_locate(human_seq_wild,unlist(result[[2]][1,"seq"])%>% str_remove_all("-")  ) %>% .[,"end"]
          
          
          dot_1 <- rep(".",start_site_1 - 1) %>% str_c(collapse = "")
          dot_2 <- rep(".",start_site_2 - end_site_1 - 1) %>% str_c(collapse = "")
          dot_3 <- rep(".",human_seq_length-end_site_2) %>% str_c(collapse = "")
          
          
          ailgn_seq_add_dot <- paste0(dot_1,alig_seq_1,dot_2,alig_seq_2,dot_3)
          
          human_seq <- human_seq_wild
          ailgn_seq <- ailgn_seq_add_dot
          
        }

        
        
        
      }
      
      
      
      
      #同上 
    }else if (length(result) == 3) {
      
      
      if (strand == -1) {
        
        output_align_reigon[1,"common_name"] <- name_table[name_table$Scientific.name == scientific_name,"Common.name"]
        output_align_reigon[1,"scientific_name"] <- scientific_name
        
        output_align_reigon[1,"seq_region_1"] <- result[[3]][2,"seq_region"] %>% unlist()
        output_align_reigon[1,"start_1"] <- result[[3]][2,"start"]  %>% unlist()
        output_align_reigon[1,"end_1"] <- result[[3]][2,"end"] %>% unlist()
        output_align_reigon[1,"strand_1"] <- result[[3]][2,"strand"] %>% unlist()
        output_align_reigon[1,"alig_seq_1"] <- result[[3]][2,"seq"] %>% unlist()
        
        
        output_align_reigon[1,"seq_region_2"] <- result[[2]][2,"seq_region"] %>% unlist()
        output_align_reigon[1,"start_2"] <- result[[2]][2,"start"]  %>% unlist()
        output_align_reigon[1,"end_2"] <- result[[2]][2,"end"] %>% unlist()
        output_align_reigon[1,"strand_2"] <- result[[2]][2,"strand"] %>% unlist()
        output_align_reigon[1,"alig_seq_2"] <- result[[2]][2,"seq"] %>% unlist()
        
        
        output_align_reigon[1,"seq_region_3"] <- result[[1]][2,"seq_region"] %>% unlist()
        output_align_reigon[1,"start_3"] <- result[[1]][2,"start"]  %>% unlist()
        output_align_reigon[1,"end_3"] <- result[[1]][2,"end"] %>% unlist()
        output_align_reigon[1,"strand_3"] <- result[[1]][2,"strand"] %>% unlist()
        output_align_reigon[1,"alig_seq_3"] <- result[[1]][2,"seq"] %>% unlist()
        
        
        
      }else{
        
        output_align_reigon[1,"common_name"] <- name_table[name_table$Scientific.name == scientific_name,"Common.name"]
        output_align_reigon[1,"scientific_name"] <- scientific_name
        
        output_align_reigon[1,"seq_region_1"] <- result[[1]][2,"seq_region"] %>% unlist()
        output_align_reigon[1,"start_1"] <- result[[1]][2,"start"]  %>% unlist()
        output_align_reigon[1,"end_1"] <- result[[1]][2,"end"] %>% unlist()
        output_align_reigon[1,"strand_1"] <- result[[1]][2,"strand"] %>% unlist()
        output_align_reigon[1,"alig_seq_1"] <- result[[1]][2,"seq"] %>% unlist()
        
        
        output_align_reigon[1,"seq_region_2"] <- result[[2]][2,"seq_region"] %>% unlist()
        output_align_reigon[1,"start_2"] <- result[[2]][2,"start"]  %>% unlist()
        output_align_reigon[1,"end_2"] <- result[[2]][2,"end"] %>% unlist()
        output_align_reigon[1,"strand"] <- result[[2]][2,"strand"] %>% unlist()
        output_align_reigon[1,"alig_seq_2"] <- result[[2]][2,"seq"] %>% unlist()
        
        output_align_reigon[1,"seq_region_3"] <- result[[3]][2,"seq_region"] %>% unlist()
        output_align_reigon[1,"start_3"] <- result[[3]][2,"start"]  %>% unlist()
        output_align_reigon[1,"end_3"] <- result[[3]][2,"end"] %>% unlist()
        output_align_reigon[1,"strand"] <- result[[3]][2,"strand"] %>% unlist()
        output_align_reigon[1,"alig_seq_3"] <- result[[3]][2,"seq"] %>% unlist()
        
        
        
      }
      
      
      
      # 没有...
      if (result[[1]][1,"start"] ==  human_reigon_start & result[[3]][1,"end"] ==  human_reigon_end & 
          (as.numeric(result[[1]][1,"end"])+1 == result[[2]][1,"start"]) & (as.numeric(result[[2]][1,"end"])+1 == result[[3]][1,"start"]) ) {
        
        
        human_seq <- if_else(strand == -1,paste0(unlist(result[[3]][1,"seq"]),unlist(result[[2]][1,"seq"]),unlist(result[[1]][1,"seq"]) ), 
                             paste0(unlist(result[[1]][1,"seq"]),unlist(result[[2]][1,"seq"]),unlist(result[[3]][1,"seq"]) )   )
        
        
        ailgn_seq <- if_else(strand == -1,paste0(unlist(result[[3]][2,"seq"]),unlist(result[[2]][2,"seq"]),unlist(result[[1]][2,"seq"]) ), 
                             paste0(unlist(result[[1]][2,"seq"]),unlist(result[[2]][2,"seq"]),unlist(result[[3]][2,"seq"]) ))
        
        # mult_alig_human_seq <- ailgn_seq
        # output_mult_alig_fasta[1,1] <- ">human"
        # output_mult_alig_fasta[1,1] <- human_seq
        # output_mult_alig_fasta[1,1] <- paste0(">",name_table[name_table$Scientific.name == scientific_name,"Common.name"] )
        # output_mult_alig_fasta[1,1] <- mult_alig_human_seq
        
        
      # 有...  暂时先stop() 停止运行 遇到了再完善
      }else {
        
        print("error: multiple reigon == 3,exit ....")
        stop()
      }
      
      
      
      # multiple reigon >3的情况
    }else{
      
      # print("error:multiple reigon >3")
      # stop()
      
      
      print(paste0(name_table[name_table$Scientific.name == scientific_name,"Common.name"],
                   ":   multiple reigon >3!,结果另存在reigon大于3的比对结果.xlsx表格中" ) )
      
      return(result)

    }
    
    
    
    
    
    output_align_fasta[1,1] <- ">human"
    output_align_fasta[2,1] <- human_seq
    # output_align_fasta[3,1] <- scientific_name
    output_align_fasta[3,1] <- paste0(">",name_table[name_table$Scientific.name == scientific_name,"Common.name"] )
    output_align_fasta[4,1] <- ailgn_seq
    
    
    print(paste0("success pairwise alignment","   ",name_table[name_table$Scientific.name == scientific_name,"Common.name"],"   ",scientific_name))
    
    
    alig_data <- list(output_align_fasta,output_align_reigon)
    # assign( "alig_data",value = alig_data,envir = .GlobalEnv)
    
    
    
    return(alig_data)
    

  }
  

}


fix_fa <- function(fa){
  
  
  for (i in 1:nrow(fa)) {
    
    
    if (i %% 4 == 1) {
      
      len <- nchar(fa[i+1,1])
      fa[i,1] <- paste0("homo_sapiens/1-",len)
      
    }
    
    if (i %% 4 == 3) {
      
      len <- nchar(fa[i+1,1])
      fa[i,1] <- paste0(fa[i,1],"/1-",len)
      
    }
    
  }
  
  return(fa)
  
}


get_mult_reigon_species_full_seq <- function(mult_reigon = mult_reigon,output_path = output_path,name_table){
  
  if (file.exists(paste0(output_path,"mult_reigon_species_full_seq.txt"))) {
    file.remove( paste0(output_path,"mult_reigon_species_full_seq.txt") )
  }
  
  # write.table(NULL,file = paste0(output_path,"mult_reigon_species_full_seq.txt"),quote = F,row.names = F,col.names = F)
  
  if (nrow(mult_reigon) >0) {
    
    for (i in 1:nrow(mult_reigon)) {
      
      if (is.null(mult_reigon[i,"seq_region_3"]) || is.na(mult_reigon[i,"seq_region_3"]) ) {
        if (mult_reigon[i,"seq_region_1"] == mult_reigon[i,"seq_region_2"]) {
          

          full_reigon <- paste0(mult_reigon[i,"seq_region_1"],":",mult_reigon[i,"start_1"],"-",mult_reigon[i,"end_2"],":",mult_reigon[i,"strand_1"])
          
          
          if (mult_reigon[i,"strand_1"] == -1) {
            
            full_reigon <- paste0(mult_reigon[i,"seq_region_1"],":",mult_reigon[i,"start_2"],"-",mult_reigon[i,"end_1"],":",mult_reigon[i,"strand_1"])
            
          }
          

          
          print(paste0(mult_reigon[i,"scientific_name"]," ",full_reigon))
          
          
          full_seq <- reigon_to_seq(reigon = full_reigon,scientific_name = mult_reigon[i,"scientific_name"],name_table = name_table)
        
          # if (str_ends(full_seq,"fail to get seq")) {
          #   print(full_seq)
          #   next
          # }
          
          
          
          full_seq_len <- nchar(full_seq)
          
          id <- paste0(">",mult_reigon[i,"common_name"],"/1-",full_seq_len )
          
          write.table(id,file = paste0(output_path,"mult_reigon_species_full_seq.txt"),quote = F,row.names = F,col.names = F,append = T)
          write.table(full_seq[2,1],file = paste0(output_path,"mult_reigon_species_full_seq.txt"),quote = F,row.names = F,col.names = F,append = T)
          
          
        }else{
          print("reigon1 reigon2 site in different chromesome")
        }
        
      }else{
        if (mult_reigon[i,"seq_region_1"] == mult_reigon[i,"seq_region_2"] & mult_reigon[i,"seq_region_2"] == mult_reigon[i,"seq_region_3"]) {
          
          full_reigon <- paste0(mult_reigon[i,"seq_region_1"],":",mult_reigon[i,"start_1"],"-",mult_reigon[i,"end_3"],":",mult_reigon[i,"strand_1"])
          
          
          if (mult_reigon[i,"strand_1"] == -1) {
            
            full_reigon <- paste0(mult_reigon[i,"seq_region_1"],":",mult_reigon[i,"start_3"],"-",mult_reigon[i,"end_1"],":",mult_reigon[i,"strand_1"])
            
          }
          
          
          
          
          # full_seq <- reigon_to_seq(reigon = full_reigon,scientific_name = mult_reigon[i,"scientific_name"])
          
          full_seq <- reigon_to_seq(reigon = full_reigon,scientific_name = mult_reigon[i,"scientific_name"],name_table = name_table)
    
          
          
          
          
          # if (str_ends(full_seq,"fail to get seq")) {
          #   print(full_seq)
          #   next
          # }
          
          
          
          full_seq_len <- nchar(full_seq)
          
          write.table(id,file = paste0(output_path,"mult_reigon_species_full_seq.txt"),quote = F,row.names = F,col.names = F,append = T)
          write.table(full_seq[2,1],file = paste0(output_path,"mult_reigon_species_full_seq.txt"),quote = F,row.names = F,col.names = F,append = T)
          
          
          
        }else{
          print("reigon1 reigon2 reigon3 site in different chromesome")
        }
        
      }
      
    }  
  }
  
}



















