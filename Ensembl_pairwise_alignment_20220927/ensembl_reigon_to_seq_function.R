library(httr)
library(jsonlite)
library(xml2)
library(tidyverse)


# Scientific_name = "Tupaia_belangeri"
# reigon = "scaffold_126628:361533-361551:-1"


name_df <<- read.csv("./all_ensembl_Species_name_no_blank.csv")


reigon_to_seq <- function(human_reigon,scientific_name,seq_id = NULL){
  
  df <- data.frame()
  server <- "https://rest.ensembl.org"
  ext <- paste0("/sequence/region/",scientific_name,"/",human_reigon,"?")
  
  r <- GET(paste(server, ext, sep = ""), content_type("text/plain"))
  r
  
    if (r$status_code == 400) {
    
    print(paste0(scientific_name,": fail to get seq"))
  }else{
    
    
    stop_for_status(r)
    
    #print(content(r))
    
    if (is.null(seq_id)) {
      # df[1,1] <- paste0(">",scientific_name)
      df[1,1] <- paste0(">",name_df[name_df$scientific.name == str_remove(df[1,1],">"),1] )
      
    }else{
      
      df[1,1] <- seq_id
    }

    df[2,1] <- content(r)
        
    return(df)
    
  }
  
  

  # return(list(Scientific_name,reigon,content(r)) )
  # return(content(r) )

}








# seq <- reigon_to_seq(human_reigon = "9:40915717-40992195:-1",scientific_name = "Homo_sapiens")
# seq2 <- reigon_to_seq(human_reigon = "10:26936871-27004115:-1",scientific_name = "Homo_sapiens")
# 
# 
# write.table(seq,file = "seq.txt",quote = F,row.names = F,col.names = F)
# write.table(seq2,file = "seq2.txt",quote = F,row.names = F,col.names = F)


# reigon_to_seq(Scientific_name = "Tupaia_belangeri",reigon = "scaffold_126628:361521-361532:-1") #reigon1
# 
# reigon_to_seq(Scientific_name = "Tupaia_belangeri",reigon = "scaffold_126628:361521-361551:-1") #reigon1的start-reigon2的end
# 
# 
# reigon_to_seq(Scientific_name = "Vulpes_vulpes",reigon = "NBDQ01000010.1:20467005-20467149:-1")


# reigon_to_seq(Scientific_name = "Tupaia_belangeri",reigon = "scaffold_126628:361533-361551:-1")
