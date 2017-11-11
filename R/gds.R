h5ls_df <- function(h5filename){



}



# gds2HDF5 <- function(gds,h5filename){
#   library(gdsfmt)
#   library(purrr)
#   all_nodes <-ls.gdsn(gds)
#   all_idx <-invoke_map(compose(objdesp.gdsn,index.gdsn),all_nodes,node=gds)
#   all_idx <-invoke_map(compose(objdesp.gdsn,index.gdsn),all_nodes,node=gds) %>% modify_depth(.depth = 2,function(x){x%||% NA}) %>% bind_rows()
#   index.gdsn()
#   array_grps <- filter(all_idx,is.array)
#   folder_grps <- filter(all_idx,!is.array)
#   objd <- objdesp.gdsn(gds,)
# }

gds_walk <- function(gds,path=""){
  gds_node <-index.gdsn(gds,path=path)
  all_nodes <-ls.gdsn(gds_node)
  if(length(all_nodes)==0){
    retl <-objdesp.gdsn(gds_node)
    return(retl)
  }else{
    return(lapply(all_nodes,gds_walk,gds=gds_node))
  }
}

gds2HDF5 <- function(gds,path="",outfile){
  gds_node <-index.gdsn(gds,path=path)
  all_nodes <-ls.gdsn(gds_node)
  if(length(all_nodes)==0){
    retl <-objdesp.gdsn(gds_node)
    if(retl$type!="Folder"){
      retd <- read.gdsn(gds_node)
    }
    return(retl)
  }else{
    return(lapply(all_nodes,gds_walk,gds=gds_node))
  }
}

