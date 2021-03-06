

dir2_hdf5 <- function(base_dir,out_file,deflate_level=8L){

  sub_files <- dir(base_dir,full.names = T)
  sub_dirs<- file.info(sub_files) %>%
    dplyr::mutate(fname=rownames(.)) %>%
    dplyr::filter(isdir) %>%
    select(fname) %>%
    mutate(groupname=basename(fname))



  load_obj <- function(f)
  {
    env <- new.env()
    nm <- load(f, env)[1]
    env[[nm]]
  }
  sload_obj <- purrr::possibly(load_obj,otherwise=NULL)
  swrite_ivec <- purrr::safely(write_ivec_h5)
  sub_dirl <- split(sub_dirs,sub_dirs$groupname)
  purrr::map_df(sub_dirl,function(tsub_dirl){
    cat(tsub_dirl$groupname,"\n")
    if(!purrr::possibly(group_exists,otherwise=F)(out_file,tsub_dirl$groupname)){
      all_files <- dir(tsub_dirl$fname,full.names = T,pattern = "*Rdata")
      gene_names <-gsub("(.+).Rdata","\\1",basename(all_files))
      prep_dataf <- data_frame(inf=all_files,h5file=out_file,dataname=gene_names,groupname=tsub_dirl$groupname) %>% mutate(input=1:n())
      tres <- group_by(prep_dataf,input) %>% do({
        iv <- sload_obj(.$inf)
        if(!is.null(iv)){
          res <- swrite_ivec(h5file=.$h5file,groupname = .$groupname,dataname = .$dataname,data = iv,deflate_level=deflate_level)
        }else{
          res <- list(error=paste0(.$inf,"Not Found"))
        }
        data_frame(filename=.$dataname,write_succeed=is.null(res$error))
      })
    }
  })




}


read_attr <- function(h5filename,datapath,attrname){
  requireNamespace("h5")
  h5f <- h5::h5file(h5filename,'r')
  datag <- h5f[datapath]
  retval <- h5::h5attr(datag,attrname)
  h5::h5close(h5f)
  return(retval)
}

read_vec <- function(h5filename,datapath){
  requireNamespace("h5")
  h5f <- h5::h5file(h5filename,'r')
  data <- h5f[datapath][]
  h5::h5close(h5f)
  return(data)
}

write_vec <- function(h5filename,groupname,dataname,data,deflate_level=0L){
  stopifnot(is.null(dim(data)))
  td <- typeof(data)
  wfun <- NULL
  if(td=="integer"){
    wfun <- write_ivec_h5
  }
  if(td=="character"){
    wfun <-write_svec_h5
  }
  if(td=="double"){
    wfun <- write_dvec_h5
  }
  if(is.null(wfun)){
    stop("vec must be of integer, character, or double type")
  }
  wfun(h5file=h5filename,groupname=groupname,dataname=dataname,data=data,deflate_level=deflate_level)
}

read_df_h5 <- function(h5filepath,groupname=NULL,subcols=NULL,filtervec=NULL){
  require(h5)
  require(dplyr)
  require(lazyeval)
  stopifnot(file.exists(h5filepath))

  hf <- h5file(h5filepath,mode = 'r')
  group <- NULL
  if(!is.null(groupname)){
    stopifnot(existsGroup(hf,groupname))
    group <- openGroup(hf,groupname)
    dsets <- list.datasets(group)
  }else{

    group <- openGroup(hf,groupname)
    dsets <- list.datasets(hf)
  }
  if(!is.null(subcols)){
    if(!is.null(groupname)){
      subcols <- paste0("/",groupname,"/",subcols)
    }else{
      subcols <- paste0("/",subcols)
    }
    dsets <- dsets[dsets %in% subcols]
  }
  stopifnot(length(dsets)>0)
  if(!is.null(groupname)){
    dsnames <- gsub(pattern = paste0("/",groupname,"/"),"",dsets)
  }else{
    dsnames <-gsub("/","",dsets)
  }
  retdf <- as_data_frame(setNames(lapply(dsets,function(x,file,fvec){
    if(is.null(fvec)){
      return(x=c(file[x][]))
    }else{
      if(fvec[length(fvec)]==(fvec[1]+length(fvec)-1)){
        return(file[x][fvec])
      }else{
        return(file[x][][fvec])
      }
    }
  },file=hf,fvec=filtervec),dsnames)
  )
  if(!is.null(group)){
    h5close(group)
  }
  h5close(hf)

  return(retdf)
}



transpose_mat <- function(infilename,ingroupname,indataname,outfilename,outgroupname,outdataname,chunksize=100000,index=NULL){
  library(BBmisc)
  library(dplyr)
  tot_col <-get_colnum_h5(infilename[1],ingroupname,indataname)
  if(is.null(index)){
    index <- 1:tot_col
  }
  if(length(infilename)>1){
    stopifnot(length(infilename)==length(index))
    tot_col<-sum(lengths(index))
    tot_row <- get_rownum_h5(infilename[1],ingroupname,indataname)
  }else{
    tot_col <- length(index)
    tot_row <- get_rownum_h5(infilename,ingroupname,indataname)
    tot_chunks <- max(c(ceiling(tot_col/chunksize),2))
  }

  create_mat_dataset_h5(h5file = outfilename,
                        groupname = outgroupname,
                        dataname = outdataname,
                        dims = as.integer(c(tot_row,tot_col)),
                        chunkdims = c(tot_row,min(c(10000,tot_col))),
                        deflate_level = 4L,
                        doTranspose = T)
  tot_colo <- get_colnum_h5(outfilename,groupname = outgroupname,dataname = outdataname)
  stopifnot(tot_colo==tot_col)
  ichunksl <- list()
  index_dfl <- list()
  if(length(infilename)>1){
    for(j in 1:length(infilename)){
      cat(j,"\n")
      #      ichunksl[[i]] <-chunk(index[[j]],chunk.size = chunksize)
      chunknum <-max(c(ceiling(length(index[[j]])/chunksize),2))
      index_dfl[[j]] <- data_frame(index=index[[j]],chrom=j) %>% mutate(data_ind=1:n(),chunks=cut(data_ind,chunknum,labels=F,include.lowest=T))
    }
    index_df <- bind_rows(index_dfl) %>% mutate(data_ind=1:n())
    ichunks <- split(index_df,list(index_df$chrom,index_df$chunks))
  }else{
    index_df <-data_frame(index=index) %>% mutate(data_ind=1:n(),chunks=cut(data_ind,tot_chunks,labels=F,include.lowest=T))
    ichunks <- split(index_df,index_df[["chunks"]])
  }

  for(i in 1:length(ichunks)){
    cat(i," of ",length(ichunks),"\n")
    chunk <-ichunks[[i]][["index"]]
    writechunk <- ichunks[[i]][["data_ind"]]
    if(length(infilename)>1){
      tinfilename <- infilename[ichunks[[i]][["chrom"]][1]]
    }else{
      tinfilename <- infilename
    }

    datam <-read_2d_index_h5(tinfilename,ingroupname,indataname,indvec = chunk)
    write_mat_chunk_h5(outfilename,outgroupname,outdataname,datam,offsets = as.integer(c(0L,writechunk[1]-1)))
  }

}

libpath <- function() {
  cat(sprintf(
    "%s/RcppEigenH5/libs/RcppEigenH5%s",
    installed.packages()["RcppEigenH5","LibPath"][1],
    .Platform$dynlib.ext
  ))
}






read_ccs_h5 <- function(h5filename,groupname,dataname="data",iname="ir",pname="jc",isSymmetric=NULL,dims=NULL){
  # require("h5")
  require("Matrix")

  stopifnot(file.exists(h5filename))
  h5gn <- list_groups_h5(h5filename)
  stopifnot(groupname %in% h5gn)

  # h5g <- h5f[groupname]
  h5_attrs <- list_attrs_h5(h5file = h5filename,base_group = groupname)
  if("Layout" %in% h5_attrs){

    layout <- read_group_attr_h5(h5file = h5filename,groupname = groupname,attr_name = "Layout")
    stopifnot(layout=="CCS")
  }

  data <- read_dvec(h5filename,groupname,dataname)
  i <- read_ivec(h5filename,groupname,iname)
  p <- read_ivec(h5filename,groupname,pname)
  if(is.null(dims)){
    rows <- NULL
    cols <- NULL
    if("rows" %in% h5_attrs){
      rows <- read_igroup_attr_h5(h5filename,groupname,"rows")
    }
    if("cols" %in% h5_attrs){
      cols <- read_igroup_attr_h5(h5filename,groupname,"cols")
    }
    if(is.null(rows)){
      rows <- cols
    }
    if(is.null(cols)){
      cols <- rows
    }
    dims <- c(rows,cols)
    if("Dimensions" %in% h5_attrs){
      dims <- read_group_iarray_attr_h5(h5filename,groupname,"Dimensions")
    }else{
      stop("Must specify rows,cols or Dimensions with sparse matrix storage")
    }
  }
  return(Matrix::sparseMatrix(i=i,p=p,x=data,index1 = F,dims = dims))
}

#
# write_data_fun(h5file,h5group,h5dataname,data,deflate_level=4,chunksize=1000L){
#   stopifnot(typeof(data)!="list")
#   if(!is.null(dim(data))){
#     stopifnot(length(dim(data))==2)
#     if(typeof(data)=="double"){
#       write_mat_h5(h5file,h5group,h5dataname,data,deflate_level = deflate_level)
#     }else{
#       stop("writing of matrices only supported for doubles(currently)")
#     }
#   }else{
#     if(typeof(data)=="double"){
#       write_dvec_h5(h5file,h5group,h5dataname,data,deflate_level = deflate_level)
#     }else{
#       if(typeof(data)=="character"){
#         stop("Writing of strings not currently supported")
#       }else{
#         if(typeof(data)=="integer"){
#           write_ivec_h5(h5file,h5group,h5dataname,data,deflate_level = deflate_level)
#         }
#       }
#     }
#   }
#
# }

write_df_h5g_alt <- function(df,grpname,deflate_level=4,chunksize=1000L){
  dfl <- nrow(df)
  if(dfl<chunksize){
    chunksize <- nrow(df)
  }

  deflate_level <- as.integer(deflate_level)
  library(h5)
  # dataname <- colnames(df)
  f <-h5file(outfile,mode = 'a')
  stopifnot(groupname!="/")
  if(existsGroup(f,groupname)){
    h5close(f)
    res <- append_df_h5(df,groupname,outfile,deflate_level)
    return(res)
  }else{
    group <- createGroup(f,groupname)
  }
  write_df_h5g(df = df,grp = group,deflate_level = deflate_level,chunksize = chunksize)
}





write_df_h5g <- function(df,grp,deflate_level=4L,chunksize=1000L){
  dataname <- colnames(df)
  dtypes <- sapply(df,typeof)
  dfl <- nrow(df)
  for(i in 1:length(dataname)){
    dsn <- dataname[i]
    td <- df[[dsn]]
    dt <- dtypes[i]
    if(dt!="list"){
      if(is.null(dim(td))){
        tdata <- createDataSet(.Object = group,
                               datasetname = dsn,
                               type = typeof(td),
                               dimensions = length(td),
                               chunksize =chunksize,
                               maxdimensions = NA_integer_,
                               compression = deflate_level)
        tdata[] <- td
      }else{
        tdata <- createDataSet(.Object = group,
                               datasetname = dsn,
                               type = typeof(td),
                               dimensions = dim(td),
                               chunksize =chunksize,
                               maxdimensions = NA_integer_,
                               compression = deflate_level)
        tdata[,] <- td
      }
    }else{
      stopifnot(length(td)==dfl)
      for(j in 1:length(td)){
        # tbasegroup <-paste(groupname,j,sep="/")
        tgrp <- createGroup(group,groupname = as.character(j))
        write_df_h5g(df = td,grp = tgrp)
        h5close(tgrp)
      }
    }
  }
  return(T)
}




write_df_h5 <- function(df,groupname,outfile,deflate_level=4L){
  # library(h5)
  if(file.exists(outfile)){
    stopifnot(!group_exists(outfile,groupname))
  }
  deflate_level <- as.integer(deflate_level)
  # require(h5)
  dataname <- colnames(df)
  # f <-h5file(outfile,mode = 'a')
  # if(existsGroup(f,groupname)){
  #   h5close(f)
  #   res <- append_df_h5(df,groupname,outfile,deflate_level)
  #   return(res)
  # }
  # group <- createGroup(f,groupname)
  for(i in 1:length(dataname)){
    dsn <- dataname[i]
    td <- df[[dsn]]
    write_vec(h5filename = outfile,groupname = groupname,dataname = dsn,data = td,deflate_level = deflate_level)
  }
}
#                       )
#     # tdata <- createDataSet(.Object = group,
#     #                        datasetname = dsn,
#     #                        type = typeof(td),
#     #                        dimensions = length(td),
#     #                        chunksize =chunksize,
#     #                        maxdimensions = NA_integer_,
#     #                        compression = deflate_level)
#     # tdata[] <- td
#   }
#   # h5close(f)
#   return(T)
# }

#
#
# write_df_h5 <- function(df,groupname,outfile,deflate_level=4L,chunksize=1000L){
#   dfl <- nrow(df)
#   if(dfl<chunksize){
#     chunksize <- nrow(df)
#   }
#
#   deflate_level <- as.integer(deflate_level)
#   library(h5)
#   # dataname <- colnames(df)
#   f <-h5file(outfile,mode = 'a')
#   stopifnot(groupname!="/")
#   if(existsGroup(f,groupname)){
#     h5close(f)
#     res <- append_df_h5(df,groupname,outfile,deflate_level)
#     return(res)
#   }else{
#     group <- createGroup(f,groupname)
#   }
#   write_df_h5g(df = df,grp = group,deflate_level = deflate_level,chunksize = chunksize)
#
#
#   h5close(f)
#   return(T)
# }
#
#

