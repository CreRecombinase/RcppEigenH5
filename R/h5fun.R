
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

write_vec <- function(h5filename,vec,datapath,deflate_level=4,chunksize=NULL){
  require(h5)
  h5f <- h5::h5file(h5filename,'a')
  if(is.null(chunksize)){
    chunksize <- length(vec)/2
  }
  data <- h5f[datapath,compression=deflate_level,chunksize=chunksize] <- vec
  h5::h5close(h5f)
}

read_df_h5 <- function(h5filepath,groupname=NULL,subcols=NULL,filtervec=NULL){
  require(h5)
  require(dplyr)
  require(lazyeval)
  stopifnot(file.exists(h5filepath))

  f <- h5file(h5filepath,mode = 'r')
  if(!is.null(groupname)){
    stopifnot(existsGroup(f,groupname))
    group <- openGroup(f,groupname)
    dsets <- list.datasets(group)
  }else{
    dsets <- list.datasets(f)
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
  return(as_data_frame(setNames(lapply(dsets,function(x,file,fvec){
    if(is.null(fvec)){
      return(x=c(file[x][]))
    }else{
      if(fvec[length(fvec)]==(fvec[1]+length(fvec)-1)){
        return(file[x][fvec])
      }else{
        return(file[x][][fvec])
      }
    }
  },file=f,fvec=filtervec),dsnames)
  ))
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






read_ccs_h5 <- function(h5filename,groupname,dataname="data",iname="ir",pname="jc",isSymmetric=NULL){
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
  if(is.null(dims)){
    if("Dimensions" %in% h5_attrs){
      dims <- read_group_iarray_attr_h5(h5filename,groupname,"Dimensions")
    }else{
      stop("Must specify rows,cols or Dimensions with sparse matrix storage")
    }
  }
  return(Matrix::sparseMatrix(i=i,p=p,x=data,index1 = F,dims = dims))
}



write_df_h5 <- function(df,groupname,outfile,deflate_level=4L,chunksize=1000L){
  if(nrow(df)<chunksize){
    chunksize <- nrow(df)
  }
  deflate_level <- as.integer(deflate_level)
  require(h5)
  dataname <- colnames(df)
  f <-h5file(outfile,mode = 'a')
  if(existsGroup(f,groupname)){
    h5close(f)
    res <- append_df_h5(df,groupname,outfile,deflate_level)
    return(res)
  }
  group <- createGroup(f,groupname)
  for(i in 1:length(dataname)){
    dsn <- dataname[i]
    td <- df[[dsn]]
    tdata <- createDataSet(.Object = group,
                           datasetname = dsn,
                           type = typeof(td),
                           dimensions = length(td),
                           chunksize =chunksize,
                           maxdimensions = NA_integer_,
                           compression = deflate_level)
    tdata[] <- td
  }
  h5close(f)
  return(T)
}



