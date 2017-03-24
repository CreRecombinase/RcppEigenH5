
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


libpath <- function() {
  cat(sprintf(
    "%s/RcppEigenH5/libs/RcppEigenH5%s",
    installed.packages()["RcppEigenH5","LibPath"][1],
    .Platform$dynlib.ext
  ))
}

