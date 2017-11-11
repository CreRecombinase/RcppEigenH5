context("h5 compatability")


test_that("'Sparse Matrix' multiply works",{

  nrep <- 30
  p <- 500

  uhmat <- matrix(rnorm(nrep*p),p,nrep)
  chunksize <- 200
  bp <- matrix(0,p,p)
  chunki <- BBmisc::chunk(1:p,chunk.size = chunksize)
  names(chunki) <- letters[1:length(chunki)]
  imatl <- list()
  gl <- names(chunki)
  dl <- names(chunki)
  tf <- tempfile()
  for(i in 1:length(chunki)){
    ti <- chunki[[i]]
    bp[ti,ti] <- runif(length(ti)*length(ti))
    imatl[[i]] <- bp[ti,ti]
    write_mat_h5(tf,gl[i],dl[i],imatl[[i]])

  }
  true_res_t <-t(bp)%*%uhmat
  cpp_res_t <- RSSp::block_mat_mul(imatl,uhmat,transpose_mat_l = T)
  expect_equal(true_res_t,cpp_res_t)
  nres_t <- sparse_matmul(tf,gl,dl,uhmat,T)

  true_res <-bp%*%uhmat
  cpp_res <- block_mat_mul(imatl,uhmat,transpose_mat_l = F)
  expect_equal(true_res,cpp_res)


  # uh_l <- lapply(indl,function(x,mat)return(mat[x,,drop=F]),mat=uhmat)
  # sm_mul <- function(imatl,tmat){
  #   res_mat <- matrix(0,nrow(tmat),ncol(tmat))
  #   toffset <- c(1,1)
  #   for(i in 1:length(imatl)){
  #     lp <- dim(imatl[[i]])
  #     res_mat[toffset] <- imatl[[i]]%*%tmat[toffset[1],]
  #     toffset <- toffset+lp
  #   }
  #
  # }
})


test_that("Checking for transpose works",{

  test_mat <-matrix(runif(9*8),9,8)
  tfile <- tempfile()
  write_mat_h5(tfile,"test","geno2",data=test_mat,doTranspose = F)
  write_mat_h5(tfile,"test","geno",data=test_mat,doTranspose = T)

  expect_equal(is_transposed_h5(tfile,"test","geno"),T)
  expect_equal(is_transposed_h5(tfile,"test","geno2"),F)
})




test_that("Matrix read/write works with boost",{

  test_mat <-matrix(runif(9*8),9,8)
  tfile <- tempfile()
  write_mat_h5(tfile,"test","geno",data=test_mat,deflate_level = 3)
  rmat <- read_2d_boost_h5(tfile,"test","geno")
  expect_equal(test_mat,rmat)



})



test_that("matrices written and read from h5 reads the same with RcppEigenH5",{
  library(h5)
  test_mat <-matrix(runif(9*8),9,8)
  tfile <- tempfile()
  th5f <- h5file(tfile,'a')
  th5f['test'] <- test_mat
  h5close(th5f)
  read_mat <- read_2d_h5(tfile,"/","test",c(0L,0L),c(9L,8L))
  sub_read_mat <- read_2d_h5(tfile,"/","test",c(1L,1L),c(8L,7L))
  expect_equal(sub_read_mat,test_mat[2:9,2:8])
  expect_equal(read_mat,test_mat)

})


test_that("matrices written by  RcppEigenH5 is read back the same (by RcppEigenH5)",{
  # library(h5)
  n <- 4L
  p <- 3L
   test_mat <-matrix(1:(n*p)+0.0,n,p)
  # tfile <- tempfile()
  # th5f <- h5file(tfile,'a')
  # th5f['test'] <- test_mat
  # h5close(th5f)
  tf2 <- tempfile()
  write_mat_h5(tf2,"/","test",test_mat)
  write_mat_h5(tf2,"/","test_T",test_mat,deflate_level = 4L,doTranspose = T)
  # read_mat <- read_2d_h5(tfile,"/","test",c(0L,0L),c(9L,8L))
  read_mat_2 <-read_2d_h5(tf2,"/","test",c(0L,0L),c(n,p))
  read_mat_t <-read_2d_h5(tf2,"/","test_T",c(0L,0L),c(n,p))
  expect_equal(test_mat,read_mat_2)
})

test_that("numeric vectors written from RcppEigenH5 are read the same (by RcppEigenH5)",{
  test_fvec <-runif(9*8)
  test_ivec <- sample(1:100,80)

  tf2 <- tempfile()
  write_dvec_h5(tf2,"/","testf",test_fvec)
  write_ivec_h5(tf2,"/","testi",test_ivec)

  write_dvec_h5(tf2,"test1","testf",test_fvec)
  write_ivec_h5(tf2,"test1","testi",test_ivec)




  read_fvec <- read_dvec(tf2,"/","testf")
  read_ivec <- read_ivec(tf2,"/","testi")
  expect_equal(read_fvec,test_fvec)
  expect_equal(read_ivec,test_ivec)
})


test_that("list.datasets(h5) works just like h5ls (RcppEigenH5)",{
  test_fvec <-runif(9*8)
  test_ivec <- sample(1:100,80)

  tf2 <- tempfile()
  write_dvec_h5(tf2,"/","testf",test_fvec)
  write_ivec_h5(tf2,"/","testi",test_ivec)

  write_dvec_h5(tf2,"test1","testf",test_fvec)
  write_ivec_h5(tf2,"test1","testi",test_ivec)

  tf <- h5::h5file(tf2,mode="a")
  ld <- h5::list.datasets(tf)
  rld <- h5ls(tf2)
  expect_equal(ld,rld)

  ld <- list.datasets(tf,"test1")
  rld <- h5ls(tf2,"test1")
  expect_equal(ld,rld)



  read_fvec <- read_dvec(tf2,"/","testf")
  read_ivec <- read_ivec(tf2,"/","testi")
  expect_equal(read_fvec,test_fvec)
  expect_equal(read_ivec,test_ivec)
})



test_that("string vectors written from h5 are read the same (by RcppEigenH5)",{
  library(h5)
  tessvec <-rep(sample(letters,sample(1:10),replace=T))
  tessvec <- paste0(tessvec,sample(tessvec),sample(tessvec))
  tfile <- tempfile()
  th5f <- h5file(tfile,'a')
  th5f['test'] <- tessvec
  h5close(th5f)
  r_tessvec <- read_svec(tfile,"/","test")
  expect_equal(r_tessvec,tessvec)
})



test_that("Numeric vectors with NAs can be read and written",{

  myn <-c(1.2,NA,3)
  tf <- tempfile()
  write_dvec_h5(tf,"test","testd",myn)
  rn <- read_dvec(tf,"test","testd")
  expect_equal(rn,myn)
})


test_that("String vectors written from RcppEigen will make up groups as they go",{

  tessvec <-rep(sample(letters,sample(1:10),replace=T))
  tessvec <- paste0(tessvec,sample(tessvec),sample(tessvec))
  tf <- tempfile()
  ## No group
  write_svec_h5(h5file = tf,groupname = "/",dataname = "test",data = tessvec,deflate_level = 4L)
  r_tessvec <- read_svec(tf,"/","test")
  expect_equal(r_tessvec,tessvec)

  ## One group
  write_svec_h5(h5file = tf,groupname = "/one",dataname = "test",data = tessvec,deflate_level = 4L)
  r_tessvec <- read_svec(tf,"one","test")
  expect_equal(r_tessvec,tessvec)
  ## Two groups
  write_svec_h5(h5file = tf,groupname = "/one/two",dataname = "test",data = tessvec,deflate_level = 4L)
  r_tessvec <- read_svec(tf,"one/two","test")
  expect_equal(r_tessvec,tessvec)


})






test_that("matrices written from RcppEigen is read back the same (by h5)",{
  library(h5)
  test_mat <-matrix(runif(9*8),9,8)
  tf2 <- tempfile()
  write_mat_h5(tf2,"/","test",test_mat)
  th5f <- h5file(tf2,'a')
  read_mat <- th5f['test'][]
  expect_equal(read_mat,test_mat)
})







test_that("numeric vectors written and read from h5 reads the same with RcppEigenH5",{
  library(h5)
  test_vec <-runif(9*8)
  tfile <- tempfile()
  th5f <- h5file(tfile,'a')
  th5f['test'] <- test_vec
  h5close(th5f)
  read_vec <- read_1d_h5(tfile,"/","test",0L,72L)
  sub_read_vec <- read_1d_h5(tfile,"/","test",1L,71L)
  expect_equal(read_vec,test_vec)
  expect_equal(sub_read_vec,test_vec[2:72])
})


test_that("integer vectors written and read from h5 reads the same with RcppEigenH5",{

  library(h5)
  test_vec <-sample(1:72,9*8,replace = T)
  tfile <- tempfile()
  th5f <- h5file(tfile,'a')
  th5f['test'] <- test_vec
  h5close(th5f)
  read_vec <- read_1i_h5(tfile,"/","test",0L,72L)
  sub_read_vec <- read_1i_h5(tfile,"/","test",1L,71L)
  expect_equal(read_vec,test_vec)
  expect_equal(sub_read_vec,test_vec[2:72])
})

test_that("subsetting matrix columns works as expected",{


  library(h5)

  test_mat <-matrix(runif(9*8),9,8)
  tfile <- tempfile()

   write_mat_h5(h5file = tfile,groupname = "/",dataname = "test",deflate_level = 0L,data = test_mat,doTranspose = T)
  # th5f <- h5file(tfile,'a')
  # th5f['test'] <- test_mat
  # h5close(th5f)
  sub_read_mat <- read_2d_index_h5(tfile,"/","test",as.integer(c(1,3,5)))
  expect_equal(sub_read_mat,test_mat[,c(1,3,5)])

  # rc_sub_mat <- RColumbo::read_ind_h5(tfile,"/","test",c(1,3,5))
  # expect_equal(sub_read_mat,rc_sub_mat)

  sub_read_mat <- read_2d_index_h5(tfile,"/","test",as.integer(c(2,3,6)))
  rc_sub_mat <- RColumbo::read_ind_h5(tfile,"/","test",c(2,3,6))
  expect_equal(sub_read_mat,test_mat[,c(2,3,6)])
})



test_that("'transposing' a matrix doesn't change how it's read",{
  library(h5)
  n <- 9
  p <- 8000
  test_mat <- matrix(1:n*p+0.0,n,p)
  tfile <- tempfile()
  write_mat_h5(tfile,"test","test",test_mat,deflate_level = 1L)
  write_mat_h5(tfile,"test","testT",test_mat,doTranspose = T)
  transpose_mat(infilename = tfile,
                ingroupname = "test",
                indataname = "testT",
                outfilename = tfile,
                outgroupname = "test",
                outdataname = "testT2",chunksize = 300)
  read_t <- read_2d_mat_h5(tfile,"test","testT")
  read_t2 <- read_2d_mat_h5(tfile,"test","testT2")
  read_r <- read_2d_mat_h5(tfile,"test","test")
  expect_equal(test_mat,read_t2)
  expect_equal(test_mat,read_t)
  expect_equal(read_r,read_t)
  expect_equal(test_mat,read_r)
})


test_that("'transposing' and subsetting a matrix doesn't change how it's read",{
  library(h5)
  n <- 9
  p <- 8000
  test_mat <- matrix(1:n*p+0.0,n,p)
  tfile <- tempfile()
  write_mat_h5(tfile,"test","test",test_mat,deflate_level = 1L)
  write_mat_h5(tfile,"test","testT",test_mat,doTranspose = T)
  transpose_mat(infilename = tfile,
                ingroupname = "test",
                indataname = "testT",
                outfilename = tfile,
                outgroupname = "test",
                outdataname = "testT2",chunksize = 300)
  test_ind <- as.integer(seq(1,p,length.out = 100))
  transpose_mat(infilename = tfile,
                ingroupname = "test",
                indataname = "testT",
                outfilename = tfile,
                outgroupname = "test",
                outdataname = "testT3",chunksize = 30,index = test_ind)


  read_t <- read_2d_mat_h5(tfile,"test","testT")
  read_t2 <- read_2d_mat_h5(tfile,"test","testT2")
  read_r <- read_2d_mat_h5(tfile,"test","test")

  sub_ind <- test_mat[,test_ind]
  read_t3 <-read_2d_mat_h5(tfile,"test","testT3")

  expect_equal(sub_ind,read_t3)
  expect_equal(test_mat,read_t2)
  expect_equal(test_mat,read_t)
  expect_equal(read_r,read_t)
  expect_equal(test_mat,read_r)
})

test_that("'transposing' and subsetting multiple matrices doesn't change how they're read",{
  n <- 9
  p <- 8000
  test_mat_1 <- matrix(runif(n*(p-1)),n,p-1)
  tfile_1 <- tempfile()
  test_mat_2 <- matrix(runif(n*(p+2)),n,p+2)
  tfile_2 <- tempfile()
  test_mat_3 <- matrix(runif(n*(p+4)),n,p+4)
  tfile_3 <- tempfile()

  write_mat_h5(tfile_1,"test","test",test_mat_1,deflate_level = 1L)
  write_mat_h5(tfile_2,"test","test",test_mat_2,deflate_level = 1L)
  write_mat_h5(tfile_3,"test","test",test_mat_3,deflate_level = 1L)
  sub_i_1 <- sort(sample(1:(p-1),4,replace = F))
  sub_i_2 <- sort(sample(1:(p+2),6,replace = F))
  sub_i_3 <- sort(sample(1:(p+4),2,replace = F))
  exp_mat <- cbind(test_mat_1[,sub_i_1],test_mat_2[,sub_i_2],test_mat_3[,sub_i_3])
  ofile <- tempfile()
  transpose_mat(infilename = c(tfile_1,tfile_2,tfile_3),
                ingroupname = "test",
                indataname = "test",
                outfilename = ofile,
                outgroupname = "test",
                outdataname = "test",chunksize = 3,index = list(sub_i_1,sub_i_2,sub_i_3))
  rmat <-read_2d_mat_h5(ofile,"test","test")
  expect_equal(exp_mat,rmat)
})
test_that("subsetting matrix columns (chunked/transposed) works as expected",{


  library(h5)
  n <- 9
  p <- 8000
  test_mat <- matrix(runif(n*p),n,p)
  tfile <- tempfile()
  write_mat_h5(tfile,"test","test",test_mat)
  write_mat_h5(tfile,"test","testT",test_mat,doTranspose = T)
  sind <- as.integer(c(1,3,7000))
  sub_read_mat <- read_2d_index_h5(tfile,"test","testT",sind)
  expect_equal(sub_read_mat,test_mat[,sind])

  csub_read_mat <- read_2d_index_chunk_h5(tfile,"test","testT",sind,chunksize = 100L)
  expect_equal(sub_read_mat,csub_read_mat)
  csub_read_mat <- read_2d_index_chunk_h5(tfile,"test","testT",sind,chunksize = 1L)
  expect_equal(sub_read_mat,csub_read_mat)
})




test_that("subsetting matrix columns (chunked) works as expected",{


  library(h5)
  n <- 9
  p <- 8000
  test_mat <- matrix(runif(n*p),n,p)
  tfile <- tempfile()
  write_mat_h5(tfile,"test","test",test_mat)
  sind <- as.integer(c(1,3,7000))
  sub_read_mat <- read_2d_index_h5(tfile,"test","test",sind)
  expect_equal(sub_read_mat,test_mat[,sind])

  csub_read_mat <- read_2d_index_chunk_h5(tfile,"test","test",sind,chunksize = 100L)
  expect_equal(sub_read_mat,csub_read_mat)
  csub_read_mat <- read_2d_index_chunk_h5(tfile,"test","test",sind,chunksize = 1L)
  expect_equal(sub_read_mat,csub_read_mat)
})


test_that("subsetting matrix columns (chunked, transposed) works as expected",{


  library(h5)
  n <- 9
  p <- 8000
  test_mat <- matrix(runif(n*p),n,p)
  tfile <- tempfile()
  write_mat_h5(tfile,"test","test",test_mat,deflate_level = 4L,doTranspose = T)
  sind <- as.integer(c(1,3,7000))
  sub_read_mat <- read_2d_index_h5(tfile,"test","test",sind)
  expect_equal(sub_read_mat,test_mat[,sind])

  csub_read_mat <- read_2d_index_chunk_h5(tfile,"test","test",sind,chunksize = 100L)
  expect_equal(sub_read_mat,csub_read_mat)
  csub_read_mat <- read_2d_index_chunk_h5(tfile,"test","test",sind,chunksize = 1L)
  expect_equal(sub_read_mat,csub_read_mat)
})




test_that("subsetting vectors works as expected",{
  library(h5)
  test_vec <-runif(9*8)
  tfile <- tempfile()
  th5f <- h5file(tfile,'a')
  th5f['test'] <- test_vec
  h5close(th5f)
  sub_vec <- read_1d_index_h5(tfile,"/","test",as.integer(c(1,3,5)))
  expect_equal(sub_vec,test_vec[c(1,3,5)])
})



test_that("Checking for group existence works",{

  n <- 10
  p <- 8000
  test_mat_1 <- matrix(runif(n*p),n,p)
  tfile <- tempfile()
  write_mat_h5(tfile,"test","testmat1",test_mat_1,deflate_level = 4L,doTranspose = T)
  expect_true(group_exists(tfile,"test"))
  expect_true(group_exists(tfile))
  expect_false(group_exists(tfile,"nottest"))
  expect_false(group_exists(tfile,"test/2"))
  expect_false(group_exists(tfile,"test/testmat1"))

})

test_that("Nested groups work",{

  n <- 10
  p <- 8000
  test_mat_1 <- matrix(runif(n*p),n,p)
  test_mat_2 <- matrix(runif(n*p),n,p)

  tfile <- tempfile()
  write_mat_h5(tfile,"test","testmat1",test_mat_1,deflate_level = 4L,doTranspose = T)
  expect_true(group_exists(tfile,"test"))
  write_mat_h5(tfile,"test/2","testmat2",test_mat_2,deflate_level = 4L,doTranspose = T)
  expect_true(group_exists(tfile,"test/2"))
  list_groups_h5(tfile,base_group = "test")
  tt_1 <- read_2d_mat_h5(tfile,"test","testmat1")
  tt_2 <- read_2d_mat_h5(tfile,"test/2","testmat2")
  expect_equal(test_mat_1,tt_1)
  expect_equal(test_mat_2,tt_2)
})



test_that("H5LS works",{

  n <- 10
  p <- 8000
  test_mat_1 <- matrix(runif(n*p),n,p)
  test_mat_2 <- matrix(runif(n*p),n,p)

  tfile <- tempfile()
  write_mat_h5(tfile,"test","testmat1",test_mat_1,deflate_level = 4L,doTranspose = T)
  write_mat_h5(tfile,"test/2","testmat2",test_mat_2,deflate_level = 4L,doTranspose = T)
  write_mat_h5(tfile,"/","testmat0",test_mat_2,deflate_level = 4L,doTranspose = T)
  list_groups_h5(tfile,base_group = "test")
  tt_1 <- read_2d_mat_h5(tfile,"test","testmat1")
  tt_2 <- read_2d_mat_h5(tfile,"test/2","testmat2")
  expect_equal(h5ls(tfile),c("/test/2/testmat2","/test/testmat1","/testmat0"))


})


# test_that("Writing nested DFs works",{
#
#   # n <- 10
#   # p <- 8000
#   # test_mat_1 <- matrix(runif(n*p),n,p)
#   # test_mat_2 <- matrix(runif(n*p),n,p)
#   #
#   # tnest_df <- beaver1 %>% nest(-activ)
#   # tfile <- tempfile()
#   # write_df_h5(tnest_df,groupname="data",outfile=tfile)
#   # list_groups_h5(tfile,base_group = "test")
#   # tt_1 <- read_2d_mat_h5(tfile,"test","testmat1")
#   # tt_2 <- read_2d_mat_h5(tfile,"test/2","testmat2")
#   # expect_equal(h5ls(tfile),c("/test/2/testmat2","/test/testmat1","/testmat0"))
#
#
# })

testthat::test_that("Writing attributes works",{

  n <- 10
  p <- 8000
  test_mat_1 <- matrix(runif(n*p),n,p)
  tfile <- tempfile()
  write_mat_h5(tfile,"test","testmat1",test_mat_1,deflate_level = 4L,doTranspose = T)
  write_data_string_attr_h5(h5filename = tfile,h5_groupname = "test",h5_dataname = "testmat1",h5_attr_name = "testv",h5_attr_value = "test_value")
  read_data_attr_h5(tfile,"test","testmat1","testv")

})



test_that("Reading and writng sparse matrices works",{

  library(Matrix)
  tfile <- tempfile()

  i <- c(1,3:8); j <- c(2,9,6:10); x <- 7.2 * (1:7)
  A <- sparseMatrix(i, j, x = x)
  ai <-A@i
  aj <-A@p
  ax <-A@x

  write_ivec_h5(tfile,"R","ir",data=ai)
  write_ivec_h5(tfile,"R","jc",data=aj)
  write_dvec_h5(tfile,"R","data",data=ax)
  tdims <- dim(A)
  res_mat <- read_ccs_h5(tfile,"R",dims = tdims)
  expect_equal(res_mat,A)
})


test_that("reduce sum operations work ",{

  m <-4
  rownum <- 5
  colnum <- 6
  tfv <- character(m)
  tml <- list()
  R_sum <- matrix(0,rownum,colnum)
  for(i in 1:m){
    tfv[i] <- tempfile()
    tmat <- matrix(runif(rownum*colnum),rownum,colnum)
    tml[[i]] <- tmat
    RcppEigenH5::write_mat_h5(h5file = tfv[i],groupname = as.character(i),dataname = as.character(i),data = tmat)
    R_sum <- R_sum+tmat
  }

  c_sum <- RcppEigenH5::sum_mats(tfv[1],as.character(1:m)[1],as.character(1:m)[1])
  expect_equal(c_sum,tml[[1]])
  c_sum <- RcppEigenH5::sum_mats(tfv,as.character(1:m),as.character(1:m))
  expect_equal(c_sum,R_sum)
})



test_that("Concatenating matrix chunks works",{
  m <-4
  rownum <- 5
  colnum <- 6
  tfv <- character(m)
  tml <- list()
  R_sum <- matrix(0,0,colnum)
  for(i in 1:m){
    tfv[i] <- tempfile()
    tmat <- matrix(runif(rownum*colnum),rownum,colnum)
    tml[[i]] <- tmat
    RcppEigenH5::write_mat_h5(h5file = tfv[i],groupname = as.character(i),dataname = as.character(i),data = tmat)
    R_sum <- rbind(R_sum,tmat)
  }

  c_sum <- concat_mat_chunks(tfv,as.character(1:m),as.character(1:m))
  expect_equal(c_sum,R_sum)

})

