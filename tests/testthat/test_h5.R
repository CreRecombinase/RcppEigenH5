context("h5 compatability")
#
test_that("calcAF works as expected",{
  test_mat <-matrix(runif(9*8),9,8)
  test_means <- colMeans(test_mat)/2
  tfile <- tempfile()
  write_mat_h5(tfile,"test","geno",data=test_mat)
  sub_i <- c(1,3,5,7)
  sub_means <- calc_af(tfile,"test","geno",index = sub_i,chunksize = 2)
  expect_equal(test_means[sub_i],sub_means)
  sub_i <- c(1:8)
  sub_means <- calc_af(tfile,"test","geno",index = sub_i,chunksize = 1)
  expect_equal(test_means[sub_i],sub_means)
  sub_means<-calc_af(tfile,"test","geno",index = sub_i,chunksize = 5)
  expect_equal(test_means[sub_i],sub_means)

})


# test_that("calcAF works as expected",{
#   test_mat <-matrix(1:(9*8)+0.0,9,8)
#
#   tfile <- tempfile()
#   write_mat_h5(tfile,"test","geno",data=test_mat)
#   write_mat_wrong_h5(tfile,"test","genoW",data=test_mat)
#
#   mat_wright_rright <- read_2d_mat_h5(tfile,"test","geno")
#   expect_equal(mat_wright_rright,test_mat)
#   mat_wwrong_rwrong <- read_mat_wrong_h5(tfile,"test","genoW",offset = as.integer(c(0,0)),chunksize = as.integer(c(9,8)))
#   sub_mat_wmat <-read_mat_wrong_h5(tfile,"test","genoW",offset = c(1,1),chunksize=as.integer(c(8,7)))
#   expect_equal(mat_wwrong_rwrong,test_mat)
#   tsm <- test_mat[2:9,2:8]
#   expect_equal(test_mat[2:9,2:8],sub_mat_wmat)
# })

#
test_that("calc_yh works as expected",{
  n <- 9
  p <- 8
  test_mat <- matrix(runif(n*p),n,p)
  tpi <- 0.5
  Z <-rbinom(n = p,size = 1,prob = tpi)
  beta <-numeric(p)
  beta[Z==1] <- rnorm(n = sum(Z))
  yh <- c(test_mat%*%beta)
  tfile <- tempfile()
  write_mat_h5(tfile,"test","test",test_mat)
  bii <- which(Z==1)
  t_yh <- calc_yh(h5file = tfile,groupname = "test",dataname = "test",index = bii,beta = beta[bii],chunksize = 2)
  expect_equal(t_yh,yh)


})


test_that("calc_yh works with 'transposed'  data",{
  n <- 9
  p <- 8
  test_mat <- matrix(runif(n*p),n,p)
  tpi <- 0.5
  Z <-rbinom(n = p,size = 1,prob = tpi)
  beta <-numeric(p)
  beta[Z==1] <- rnorm(n = sum(Z))
  yh <- c(test_mat%*%beta)
  tfile <- tempfile()
  write_mat_h5(tfile,"test","test",test_mat,doTranspose = T)
  bii <- which(Z==1)
  t_yh <- calc_yh(h5file = tfile,groupname = "test",dataname = "test",index = bii,beta = beta[bii],chunksize = 2)
  expect_equal(t_yh,yh)
})



test_that("calc_yh works with 'transposed' and subset  data",{
  n <- 9
  p <- 8
  test_mat <- matrix(runif(n*p),n,p)
  tpi <- 0.5
  Z <-rbinom(n = p,size = 1,prob = tpi)
  beta <-numeric(p)
  beta[Z==1] <- rnorm(n = sum(Z))
  yh <- c(test_mat%*%beta)
  tfile <- tempfile()
  write_mat_h5(tfile,"test","test",test_mat,doTranspose = T)
  bii <- which(Z==1)
  t_yh <- calc_yh(h5file = tfile,groupname = "test",dataname = "test",index = bii,beta = beta[bii],chunksize = 2)
  expect_equal(t_yh,yh)
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


test_that("matrices written from h5 and RcppEigenH5  read the same (by RcppEigenH5)",{
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
  write_mat_h5(tf2,"/","test_T",test_mat,deflate_level = 4,doTranspose = T)
  # read_mat <- read_2d_h5(tfile,"/","test",c(0L,0L),c(9L,8L))
  read_mat_2 <-read_2d_h5(tf2,"/","test",c(0L,0L),c(n,p))
  read_mat_t <-read_2d_h5(tf2,"/","test_T",c(0L,0L),c(n,p))
  expect_equal(test_mat,read_mat_2)
})

test_that("numeric vectors written from h5 and RcppEigenH5  are read the same (by RcppEigenH5)",{
  library(h5)
  test_fvec <-runif(9*8)
  test_ivec <- sample(1:100,80)
  # tfile <- tempfile()
  # th5f <- h5file(tfile,'a')
  # th5f['test'] <- test_mat
  # h5close(th5f)
  tf2 <- tempfile()
  write_dvec_h5(tf2,"/","testf",test_fvec)
  write_ivec_h5(tf2,"/","testi",test_ivec)

  read_fvec <- read_dvec(tf2,"/","testf")
  read_ivec <- read_ivec(tf2,"/","testi")
  expect_equal(read_fvec,test_fvec)
  expect_equal(read_ivec,test_ivec)
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


test_that("Attributes are read and written correctly",{

  tfile <- tempfile()
  # th5f <- h5file(tfile,'a')
  Rattr <- "SNAKE"
  # RColumbo::set_attr(tfile,groupname = "grp",attrname = "data",attrvalue = Rattr)
  RColumbo::set_attr(tfile,groupname = "grp",attrname = "data",attrvalue = Rattr)
  # RColumbo::set_attr(tfile,groupname = "grp/data",attrname = "data",attrvalue = Rattr)
  mattr <- read_group_attr_h5(h5file = tfile,groupname = "grp",attr_name = "data")
  # dattr <- read_data_attr_h5(tfile,"grp","data","data")
  expect_equal(Rattr,mattr)

})
