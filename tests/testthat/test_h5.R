context("h5 compatability")




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


test_that("matrices written from h5 and RcppEIgen  read the same (by RcppEigenH5)",{
  # library(h5)
   test_mat <-matrix(runif(9*8),9,8)
  # tfile <- tempfile()
  # th5f <- h5file(tfile,'a')
  # th5f['test'] <- test_mat
  # h5close(th5f)
  tf2 <- tempfile()
  write_mat_h5(tf2,"/","test",test_mat)
  # read_mat <- read_2d_h5(tfile,"/","test",c(0L,0L),c(9L,8L))
  read_mat_2 <-read_2d_h5(tf2,"/","test",c(0L,0L),c(9L,8L))
  expect_equal(test_mat,read_mat)
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
  th5f <- h5file(tfile,'a')
  th5f['test'] <- test_mat
  h5close(th5f)
  sub_read_mat <- read_2d_index_h5(tfile,"/","test",as.integer(c(1,3,5)))
  expect_equal(sub_read_mat,test_mat[,c(1,3,5)])

  rc_sub_mat <- RColumbo::read_ind_h5(tfile,"/","test",c(1,3,5))
  expect_equal(sub_read_mat,rc_sub_mat)

  sub_read_mat <- read_2d_index_h5(tfile,"/","test",as.integer(c(2,3,6)))
  rc_sub_mat <- RColumbo::read_ind_h5(tfile,"/","test",c(2,3,6))
  expect_equal(sub_read_mat,test_mat[,c(2,3,6)])
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
