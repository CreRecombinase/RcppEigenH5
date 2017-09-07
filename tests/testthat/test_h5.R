context("h5 compatability")
# #
test_that("calcAF works as expected",{
  test_mat <-matrix(runif(9*8),9,8)
  test_means <- colMeans(test_mat)/2
  tfile <- tempfile()
  write_mat_h5(tfile,"test","geno",data=test_mat)
  sub_i <- c(1,3,5,7)
  sub_means <- calc_af(tfile,"test","geno",index = sub_i,chunksize = 2,check_dup = F)
  expect_equal(test_means[sub_i],sub_means)
  sub_i <- c(1:8)
  sub_means <- calc_af(tfile,"test","geno",index = sub_i,chunksize = 1,check_dup = T)
  tmi <- test_mat
  #  tmi[] <- sample(0:1,9*8,replace = T)
  #  tmd <- duplicated(tmi,MARGIN = 2)
  #
  # # test_means[tmd] <- 0
  # # expect_equal(test_means[sub_i],sub_means)
  # # sub_means<-calc_af(tfile,"test","geno",index = sub_i,chunksize = 1,check_dup = T)
  # # expect_equal(test_means[sub_i],sub_means)

})







test_that("calcvar works as expected",{
  test_mat <-matrix(runif(9*8),9,8)
  test_vars <- apply(test_mat,2,var)
  tfile <- tempfile()
  write_mat_h5(tfile,"test","geno",data=test_mat)
  sub_i <- c(1,3,5,7)
  sub_vars <- calc_var(tfile,"test","geno",index = sub_i,chunksize = 2)
  expect_equal(test_vars[sub_i],sub_vars)
  sub_i <- c(1:8)
  sub_vars <- calc_var(tfile,"test","geno",index = sub_i,chunksize = 1)
  expect_equal(test_vars[sub_i],sub_vars)

})

test_that("calc_yh works as expected",{
  n <- 9
  p <- 8
  test_mat <- matrix(runif(n*p),n,p)
  tpi <- 0.5
  Z <-rbinom(n = p,size = 1,prob = tpi)
  beta <-numeric(p)
  beta[Z==1] <- rnorm(n = sum(Z))
  yh <- c(scale(test_mat,center = T,scale = F)%*%beta)
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
  yh <- c(scale(test_mat,center=T,scale=F)%*%beta)
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
  yh <- c(scale(test_mat,center=T,scale=F)%*%beta)
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

