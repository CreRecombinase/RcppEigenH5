context("text processing/boost")

test_that("line counting works",{


  my_mat <-matrix(runif(9*8),9,8)
  tf <- tempfile()
  write.table(my_mat,tf,sep=" ",col.names=F,row.names=F,quote=F)
  rc <- file_rownum(tf)
  expect_equal(rc,nrow(my_mat))
  file.remove(tf)
})


test_that("col counting works",{


  my_mat <-matrix(sample(1:(9*8)),9,8)
  tf <- tempfile()
  write.table(my_mat,tf,sep=" ",col.names=F,row.names=F,quote=F)
  rc <- file_colnum(tf)
  expect_equal(rc,ncol(my_mat))
  file.remove(tf)
})

test_that("reading haplotypes works",{

  my_mat <-matrix(sample(0:1,(9*8),replace=T),9,8)
  tf <- tempfile()
  write.table(my_mat,tf,sep=" ",col.names=F,row.names=F,quote=F)
  read_dat <-read_haps(tf)
  file.remove(tf)
  expect_equal(my_mat,read_dat)

})


test_that("writing (and reading) haplotypes works",{

  my_mat <-matrix(sample(0:1,(9*8),replace=T),9,8)
  tf <- tempfile()
  write.table(my_mat,tf,sep=" ",col.names=F,row.names=F,quote=F)
  ntf <- tempfile()

  write_haps_h5(filename = tf,h5file = ntf,groupname = "SNPdata",dataname = "haplotype",deflate_level = 4)
  read_mat <- read_2d_mat_h5(ntf,"SNPdata","haplotype")
  file.remove(c(tf,ntf))
  expect_equal(read_mat,my_mat)

  big_hapmat <- read_haps("/media/nwknoblauch/Data/1kg/Simulation/chr1_1000.cases.haps")

})


