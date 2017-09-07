context("Data Frames")

test_that("We can read and write integer dataframes",{
  library(dplyr)
  tf <- tempfile()
  tdf <-data_frame(a=1L:10L)
  write_df_h5(tdf,"test_df",tf)
  rdf <- read_df_h5(tf,"test_df")
  expect_equal(tdf,rdf)
  tdf <- mutate(tdf,b=2L:11L)
  tf2 <- tempfile()

  write_df_h5(tdf,"test_df3",tf2)
  rdf <- read_df_h5(tf2,"test_df3")
  expect_equal(tdf,rdf)
})

test_that("We can read and write string dataframes",{
  library(dplyr)
  tf <- tempfile()
  tdf <-data_frame(a=letters)
  write_df_h5(tdf,"test_df",tf)
  rdf <- read_df_h5(tf,"test_df")
  expect_equal(tdf,rdf)
  tdf <- mutate(tdf,b=1:n())
  tf2 <- tempfile()

  write_df_h5(tdf,"test_df3",tf2)
  rdf <- read_df_h5(tf2,"test_df3")
  expect_equal(tdf,rdf)
})
