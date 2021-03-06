# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

read_2d_boost_h5 <- function(h5file, groupname, dataname) {
    .Call('_RcppEigenH5_read_2d_boost_h5', PACKAGE = 'RcppEigenH5', h5file, groupname, dataname)
}

match_sorted <- function(query, target) {
    .Call('_RcppEigenH5_match_sorted_exp', PACKAGE = 'RcppEigenH5', query, target)
}

rcpsplit <- function(s, delim) {
    .Call('_RcppEigenH5_rcpsplit', PACKAGE = 'RcppEigenH5', s, delim)
}

read_2d_h5 <- function(h5file, groupname, dataname, offset, chunksize) {
    .Call('_RcppEigenH5_read_2d_h5_exp', PACKAGE = 'RcppEigenH5', h5file, groupname, dataname, offset, chunksize)
}

read_1d_h5 <- function(h5file, groupname, dataname, offset, chunksize) {
    .Call('_RcppEigenH5_read_1d_h5', PACKAGE = 'RcppEigenH5', h5file, groupname, dataname, offset, chunksize)
}

read_svec <- function(h5file, groupname, dataname) {
    .Call('_RcppEigenH5_read_svec_exp', PACKAGE = 'RcppEigenH5', h5file, groupname, dataname)
}

read_svec_chunk <- function(h5file, groupname, dataname, offset, chunksize) {
    .Call('_RcppEigenH5_read_svec_chunk_exp', PACKAGE = 'RcppEigenH5', h5file, groupname, dataname, offset, chunksize)
}

read_1d_sindex_h5 <- function(h5file, groupname, dataname, indvec) {
    .Call('_RcppEigenH5_read_1d_sindex_h5', PACKAGE = 'RcppEigenH5', h5file, groupname, dataname, indvec)
}

read_dvec <- function(h5file, groupname, dataname) {
    .Call('_RcppEigenH5_read_dvec_exp', PACKAGE = 'RcppEigenH5', h5file, groupname, dataname)
}

read_ivec <- function(h5file, groupname, dataname) {
    .Call('_RcppEigenH5_read_ivec_exp', PACKAGE = 'RcppEigenH5', h5file, groupname, dataname)
}

read_uivec <- function(h5file, groupname, dataname) {
    .Call('_RcppEigenH5_read_uivec_exp', PACKAGE = 'RcppEigenH5', h5file, groupname, dataname)
}

read_1i_h5 <- function(h5file, groupname, dataname, offset, chunksize) {
    .Call('_RcppEigenH5_read_1i_h5', PACKAGE = 'RcppEigenH5', h5file, groupname, dataname, offset, chunksize)
}

read_2d_mat_h5 <- function(h5file, groupname, dataname) {
    .Call('_RcppEigenH5_read_2d_mat_h5', PACKAGE = 'RcppEigenH5', h5file, groupname, dataname)
}

read_2d_index_h5 <- function(h5file, groupname, dataname, indvec) {
    .Call('_RcppEigenH5_read_2d_index_h5', PACKAGE = 'RcppEigenH5', h5file, groupname, dataname, indvec)
}

read_2d_index_chunk_h5 <- function(h5file, groupname, dataname, indvec, chunksize) {
    .Call('_RcppEigenH5_read_2d_index_chunk_h5', PACKAGE = 'RcppEigenH5', h5file, groupname, dataname, indvec, chunksize)
}

read_1d_index_h5 <- function(h5file, groupname, dataname, indvec) {
    .Call('_RcppEigenH5_read_1d_index_h5', PACKAGE = 'RcppEigenH5', h5file, groupname, dataname, indvec)
}

#'Some documentation
#' @param
#' @export
is_transposed_h5 <- function(h5file, groupname, dataname) {
    .Call('_RcppEigenH5_is_transposed_h5', PACKAGE = 'RcppEigenH5', h5file, groupname, dataname)
}

get_selfpath <- function() {
    .Call('_RcppEigenH5_get_selfpath', PACKAGE = 'RcppEigenH5')
}

h5ls <- function(h5file, base_groupname = "/") {
    .Call('_RcppEigenH5_h5ls', PACKAGE = 'RcppEigenH5', h5file, base_groupname)
}

get_rownum_h5 <- function(h5file, groupname, dataname) {
    .Call('_RcppEigenH5_h5_rownum', PACKAGE = 'RcppEigenH5', h5file, groupname, dataname)
}

get_colnum_h5 <- function(h5file, groupname, dataname) {
    .Call('_RcppEigenH5_h5_colnum', PACKAGE = 'RcppEigenH5', h5file, groupname, dataname)
}

calc_af <- function(h5file, groupname, dataname, index, chunksize, display_progress = TRUE, check_dup = TRUE) {
    .Call('_RcppEigenH5_calc_af', PACKAGE = 'RcppEigenH5', h5file, groupname, dataname, index, chunksize, display_progress, check_dup)
}

calc_summ_h5 <- function(h5file, groupname, dataname, index, chunksize, display_progress = TRUE, check_dup = TRUE) {
    .Call('_RcppEigenH5_calc_summ_h5', PACKAGE = 'RcppEigenH5', h5file, groupname, dataname, index, chunksize, display_progress, check_dup)
}

calc_var <- function(h5file, groupname, dataname, index, chunksize, display_progress = TRUE) {
    .Call('_RcppEigenH5_calc_var', PACKAGE = 'RcppEigenH5', h5file, groupname, dataname, index, chunksize, display_progress)
}

calc_yh <- function(h5file, groupname, dataname, index, beta, chunksize, display_progress = TRUE) {
    .Call('_RcppEigenH5_calc_yh', PACKAGE = 'RcppEigenH5', h5file, groupname, dataname, index, beta, chunksize, display_progress)
}

list_groups_h5 <- function(h5file, base_group = "/") {
    .Call('_RcppEigenH5_h5ls_grp_exp', PACKAGE = 'RcppEigenH5', h5file, base_group)
}

group_exists <- function(h5file, base_group = "/") {
    .Call('_RcppEigenH5_group_exists', PACKAGE = 'RcppEigenH5', h5file, base_group)
}

data_exists <- function(h5file, data_name, base_group = "/") {
    .Call('_RcppEigenH5_data_exists', PACKAGE = 'RcppEigenH5', h5file, data_name, base_group)
}

list_attrs_h5 <- function(h5file, base_group = "/") {
    .Call('_RcppEigenH5_h5ls_attr_exp', PACKAGE = 'RcppEigenH5', h5file, base_group)
}

get_h5_version <- function() {
    .Call('_RcppEigenH5_get_h5_version_exp', PACKAGE = 'RcppEigenH5')
}

read_data_attr_h5 <- function(h5filename, h5_groupname, h5_dataname, h5_attr_name) {
    .Call('_RcppEigenH5_read_data_attr_h5_exp', PACKAGE = 'RcppEigenH5', h5filename, h5_groupname, h5_dataname, h5_attr_name)
}

read_group_attr_h5 <- function(h5file, groupname, attr_name) {
    .Call('_RcppEigenH5_read_group_attr_h5_exp', PACKAGE = 'RcppEigenH5', h5file, groupname, attr_name)
}

read_group_iarray_attr_h5 <- function(h5file, groupname, attr_name) {
    .Call('_RcppEigenH5_read_group_iarray_attr_h5_exp', PACKAGE = 'RcppEigenH5', h5file, groupname, attr_name)
}

read_data_iarray_attr_h5 <- function(h5file, groupname, dataname, attr_name) {
    .Call('_RcppEigenH5_read_data_iarray_attr_h5_exp', PACKAGE = 'RcppEigenH5', h5file, groupname, dataname, attr_name)
}

read_igroup_attr_h5 <- function(h5file, groupname, attr_name) {
    .Call('_RcppEigenH5_read_igroup_attr_h5_exp', PACKAGE = 'RcppEigenH5', h5file, groupname, attr_name)
}

read_idata_attr_h5 <- function(h5file, groupname, dataname, attr_name) {
    .Call('_RcppEigenH5_read_idata_attr_h5_exp', PACKAGE = 'RcppEigenH5', h5file, groupname, dataname, attr_name)
}

write_data_string_attr_h5 <- function(h5filename, h5_groupname, h5_dataname, h5_attr_name, h5_attr_value) {
    invisible(.Call('_RcppEigenH5_write_data_string_attr_h5_exp', PACKAGE = 'RcppEigenH5', h5filename, h5_groupname, h5_dataname, h5_attr_name, h5_attr_value))
}

write_group_string_attr_h5 <- function(h5filename, h5_groupname, h5_attr_name, h5_attr_value) {
    invisible(.Call('_RcppEigenH5_write_group_string_attr_h5_exp', PACKAGE = 'RcppEigenH5', h5filename, h5_groupname, h5_attr_name, h5_attr_value))
}

write_data_int_attr_h5 <- function(h5filename, h5_groupname, h5_dataname, h5_attr_name, h5_attr_value) {
    invisible(.Call('_RcppEigenH5_write_data_int_attr_h5_exp', PACKAGE = 'RcppEigenH5', h5filename, h5_groupname, h5_dataname, h5_attr_name, h5_attr_value))
}

write_group_int_attr_h5 <- function(h5filename, h5_groupname, h5_attr_name, h5_attr_value) {
    invisible(.Call('_RcppEigenH5_write_group_int_attr_h5_exp', PACKAGE = 'RcppEigenH5', h5filename, h5_groupname, h5_attr_name, h5_attr_value))
}

write_mat_chunk_h5 <- function(h5file, groupname, dataname, data, offsets) {
    invisible(.Call('_RcppEigenH5_write_mat_chunk_h5_exp', PACKAGE = 'RcppEigenH5', h5file, groupname, dataname, data, offsets))
}

create_mat_dataset_h5 <- function(h5file, groupname, dataname, dims, chunkdims, deflate_level = as.integer( c(0)), doTranspose = FALSE) {
    invisible(.Call('_RcppEigenH5_create_mat_dataset_h5_exp', PACKAGE = 'RcppEigenH5', h5file, groupname, dataname, dims, chunkdims, deflate_level, doTranspose))
}

write_mat_h5 <- function(h5file, groupname, dataname, data, deflate_level = as.integer( c(0)), doTranspose = FALSE) {
    invisible(.Call('_RcppEigenH5_write_mat_h5_exp', PACKAGE = 'RcppEigenH5', h5file, groupname, dataname, data, deflate_level, doTranspose))
}

write_dvec_h5 <- function(h5file, groupname, dataname, data, deflate_level = as.integer( c(0))) {
    invisible(.Call('_RcppEigenH5_write_dvec_h5_exp', PACKAGE = 'RcppEigenH5', h5file, groupname, dataname, data, deflate_level))
}

write_ivec_h5 <- function(h5file, groupname, dataname, data, deflate_level = as.integer( c(0))) {
    invisible(.Call('_RcppEigenH5_write_ivec_h5_exp', PACKAGE = 'RcppEigenH5', h5file, groupname, dataname, data, deflate_level))
}

write_svec_h5 <- function(h5file, groupname, dataname, data, deflate_level = as.integer( c(0))) {
    invisible(.Call('_RcppEigenH5_write_svec_h5_exp', PACKAGE = 'RcppEigenH5', h5file, groupname, dataname, data, deflate_level))
}

file_rownum <- function(filename) {
    .Call('_RcppEigenH5_file_rownum', PACKAGE = 'RcppEigenH5', filename)
}

file_colnum <- function(filename) {
    .Call('_RcppEigenH5_file_colnum', PACKAGE = 'RcppEigenH5', filename)
}

read_haps <- function(filename, rownum = 0L, colnum = 0L) {
    .Call('_RcppEigenH5_read_haps_exp', PACKAGE = 'RcppEigenH5', filename, rownum, colnum)
}

write_haps_h5 <- function(filename, h5file, groupname, dataname, deflate_level, rownum = 0L, colnum = 0L) {
    invisible(.Call('_RcppEigenH5_write_haps_h5', PACKAGE = 'RcppEigenH5', filename, h5file, groupname, dataname, deflate_level, rownum, colnum))
}

sum_mats <- function(h5files, groupnames, datanames) {
    .Call('_RcppEigenH5_sum_mats', PACKAGE = 'RcppEigenH5', h5files, groupnames, datanames)
}

sparse_matmul <- function(h5files, groupnames, datanames, xmat, do_transpose = TRUE) {
    .Call('_RcppEigenH5_sparse_matmul', PACKAGE = 'RcppEigenH5', h5files, groupnames, datanames, xmat, do_transpose)
}

concat_mat_chunks <- function(h5files, groupnames, datanames) {
    .Call('_RcppEigenH5_concat_mat_chunks', PACKAGE = 'RcppEigenH5', h5files, groupnames, datanames)
}

