// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/RcppEigenH5.h"
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// read_2d_h5_exp
Eigen::MatrixXd read_2d_h5_exp(const StringVector h5file, const StringVector groupname, const StringVector dataname, const IntegerVector offset, const IntegerVector chunksize);
RcppExport SEXP RcppEigenH5_read_2d_h5_exp(SEXP h5fileSEXP, SEXP groupnameSEXP, SEXP datanameSEXP, SEXP offsetSEXP, SEXP chunksizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const StringVector >::type h5file(h5fileSEXP);
    Rcpp::traits::input_parameter< const StringVector >::type groupname(groupnameSEXP);
    Rcpp::traits::input_parameter< const StringVector >::type dataname(datanameSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type offset(offsetSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type chunksize(chunksizeSEXP);
    rcpp_result_gen = Rcpp::wrap(read_2d_h5_exp(h5file, groupname, dataname, offset, chunksize));
    return rcpp_result_gen;
END_RCPP
}
// read_1d_h5
Eigen::ArrayXd read_1d_h5(const StringVector h5file, const StringVector groupname, const StringVector dataname, const IntegerVector offset, const IntegerVector chunksize);
RcppExport SEXP RcppEigenH5_read_1d_h5(SEXP h5fileSEXP, SEXP groupnameSEXP, SEXP datanameSEXP, SEXP offsetSEXP, SEXP chunksizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const StringVector >::type h5file(h5fileSEXP);
    Rcpp::traits::input_parameter< const StringVector >::type groupname(groupnameSEXP);
    Rcpp::traits::input_parameter< const StringVector >::type dataname(datanameSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type offset(offsetSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type chunksize(chunksizeSEXP);
    rcpp_result_gen = Rcpp::wrap(read_1d_h5(h5file, groupname, dataname, offset, chunksize));
    return rcpp_result_gen;
END_RCPP
}
// read_svec_exp
Rcpp::StringVector read_svec_exp(const StringVector h5file, const StringVector groupname, const StringVector dataname);
RcppExport SEXP RcppEigenH5_read_svec_exp(SEXP h5fileSEXP, SEXP groupnameSEXP, SEXP datanameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const StringVector >::type h5file(h5fileSEXP);
    Rcpp::traits::input_parameter< const StringVector >::type groupname(groupnameSEXP);
    Rcpp::traits::input_parameter< const StringVector >::type dataname(datanameSEXP);
    rcpp_result_gen = Rcpp::wrap(read_svec_exp(h5file, groupname, dataname));
    return rcpp_result_gen;
END_RCPP
}
// read_svec_chunk_exp
Rcpp::StringVector read_svec_chunk_exp(const StringVector h5file, const StringVector groupname, const StringVector dataname, const int offset, const int chunksize);
RcppExport SEXP RcppEigenH5_read_svec_chunk_exp(SEXP h5fileSEXP, SEXP groupnameSEXP, SEXP datanameSEXP, SEXP offsetSEXP, SEXP chunksizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const StringVector >::type h5file(h5fileSEXP);
    Rcpp::traits::input_parameter< const StringVector >::type groupname(groupnameSEXP);
    Rcpp::traits::input_parameter< const StringVector >::type dataname(datanameSEXP);
    Rcpp::traits::input_parameter< const int >::type offset(offsetSEXP);
    Rcpp::traits::input_parameter< const int >::type chunksize(chunksizeSEXP);
    rcpp_result_gen = Rcpp::wrap(read_svec_chunk_exp(h5file, groupname, dataname, offset, chunksize));
    return rcpp_result_gen;
END_RCPP
}
// read_1d_sindex_h5
Rcpp::StringVector read_1d_sindex_h5(const StringVector h5file, const StringVector groupname, const StringVector dataname, const IntegerVector indvec);
RcppExport SEXP RcppEigenH5_read_1d_sindex_h5(SEXP h5fileSEXP, SEXP groupnameSEXP, SEXP datanameSEXP, SEXP indvecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const StringVector >::type h5file(h5fileSEXP);
    Rcpp::traits::input_parameter< const StringVector >::type groupname(groupnameSEXP);
    Rcpp::traits::input_parameter< const StringVector >::type dataname(datanameSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type indvec(indvecSEXP);
    rcpp_result_gen = Rcpp::wrap(read_1d_sindex_h5(h5file, groupname, dataname, indvec));
    return rcpp_result_gen;
END_RCPP
}
// read_dvec_exp
Eigen::ArrayXd read_dvec_exp(const StringVector h5file, const StringVector groupname, const StringVector dataname);
RcppExport SEXP RcppEigenH5_read_dvec_exp(SEXP h5fileSEXP, SEXP groupnameSEXP, SEXP datanameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const StringVector >::type h5file(h5fileSEXP);
    Rcpp::traits::input_parameter< const StringVector >::type groupname(groupnameSEXP);
    Rcpp::traits::input_parameter< const StringVector >::type dataname(datanameSEXP);
    rcpp_result_gen = Rcpp::wrap(read_dvec_exp(h5file, groupname, dataname));
    return rcpp_result_gen;
END_RCPP
}
// read_ivec_exp
Eigen::ArrayXi read_ivec_exp(const StringVector h5file, const StringVector groupname, const StringVector dataname);
RcppExport SEXP RcppEigenH5_read_ivec_exp(SEXP h5fileSEXP, SEXP groupnameSEXP, SEXP datanameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const StringVector >::type h5file(h5fileSEXP);
    Rcpp::traits::input_parameter< const StringVector >::type groupname(groupnameSEXP);
    Rcpp::traits::input_parameter< const StringVector >::type dataname(datanameSEXP);
    rcpp_result_gen = Rcpp::wrap(read_ivec_exp(h5file, groupname, dataname));
    return rcpp_result_gen;
END_RCPP
}
// read_1i_h5
Eigen::ArrayXi read_1i_h5(const StringVector h5file, const StringVector groupname, const StringVector dataname, const IntegerVector offset, const IntegerVector chunksize);
RcppExport SEXP RcppEigenH5_read_1i_h5(SEXP h5fileSEXP, SEXP groupnameSEXP, SEXP datanameSEXP, SEXP offsetSEXP, SEXP chunksizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const StringVector >::type h5file(h5fileSEXP);
    Rcpp::traits::input_parameter< const StringVector >::type groupname(groupnameSEXP);
    Rcpp::traits::input_parameter< const StringVector >::type dataname(datanameSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type offset(offsetSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type chunksize(chunksizeSEXP);
    rcpp_result_gen = Rcpp::wrap(read_1i_h5(h5file, groupname, dataname, offset, chunksize));
    return rcpp_result_gen;
END_RCPP
}
// read_2d_mat_h5
Eigen::MatrixXd read_2d_mat_h5(const StringVector h5file, const StringVector groupname, const StringVector dataname);
RcppExport SEXP RcppEigenH5_read_2d_mat_h5(SEXP h5fileSEXP, SEXP groupnameSEXP, SEXP datanameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const StringVector >::type h5file(h5fileSEXP);
    Rcpp::traits::input_parameter< const StringVector >::type groupname(groupnameSEXP);
    Rcpp::traits::input_parameter< const StringVector >::type dataname(datanameSEXP);
    rcpp_result_gen = Rcpp::wrap(read_2d_mat_h5(h5file, groupname, dataname));
    return rcpp_result_gen;
END_RCPP
}
// read_2d_index_h5
Eigen::MatrixXd read_2d_index_h5(const StringVector h5file, const StringVector groupname, const StringVector dataname, const IntegerVector indvec);
RcppExport SEXP RcppEigenH5_read_2d_index_h5(SEXP h5fileSEXP, SEXP groupnameSEXP, SEXP datanameSEXP, SEXP indvecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const StringVector >::type h5file(h5fileSEXP);
    Rcpp::traits::input_parameter< const StringVector >::type groupname(groupnameSEXP);
    Rcpp::traits::input_parameter< const StringVector >::type dataname(datanameSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type indvec(indvecSEXP);
    rcpp_result_gen = Rcpp::wrap(read_2d_index_h5(h5file, groupname, dataname, indvec));
    return rcpp_result_gen;
END_RCPP
}
// read_2d_index_chunk_h5
Eigen::MatrixXd read_2d_index_chunk_h5(const StringVector h5file, const StringVector groupname, const StringVector dataname, const IntegerVector indvec, const IntegerVector chunksize);
RcppExport SEXP RcppEigenH5_read_2d_index_chunk_h5(SEXP h5fileSEXP, SEXP groupnameSEXP, SEXP datanameSEXP, SEXP indvecSEXP, SEXP chunksizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const StringVector >::type h5file(h5fileSEXP);
    Rcpp::traits::input_parameter< const StringVector >::type groupname(groupnameSEXP);
    Rcpp::traits::input_parameter< const StringVector >::type dataname(datanameSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type indvec(indvecSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type chunksize(chunksizeSEXP);
    rcpp_result_gen = Rcpp::wrap(read_2d_index_chunk_h5(h5file, groupname, dataname, indvec, chunksize));
    return rcpp_result_gen;
END_RCPP
}
// read_1d_index_h5
Eigen::ArrayXd read_1d_index_h5(const StringVector h5file, const StringVector groupname, const StringVector dataname, const IntegerVector indvec);
RcppExport SEXP RcppEigenH5_read_1d_index_h5(SEXP h5fileSEXP, SEXP groupnameSEXP, SEXP datanameSEXP, SEXP indvecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const StringVector >::type h5file(h5fileSEXP);
    Rcpp::traits::input_parameter< const StringVector >::type groupname(groupnameSEXP);
    Rcpp::traits::input_parameter< const StringVector >::type dataname(datanameSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type indvec(indvecSEXP);
    rcpp_result_gen = Rcpp::wrap(read_1d_index_h5(h5file, groupname, dataname, indvec));
    return rcpp_result_gen;
END_RCPP
}
// h5_rownum
Rcpp::IntegerVector h5_rownum(const std::string h5file, const std::string groupname, const std::string dataname);
RcppExport SEXP RcppEigenH5_h5_rownum(SEXP h5fileSEXP, SEXP groupnameSEXP, SEXP datanameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string >::type h5file(h5fileSEXP);
    Rcpp::traits::input_parameter< const std::string >::type groupname(groupnameSEXP);
    Rcpp::traits::input_parameter< const std::string >::type dataname(datanameSEXP);
    rcpp_result_gen = Rcpp::wrap(h5_rownum(h5file, groupname, dataname));
    return rcpp_result_gen;
END_RCPP
}
// h5_colnum
Rcpp::IntegerVector h5_colnum(const std::string h5file, const std::string groupname, const std::string dataname);
RcppExport SEXP RcppEigenH5_h5_colnum(SEXP h5fileSEXP, SEXP groupnameSEXP, SEXP datanameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string >::type h5file(h5fileSEXP);
    Rcpp::traits::input_parameter< const std::string >::type groupname(groupnameSEXP);
    Rcpp::traits::input_parameter< const std::string >::type dataname(datanameSEXP);
    rcpp_result_gen = Rcpp::wrap(h5_colnum(h5file, groupname, dataname));
    return rcpp_result_gen;
END_RCPP
}
// calc_af
Rcpp::NumericVector calc_af(const std::string h5file, const std::string groupname, const std::string dataname, const Eigen::ArrayXi index, const Rcpp::IntegerVector chunksize, bool display_progress);
RcppExport SEXP RcppEigenH5_calc_af(SEXP h5fileSEXP, SEXP groupnameSEXP, SEXP datanameSEXP, SEXP indexSEXP, SEXP chunksizeSEXP, SEXP display_progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string >::type h5file(h5fileSEXP);
    Rcpp::traits::input_parameter< const std::string >::type groupname(groupnameSEXP);
    Rcpp::traits::input_parameter< const std::string >::type dataname(datanameSEXP);
    Rcpp::traits::input_parameter< const Eigen::ArrayXi >::type index(indexSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector >::type chunksize(chunksizeSEXP);
    Rcpp::traits::input_parameter< bool >::type display_progress(display_progressSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_af(h5file, groupname, dataname, index, chunksize, display_progress));
    return rcpp_result_gen;
END_RCPP
}
// calc_var
Rcpp::NumericVector calc_var(const std::string h5file, const std::string groupname, const std::string dataname, const Eigen::ArrayXi index, const Rcpp::IntegerVector chunksize, bool display_progress);
RcppExport SEXP RcppEigenH5_calc_var(SEXP h5fileSEXP, SEXP groupnameSEXP, SEXP datanameSEXP, SEXP indexSEXP, SEXP chunksizeSEXP, SEXP display_progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string >::type h5file(h5fileSEXP);
    Rcpp::traits::input_parameter< const std::string >::type groupname(groupnameSEXP);
    Rcpp::traits::input_parameter< const std::string >::type dataname(datanameSEXP);
    Rcpp::traits::input_parameter< const Eigen::ArrayXi >::type index(indexSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector >::type chunksize(chunksizeSEXP);
    Rcpp::traits::input_parameter< bool >::type display_progress(display_progressSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_var(h5file, groupname, dataname, index, chunksize, display_progress));
    return rcpp_result_gen;
END_RCPP
}
// calc_yh
Rcpp::NumericVector calc_yh(const std::string h5file, const std::string groupname, const std::string dataname, const Eigen::ArrayXi index, const Eigen::VectorXd beta, const Rcpp::IntegerVector chunksize, bool display_progress);
RcppExport SEXP RcppEigenH5_calc_yh(SEXP h5fileSEXP, SEXP groupnameSEXP, SEXP datanameSEXP, SEXP indexSEXP, SEXP betaSEXP, SEXP chunksizeSEXP, SEXP display_progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string >::type h5file(h5fileSEXP);
    Rcpp::traits::input_parameter< const std::string >::type groupname(groupnameSEXP);
    Rcpp::traits::input_parameter< const std::string >::type dataname(datanameSEXP);
    Rcpp::traits::input_parameter< const Eigen::ArrayXi >::type index(indexSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector >::type chunksize(chunksizeSEXP);
    Rcpp::traits::input_parameter< bool >::type display_progress(display_progressSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_yh(h5file, groupname, dataname, index, beta, chunksize, display_progress));
    return rcpp_result_gen;
END_RCPP
}
// h5ls_grp_exp
std::vector<std::string> h5ls_grp_exp(std::string h5file, std::string base_group);
RcppExport SEXP RcppEigenH5_h5ls_grp_exp(SEXP h5fileSEXP, SEXP base_groupSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type h5file(h5fileSEXP);
    Rcpp::traits::input_parameter< std::string >::type base_group(base_groupSEXP);
    rcpp_result_gen = Rcpp::wrap(h5ls_grp_exp(h5file, base_group));
    return rcpp_result_gen;
END_RCPP
}
// h5ls_attr_exp
std::vector<std::string> h5ls_attr_exp(std::string h5file, std::string base_group);
RcppExport SEXP RcppEigenH5_h5ls_attr_exp(SEXP h5fileSEXP, SEXP base_groupSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type h5file(h5fileSEXP);
    Rcpp::traits::input_parameter< std::string >::type base_group(base_groupSEXP);
    rcpp_result_gen = Rcpp::wrap(h5ls_attr_exp(h5file, base_group));
    return rcpp_result_gen;
END_RCPP
}
// get_h5_version_exp
StringVector get_h5_version_exp();
RcppExport SEXP RcppEigenH5_get_h5_version_exp() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(get_h5_version_exp());
    return rcpp_result_gen;
END_RCPP
}
// read_data_attr_h5_exp
Rcpp::CharacterVector read_data_attr_h5_exp(const StringVector h5filename, const StringVector h5_groupname, const StringVector h5_dataname, const StringVector h5_attr_name);
RcppExport SEXP RcppEigenH5_read_data_attr_h5_exp(SEXP h5filenameSEXP, SEXP h5_groupnameSEXP, SEXP h5_datanameSEXP, SEXP h5_attr_nameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const StringVector >::type h5filename(h5filenameSEXP);
    Rcpp::traits::input_parameter< const StringVector >::type h5_groupname(h5_groupnameSEXP);
    Rcpp::traits::input_parameter< const StringVector >::type h5_dataname(h5_datanameSEXP);
    Rcpp::traits::input_parameter< const StringVector >::type h5_attr_name(h5_attr_nameSEXP);
    rcpp_result_gen = Rcpp::wrap(read_data_attr_h5_exp(h5filename, h5_groupname, h5_dataname, h5_attr_name));
    return rcpp_result_gen;
END_RCPP
}
// read_group_attr_h5_exp
Rcpp::CharacterVector read_group_attr_h5_exp(const std::string h5file, const std::string groupname, const std::string attr_name);
RcppExport SEXP RcppEigenH5_read_group_attr_h5_exp(SEXP h5fileSEXP, SEXP groupnameSEXP, SEXP attr_nameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string >::type h5file(h5fileSEXP);
    Rcpp::traits::input_parameter< const std::string >::type groupname(groupnameSEXP);
    Rcpp::traits::input_parameter< const std::string >::type attr_name(attr_nameSEXP);
    rcpp_result_gen = Rcpp::wrap(read_group_attr_h5_exp(h5file, groupname, attr_name));
    return rcpp_result_gen;
END_RCPP
}
// read_group_iarray_attr_h5_exp
Rcpp::IntegerVector read_group_iarray_attr_h5_exp(const std::string h5file, const std::string groupname, const std::string attr_name);
RcppExport SEXP RcppEigenH5_read_group_iarray_attr_h5_exp(SEXP h5fileSEXP, SEXP groupnameSEXP, SEXP attr_nameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string >::type h5file(h5fileSEXP);
    Rcpp::traits::input_parameter< const std::string >::type groupname(groupnameSEXP);
    Rcpp::traits::input_parameter< const std::string >::type attr_name(attr_nameSEXP);
    rcpp_result_gen = Rcpp::wrap(read_group_iarray_attr_h5_exp(h5file, groupname, attr_name));
    return rcpp_result_gen;
END_RCPP
}
// read_data_iarray_attr_h5_exp
Rcpp::IntegerVector read_data_iarray_attr_h5_exp(const std::string h5file, const std::string groupname, const std::string dataname, const std::string attr_name);
RcppExport SEXP RcppEigenH5_read_data_iarray_attr_h5_exp(SEXP h5fileSEXP, SEXP groupnameSEXP, SEXP datanameSEXP, SEXP attr_nameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string >::type h5file(h5fileSEXP);
    Rcpp::traits::input_parameter< const std::string >::type groupname(groupnameSEXP);
    Rcpp::traits::input_parameter< const std::string >::type dataname(datanameSEXP);
    Rcpp::traits::input_parameter< const std::string >::type attr_name(attr_nameSEXP);
    rcpp_result_gen = Rcpp::wrap(read_data_iarray_attr_h5_exp(h5file, groupname, dataname, attr_name));
    return rcpp_result_gen;
END_RCPP
}
// read_igroup_attr_h5_exp
Rcpp::IntegerVector read_igroup_attr_h5_exp(const std::string h5file, const std::string groupname, const std::string attr_name);
RcppExport SEXP RcppEigenH5_read_igroup_attr_h5_exp(SEXP h5fileSEXP, SEXP groupnameSEXP, SEXP attr_nameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string >::type h5file(h5fileSEXP);
    Rcpp::traits::input_parameter< const std::string >::type groupname(groupnameSEXP);
    Rcpp::traits::input_parameter< const std::string >::type attr_name(attr_nameSEXP);
    rcpp_result_gen = Rcpp::wrap(read_igroup_attr_h5_exp(h5file, groupname, attr_name));
    return rcpp_result_gen;
END_RCPP
}
// read_idata_attr_h5_exp
Rcpp::IntegerVector read_idata_attr_h5_exp(const std::string h5file, const std::string groupname, const std::string dataname, const std::string attr_name);
RcppExport SEXP RcppEigenH5_read_idata_attr_h5_exp(SEXP h5fileSEXP, SEXP groupnameSEXP, SEXP datanameSEXP, SEXP attr_nameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string >::type h5file(h5fileSEXP);
    Rcpp::traits::input_parameter< const std::string >::type groupname(groupnameSEXP);
    Rcpp::traits::input_parameter< const std::string >::type dataname(datanameSEXP);
    Rcpp::traits::input_parameter< const std::string >::type attr_name(attr_nameSEXP);
    rcpp_result_gen = Rcpp::wrap(read_idata_attr_h5_exp(h5file, groupname, dataname, attr_name));
    return rcpp_result_gen;
END_RCPP
}
// write_data_string_attr_h5_exp
void write_data_string_attr_h5_exp(const StringVector h5filename, const StringVector h5_groupname, const StringVector h5_dataname, const StringVector h5_attr_name, const StringVector h5_attr_value);
RcppExport SEXP RcppEigenH5_write_data_string_attr_h5_exp(SEXP h5filenameSEXP, SEXP h5_groupnameSEXP, SEXP h5_datanameSEXP, SEXP h5_attr_nameSEXP, SEXP h5_attr_valueSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const StringVector >::type h5filename(h5filenameSEXP);
    Rcpp::traits::input_parameter< const StringVector >::type h5_groupname(h5_groupnameSEXP);
    Rcpp::traits::input_parameter< const StringVector >::type h5_dataname(h5_datanameSEXP);
    Rcpp::traits::input_parameter< const StringVector >::type h5_attr_name(h5_attr_nameSEXP);
    Rcpp::traits::input_parameter< const StringVector >::type h5_attr_value(h5_attr_valueSEXP);
    write_data_string_attr_h5_exp(h5filename, h5_groupname, h5_dataname, h5_attr_name, h5_attr_value);
    return R_NilValue;
END_RCPP
}
// write_group_string_attr_h5_exp
void write_group_string_attr_h5_exp(const StringVector h5filename, const StringVector h5_groupname, const StringVector h5_attr_name, const StringVector h5_attr_value);
RcppExport SEXP RcppEigenH5_write_group_string_attr_h5_exp(SEXP h5filenameSEXP, SEXP h5_groupnameSEXP, SEXP h5_attr_nameSEXP, SEXP h5_attr_valueSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const StringVector >::type h5filename(h5filenameSEXP);
    Rcpp::traits::input_parameter< const StringVector >::type h5_groupname(h5_groupnameSEXP);
    Rcpp::traits::input_parameter< const StringVector >::type h5_attr_name(h5_attr_nameSEXP);
    Rcpp::traits::input_parameter< const StringVector >::type h5_attr_value(h5_attr_valueSEXP);
    write_group_string_attr_h5_exp(h5filename, h5_groupname, h5_attr_name, h5_attr_value);
    return R_NilValue;
END_RCPP
}
// write_data_int_attr_h5_exp
void write_data_int_attr_h5_exp(const StringVector h5filename, const StringVector h5_groupname, const StringVector h5_dataname, const StringVector h5_attr_name, const IntegerVector h5_attr_value);
RcppExport SEXP RcppEigenH5_write_data_int_attr_h5_exp(SEXP h5filenameSEXP, SEXP h5_groupnameSEXP, SEXP h5_datanameSEXP, SEXP h5_attr_nameSEXP, SEXP h5_attr_valueSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const StringVector >::type h5filename(h5filenameSEXP);
    Rcpp::traits::input_parameter< const StringVector >::type h5_groupname(h5_groupnameSEXP);
    Rcpp::traits::input_parameter< const StringVector >::type h5_dataname(h5_datanameSEXP);
    Rcpp::traits::input_parameter< const StringVector >::type h5_attr_name(h5_attr_nameSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type h5_attr_value(h5_attr_valueSEXP);
    write_data_int_attr_h5_exp(h5filename, h5_groupname, h5_dataname, h5_attr_name, h5_attr_value);
    return R_NilValue;
END_RCPP
}
// write_group_int_attr_h5_exp
void write_group_int_attr_h5_exp(const StringVector h5filename, const StringVector h5_groupname, const StringVector h5_attr_name, const IntegerVector h5_attr_value);
RcppExport SEXP RcppEigenH5_write_group_int_attr_h5_exp(SEXP h5filenameSEXP, SEXP h5_groupnameSEXP, SEXP h5_attr_nameSEXP, SEXP h5_attr_valueSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const StringVector >::type h5filename(h5filenameSEXP);
    Rcpp::traits::input_parameter< const StringVector >::type h5_groupname(h5_groupnameSEXP);
    Rcpp::traits::input_parameter< const StringVector >::type h5_attr_name(h5_attr_nameSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type h5_attr_value(h5_attr_valueSEXP);
    write_group_int_attr_h5_exp(h5filename, h5_groupname, h5_attr_name, h5_attr_value);
    return R_NilValue;
END_RCPP
}
// write_mat_chunk_h5_exp
void write_mat_chunk_h5_exp(const StringVector h5file, const StringVector groupname, const StringVector dataname, const Matrix_external data, const IntegerVector offsets);
RcppExport SEXP RcppEigenH5_write_mat_chunk_h5_exp(SEXP h5fileSEXP, SEXP groupnameSEXP, SEXP datanameSEXP, SEXP dataSEXP, SEXP offsetsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const StringVector >::type h5file(h5fileSEXP);
    Rcpp::traits::input_parameter< const StringVector >::type groupname(groupnameSEXP);
    Rcpp::traits::input_parameter< const StringVector >::type dataname(datanameSEXP);
    Rcpp::traits::input_parameter< const Matrix_external >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type offsets(offsetsSEXP);
    write_mat_chunk_h5_exp(h5file, groupname, dataname, data, offsets);
    return R_NilValue;
END_RCPP
}
// create_mat_dataset_h5_exp
void create_mat_dataset_h5_exp(const StringVector h5file, const StringVector groupname, const StringVector dataname, const IntegerVector dims, const IntegerVector chunkdims, const IntegerVector deflate_level, const bool doTranspose);
RcppExport SEXP RcppEigenH5_create_mat_dataset_h5_exp(SEXP h5fileSEXP, SEXP groupnameSEXP, SEXP datanameSEXP, SEXP dimsSEXP, SEXP chunkdimsSEXP, SEXP deflate_levelSEXP, SEXP doTransposeSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const StringVector >::type h5file(h5fileSEXP);
    Rcpp::traits::input_parameter< const StringVector >::type groupname(groupnameSEXP);
    Rcpp::traits::input_parameter< const StringVector >::type dataname(datanameSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type dims(dimsSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type chunkdims(chunkdimsSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type deflate_level(deflate_levelSEXP);
    Rcpp::traits::input_parameter< const bool >::type doTranspose(doTransposeSEXP);
    create_mat_dataset_h5_exp(h5file, groupname, dataname, dims, chunkdims, deflate_level, doTranspose);
    return R_NilValue;
END_RCPP
}
// write_mat_h5_exp
void write_mat_h5_exp(const StringVector h5file, const StringVector groupname, const StringVector dataname, const Matrix_external data, const IntegerVector deflate_level, const bool doTranspose);
RcppExport SEXP RcppEigenH5_write_mat_h5_exp(SEXP h5fileSEXP, SEXP groupnameSEXP, SEXP datanameSEXP, SEXP dataSEXP, SEXP deflate_levelSEXP, SEXP doTransposeSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const StringVector >::type h5file(h5fileSEXP);
    Rcpp::traits::input_parameter< const StringVector >::type groupname(groupnameSEXP);
    Rcpp::traits::input_parameter< const StringVector >::type dataname(datanameSEXP);
    Rcpp::traits::input_parameter< const Matrix_external >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type deflate_level(deflate_levelSEXP);
    Rcpp::traits::input_parameter< const bool >::type doTranspose(doTransposeSEXP);
    write_mat_h5_exp(h5file, groupname, dataname, data, deflate_level, doTranspose);
    return R_NilValue;
END_RCPP
}
// write_dvec_h5_exp
void write_dvec_h5_exp(const StringVector h5file, const StringVector groupname, const StringVector dataname, const arrayxd_external data, const IntegerVector deflate_level);
RcppExport SEXP RcppEigenH5_write_dvec_h5_exp(SEXP h5fileSEXP, SEXP groupnameSEXP, SEXP datanameSEXP, SEXP dataSEXP, SEXP deflate_levelSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const StringVector >::type h5file(h5fileSEXP);
    Rcpp::traits::input_parameter< const StringVector >::type groupname(groupnameSEXP);
    Rcpp::traits::input_parameter< const StringVector >::type dataname(datanameSEXP);
    Rcpp::traits::input_parameter< const arrayxd_external >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type deflate_level(deflate_levelSEXP);
    write_dvec_h5_exp(h5file, groupname, dataname, data, deflate_level);
    return R_NilValue;
END_RCPP
}
// write_ivec_h5_exp
void write_ivec_h5_exp(const StringVector h5file, const StringVector groupname, const StringVector dataname, const arrayxi_external data, const IntegerVector deflate_level);
RcppExport SEXP RcppEigenH5_write_ivec_h5_exp(SEXP h5fileSEXP, SEXP groupnameSEXP, SEXP datanameSEXP, SEXP dataSEXP, SEXP deflate_levelSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const StringVector >::type h5file(h5fileSEXP);
    Rcpp::traits::input_parameter< const StringVector >::type groupname(groupnameSEXP);
    Rcpp::traits::input_parameter< const StringVector >::type dataname(datanameSEXP);
    Rcpp::traits::input_parameter< const arrayxi_external >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type deflate_level(deflate_levelSEXP);
    write_ivec_h5_exp(h5file, groupname, dataname, data, deflate_level);
    return R_NilValue;
END_RCPP
}
