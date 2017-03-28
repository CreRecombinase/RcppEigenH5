#include "RcppEigenH5.h"


//[[Rcpp::depends(RcppEigen)]]
using namespace Eigen;


//[[Rcpp::export(name="write_mat_h5")]]
void write_mat_h5_exp(const StringVector h5file, const StringVector groupname, const StringVector dataname,const Matrix_external data, const IntegerVector deflate_level=IntegerVector::create(0)){

  std::string th5file(h5file[0]);
  std::string tgroupname(groupname[0]);
  std::string tdataname(dataname[0]);
  int deflate=deflate_level[0];
  write_mat_h5(th5file,tgroupname,tdataname,data,deflate);
}

//[[Rcpp::export(name="write_dvec_h5")]]
void write_dvec_h5_exp(const StringVector h5file, const StringVector groupname, const StringVector dataname,const arrayxd_external data, const IntegerVector deflate_level=IntegerVector::create(0)){

  std::string th5file(h5file[0]);
  std::string tgroupname(groupname[0]);
  std::string tdataname(dataname[0]);
  int deflate=deflate_level[0];
  write_dvec_h5(th5file,tgroupname,tdataname,data,deflate);
}


//[[Rcpp::export(name="write_ivec_h5")]]
void write_ivec_h5_exp(const StringVector h5file, const StringVector groupname, const StringVector dataname,const arrayxi_external data, const IntegerVector deflate_level=IntegerVector::create(0)){

  std::string th5file(h5file[0]);
  std::string tgroupname(groupname[0]);
  std::string tdataname(dataname[0]);
  int deflate=deflate_level[0];
  write_ivec_h5(th5file,tgroupname,tdataname,data,deflate);
}



