#include <RcppEigenH5.h>
//[[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;



//[[Rcpp::export(name="read_2d_h5")]]
Eigen::MatrixXd read_2d_h5_exp(const std::string h5file, const std::string groupname, const std::string dataname,const  IntegerVector offset ,const  IntegerVector chunksize){

  return(read_2d_h5(h5file,groupname,dataname,as<Eigen::ArrayXi>(offset),as<Eigen::ArrayXi>(chunksize)));
}
//[[Rcpp::export]]
Eigen::ArrayXd read_1d_h5(const std::string h5file, const std::string groupname, const std::string dataname,const  IntegerVector offset ,const  IntegerVector chunksize){

  size_t v_offset=offset[0];
  size_t v_chunksize=chunksize[0];
  return(read_dvec_h5(h5file,groupname,dataname,v_offset,v_chunksize));
}

//[[Rcpp::export]]
Eigen::ArrayXi read_1i_h5(const std::string h5file, const std::string groupname, const std::string dataname,const  IntegerVector offset ,const  IntegerVector chunksize){



  size_t v_offset=offset[0];
  size_t v_chunksize=chunksize[0];
  Eigen::ArrayXi retvec(v_chunksize);
  read_ivec_h5(h5file,groupname,dataname,v_offset,v_chunksize,retvec.data());
  return(retvec);
}

//[[Rcpp::export]]
Eigen::MatrixXd read_2d_mat_h5(const std::string h5file, const std::string groupname, const std::string dataname){

  Rcpp::IntegerVector offset(2);
  offset[0]=0;
  offset[1]=0;
  Rcpp::IntegerVector chunksize(2);


  int rownum=get_rownum_h5(h5file,groupname,dataname);
  int colnum=get_colnum_h5(h5file,groupname,dataname);
  chunksize[0]=rownum;
  chunksize[1]=colnum;
  return(read_2d_h5_exp(h5file,groupname,dataname,offset,chunksize));
}


//[[Rcpp::export]]
Eigen::MatrixXd read_2d_index_h5(const std::string h5file,const std::string groupname, const std::string dataname, const  IntegerVector indvec){
  Eigen::ArrayXi ivec(as<Eigen::ArrayXi>(indvec));

  return(read_2d_cindex_h5(h5file,groupname,dataname,ivec));
}

//[[Rcpp::export]]
Eigen::ArrayXd read_1d_index_h5(const std::string h5file,const std::string groupname, const std::string dataname, const  IntegerVector indvec){

  Eigen::ArrayXi ivec(as<Eigen::ArrayXi>(indvec));
  return(read_2d_cindex_h5(h5file,groupname,dataname,ivec));

}
