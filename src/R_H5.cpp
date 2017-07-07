#include <RcppEigenH5.h>
//[[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;



//[[Rcpp::export(name="read_2d_h5")]]
Eigen::MatrixXd read_2d_h5_exp(const StringVector h5file, const StringVector groupname, const StringVector dataname,const  IntegerVector offset ,const  IntegerVector chunksize){

  std::string th5file(h5file[0]);
  std::string tgroupname(groupname[0]);
  std::string tdataname(dataname[0]);

  return(read_2d_h5(th5file,tgroupname,tdataname,as<Eigen::ArrayXi>(offset),as<Eigen::ArrayXi>(chunksize)));
}


//[[Rcpp::export]]
Eigen::ArrayXd read_1d_h5(const StringVector h5file, const StringVector groupname, const StringVector dataname,const  IntegerVector offset ,const  IntegerVector chunksize){

  std::string th5file(h5file[0]);
  std::string tgroupname(groupname[0]);
  std::string tdataname(dataname[0]);
  size_t v_offset=offset[0];
  size_t v_chunksize=chunksize[0];
  return(read_dvec_h5(th5file,tgroupname,tdataname,v_offset,v_chunksize));
}

//[[Rcpp::export(name="read_svec")]]
Rcpp::StringVector read_svec_exp(const StringVector h5file, const StringVector groupname, const StringVector dataname){
  std::string th5file(h5file[0]);
  std::string tgroupname(groupname[0]);
  std::string tdataname(dataname[0]);
  return(Rcpp::wrap(read_svec_h5(th5file,tgroupname,tdataname)));
}

//[[Rcpp::export(name="read_svec_chunk")]]
Rcpp::StringVector read_svec_chunk_exp(const StringVector h5file, const StringVector groupname, const StringVector dataname,const int offset,const int chunksize){
  std::string th5file(h5file[0]);
  std::string tgroupname(groupname[0]);
  std::string tdataname(dataname[0]);
  return(Rcpp::wrap(read_svec_h5(th5file,tgroupname,tdataname,offset,chunksize)));
}



//[[Rcpp::export]]
Rcpp::StringVector read_1d_sindex_h5(const StringVector h5file,const StringVector groupname, const StringVector dataname, const  IntegerVector indvec){

  std::string th5file(h5file[0]);
  std::string tgroupname(groupname[0]);

  std::string tdataname(dataname[0]);


  Eigen::ArrayXi ivec(as<Eigen::ArrayXi>(indvec));
  return(Rcpp::wrap(read_s1d_index_h5(th5file,tgroupname,tdataname,ivec)));

}

//[[Rcpp::export(name="read_dvec")]]
Eigen::ArrayXd read_dvec_exp(const StringVector h5file, const StringVector groupname, const StringVector dataname){

  std::string th5file(h5file[0]);
  std::string tgroupname(groupname[0]);
  std::string tdataname(dataname[0]);
  return(read_dvec_h5(th5file,tgroupname,tdataname));
}

//[[Rcpp::export(name="read_ivec")]]
Eigen::ArrayXi read_ivec_exp(const StringVector h5file, const StringVector groupname, const StringVector dataname){

  std::string th5file(h5file[0]);
  std::string tgroupname(groupname[0]);
  std::string tdataname(dataname[0]);
  return(read_ivec_h5(th5file,tgroupname,tdataname));
}


//[[Rcpp::export(name="read_uivec")]]
Eigen::ArrayXi read_uivec_exp(const StringVector h5file, const StringVector groupname, const StringVector dataname){

  std::string th5file(h5file[0]);
  std::string tgroupname(groupname[0]);
  std::string tdataname(dataname[0]);
  return(read_uivec_h5(th5file,tgroupname,tdataname));
}



//[[Rcpp::export]]
Eigen::ArrayXi read_1i_h5(const StringVector h5file, const StringVector groupname, const StringVector dataname,const  IntegerVector offset ,const  IntegerVector chunksize){

  std::string th5file(h5file[0]);
  std::string tgroupname(groupname[0]);
  std::string tdataname(dataname[0]);

  size_t v_offset=offset[0];
  size_t v_chunksize=chunksize[0];
  Eigen::ArrayXi retvec(v_chunksize);
  read_ivec_h5(th5file,tgroupname,tdataname,v_offset,v_chunksize,retvec.data());
  return(retvec);
}

//[[Rcpp::export]]
Eigen::MatrixXd read_2d_mat_h5(const StringVector h5file, const StringVector groupname, const StringVector dataname){

  std::string th5file(h5file[0]);
  std::string tgroupname(groupname[0]);
  std::string tdataname(dataname[0]);

  Rcpp::IntegerVector offset(2);
  offset[0]=0;
  offset[1]=0;
  Rcpp::IntegerVector chunksize(2);


  int rownum=get_rownum_h5(th5file,tgroupname,tdataname);
  int colnum=get_colnum_h5(th5file,tgroupname,tdataname);
  chunksize[0]=rownum;
  chunksize[1]=colnum;
  return(read_2d_h5_exp(th5file,tgroupname,tdataname,offset,chunksize));
}


//[[Rcpp::export]]
Eigen::MatrixXd read_2d_index_h5(const StringVector h5file,const StringVector groupname, const StringVector dataname,const IntegerVector indvec){

  std::string th5file(h5file[0]);
  std::string tgroupname(groupname[0]);
  std::string tdataname(dataname[0]);

  Eigen::ArrayXi ivec(indvec.size());
  for(int i=0; i<indvec.size(); i++){
    ivec.coeffRef(i)=indvec[i];
  }

  //  Eigen::ArrayXi ivec(as<Eigen::ArrayXi>(indvec));
  const  int i_colnum(indvec.size());
  const  int i_rownum (get_rownum_h5(th5file,tgroupname,tdataname));

  //  Rcpp::Rcout<<"retmat should be "<<i_rownum<<"x"<<i_colnum<<std::endl;
  Eigen::MatrixXd retmat(i_rownum,i_colnum);

  retmat.setZero();
  //  Rcpp::Rcout<<"";
  //   Rcpp::Rcout<<"retmat is "<<retmat.rows()<<"x"<<retmat.cols()<<std::endl;
  read_2d_cindex_h5(th5file,tgroupname,tdataname,ivec,retmat);

  return(retmat);
}


//[[Rcpp::export]]
Eigen::MatrixXd read_2d_index_chunk_h5(const StringVector h5file,const StringVector groupname, const StringVector dataname,const IntegerVector indvec,const IntegerVector chunksize){

  size_t csize=chunksize[0];
  std::string th5file(h5file[0]);
  std::string tgroupname(groupname[0]);
  std::string tdataname(dataname[0]);

  Eigen::ArrayXi ivec(indvec.size());
  for(int i=0; i<indvec.size(); i++){
    ivec.coeffRef(i)=indvec[i];
  }

  //  Eigen::ArrayXi ivec(as<Eigen::ArrayXi>(indvec));
  const  int i_colnum(indvec.size());
  const  int i_rownum (get_rownum_h5(th5file,tgroupname,tdataname));

  //  Rcpp::Rcout<<"retmat should be "<<i_rownum<<"x"<<i_colnum<<std::endl;
  Eigen::MatrixXd retmat(i_rownum,i_colnum);

  retmat.setZero();
  //  Rcpp::Rcout<<"";
  //   Rcpp::Rcout<<"retmat is "<<retmat.rows()<<"x"<<retmat.cols()<<std::endl;
  read_2d_cindex_chunk_h5(th5file,tgroupname,tdataname,ivec,retmat,csize);

  return(retmat);
}



//[[Rcpp::export]]
Eigen::ArrayXd read_1d_index_h5(const StringVector h5file,const StringVector groupname, const StringVector dataname, const  IntegerVector indvec){

  std::string th5file(h5file[0]);
  std::string tgroupname(groupname[0]);

  std::string tdataname(dataname[0]);


  Eigen::ArrayXi ivec(as<Eigen::ArrayXi>(indvec));
  return(read_1d_cindex_h5(th5file,tgroupname,tdataname,ivec));

}




