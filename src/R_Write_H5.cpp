#include "RcppEigenH5.h"


//[[Rcpp::depends(RcppEigen)]]
using namespace Eigen;

//[[Rcpp::export(name="write_mat_chunk_h5")]]
void write_mat_chunk_h5_exp(const StringVector h5file, const StringVector groupname, const StringVector dataname,const Matrix_external data,const IntegerVector offsets){

  Eigen::ArrayXi offset=Rcpp::as<Eigen::ArrayXi>(offsets);
  std::string h5fn(h5file[0]);
  std::string gname(groupname[0]);
  std::string dname(dataname[0]);
  write_mat_chunk_h5(h5fn,gname,dname,data,offset);
}





//[[Rcpp::export(name="create_mat_dataset_h5")]]
void create_mat_dataset_h5_exp(const StringVector h5file, const StringVector groupname, const StringVector dataname,const IntegerVector dims, const IntegerVector chunkdims, const IntegerVector deflate_level=IntegerVector::create(0),const bool doTranspose=false){

  std::string filename(h5file[0]);
  std::string gname(groupname[0]);
  std::string dname(dataname[0]);


  int d_level =deflate_level[0];





  size_t rownum =dims[0];
  size_t colnum=dims[1];

  size_t rowchunk=chunkdims[0];
  size_t colchunk=chunkdims[1];

  if(doTranspose){
    std::swap(rownum,colnum);
    std::swap(rowchunk,colchunk);

  }
  H5FilePtr file =create_or_open_file(filename);

  H5GroupPtr group =create_or_open_group(file,gname);
  FloatType ftypew(PredType::NATIVE_DOUBLE);

  std::vector<hsize_t> cumdim{rownum,colnum};
  std::vector<hsize_t> maxdim{rownum,colnum};
  std::vector<hsize_t> chunkdim{rowchunk,colchunk};

  H5DataSetPtr dataset =create_or_open_dataset(group,dname,ftypew,cumdim,maxdim,chunkdim,d_level);
  write_transpose(dataset,doTranspose);
  dataset->close();
  group->close();
  file->close();

}


// //[[Rcpp::export]]
// Eigen::MatrixXd test_rowcol(const Matrix_external data){
//   Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>
//  }




//[[Rcpp::export(name="write_mat_h5")]]
void write_mat_h5_exp(const StringVector h5file, const StringVector groupname, const StringVector dataname,const Matrix_external data, const IntegerVector deflate_level=IntegerVector::create(0),const bool doTranspose=false){

  std::string th5file(h5file[0]);
  std::string tgroupname(groupname[0]);
  std::string tdataname(dataname[0]);
  int deflate=deflate_level[0];

  Eigen::ArrayXi chunksize(2);
  chunksize[0]=data.rows();
  chunksize[1]=std::min(10000,(int)data.cols());
//
//   if(bdoTranspose){
//     std::swap(chunksize[0],chunksize[1]);
//   }

  write_mat_h5(th5file,tgroupname,tdataname,data,chunksize,deflate,doTranspose);
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



