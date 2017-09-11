#include "RcppEigenH5.h"
// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
#include <progress.hpp>

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

  size_t groupnum=groupname.size();
  std::vector<std::string> groupnames(groupnum);
  for(size_t i=0;i<groupnum;i++){
    groupnames[i]=groupname[i];
  }




//  std::string gname(groupname[0]);
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
  FloatType ftypew(PredType::NATIVE_DOUBLE);

  std::vector<hsize_t> cumdim{rownum,colnum};
  std::vector<hsize_t> maxdim{rownum,colnum};
  std::vector<hsize_t> chunkdim{rowchunk,colchunk};
  for(size_t i=0; i<groupnum;i++){
    std::string gname=groupnames[i];
    H5GroupPtr group =create_or_open_group(file,gname);

    H5DataSetPtr dataset =create_or_open_dataset(group,dname,ftypew,cumdim,maxdim,chunkdim,d_level);
    write_transpose(dataset,doTranspose);
    dataset->close();
    group->close();
  }
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





//[[Rcpp::export(name="write_svec_h5")]]
void write_svec_h5_exp(const StringVector h5file, const StringVector groupname, const StringVector dataname,const StringVector data, const IntegerVector deflate_level=IntegerVector::create(0)){

  // using namespace HighFive
  std::string th5file(h5file[0]);
  std::string tgroupname(groupname[0]);
  std::string tdataname(dataname[0]);
  std::vector<std::string> sdata=Rcpp::as<std::vector<std::string> >(data);
  const bool hfile_exists =f_exists(th5file);

  auto writers_type = hfile_exists ? HighFive::File::ReadWrite : HighFive::File::ReadWrite |  HighFive::File::Create;

  HighFive::File file(th5file,writers_type);


  if(!file.exist(tgroupname)){
    file.createGroup(tgroupname);
  }

  // HighFive::Group g =file.open_group(tgroupname);

  HighFive::DataSet dset =file.createDataSet<std::string>(tgroupname+"/"+tdataname,HighFive::DataSpace::From(sdata));
  dset.write(sdata);

}

//
//   // std::vector<std::string> read_svec_h5(const std::string h5file, const std::string groupname, const std::string dataname){
//
//     size_t offset=0;
//
//     char** datasa;
//
//     H5FilePtr file =create_or_open_file(h5file);
//
//
//     size_t chunksize=1000;
//     if(chunksize>datasize){
//       chunksize=datasize;
//     }
//
//     H5GroupPtr group =create_or_open_group(file,groupname);
//     std::vector<hsize_t> cumdim{datasize};
//     std::vector<hsize_t> maxdim{datasize};
//     std::vector<hsize_t> chunkdim{chunksize};
//
//     H5DataSetPtr dataset =create_or_open_dataset(group,dataname,ftypew,cumdim,maxdim,chunkdim,deflate_level);
//
//     DataSpace* fdataspace;
//     H5FilePtr file=create_or_open_file(Rcpp::as<std::string>(h5file[0]));
//     H5GroupPtr group= create_or_open_group(file,Rcpp::as<std::string>(groupname[0]));
//     H5DataSetPtr dataset = create_or_open_dataset(group,Rcpp::as<std::string>(dataname[0]));
//     DSetCreatPropList cparms= dataset->getCreatePlist();
//
//
//     DataType dt= dataset->getDataType();
//     hsize_t datadims[]={0};
//     DataSpace fspace=dataset->getSpace();
//     fspace.getSimpleExtentDims(datadims,NULL);
//     //  std::cout<<"Full data is of dimensions"<<datadims[0]<<"x"<<datadims[1]<<std::endl;
//     hsize_t vec_dims[1];
//     size_t  chunksize =datadims[0];
//     if(offset+chunksize>datadims[0]){
//       Rcpp::stop("offset+chunksize>datadims[0]!");
//     }
//
//     datasa = new char*[chunksize];
//     vec_dims[0]=chunksize;
//     datadims[0]=chunksize;
//
//     // std::cout<<"Full data is of size "<<datadim[0]<<std::endl;
//     hsize_t offseta[1];
//     offseta[0]=offset;
//
//     fspace.selectHyperslab(H5S_SELECT_SET,datadims,offseta);
//     //  std::cout<<"Allocating matrix of size:"<<matrix_dims[0]<<"x"<<matrix_dims[1]<<std::endl;
//     //  std::cout<<"Matrix starts at"<<row_offset<<"x"<<col_offset<<std::endl;
//     DataSpace memspace(1,vec_dims);
//     //  std::cout<<"Reading data"<<std::endl;
//     dataset->read((void*)datasa,dt,memspace,fspace);
//     //  std::cout<<"Read complete!"<<std::endl;
//     dt.close();
//     memspace.close();
//     fspace.close();
//     group->close();
//     dataset->close();
//     file->close();
//     std::vector<std::string> retvec(chunksize);
//     for(int i=0;i<chunksize;i++){
//       retvec[i]=std::string(datasa[i]);
//     }
//
//
//     return(retvec);
//   }



