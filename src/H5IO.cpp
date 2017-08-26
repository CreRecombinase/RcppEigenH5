#include "RcppEigenH5.h"
#include <type_traits>
//[[Rcpp::depends(RcppEigen)]]
using namespace Eigen;







void read_2dfmat_h5(const std::string h5file, const std::string groupname, const std::string dataname, int row_offset, int col_offset, int row_chunksize, int col_chunksize,  float* data){

  //Try breaking up reads in to chunks
  H5FilePtr file=open_file(h5file);
  H5GroupPtr group= open_group(file,groupname);
  H5DataSetPtr dataset = open_dataset(group,dataname);
  bool isTranspose=check_transpose(dataset);
  if(isTranspose){
    // Rcpp::Rcout<<"row_chunksize is :"<<row_chunksize<<std::endl;
    // Rcpp::Rcout<<"col_chunksize is :"<<col_chunksize<<std::endl;
    std::swap(row_offset,col_offset);
    std::swap(row_chunksize,col_chunksize);
    // Rcpp::Rcout<<"row_chunksize is :"<<row_chunksize<<std::endl;
    // Rcpp::Rcout<<"col_chunksize is :"<<col_chunksize<<std::endl;
  }

//  DataType dt= dataset->getDataType();
FloatType dt(PredType::NATIVE_FLOAT);
  hsize_t datadims[]={0,0};
  DataSpace fspace=dataset->getSpace();
  fspace.getSimpleExtentDims(datadims,NULL);
  //  std::cout<<"Full data is of dimensions"<<datadims[0]<<"x"<<datadims[1]<<std::endl;
  hsize_t matrix_dims[2];
  if(row_offset+row_chunksize>datadims[0]){
    Rcpp::stop("row_offset ("+std::to_string(static_cast<long long>(row_offset))+") + row_chunksize ("+std::to_string(static_cast<long long>(row_chunksize))+") >datadims[0] ("+std::to_string(static_cast<long long>(datadims[0]))+")");
  }
  if(col_offset+col_chunksize>datadims[1]){
    Rcpp::stop("col_offset ("+std::to_string(static_cast<long long>(col_offset))+") + col_chunksize ("+std::to_string(static_cast<long long>(col_chunksize))+") >datadims[1] ("+std::to_string(static_cast<long long>(datadims[1]))+")");
  }

  // std::cout<<"read consists of"<<colchunknum*rowchunknum<<"chunks"<<std::endl;
  // long slab_bytes= rchunksize*cchunksize*sizeof(float);
  // std::cout<<"Slab is "<<slab_bytes<<"bytes"<<std::endl;
  // plist.setCache(nmdc,cache_elem,slab_bytes,1);



  matrix_dims[0]=row_chunksize;
  matrix_dims[1]=col_chunksize;
  datadims[0]=row_chunksize;
  datadims[1]=col_chunksize;

  // std::cout<<"Full data is of size "<<datadim[0]<<std::endl;
  hsize_t offseta[2];
  offseta[0]=row_offset;
  offseta[1]=col_offset;

  fspace.selectHyperslab(H5S_SELECT_SET,datadims,offseta);
  //  std::cout<<"Allocating matrix of size:"<<matrix_dims[0]<<"x"<<matrix_dims[1]<<std::endl;
  //  std::cout<<"Matrix starts at"<<row_offset<<"x"<<col_offset<<std::endl;
  DataSpace memspace(2,matrix_dims);
  //  std::cout<<"Reading data"<<std::endl;
  if(isTranspose){

    Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> readmat(row_chunksize,col_chunksize);
    dataset->read(readmat.data(),dt,memspace,fspace);

    // Rcout<<"readmat was "<<print_dims(readmat)<<std::endl;
    readmat=readmat.transpose().eval();
    // Rcout<<"readmat is "<<print_dims(readmat)<<std::endl;
    std::copy(readmat.data(),readmat.data()+readmat.size(),data);
  }else{
    dataset->read(data,dt,memspace,fspace);
  }
  //  std::cout<<"Read complete!"<<std::endl;

  dt.close();
  memspace.close();
  fspace.close();
  group->close();
  dataset->close();
  file->close();
}






void read_2ddmat_h5(const std::string h5file, const std::string groupname, const std::string dataname, int row_offset, int col_offset, int row_chunksize, int col_chunksize,  double* data){

  //Try breaking up reads in to chunks
  H5FilePtr file=open_file(h5file);
  H5GroupPtr group= open_group(file,groupname);
  H5DataSetPtr dataset = open_dataset(group,dataname);
  bool isTranspose=check_transpose(dataset);
  if(isTranspose){
    // Rcpp::Rcout<<"row_chunksize is :"<<row_chunksize<<std::endl;
    // Rcpp::Rcout<<"col_chunksize is :"<<col_chunksize<<std::endl;
    std::swap(row_offset,col_offset);
    std::swap(row_chunksize,col_chunksize);
    // Rcpp::Rcout<<"row_chunksize is :"<<row_chunksize<<std::endl;
    // Rcpp::Rcout<<"col_chunksize is :"<<col_chunksize<<std::endl;
  }

  DataType dt= dataset->getDataType();
  hsize_t datadims[]={0,0};
  DataSpace fspace=dataset->getSpace();
  fspace.getSimpleExtentDims(datadims,NULL);
  //  std::cout<<"Full data is of dimensions"<<datadims[0]<<"x"<<datadims[1]<<std::endl;
  hsize_t matrix_dims[2];
  if(row_offset+row_chunksize>datadims[0]){
    Rcpp::stop("row_offset ("+std::to_string(static_cast<long long>(row_offset))+") + row_chunksize ("+std::to_string(static_cast<long long>(row_chunksize))+") >datadims[0] ("+std::to_string(static_cast<long long>(datadims[0]))+")");
  }
  if(col_offset+col_chunksize>datadims[1]){
    Rcpp::stop("col_offset ("+std::to_string(static_cast<long long>(col_offset))+") + col_chunksize ("+std::to_string(static_cast<long long>(col_chunksize))+") >datadims[1] ("+std::to_string(static_cast<long long>(datadims[1]))+")");
  }

  // std::cout<<"read consists of"<<colchunknum*rowchunknum<<"chunks"<<std::endl;
  // long slab_bytes= rchunksize*cchunksize*sizeof(float);
  // std::cout<<"Slab is "<<slab_bytes<<"bytes"<<std::endl;
  // plist.setCache(nmdc,cache_elem,slab_bytes,1);



  matrix_dims[0]=row_chunksize;
  matrix_dims[1]=col_chunksize;
  datadims[0]=row_chunksize;
  datadims[1]=col_chunksize;

  // std::cout<<"Full data is of size "<<datadim[0]<<std::endl;
  hsize_t offseta[2];
  offseta[0]=row_offset;
  offseta[1]=col_offset;

  fspace.selectHyperslab(H5S_SELECT_SET,datadims,offseta);
  //  std::cout<<"Allocating matrix of size:"<<matrix_dims[0]<<"x"<<matrix_dims[1]<<std::endl;
  //  std::cout<<"Matrix starts at"<<row_offset<<"x"<<col_offset<<std::endl;
  DataSpace memspace(2,matrix_dims);
  //  std::cout<<"Reading data"<<std::endl;
  if(isTranspose){

    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> readmat(row_chunksize,col_chunksize);
    dataset->read(readmat.data(),dt,memspace,fspace);

    // Rcout<<"readmat was "<<print_dims(readmat)<<std::endl;
    readmat=readmat.transpose().eval();
    // Rcout<<"readmat is "<<print_dims(readmat)<<std::endl;
    std::copy(readmat.data(),readmat.data()+readmat.size(),data);
  }else{
    dataset->read(data,dt,memspace,fspace);
  }
  //  std::cout<<"Read complete!"<<std::endl;

  dt.close();
  memspace.close();
  fspace.close();
  group->close();
  dataset->close();
  file->close();


}





Eigen::MatrixXd read_2d_h5(const std::string h5file, const std::string groupname, const std::string dataname,const  Eigen::ArrayXi offset ,const  Eigen::ArrayXi chunksize){

  int row_offset=offset[0];
  int col_offset=offset[1];

  int row_chunksize=chunksize[0];
  int col_chunksize=chunksize[1];

  int tot_size = row_chunksize*col_chunksize;
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> readmat(row_chunksize,col_chunksize);
  read_2ddmat_h5(h5file,groupname,dataname,row_offset,col_offset,row_chunksize,col_chunksize,readmat.data());

  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> retmat=readmat;
  return(retmat);
}

Eigen::MatrixXd read_2d_h5(const std::string h5file, const std::string groupname, const std::string dataname){

  int row_offset=0;
  int col_offset=0;

  int row_chunksize=get_rownum_h5(h5file,groupname,dataname);
  int col_chunksize=get_colnum_h5(h5file,groupname,dataname);

  int tot_size = row_chunksize*col_chunksize;
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> readmat(row_chunksize,col_chunksize);
  read_2ddmat_h5(h5file,groupname,dataname,row_offset,col_offset,row_chunksize,col_chunksize,readmat.data());
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> retmat=readmat;
  return(retmat);
}




Eigen::ArrayXd read_dvec_h5(const std::string h5file, const std::string groupname, const std::string dataname,const int offset,const int chunksize){

  H5FilePtr file=open_file(h5file);
  H5GroupPtr group= open_group(file,groupname);
  H5DataSetPtr dataset = open_dataset(group,dataname);
  DSetCreatPropList cparms= dataset->getCreatePlist();


  DataType dt= dataset->getDataType();
  hsize_t datadims[]={0};
  DataSpace fspace=dataset->getSpace();
  fspace.getSimpleExtentDims(datadims,NULL);
  //  std::cout<<"Full data is of dimensions"<<datadims[0]<<"x"<<datadims[1]<<std::endl;
  hsize_t vec_dims[1];
  if(offset+chunksize>datadims[0]){
    Rcpp::stop("offset+chunksize>datadims[0]!");
  }
  vec_dims[0]=chunksize;
  datadims[0]=chunksize;

  // std::cout<<"Full data is of size "<<datadim[0]<<std::endl;
  hsize_t offseta[1];
  offseta[0]=offset;
  Eigen::ArrayXd dataxd(chunksize);
  fspace.selectHyperslab(H5S_SELECT_SET,datadims,offseta);
  //  std::cout<<"Allocating matrix of size:"<<matrix_dims[0]<<"x"<<matrix_dims[1]<<std::endl;
  //  std::cout<<"Matrix starts at"<<row_offset<<"x"<<col_offset<<std::endl;
  DataSpace memspace(1,vec_dims);
  //  std::cout<<"Reading data"<<std::endl;
  dataset->read(dataxd.data(),dt,memspace,fspace);
  //  std::cout<<"Read complete!"<<std::endl;
  dt.close();
  memspace.close();
  fspace.close();
  group->close();
  dataset->close();
  file->close();
  return(dataxd);
}


Eigen::ArrayXf read_dfvec_h5(const std::string h5file, const std::string groupname, const std::string dataname,const int offset,const int chunksize){

  H5FilePtr file=open_file(h5file);
  H5GroupPtr group= open_group(file,groupname);
  H5DataSetPtr dataset = open_dataset(group,dataname);
  DSetCreatPropList cparms= dataset->getCreatePlist();


  FloatType dt(PredType::NATIVE_FLOAT);
  hsize_t datadims[]={0};
  DataSpace fspace=dataset->getSpace();
  fspace.getSimpleExtentDims(datadims,NULL);
  //  std::cout<<"Full data is of dimensions"<<datadims[0]<<"x"<<datadims[1]<<std::endl;
  hsize_t vec_dims[1];
  if(offset+chunksize>datadims[0]){
    Rcpp::stop("offset+chunksize>datadims[0]!");
  }
  vec_dims[0]=chunksize;
  datadims[0]=chunksize;

  // std::cout<<"Full data is of size "<<datadim[0]<<std::endl;
  hsize_t offseta[1];
  offseta[0]=offset;
  Eigen::ArrayXf dataxd(chunksize);
  fspace.selectHyperslab(H5S_SELECT_SET,datadims,offseta);
  //  std::cout<<"Allocating matrix of size:"<<matrix_dims[0]<<"x"<<matrix_dims[1]<<std::endl;
  //  std::cout<<"Matrix starts at"<<row_offset<<"x"<<col_offset<<std::endl;
  DataSpace memspace(1,vec_dims);
  //  std::cout<<"Reading data"<<std::endl;
  dataset->read(dataxd.data(),dt,memspace,fspace);
  //  std::cout<<"Read complete!"<<std::endl;
  dt.close();
  memspace.close();
  fspace.close();
  group->close();
  dataset->close();
  file->close();
  return(dataxd);
}




void read_dvec_h5(const std::string h5file, const std::string groupname, const std::string dataname, double* data){

  H5FilePtr file=open_file(h5file);
  H5GroupPtr group= open_group(file,groupname);
  H5DataSetPtr dataset = open_dataset(group,dataname);
  DSetCreatPropList cparms= dataset->getCreatePlist();


  DataType dt= dataset->getDataType();
  hsize_t datadims[]={0};
  DataSpace fspace=dataset->getSpace();
  fspace.getSimpleExtentDims(datadims,NULL);
  //  std::cout<<"Full data is of dimensions"<<datadims[0]<<"x"<<datadims[1]<<std::endl;
  hsize_t vec_dims[1];
  int offset=0;
  int chunksize=datadims[0];

  if(offset+chunksize>datadims[0]){
    Rcpp::stop("offset+chunksize>datadims[0]!");
  }
  vec_dims[0]=chunksize;
  datadims[0]=chunksize;

  // std::cout<<"Full data is of size "<<datadim[0]<<std::endl;
  hsize_t offseta[1];
  offseta[0]=offset;
  fspace.selectHyperslab(H5S_SELECT_SET,datadims,offseta);
  //  std::cout<<"Allocating matrix of size:"<<matrix_dims[0]<<"x"<<matrix_dims[1]<<std::endl;
  //  std::cout<<"Matrix starts at"<<row_offset<<"x"<<col_offset<<std::endl;
  DataSpace memspace(1,vec_dims);
  //  std::cout<<"Reading data"<<std::endl;
  dataset->read(data,dt,memspace,fspace);
  //  std::cout<<"Read complete!"<<std::endl;
  dt.close();
  memspace.close();
  fspace.close();
  group->close();
  dataset->close();
  file->close();

}




Eigen::ArrayXd read_dvec_h5(const std::string h5file, const std::string groupname, const std::string dataname){

  H5FilePtr file=open_file(h5file);
  H5GroupPtr group= open_group(file,groupname);
  H5DataSetPtr dataset = open_dataset(group,dataname);
  DSetCreatPropList cparms= dataset->getCreatePlist();


  DataType dt= dataset->getDataType();
  hsize_t datadims[]={0};
  DataSpace fspace=dataset->getSpace();
  fspace.getSimpleExtentDims(datadims,NULL);
  //  std::cout<<"Full data is of dimensions"<<datadims[0]<<"x"<<datadims[1]<<std::endl;
  hsize_t vec_dims[1];
  int offset=0;
  int chunksize=datadims[0];
  Eigen::ArrayXd dataxd(chunksize);

  if(offset+chunksize>datadims[0]){
    Rcpp::stop("offset+chunksize>datadims[0]!");
  }
  vec_dims[0]=chunksize;
  datadims[0]=chunksize;

  // std::cout<<"Full data is of size "<<datadim[0]<<std::endl;
  hsize_t offseta[1];
  offseta[0]=offset;
  fspace.selectHyperslab(H5S_SELECT_SET,datadims,offseta);
  //  std::cout<<"Allocating matrix of size:"<<matrix_dims[0]<<"x"<<matrix_dims[1]<<std::endl;
  //  std::cout<<"Matrix starts at"<<row_offset<<"x"<<col_offset<<std::endl;
  DataSpace memspace(1,vec_dims);
  //  std::cout<<"Reading data"<<std::endl;
  dataset->read(dataxd.data(),dt,memspace,fspace);
  //  std::cout<<"Read complete!"<<std::endl;
  dt.close();
  memspace.close();
  fspace.close();
  group->close();
  dataset->close();
  file->close();
  return(dataxd);
}


std::vector<std::string> read_svec_h5(const std::string h5file, const std::string groupname, const std::string dataname,const int offset,const int chunksize){



  char** datasa;

  H5FilePtr file=open_file(h5file);
  H5GroupPtr group= open_group(file,groupname);
  H5DataSetPtr dataset = open_dataset(group,dataname);
  DSetCreatPropList cparms= dataset->getCreatePlist();


  DataType dt= dataset->getDataType();
  hsize_t datadims[]={0};
  DataSpace fspace=dataset->getSpace();
  fspace.getSimpleExtentDims(datadims,NULL);
  //  std::cout<<"Full data is of dimensions"<<datadims[0]<<"x"<<datadims[1]<<std::endl;
  hsize_t vec_dims[1];
  if(offset+chunksize>datadims[0]){
    Rcpp::stop("offset+chunksize>datadims[0]!");
  }

  datasa = new char*[chunksize];
  vec_dims[0]=chunksize;
  datadims[0]=chunksize;

  // std::cout<<"Full data is of size "<<datadim[0]<<std::endl;
  hsize_t offseta[1];
  offseta[0]=offset;

  fspace.selectHyperslab(H5S_SELECT_SET,datadims,offseta);
  //  std::cout<<"Allocating matrix of size:"<<matrix_dims[0]<<"x"<<matrix_dims[1]<<std::endl;
  //  std::cout<<"Matrix starts at"<<row_offset<<"x"<<col_offset<<std::endl;
  DataSpace memspace(1,vec_dims);
  //  std::cout<<"Reading data"<<std::endl;
  dataset->read((void*)datasa,dt,memspace,fspace);
  //  std::cout<<"Read complete!"<<std::endl;
  dt.close();
  memspace.close();
  fspace.close();
  group->close();
  dataset->close();
  file->close();
  std::vector<std::string> retvec(chunksize);
  for(int i=0;i<chunksize;i++){
    retvec[i]=std::string(datasa[i]);
  }


  return(retvec);
}


std::vector<std::string> read_svec_h5(const std::string h5file, const std::string groupname, const std::string dataname){

  size_t offset=0;

  char** datasa;

  H5FilePtr file=open_file(h5file);
  H5GroupPtr group= open_group(file,groupname);
  H5DataSetPtr dataset = open_dataset(group,dataname);
  DSetCreatPropList cparms= dataset->getCreatePlist();


  DataType dt= dataset->getDataType();
  hsize_t datadims[]={0};
  DataSpace fspace=dataset->getSpace();
  fspace.getSimpleExtentDims(datadims,NULL);
  //  std::cout<<"Full data is of dimensions"<<datadims[0]<<"x"<<datadims[1]<<std::endl;
  hsize_t vec_dims[1];
  size_t  chunksize =datadims[0];
  if(offset+chunksize>datadims[0]){
    Rcpp::stop("offset+chunksize>datadims[0]!");
  }

  datasa = new char*[chunksize];
  vec_dims[0]=chunksize;
  datadims[0]=chunksize;

  // std::cout<<"Full data is of size "<<datadim[0]<<std::endl;
  hsize_t offseta[1];
  offseta[0]=offset;

  fspace.selectHyperslab(H5S_SELECT_SET,datadims,offseta);
  //  std::cout<<"Allocating matrix of size:"<<matrix_dims[0]<<"x"<<matrix_dims[1]<<std::endl;
  //  std::cout<<"Matrix starts at"<<row_offset<<"x"<<col_offset<<std::endl;
  DataSpace memspace(1,vec_dims);
  //  std::cout<<"Reading data"<<std::endl;
  dataset->read((void*)datasa,dt,memspace,fspace);
  //  std::cout<<"Read complete!"<<std::endl;
  dt.close();
  memspace.close();
  fspace.close();
  group->close();
  dataset->close();
  file->close();
  std::vector<std::string> retvec(chunksize);
  for(int i=0;i<chunksize;i++){
    retvec[i]=std::string(datasa[i]);
  }


  return(retvec);
}




Eigen::ArrayXi read_ivec_h5(const std::string h5file, const std::string groupname, const std::string dataname,const int offset,const int chunksize){

  H5FilePtr file=open_file(h5file);
  H5GroupPtr group= open_group(file,groupname);
  H5DataSetPtr dataset = open_dataset(group,dataname);
  DSetCreatPropList cparms= dataset->getCreatePlist();


  DataType dt= PredType::NATIVE_INT;
  std::vector<hsize_t> datadims;
  DataSpace fspace=dataset->getSpace();
  hsize_t nd = fspace.getSimpleExtentNdims();
  datadims.resize(nd);
  fspace.getSimpleExtentDims(datadims.data(),NULL);

  //  std::cout<<"Full data is of dimensions"<<datadims[0]<<"x"<<datadims[1]<<std::endl;
  hsize_t vec_dims[1];
  if(offset+chunksize>datadims[0]){
    Rcpp::stop("offset+chunksize>datadims[0]!");
  }
  vec_dims[0]=chunksize;
  datadims[0]=chunksize;

  // std::cout<<"Full data is of size "<<datadim[0]<<std::endl;
  hsize_t offseta[1];
  offseta[0]=offset;
  Eigen::ArrayXi dataxi(chunksize);
  fspace.selectHyperslab(H5S_SELECT_SET,datadims.data(),offseta);
  //  std::cout<<"Allocating matrix of size:"<<matrix_dims[0]<<"x"<<matrix_dims[1]<<std::endl;
  //  std::cout<<"Matrix starts at"<<row_offset<<"x"<<col_offset<<std::endl;
  DataSpace memspace(1,vec_dims);
  //  std::cout<<"Reading data"<<std::endl;
  dataset->read(dataxi.data(),dt,memspace,fspace);
  //  std::cout<<"Read complete!"<<std::endl;
  dt.close();
  memspace.close();
  fspace.close();
  group->close();
  dataset->close();
  file->close();
  return(dataxi);
}









void read_ivec_h5(const std::string h5file, const std::string groupname, const std::string dataname,const int offset,const int chunksize, int* data){

  H5FilePtr file=open_file(h5file);
  H5GroupPtr group= open_group(file,groupname);
  H5DataSetPtr dataset = open_dataset(group,dataname);
  DSetCreatPropList cparms= dataset->getCreatePlist();


  DataType dt= PredType::NATIVE_INT;
  hsize_t datadims[]={0};
  DataSpace fspace=dataset->getSpace();
  fspace.getSimpleExtentDims(datadims,NULL);
  //  std::cout<<"Full data is of dimensions"<<datadims[0]<<"x"<<datadims[1]<<std::endl;
  hsize_t vec_dims[1];
  if(offset+chunksize>datadims[0]){
    Rcpp::stop("offset+chunksize>datadims[0]!");
  }
  vec_dims[0]=chunksize;
  datadims[0]=chunksize;

  // std::cout<<"Full data is of size "<<datadim[0]<<std::endl;
  hsize_t offseta[1];
  offseta[0]=offset;
  fspace.selectHyperslab(H5S_SELECT_SET,datadims,offseta);
  //  std::cout<<"Allocating matrix of size:"<<matrix_dims[0]<<"x"<<matrix_dims[1]<<std::endl;
  //  std::cout<<"Matrix starts at"<<row_offset<<"x"<<col_offset<<std::endl;
  DataSpace memspace(1,vec_dims);
  //  std::cout<<"Reading data"<<std::endl;
  dataset->read(data,dt,memspace,fspace);
  //  std::cout<<"Read complete!"<<std::endl;
  dt.close();
  memspace.close();
  fspace.close();
  group->close();
  dataset->close();
  file->close();

}


Eigen::ArrayXi read_ivec_h5(const std::string h5file, const std::string groupname, const std::string dataname){

  H5FilePtr file=open_file(h5file);
  H5GroupPtr group= open_group(file,groupname);
  H5DataSetPtr dataset = open_dataset(group,dataname);
  DSetCreatPropList cparms= dataset->getCreatePlist();


  DataType dt= PredType::NATIVE_INT;
  //  hsize_t datadims[]={0};

  DataSpace fspace=dataset->getSpace();
  size_t nd = fspace.getSimpleExtentNdims();
  std::vector<hsize_t> datadims(nd);
  fspace.getSimpleExtentDims(datadims.data(),NULL);
  //  std::cout<<"Full data is of dimensions"<<datadims[0]<<"x"<<datadims[1]<<std::endl;
  //  std::vector<hsize_t> chunksizes=datadims;
  std::vector<hsize_t> offset(nd);
  std::fill(offset.begin(),offset.end(),0);
  //  int chunksize=datadims[0];
  hsize_t tot_size=1;
  tot_size=std::accumulate(datadims.begin(),datadims.end(),tot_size,std::multiplies<hsize_t>());

  // std::cout<<"tot_size is of returned array is: "<<tot_size<<std::endl;
  if(tot_size==0){

    Rcpp::stop("total size of returned array shouldn't be zero");
  }
  Eigen::ArrayXi arrayxi(tot_size);
  arrayxi.setZero();

  // std::cout<<"Full data is of size: ";
  // for(auto veci:datadims){
  //   std::cout<<veci<<" ";
  // }
  // std::cout<<std::endl;
  // hsize_t offseta[1];
  // offseta[0]=offset;
  fspace.selectHyperslab(H5S_SELECT_SET,datadims.data(),offset.data());
  //  std::cout<<"Allocating matrix of size:"<<matrix_dims[0]<<"x"<<matrix_dims[1]<<std::endl;
  //  std::cout<<"Matrix starts at"<<row_offset<<"x"<<col_offset<<std::endl;
  DataSpace memspace(nd,datadims.data());
  //  std::cout<<"Reading data"<<std::endl;
  try{
    dataset->read(arrayxi.data(),dt,memspace,fspace);
  }catch(DataSetIException error){
    error.printError();
    Rcpp::stop("Error reading column!");
  }
  dt.close();
  memspace.close();
  fspace.close();
  group->close();
  dataset->close();
  file->close();
  return(arrayxi);
}


Eigen::ArrayXi read_uivec_h5(const std::string h5file, const std::string groupname, const std::string dataname){

  H5FilePtr file=open_file(h5file);
  H5GroupPtr group= open_group(file,groupname);
  H5DataSetPtr dataset = open_dataset(group,dataname);
  DSetCreatPropList cparms= dataset->getCreatePlist();

  DataType dt= PredType::NATIVE_INT;
  // DataType dt= dataset->getDataType();
  hsize_t datadims[]={0};
  DataSpace fspace=dataset->getSpace();
  fspace.getSimpleExtentDims(datadims,NULL);
  //  std::cout<<"Full data is of dimensions"<<datadims[0]<<"x"<<datadims[1]<<std::endl;
  hsize_t vec_dims[1];
  int offset=0;
  int chunksize=datadims[0];
  std::vector<unsigned int> tai(chunksize);
  Eigen::ArrayXi arrayxi(chunksize);
  arrayxi.setZero();


  if(offset+chunksize>datadims[0]){
    Rcpp::stop("offset+chunksize>datadims[0]!");
  }
  vec_dims[0]=chunksize;
  datadims[0]=chunksize;

  // std::cout<<"Full data is of size "<<datadim[0]<<std::endl;
  hsize_t offseta[1];
  offseta[0]=offset;
  fspace.selectHyperslab(H5S_SELECT_SET,datadims,offseta);
  //  std::cout<<"Allocating matrix of size:"<<matrix_dims[0]<<"x"<<matrix_dims[1]<<std::endl;
  //  std::cout<<"Matrix starts at"<<row_offset<<"x"<<col_offset<<std::endl;
  DataSpace memspace(1,vec_dims);
  //  std::cout<<"Reading data"<<std::endl;
  try{
    dataset->read(&tai[0],dt,memspace,fspace);
  }catch(DataSetIException error){
    error.printError();
    Rcpp::stop("Error reading column!");
  }



  dt.close();
  memspace.close();
  fspace.close();
  group->close();
  dataset->close();
  file->close();
  for(size_t i=0;i<tai.size();i++){
    arrayxi[i]=(int)tai[i];
  }
  return(arrayxi);
}





void read_2d_cindex_h5(const std::string h5file,const std::string groupname, const std::string dataname, const c_arrayxi_internal indvec, Matrix_internal retmat){

  // Rcpp::Rcout<<"Starting!"<<std::endl;

  const  int minind (indvec.minCoeff());
  const  int maxind(indvec.maxCoeff());
  const  int t_chunksize((maxind-minind)+1);
  // Rcpp::Rcout<<"t_chunksize is :"<<t_chunksize<<std::endl;
  // Rcpp::Rcout<<"maxind-minind+1 is "<<maxind-minind+1<<std::endl;


  // Rcpp::Rcout<<"indvec is :"<<std::endl<<indvec<<std::endl;
  // Rcpp::Rcout<<"minind :"<<minind<<std::endl;
  // Rcpp::Rcout<<"maxind :"<<maxind<<std::endl;

  const int rownum(get_rownum_h5(h5file,groupname,dataname));
  retmat.resize(rownum,indvec.size());
  Eigen::ArrayXi rowind(rownum);

  rowind.setLinSpaced(rownum,0,rownum-1);

  const Eigen::ArrayXi newind=indvec-minind;


  //Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> readmat(rownum,t_chunksize);
  double *x = new double[rownum*t_chunksize];
  Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> > readmat((double*)x,rownum,t_chunksize);
  // Rcpp::Rcout<<"Sum of x is "<<readmat.sum()<<std::endl;
  memset((double* )x,0,rownum*t_chunksize*sizeof(double));
  // Rcpp::Rcout<<"Sum of x is "<<readmat.sum()<<std::endl;
  //  Rcpp::Rcout<<"chunksize is :"<<t_chunksize<<std::endl;
  read_2ddmat_h5(h5file,groupname,dataname,0,minind-1,rownum,t_chunksize,x);

  //    Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> > readmat(x,rownum,t_chunksize);

  // Rcpp::Rcout<<"Sum of x is "<<readmat.sum()<<std::endl;

  //Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> retmat(rownum,t_chunksize);
  int indsize=indvec.size();
  for(int i=0; i<indsize; i++){
    retmat.col(i)=readmat.col(newind(i));
  }
  delete[] x;
}

void read_2df_cindex_h5(const std::string h5file,const std::string groupname, const std::string dataname, const c_arrayxi_internal indvec, Matrix_finternal retmat){

  // Rcpp::Rcout<<"Starting!"<<std::endl;

  const  int minind (indvec.minCoeff());
  const  int maxind(indvec.maxCoeff());
  const  int t_chunksize((maxind-minind)+1);
  // Rcpp::Rcout<<"t_chunksize is :"<<t_chunksize<<std::endl;
  // Rcpp::Rcout<<"maxind-minind+1 is "<<maxind-minind+1<<std::endl;


  // Rcpp::Rcout<<"indvec is :"<<std::endl<<indvec<<std::endl;
  // Rcpp::Rcout<<"minind :"<<minind<<std::endl;
  // Rcpp::Rcout<<"maxind :"<<maxind<<std::endl;

  const int rownum(get_rownum_h5(h5file,groupname,dataname));
  retmat.resize(rownum,indvec.size());
  Eigen::ArrayXi rowind(rownum);

  rowind.setLinSpaced(rownum,0,rownum-1);

  const Eigen::ArrayXi newind=indvec-minind;


  //Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> readmat(rownum,t_chunksize);

  float *x = new float[rownum*t_chunksize];
  Eigen::Map<Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> > readmat((float*)x,rownum,t_chunksize);
  // Rcpp::Rcout<<"Sum of x is "<<readmat.sum()<<std::endl;
  memset((float* )x,0,rownum*t_chunksize*sizeof(float));
  // Rcpp::Rcout<<"Sum of x is "<<readmat.sum()<<std::endl;
  //  Rcpp::Rcout<<"chunksize is :"<<t_chunksize<<std::endl;
  read_2dfmat_h5(h5file,groupname,dataname,0,minind-1,rownum,t_chunksize,x);

  //    Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> > readmat(x,rownum,t_chunksize);

  // Rcpp::Rcout<<"Sum of x is "<<readmat.sum()<<std::endl;

  //Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> retmat(rownum,t_chunksize);
  int indsize=indvec.size();
  for(int i=0; i<indsize; i++){
    retmat.col(i)=readmat.col(newind(i));
  }
  delete[] x;
}








void read_2d_cindex_chunk_h5(const std::string h5file,const std::string groupname, const std::string dataname, const c_arrayxi_internal indvec, Matrix_internal retmat,const size_t chunksize){

  // Rcpp::Rcout<<"Starting!"<<std::endl;

  const  int minind (indvec.minCoeff());
  const  int maxind(indvec.maxCoeff());
  const  int t_chunksize((maxind-minind)+1);
  // Rcpp::Rcout<<"t_chunksize is :"<<t_chunksize<<std::endl;
  // Rcpp::Rcout<<"maxind-minind+1 is "<<maxind-minind+1<<std::endl;

  size_t p=indvec.size();

  if(t_chunksize>chunksize){
    size_t midp= p/2;
    read_2d_cindex_chunk_h5(h5file,groupname,dataname,indvec.head(midp),retmat.leftCols(midp));
    read_2d_cindex_chunk_h5(h5file,groupname,dataname,indvec.tail(p-midp),retmat.rightCols(p-midp));
  }else{
    read_2d_cindex_h5(h5file,groupname,dataname,indvec,retmat);
  }

}



Eigen::ArrayXd read_1d_cindex_h5(const std::string h5file,const std::string groupname, const std::string dataname, const  c_arrayxi_internal indvec){

  int minind = indvec.minCoeff();
  int maxind = indvec.maxCoeff();
  int chunksize=maxind-minind+1;
  //  Rcpp::Rcout<<"rowind is :"<<std::endl<<rowind<<std::endl;
  ArrayXi newind=indvec-minind;
  ArrayXi colind(1);
  colind[0]=0;
  //  Rcpp::Rcout<<"colind is :"<<std::endl<<newind<<std::endl;
  Eigen::ArrayXd readvec(chunksize,1);
  //  Rcpp::Rcout<<"About to read  matrix"<<std::endl;
  readvec=read_dvec_h5(h5file,groupname,dataname,minind-1,chunksize);
  //  Rcpp::Rcout<<"Finished reading matrix"<<std::endl;
  //  Rcpp::Rcout<<"readmat has "<<readmat.rows()<<" rows and "<<readmat.cols()<<" cols"<<std::endl;
  //  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> readmat(rownum,chunksize);
  //  Eigen::MatrixXd tmat=indexing(readmat,newind,newind);
  int indsize=indvec.size();
  Eigen::ArrayXd retvec(indsize);

  for(int i=0; i<indsize; i++){
    retvec.coeffRef(i)=readvec.coeffRef(newind.coeff(i));
  }
  return(retvec);
}

Eigen::ArrayXf read_1df_cindex_h5(const std::string h5file,const std::string groupname, const std::string dataname, const  c_arrayxi_internal indvec){

  int minind = indvec.minCoeff();
  int maxind = indvec.maxCoeff();
  int chunksize=maxind-minind+1;
  //  Rcpp::Rcout<<"rowind is :"<<std::endl<<rowind<<std::endl;
  ArrayXi newind=indvec-minind;
  ArrayXi colind(1);
  colind[0]=0;
  //  Rcpp::Rcout<<"colind is :"<<std::endl<<newind<<std::endl;
  Eigen::ArrayXf readvec(chunksize,1);
  //  Rcpp::Rcout<<"About to read  matrix"<<std::endl;
  readvec=read_dfvec_h5(h5file,groupname,dataname,minind-1,chunksize);
  //  Rcpp::Rcout<<"Finished reading matrix"<<std::endl;
  //  Rcpp::Rcout<<"readmat has "<<readmat.rows()<<" rows and "<<readmat.cols()<<" cols"<<std::endl;
  //  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> readmat(rownum,chunksize);
  //  Eigen::MatrixXd tmat=indexing(readmat,newind,newind);
  int indsize=indvec.size();
  Eigen::ArrayXf retvec(indsize);

  for(int i=0; i<indsize; i++){
    retvec.coeffRef(i)=readvec.coeffRef(newind.coeff(i));
  }
  return(retvec);
}


std::vector<std::string> read_s1d_index_h5(const std::string h5file,const std::string groupname, const std::string dataname, const  c_arrayxi_internal indvec){

  int minind = indvec.minCoeff();
  int maxind = indvec.maxCoeff();
  int chunksize=maxind-minind+1;
  //  Rcpp::Rcout<<"rowind is :"<<std::endl<<rowind<<std::endl;
  ArrayXi newind=indvec-minind;
  ArrayXi colind(1);
  colind[0]=0;
  //  Rcpp::Rcout<<"colind is :"<<std::endl<<newind<<std::endl;
  std::vector<std::string> readvec(chunksize);
  //  Rcpp::Rcout<<"About to read  matrix"<<std::endl;
  readvec=read_svec_h5(h5file,groupname,dataname,minind-1,chunksize);
  int indsize=indvec.size();
  std::vector<std::string> retvec(indsize);
  for(int i=0;i<indsize;i++){
    retvec[i]=readvec[newind[i]];
  }
  return(retvec);
}






