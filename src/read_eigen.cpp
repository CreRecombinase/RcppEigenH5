// #include <RcppEigen.h>
// #include "RcppEigenH5.h"
//
//
// template <typename T>
// struct DatatypeSpecialization;
//
// // floating-point types
// template <>
// struct DatatypeSpecialization<float>
// {
//   static inline const H5::DataType * get (void)
//   {
//     return &H5::PredType::NATIVE_FLOAT;
//   }
// };
//
// template <>
// struct DatatypeSpecialization<double>
// {
//   static inline const H5::DataType * get (void)
//   {
//     return &H5::PredType::NATIVE_DOUBLE;
//   }
// };
//
// template <>
// struct DatatypeSpecialization<long double>
// {
//   static inline const H5::DataType * get (void)
//   {
//     return &H5::PredType::NATIVE_LDOUBLE;
//   }
// };
//
// // integer types
//
// template <>
// struct DatatypeSpecialization<short>
// {
//   static inline const H5::DataType * get (void)
//   {
//     return &H5::PredType::NATIVE_SHORT;
//   }
// };
//
// template <>
// struct DatatypeSpecialization<unsigned short>
// {
//   static inline const H5::DataType * get (void)
//   {
//     return &H5::PredType::NATIVE_USHORT;
//   }
// };
//
// template <>
// struct DatatypeSpecialization<int>
// {
//   static inline const H5::DataType * get (void)
//   {
//     return &H5::PredType::NATIVE_INT;
//   }
// };
//
// template <>
// struct DatatypeSpecialization<unsigned int>
// {
//   static inline const H5::DataType * get (void)
//   {
//     return &H5::PredType::NATIVE_UINT;
//   }
// };
//
// template <>
// struct DatatypeSpecialization<long>
// {
//   static inline const H5::DataType * get (void)
//   {
//     return &H5::PredType::NATIVE_LONG;
//   }
// };
//
// template <>
// struct DatatypeSpecialization<unsigned long>
// {
//   static inline const H5::DataType * get (void)
//   {
//     return &H5::PredType::NATIVE_ULONG;
//   }
// };
//
// template <>
// struct DatatypeSpecialization<long long>
// {
//   static inline const H5::DataType * get (void)
//   {
//     return &H5::PredType::NATIVE_LLONG;
//   }
// };
//
// template <>
// struct DatatypeSpecialization<unsigned long long>
// {
//   static inline const H5::DataType * get (void)
//   {
//     return &H5::PredType::NATIVE_ULLONG;
//   }
// };
//
//
// // string types, to be used mainly for attributes
//
// template <>
// struct DatatypeSpecialization<const char *>
// {
//   static inline const H5::DataType * get (void)
//   {
//     static const H5::StrType strtype(0, H5T_VARIABLE);
//     return &strtype;
//   }
// };
//
// template <>
// struct DatatypeSpecialization<char *>
// {
//   static inline const H5::DataType * get (void)
//   {
//     static const H5::StrType strtype(0, H5T_VARIABLE);
//     return &strtype;
//   }
// };
//
// // XXX: for some unknown reason the following two functions segfault if
// // H5T_VARIABLE is used.  The passed strings should still be null-terminated,
// // so this is a bit worrisome.
//
// template <std::size_t N>
// struct DatatypeSpecialization<const char [N]>
// {
//   static inline const H5::DataType * get (void)
//   {
//     static const H5::StrType strtype(0, N);
//     return &strtype;
//   }
// };
//
// template <std::size_t N>
// struct DatatypeSpecialization<char [N]>
// {
//   static inline const H5::DataType * get (void)
//   {
//     static const H5::StrType strtype(0, N);
//     return &strtype;
//   }
// };

// //
// template<typename Derived> Eigen::EigenBase<Derived> read_tvec_h5(const std::string h5file,
//                                      const std::string groupname,
//                                      const std::string dataname){
//
//   H5FilePtr file=open_file(h5file);
//   H5GroupPtr group= open_group(file,groupname);
//   H5DataSetPtr dataset = open_dataset(group,dataname);
//   DSetCreatPropList cparms= dataset->getCreatePlist();
//
// DataType
//   DataType dt= PredType::NATIVE_INT;
//   //  hsize_t datadims[]={0};
//
//   DataSpace fspace=dataset->getSpace();
//   size_t nd = fspace.getSimpleExtentNdims();
//   std::vector<hsize_t> datadims(nd);
//   fspace.getSimpleExtentDims(datadims.data(),NULL);
//   //  std::cout<<"Full data is of dimensions"<<datadims[0]<<"x"<<datadims[1]<<std::endl;
//   //  std::vector<hsize_t> chunksizes=datadims;
//   std::vector<hsize_t> offset(nd);
//   std::fill(offset.begin(),offset.end(),0);
//   //  int chunksize=datadims[0];
//   hsize_t tot_size=1;
//   tot_size=std::accumulate(datadims.begin(),datadims.end(),tot_size,std::multiplies<hsize_t>());
//
//   // std::cout<<"tot_size is of returned array is: "<<tot_size<<std::endl;
//   if(tot_size==0){
//
//     Rcpp::stop("total size of returned array shouldn't be zero");
//   }
//   Eigen::ArrayXi arrayxi(tot_size);
//   arrayxi.setZero();
//
//   // std::cout<<"Full data is of size: ";
//   for(auto veci:datadims){
//     std::cout<<veci<<" ";
//   }
//   std::cout<<std::endl;
//   // hsize_t offseta[1];
//   // offseta[0]=offset;
//   fspace.selectHyperslab(H5S_SELECT_SET,datadims.data(),offset.data());
//   //  std::cout<<"Allocating matrix of size:"<<matrix_dims[0]<<"x"<<matrix_dims[1]<<std::endl;
//   //  std::cout<<"Matrix starts at"<<row_offset<<"x"<<col_offset<<std::endl;
//   DataSpace memspace(nd,datadims.data());
//   //  std::cout<<"Reading data"<<std::endl;
//   try{
//     dataset->read(arrayxi.data(),dt,memspace,fspace);
//   }catch(DataSetIException error){
//     error.printError();
//     Rcpp::stop("Error reading column!");
//   }
//   dt.close();
//   memspace.close();
//   fspace.close();
//   group->close();
//   dataset->close();
//   file->close();
//   return(arrayxi);
//
// }
//







// template <typename Type> using ArrayX = Eigen::Array< Type,Eigen::Dynamic,1,Eigen::Dynamic>;
// template <typename Type> using MatrixX = Eigen::Array< Type,Eigen::Dynamic,Eigen::Dynamic,Eigen::Dynamic>;

// template<typename T>
// struct data_rank{
//   static const int rank=2;
// };
//
// template<ArrayX<T>>
// struct data_rank{
//   static const int rank=2;
// };

//
//
// template <typename D> class vec_chunk_ind{
//
//   const size_t offset;
//   const size_t chunksize;
// public:
//   vec_chunk_ind(const size_t offset_,const size_t chunksize_):offset(offset_),chunksize(chunksize_){};
//
//   void set_fspace(H5DataSpacePtr p_fspace){
//     hsize_t datadims[]={0};
//     const int datarank=p_fspace->getSimpleExtentNdims();
//     if(datarank>1){
//       Rcpp::stop("Can't use vector reading method on data with more than one dimeension!");
//     }
//     //    DataSpace fspace=dataset->getSpace();
//     p_fspace->getSimpleExtentDims(datadims,NULL);
//     //  std::cout<<"Full data is of dimensions"<<datadims[0]<<"x"<<datadims[1]<<std::endl;
//    // hsize_t vec_dims[1];
//     if(offset+chunksize>datadims[0]){
//       Rcpp::stop("offset+chunksize>datadims[0]!");
//     }
// //    vec_dims[0]=chunksize;
//     datadims[0]=chunksize;
//
//     // std::cout<<"Full data is of size "<<datadim[0]<<std::endl;
//     hsize_t offseta[1];
//     offseta[0]=offset;
//     p_fspace->selectHyperslab(H5S_SELECT_SET,datadims,offseta);
//   }
//
//
// };
/*

template<typename Ind_p, typename Data_p> void read_vec_h5(const std::string h5file, const std::string groupname, const std::string dataname,const Ind_p index_method, Data_p data_method){
  H5FilePtr file=open_file(h5file);
  H5GroupPtr group= open_group(file,groupname);
  H5DataSetPtr dataset = open_dataset(group,dataname);
  DSetCreatPropList cparms= dataset->getCreatePlist();


  DataType dt= dataset->getDataType();
  hsize_t datadims[]={0};
  H5DataSpacePtr fspace=std::make_shared<DataSpace>(dataset->getSpace());

  // fspace.getSimpleExtentDims(datadims,NULL);
  // //  std::cout<<"Full data is of dimensions"<<datadims[0]<<"x"<<datadims[1]<<std::endl;
  // hsize_t vec_dims[1];
  // if(offset+chunksize>datadims[0]){
  //   Rcpp::stop("offset+chunksize>datadims[0]!");
  // }
  // vec_dims[0]=chunksize;
  // datadims[0]=chunksize;

  // std::cout<<"Full data is of size "<<datadim[0]<<std::endl;
  // hsize_t offseta[1];
  // offseta[0]=offset;
  // Eigen::ArrayXd dataxd(chunksize);
  // fspace.selectHyperslab(H5S_SELECT_SET,datadims,offseta);
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
  Eigen::ArrayXi dataxi(chunksize);
  fspace.selectHyperslab(H5S_SELECT_SET,datadims,offseta);
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


  DataType dt= dataset->getDataType();
  hsize_t datadims[]={0};
  DataSpace fspace=dataset->getSpace();
  fspace.getSimpleExtentDims(datadims,NULL);
  //  std::cout<<"Full data is of dimensions"<<datadims[0]<<"x"<<datadims[1]<<std::endl;
  hsize_t vec_dims[1];
  int offset=0;
  int chunksize=datadims[0];
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


  DataType dt= dataset->getDataType();
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

 */
