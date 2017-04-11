#include "RcppEigenH5.h"
#include <type_traits>
//[[Rcpp::depends(RcppEigen)]]
using namespace Eigen;



void write_mat_chunk_h5(const std::string h5file, const std::string groupname, const std::string dataname, Matrix_internal data, c_arrayxi_internal offsets){



  FloatType ftypew(PredType::NATIVE_DOUBLE);
  H5FilePtr file =create_or_open_file(h5file);

  size_t row_offset = offsets[0];
  size_t col_offset = offsets[1];



  H5GroupPtr group =create_or_open_group(file,groupname);


  H5DataSetPtr dataset =open_dataset(group,dataname);




  Eigen::MatrixXd tdata = data;

  bool doTranspose = check_transpose(dataset);

  if(doTranspose){
    //    Rcout<<"Transposing!"<<std::endl;
    // Rcout<<"data  is "<<print_dims(tdata)<<std::endl;
    std::swap(row_offset,col_offset);
    tdata.transposeInPlace();
    // Rcout<<"data  is "<<print_dims(tdata)<<std::endl;
    // Rcout<<data<<std::endl;
  }



  size_t rowsize=tdata.rows();
  size_t colsize=tdata.cols();

  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> writemat = tdata;







  hsize_t datadims[]={0,0};
  DataSpace fspace=dataset->getSpace();
  fspace.getSimpleExtentDims(datadims,NULL);
  //  std::cout<<"Full data is of dimensions"<<datadims[0]<<"x"<<datadims[1]<<std::endl;
  hsize_t matrix_dims[2];
  if(row_offset+rowsize>datadims[0]){
    Rcpp::stop("row_offset ("+std::to_string(static_cast<long long>(row_offset))+") + row_chunksize ("+std::to_string(static_cast<long long>(rowsize))+") >datadims[0] ("+std::to_string(static_cast<long long>(datadims[0]))+")");
  }
  if(col_offset+colsize>datadims[1]){
    Rcpp::stop("col_offset ("+std::to_string(static_cast<long long>(col_offset))+") + col_chunksize ("+std::to_string(static_cast<long long>(colsize))+") >datadims[1] ("+std::to_string(static_cast<long long>(datadims[1]))+")");
  }





  DataSpace* fdataspace;
  try{
    fdataspace= new DataSpace(dataset->getSpace());
  }catch(DataSpaceIException error){
    error.printError();
    Rcpp::stop("Error creating memory dataspace ");
  }

  hsize_t datadim[2];
  fdataspace->getSimpleExtentDims(datadim,NULL);

  hsize_t memdim[]={rowsize,colsize};

  DataSpace *mspace;
  try{
    mspace= new DataSpace(2,memdim); //Size of first dataset (in memory, can be bigger or smaller than size on disk, depending on how much you're writing)
  }catch(DataSpaceIException error){
    error.printError();
    Rcpp::stop("Error creating memory dataspace ");
  }
  hsize_t odim[]={row_offset,col_offset};//dimension of each offset (current_chunk*chunksize)
  hsize_t stridea[]={1,1};
  hsize_t blocka[]={1,1};


  fdataspace->selectHyperslab( H5S_SELECT_SET, memdim,odim,stridea,blocka);


  //  std::cout<<"Starting to write data"<<std::endl;
  try{
    dataset->write(writemat.data(),PredType::NATIVE_DOUBLE,*mspace,*fdataspace);
  }  catch( DataSetIException error )
  {
    error.printError();
    Rcpp::stop("Error writing file");
  }
  //  std::cout<<"Data sucessfully written"<<std::endl;
  try{
    file->flush(H5F_SCOPE_GLOBAL);
  }catch(FileIException error)
  {
    error.printError();
    Rcpp::stop("Error flushing file");
  }
  //  std::cout<<"File flushed"<<std::endl;
  dataset->close();
  fdataspace->close();
  mspace->close();
  group->close();
  file->close();
  fdataspace->close();
  mspace->close();
}



void write_mat_h5(const std::string h5file, const std::string groupname, const std::string dataname, Matrix_internal data, const int deflate_level){

  size_t rowsize=data.rows();
  size_t colsize=data.cols();

  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> writemat(data);

  FloatType ftypew(PredType::NATIVE_DOUBLE);
  H5FilePtr file =create_or_open_file(h5file);


  size_t rchunksize=rowsize;
  size_t cchunksize=1000;
  if(cchunksize>colsize){
    cchunksize=colsize;
  }

  H5GroupPtr group =create_or_open_group(file,groupname);
  std::vector<hsize_t> cumdim{rowsize,colsize};
  std::vector<hsize_t> maxdim{rowsize,colsize};
  std::vector<hsize_t> chunkdim{rchunksize,cchunksize};

  H5DataSetPtr dataset =create_or_open_dataset(group,dataname,ftypew,cumdim,maxdim,chunkdim,deflate_level);
  write_transpose(dataset,false);
  DataSpace* fdataspace;
  try{
    fdataspace= new DataSpace(dataset->getSpace());
  }catch(DataSpaceIException error){
    error.printError();
    Rcpp::stop("Error creating memory dataspace ");
  }

  hsize_t datadim[2];
  fdataspace->getSimpleExtentDims(datadim,NULL);

  hsize_t memdim[]={rowsize,colsize};

  DataSpace *mspace;
  try{
    mspace= new DataSpace(2,memdim); //Size of first dataset (in memory, can be bigger or smaller than size on disk, depending on how much you're writing)
  }catch(DataSpaceIException error){
    error.printError();
    Rcpp::stop("Error creating memory dataspace ");
  }
  hsize_t odim[]={0,0};//dimension of each offset (current_chunk*chunksize)
  hsize_t stridea[]={1,1};
  hsize_t blocka[]={1,1};


  fdataspace->selectHyperslab( H5S_SELECT_SET, memdim,odim,stridea,blocka);


  //  std::cout<<"Starting to write data"<<std::endl;
  try{
    dataset->write(writemat.data(),PredType::NATIVE_DOUBLE,*mspace,*fdataspace);
  }  catch( DataSetIException error )
  {
    error.printError();
    Rcpp::stop("Error writing file");
  }
  //  std::cout<<"Data sucessfully written"<<std::endl;
  try{
    file->flush(H5F_SCOPE_GLOBAL);
  }catch(FileIException error)
  {
    error.printError();
    Rcpp::stop("Error flushing file");
  }
  //  std::cout<<"File flushed"<<std::endl;
  dataset->close();
  fdataspace->close();
  mspace->close();
  group->close();
  file->close();
  fdataspace->close();
  mspace->close();

}







void write_mat_h5(const std::string h5file, const std::string groupname, const std::string dataname, Matrix_internal data, Eigen::ArrayXi chunksize, const int deflate_level,bool doTranspose){

  if(doTranspose){

    // data=data.transposeInPlace();
    std::swap(chunksize.coeffRef(0),chunksize.coeffRef(1));
    // int tsize=chunksize.coeff(0);
    // chunksize.coeffRef(0)=chunksize.coeff(1);
    // chunksize.coeffRef(1)=tsize;
  }
  Eigen::MatrixXd tdata = data;

  // if(tdata.coeff(0,0)!=data.coeff(0,0)){
  //   Rcpp::stop("copy didn't work!");
  // }


  if(doTranspose){
    //    Rcout<<"Transposing!"<<std::endl;
    // Rcout<<"data  is "<<print_dims(tdata)<<std::endl;
    tdata.transposeInPlace();
    // Rcout<<"data  is "<<print_dims(tdata)<<std::endl;
    // Rcout<<data<<std::endl;
  }else{

  }
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> writemat(tdata);
  if(doTranspose){
    if(writemat.rows()!=tdata.rows()){
      Rcpp::stop("Transposition didn't happen!(doTranspose==TRUE)");
    }
  }else{
    if(writemat.rows()!=data.rows()){
      Rcout<<"writemat is "<<print_dims(writemat)<<std::endl;
//      Rcout<<writemat<<std::endl;
      Rcout<<"data  is "<<print_dims(data)<<std::endl;
//      Rcout<<data<<std::endl;
      Rcpp::stop("Transposition happend!(doTranspose==FALSE)");
    }
  }
//  Rcout<<"Writemat is "<<print_dims(writemat)<<std::endl;



  size_t rowsize=writemat.rows();
  size_t colsize=writemat.cols();

  FloatType ftypew(PredType::NATIVE_DOUBLE);
  H5FilePtr file =create_or_open_file(h5file);


  size_t rchunksize=chunksize[0];
  size_t cchunksize=chunksize[1];
  if(cchunksize>colsize){
    cchunksize=colsize;
  }

  H5GroupPtr group =create_or_open_group(file,groupname);
  std::vector<hsize_t> cumdim{rowsize,colsize};
  std::vector<hsize_t> maxdim{rowsize,colsize};
  std::vector<hsize_t> chunkdim{rchunksize,cchunksize};

  H5DataSetPtr dataset =create_or_open_dataset(group,dataname,ftypew,cumdim,maxdim,chunkdim,deflate_level);
  write_transpose(dataset,doTranspose);


  DataSpace* fdataspace;
  try{
    fdataspace= new DataSpace(dataset->getSpace());
  }catch(DataSpaceIException error){
    error.printError();
    Rcpp::stop("Error creating memory dataspace ");
  }

  hsize_t datadim[2];
  fdataspace->getSimpleExtentDims(datadim,NULL);

  hsize_t memdim[]={rowsize,colsize};

  DataSpace *mspace;
  try{
    mspace= new DataSpace(2,memdim); //Size of first dataset (in memory, can be bigger or smaller than size on disk, depending on how much you're writing)
  }catch(DataSpaceIException error){
    error.printError();
    Rcpp::stop("Error creating memory dataspace ");
  }
  hsize_t odim[]={0,0};//dimension of each offset (current_chunk*chunksize)
  hsize_t stridea[]={1,1};
  hsize_t blocka[]={1,1};


  fdataspace->selectHyperslab( H5S_SELECT_SET, memdim,odim,stridea,blocka);


  //  std::cout<<"Starting to write data"<<std::endl;
  try{
    dataset->write(writemat.data(),PredType::NATIVE_DOUBLE,*mspace,*fdataspace);
  }  catch( DataSetIException error )
  {
    error.printError();
    Rcpp::stop("Error writing file");
  }
  //  std::cout<<"Data sucessfully written"<<std::endl;
  try{
    file->flush(H5F_SCOPE_GLOBAL);
  }catch(FileIException error)
  {
    error.printError();
    Rcpp::stop("Error flushing file");
  }
//    std::cout<<"File flushed"<<std::endl;
  dataset->close();
  fdataspace->close();
  mspace->close();
  group->close();
  file->close();
  fdataspace->close();
  mspace->close();

}




void write_dvec_h5(const std::string h5file, const std::string groupname, const std::string dataname, c_arrayxd_internal data, const int deflate_level){

  size_t datasize=data.size();

  FloatType ftypew(PredType::NATIVE_DOUBLE);
  H5FilePtr file =create_or_open_file(h5file);


  size_t chunksize=1000;
  if(chunksize>datasize){
    chunksize=datasize;
  }

  H5GroupPtr group =create_or_open_group(file,groupname);
  std::vector<hsize_t> cumdim{datasize};
  std::vector<hsize_t> maxdim{datasize};
  std::vector<hsize_t> chunkdim{chunksize};

  H5DataSetPtr dataset =create_or_open_dataset(group,dataname,ftypew,cumdim,maxdim,chunkdim,deflate_level);

  DataSpace* fdataspace;
  try{
    fdataspace= new DataSpace(dataset->getSpace());
  }catch(DataSpaceIException error){
    error.printError();
    Rcpp::stop("Error creating memory dataspace ");
  }

  hsize_t datadim[1];
  fdataspace->getSimpleExtentDims(datadim,NULL);

  hsize_t memdim[]={datasize};

  DataSpace *mspace;
  try{
    mspace= new DataSpace(1,memdim); //Size of first dataset (in memory, can be bigger or smaller than size on disk, depending on how much you're writing)
  }catch(DataSpaceIException error){
    error.printError();
    Rcpp::stop("Error creating memory dataspace ");
  }
  hsize_t odim[]={0};//dimension of each offset (current_chunk*chunksize)
  hsize_t stridea[]={1};
  hsize_t blocka[]={1};


  fdataspace->selectHyperslab( H5S_SELECT_SET, memdim,odim,stridea,blocka);


  //  std::cout<<"Starting to write data"<<std::endl;
  try{
    dataset->write(data.data(),PredType::NATIVE_DOUBLE,*mspace,*fdataspace);
  }  catch( DataSetIException error )
  {
    error.printError();
    Rcpp::stop("Error writing file");
  }
//  std::cout<<"Data sucessfully written"<<std::endl;
  try{
    file->flush(H5F_SCOPE_GLOBAL);
  }catch(FileIException error)
  {
    error.printError();
    Rcpp::stop("Error flushing file");
  }
//  std::cout<<"File flushed"<<std::endl;
  dataset->close();
  fdataspace->close();
  mspace->close();
  group->close();
  file->close();
  fdataspace->close();
  mspace->close();

}


void write_ivec_h5(const std::string h5file, const std::string groupname, const std::string dataname, c_arrayxi_internal data, const int deflate_level){

  size_t datasize=data.size();

  IntType ftypew(PredType::NATIVE_INT);
  H5FilePtr file =create_or_open_file(h5file);

  size_t chunksize=1000;
  if(chunksize>datasize){
    chunksize=datasize;
  }

  H5GroupPtr group =create_or_open_group(file,groupname);
  std::vector<hsize_t> cumdim{datasize};
  std::vector<hsize_t> maxdim{datasize};
  std::vector<hsize_t> chunkdim{chunksize};

  H5DataSetPtr dataset =create_or_open_dataset(group,dataname,ftypew,cumdim,maxdim,chunkdim,deflate_level);

  DataSpace* fdataspace;
  try{
    fdataspace= new DataSpace(dataset->getSpace());
  }catch(DataSpaceIException error){
    error.printError();
    Rcpp::stop("Error creating memory dataspace ");
  }

  hsize_t datadim[1];
  fdataspace->getSimpleExtentDims(datadim,NULL);

  hsize_t memdim[]={datasize};

  DataSpace *mspace;
  try{
    mspace= new DataSpace(1,memdim); //Size of first dataset (in memory, can be bigger or smaller than size on disk, depending on how much you're writing)
  }catch(DataSpaceIException error){
    error.printError();
    Rcpp::stop("Error creating memory dataspace ");
  }
  hsize_t odim[]={0};//dimension of each offset (current_chunk*chunksize)
  hsize_t stridea[]={1};
  hsize_t blocka[]={1};


  fdataspace->selectHyperslab( H5S_SELECT_SET, memdim,odim,stridea,blocka);


  //  std::cout<<"Starting to write data"<<std::endl;
  try{
    dataset->write(data.data(),PredType::NATIVE_INT,*mspace,*fdataspace);
  }  catch( DataSetIException error )
  {
    error.printError();
    Rcpp::stop("Error writing file");
  }
//  std::cout<<"Data sucessfully written"<<std::endl;
  try{
    file->flush(H5F_SCOPE_GLOBAL);
  }catch(FileIException error)
  {
    error.printError();
    Rcpp::stop("Error flushing file");
  }
//  std::cout<<"File flushed"<<std::endl;
  dataset->close();
  fdataspace->close();
  mspace->close();
  group->close();
  file->close();
  fdataspace->close();
  mspace->close();

}


