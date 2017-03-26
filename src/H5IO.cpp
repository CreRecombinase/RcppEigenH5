#include <RcppEigenH5.h>
#include <valarray>
//[[Rcpp::depends(RcppEigen)]]
using namespace Eigen;

void read_2ddmat_h5(const std::string h5file, const std::string groupname, const std::string dataname,const size_t row_offset,const size_t col_offset,const size_t row_chunksize,const size_t col_chunksize,  double* data){

  //Try breaking up reads in to chunks
  H5FilePtr file=open_file(h5file);
  H5GroupPtr group= open_group(file,groupname);
  H5DataSetPtr dataset = open_dataset(group,dataname);


   // Rcpp::Rcout<<"row_chunksize is :"<<row_chunksize<<std::endl;
   // Rcpp::Rcout<<"col_chunksize is :"<<col_chunksize<<std::endl;
  DataType dt= dataset->getDataType();
  hsize_t datadims[]={0,0};
  DataSpace fspace=dataset->getSpace();
  fspace.getSimpleExtentDims(datadims,NULL);
  //  std::cout<<"Full data is of dimensions"<<datadims[0]<<"x"<<datadims[1]<<std::endl;
  hsize_t matrix_dims[2];
  if(row_offset+row_chunksize>datadims[0]){
    Rcpp::stop("row_offset ("+std::to_string(row_offset)+") + row_chunksize ("+std::to_string(row_chunksize)+") >datadims[0] ("+std::to_string(datadims[0])+")");
  }
  if(col_offset+col_chunksize>datadims[1]){
    Rcpp::stop("col_offset ("+std::to_string(col_offset)+") + col_chunksize ("+std::to_string(col_chunksize)+") >datadims[1] ("+std::to_string(datadims[1])+")");
  }

  // std::cout<<"read consists of"<<colchunknum*rowchunknum<<"chunks"<<std::endl;
  // size_t slab_bytes= rchunksize*cchunksize*sizeof(float);
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
  dataset->read(data,dt,memspace,fspace);
  //  std::cout<<"Read complete!"<<std::endl;

  dt.close();
  memspace.close();
  fspace.close();
  group->close();
  dataset->close();
  file->close();
}



Eigen::MatrixXd read_2d_h5(const std::string h5file, const std::string groupname, const std::string dataname,const  Eigen::ArrayXi offset ,const  Eigen::ArrayXi chunksize){

  size_t row_offset=offset[0];
  size_t col_offset=offset[1];

  size_t row_chunksize=chunksize[0];
  size_t col_chunksize=chunksize[1];

  size_t tot_size = row_chunksize*col_chunksize;
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> readmat(row_chunksize,col_chunksize);
  read_2ddmat_h5(h5file,groupname,dataname,row_offset,col_offset,row_chunksize,col_chunksize,readmat.data());
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> retmat=readmat;
  return(retmat);
}

Eigen::MatrixXd read_2d_h5(const std::string h5file, const std::string groupname, const std::string dataname){

  size_t row_offset=0;
  size_t col_offset=0;

  size_t row_chunksize=get_rownum_h5(h5file,groupname,dataname);
  size_t col_chunksize=get_colnum_h5(h5file,groupname,dataname);

  size_t tot_size = row_chunksize*col_chunksize;
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> readmat(row_chunksize,col_chunksize);
  read_2ddmat_h5(h5file,groupname,dataname,row_offset,col_offset,row_chunksize,col_chunksize,readmat.data());
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> retmat=readmat;
  return(retmat);
}




Eigen::ArrayXd read_dvec_h5(const std::string h5file, const std::string groupname, const std::string dataname,const size_t offset,const size_t chunksize){

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
  size_t offset=0;
  size_t chunksize=datadims[0];

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
  size_t offset=0;
  size_t chunksize=datadims[0];
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


Eigen::ArrayXi read_ivec_h5(const std::string h5file, const std::string groupname, const std::string dataname,const size_t offset,const size_t chunksize){

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









void read_ivec_h5(const std::string h5file, const std::string groupname, const std::string dataname,const size_t offset,const size_t chunksize, int* data){

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
  size_t offset=0;
  size_t chunksize=datadims[0];
  Eigen::ArrayXi arrayxi(chunksize);

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
  dataset->read(arrayxi.data(),dt,memspace,fspace);
  //  std::cout<<"Read complete!"<<std::endl;
  dt.close();
  memspace.close();
  fspace.close();
  group->close();
  dataset->close();
  file->close();
  return(arrayxi);
}





void read_2d_cindex_h5(const std::string h5file,const std::string groupname, const std::string dataname, const  Eigen::ArrayXi indvec, Matrix_internal retmat){

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

  Eigen::ArrayXi rowind(rownum);

  rowind.setLinSpaced(rownum,0,rownum-1);

  const Eigen::ArrayXi newind=indvec-minind;


  //Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> readmat(rownum,t_chunksize);
  double *x =(double* )malloc(rownum*t_chunksize*sizeof(double));
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
//  retmat=indexing(readmat,rowind,newind);
  //    free(x);
  // Rcpp::Rcout<<"Sum of x is "<<retmat.sum()<<std::endl;

  //    return(retmat);
}



Eigen::ArrayXd read_1d_cindex_h5(const std::string h5file,const std::string groupname, const std::string dataname, const  Eigen::ArrayXi indvec){

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





