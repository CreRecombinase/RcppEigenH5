#include "RcppEigenH5.h"




H5DataSetPtr create_or_open_dataset(H5GroupPtr &group,const std::string &dataname, const DataType &data_type,std::vector<hsize_t> &cdatadim,std::vector<hsize_t> &mdatadim,std::vector<hsize_t> &chunkdim,const int deflate_level)
{
  H5::Exception::dontPrint();
  DataSet* dataset;
  DataSpace* fdataspace;
  hsize_t objc= group->getNumObjs();



  char* version;
  char* date;
  int r=0;
  unsigned int cd_values[7];
  if(deflate_level>0){
    r = register_blosc(&version,&date);

    //    printf("Blosc version info: %s (%s)\n", version, date);

    cd_values[4]=deflate_level;
    cd_values[5]=1;
    cd_values[6]=BLOSC_BLOSCLZ;
  }
  bool fdat=false;
  if(objc!=0){
    for(hsize_t i=0; i<objc;i++){
      std::string tst=group->getObjnameByIdx(i);
      if(tst==dataname){
        fdat = true;
      }
    }
  }
  if(fdat){
    try{
      dataset = new DataSet(group->openDataSet(dataname));
    }  catch( DataSetIException error )
    {
      error.printError();
      Rcpp::stop("Error opening dataset");
    }
  }else{
    hsize_t *cumdima = new hsize_t[cdatadim.size()];
    hsize_t *chunkdima = new hsize_t[chunkdim.size()];
    hsize_t *chunksizea = new hsize_t[chunkdim.size()];
    hsize_t*mdima = new hsize_t[mdatadim.size()];

    for(int i=0; i<cdatadim.size();i++){
      cumdima[i]=cdatadim[i];
      mdima[i]=mdatadim[i];
      chunkdima[i]=chunkdim[i];
      hsize_t tchunk=std::min((int)chunkdima[i],1000);
      chunksizea[i]=tchunk;
    }
    try{
      fdataspace= new DataSpace(cdatadim.size(),cumdima,mdima); //Create dataspace for dataset (on disk)
    }catch(DataSpaceIException error)
    {
      error.printError();
      Rcpp::stop("Error creating file dataspace");
    }
    DSetCreatPropList cparms; //Create chunksize file parameters
    cparms.setChunk(chunkdim.size(),chunkdima); //Set chunksize
    if(deflate_level>0){
      cparms.setFilter(FILTER_BLOSC,H5Z_FLAG_OPTIONAL,7,cd_values);
    }else{
      //    std::cout<<"Using DEFLATE for compression"<<std::endl;
      cparms.setDeflate(-deflate_level);
      //    cparms.setFilter(FILTER_BLOSC,H5Z_FLAG_OPTIONAL,7,cd_values);
    }
    try{
      dataset = new DataSet(group->createDataSet(dataname,data_type,*fdataspace,cparms));
    }  catch( DataSetIException error )
    {
      error.printError();
      Rcpp::stop("Error creating dataset");
    }
  }
  return H5DataSetPtr(dataset);
}


H5DataSetPtr open_dataset(H5GroupPtr &group,const std::string &dataname)
{
  H5::Exception::dontPrint();
  DataSet* dataset;
  DataSpace* fdataspace;
  hsize_t objc= group->getNumObjs();
  char* version;
  char* date;
  int r=0;
  //  std::cout<<"registering blosc"<<std::endl;
  r = register_blosc(&version,&date);
  unsigned int cd_values[7];
  //  printf("Blosc version info: %s (%s)\n", version, date);
  // cd_values[4]=deflate_level;
  // cd_values[5]=1;
  // cd_values[6]=BLOSC_BLOSCLZ;
  bool fdat=false;
  if(objc!=0){
    for(hsize_t i=0; i<objc;i++){
      std::string tst=group->getObjnameByIdx(i);
      if(tst==dataname){
        fdat = true;
      }
    }
  }
  if(fdat){
    try{
      dataset = new DataSet(group->openDataSet(dataname));
    }  catch( DataSetIException error )
    {
      error.printError();
      Rcpp::stop("Error opening dataset");
    }
  }else{
    Rcpp::Rcerr<<"Can't find "<<dataname<<std::endl;
    Rcpp::stop("Dataset not found!");
  }
  return H5DataSetPtr(dataset);
}

