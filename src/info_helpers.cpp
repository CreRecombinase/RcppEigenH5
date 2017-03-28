#include "RcppEigenH5.h"


int get_rownum_h5(const std::string h5file,const std::string groupname, const std::string dataname){
  H5FilePtr file = open_file(h5file);
  H5GroupPtr group= open_group(file,groupname);
  H5DataSetPtr dataset = open_dataset(group,dataname);
  DataSpace fspace =dataset->getSpace();
  hsize_t datadim[]={0};
  fspace.getSimpleExtentDims(datadim,NULL);
  fspace.close();
  dataset->close();
  file->close();
  int retnum=datadim[0];
  return(retnum);
}


int get_colnum_h5(const std::string h5file,const std::string groupname, const std::string dataname){
  H5FilePtr file = open_file(h5file);
  H5GroupPtr group= open_group(file,groupname);
  H5DataSetPtr dataset = open_dataset(group,dataname);
  DataSpace fspace =dataset->getSpace();
  hsize_t datadim[]={0,0};
  fspace.getSimpleExtentDims(datadim,NULL);
  fspace.close();
  dataset->close();
  file->close();
  int retnum=datadim[1];
  return(retnum);
}




std::vector<std::string> getGroups(const std::string h5file){

  H5File* file;

  try{
    file= new H5File( h5file.c_str(), H5F_ACC_RDONLY);
  }
  catch( FileIException error )
  {
    error.printError();
  }
  // file->getNumObjs()
  // Group group = file->openGroup(groupname);
  hsize_t objc= file->getNumObjs();
  std::vector<std::string> retvec(objc);
  bool fdat=false;
  if(objc!=0){
    for(hsize_t i=0; i<objc;i++){
      std::string tst=file->getObjnameByIdx(i);
      retvec[i]=tst;
    }
  }
  file->close();
  return(retvec);
}


