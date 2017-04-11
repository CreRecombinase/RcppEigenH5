#include "RcppEigenH5.h"


std::string print_dims(c_Matrix_internal mat){
  return(std::to_string(static_cast<long long>(mat.rows()))+"x"+std::to_string(static_cast<long long>(mat.cols())));
}


int get_rownum_h5(const std::string h5file,const std::string groupname, const std::string dataname){
  H5FilePtr file = open_file(h5file);
  H5GroupPtr group= open_group(file,groupname);
  H5DataSetPtr dataset = open_dataset(group,dataname);

  DataSpace fspace =dataset->getSpace();
  hsize_t datadim[]={0,0};
  fspace.getSimpleExtentDims(datadim,NULL);
  int retnum=0;
  if(check_transpose(dataset)){
    retnum=datadim[1];
  }else{
    retnum=datadim[0];
  }
  fspace.close();
  dataset->close();
  file->close();


  return(retnum);
}


int get_colnum_h5(const std::string h5file,const std::string groupname, const std::string dataname){
  H5FilePtr file = open_file(h5file);
  H5GroupPtr group= open_group(file,groupname);
  H5DataSetPtr dataset = open_dataset(group,dataname);



  DataSpace fspace =dataset->getSpace();
  hsize_t datadim[]={0,0};
  fspace.getSimpleExtentDims(datadim,NULL);
  int retnum=0;
  if(check_transpose(dataset)){
    retnum=datadim[0];
  }else{
    retnum=datadim[1];
  }


  fspace.close();
  dataset->close();
  file->close();
//  int retnum=datadim[1];
  return(retnum);
}




std::vector<std::string> getGroups(const std::string h5file){


  H5FilePtr file = open_file(h5file);
  hsize_t objc= file->getNumObjs();
  std::vector<std::string> retvec(objc);

  if(objc!=0){
    for(hsize_t i=0; i<objc;i++){
      std::string tst=file->getObjnameByIdx(i);
      retvec[i]=tst;
    }
  }
  file->close();
  return(retvec);
}



