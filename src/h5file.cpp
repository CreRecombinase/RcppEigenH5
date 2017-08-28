#include "RcppEigenH5.h"



using namespace H5;

bool f_exists (const std::string name) {
  if (FILE *file = fopen(name.c_str(), "r")) {
    fclose(file);
    return true;
  } else {
    return false;
  }
}


H5FilePtr create_or_open_file(const std::string fname)
{
  H5::H5File* file;

#ifdef BLOSC_BLOSCLZ
  char* version;
  char* date;  int r=0;
  r = register_blosc(&version,&date);
#endif

//  bool fe = file_exists(fname);
  if(!f_exists(fname)){
    try{
      file = new H5::H5File(fname.c_str(), H5F_ACC_EXCL);
    }catch(FileIException error){
      error.printError();
      Rcpp::stop("Error creating file!");
    }
  }
  else{
    try{
      file = new H5::H5File(fname.c_str(), H5F_ACC_RDWR);
    }catch(FileIException error){
      error.printError();
      Rcpp::stop("Error opening file!");
    }
  }
  return H5FilePtr(file);
}

H5FilePtr open_file(const std::string fname)
{
  H5::H5File* file;

#ifdef BLOSC_BLOSCLZ
  char* version;
  char* date;  int r=0;
  r = register_blosc(&version,&date);
#endif


  if(!f_exists(fname)){
    Rcpp::stop("File does not exist!");
  }
  else{
    try{
      file = new H5::H5File(fname.c_str(), H5F_ACC_RDONLY);
    }catch(FileIException error){
      error.printError();
      Rcpp::stop("Error opening file!");
    }
  }
  return H5FilePtr(file);
}
