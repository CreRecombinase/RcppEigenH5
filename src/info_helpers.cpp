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






//[[Rcpp::export]]
std::vector<std::string> rcpsplit(const std::string s, const std::string delim) {
  std::vector<std::string> elems;
  const char tdelim=delim.c_str()[0];
  return split_string(s, tdelim);
}
//

//Starting at some base group, recursively(?) list all datasets in the group (and subgroups)
std::vector<std::string>h5_ls_dset(const H5GroupPtr grp,const std::string basename){

  size_t objc =grp->getNumObjs();
  std::vector<std::string> dsets(0);
  for(hsize_t i=0; i<objc;i++){
    if(grp->getObjTypeByIdx(i)==H5G_GROUP){
      std::string grpname=grp->getObjnameByIdx(i);
      std::string nbasename=basename+grpname+"/";
      H5GroupPtr trg=std::make_shared<Group>(grp->openGroup(grpname));
      std::vector<std::string> tdsets=h5_ls_dset(trg,nbasename);
      dsets.insert(dsets.end(),tdsets.begin(),tdsets.end());
      trg->close();
    }else{
      if(grp->getObjTypeByIdx(i)==H5G_DATASET){
        std::string dsname=grp->getObjnameByIdx(i);
        dsets.push_back(basename+dsname);
      }
    }
  }
  return(dsets);
}

std::vector<std::string> list_subgroups(const H5FilePtr file,const std::string base){


//  H5GroupPtr trg=get_groups(file,base);
  std::vector<std::string> groupvec=split_string(base,'/');
  size_t groupsize=groupvec.size();
  std::vector<std::string> retvec;

  if(grp_path_exists(file,base)){
    H5GroupPtr trg = open_group(file,base);
    return(subgrp_grp(trg));
  }
  return(retvec);
}

using namespace H5;

enum DTYPE { T_DOUBLE, T_INTEGER, T_LOGICAL, T_CHARACTER, T_VLEN_FLOAT,
             T_VLEN_DOUBLE, T_VLEN_INTEGER, T_VLEN_LOGICAL, T_COMPOUND, T_DATETIME, T_ENUM};

#define CPTR(VAR,CONST) ((VAR)=(CONST),&(VAR))

DataType GetDataType(const DTYPE datatype, int size = -1) {
  switch(datatype){
  case T_DOUBLE: return PredType::NATIVE_DOUBLE;
  case T_INTEGER: return PredType::NATIVE_INT32;
  case T_LOGICAL: {
    int val;
    EnumType boolenumtype = EnumType(sizeof(char));
    boolenumtype.insert("FALSE", CPTR(val, FALSE));
    boolenumtype.insert("TRUE", CPTR(val, TRUE));
    boolenumtype.insert("NA", CPTR(val, -1));
    return boolenumtype;
  }
  case T_CHARACTER: {
    if ( (size_t)size == H5T_VARIABLE ) { // Assume Variable string size
    StrType datatype(0, H5T_VARIABLE);
    return datatype;
  } else {
    PredType datatype = PredType::C_S1;
    datatype.setSize(size);
    return datatype;
  }
  }
  case T_VLEN_FLOAT: {
    DataType type = PredType::NATIVE_DOUBLE;
    return VarLenType(&type);
  }
  case T_VLEN_DOUBLE: {
    DataType type = PredType::NATIVE_DOUBLE;
    return VarLenType(&type);
  }
  case T_VLEN_INTEGER: {
    DataType type = PredType::NATIVE_INT32;
    return VarLenType(&type);
  }
  case T_VLEN_LOGICAL: {
    DataType type = GetDataType(T_VLEN_LOGICAL);
    return VarLenType(&type);
  }
  case T_COMPOUND:
    throw Rcpp::exception("Writing of compound datatypes is not yet supported.");
  case T_DATETIME:
    throw Rcpp::exception("Writing of date/time datatypes is not yet supported.");
  case T_ENUM:
    throw Rcpp::exception("Writing of enum datatypes is not yet supported.");
  default: throw Rcpp::exception("Unknown data type.");
  }
}

struct cmpDataType {
  bool operator()(const hid_t& a, const hid_t& b) const {
    return H5Tequal(a, b);
  }
};

DTYPE GetTypechar(const DataType &dtype) {
  if ( (dtype == PredType::NATIVE_FLOAT) ||
       (dtype == PredType::NATIVE_DOUBLE) ||
       (dtype == PredType::NATIVE_INT64) ||
       (dtype == PredType::NATIVE_UINT32) ||
       (dtype == PredType::NATIVE_UINT64) ||
       (dtype == PredType::IEEE_F32BE) ||
       (dtype == PredType::IEEE_F32LE) ||
       (dtype == PredType::IEEE_F64BE) ||
       (dtype == PredType::IEEE_F64LE)
  ) {
    return T_DOUBLE;
  }

  if( (dtype == PredType::NATIVE_INT) ||
      (dtype == PredType::NATIVE_INT8) ||
      (dtype == PredType::NATIVE_INT16) ||
      (dtype == PredType::NATIVE_INT32) ||
      (dtype == PredType::NATIVE_UINT8) ||
      (dtype == PredType::NATIVE_UINT16) ||
      (dtype == PredType::STD_U8BE) ||
      (dtype == PredType::STD_U8LE)) {
    return T_INTEGER;
  }
  if (dtype == PredType::C_S1 || dtype.getClass() == H5T_STRING) {
    return T_CHARACTER;
  }
  if (dtype == GetDataType(T_LOGICAL)) {
    return T_LOGICAL;
  }
  if ( (dtype == VarLenType(&PredType::NATIVE_FLOAT)) ||
       (dtype == VarLenType(&PredType::NATIVE_DOUBLE)) ||
       (dtype == VarLenType(&PredType::NATIVE_INT64)) ||
       (dtype == VarLenType(&PredType::NATIVE_UINT32)) ||
       (dtype == VarLenType(&PredType::NATIVE_UINT64))) {
    return T_VLEN_DOUBLE;
  }
  if ( (dtype == VarLenType(&PredType::NATIVE_INT)) ||
       (dtype == VarLenType(&PredType::NATIVE_INT8)) ||
       (dtype == VarLenType(&PredType::NATIVE_INT16)) ||
       (dtype == VarLenType(&PredType::NATIVE_INT32)) ||
       (dtype == VarLenType(&PredType::NATIVE_UINT8)) ||
       (dtype == VarLenType(&PredType::NATIVE_UINT16)) ) {
    return T_VLEN_INTEGER;
  }
  if (dtype.getClass() == H5T_COMPOUND) {
    return T_COMPOUND;
  } else if (dtype.getClass() == H5T_TIME) {
    return T_DATETIME;
  } else if (dtype.getClass() == H5T_ENUM) {
    return T_ENUM;
  }

  /*
   if (dtype == GetDataType(T_VLEN_LOGICAL)) {
   return T_VLEN_LOGICAL;
   } */

  throw Rcpp::exception("Datatype unknown.");
}

DTYPE GetTypechar(char typechar) {
  switch(typechar) {
  case 'd': return T_DOUBLE;
  case 'i': return T_INTEGER;
  case 'l': return T_LOGICAL;
  case 'c': return T_CHARACTER;
  case 'x': return T_VLEN_DOUBLE;
  case 'y': return T_VLEN_INTEGER;
  case 'z': return T_VLEN_LOGICAL;
  case 't': return T_COMPOUND;
  case 'm': return T_DATETIME;
  case 'f': return T_ENUM;
  default: throw new Exception("Typechar unknown");
  }
}

char GetTypechar(DTYPE typechar) {
  switch(typechar) {
  case T_DOUBLE: return 'd';
  case T_INTEGER: return 'i';
  case T_LOGICAL: return 'l';
  case T_CHARACTER: return 'c';
  case T_VLEN_DOUBLE: return 'x';
  case T_VLEN_INTEGER: return 'y';
  case T_VLEN_LOGICAL: return 'z';
  case T_COMPOUND: return 't';
  case T_DATETIME: return 'm';
  case T_ENUM: return 'f';
  default: throw new Exception("Typechar unknown");
  }
}


