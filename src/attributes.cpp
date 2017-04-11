#include "RcppEigenH5.h"


int write_transpose(H5DataSetPtr dataset, bool doTranspose){
  IntType intdatatype(PredType::NATIVE_INT);
  DataSpace att_space(H5S_SCALAR);
  int transpose_attr=0;
  if(doTranspose){
    transpose_attr=1;
  }
  Attribute attr = dataset->createAttribute("doTranspose",intdatatype,att_space);
  attr.write(intdatatype,&transpose_attr);
  attr.close();
  att_space.close();
  return(transpose_attr);
}


int get_int_attr(H5DataSetPtr dataset,std::string attr_name){
  Attribute attr = dataset->openAttribute(attr_name);
  int ireadbuf=0;
  DataType rdatatype = attr.getDataType();
  attr.read(rdatatype, &ireadbuf);
  attr.close();
  return(ireadbuf);
}


bool check_transpose(H5DataSetPtr dataset){

  if(!dataset->attrExists("doTranspose")){
    return(false);
  }
  return(get_int_attr(dataset,"doTranspose")==1);
}




std::string read_data_attr_h5(const std::string h5file, const std::string groupname, const std::string dataname,const std::string attr_name){
  using namespace H5;
  H5std_string strreadbuf ("");

  H5FilePtr file=open_file(h5file);
  H5GroupPtr group= open_group(file,groupname);
  H5DataSetPtr dataset = open_dataset(group,dataname);
  Attribute attr = dataset->openAttribute(attr_name);


  DataType strdatatype = attr.getDataType();
  attr.read(strdatatype, strreadbuf);
  dataset->close();
  group->close();
  file->close();

  return(strreadbuf);
}

int read_idata_attr_h5(const std::string h5file, const std::string groupname, const std::string dataname,const std::string attr_name){
  using namespace H5;
  int ireadbuf=0;

  H5FilePtr file=open_file(h5file);
  H5GroupPtr group= open_group(file,groupname);
  H5DataSetPtr dataset = open_dataset(group,dataname);
  Attribute attr = dataset->openAttribute(attr_name);


  DataType rdatatype = attr.getDataType();
  attr.read(rdatatype, &ireadbuf);
  dataset->close();
  group->close();
  file->close();
  return(ireadbuf);
}


int read_igroup_attr_h5(const std::string h5file, const std::string groupname,const std::string attr_name){
  using namespace H5;
  int ireadbuf=0;

  H5FilePtr file=open_file(h5file);
  H5GroupPtr group= open_group(file,groupname);
  Attribute attr = group->openAttribute(attr_name);


  DataType rdatatype = attr.getDataType();
  attr.read(rdatatype, &ireadbuf);
  group->close();
  file->close();
  return(ireadbuf);
}


std::string read_group_attr_h5(const std::string h5file, const std::string groupname,const std::string attr_name){
  using namespace H5;
  H5std_string strreadbuf ("");

  H5FilePtr file=open_file(h5file);
  H5GroupPtr group= open_group(file,groupname);

  Attribute attr = group->openAttribute(attr_name);

  DataType strdatatype = attr.getDataType();
  attr.read(strdatatype, strreadbuf);
  group->close();
  file->close();
  return(strreadbuf);
}

std::vector<int> read_data_iarray_attr_h5(const std::string h5file, const std::string groupname, const std::string dataname,const std::string attrname){
  using namespace H5;

  H5FilePtr file=open_file(h5file);
  H5GroupPtr group= open_group(file,groupname);
  H5DataSetPtr dataset = open_dataset(group,dataname);

  Attribute attr = dataset->openAttribute(attrname);
  DataSpace fspace=attr.getSpace();
  hsize_t datadims[]={0};
  fspace.getSimpleExtentDims(datadims,NULL);
  //  std::cout<<"Full data is of dimensions"<<datadims[0]<<"x"<<datadims[1]<<std::endl;
  DataType adatat = attr.getDataType();

  std::vector<int> retvec(datadims[0]);
  attr.read(adatat, &retvec[0]);
  fspace.close();
  group->close();
  file->close();
  return(retvec);


}

std::vector<int> read_group_iarray_attr_h5(const std::string h5file, const std::string groupname,const std::string attrname){
  using namespace H5;

  H5FilePtr file=open_file(h5file);
  H5GroupPtr group= open_group(file,groupname);

  Attribute attr = group->openAttribute(attrname);

  DataSpace fspace=attr.getSpace();
  hsize_t datadims[]={0};
  fspace.getSimpleExtentDims(datadims,NULL);
  //  std::cout<<"Full data is of dimensions"<<datadims[0]<<"x"<<datadims[1]<<std::endl;
  DataType adatat = attr.getDataType();

  std::vector<int> retvec(datadims[0]);
  attr.read(adatat, &retvec[0]);
  fspace.close();
  group->close();
  file->close();
  return(retvec);
}








void write_data_string_attr_h5(const std::string h5file, const std::string groupname, const std::string dataname,const std::string attr_name,const std::string attr_value){
  using namespace H5;
  H5std_string strreadbuf ("");
  H5FilePtr file=open_file(h5file);
  H5GroupPtr group= open_group(file,groupname);
  H5DataSetPtr dataset = open_dataset(group,dataname);
  StrType strdatatype(0,H5T_VARIABLE);
  DataSpace att_space(H5S_SCALAR);
  Attribute attr = dataset->createAttribute(attr_name,strdatatype,att_space);
  attr.write(strdatatype,attr_value);
  dataset->close();
  group->close();
  file->close();
}

void write_group_string_attr_h5(const std::string h5file, const std::string groupname,const std::string attr_name,const std::string attr_value){
  using namespace H5;
  H5std_string strreadbuf ("");
  H5FilePtr file=open_file(h5file);
  H5GroupPtr group= open_group(file,groupname);
  StrType strdatatype(0,H5T_VARIABLE);
  DataSpace att_space(H5S_SCALAR);
  Attribute attr = group->createAttribute(attr_name,strdatatype,att_space);
  attr.write(strdatatype,attr_value);
  group->close();
  file->close();
}


void write_group_int_attr_h5(const std::string h5file, const std::string groupname,const std::string attr_name,const int attr_value){
  using namespace H5;
  H5std_string strreadbuf ("");
  H5FilePtr file=open_file(h5file);
  H5GroupPtr group= open_group(file,groupname);
  IntType intdatatype(PredType::NATIVE_INT);
  DataSpace att_space(H5S_SCALAR);
  Attribute attr = group->createAttribute(attr_name,intdatatype,att_space);
  attr.write(intdatatype,&attr_value);
  group->close();
  file->close();
}

void write_data_int_attr_h5(const std::string h5file, const std::string groupname, const std::string dataname,const std::string attr_name,const int attr_value){
  using namespace H5;
  H5std_string strreadbuf ("");
  H5FilePtr file=open_file(h5file);
  H5GroupPtr group= open_group(file,groupname);
  H5DataSetPtr dataset = open_dataset(group,dataname);
  IntType intdatatype(PredType::NATIVE_INT);
  DataSpace att_space(H5S_SCALAR);
  Attribute attr = dataset->createAttribute(attr_name,intdatatype,att_space);
  attr.write(intdatatype,&attr_value);
  dataset->close();
  group->close();
  file->close();
}








