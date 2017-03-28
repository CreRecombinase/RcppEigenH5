#include "RcppEigenH5.h"
#include "H5Cpp.h"

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


//[[Rcpp::export(name="read_data_attr_h5")]]
Rcpp::CharacterVector read_data_attr_h5_exp(const StringVector h5filename, const StringVector h5_groupname, const StringVector h5_dataname,const StringVector h5_attr_name){
  std::string h5file(h5filename[0]);
  std::string groupname(h5_groupname[0]);
  std::string dataname(h5_dataname[0]);
  std::string attr_name(h5_attr_name[0]);

  return(Rcpp::wrap(read_data_attr_h5(h5file,groupname,dataname,attr_name)));
}


//[[Rcpp::export(name="write_data_string_attr_h5")]]
void write_data_string_attr_h5_exp(const StringVector h5filename, const StringVector h5_groupname, const StringVector h5_dataname,const StringVector h5_attr_name, const StringVector h5_attr_value){
  std::string h5file(h5filename[0]);
  std::string groupname(h5_groupname[0]);
  std::string dataname(h5_dataname[0]);
  std::string attr_name(h5_attr_name[0]);
  std::string attr_value(h5_attr_value[0]);

  write_data_string_attr_h5(h5file,groupname,dataname,attr_name,attr_value);
}

//[[Rcpp::export(name="write_group_string_attr_h5")]]
void write_group_string_attr_h5_exp(const StringVector h5filename, const StringVector h5_groupname,const StringVector h5_attr_name, const StringVector h5_attr_value){
  std::string h5file(h5filename[0]);
  std::string groupname(h5_groupname[0]);
  std::string attr_name(h5_attr_name[0]);
  std::string attr_value(h5_attr_value[0]);

  write_group_string_attr_h5(h5file,groupname,attr_name,attr_value);

}


//[[Rcpp::export(name="write_data_int_attr_h5")]]
void write_data_int_attr_h5_exp(const StringVector h5filename, const StringVector h5_groupname, const StringVector h5_dataname,const StringVector h5_attr_name, const IntegerVector h5_attr_value){
  std::string h5file(h5filename[0]);
  std::string groupname(h5_groupname[0]);
  std::string dataname(h5_dataname[0]);
  std::string attr_name(h5_attr_name[0]);
  int attr_value(h5_attr_value[0]);

  write_data_int_attr_h5(h5file,groupname,dataname,attr_name,attr_value);
}

//[[Rcpp::export(name="write_group_int_attr_h5")]]
void write_group_int_attr_h5_exp(const StringVector h5filename, const StringVector h5_groupname,const StringVector h5_attr_name, const IntegerVector h5_attr_value){
  std::string h5file(h5filename[0]);
  std::string groupname(h5_groupname[0]);
  std::string attr_name(h5_attr_name[0]);
  int attr_value(h5_attr_value[0]);

  write_group_int_attr_h5(h5file,groupname,attr_name,attr_value);

}


//[[Rcpp::export(name="read_group_attr_h5")]]
Rcpp::CharacterVector read_group_attr_h5_exp(const std::string h5file, const std::string groupname,const std::string attr_name){
  return(Rcpp::wrap(read_group_attr_h5(h5file,groupname,attr_name)));
}


