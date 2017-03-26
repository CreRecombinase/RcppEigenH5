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
  return(strreadbuf);
}


std::string read_group_attr_h5(const std::string h5file, const std::string groupname,const std::string attr_name){
  using namespace H5;
  H5std_string strreadbuf ("");

  H5FilePtr file=open_file(h5file);
  H5GroupPtr group= open_group(file,groupname);

  Attribute attr = group->openAttribute(attr_name);

  DataType strdatatype = attr.getDataType();
  attr.read(strdatatype, strreadbuf);
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

//[[Rcpp::export(name="read_group_attr_h5")]]
Rcpp::CharacterVector read_group_attr_h5_exp(const std::string h5file, const std::string groupname,const std::string attr_name){
  return(Rcpp::wrap(read_group_attr_h5(h5file,groupname,attr_name)));
}


