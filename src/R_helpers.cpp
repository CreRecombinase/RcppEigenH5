#include "RcppEigenH5.h"
#include <set>
using namespace Rcpp;
//[[Rcpp::export(name="get_rownum_h5")]]
Rcpp::IntegerVector h5_rownum(const std::string h5file,const std::string groupname, const std::string dataname){

  Rcpp::IntegerVector retvec(1);
  retvec[0]=get_rownum_h5(h5file,groupname,dataname);
  return(retvec);
}

//[[Rcpp::export(name="get_colnum_h5")]]
Rcpp::IntegerVector h5_colnum(const std::string h5file,const std::string groupname, const std::string dataname){

  Rcpp::IntegerVector retvec(1);
  retvec[0]=get_colnum_h5(h5file,groupname,dataname);
  return(retvec);
}

//[[Rcpp::export(name="list_groups_h5")]]
std::vector<std::string> h5ls_grp_exp(std::string h5file,std::string base_group="/")
{
  H5FilePtr file=create_or_open_file(h5file);
  Group* group;
  Group *rg = new Group(file->openGroup(base_group));
  size_t num_grp_attrs=0;
  std::set<std::string> attr_names;
  hsize_t objc= rg->getNumObjs();
  std::vector<std::string> groupNames(objc);
  if(objc!=0){
    for(hsize_t i=0; i<objc;i++){
      std::string tst=rg->getObjnameByIdx(i);
      groupNames[i]=tst;
    }
  }
  file->close();
  return(groupNames);
}

//[[Rcpp::export(name="list_attrs_h5")]]
std::vector<std::string> h5ls_attr_exp(std::string h5file,std::string base_group="/")
{
  H5FilePtr file=create_or_open_file(h5file);
  Group* group;
  Group *rg = new Group(file->openGroup(base_group));
  size_t num_grp_attrs=0;
  std::set<std::string> attr_names;
  hsize_t objc= rg->getNumAttrs();
  std::vector<std::string> attrNames(objc);
  if(objc!=0){
    for(hsize_t i=0; i<objc;i++){
      Attribute ta =rg->openAttribute(i);
      std::string tst=  ta.getName();
      attrNames[i]=tst;
    }
  }
  file->close();
  return(attrNames);
}

//[[Rcpp::export(name="get_h5_version")]]
StringVector get_h5_version_exp(){

  StringVector ret(1);
  unsigned int majnum=0;
  unsigned int minnum=0;
  unsigned int relnum=0;
  // H5::H5Library::initH5cpp();

  H5::H5Library::getLibVersion(majnum,minnum,relnum);

  std::string ret_string=std::to_string(static_cast<long long>(majnum))+"."+std::to_string(static_cast<long long>(minnum))+"."+std::to_string(static_cast<long long>(relnum));
  ret[0]=ret_string;
  return(ret);
}





