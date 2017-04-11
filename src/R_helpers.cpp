#include "RcppEigenH5.h"
using namespace Rcpp;
#include <algorithm>
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>








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



//[[Rcpp::export]]
Rcpp::NumericVector calc_af(const std::string h5file,const std::string groupname, const std::string dataname, const Eigen::ArrayXi index, const Rcpp::IntegerVector chunksize,bool display_progress=true){

  size_t p=index.size();
  Eigen::ArrayXd retvec(p);
  size_t csize= chunksize[0];
  size_t totchunks=ceil((double) p / (double) csize);
  Eigen::ArrayXi subi;
  size_t rownum =get_rownum_h5(h5file,groupname,dataname);
  Eigen::MatrixXd temp(rownum,csize);
  // Rcpp::Rcout<<"totchunks: "<<totchunks<<std::endl;

  Progress pp(totchunks, display_progress);
  for(size_t i=0; i<totchunks;i++){
    size_t chunkstart =i*csize;
    size_t chunkstop =std::min((p-1),((i+1)*csize)-1);
    if (Progress::check_abort() )
      return Rcpp::wrap(retvec);

    pp.increment();
    size_t tchunksize= chunkstop-chunkstart+1;
    // Rcpp::Rcout<<"Chunk: "<<i<<"of size: "<<tchunksize<<std::endl;
    // Rcpp::Rcout<<index.segment(chunkstart,tchunksize)<<std::endl;

    read_2d_cindex_h5(h5file,groupname,dataname,index.segment(chunkstart,tchunksize),temp);
    retvec.segment(chunkstart,tchunksize)=temp.colwise().mean()/2;

  }
  return(Rcpp::wrap(retvec));
}



//[[Rcpp::export]]
Rcpp::NumericVector calc_yh(const std::string h5file,const std::string groupname, const std::string dataname, const Eigen::ArrayXi index, const Eigen::VectorXd beta,const Rcpp::IntegerVector chunksize,bool display_progress=true){

  size_t p=index.size();
  if(beta.size()!=p){
    Rcpp::stop("p!=beta.size()");
  }

  size_t csize= chunksize[0];
  size_t totchunks=ceil((double) p / (double) csize);

  size_t rownum =get_rownum_h5(h5file,groupname,dataname);
  Eigen::ArrayXd retvec(rownum);
  retvec.setZero();
  Eigen::MatrixXd temp(rownum,p);
  // Rcpp::Rcout<<"totchunks: "<<totchunks<<std::endl;

  // Progress pp(totchunks, display_progress);
  // for(size_t i=0; i<totchunks;i++){
    // size_t chunkstart =i*csize;
    // size_t chunkstop =std::min((p-1),((i+1)*csize)-1);
    // if (Progress::check_abort() )
    //   return Rcpp::wrap(retvec);
    //
    // pp.increment();
    // size_t tchunksize= chunkstop-chunkstart+1;
    // Rcpp::Rcout<<"Chunk: "<<i<<"of size: "<<tchunksize<<std::endl;
    // Rcpp::Rcout<<index.segment(chunkstart,tchunksize)<<std::endl;

    read_2d_cindex_chunk_h5(h5file,groupname,dataname,index,temp,csize);
    retvec=retvec.matrix()+temp*beta;

  return(Rcpp::wrap(retvec));
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

//[[Rcpp::export(name="read_group_iarray_attr_h5")]]
Rcpp::IntegerVector read_group_iarray_attr_h5_exp(const std::string h5file, const std::string groupname,const std::string attr_name){
  return(Rcpp::wrap(read_group_iarray_attr_h5(h5file,groupname,attr_name)));
}

//[[Rcpp::export(name="read_data_iarray_attr_h5")]]
Rcpp::IntegerVector read_data_iarray_attr_h5_exp(const std::string h5file, const std::string groupname,const std::string dataname, const std::string attr_name){
  return(Rcpp::wrap(read_data_iarray_attr_h5(h5file,groupname,dataname,attr_name)));
}



//[[Rcpp::export(name="read_igroup_attr_h5")]]
Rcpp::IntegerVector read_igroup_attr_h5_exp(const std::string h5file, const std::string groupname,const std::string attr_name){
  return(Rcpp::wrap(read_igroup_attr_h5(h5file,groupname,attr_name)));
}

//[[Rcpp::export(name="read_idata_attr_h5")]]
Rcpp::IntegerVector read_idata_attr_h5_exp(const std::string h5file, const std::string groupname,const std::string dataname, const std::string attr_name){
  return(Rcpp::wrap(read_idata_attr_h5(h5file,groupname,dataname,attr_name)));
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



