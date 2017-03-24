#include "RcppEigenH5.h"



H5GroupPtr create_or_open_group(H5FilePtr &file,const std::string &groupname)
{
  Group* group;
  Group *rg = new Group(file->openGroup("/"));
  if(groupname==""||groupname=="/"){
    return H5GroupPtr(rg);
  }
  hsize_t objc= rg->getNumObjs();
  bool fgroup=false;
  if(objc!=0){
    for(hsize_t i=0; i<objc;i++){
      std::string tst=rg->getObjnameByIdx(i);
      if(tst==groupname){
        fgroup = true;
      }
    }
  }
  if(fgroup){
    try{
      group = new Group(file->openGroup(groupname));
    }catch(GroupIException error){
      error.printError();
      Rcpp::stop("Error opening group");
    }
  }else{
    try{
      group = new Group(file->createGroup("/"+groupname));
    }catch(GroupIException error){
      error.printError();
      Rcpp::stop("Error creating group");
    }
  }
  return H5GroupPtr(group);
}

H5GroupPtr open_group(H5FilePtr &file,const std::string &groupname)
{
  Group* group;
  Group *rg;
  try{
    rg= new Group(file->openGroup("/"));
  }catch(GroupIException error){
    error.printError();
    Rcpp::stop("Error opening group");
  }
  if(groupname==""||groupname=="/"){
    return H5GroupPtr(rg);
  }
  hsize_t objc= rg->getNumObjs();
  bool fgroup=false;
  if(objc!=0){
    for(hsize_t i=0; i<objc;i++){
      std::string tst=rg->getObjnameByIdx(i);
      if(tst==groupname){
        fgroup = true;
      }
    }
  }
  if(fgroup){
    try{
      group = new Group(file->openGroup(groupname));
    }catch(GroupIException error){
      error.printError();
      Rcpp::stop("Error opening group");
    }
  }else{
    Rcpp::stop("Group does not exist!");
  }
  return H5GroupPtr(group);
}




