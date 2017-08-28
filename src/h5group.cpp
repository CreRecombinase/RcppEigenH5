#include "RcppEigenH5.h"
#include <string>
#include <sstream>
#include <vector>
#include <iterator>
#include <algorithm>



template<typename Out>
void split(const std::string &s, char delim, Out result) {
  std::stringstream ss;
  ss.str(s);
  std::string item;
  while (std::getline(ss, item, delim)) {
    *(result++) = item;
  }
}


std::vector<std::string> split_string(const std::string &s, char delim) {
  std::vector<std::string> elems;
  split(s, delim, std::back_inserter(elems));
  return elems;
}




bool grp_exists(const H5GroupPtr grp, const std::string groupname){
  hsize_t objc = grp->getNumObjs();
  for(hsize_t i=0; i<objc;i++){
    std::string tst=grp->getObjnameByIdx(i);
    if(tst==groupname){
      if(grp->getObjTypeByIdx(i)==H5G_GROUP){
        return(true);
      }else{
        return(false);
      }
    }
  }
  return(false);
}


std::vector<std::string> subgrp_grp(const H5GroupPtr trg){
  size_t objc =trg->getNumObjs();
  std::vector<std::string> retvec;
  retvec.reserve(objc);
  for(hsize_t i=0; i<objc;i++){
    if(trg->getObjTypeByIdx(i)==H5G_GROUP){
      retvec.push_back(trg->getObjnameByIdx(i));
    }
  }
  return(retvec);
}




bool grp_path_exists(const H5FilePtr file, const std::string groupname){
  H5GroupPtr trg =std::make_shared<Group>(file->openGroup("/"));
  if(groupname=="/"){
    trg->close();
    return true;
  }
  std::vector<std::string> groupvec=split_string(groupname,'/');
  int grpsize=groupvec.size();
  size_t i=0;

  for( auto subname : groupvec){
    if(!grp_exists(trg,subname)){
      trg->close();
      return false;
    }else{
      trg=std::make_shared<Group>(trg->openGroup(subname));
    }
  }
  trg->close();
  return true;
}


H5GroupPtr create_or_open_group(H5FilePtr file,const std::string groupname)
{

  H5GroupPtr trg =std::make_shared<Group>(file->openGroup("/"));
  if(groupname=="/"){
    return trg;
  }
  std::vector<std::string> groupvec=split_string(groupname,'/');
  int grpsize=groupvec.size();
  size_t i=0;
  for( auto subname : groupvec){

    if(!grp_exists(trg,subname)){
      try{
        trg = std::make_shared<Group>(trg->createGroup(subname));
        //        trg = H5GroupPtr(group);
      }catch(GroupIException error){
        error.printError();
        Rcpp::stop("Error creating group");
      }
    }else{
      try{
        trg = std::make_shared<Group>(trg->openGroup(subname));
        //        trg = H5GroupPtr(group);
      }catch(GroupIException error){
        error.printError();
        Rcpp::stop("Error opening group");
      }
    }
  }
  return(trg);
}







H5GroupPtr open_group(const H5FilePtr file,const std::string groupname)
{
  H5GroupPtr trg =std::make_shared<Group>(file->openGroup("/"));
  if(groupname=="/"){
    return trg;
  }
  std::vector<std::string> groupvec=split_string(groupname,'/');
  int grpsize=groupvec.size();
  size_t i=0;

  for( auto subname : groupvec){
    std::vector<std::string> subgrps= subgrp_grp(trg);
    if((std::find(subgrps.begin(), subgrps.end(), subname) == subgrps.end())||subgrps.size()==0){
      Rcpp::Rcerr<<"Unable to find group "<<groupname<<std::endl;
      Rcpp::stop("Group doesn't exist!");
    }else{
      try{
        trg = std::make_shared<Group>(trg->openGroup(subname));
        //        trg = H5GroupPtr(group);
      }catch(GroupIException error){
        error.printError();
        Rcpp::Rcerr<<"Unable to open group "<<groupname<<std::endl;
        Rcpp::stop("Error opening group");
      }
    }


  }
  return trg;
}






