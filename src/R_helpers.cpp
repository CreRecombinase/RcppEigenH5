#include "RcppEigenH5.h"
using namespace Rcpp;
#include <algorithm>
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <unordered_set>
#include <unordered_map>
#include <unistd.h>
#include <limits.h>

//'Some documentation
//' @param
//' @export
//[[Rcpp::export]]
bool is_transposed_h5(std::string h5file,std::string groupname,std::string dataname){
  H5FilePtr tf=open_file(h5file);
  H5GroupPtr grp = open_group(tf,groupname);
  H5DataSetPtr dst = open_dataset(grp,dataname);
  const bool ret =check_transpose(dst);
  dst->close();
  grp->close();
  tf->close();
  return(ret);
}


//[[Rcpp::export]]
std::string get_selfpath() {
  char buff[PATH_MAX];
  ssize_t len = ::readlink("/proc/self/exe", buff, sizeof(buff)-1);
  if (len != -1) {
    buff[len] = '\0';
    return std::string(buff);
  }
  /* handle error condition */
}

//[[Rcpp::export]]
Rcpp::StringVector h5ls(std::string h5file,std::string base_groupname="/"){

//  H5FilePtr open_file(const std::string fname)
H5FilePtr tf=open_file(h5file);

  H5GroupPtr grp = open_group(tf,base_groupname);
  if(base_groupname!="/"){
    base_groupname=base_groupname+"/";
  }
  Rcpp::StringVector retstring = Rcpp::wrap(h5_ls_dset(grp,base_groupname));

  return(retstring);
}






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


template<typename T>
struct matrix_hash : std::unary_function<T, size_t> {
  std::size_t operator()(T const& matrix) const {
    // Note that it is oblivious to the storage order of Eigen matrix (column- or
    // row-major). It will give you the same hash value for two different matrices if they
    // are the transpose of each other in different storage order.
    size_t seed = 0;
    for (size_t i = 0; i < matrix.size(); ++i) {
      auto elem = *(matrix.data() + i);
      seed ^= std::hash<typename T::Scalar>()(elem) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }
    return seed;
  }
};


//[[Rcpp::export]]
Rcpp::NumericVector calc_af(const std::string h5file,const std::string groupname, const std::string dataname, const Eigen::ArrayXi index, const Rcpp::IntegerVector chunksize,bool display_progress=true,bool check_dup=true){
  using namespace Eigen;
  size_t p=index.size();
  Eigen::ArrayXd retvec(p);
  size_t csize= chunksize[0];
  size_t totchunks=ceil((double) p / (double) csize);
  Eigen::ArrayXi subi;
  size_t rownum =get_rownum_h5(h5file,groupname,dataname);
  Eigen::MatrixXd temp(rownum,csize);
  // Eigen::Matrix<char,Dynamic,Dynamic> tempi(rownum,csize);

  // Rcpp::Rcout<<"totchunks: "<<totchunks<<std::endl;
  std::vector<bool> isDup(p);
  std::unordered_map<ColumnMatrixi,int,matrix_hash<ColumnMatrixi> > my_map;
  // if(check_dup){
  // my_map.reserve(p);
  // }
  // std::unordered_map<Array<char,Dynamic,1>,int >::iterator my_map_it;

  ColumnMatrixi ti=temp.col(0).cast<char>();
  int tti=0;
  Progress pp(totchunks, display_progress);
  for(size_t i=0; i<totchunks;i++){
    size_t chunkstart =i*csize;
    size_t chunkstop =std::min((p-1),((i+1)*csize)-1);
    if (Progress::check_abort() )
      return Rcpp::wrap(retvec);


    size_t tchunksize= chunkstop-chunkstart+1;
    // Rcpp::Rcout<<"Chunk: "<<i<<"of size: "<<tchunksize<<std::endl;
    // Rcpp::Rcout<<index.segment(chunkstart,tchunksize)<<std::endl;

    read_2d_cindex_h5(h5file,groupname,dataname,index.segment(chunkstart,tchunksize),temp);
    // tempi =temp.cast<char>();
    retvec.segment(chunkstart,tchunksize)=temp.colwise().mean()/2;
    if(check_dup){
      for(int c=chunkstart;c<(chunkstart+tchunksize);c++){
        ti = round(temp.col(c-chunkstart).array()).cast<char>();
        tti=my_map[ti]++;
        if(tti!=0){
          retvec(c)=0;
        }
      }
    }
    pp.increment();
  }
  return(Rcpp::wrap(retvec));
}


//[[Rcpp::export]]
Rcpp::DataFrame calc_summ_h5(const std::string h5file,const std::string groupname, const std::string dataname, const Eigen::ArrayXi index, const Rcpp::IntegerVector chunksize,bool display_progress=true,bool check_dup=true){
  using namespace Eigen;
  using namespace Rcpp;
  size_t p=index.size();
  Eigen::ArrayXd retvec(p);
  size_t csize= chunksize[0];
  size_t totchunks=ceil((double) p / (double) csize);
  Eigen::ArrayXi subi;
  Rcpp::IntegerVector rindex(Rcpp::wrap(index));
  size_t rownum =get_rownum_h5(h5file,groupname,dataname);

  Eigen::ArrayXd dat_min(p);
  Eigen::ArrayXd dat_max(p);
  Eigen::ArrayXd dat_af(p);
  Rcpp::LogicalVector dat_isdup(p);


  Eigen::MatrixXd temp(rownum,csize);
  // Eigen::Matrix<char,Dynamic,Dynamic> tempi(rownum,csize);

  // Rcpp::Rcout<<"totchunks: "<<totchunks<<std::endl;
  // std::vector<bool> isDup(p);
  std::unordered_map<ColumnMatrixi,int,matrix_hash<ColumnMatrixi> > my_map;
  // if(check_dup){
  // my_map.reserve(p);
  // }
  // std::unordered_map<Array<char,Dynamic,1>,int >::iterator my_map_it;

  ColumnMatrixi ti=temp.col(0).cast<char>();
  int tti=0;
  Progress pp(totchunks, display_progress);
  for(size_t i=0; i<totchunks;i++){
    size_t chunkstart =i*csize;
    size_t chunkstop =std::min((p-1),((i+1)*csize)-1);
    if (Progress::check_abort() )
      return Rcpp::wrap(retvec);


    size_t tchunksize= chunkstop-chunkstart+1;
    // Rcpp::Rcout<<"Chunk: "<<i<<"of size: "<<tchunksize<<std::endl;
    // Rcpp::Rcout<<index.segment(chunkstart,tchunksize)<<std::endl;

    read_2d_cindex_h5(h5file,groupname,dataname,index.segment(chunkstart,tchunksize),temp);
    // tempi =temp.cast<char>();
    dat_af.segment(chunkstart,tchunksize)=temp.colwise().mean()/2;
    dat_min.segment(chunkstart,tchunksize)=temp.colwise().minCoeff();
    dat_max.segment(chunkstart,tchunksize)=temp.colwise().maxCoeff();


    if(check_dup){
      for(int c=chunkstart;c<(chunkstart+tchunksize);c++){
        ti = temp.col(c-chunkstart).cast<char>();
        tti=my_map[ti]++;
        dat_isdup(c)=(tti!=0);

        // if(tti!=0){
        //
        //   retvec(c)=0;
        // }
      }
    }
    pp.increment();
  }
  return(Rcpp::DataFrame::create(_["af"]=Rcpp::wrap(dat_af),
                                 _["min"]=Rcpp::wrap(dat_min),
                                 _["max"]=Rcpp::wrap(dat_max),
                                 _["isDup"]=dat_isdup,
                                 _["index"]=rindex));

}






Eigen::ArrayXd calc_variance(c_Matrix_internal mat){
  int n=mat.rows();
  return(((mat.rowwise()-(mat.colwise().mean())).array().square().colwise().sum())/(n-1));
}


//[[Rcpp::export]]
Rcpp::NumericVector calc_var(const std::string h5file,const std::string groupname, const std::string dataname, const Eigen::ArrayXi index, const Rcpp::IntegerVector chunksize,bool display_progress=true){

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

    retvec.segment(chunkstart,tchunksize)=calc_variance(temp);
  }
  return(Rcpp::wrap(retvec));
}


Eigen::MatrixXd scale_mat(Matrix_internal mat){
  return((mat.rowwise()-(mat.colwise().mean())));
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

  read_2d_cindex_chunk_h5(h5file,groupname,dataname,index,temp,csize);

  retvec=retvec.matrix()+(temp.rowwise()-(temp.colwise().mean()))*beta;

  return(Rcpp::wrap(retvec));
}

//[[Rcpp::export(name="list_groups_h5")]]
std::vector<std::string> h5ls_grp_exp(const std::string h5file,const std::string base_group="/")
{
  H5FilePtr file=open_file(h5file);
std::vector<std::string> groupNames=list_subgroups(file, base_group);
  file->close();
  return(groupNames);
}



//[[Rcpp::export]]
Rcpp::LogicalVector group_exists(std::string h5file,std::string base_group="/"){
Rcpp::LogicalVector res(1);
  H5FilePtr file=open_file(h5file);
  res(0)=grp_path_exists(file,base_group);
  file->close();
  return(res);
}


//[[Rcpp::export(name="list_attrs_h5")]]
std::vector<std::string> h5ls_attr_exp(std::string h5file,std::string base_group="/")
{
  H5FilePtr file=open_file(h5file);
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



