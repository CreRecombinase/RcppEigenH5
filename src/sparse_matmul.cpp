#include "RcppEigenH5.h"

// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::depends(RcppParallel)]]
//[[Rcpp::depends(BH)]]
#include <boost/multi_array.hpp>
#include <RcppParallel.h>
#include <progress.hpp>

//[[Rcpp::depends(RcppEigen)]]
using namespace Eigen;


//[[Rcpp::export]]
Rcpp::NumericMatrix sparse_matmul(StringVector h5files, StringVector groupnames,StringVector datanames,Rcpp::NumericMatrix &xmat,bool do_transpose=true){



  std::vector<std::string> tfilenames=Rcpp::as<std::vector<std::string> >(h5files);
  std::vector<std::string> tgroupnames=Rcpp::as<std::vector<std::string> >(groupnames);
  std::vector<std::string> tdatanames=Rcpp::as<std::vector<std::string> >(datanames);
  for(size_t i=0;i<tfilenames.size();i++){
    const bool hfile_exists =f_exists(tfilenames[i]);
    if(!hfile_exists){
      Rcpp::Rcerr<<"Missing file: "<<tfilenames[i]<<std::endl;
      Rcpp::stop("file does not exist!");
    }
  }
  Rcpp::NumericMatrix retmat(xmat.nrow(),xmat.ncol());
  Eigen::Map<Eigen::MatrixXd> tretmat(&retmat[0],xmat.nrow(),xmat.ncol());
  Eigen::Map<Eigen::MatrixXd> txmat(&xmat[0],xmat.nrow(),xmat.ncol());
  size_t n_groups=tgroupnames.size();
  std::vector<double> pre_vec(10);

  Progress p(n_groups, true);

  size_t i_offset=0;

  for(size_t i=0;i<n_groups;i++){
    const std::string th5file = tfilenames[i];
    const std::string groupname=tgroupnames[i];
    const std::string dataname=tdatanames[i];
    int row_offset=0;
    int col_offset=0;
    if (Progress::check_abort() )
      Rcpp::stop("Interrupted");
    int row_chunksize=get_rownum_h5(th5file,groupname,dataname);
    int col_chunksize=get_colnum_h5(th5file,groupname,dataname);

    int tot_size = row_chunksize*col_chunksize;
    // Rcpp::Rcerr<<"i: "<<i<<" row_chunksize: "<<row_chunksize<<" col_chunksize: "<<col_chunksize<<" offset: "<<i_offset<<std::endl;
    if(i_offset+row_chunksize>tretmat.rows()){

      Rcpp::stop("This is not going to work");
    }
    pre_vec.resize(tot_size);
    Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> >readmat(pre_vec.data(),row_chunksize,col_chunksize);
    read_2ddmat_h5(th5file,groupname,dataname,row_offset,col_offset,row_chunksize,col_chunksize,readmat.data());
    if(do_transpose){
      tretmat.block(i_offset,0,row_chunksize,tretmat.cols())=readmat.transpose()*txmat.block(i_offset,0,row_chunksize,tretmat.cols());
    }else{
      tretmat.block(i_offset,0,row_chunksize,tretmat.cols())=readmat*txmat.block(i_offset,0,row_chunksize,tretmat.cols());
    }
    i_offset+=row_chunksize;
    p.increment();
  }
  return(retmat);
}




  //[[Rcpp::export]]
Rcpp::NumericMatrix concat_mat_chunks(StringVector h5files, StringVector groupnames,StringVector datanames){



    std::vector<std::string> tfilenames=Rcpp::as<std::vector<std::string> >(h5files);
    std::vector<std::string> tgroupnames=Rcpp::as<std::vector<std::string> >(groupnames);
    std::vector<std::string> tdatanames=Rcpp::as<std::vector<std::string> >(datanames);
    int nfiles =tfilenames.size();
    std::vector<int> rowsize_vec(nfiles);
    std::vector<int> colsize_vec(nfiles);
    for(size_t i=0;i<tfilenames.size();i++){
      const bool hfile_exists =f_exists(tfilenames[i]);
      rowsize_vec[i]=get_rownum_h5(tfilenames[i],tgroupnames[i],tdatanames[i]);
      colsize_vec[i]=get_colnum_h5(tfilenames[i],tgroupnames[i],tdatanames[i]);
      if(!hfile_exists){
        Rcpp::Rcerr<<"Missing file: "<<tfilenames[i]<<std::endl;
        Rcpp::stop("file does not exist!");
      }
    }
    size_t rows= std::accumulate(rowsize_vec.begin(),rowsize_vec.end(),0);
    size_t cols=colsize_vec[0];
    if (std::adjacent_find( colsize_vec.begin(), colsize_vec.end(), std::not_equal_to<int>() ) == colsize_vec.end() )
    {
      Rcpp::stop("All matrices do not have the same number of columns");
    }

    Rcpp::NumericMatrix retmat(rows,cols);
    Eigen::Map<Eigen::MatrixXd> tretmat(&retmat[0],retmat.nrow(),retmat.ncol());
    // Eigen::Map<Eigen::MatrixXd> txmat(&xmat[0],xmat.nrow(),xmat.ncol());
    size_t n_groups=tgroupnames.size();
    std::vector<double> pre_vec(rowsize_vec[0]*cols);

    Progress p(n_groups, true);

    size_t i_offset=0;

    for(size_t i=0;i<n_groups;i++){
      const std::string th5file = tfilenames[i];
      const std::string groupname=tgroupnames[i];
      const std::string dataname=tdatanames[i];
      int row_offset=0;
      int col_offset=0;
      if (Progress::check_abort() )
        Rcpp::stop("Interrupted");
      int row_chunksize=rowsize_vec[i];
      int col_chunksize=cols;

      int tot_size = row_chunksize*col_chunksize;
      if(col_chunksize!=tretmat.cols()){
        Rcpp::stop("There's a matrix with the wrong dimensions!");
      }
      // Rcpp::Rcerr<<"i: "<<i<<" row_chunksize: "<<row_chunksize<<" col_chunksize: "<<col_chunksize<<" offset: "<<i_offset<<std::endl;
      if(i_offset+row_chunksize>tretmat.rows()){

        Rcpp::stop("This is not going to work");
      }
      pre_vec.resize(tot_size);
      Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> >readmat(pre_vec.data(),row_chunksize,col_chunksize);
      read_2ddmat_h5(th5file,groupname,dataname,row_offset,col_offset,row_chunksize,col_chunksize,readmat.data());

      tretmat.block(i_offset,0,row_chunksize,tretmat.cols())=readmat;

      i_offset+=row_chunksize;
      p.increment();
    }
    return(retmat);
  }
