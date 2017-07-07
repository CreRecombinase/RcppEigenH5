#include "RcppEigenH5.h"
//[[Rcpp::depends(BH,RcppEigen)]]
#include <boost/algorithm/string.hpp>
//#include <boost/iostreams/filtering_stream.hpp>
//#include <boost/iostreams/filter/gzip.hpp>
//#include <boost/iostreams/copy.hpp>
#include <boost/utility/string_ref.hpp>
#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/qi_match.hpp>
#include <boost/spirit/include/qi_uint.hpp>

//[[Rcpp::export]]
int file_rownum(const std::string filename){
  std::ifstream myfile(filename);

  // new lines will be skipped unless we stop it from happening:
  myfile.unsetf(std::ios_base::skipws);

  // count the newlines with an algorithm specialized for counting:
  int line_count = std::count(
    std::istream_iterator<char>(myfile),
    std::istream_iterator<char>(),
    '\n');
  myfile.close();
  return(line_count);
}

//[[Rcpp::export]]
int file_colnum(const std::string filename){


  namespace qi = boost::spirit::qi;
  namespace ascii = boost::spirit::ascii;
  namespace phoenix = boost::phoenix;

  using qi::int_;
  using qi::_1;
  using ascii::space;
  using phoenix::ref;
  using namespace std;
  int n=0;

  ifstream myfile(filename);
  string iline;
  getline(myfile,iline);
  myfile.close();
//  Rcpp::Rcout<<iline<<std::endl;
bool r= qi::parse(iline.begin(),iline.end(),
                         (
int_[ref(n)+=1]>>*(' ' >> int_[ref(n)+=1])
                         ));

  return(n);
}

template <typename Iterator> void parse_haps_int(Iterator first, Iterator last,std::vector<int> &datavec){

  namespace qi = boost::spirit::qi;
  namespace ascii = boost::spirit::ascii;
  namespace phoenix = boost::phoenix;

  using qi::int_;
  using qi::_1;
  using ascii::space;
  using phoenix::ref;
  using namespace std;
  using phoenix::push_back;
  qi::phrase_parse(first,last,

                   //  Begin grammar
                   (
                       *int_[push_back(phoenix::ref(datavec), _1)]
                   )
                     ,
                     //  End grammar
                     qi::space);

}

template <typename Iterator> void parse_haps_double(Iterator first, Iterator last,std::vector<double> &datavec){

  namespace qi = boost::spirit::qi;
  namespace ascii = boost::spirit::ascii;
  namespace phoenix = boost::phoenix;

  using qi::double_;
  using qi::_1;
  using ascii::space;
  using phoenix::ref;
  using namespace std;
  using phoenix::push_back;
  qi::phrase_parse(first,last,

                   //  Begin grammar
                   (
                       *double_[push_back(phoenix::ref(datavec), _1)]
                   )
                     ,
                     //  End grammar
                     qi::space);

}


//[[Rcpp::export(name="read_haps")]]
Rcpp::IntegerMatrix read_haps_exp(const std::string filename, int rownum=0, int colnum=0){
  using namespace Rcpp;

  using namespace Eigen;

  if(rownum==0){
    rownum=file_rownum(filename);
  }
  if(colnum==0){
    colnum=file_colnum(filename);
  }

  boost::iostreams::mapped_file_source file(filename);
  size_t matsize=rownum*colnum;
  std::vector<int> datavec;
  datavec.reserve(matsize);
  parse_haps_int(file.data(),file.data()+file.size(),datavec);


  //  Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> > readmat((double*)x,rownum,t_chunksize);

  Eigen::Map<Eigen::Matrix<int,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> > M(datavec.data(),rownum,colnum);
  Rcpp::IntegerMatrix retmat(rownum,colnum);
  Eigen::Map<Eigen::Matrix<int,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> > rM(&retmat[0],rownum,colnum);
  rM=M;
  return(retmat);
}


//[[Rcpp::export]]
void write_haps_h5(const std::string filename,const std::string h5file, const std::string groupname, const std::string dataname, const int deflate_level, int rownum=0,int colnum=0){
  using namespace Rcpp;

  using namespace std;
  using namespace Eigen;
  if(rownum==0){
    rownum=file_rownum(filename);
  }
  if(colnum==0){
    colnum=file_colnum(filename);
  }

  boost::iostreams::mapped_file_source file(filename);
  size_t matsize=rownum*colnum;
  std::vector<double> datavec;
  datavec.reserve(matsize);
  parse_haps_double(file.data(),file.data()+file.size(),datavec);
  hnames names(h5file,groupname,dataname);
  matdim dims(rownum,colnum);
  write_dmat_h5(names,datavec.data(),deflate_level,dims);
  //  Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> > readmat((double*)x,rownum,t_chunksize);

}

