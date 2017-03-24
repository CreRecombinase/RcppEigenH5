#ifndef H5IO_H
#define H5IO_H
#include <H5Cpp.h>
#include <RcppEigen.h>

using namespace Eigen;

typedef std::ptrdiff_t Index;


typedef Eigen::Ref<Eigen::MatrixXd > Matrix_internal;

template<class ArgType, class RowIndexType, class ColIndexType>
class indexing_functor {
  const ArgType &m_arg;
  const RowIndexType &m_rowIndices;
  const ColIndexType &m_colIndices;
public:
  typedef Matrix<typename ArgType::Scalar,
                 RowIndexType::SizeAtCompileTime,
                 ColIndexType::SizeAtCompileTime,
                 ArgType::Flags&RowMajorBit?RowMajor:ColMajor,
                 RowIndexType::MaxSizeAtCompileTime,
                 ColIndexType::MaxSizeAtCompileTime> MatrixType;
  indexing_functor(const ArgType& arg, const RowIndexType& row_indices, const ColIndexType& col_indices)
    : m_arg(arg), m_rowIndices(row_indices), m_colIndices(col_indices)
  {}
  const typename ArgType::Scalar& operator() (Index row, Index col) const {
    return m_arg(m_rowIndices[row], m_colIndices[col]);
  }
};

template <class ArgType, class RowIndexType, class ColIndexType>
CwiseNullaryOp<indexing_functor<ArgType,RowIndexType,ColIndexType>, typename indexing_functor<ArgType,RowIndexType,ColIndexType>::MatrixType>
indexing(const Eigen::MatrixBase<ArgType>& arg, const RowIndexType& row_indices, const ColIndexType& col_indices)
{
  typedef indexing_functor<ArgType,RowIndexType,ColIndexType> Func;
  typedef typename Func::MatrixType MatrixType;
  return MatrixType::NullaryExpr(row_indices.size(), col_indices.size(), Func(arg.derived(), row_indices, col_indices));
}






void read_2ddmat_h5(const std::string h5file, const std::string groupname, const std::string dataname,const size_t row_offset,const size_t col_offset,const size_t row_chunksize,const size_t col_chunksize, double* data);

Eigen::MatrixXd read_2d_h5(const std::string h5file, const std::string groupname, const std::string dataname,const  Eigen::ArrayXi offset ,const  Eigen::ArrayXi chunksize);
Eigen::MatrixXd read_2d_h5(const std::string h5file, const std::string groupname, const std::string dataname);

void read_dvec_h5(const std::string h5file, const std::string groupname, const std::string dataname,const size_t offset,const size_t chunksize, double* data);
Eigen::ArrayXd read_dvec_h5(const std::string h5file, const std::string groupname, const std::string dataname,const size_t offset,const size_t chunksize);
Eigen::ArrayXd read_dvec_h5(const std::string h5file, const std::string groupname, const std::string dataname);


void read_ivec_h5(const std::string h5file, const std::string groupname, const std::string dataname,const size_t offset,const size_t chunksize, int* data);
Eigen::ArrayXi read_ivec_h5(const std::string h5file, const std::string groupname, const std::string dataname,const size_t offset,const size_t chunksize);
Eigen::ArrayXi read_ivec_h5(const std::string h5file, const std::string groupname, const std::string dataname);


Eigen::MatrixXd read_2d_cindex_h5(const std::string h5file,const std::string groupname, const std::string dataname, const  Eigen::ArrayXi indvec);
Eigen::ArrayXd read_1d_cindex_h5(const std::string h5file,const std::string groupname, const std::string dataname, const  Eigen::ArrayXi indvec);
#endif
