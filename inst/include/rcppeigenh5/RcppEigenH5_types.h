#ifndef RCPPEIGENH5_TYPES_H
#define RCPPEIGENH5_TYPES_H

#include <RcppEigen.h>
#include <H5Cpp.h>
#include <memory>

using namespace H5;

typedef std::shared_ptr<Group> H5GroupPtr;
typedef std::shared_ptr<H5File> H5FilePtr;
typedef std::shared_ptr<DataSet> H5DataSetPtr;
typedef std::shared_ptr<DataSpace> H5DataSpacePtr;

typedef Eigen::Array<double, Eigen::Dynamic, 1> ColumnArray;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> ColumnMatrix;
typedef Eigen::Matrix<char, Eigen::Dynamic, 1> ColumnMatrixi;
typedef Eigen::Array<double, 1,Eigen::Dynamic > RowArray;


typedef Eigen::Ref<const ColumnArray> c_column_internal;
typedef Eigen::Ref<const RowArray> c_row_internal;
typedef Eigen::Ref<const Eigen::ArrayXd> c_arrayxd_internal;
typedef Eigen::Ref<const Eigen::ArrayXi> c_arrayxi_internal;
typedef Eigen::Ref<const Eigen::ArrayXXd> c_arrayxxd_internal;
typedef Eigen::Ref<const Eigen::VectorXd> c_vectorxd_internal;
typedef Eigen::MappedSparseMatrix<double> c_sparseMatrix_internal;

typedef Eigen::Ref<const Eigen::MatrixXd > c_Matrix_internal;
typedef Eigen::Ref<Eigen::MatrixXd > Matrix_internal;
typedef Eigen::Ref<Eigen::MatrixXf > Matrix_finternal;

typedef Eigen::Ref<Eigen::ArrayXd> arrayxd_internal;
typedef Eigen::Ref<Eigen::ArrayXi> arrayxi_internal;

typedef Eigen::Ref<Eigen::ArrayXXd> arrayxxd_internal;
typedef Eigen::Ref<ColumnArray> column_internal;
typedef Eigen::Ref<RowArray> row_internal;

typedef Eigen::Map<Eigen::ArrayXd> arrayxd_external;
typedef Eigen::Map<Eigen::ArrayXi> arrayxi_external;
typedef Eigen::MappedSparseMatrix<double> sparseMatrix_external;
typedef Eigen::Map<Eigen::MatrixXd > Matrix_external;
typedef Eigen::Map<Eigen::VectorXd> vectorxd_external;

typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> rowmat;
typedef Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> rowfmat;
typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> colmat;

typedef Eigen::Map<rowmat> mrowmat;
typedef Eigen::Map<colmat> mmat;

typedef std::tuple<std::string, std::string,std::string> hnames;
typedef std::pair<size_t,size_t> matdim;


namespace Eigen
{

#define EIGEN_MAKE_TYPEDEFS(Size, SizeSuffix)                     \
/** \ingroup matrixtypedefs */                                    \
template <typename Type>                                          \
using Matrix##SizeSuffix = Matrix<Type, Size, Size>;              \
/** \ingroup matrixtypedefs */                                    \
template <typename Type>                                          \
using Vector##SizeSuffix = Matrix<Type, Size, 1>;                 \
/** \ingroup matrixtypedefs */                                    \
template <typename Type>                                          \
using RowVector##SizeSuffix = Matrix<Type, 1, Size>;

#define EIGEN_MAKE_FIXED_TYPEDEFS(Size)                           \
/** \ingroup matrixtypedefs */                                    \
template <typename Type>                                          \
using Matrix##Size##X = Matrix<Type, Size, Dynamic>;              \
/** \ingroup matrixtypedefs */                                    \
template <typename Type>                                          \
using Matrix##X##Size = Matrix<Type, Dynamic, Size>;

EIGEN_MAKE_TYPEDEFS(2, 2)
  EIGEN_MAKE_TYPEDEFS(3, 3)
  EIGEN_MAKE_TYPEDEFS(4, 4)
  EIGEN_MAKE_TYPEDEFS(Dynamic, X)
  EIGEN_MAKE_FIXED_TYPEDEFS(2)
  EIGEN_MAKE_FIXED_TYPEDEFS(3)
  EIGEN_MAKE_FIXED_TYPEDEFS(4)

#undef EIGEN_MAKE_TYPEDEFS
#undef EIGEN_MAKE_FIXED_TYPEDEFS

}


#endif
