#ifndef WRITE_DATA_H
#define WRITE_DATA_H
#include "RcppEigenH5_types.h"

void write_mat_chunk_h5(const std::string h5file, const std::string groupname, const std::string dataname, Matrix_internal data, c_arrayxi_internal offsets);

void write_mat_h5(const std::string h5file, const std::string groupname, const std::string dataname, Matrix_internal data, const int deflate_level=0);
void write_mat_h5(const std::string h5file, const std::string groupname, const std::string dataname, Matrix_internal data, Eigen::ArrayXi chunksize, const int deflate_level,bool doTranspose);

void write_dvec_h5(const std::string h5file, const std::string groupname, const std::string dataname, c_arrayxd_internal data, const int deflate_level);

void write_ivec_h5(const std::string h5file, const std::string groupname, const std::string dataname, c_arrayxi_internal data, const int deflate_level);
#endif
