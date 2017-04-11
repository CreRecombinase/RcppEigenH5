#ifndef INFO_HELPERS_H
#define INFO_HELPERS_H

#include"RcppEigenH5_types.h"

std::string print_dims(c_Matrix_internal mat);
int get_rownum_h5(const std::string h5file,const std::string groupname, const std::string dataname);

int get_colnum_h5(const std::string h5file,const std::string groupname, const std::string dataname);

#endif
