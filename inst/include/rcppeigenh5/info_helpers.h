#ifndef INFO_HELPERS_H
#define INFO_HELPERS_H

#include"RcppEigenH5_types.h"

std::string print_dims(c_Matrix_internal mat);
int get_rownum_h5(const std::string h5file,const std::string groupname, const std::string dataname);

int get_colnum_h5(const std::string h5file,const std::string groupname, const std::string dataname);

std::vector<std::string> list_subgroups(const H5FilePtr file,const std::string base);
std::vector<std::string>h5_ls_dset(const H5GroupPtr grp,const std::string basename);
#endif
