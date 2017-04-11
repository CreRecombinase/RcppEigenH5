#ifndef H5GROUP_H
#define H5GROUP_H

#include "RcppEigenH5_types.h"
#include "h5file.h"




H5GroupPtr create_or_open_group(H5FilePtr &file,const std::string groupname);

H5GroupPtr open_group(H5FilePtr &file,const std::string groupname);

#endif
