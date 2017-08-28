#ifndef H5GROUP_H
#define H5GROUP_H

#include "RcppEigenH5_types.h"
#include "h5file.h"




H5GroupPtr create_or_open_group(const H5FilePtr file,const std::string groupname);

H5GroupPtr open_group(const H5FilePtr file,const std::string groupname);
std::vector<std::string> subgrp_grp(const H5GroupPtr trg);
std::vector<std::string> split_string(const std::string &s, const char delim);
bool grp_path_exists(const H5FilePtr file, const std::string groupname);
#endif
