#ifndef ATTRIBUTES_H
#define ATTRIBUTES_H

#include <H5Cpp.h>
#include "RcppEigen.h"
using namespace Rcpp;

std::string read_data_attr_h5(const StringVector h5filename, const StringVector h5_groupname, const StringVector h5_dataname,const StringVector h5_attr_name);
std::string read_data_attr_h5(const std::string h5file, const std::string groupname, const std::string dataname,const std::string attrname);
std::string read_group_attr_h5(const std::string h5file, const std::string groupname, const std::string attrname);

#endif
