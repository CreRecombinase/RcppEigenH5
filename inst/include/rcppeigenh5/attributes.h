#ifndef ATTRIBUTES_H
#define ATTRIBUTES_H

#include <H5Cpp.h>

std::string read_data_attr_h5(const std::string h5file, const std::string groupname, const std::string dataname,const std::string attrname);
std::string read_group_attr_h5(const std::string h5file, const std::string groupname, const std::string attrname);

#endif
