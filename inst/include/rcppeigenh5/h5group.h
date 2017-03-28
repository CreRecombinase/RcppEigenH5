#ifndef H5GROUP_H
#define H5GROUP_H
#include <H5Cpp.h>
#include "h5file.h"
#include <memory>

typedef std::shared_ptr<H5::Group> H5GroupPtr;

H5GroupPtr create_or_open_group(H5FilePtr &file,const std::string groupname);

H5GroupPtr open_group(H5FilePtr &file,const std::string groupname);

#endif
