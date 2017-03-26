#ifndef H5FILE_H
#define H5FILE_H

#include <H5Cpp.h>
#include <memory>

using namespace H5;

 bool f_exists (const std::string name);

typedef std::shared_ptr<H5::H5File> H5FilePtr;


H5FilePtr create_or_open_file(const std::string fname);

H5FilePtr open_file(const std::string fname);
#endif
