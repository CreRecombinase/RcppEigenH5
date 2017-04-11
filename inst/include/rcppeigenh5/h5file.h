#ifndef H5FILE_H
#define H5FILE_H

#include <H5Cpp.h>
#include "RcppEigenH5_types.h"
#include <memory>

using namespace H5;

bool f_exists (const std::string name);


H5FilePtr create_or_open_file(const std::string fname);

H5FilePtr open_file(const std::string fname);
#endif
