#ifndef ATTRIBUTES_H
#define ATTRIBUTES_H

#include <H5Cpp.h>
#include "RcppEigen.h"
#include "RcppEigenH5_types.h"
#include "h5file.h"
#include "h5group.h"



using namespace Rcpp;
int write_transpose(H5DataSetPtr dataset, bool doTranspose);
bool check_transpose(H5DataSetPtr dataset);
int get_int_attr(H5DataSetPtr dataset,std::string attr_name);
std::string read_data_attr_h5(const StringVector h5filename, const StringVector h5_groupname, const StringVector h5_dataname,const StringVector h5_attr_name);
std::string read_data_attr_h5(const std::string h5file, const std::string groupname, const std::string dataname,const std::string attrname);
std::string read_group_attr_h5(const std::string h5file, const std::string groupname, const std::string attrname);


void write_group_string_attr_h5(const std::string h5file, const std::string groupname,const std::string attr_name,const std::string attr_value);
void write_data_string_attr_h5(const std::string h5file, const std::string groupname, const std::string dataname,const std::string attr_name,const std::string attr_value);
void write_group_int_attr_h5(const std::string h5file, const std::string groupname,const std::string attr_name,const int attr_value);
void write_data_int_attr_h5(const std::string h5file, const std::string groupname, const std::string dataname,const std::string attr_name,const int attr_value);

int read_igroup_attr_h5(const std::string h5file, const std::string groupname,const std::string attr_name);
int read_idata_attr_h5(const std::string h5file, const std::string groupname, const std::string dataname,const std::string attr_name);

std::vector<int> read_data_iarray_attr_h5(const std::string h5file, const std::string groupname, const std::string dataname,const std::string attrname);
std::vector<int> read_group_iarray_attr_h5(const std::string h5file, const std::string groupname, const std::string attrname);

#endif
