#ifndef H5DATASET_H
#define H5DATASET_H
#include "h5dataset.h"
#include "h5group.h"
#include <vector>
#include <H5Cpp.h>
#include "RcppEigenH5_types.h"


using namespace H5;


H5DataSetPtr create_or_open_dataset(H5GroupPtr &group,
                                    const std::string &dataname,
                                    const DataType &data_type,
                                    const std::vector<hsize_t> &cdatadim,const std::vector<hsize_t> &mdatadim,
                                    const std::vector<hsize_t> &chunkdim,const int deflate_level,const bool prealloc=true);
// H5DataSetPtr create_or_open_dataset(H5GroupPtr &group,const std::string &dataname, const DataType &data_type,std::vector<hsize_t> &cdatadim,std::vector<hsize_t> &mdatadim,std::vector<hsize_t> &chunkdim,const int deflate_level);
H5DataSetPtr open_dataset(H5GroupPtr &group,const std::string &dataname);


bool data_exists(const H5GroupPtr grp, const std::string dataname);

#endif
