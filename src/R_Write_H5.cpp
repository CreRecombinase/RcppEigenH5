#include "RcppEigenH5.h"
// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
#include <progress.hpp>

//[[Rcpp::depends(RcppEigen)]]
using namespace Eigen;

//[[Rcpp::export(name="write_mat_chunk_h5")]]
void write_mat_chunk_h5_exp(const StringVector h5file, const StringVector groupname, const StringVector dataname,const Matrix_external data,const IntegerVector offsets){

  Eigen::ArrayXi offset=Rcpp::as<Eigen::ArrayXi>(offsets);
  std::string h5fn(h5file[0]);
  std::string gname(groupname[0]);
  std::string dname(dataname[0]);
  write_mat_chunk_h5(h5fn,gname,dname,data,offset);
}




class PReader: public Worker
{
  std::vector<RMatrix<double>> input;
  std::vector<int> snp_offsets;
  std::vector<int> gene_offsets;
  const std::vector<std::vector<std::string>> *in_filenames_l;
  const std::vector<std::vector<size_t>> *chunk_rows; 
  const std::vector<std::string> *out_filenames;
  const std::vector<std::string> *out_datanames;
  const size_t num_snps;
  const size_t num_genes;
  const size_t row_offset;
private:
  std::vector<RMatrix<double>> createInput(std::vector<NumericMatrix> &inp){
    const size_t num_el=inp.size();
    std::vector<RMatrix<double>> retvec(num_el);
    
    for(size_t i=0; i<num_el;i++){
      retvec[i]=RMatrix<double>(inp[i]);
    }
    return(retvec);
  }
  PWriter(const NumericMatrix &input_,
	  const std::vector<std::vector<std::string>> &out_groupnames_l_,
	  const std::vector<std::vector<size_t>> &chunk_rows_; 
	  const std::vector<std::string> &out_filenames_,
	  const std::vector<std::string> &out_datanames_,
	  const size_t row_offset_,
	  const size_t num_snps_,
	  const size_t num_genes_): input(createInput(input_)),
				    out_groupnames_l(out_groupnames_l_),
				    chunk_rows(chunk_rows_),
				    out_filenames(out_filenames_),
				    out_datanames(out_datanames),
				    row_offset(row_offset_),
				    num_snps(num_snps_),
				    num_genes(num_genes_){}
  
  void operator()(std::size_t begin, std::size_t end){
    const size_t num_data=out_datanames.size();
    const FloatType ftypew(PredType::NATIVE_DOUBLE);
    for(size_t i=begin; i<end;i++){
      
      const std::vector<size_t> chunk_r=(*chunk_rows)[i];
      const std::vector<std::string> chunk_groupnames= (*out_groupnames_l)[i];
      const std::string o_filename=(*out_filenames)[i];
      const size_t num_groups=chunk_groupnames.size();
      
      H5FilePtr file_o =create_or_open_file(o_filename);
      for(size_t j=0;j<num_groups;j++){
	const std::string gname_o=chunk_groupnames[j];
	const size_t trownum=chunk_r[j];
	H5GroupPtr group_o = open_group(file_o,i_groupname);
	for(size_t k=0; k<num_data;k++){
	  const std::string o_dataname=o_datanames[k];
	  RMatrix<double>::Col tcol=(*input)[k].col(trownum);
	  const double *datar=&tcol.begin();
	  H5DataSetPtr dataset_o = open_dataset(group_o,o_dataname);
	  
	  DataSpace ofdataspace(dataset_o->getSpace());
	  
          hsize_t datadim[1];
          ofdataspace.getSimpleExtentDims(datadim,NULL);
          const hsize_t memdim[]={num_snps};
          DataSpace mspace(1,memdim);
          if(row_offset+num_snps>datadim[0]){
            Rcpp::stop("Trying to write off the end of the file!");
          }
          //        Rcpp::Rcerr<<"row_offset is :"<<row_offset_i<<std::endl;
          const hsize_t odim[]={row_offset};//dimension of each offset (current_chunk*chunksize)
          //      Rcpp::Rcerr<<"odim is :"<<odim[0]<<"x"<<odim[1]<<std::endl;
	  
          ofdataspace.selectHyperslab( H5S_SELECT_SET,memdim,odim);
          try{
            dataset_o->write(datar,ftypew,mspace,ofdataspace);
          }  catch( DataSetIException error )
	    {
	      error.printError();
	      Rcpp::stop("Error writing file");
	    }
	  //   file_o->flush(H5F_SCOPE_GLOBAL);
          dataset_o->close();
	  ofdataspace->close();
	  mspace->close();
	}
	group_o->close();
      }
      file_o.close();
    }
  }
}



using namespace RcppParallel;

class PWriter: public Worker
{
  const std::vector<RMatrix<double>> input;
  const std::vector<std::vector<std::string>> *out_groupnames_l;
  const std::vector<std::vector<size_t>> *chunk_rows; 
  const std::vector<std::string> *out_filenames;
  const std::vector<std::string> *out_datanames;
  const size_t num_snps;
  const size_t num_genes;
  const size_t row_offset;
private:
  std::vector<RMatrix<double>> createInput(std::vector<NumericMatrix> &inp){
    const size_t num_el=inp.size();
    std::vector<RMatrix<double>> retvec(num_el);
    
    for(size_t i=0; i<num_el;i++){
      retvec[i]=RMatrix<double>(inp[i]);
    }
    return(retvec);
  }
  PWriter(const NumericMatrix &input_,
	  const std::vector<std::vector<std::string>> &out_groupnames_l_,
	  const std::vector<std::vector<size_t>> &chunk_rows_; 
	  const std::vector<std::string> &out_filenames_,
	  const std::vector<std::string> &out_datanames_,
	  const size_t row_offset_,
	  const size_t num_snps_,
	  const size_t num_genes_): input(createInput(input_)),
				    out_groupnames_l(out_groupnames_l_),
				    chunk_rows(chunk_rows_),
				    out_filenames(out_filenames_),
				    out_datanames(out_datanames),
				    row_offset(row_offset_),
				    num_snps(num_snps_),
				    num_genes(num_genes_){}
  
  void operator()(std::size_t begin, std::size_t end){
    const size_t num_data=out_datanames.size();
    const FloatType ftypew(PredType::NATIVE_DOUBLE);
    for(size_t i=begin; i<end;i++){
      
      const std::vector<size_t> chunk_r=(*chunk_rows)[i];
      const std::vector<std::string> chunk_groupnames= (*out_groupnames_l)[i];
      const std::string o_filename=(*out_filenames)[i];
      const size_t num_groups=chunk_groupnames.size();
      
      H5FilePtr file_o =create_or_open_file(o_filename);
      for(size_t j=0;j<num_groups;j++){
	const std::string gname_o=chunk_groupnames[j];
	const size_t trownum=chunk_r[j];
	H5GroupPtr group_o = open_group(file_o,i_groupname);
	for(size_t k=0; k<num_data;k++){
	  const std::string o_dataname=o_datanames[k];
	  RMatrix<double>::Col tcol=(*input)[k].col(trownum);
	  const double *datar=&tcol.begin();
	  H5DataSetPtr dataset_o = open_dataset(group_o,o_dataname);
	  
	  DataSpace ofdataspace(dataset_o->getSpace());
	  
          hsize_t datadim[1];
          ofdataspace.getSimpleExtentDims(datadim,NULL);
          const hsize_t memdim[]={num_snps};
          DataSpace mspace(1,memdim);
          if(row_offset+num_snps>datadim[0]){
            Rcpp::stop("Trying to write off the end of the file!");
          }
          //        Rcpp::Rcerr<<"row_offset is :"<<row_offset_i<<std::endl;
          const hsize_t odim[]={row_offset};//dimension of each offset (current_chunk*chunksize)
          //      Rcpp::Rcerr<<"odim is :"<<odim[0]<<"x"<<odim[1]<<std::endl;
	  
          ofdataspace.selectHyperslab( H5S_SELECT_SET,memdim,odim);
          try{
            dataset_o->write(datar,ftypew,mspace,ofdataspace);
          }  catch( DataSetIException error )
	    {
	      error.printError();
	      Rcpp::stop("Error writing file");
	    }
	  //   file_o->flush(H5F_SCOPE_GLOBAL);
          dataset_o->close();
	  ofdataspace->close();
	  mspace->close();
	}
	group_o->close();
      }
      file_o.close();
    }
  }
}
	  
  

//[[Rcpp::export]]
void create_groups_rows_split_cols_h5(const StringVector in_h5files,
                                      const StringVector in_groupname,
                                      const StringVector in_datanames,
                                      const List out_groupname_list){


  const std::vector<std::string> i_datanames=Rcpp::as<std::vector<std::string>>(in_datanames);
  const size_t num_data=i_datanames.size();
  const std::vector<std::string> out_filenames=Rcpp::as<std::vector<std::string> >(out_groupname_list.attr("names"));
  const size_t num_out_files=out_filenames.size();
  std::vector<std::vector<std::string> > out_groupnames_l(num_out_files);
  size_t tot_out_cols=0;
  for(size_t i=0; i<num_out_files;i++){
    out_groupnames_l[i]=Rcpp::as<std::vector<std::string> >(out_groupname_list[i]);
    tot_out_cols+=out_groupnames_l[i].size();
  }

  const size_t n_in_files=in_h5files.size();
  size_t tot_rows=0;

  const std::string i_groupname(in_groupname[0]);


  std::vector<size_t> tot_cols_v(n_in_files);
  std::vector<size_t> tot_rows_v(n_in_files);

  size_t max_rows=0;
  for(size_t i=0;i<n_in_files;i++){
    size_t temp_rows=get_rownum_h5(std::string(in_h5files[i]),i_groupname,i_datanames[0]);
    tot_cols_v[i]= get_colnum_h5(std::string(in_h5files[i]),i_groupname,i_datanames[0]);
    for(size_t j=0; j<num_data;j++){
      if((get_rownum_h5(std::string(in_h5files[i]),i_groupname,i_datanames[j])!=temp_rows)||
         (get_colnum_h5(std::string(in_h5files[i]),i_groupname,i_datanames[j])!=tot_cols_v[i])){
        Rcpp::stop("all data must have the same dimensions");
      }
    }

    tot_rows+=temp_rows;
    tot_rows_v[i]=temp_rows;
    if(temp_rows>max_rows){
      max_rows=temp_rows;
    }

    if(i>0){
      if(tot_cols_v[i-1]!=tot_cols_v[i]){
        Rcpp::stop("All files must have equal number of columns");
      }
    }
  }
  std::vector<size_t> row_offset(n_in_files);
  row_offset[0]=0;
  std::partial_sum(tot_rows_v.begin(),tot_rows_v.end()-1,row_offset.begin()+1);


  const size_t tot_cols=tot_cols_v[0];
  if(tot_cols!=tot_out_cols){
    Rcpp::stop("Groupnames not the same size as number of columns!");
  }

  const FloatType ftypew(PredType::NATIVE_DOUBLE);

  const std::vector<hsize_t> cumdim{tot_rows};
  const std::vector<hsize_t> maxdim{tot_rows};
  const std::vector<hsize_t> chunkdim{tot_rows/2};


  Progress pp(tot_cols, true);
  for(size_t i=0; i<num_out_files;i++){
    const std::vector<std::string> *out_groupnames=&out_groupnames_l[i];
    const std::string o_filename=out_filenames[i];
    H5FilePtr file_o =create_or_open_file(o_filename);
    const size_t num_cols_chunk=out_groupnames->size();
    for(size_t i=0; i<num_cols_chunk;i++){
      if(Progress::check_abort()){
        Rcpp::stop("Process Aborted");
      }
      const std::string gname((*out_groupnames)[i]);
      H5GroupPtr group =create_or_open_group(file_o,gname);
      for(size_t j=0; j<num_data;j++){
        const std::string o_dataname=i_datanames[j];
        H5DataSetPtr dataset =create_or_open_dataset(group,o_dataname,ftypew,cumdim,maxdim,chunkdim,1);
        write_transpose(dataset,false);
        dataset->close();
      }
      group->close();
      pp.increment();
    }
    file_o->flush(H5F_SCOPE_GLOBAL);
    file_o->close();
  }
}



//[[Rcpp::export]]
void concat_rows_split_cols_h5(const StringVector in_h5files,
                               const StringVector in_groupname,
                               const StringVector in_datanames,
                               const List out_groupname_list){

  const std::vector<std::string> i_datanames=Rcpp::as<std::vector<std::string>>(in_datanames);
  const size_t num_data=i_datanames.size();
  const std::vector<std::string> out_filenames=Rcpp::as<std::vector<std::string> >(out_groupname_list.attr("names"));
  const size_t num_out_files=out_filenames.size();
  std::vector<std::vector<std::string> > out_groupnames_l(num_out_files);
  size_t tot_out_cols=0;
  for(size_t i=0; i<num_out_files;i++){
    out_groupnames_l[i]=Rcpp::as<std::vector<std::string> >(out_groupname_list[i]);
    tot_out_cols+=out_groupnames_l[i].size();
  }

  const size_t n_in_files=in_h5files.size();
  size_t tot_rows=0;

  const std::string i_groupname(in_groupname[0]);


  std::vector<size_t> tot_cols_v(n_in_files);
  std::vector<size_t> tot_rows_v(n_in_files);

  size_t max_rows=0;
  for(size_t i=0;i<n_in_files;i++){
    size_t temp_rows=get_rownum_h5(std::string(in_h5files[i]),i_groupname,i_datanames[0]);
    tot_cols_v[i]= get_colnum_h5(std::string(in_h5files[i]),i_groupname,i_datanames[0]);
    for(size_t j=0; j<num_data;j++){
      if((get_rownum_h5(std::string(in_h5files[i]),i_groupname,i_datanames[j])!=temp_rows)||
         (get_colnum_h5(std::string(in_h5files[i]),i_groupname,i_datanames[j])!=tot_cols_v[i])){
        Rcpp::stop("all data must have the same dimensions");
      }
    }

    tot_rows+=temp_rows;
    tot_rows_v[i]=temp_rows;
    if(temp_rows>max_rows){
      max_rows=temp_rows;
    }

    if(i>0){
      if(tot_cols_v[i-1]!=tot_cols_v[i]){
        Rcpp::stop("All files must have equal number of columns");
      }
    }
  }
  std::vector<size_t> row_offset(n_in_files);
  row_offset[0]=0;
  std::partial_sum(tot_rows_v.begin(),tot_rows_v.end()-1,row_offset.begin()+1);


  const size_t tot_cols=tot_cols_v[0];
  if(tot_cols!=tot_out_cols){
    Rcpp::stop("Groupnames not the same size as number of columns!");
  }

  const FloatType ftypew(PredType::NATIVE_DOUBLE);
  std::vector<Rcpp::NumericMatrix> mat_vec(num_data);
  for(size_t i=0; i<num_data;i++){
    mat_vec[i]=Rcpp::NumericMatrix(tot_cols,max_rows);
  }
  


  Progress pp(tot_cols*tot_rows, true);
  //  file_o->close();

  // file_o =create_or_open_file(std::string(out_h5file[0]));
  DataSpace* ofdataspace;
  for(size_t i=0; i<n_in_files;i++){

//The  layout (abstractly) is snps as rows, genes as columns.  The data is "transposed" though (also stored with a row-major storage order),
//so that means the number of genes is the first dimension (on disk)
    const size_t num_snps=tot_rows_v[i];
    const size_t num_genes=tot_cols;
    //    dbuffer.resize(num_snps*num_genes);
    H5FilePtr file_i =open_file(std::string(in_h5files[i]));
    H5GroupPtr group_i = open_group(file_i,i_groupname);
    const hsize_t matrix_dims[]={num_genes,num_snps};
    const hsize_t offseta[]={0,0};
    for(size_t k=0;k<num_data;k++){
      H5DataSetPtr dataset_i = open_dataset(group_i,i_datanames[k]);
      bool isTranspose=check_transpose(dataset_i);
      if(!isTranspose){
        Rcpp::stop("Not transposed data interface hasn't been implemented yet(sorry)");
      }
      double *dbuffer=(mat_vec[k]).begin();

      DataType dt= dataset_i->getDataType();
      hsize_t datadims[]={0,0};
      DataSpace fspace=dataset_i->getSpace();
      fspace.getSimpleExtentDims(datadims,NULL);
      //  std::cout<<"Full data is of dimensions"<<datadims[0]<<"x"<<datadims[1]<<std::endl;

      if(num_genes>datadims[0]){
        Rcpp::stop("You can't read beyond extent in first dimension (rows)");
      }
      if(num_snps>datadims[1]){
        Rcpp::stop("You can't read beyond extent in the second dimension (cols)");
      }
      //Dimensions in memory (storing data in row major order for now (because why not)


      fspace.selectHyperslab(H5S_SELECT_SET,datadims,offseta);

      DataSpace memspace(2,matrix_dims);
      dataset_i->read(dbuffer.data(),dt,memspace,fspace);

      Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> > datmap(dbuffer.data(),num_genes,num_snps);


      size_t gene_iter=0;
      for(size_t l=0;l<num_out_files;l++){
        const std::vector<std::string> *out_groupnames=&out_groupnames_l[l];
        const std::string o_filename=out_filenames[l];
        H5FilePtr file_o =create_or_open_file(o_filename);
        const size_t num_cols_chunk=out_groupnames->size();

        for(size_t j=0; j<num_cols_chunk;j++){
          //Remember that data is transposed
          const std::string gname((*out_groupnames)[j]);
          const double *trow=&datmap.coeffRef(gene_iter,0);
          gene_iter++;
          if(Progress::check_abort()){
            Rcpp::stop("Process Aborted");
          }

          //        Rcpp::Rcout<<std::endl<<"Here's a peek at tcol(length is "<<num_snps<<"):"<<std::endl<<Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>>(trow,1,num_snps)<<std::endl;
          H5GroupPtr group_o = open_group(file_o,gname);

          H5DataSetPtr dataset_o = open_dataset(group_o,i_datanames[k]);

          try{
            ofdataspace= new DataSpace(dataset_o->getSpace());
          }catch(DataSpaceIException error){
            error.printError();
            Rcpp::stop("Error creating memory dataspace ");
          }
          hsize_t datadim[1];
          ofdataspace->getSimpleExtentDims(datadim,NULL);
          const hsize_t memdim[]={num_snps};
          DataSpace *mspace;
          try{
            mspace= new DataSpace(1,memdim); //Size of first dataset (in memory, can be bigger or smaller than size on disk, depending on how much you're writing)
          }catch(DataSpaceIException error){
            error.printError();
            Rcpp::stop("Error creating memory dataspace ");
          }
          hsize_t row_offset_i=row_offset[i];

          if(row_offset_i+num_snps>datadim[0]){
            Rcpp::stop("Trying to write off the end of the file!");
          }
          //        Rcpp::Rcerr<<"row_offset is :"<<row_offset_i<<std::endl;
          const hsize_t odim[]={row_offset_i};//dimension of each offset (current_chunk*chunksize)
          //      Rcpp::Rcerr<<"odim is :"<<odim[0]<<"x"<<odim[1]<<std::endl;

          ofdataspace->selectHyperslab( H5S_SELECT_SET,memdim,odim);
          try{
            dataset_o->write(trow,ftypew,*mspace,*ofdataspace);
          }  catch( DataSetIException error )
          {
            error.printError();
            Rcpp::stop("Error writing file");
          }
          file_o->flush(H5F_SCOPE_GLOBAL);

          //  std::cout<<"File flushed"<<std::endl;
          dataset_o->close();
          ofdataspace->close();
          mspace->close();
          delete mspace;
          group_o->close();
          delete ofdataspace;
          pp.increment();
        }
        file_o->flush(H5F_SCOPE_GLOBAL);
        file_o->close();
      }

      dataset_i->close();
      dt.close();
      memspace.close();
      fspace.close();
    }
    group_i->close();
    file_i->close();
  }
  // file_o->flush(H5F_SCOPE_GLOBAL);
  // file_o->close();




}




//[[Rcpp::export(name="create_mat_dataset_h5")]]
void create_mat_dataset_h5_exp(const StringVector h5file, const StringVector groupname, const StringVector dataname,const IntegerVector dims, const IntegerVector chunkdims, const IntegerVector deflate_level=IntegerVector::create(0),const bool doTranspose=false){

  std::string filename(h5file[0]);

  size_t groupnum=groupname.size();
  std::vector<std::string> groupnames(groupnum);
  for(size_t i=0;i<groupnum;i++){
    groupnames[i]=groupname[i];
  }




//  std::string gname(groupname[0]);
  std::string dname(dataname[0]);


  int d_level =deflate_level[0];





  size_t rownum =dims[0];
  size_t colnum=dims[1];

  size_t rowchunk=chunkdims[0];
  size_t colchunk=chunkdims[1];

  if(doTranspose){
    std::swap(rownum,colnum);
    std::swap(rowchunk,colchunk);

  }
  H5FilePtr file =create_or_open_file(filename);
  FloatType ftypew(PredType::NATIVE_DOUBLE);

  std::vector<hsize_t> cumdim{rownum,colnum};
  std::vector<hsize_t> maxdim{rownum,colnum};
  std::vector<hsize_t> chunkdim{rowchunk,colchunk};
  for(size_t i=0; i<groupnum;i++){
    std::string gname=groupnames[i];
    H5GroupPtr group =create_or_open_group(file,gname);

    H5DataSetPtr dataset =create_or_open_dataset(group,dname,ftypew,cumdim,maxdim,chunkdim,d_level);
    write_transpose(dataset,doTranspose);
    dataset->close();
    group->close();
  }
  file->close();

}


// //[[Rcpp::export]]
// Eigen::MatrixXd test_rowcol(const Matrix_external data){
//   Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>
//  }




//[[Rcpp::export(name="write_mat_h5")]]
void write_mat_h5_exp(const StringVector h5file, const StringVector groupname, const StringVector dataname,const Matrix_external data, const IntegerVector deflate_level=IntegerVector::create(0),const bool doTranspose=false){

  std::string th5file(h5file[0]);
  std::string tgroupname(groupname[0]);
  std::string tdataname(dataname[0]);
  int deflate=deflate_level[0];

  Eigen::ArrayXi chunksize(2);
  chunksize[0]=data.rows();
  chunksize[1]=std::min(10000,(int)data.cols());
//
//   if(bdoTranspose){
//     std::swap(chunksize[0],chunksize[1]);
//   }

  write_mat_h5(th5file,tgroupname,tdataname,data,chunksize,deflate,doTranspose);
}


//[[Rcpp::export(name="write_dvec_h5")]]
void write_dvec_h5_exp(const StringVector h5file, const StringVector groupname, const StringVector dataname,const arrayxd_external data, const IntegerVector deflate_level=IntegerVector::create(0)){

  std::string th5file(h5file[0]);
  std::string tgroupname(groupname[0]);
  std::string tdataname(dataname[0]);
  int deflate=deflate_level[0];
  write_dvec_h5(th5file,tgroupname,tdataname,data,deflate);
}


//[[Rcpp::export(name="write_ivec_h5")]]
void write_ivec_h5_exp(const StringVector h5file, const StringVector groupname, const StringVector dataname,const arrayxi_external data, const IntegerVector deflate_level=IntegerVector::create(0)){

  std::string th5file(h5file[0]);
  std::string tgroupname(groupname[0]);
  std::string tdataname(dataname[0]);
  int deflate=deflate_level[0];
  write_ivec_h5(th5file,tgroupname,tdataname,data,deflate);
}



