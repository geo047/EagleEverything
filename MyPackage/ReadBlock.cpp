// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
// #include  "CHDF5File.h"
#ifdef _OPENMP
#include <omp.h>
//   [[Rcpp::plugins(openmp)]]
#endif






// [[Rcpp::export]]
Eigen::MatrixXd  ReadBlock(std::string asciifname,
                           long start_row,
                           long numcols,
                           long numrows_in_block)

{
 // reads in data from ASCII file 
 // to form M Eign double matrix 
std::ostringstream
      os;
std::string
   line;


long
  // coli = 0, 
  rowi = 0;


Eigen::MatrixXd
      M(numrows_in_block, numcols) ;


// Open no-space ASCII file
   std::ifstream fileIN(asciifname.c_str(), std::ios::in );

    if(!fileIN.good()) {
      os << "ERROR: Could not open  " << asciifname << std::endl;
      Rcpp::stop(os.str() );
     }

   for(long rr=0; rr < (start_row + numrows_in_block) ; rr++){
      // read a line of data from ASCII file
      getline(fileIN, line);
      if(rr >= start_row){
          std::istringstream streamA(line);
          for(long ii=0; ii < numcols  ; ii++){
            int tmp  = line[ii] - '0'; // trick to removes ASCII character offset for numbers
            M(rowi, ii) = (double) tmp - 1;   // converting data to -1, 0, 1 
          }
          rowi++;
      } // end if rr
   } // end for(rr



// Close the ascii file
   fileIN.close();


 return M;

}


/*
//'
//'
// [[Rcpp::export]]
Eigen::MatrixXd  ReadBlockWriteHDF(std::string asciifname,
                           long numcols,
                           long numrows)

{
    // reads in data from ASCII file 
    // and writes to a hdf5 file
    std::ostringstream os;
    std::string line;

    // replace the ascii filetype
    std::size_t lastindex = asciifname.find_last_of("."); 
    std::string rawname = asciifname.substr(0, lastindex); 
    std::string hdffilename(rawname+std::string(".h5")) ;

    CHDF5File<char>* tempHDF5File = (CHDF5File<char>*) new CHDF5File<char>(hdffilename , H5F_ACC_RDWR ) ;  // H5F_ACC_RDONLY
    hid_t datatype_mem = H5T_NATIVE_SCHAR;  //  H5T_NATIVE_DOUBLE ;
    hid_t datatype_disk = H5T_NATIVE_SCHAR;
    tempHDF5File->CreateDatasetNoCompression(std::string("/eagle_array"), NULL ,numrows, numrows*numcols, datatype_mem, datatype_disk) ;

    signed char * char_line = new signed char[numcols] ;
    // Open no-space ASCII file
    std::ifstream fileIN(asciifname.c_str(), std::ios::in );

    if(!fileIN.good()) {
      os << "ERROR: Could not open  " << asciifname << std::endl;
      Rcpp::stop(os.str() );
     }
     
   std::size_t rowi = 0 ;
   for(long rr=0; rr < (numrows) ; rr++){
      // read a line of data from ASCII file
      getline(fileIN, line);
      std::istringstream streamA(line);
      for(long ii=0; ii < numcols  ; ii++){
            int tmp  = line[ii] - '0'; // trick to removes ASCII character offset for numbers
            char_line[ii] = (signed char) tmp - 1;   // converting data to -1, 0, 1 
     }
     
     tempHDF5File->WriteRowsToExistingDataset(std::string("/eagle_array"), (void *) char_line, 1 , numcols, rowi, datatype_mem ) ;
     rowi++;
   } // end for(rr

    // Close the ascii file
    fileIN.close();
    // Close the hdf5 file
    delete tempHDF5File ;
    delete [] char_line ;

 return M;

}



// [[Rcpp::export]]
Eigen::MatrixXd  ReadHDF(std::string hdffname,
                           long start_row,
                           long numcols,
                           long numrows_in_block)

{
 // reads in data from HDF5 file 
 // to form M Eign double matrix 

    Eigen::MatrixXd  M(numrows_in_block, numcols) ;


    CHDF5File<char>* tempHDF5File = (CHDF5File<char>*) new CHDF5File<char>(hdffilename , H5F_ACC_RDWR ) ;  // H5F_ACC_RDONLY
    hid_t datatype_mem = H5T_NATIVE_DOUBLE;  //  H5T_NATIVE_DOUBLE ;
    hid_t datatype_disk = H5T_NATIVE_SCHAR;
    tempHDF5File->GetSectionOfDatasetReturnInPtr( std::string("/eagle_array"),   T *&  memStreamIn, start_row, numrows_in_block,  datatype_mem) ;
    // Close the hdf5 file
    delete tempHDF5File ;
    
    
 return M;

}

*/


