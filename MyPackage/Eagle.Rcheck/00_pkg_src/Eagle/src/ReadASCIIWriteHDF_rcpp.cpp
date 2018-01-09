// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include  "CHDF5File.h"
#ifdef _OPENMP
#include <omp.h>
//   [[Rcpp::plugins(openmp)]]
#endif


// This function is called from ReadMarker() twise, once for M and again for Mt
// [[Rcpp::export]]
bool  ReadASCIIWriteHDF_rcpp(std::string asciifname,
                  std::vector <long> dims,
                  bool transpose                 )

{
    int numrows, numcols ;
    if (transpose == false) {
         numcols = dims[1] ;
         numrows = dims[0] ;
    } else {
         numcols = dims[0] ;
         numrows = dims[1] ;  
    }
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
     
   Rcpp::Rcout << " ReadASCIIWriteHDF_rcpp : numrows=" << numrows  << std::endl ;
   Rcpp::Rcout << " ReadASCIIWriteHDF_rcpp : numcols="  << numcols << std::endl ;
   
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

    return true;

}



