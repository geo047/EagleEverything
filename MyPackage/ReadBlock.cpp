// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#if defined (_HASHDF5)
#include  "CHDF5File.h"
#endif  // defined (_HASHDF5)
#ifdef _OPENMP
#include <omp.h>
//   [[Rcpp::plugins(openmp)]]
#endif

void PrintEigenRowCols(std::string mat_filename, Eigen::MatrixXd Min) ;


#if defined (_HASHDF5)
// [[Rcpp::export]]
Eigen::MatrixXd  ReadBlock(std::string asciifname,
                           long start_row,
                           long numcols,
                           long numrows_in_block)

{
    // Use this to get the transposed version so that Eigen matrix class is still Column Major
    std::size_t dir_index = asciifname.find_last_of(SYSTEMDIRDELIM) ;
    std::string rawdir = asciifname.substr(0, dir_index); 
    
    std::size_t Mt_index =  asciifname.rfind("Mt.ascii");
    std::size_t M_index =  asciifname.rfind("M.ascii");
    std::string hdffilename ;
    if (Mt_index == std::string::npos)
    {
        hdffilename = rawdir+std::string(SYSTEMDIRDELIM)+std::string("Mt.h5") ;
    } else {
        hdffilename = rawdir+std::string(SYSTEMDIRDELIM)+std::string("M.h5") ;
    }
    Rcpp::Rcout << " Mt_index: " << Mt_index << std::endl ;
    Rcpp::Rcout << " M_index: "  << M_index << std::endl ;
    Rcpp::Rcout << " start_row: "  << start_row << std::endl ;
    Rcpp::Rcout << " numrows_in_block: "  << numrows_in_block << std::endl ;
    Rcpp::Rcout << " numcols: "  << numcols << std::endl ;
    Rcpp::Rcout << " asciifname  file: " << asciifname << std::endl ;
    Rcpp::Rcout << " hdffilename  file: " << hdffilename << std::endl ;

    
    // reads in data from HDF5 file 
    // to form M Eign double matrix 
    Eigen::MatrixXd  M(numrows_in_block, numcols) ; // This may need to be row major

    CHDF5File<double>* tempHDF5File = (CHDF5File<double>*) new CHDF5File<double>(hdffilename , H5F_ACC_RDWR ) ;  // H5F_ACC_RDONLY
    hid_t datatype_mem = H5T_NATIVE_DOUBLE;  //  H5T_NATIVE_DOUBLE ;  // This is the data type we want returned
    // hid_t datatype_disk = H5T_NATIVE_SCHAR;
    // double * res_data_ptr = NULL ;
    double * M_ptr = &M(0) ;
    size_t num_elements = tempHDF5File->GetSectionOfDatasetReturnInPreallocatedPtr( std::string("/eagle_array"),  M_ptr, start_row, numcols,  datatype_mem) ;    // if we did not swap the M.h5 for the Mt.h5 then numcols would be numrows_in_block
    
   // double * M_ptr = &M(0) ;
    // Close the hdf5 file
    delete tempHDF5File ; 

   // PrintEigenRowCols( hdffilename,  M) ;
    
    return M;

}
#else

Eigen::MatrixXd  ReadBlock(std::string asciifname,
                           long start_row,
                           long numcols,
                           long numrows_in_block)
{
    // reads in data from ASCII file 
    // to form M Eign double matrix 
    std::ostringstream  os;
    std::string line;

    //long coli = 0 ;
    long  rowi = 0;

    // Can we try #define EIGEN_DEFAULT_TO_ROW_MAJOR
    // This does not compile:
    // Eigen::MatrixXd  M(numrows_in_block, numcols, Eigen::RowMajor) ;
    Eigen::MatrixXd  M(numrows_in_block, numcols) ;
    
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
           // int tmp  = line[ii] - '0'; // trick to removes ASCII character offset for numbers
           // M(rowi, ii) = (double) tmp - 1;   // converting data to -1, 0, 1 
            M(rowi, ii) = (double) line[ii] - 49.0;   // converting data to -1, 0, 1 
          }
          rowi++;
      } // end if rr
    } // end for(rr
    
 
    // Close the ascii file
    fileIN.close();

    return M;
}
#endif 


/*
ReadASCIIWriteHDF_rcpp.cpp
*/


void PrintEigenRowCols(std::string mat_filename, Eigen::MatrixXd Min)
{
    int nrowsp =  5;
    int ncolsp = 12;
    int counter = 0 ;
    Rcpp::Rcout << std::endl ;
    Rcpp::Rcout << " First " << nrowsp << " lines and " << ncolsp << " columns of the M matrix  file: " << mat_filename << std::endl ;
   // Rcpp::Rcout << "data size (bytes) = " << num_elements << std::endl ; 

    while(counter < nrowsp)
    {   
       for(int i=0; i < ncolsp ; i++){
           Rcpp::Rcout << Min(counter,i)+1 << " " ;
        }
         Rcpp::Rcout << std::endl ;
         Rcpp::Rcout << std::flush ;
        counter++;
      }  // end  while(getline(fileIN, line ))
    
}


