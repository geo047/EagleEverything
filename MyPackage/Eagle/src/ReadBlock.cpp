// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#ifdef _OPENMP
#include <omp.h>
//   [[Rcpp::plugins(openmp)]]
#endif

void PrintEigenRowCols(std::string mat_filename, Eigen::MatrixXd Min) ;



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


