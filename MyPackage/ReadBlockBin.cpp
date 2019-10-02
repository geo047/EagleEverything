// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

#ifdef _OPENMP
#include <omp.h>
//   [[Rcpp::plugins(openmp)]]
#endif






// [[Rcpp::export]]
Eigen::MatrixXd  ReadBlockBin(std::string binfname,
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
  rowi = 0;

Eigen::MatrixXd
      M(numrows_in_block, numcols) ;

long   numelem = numrows_in_block * numcols;

char * buffer ;
buffer = new char[numelem];


// Binary file
   std::ifstream fileIN(binfname.c_str(), std::ios::in | std::ios::binary);

   if(!fileIN.good()) {
      os << "ERROR: Could not open  " << binfname << std::endl;
      Rcpp::stop(os.str() );
   }

 // set pointer to beginning of block of binary data
 long begofblock = start_row *  numcols;  // number of characters to skip 
 fileIN.seekg(begofblock, std::ios::beg);

 // read block into buffer
 // writing vector to binary file
   fileIN.read( buffer, numelem );
 
  // check 
  Rcpp::Rcout << " this is what is in buffer before doing anything " << std::endl;
  Rcpp::Rcout << buffer[0] << std::endl;
  Rcpp::Rcout << buffer[1] << std::endl;
  Rcpp::Rcout << buffer[2] << std::endl;
  Rcpp::Rcout << buffer[3] << std::endl;
  Rcpp::Rcout << buffer[4] << std::endl;
  Rcpp::Rcout << buffer[5] << std::endl;
  Rcpp::Rcout << buffer[6] << std::endl;
  Rcpp::Rcout << buffer[7] << std::endl;
  Rcpp::Rcout << buffer[8] << std::endl;
  Rcpp::Rcout << buffer[9] << std::endl;

 for(long ii=0; ii < numrows_in_block; ii++){
  for(long jj=0; jj < numcols ; jj++){
    int tmp = buffer[ii] - '0';
    M(start_row + ii , jj) = (double) ((buffer[ ii*numcols + jj ] - '0') - 1);  // converting data into -1, 0, 1 
  } 
 }
  Rcpp::Rcout << M(0,0) << " " << M(0, 1) << " " << M(0,2) << " " << M(0, 3) << " " << M(0,4) << std::endl;
  Rcpp::Rcout << M(1,0) << " " << M(1, 1) << " " << M(1,2) << " " << M(1, 3) << " " << M(1,4) << std::endl;
  Rcpp::Rcout << M(2,0) << " " << M(2, 1) << " " << M(2,2) << " " << M(2, 3) << " " << M(2,4) << std::endl;
  Rcpp::Rcout << M(3,0) << " " << M(3, 1) << " " << M(3,2) << " " << M(3, 3) << " " << M(3,4) << std::endl;
  Rcpp::Rcout << M(4,0) << " " << M(4, 1) << " " << M(4,2) << " " << M(4, 3) << " " << M(4,4) << std::endl;

//   for(long rr=0; rr < (start_row + numrows_in_block) ; rr++){
//      // read a line of data from ASCII file
//      getline(fileIN, line);
//      if(rr >= start_row){
//          std::istringstream streamA(line);
//          for(long ii=0; ii < numcols  ; ii++){
//            int tmp  = line[ii] - '0'; // trick to removes ASCII character offset for numbers
//            M(rowi, ii) = (double) tmp - 1;   // converting data to -1, 0, 1 
//          }
//          rowi++;
//      } // end if rr
//   } // end for(rr



// Close the bin file
   fileIN.close();

delete [] buffer;
buffer = NULL;

 return M;

}



