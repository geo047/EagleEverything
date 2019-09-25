#include <Rcpp.h>
#include<magma_v2.h>


/* Author: Andrew W. George
   Date:   16 Sep, 2019
   Purpose: to implement a multi-GPU version of the R function qr() which is very slow. 
*/




//--------------------------------------------
// Magma test code 
//--------------------------------------------
// [[Rcpp::export]]
 int  magma_solve(Rcpp::NumericMatrix  X , int numgpus, bool printInfo, std::string fname   )
{

  /*
  - A square symmetric data matrix is assumed. 
  - The data is read in from R and a copy is made otherwise the contents are overridden. 
  - The interface between R and Magma is 32 bits, even when Magma and R are built as 64 bit code. 
    Josh got around this by designing a server to run on shared memory. However, shared memory is limited to 
    matrices of around 46000 rows/cols. 
    I got around the problem by writing the inverse to disk, as a binary file, and reading this binary file back into R (very fast).   
    A simple solution to a complicated problem. 
  */

   magma_init();


  if (printInfo){
   magma_print_environment();
  }






  magma_int_t n, n2, lwork, info = 0;  // define MAGMA_ILP64 to get these as 64 bit integers
   n = X.rows() ;   // this is X.rows() but passed in from command line
   n2     = n*n;


  double * h_X;
  if (  MAGMA_SUCCESS !=   magma_dmalloc_cpu(&h_X, n2) )
  {
    std::cout << " Error: magma_eigen() magma_dmalloc_cpu failed for h_X. Need more memory" << std::endl;
    return -1;
  }


 // Assign data to CPU 
  for (int j=0; j< n; j++){
    for (int i=0; i < n; i++){
        h_X[i+n*j] = X(i,j);
    }
  }





  if (printInfo){
     std::cout << " Prior to magma_dpotrf: First 5 rows and columns of the  X matrix " << std::endl;
     magma_dprint(5,5, h_X, n);
  }

  if (numgpus == 1){
      magma_dpotrf( MagmaLower, n, h_X, n, &info );
  
      if ( info != 0){
        std::cout << "Error: magma_dpotrf has returned a non-zero info value of " << info << std::endl;
        std::cout << "       This error is most likely due to the matrix beig too large for the GPU. " << std::endl;
        std::cout << "        Need to use the CPU-based solve function in R instead. " << std::endl;         
        std::cout << "        Select this option in the AM function.      " << std::endl;
        return info;
      }
   } else {
       magma_dpotrf_m (numgpus,  MagmaLower, n, h_X, n, &info ); 

      if ( info != 0){
        std::cout << "Error: magma_dpotrf has returned a non-zero info value of " << info << std::endl;
        std::cout << "       This error is most likely due to the matrix beig too large for the GPU. " << std::endl;
        std::cout << "        Need to use the CPU-based solve function in R instead. " << std::endl;
        std::cout << "        Select this option in the AM function.      " << std::endl;
        return info;
      }


 
   }


  if (printInfo){
     std::cout << " After magma_dpotrf but prior to magma_dpotri: First 5 rows and columns of the  X matrix " << std::endl;
     magma_dprint(5,5, h_X, n);
  }


   magma_dpotri( MagmaLower, n, h_X, n, &info );

      if ( info != 0){
        std::cout << "Error: magma_dpotri has returned a non-zero info value of " << info << std::endl;
        std::cout << "       This error is most likely due to the matrix beig too large for the GPU. " << std::endl;
        std::cout << "        Need to use the CPU-based solve function in R instead. " << std::endl;
        std::cout << "        Select this option in the AM function.      " << std::endl;
        return info;
      }



  if (printInfo){
     std::cout << " After magma calls: First 5 rows and columns of the INVERSE  matrix " << std::endl;
     magma_dprint(5,5, h_X, n);
  }

 //Move lower to upper triangle
  for (int j=0; j< n; j++){
    for (int i=j; i < n; i++){
        // h_X[i+n*j] = h_X[j+n*i];
        h_X[j+n*i] = h_X[i+n*j];
    }
  }



   if( printInfo){
       std::cout << " About to write inverse to binary file. " << std::endl;
    }



    if (fname.empty()){
      std::cout << " Error: in magma_eigen. Name of binary file for the eigenvectors has not been specified. \n" << std::endl;
      return -1;
    }

    FILE* file = fopen(fname.c_str() , "wb");
    if (file == NULL){
       std::cout << " Error: in magma_eigen. Binary file for eigenvectors has failed to open." << std::endl;
       return -1;
    } else {
       fwrite(&h_X[0], 1 , n*n*sizeof(double) , file);
       fclose(file);
    }
 

  // tidying up 
  magma_free_cpu (h_X)  ;

 magma_finalize();
  
 return 0;


}
                       
  


