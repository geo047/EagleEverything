#include <Rcpp.h>
#include<magma_v2.h>
#include<magma_lapack.h>
#include<magma_operators.h>

#include <fstream>
#include <chrono>
#include <vector>
#include <cstdint>
#include <numeric>
#include <random>
#include <algorithm>
#include <iostream>
#include <cassert>




using namespace Rcpp;

// Potential issue
// First, you are accessing the matrix row-wise, whereas MAGMA and LAPACK use column-wise ordering. That is,
//
// A[ i + j*size ]


/***************************************************************************//**
 * Macros 
 */

#define PRINT(x)  \
     std::cout << x << std::endl

#define CHECK_MALLOC( err )                                                                        \
    if (err != MAGMA_SUCCESS){                                                                     \
        Rcpp::Rcout << "\nError: Eigenvalue calculation has failed due to failure of memory allocation.\n" << std::endl;  \
     }

#define CHECK_GPU( err )   \
   if (err != MAGMA_SUCCESS ){                         \
       Rcpp::Rcout <<"\n" << std::endl;    \
       Rcpp::Rcout << "\nError:    " << std::endl;    \
       Rcpp::Rcout << "  Eigenvalue calculation has failed. " << std::endl;  \
       Rcpp::Rcout << "  Things to try to solve this problem are:            " << std::endl;  \
       Rcpp::Rcout << "    1. Reduce the number of individuals to see if this is a GPU memory issue. " << std::endl; \
       Rcpp::Rcout << "    2. Ensure that you do not have perfectly correlated fixed effects in the model. " << std::endl; \
       Rcpp::Rcout << "       This will cause collinearity between columns in the model matrix. \n\n " << std::endl;     \
   }


#define SHUTDOWN   \
     magma_free_cpu ( work ); \
     magma_free_cpu ( h_R ); \
     magma_free_cpu (tau); \
     magma_queue_destroy ( queue ); \
     magma_finalize (); 





// [[Rcpp::export]]
Rcpp::List   gpuQR_magma(const Rcpp::NumericMatrix  X)
{


   // convert to C++ type
   double const* d_X = X.begin();   // this is a column-wise vector which magma likes
  


  // Initialize magma and cublas
   magma_init();

   magma_print_environment();

   // Initialize the queue
   magma_queue_t queue = NULL ;
   magma_int_t dev =0;
   magma_queue_create (dev ,& queue );


   double  *work, *h_R, *tau ;
   magma_int_t M, N, n2, lda, ldda,  info,  nb;
   magma_int_t ngpu=4;


   M = X.rows();
   N = X.cols();
   lda    = M;
   n2     = lda*N;
   ldda = ((M +31)/32)*32; // ldda = m if 32 divides m
   nb     = magma_get_dgeqrf_nb( M, N );
   magma_int_t lwork = N*nb;


   // mallocs
   magma_dmalloc_cpu(&work, lwork);
   magma_dmalloc_cpu(&tau, N);
   magma_dmalloc_cpu(&h_R, n2);


 // Assign data to CPU 
  for (int j=0; j< N; j++){
    for (int i=0; i < N; i++){
        h_R[i+N*j] = X(i,j);
    }
  }


// magma_dprint(5,5, h_R, N); 


  /* ====================================================================
               Performs operation using MAGMA
  =================================================================== */

  magma_dgeqrf_m( ngpu, M, N, h_R , N   , tau, work, lwork, &info );
  if (info != 0) {
       printf("magma_dgeqrf_m returned error %lld: %s.\n",
                       (long long) info, magma_strerror( info ));
       }

    
// create R list structure
NumericVector values =  NumericVector(tau,tau + N  );
NumericVector vectors = NumericVector(h_R, h_R + n2);
vectors.attr("dim") = Dimension(N,N);


for(int i=0; i< 3; i++){
  for(int j=0; j< 3; j++){
    std::cout << h_R[i+N*j] << " " ;
  }
  std::cout << std::endl;
}



// AWG



std::cout << "Staring writing fwrite ... " << std::endl;
auto startTime2 = std::chrono::high_resolution_clock::now();
FILE* file = fopen("/scratch1/geo047/test.bin", "wb");
    fwrite(&h_R[0], 1 , N*N*sizeof(double) , file);
    fclose(file);
auto endTime2 = std::chrono::high_resolution_clock::now();
std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(endTime2 - startTime2).count() << std::endl;



std::cout << "Staring writing myfile <<  ... " << std::endl;
auto startTime = std::chrono::high_resolution_clock::now();
std::ofstream myFile ("/scratch1/geo047/test.dat");
if (myFile.is_open()){
 for(int ii=0; ii < N; ii++){
   for(int jj=0; jj < N; jj++){
      // myFile << h_R[ii+N*jj] << std::setprecision(4) << " ";
      myFile << h_R[ii+N*jj] << " ";
   }
   myFile << std::endl;
 }
}
myFile.close();
auto endTime = std::chrono::high_resolution_clock::now();
std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count() << std::endl;



// clean up and shut down 
 SHUTDOWN;

// Return list structure 
return Rcpp::List::create(
       Rcpp::Named("tau") =   values,
      Rcpp::Named("A") =  vectors );




}


