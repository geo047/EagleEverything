#include <Rcpp.h>
#include <stdio.h>
#include <stdlib.h>
#include "cublas_v2.h"
#include "magma_v2.h"      // also includes cublas_v2.h
#include "magma_lapack.h"  // if you need BLAS & LAPACK
#include<magma_operators.h>

// set to 1 if debug information is desired
#define DEBUG 1


using namespace Rcpp;

/***************************************************************************//**
 * Macros 
 */

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
     magma_free_cpu (w1 ); \
     magma_free_cpu (r); \
     magma_free_cpu ( h_work ); \
     magma_free_cpu ( h_r ); \
     magma_free (d_r ); \
     magma_queue_destroy ( queue ); \
     magma_finalize (); \
     delete []  iwork






// [[Rcpp::export]]
Rcpp::List   gpuEigen_magma(const Rcpp::NumericMatrix&  X)
{

// Initialize Magma and print GPU environment
magma_init (); 
if (DEBUG)
  magma_print_environment();


// Initialize the queue
magma_queue_t queue = NULL ;
magma_int_t dev =0;
magma_queue_create (dev ,& queue );


// Variables 
double gpu_time , cpu_time ;
magma_int_t n=X.nrow() , n2=n*n;
magma_int_t  ldda = magma_roundup( n, 32 );


double *r; //  r - nxn matrices on the host
double *h_r ;   // nxn matrix with eigenvector on host
double *d_r ; // nxn matrix on the device
double * h_work ; // workspace
magma_int_t lwork ; // h_work size
magma_int_t * iwork ; // workspace
magma_int_t liwork ; // iwork size
double *w1  ; // w1  - vectors of eigenvalues
double error , work [1]; // used in difference computations
magma_int_t ione = 1, info ;
double mione = -1.0;
magma_int_t incr = 1;
magma_int_t ISEED [4] = {0 ,0 ,0 ,1}; // seed
CHECK_MALLOC(magma_dmalloc_cpu (&w1 ,n)); // host memory for real
CHECK_MALLOC(magma_dmalloc_cpu (&r,n2 )); // host memory for r
CHECK_MALLOC(magma_dmalloc_cpu (&h_r,n2 )); // host memory for r
CHECK_MALLOC(magma_dmalloc (& d_r ,n2 )); // device memory for d_r


// Query for workspace sizes
double aux_work [1];
magma_int_t aux_iwork [1];
CHECK_GPU(magma_dsyevd_gpu ( MagmaVec , MagmaLower ,n,NULL ,n,NULL ,NULL,n, aux_work ,
               -1, aux_iwork ,-1,& info ));
if (info != MAGMA_SUCCESS){
     // clean up and shut down 
     SHUTDOWN;
     return 0;
};



lwork = ( magma_int_t ) aux_work [0];
liwork = aux_iwork [0];
iwork =( magma_int_t *) malloc ( liwork * sizeof ( magma_int_t ));
magma_dmalloc_cpu (& h_work , lwork ); // memory for workspace



// Convert X into a form that Magma can deal with of type double *
// Have to do it in stages
// double const *xx = &X(0,0);
double const *xx = X.begin();

// Copy contents of X that sit in CPU memory to GPU memory : xx -> d_r  
magma_dsetmatrix ( n, n, xx  , n, d_r ,n , queue ); 
if (DEBUG){
         Rcpp::Rcout << " Contents of first 10 rows and columns of data on GPU  " << std::endl;
         magma_dprint_gpu(10,10, d_r, n , queue);
};



// Doing eigenvalue calculation in GPU land
CHECK_GPU( magma_dsyevd_gpu(MagmaVec,MagmaLower,n,d_r,n,w1,r,n,h_work,
 lwork,iwork,liwork,&info));
if (info != MAGMA_SUCCESS){
    // clean up and shut down 
     SHUTDOWN;
     return 0;
};
if (DEBUG){
     Rcpp::Rcout << " First 10 eigenvalues (descending order)... " << std::endl;
     for(int i=0; i<10; i++){
        Rcpp::Rcout << w1[i] << " " ;
     }
     Rcpp::Rcout << "\n\n" << std::endl;

     Rcpp::Rcout << " First 10 rows and columns of the eigenvector matrix (residing in GPU memory) " << std::endl;
     magma_dprint_gpu(10,10,d_r,n, queue);
};





// Move eigenvector matrix from GPU (d_r) to CPU memory (h_r)
magma_dgetmatrix(n,n,d_r, n, h_r, n, queue);
if (DEBUG){
     Rcpp::Rcout << " First 10 rows and columns of the eigenvector matrix (residing in CPU memory) " << std::endl;
     Rcpp::Rcout << " Should be the same as what is in GPU memory (above) " << std::endl;
     magma_dprint(10,10,h_r,n) ;
};




// create R list structure
NumericVector values =  NumericVector(w1,w1 + n  );
NumericVector vectors = NumericVector(h_r, h_r + n2);
vectors.attr("dim") = Dimension(n,n);


// clean up and shut down 
SHUTDOWN;




// Return list structure 
return Rcpp::List::create( 
       Rcpp::Named("values") =   values, 
       Rcpp::Named("vectors") =  vectors );

} 

